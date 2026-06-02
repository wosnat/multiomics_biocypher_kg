#!/usr/bin/env python3
"""
release_kg.py — Agnostic-core release orchestrator for the multi-omics KG.

Six phases (all idempotent):
  1. Preflight       — version regex, tooling, gh auth, git state, .env present
  2. CHANGELOG cut   — rename [Unreleased] → [<version>] - <date>; PAUSE for polish
  3. Commit/tag/push — chore(release): kg-<version> + annotated tag + push --follow-tags
  4. Clean clone     — git clone --branch kg-<version> --depth 1 → /tmp/...
  5. Build + verify  — docker compose -p kg-release-staging up --build -d; query Schema_info
  6. Deploy + publish — `--target staging` keeps staging up; gh release create + metadata.json

`--dry-run` exercises every phase without mutating anything.
See SKILL.md and references/PHASES.md for the full design.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

# ─── Constants ──────────────────────────────────────────────────────────────
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+(-(alpha|beta|rc)\.\d+)?$")
DEFAULT_MCP_MIN = "0.1.0"
DEFAULT_TARGET = "staging"

STAGING_PROJECT = "kg-release-staging"
STAGING_HTTP_BIND = "127.0.0.1:27474"
STAGING_BOLT_BIND = "127.0.0.1:27687"
STAGING_BOLT_URL = "bolt://localhost:27687"
STAGING_DEPLOY_CONTAINER = "staging-deploy"  # set by docker-compose.staging.yml
STAGING_COMPOSE_OVERRIDE = "docker-compose.staging.yml"
# Phase 6 queries the staging deploy *from inside the container* (`docker exec
# staging-deploy cypher-shell …`), so we connect to Bolt over the container's
# own loopback on port 7687 — not the host-published 27687. This removes the
# host-`cypher-shell` dependency.
STAGING_INNER_BOLT_URL = "bolt://localhost:7687"

CLONE_PARENT = Path("/tmp")
COMMIT_PREFIX = "chore(release): kg-"
TAG_PREFIX = "kg-"


# ─── Logging ────────────────────────────────────────────────────────────────
def log(msg: str = "") -> None:
    print(msg, flush=True)


def section(title: str) -> None:
    log(f"\n=== {title} ===")


def dry_msg(action: str) -> None:
    log(f"[dry-run] would {action}")


def die(msg: str) -> None:
    log(f"\nFATAL: {msg}")
    sys.exit(1)


# ─── Context ────────────────────────────────────────────────────────────────
@dataclass
class Context:
    version: str
    target: str
    mcp_min: str
    allow_dirty: bool
    draft: bool
    dry_run: bool
    resume: bool
    # Filled by preflight
    git_sha: str = ""
    git_sha_short: str = ""
    git_branch: str = ""
    git_dirty: bool = False
    # Filled by clean_clone
    clone_dir: Optional[Path] = None
    # Filled by verify
    schema_info: dict = field(default_factory=dict)


# ─── Subprocess helpers ─────────────────────────────────────────────────────
def run(cmd: list[str], *, check: bool = True, capture: bool = False,
        cwd: Optional[Path] = None, env: Optional[dict] = None,
        show: bool = True) -> subprocess.CompletedProcess:
    if show:
        log(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd, check=check, capture_output=capture, text=True, cwd=cwd, env=env)


def run_or_dry(ctx: Context, action: str, cmd: list[str], **kwargs):
    if ctx.dry_run:
        dry_msg(action)
        log(f"  [would run] {' '.join(cmd)}")
        return None
    return run(cmd, **kwargs)


def git_out(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True).stdout.strip()


# ─── Phase 1: Preflight ─────────────────────────────────────────────────────
def phase_preflight(ctx: Context) -> None:
    section("Phase 1: Preflight")

    if not VERSION_RE.match(ctx.version):
        die(f"version {ctx.version!r} does not match X.Y.Z[-(alpha|beta|rc).N]")
    log(f"  version: {ctx.version} ✓")

    # `git` is needed even in --dry-run (we capture SHA/branch via subprocess).
    # The other three are only invoked in a real run, so their checks soften under --dry-run.
    if not shutil.which("git"):
        die("required tool not on PATH: git")
    if ctx.dry_run:
        log(f"  tooling: git ✓  (dry-run: skipping docker/gh availability + gh auth checks)")
    else:
        # cypher-shell is NOT a host requirement — Phase 6 invokes it via
        # `docker exec staging-deploy cypher-shell …`, which the neo4j image
        # ships at /var/lib/neo4j/bin/cypher-shell.
        for tool in ["docker", "gh"]:
            if not shutil.which(tool):
                die(f"required tool not on PATH: {tool} (only `git` is required for --dry-run)")
        log(f"  tooling: docker / git / gh ✓")

        auth = subprocess.run(["gh", "auth", "status"], capture_output=True, text=True)
        if auth.returncode != 0:
            die(f"gh not authenticated:\n{auth.stderr.strip() or auth.stdout.strip()}")
        log(f"  gh auth: ✓")

    if not Path(".env").exists():
        die(".env missing in repo root — build.sh hardcodes `cp /src/.env .` and will abort")
    log(f"  .env: present ✓")

    branch = git_out("rev-parse", "--abbrev-ref", "HEAD")
    sha = git_out("rev-parse", "HEAD")
    sha_short = git_out("rev-parse", "--short", "HEAD")
    if not branch or not sha:
        die("not in a git repo, or HEAD is detached")

    status = git_out("status", "--porcelain")
    dirty = bool(status)
    # `--resume` enters with CHANGELOG.md modified by Phase 2's cut — that's the
    # intentional handoff state, not real drift. Allow exactly that single-file
    # diff; anything else means the operator polished beyond CHANGELOG (or there
    # is unrelated drift) and we should still surface it.
    # NB: parse with split(maxsplit=1) because git_out's .strip() drops the
    # porcelain status's leading whitespace column for unstaged changes, so
    # the file path is no longer at a fixed offset.
    modified_files: set[str] = set()
    for raw_line in status.splitlines():
        parts = raw_line.strip().split(maxsplit=1)
        if len(parts) == 2:
            modified_files.add(parts[1])
    only_changelog_cut = ctx.resume and modified_files == {"CHANGELOG.md"}
    if dirty and not ctx.allow_dirty and not only_changelog_cut:
        die(f"working tree dirty (re-run with --allow-dirty to override):\n{status}")

    # In sync with origin (best-effort — don't fail if offline)
    subprocess.run(["git", "fetch", "origin", branch],
                   check=False, capture_output=True, text=True)
    behind_raw = git_out("rev-list", "--count", f"HEAD..origin/{branch}")
    behind = int(behind_raw) if behind_raw.isdigit() else 0
    if behind > 0 and not ctx.allow_dirty:
        die(f"branch {branch} is {behind} commits behind origin/{branch}; "
            f"pull first or use --allow-dirty")

    ctx.git_sha = sha
    ctx.git_sha_short = sha_short
    ctx.git_branch = branch
    # On --resume, the CHANGELOG-only diff is the cut's intentional product,
    # not real drift, so it does not flip the release's git_dirty stamp.
    ctx.git_dirty = dirty and not only_changelog_cut
    log(f"  branch: {branch} @ {sha_short} (dirty={ctx.git_dirty}, behind=0)")


# ─── CHANGELOG cut helper ───────────────────────────────────────────────────
def cut_changelog(content: str, version: str, today: str) -> tuple[str, bool]:
    """
    Cut [Unreleased] → [version] - today. Returns (new_content, performed).
    Idempotent: if [version] already present, returns (content, False).
    """
    lines = content.splitlines()
    version_header = f"## [{version}]"
    if any(l.startswith(version_header) for l in lines):
        return content, False

    try:
        unreleased_idx = next(i for i, l in enumerate(lines) if l.strip() == "## [Unreleased]")
    except StopIteration:
        raise ValueError("CHANGELOG.md missing `## [Unreleased]` section")

    end_idx = len(lines)
    for i in range(unreleased_idx + 1, len(lines)):
        if lines[i].startswith("## ["):
            end_idx = i
            break

    unreleased_body = "\n".join(lines[unreleased_idx + 1:end_idx]).strip()
    has_entries = any(l.lstrip().startswith("-") for l in lines[unreleased_idx + 1:end_idx])

    new_unreleased = [
        "## [Unreleased]",
        "",
        "### Added",
        "",
        "### Changed",
        "",
        "### Fixed",
        "",
    ]
    new_version_section = [
        f"## [{version}] - {today}",
        "",
        unreleased_body if has_entries else "_No entries — placeholder._",
        "",
    ]

    new_lines = lines[:unreleased_idx] + new_unreleased + new_version_section + lines[end_idx:]
    new_content = "\n".join(new_lines).rstrip() + "\n"
    return new_content, True


# ─── Phase 2: CHANGELOG cut ─────────────────────────────────────────────────
def phase_changelog_cut(ctx: Context) -> bool:
    section("Phase 2: CHANGELOG cut")

    path = Path("CHANGELOG.md")
    if not path.exists():
        die("CHANGELOG.md missing in repo root")

    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    try:
        new_content, performed = cut_changelog(path.read_text(), ctx.version, today)
    except ValueError as e:
        die(str(e))

    if not performed:
        log(f"  [{ctx.version}] section already exists — cut already done, skipping")
        return False

    if ctx.dry_run:
        dry_msg(f"rename [Unreleased] → [{ctx.version}] - {today} in CHANGELOG.md")
        return True

    path.write_text(new_content)
    log(f"  cut: [Unreleased] → [{ctx.version}] - {today}")
    return True


# ─── Phase 3: Commit + tag + push ───────────────────────────────────────────
def phase_commit_tag_push(ctx: Context) -> None:
    section("Phase 3: Commit + tag + push")
    commit_subject = f"{COMMIT_PREFIX}{ctx.version}"
    tag_name = f"{TAG_PREFIX}{ctx.version}"

    head_subject = git_out("log", "-1", "--pretty=%s")
    if head_subject == commit_subject:
        log(f"  HEAD already at {commit_subject!r} — skipping commit")
    else:
        run_or_dry(ctx, f"stage CHANGELOG.md", ["git", "add", "CHANGELOG.md"])
        run_or_dry(ctx, f"commit {commit_subject!r}",
                   ["git", "commit", "-m", commit_subject])

    existing = git_out("tag", "-l", tag_name)
    if existing:
        log(f"  tag {tag_name} already exists — skipping")
    else:
        run_or_dry(ctx, f"create annotated tag {tag_name}",
                   ["git", "tag", "-a", tag_name, "-m", f"KG release {ctx.version}"])

    run_or_dry(ctx, "push branch + tag",
               ["git", "push", "origin", ctx.git_branch, "--follow-tags"])


# ─── Phase 4: Clean clone of the tag ────────────────────────────────────────
def phase_clean_clone(ctx: Context) -> None:
    section(f"Phase 4: Clean clone of tag {TAG_PREFIX}{ctx.version}")
    clone_dir = CLONE_PARENT / f"kg-release-{ctx.version}"
    ctx.clone_dir = clone_dir

    if clone_dir.exists():
        log(f"  {clone_dir} exists — removing for clean clone")
        if not ctx.dry_run:
            shutil.rmtree(clone_dir)

    origin_url = git_out("remote", "get-url", "origin")
    if not origin_url:
        die("could not determine origin URL via `git remote get-url origin`")

    tag = f"{TAG_PREFIX}{ctx.version}"
    run_or_dry(ctx, f"clone tag {tag} → {clone_dir}",
               ["git", "clone", "--branch", tag, "--depth", "1", origin_url, str(clone_dir)])

    if not ctx.dry_run:
        shutil.copy(Path(".env"), clone_dir / ".env")
        log(f"  copied .env → {clone_dir}/.env (gitignored; build.sh needs it)")
    else:
        dry_msg(f"copy .env → {clone_dir}/.env")


# ─── Phase 5: Build into staging + verify ───────────────────────────────────
def phase_build_and_verify(ctx: Context) -> None:
    section(f"Phase 5: Build into staging ({STAGING_PROJECT}) + verify")

    release_env = os.environ.copy()
    release_env.update({
        "KG_RELEASE_VERSION": ctx.version,
        "KG_GIT_SHA": ctx.git_sha,
        "KG_GIT_SHA_SHORT": ctx.git_sha_short,
        "KG_GIT_BRANCH": ctx.git_branch,
        "KG_GIT_DIRTY": "true" if ctx.git_dirty else "false",
        "KG_MCP_MIN_VERSION": ctx.mcp_min,
        "KG_DEPLOY_HTTP_BIND": STAGING_HTTP_BIND,
        "KG_DEPLOY_BOLT_BIND": STAGING_BOLT_BIND,
    })

    compose_up = [
        "docker", "compose",
        "-p", STAGING_PROJECT,
        "-f", "docker-compose.yml",
        "-f", STAGING_COMPOSE_OVERRIDE,  # renames containers staging-* so dev's `deploy` etc. don't collide
        "up", "--build", "-d",
        "build", "import", "post-process", "deploy",
    ]
    if ctx.dry_run:
        dry_msg(f"`docker compose -p {STAGING_PROJECT} up --build` from {ctx.clone_dir}")
        log(f"  [would run] {' '.join(compose_up)}")
        log(f"  [would set] KG_RELEASE_VERSION={ctx.version} "
            f"KG_GIT_SHA_SHORT={ctx.git_sha_short} "
            f"KG_DEPLOY_BOLT_BIND={STAGING_BOLT_BIND}")
        dry_msg(f"query Schema_info via `docker exec {STAGING_DEPLOY_CONTAINER} "
                f"cypher-shell -a {STAGING_INNER_BOLT_URL}`; assert version == {ctx.version}")
        ctx.schema_info = {"version": ctx.version, "git_sha": ctx.git_sha,
                           "papers": -1, "experiments": -1, "genes": -1,
                           "organisms": -1, "expr_edges": -1,
                           "mcp_min": ctx.mcp_min, "built_at": "<dry-run>"}
        return

    run(compose_up, cwd=ctx.clone_dir, env=release_env)
    log(f"  staging stack up; waiting for {STAGING_DEPLOY_CONTAINER} to accept Bolt "
        f"(via docker exec) …")

    # `docker exec` requires the container to be running; depends_on chains
    # guarantee it exists by the time compose-up returns, but Bolt warmup
    # takes a few seconds. Poll up to 120s.
    cypher_shell_in_container = [
        "docker", "exec", STAGING_DEPLOY_CONTAINER,
        "cypher-shell", "-a", STAGING_INNER_BOLT_URL,
    ]
    for i in range(120):
        ready = subprocess.run(cypher_shell_in_container + ["RETURN 1;"],
                               capture_output=True, text=True)
        if ready.returncode == 0:
            log(f"  Bolt ready after {i}s")
            break
        time.sleep(1)
    else:
        die(f"staging deploy not reachable via `docker exec {STAGING_DEPLOY_CONTAINER} "
            f"cypher-shell -a {STAGING_INNER_BOLT_URL}` after 120s")

    cypher = (
        "MATCH (s:Schema_info {id:'schema_info'}) "
        "RETURN s.version AS version, s.git_sha AS git_sha, "
        "s.paper_count AS papers, s.experiment_count AS experiments, "
        "s.gene_count AS genes, s.organism_count AS organisms, "
        "s.expression_edge_count AS expr_edges, "
        "s.mcp_min_version AS mcp_min, s.built_at AS built_at"
    )
    result = run(cypher_shell_in_container + ["--format", "plain", cypher],
                 capture=True, show=False)
    ctx.schema_info = _parse_plain_row(result.stdout)

    if ctx.schema_info.get("version") != ctx.version:
        die(f"Schema_info.version = {ctx.schema_info.get('version')!r} but expected {ctx.version!r}")
    log(f"  version: {ctx.schema_info['version']} ✓ (matches tag)")
    log(f"  git_sha: {ctx.schema_info['git_sha']}")
    log(f"  built_at: {ctx.schema_info['built_at']}")
    log(f"  counts: {ctx.schema_info['papers']} papers · "
        f"{ctx.schema_info['experiments']} experiments · "
        f"{ctx.schema_info['genes']} genes · "
        f"{ctx.schema_info['organisms']} organisms · "
        f"{ctx.schema_info['expr_edges']} expression edges")
    log(f"  explorer smoke test: out-of-scope here (explorer-repo work; verify manually if needed)")


def _parse_plain_row(stdout: str) -> dict:
    """Parse cypher-shell --format plain: header line + value line, comma-separated."""
    lines = [l for l in stdout.strip().splitlines() if l.strip()]
    if len(lines) < 2:
        die(f"Schema_info query returned no data:\n{stdout}")
    headers = [h.strip() for h in lines[0].split(",")]
    raw_values = [v.strip() for v in lines[1].split(",")]
    out: dict = {}
    for h, v in zip(headers, raw_values):
        v = v.strip()
        if v.startswith('"') and v.endswith('"'):
            v = v[1:-1]
        if v.lstrip("-").isdigit():
            out[h] = int(v)
        else:
            out[h] = v
    return out


# ─── Phase 6: Deploy ────────────────────────────────────────────────────────
def deploy_staging(ctx: Context) -> None:
    log(f"  --target staging: leaving the staging stack up on {STAGING_BOLT_URL}")
    log(f"  to tear down later: docker compose -p {STAGING_PROJECT} down")


def deploy_local_stub(ctx: Context) -> None:
    raise NotImplementedError(
        "--target local is not yet implemented. Track A (local-box alpha) is gated on "
        "the hosting decision (local box vs Aura). See plans/alpha_release.md §2.2-2.6."
    )


def deploy_aura_stub(ctx: Context) -> None:
    raise NotImplementedError(
        "--target aura is not yet implemented. Track B (Aura) is gated on the hosting "
        "decision (local box vs Aura). See plans/alpha_release.md §7.3."
    )


DEPLOY_BACKENDS = {
    "staging": deploy_staging,
    "local": deploy_local_stub,
    "aura": deploy_aura_stub,
}


def phase_deploy(ctx: Context) -> None:
    section(f"Phase 6: Deploy (--target {ctx.target})")
    DEPLOY_BACKENDS[ctx.target](ctx)


# ─── Phase 7: Publish ───────────────────────────────────────────────────────
def build_metadata(ctx: Context) -> dict:
    return {
        "version": ctx.version,
        "tag": f"{TAG_PREFIX}{ctx.version}",
        "git_sha": ctx.git_sha,
        "git_sha_short": ctx.git_sha_short,
        "git_branch": ctx.git_branch,
        "git_dirty": ctx.git_dirty,
        "mcp_min_version": ctx.mcp_min,
        "target": ctx.target,
        "built_at": ctx.schema_info.get("built_at", ""),
        "counts": {
            "papers": ctx.schema_info.get("papers"),
            "experiments": ctx.schema_info.get("experiments"),
            "genes": ctx.schema_info.get("genes"),
            "organisms": ctx.schema_info.get("organisms"),
            "expression_edges": ctx.schema_info.get("expr_edges"),
        },
        "stamped_at": datetime.now(timezone.utc).isoformat(),
    }


def extract_changelog_fragment(path: Path, version: str) -> str:
    if not path.exists():
        die(f"CHANGELOG.md missing at {path}")
    lines = path.read_text().splitlines()
    version_header = f"## [{version}]"
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith(version_header))
    except StopIteration:
        die(f"CHANGELOG.md missing `[{version}]` section — Phase 2 cut not done?")
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith("## ["):
            end = i
            break
    return "\n".join(lines[start:end]).rstrip() + "\n"


def phase_publish(ctx: Context) -> None:
    section("Phase 7: Publish GitHub Release")
    base = ctx.clone_dir or Path(".")
    metadata = build_metadata(ctx)
    metadata_path = base / "metadata.json"
    fragment_path = base / f".release-notes-{ctx.version}.md"

    if ctx.dry_run:
        dry_msg(f"write {metadata_path}")
        log(f"  [would write metadata.json]")
        for line in json.dumps(metadata, indent=2).splitlines():
            log(f"    {line}")
        dry_msg(f"extract CHANGELOG [{ctx.version}] section → {fragment_path}")
        gh_cmd = ["gh", "release", "create", f"{TAG_PREFIX}{ctx.version}",
                  "--notes-file", str(fragment_path), "--prerelease"]
        if ctx.draft:
            gh_cmd.append("--draft")
        dry_msg(f"create GitHub Release {TAG_PREFIX}{ctx.version}")
        log(f"  [would run] {' '.join(gh_cmd)}")
        dry_msg("upload metadata.json to the release")
    else:
        metadata_path.write_text(json.dumps(metadata, indent=2))
        log(f"  metadata.json: {metadata_path}")
        fragment = extract_changelog_fragment(base / "CHANGELOG.md", ctx.version)
        fragment_path.write_text(fragment)
        log(f"  release-notes fragment: {fragment_path}")
        gh_cmd = ["gh", "release", "create", f"{TAG_PREFIX}{ctx.version}",
                  "--notes-file", str(fragment_path), "--prerelease"]
        if ctx.draft:
            gh_cmd.append("--draft")
        run(gh_cmd, cwd=base)
        run(["gh", "release", "upload", f"{TAG_PREFIX}{ctx.version}", str(metadata_path)],
            cwd=base)

    log("")
    log("  Operator checklist:")
    log(f"    - GitHub Release: gh release view {TAG_PREFIX}{ctx.version}")
    log(f"    - Staging URI:    {STAGING_BOLT_URL}  (--target staging keeps it up)")
    log(f"    - Tear down staging:  docker compose -p {STAGING_PROJECT} down")
    if ctx.clone_dir:
        log(f"    - Clean up clone: rm -rf {ctx.clone_dir}")


# ─── Entry point ────────────────────────────────────────────────────────────
def parse_args(argv: Optional[list[str]] = None) -> Context:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("version", help="X.Y.Z[-(alpha|beta|rc).N]")
    ap.add_argument("--target", default=DEFAULT_TARGET,
                    choices=list(DEPLOY_BACKENDS.keys()))
    ap.add_argument("--mcp-min", default=DEFAULT_MCP_MIN)
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--draft", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="Skip the post-CHANGELOG-cut pause (use on the 2nd invocation)")
    args = ap.parse_args(argv)
    return Context(
        version=args.version,
        target=args.target,
        mcp_min=args.mcp_min,
        allow_dirty=args.allow_dirty,
        draft=args.draft,
        dry_run=args.dry_run,
        resume=args.resume,
    )


def main(argv: Optional[list[str]] = None) -> int:
    ctx = parse_args(argv)
    suffix = ", --dry-run" if ctx.dry_run else ""
    log(f"=== release-kg: {ctx.version} (--target {ctx.target}{suffix}) ===")

    phase_preflight(ctx)
    cut_performed = phase_changelog_cut(ctx)

    if cut_performed and not ctx.resume and not ctx.dry_run:
        log("")
        log("=== PAUSE: CHANGELOG cut done ===")
        log("Review and polish CHANGELOG.md as needed, then re-run with --resume to continue.")
        return 0

    phase_commit_tag_push(ctx)
    phase_clean_clone(ctx)
    phase_build_and_verify(ctx)
    phase_deploy(ctx)
    phase_publish(ctx)

    log(f"\n=== release-kg: {TAG_PREFIX}{ctx.version} done ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
