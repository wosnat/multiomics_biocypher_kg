#!/usr/bin/env python3
"""Stamp a release's CHANGELOG notes onto a live Schema_info without a rebuild.

Sets Schema_info.release_highlights / release_breaking (from the TAGGED
CHANGELOG's `### Highlights` / `### Breaking` of <version>) and
release_notes_url (the GitHub Release URL), on the deploy container serving
that version. Refuses if the live Schema_info.version is not <version>.

Written for kg-0.1.0-alpha.7 / alpha.8, which shipped with both properties null
because no compose file forwarded KG_RELEASE_HIGHLIGHTS / KG_RELEASE_BREAKING
into post-process (fixed for later releases).

    uv run python .claude/skills/release-kg/restamp_release_notes.py 0.1.0-alpha.8
        [--container alpha-deploy] [--env-file .env.alpha] [--dry-run]

Credentials: NEO4J_AUTH (user/password) from --env-file.
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]


def _load_release_kg():
    spec = importlib.util.spec_from_file_location("release_kg", HERE / "release_kg.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["release_kg"] = mod
    spec.loader.exec_module(mod)
    return mod


def cypher_str(value: str) -> str:
    """Body of a Cypher single-quoted literal, same rules as post-import.sh:
    backslashes, then single quotes, then newlines as the escape \\n
    (cypher-shell -P rejects a raw newline)."""
    return value.replace("\\", "\\\\").replace("'", "\\'").replace("\n", "\\n")


def read_env(path: Path) -> dict:
    env = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, _, v = line.partition("=")
            env[k.strip()] = v.strip().strip('"').strip("'")
    return env


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("version")
    ap.add_argument("--container", default="alpha-deploy")
    ap.add_argument("--env-file", default=str(REPO_ROOT / ".env.alpha"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rkg = _load_release_kg()
    tag = f"{rkg.TAG_PREFIX}{args.version}"
    changelog = subprocess.run(["git", "-C", str(REPO_ROOT), "show", f"{tag}:CHANGELOG.md"],
                               capture_output=True, text=True, check=True).stdout
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as fh:
        fh.write(changelog)
    highlights = rkg.extract_preflight_subsection(Path(fh.name), args.version, "Highlights") or ""
    breaking = rkg.extract_preflight_subsection(Path(fh.name), args.version, "Breaking") or ""
    Path(fh.name).unlink()
    url = rkg.release_notes_url(args.version)
    print(f"{tag}: highlights {len(highlights)} chars, breaking {len(breaking)} chars, url {url or '<none>'}")

    user, _, password = read_env(Path(args.env_file))["NEO4J_AUTH"].partition("/")
    cypher = (
        "MATCH (s:Schema_info {id:'schema_info'}) WHERE s.version = $version "
        "SET s.release_highlights = CASE WHEN $h = '' THEN null ELSE $h END, "
        "    s.release_breaking   = CASE WHEN $b = '' THEN null ELSE $b END, "
        "    s.release_notes_url  = $u "
        "RETURN s.version AS version, size(coalesce(s.release_highlights,'')) AS highlights_chars, "
        "size(coalesce(s.release_breaking,'')) AS breaking_chars, s.release_notes_url AS url"
    )
    cmd = ["docker", "exec", "-e", "LANG=C.UTF-8", "-e", "LC_ALL=C.UTF-8", args.container,
           "cypher-shell", "-a", "bolt://localhost:7687", "-u", user, "-p", password,
           "--format", "plain",
           "-P", f"version => '{cypher_str(args.version)}'",
           "-P", f"h => '{cypher_str(highlights)}'",
           "-P", f"b => '{cypher_str(breaking)}'",
           "-P", f"u => '{cypher_str(url)}'",
           cypher]
    if args.dry_run:
        print(f"[dry-run] would SET release_* on {args.container} (Schema_info.version = {args.version})")
        return
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"cypher-shell failed: {out.stderr.strip()}")
    rows = out.stdout.strip().splitlines()
    if len(rows) < 2:
        sys.exit(f"no Schema_info with version {args.version!r} on {args.container}; nothing stamped")
    print(rows[-1])
    rkg.assert_subsections_stamped(rkg._parse_plain_row(out.stdout), highlights, breaking,
                                   args.container)
    print("stamped ✓")


if __name__ == "__main__":
    main()
