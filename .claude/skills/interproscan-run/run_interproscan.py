#!/usr/bin/env python3
"""Run InterProScan 5 (via the interpro/interproscan Docker image) on each
strain's protein.faa, emitting per-protein domain calls (Phase 1).

Usage:
    # One-time install
    uv run python .claude/skills/interproscan-run/run_interproscan.py --prepare-image
    uv run python .claude/skills/interproscan-run/run_interproscan.py --refresh-data

    # Run strains
    uv run python .claude/skills/interproscan-run/run_interproscan.py                  # all
    uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4
    uv run python .claude/skills/interproscan-run/run_interproscan.py --force
    uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4 --limit 100

See .claude/skills/interproscan-run/SKILL.md.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from multiomics_kg.download.utils.cli import load_genome_rows  # noqa: E402
from multiomics_kg.utils import tool_calls_io as tcio  # noqa: E402
from multiomics_kg.utils.interproscan import (  # noqa: E402
    parse_interproscan_json,
    summarize,
)

TOOL = "interproscan"
IPS_VERSION = "5.78-109.0"
DOCKER_IMAGE = f"interpro/interproscan:{IPS_VERSION}"
DATA_TARBALL = f"interproscan-data-{IPS_VERSION}.tar.gz"
FTP_BASE = f"https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{IPS_VERSION}/alt"


# ---------------------------------------------------------------------------
# Install / setup
# ---------------------------------------------------------------------------

def default_data_root() -> Path:
    """`$INTERPROSCAN_DATA_DIR` or `~/tools/InterProScan` (outside the repo)."""
    env = os.environ.get("INTERPROSCAN_DATA_DIR")
    return Path(env).expanduser() if env else Path.home() / "tools" / "InterProScan"


def resolve_data_dir(data_root: Path) -> Path | None:
    """Find the extracted `data/` dir; return None if not present."""
    primary = data_root / f"interproscan-{IPS_VERSION}" / "data"
    if primary.is_dir():
        return primary
    fallback = data_root / "data"
    if fallback.is_dir():
        return fallback
    return None


def image_digest() -> str | None:
    r = subprocess.run(
        ["docker", "image", "inspect", "--format", "{{index .RepoDigests 0}}", DOCKER_IMAGE],
        capture_output=True, text=True,
    )
    out = r.stdout.strip()
    return out or None


def check_image() -> bool:
    return subprocess.run(
        ["docker", "image", "inspect", DOCKER_IMAGE], capture_output=True
    ).returncode == 0


def prepare_image() -> int:
    print(f"Pulling {DOCKER_IMAGE} ...")
    rc = subprocess.run(["docker", "pull", DOCKER_IMAGE]).returncode
    if rc != 0:
        print("ERROR: docker pull failed.", file=sys.stderr)
        return rc
    print(f"\nImage ready. RepoDigest: {image_digest()}")
    return 0


def _download(url: str, dest: Path) -> None:
    print(f"  downloading {url}")
    with urllib.request.urlopen(url) as resp, dest.open("wb") as f:  # noqa: S310
        while chunk := resp.read(8 << 20):
            f.write(chunk)


def refresh_data(data_root: Path, force: bool) -> int:
    """Download + md5-verify + extract the InterProScan data tarball."""
    if resolve_data_dir(data_root) is not None and not force:
        print(f"Data dir already present under {data_root} (use --force to re-download).")
        return 0
    data_root.mkdir(parents=True, exist_ok=True)
    tarball = data_root / DATA_TARBALL
    md5_file = data_root / f"{DATA_TARBALL}.md5"

    _download(f"{FTP_BASE}/{DATA_TARBALL}", tarball)
    _download(f"{FTP_BASE}/{DATA_TARBALL}.md5", md5_file)

    expected = md5_file.read_text().split()[0].strip()
    print("  verifying md5 ...")
    h = hashlib.md5()  # noqa: S324 — checksum matches upstream .md5, not security
    with tarball.open("rb") as f:
        while chunk := f.read(8 << 20):
            h.update(chunk)
    if h.hexdigest() != expected:
        print(f"ERROR: md5 mismatch ({h.hexdigest()} != {expected}).", file=sys.stderr)
        return 1

    print(f"  extracting {tarball.name} ...")
    with tarfile.open(tarball) as tf:
        tf.extractall(data_root)  # noqa: S202 — trusted EBI tarball
    tarball.unlink()
    dd = resolve_data_dir(data_root)
    if dd is None:
        print("ERROR: extraction did not produce a data/ dir.", file=sys.stderr)
        return 1
    print(f"\nData ready at {dd}")
    print("Member-DB subdirs:", ", ".join(sorted(p.name for p in dd.iterdir() if p.is_dir())))
    return 0


# ---------------------------------------------------------------------------
# Per-strain run
# ---------------------------------------------------------------------------

def _write_limited_fasta(src: Path, dest: Path, limit: int) -> int:
    """Write the first *limit* sequences of *src* to *dest*. Returns count."""
    count = 0
    with src.open() as fin, dest.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                if count >= limit:
                    break
                count += 1
            fout.write(line)
    return count


def _count_seqs(faa: Path) -> int:
    with faa.open() as f:
        return sum(1 for line in f if line.startswith(">"))


def run_strain(strain: str, genome_dir: Path, data_dir: Path, args) -> tuple[str, str]:
    """Run InterProScan on one strain. Returns (status, message)."""
    faa = genome_dir / "protein.faa"
    if not faa.exists():
        return "MISSING_INPUT", "protein.faa not found"

    out_dir = genome_dir / TOOL
    calls_p = tcio.calls_path(genome_dir, TOOL, strain, limited=args.limit)
    if calls_p.exists() and not args.force:
        return "SKIPPED", f"{calls_p.name} exists"

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "temp").mkdir(exist_ok=True)

    # Choose input FASTA (full or truncated for --limit).
    infix = f".limited_{args.limit}" if args.limit is not None else ""
    if args.limit is not None:
        in_faa = out_dir / f"{strain}{infix}.faa"
        input_proteins = _write_limited_fasta(faa, in_faa, args.limit)
        container_input = f"/input/{TOOL}/{in_faa.name}"
    else:
        input_proteins = _count_seqs(faa)
        container_input = "/input/protein.faa"

    out_base = f"{strain}.{TOOL}{infix}.raw"
    raw_json = out_dir / f"{out_base}.json"

    # InterProScan REFUSES to overwrite an existing output file — it silently
    # writes `<base>_1.json`, `<base>_2.json`, … instead. Without this cleanup a
    # `--force` re-run would re-parse the *stale* raw JSON from the previous run
    # and report success, so e.g. newly-enabled --goterms output would never
    # appear. Clear the canonical file and any suffixed leftovers first.
    for stale in [raw_json, *out_dir.glob(f"{out_base}_[0-9]*.json")]:
        stale.unlink(missing_ok=True)

    cmd = [
        "docker", "run", "--rm",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{genome_dir.resolve()}:/input",
        "-v", f"{data_dir.resolve()}:/opt/interproscan/data:ro",
        DOCKER_IMAGE,
        "--input", container_input,
        "--output-file-base", f"/input/{TOOL}/{out_base}",
        "--formats", "JSON",
        "--disable-precalc",
        "--tempdir", f"/input/{TOOL}/temp",
        "--cpu", str(args.threads),
    ]
    if args.applications:
        cmd += ["--applications", args.applications]
    # GO + pathway xrefs on integrated InterPro entries. Both are local lookups
    # against the bundled entry data (no network, compatible with
    # --disable-precalc) and each IMPLIES -iprlookup. Without them IPS emits
    # `goXRefs: []` / `pathwayXRefs: []` on every entry.
    if not args.no_xrefs:
        cmd += ["--goterms", "--pathways"]

    log_dir = REPO_ROOT / "logs" / TOOL
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{strain}{infix}.log"

    print(f"\n{'='*60}\nInterProScan {IPS_VERSION} — {strain} "
          f"({input_proteins} proteins)\n  log: {log_file}\n{'='*60}")
    t0 = time.time()
    with log_file.open("w") as lf:
        rc = subprocess.run(cmd, cwd=str(REPO_ROOT), stdout=lf, stderr=subprocess.STDOUT).returncode
    wallclock = time.time() - t0

    if rc != 0:
        return "FAILED", f"docker exit {rc}; see {log_file}"
    if not raw_json.exists():
        return "FAILED", f"no raw JSON at {raw_json.name}; see {log_file}"

    import json
    with raw_json.open() as f:
        data = json.load(f)
    entry_xrefs: dict[str, dict[str, list[str]]] = {}
    calls = parse_interproscan_json(data, entry_xrefs)
    # Per-entry GO/pathway detail, normalized out of the match records.
    xrefs_p = calls_p.with_name(calls_p.name.replace(".calls.json", ".entry_xrefs.json"))
    with xrefs_p.open("w") as f:
        json.dump(dict(sorted(entry_xrefs.items())), f, indent=1, sort_keys=True)
        f.write("\n")

    summary = summarize(
        calls, strain=strain, input_proteins=input_proteins,
        tool_version=IPS_VERSION, applications=(args.applications or "ALL_DEFAULT"),
        image_digest=image_digest(), wallclock_s=wallclock,
        xrefs_requested=not args.no_xrefs,
    )
    tcio.save_calls(genome_dir, TOOL, strain, calls, limited=args.limit)
    tcio.save_skill_summary(genome_dir, TOOL, strain, summary, limited=args.limit)

    return "OK", (f"{summary['calls_made']} calls, {summary['total_matches']} matches, "
                  f"{summary['proteins_no_match']} no-match, {wallclock:.0f}s")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description="Run InterProScan 5 per strain (Phase 1).")
    p.add_argument("--strains", nargs="+", metavar="STRAIN", help="Only these strains")
    p.add_argument("--force", action="store_true", help="Re-run even if calls.json exists")
    p.add_argument("--limit", type=int, metavar="N", help="Only first N proteins (smoke test)")
    p.add_argument("--threads", type=int, default=10, help="CPUs for IPS (--cpu); default 10")
    p.add_argument("--applications", metavar="APP,APP",
                   help="Restrict to these member DBs (default: all activated apps)")
    p.add_argument("--no-xrefs", action="store_true",
                   help="Omit --goterms/--pathways (GO + pathway xrefs are on by default)")
    p.add_argument("--data-dir", metavar="PATH",
                   help="InterProScan install root (default $INTERPROSCAN_DATA_DIR or ~/tools/InterProScan)")
    p.add_argument("--prepare-image", action="store_true", help="docker pull the image, then exit")
    p.add_argument("--refresh-data", action="store_true", help="download + extract the data tarball, then exit")
    args = p.parse_args()

    data_root = Path(args.data_dir).expanduser() if args.data_dir else default_data_root()

    if args.prepare_image:
        sys.exit(prepare_image())
    if args.refresh_data:
        sys.exit(refresh_data(data_root, args.force))

    # Normal run — refuse to start if the install is incomplete.
    if not check_image():
        print(f"ERROR: Docker image '{DOCKER_IMAGE}' not found.\n"
              f"Run: uv run python {Path(__file__).relative_to(REPO_ROOT)} --prepare-image",
              file=sys.stderr)
        sys.exit(1)
    data_dir = resolve_data_dir(data_root)
    if data_dir is None:
        print(f"ERROR: InterProScan data dir not found under {data_root}.\n"
              f"Run: uv run python {Path(__file__).relative_to(REPO_ROOT)} --refresh-data",
              file=sys.stderr)
        sys.exit(1)

    rows = load_genome_rows(strains=args.strains)
    results: list[tuple[str, str, str]] = []
    for row in rows:
        strain = row["strain_name"]
        genome_dir = REPO_ROOT / row["data_dir"].rstrip("/")
        status, msg = run_strain(strain, genome_dir, data_dir, args)
        results.append((strain, status, msg))

    print(f"\n{'='*74}\n{'Strain':<12} {'Status':<14} {'Info'}\n{'-'*12} {'-'*14} {'-'*44}")
    for strain, status, msg in results:
        print(f"{strain:<12} {status:<14} {msg}")
    print(f"{'='*74}")

    failed = [s for s, st, _ in results if st == "FAILED"]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
