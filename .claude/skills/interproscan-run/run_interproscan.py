#!/usr/bin/env python3
"""Run InterProScan 6 (the ebi-pf-team/interproscan6 Nextflow pipeline, Docker
profile) on each strain's protein.faa, emitting per-protein domain calls
(Phase 1). Pinned to an InterProScan release and an InterPro data release.

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
import gzip
import shutil
import os
import subprocess
import sys
import tempfile
import time
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
PIPELINE = "ebi-pf-team/interproscan6"
IPS_VERSION = "6.0.3.0"      # pipeline release (git tag)
INTERPRO_RELEASE = "111.0"   # pins every member-DB data version
TOOL_VERSION = f"{IPS_VERSION}-{INTERPRO_RELEASE}"


# ---------------------------------------------------------------------------
# Install / setup
# ---------------------------------------------------------------------------

def default_data_root() -> Path:
    """`$INTERPROSCAN_DATA_DIR` or `~/tools/InterProScan` (outside the repo)."""
    env = os.environ.get("INTERPROSCAN_DATA_DIR")
    return Path(env).expanduser() if env else Path.home() / "tools" / "InterProScan"


def v6_data_dir(data_root: Path) -> Path:
    """IPS6 member-DB data, one subdir per database and version
    (`pfam/39.0`, `interpro/111.0`, ...) — downloaded by the pipeline itself."""
    return data_root / "v6" / "data"


def nextflow_bin() -> str | None:
    """`$NEXTFLOW`, else `~/tools/nextflow/nextflow`, else `nextflow` on PATH."""
    env = os.environ.get("NEXTFLOW")
    if env:
        return str(Path(env).expanduser())
    local = Path.home() / "tools" / "nextflow" / "nextflow"
    if local.exists():
        return str(local)
    return shutil.which("nextflow")


def data_ready(data_dir: Path) -> bool:
    return (data_dir / "interpro" / INTERPRO_RELEASE).is_dir()


def nextflow_cmd(nf: str, faa: Path, data_dir: Path, outdir: Path, outprefix: str,
                 work_dir: Path, *, cpus: int, applications: str | None,
                 xrefs: bool) -> list[str]:
    cmd = [
        nf, "run", PIPELINE, "-r", IPS_VERSION, "-profile", "docker",
        "-work-dir", str(work_dir),
        "--input", str(faa),
        "--datadir", str(data_dir),
        "--interpro", INTERPRO_RELEASE,
        # Compute every match locally instead of fetching precalculated ones
        # from the EBI matches API (IPS5's --disable-precalc): reproducible,
        # no network dependency during the scan.
        "--noMatchesApi",
        "--formats", "json",
        "--cpus", str(cpus),
        "--outdir", str(outdir),
        "--outprefix", outprefix,
    ]
    if applications:
        cmd += ["--applications", applications]
    # GO + pathway xrefs on integrated InterPro entries, looked up locally in
    # the interpro/<release> data. Pathways are dropped by the parser.
    if xrefs:
        cmd += ["--goterms", "--pathways"]
    return cmd


def prepare_pipeline(nf: str) -> int:
    print(f"Pulling {PIPELINE} {IPS_VERSION} ...")
    return subprocess.run([nf, "pull", PIPELINE, "-r", IPS_VERSION]).returncode


def refresh_data(nf: str, data_root: Path) -> int:
    """Download the pinned InterPro release's member-DB data (~6.6 GB, ~35 GB
    extracted) by running the pipeline on one protein: IPS6 fetches any
    missing database into --datadir before scanning."""
    data_dir = v6_data_dir(data_root)
    data_dir.mkdir(parents=True, exist_ok=True)
    faa = REPO_ROOT / "cache/data/Prochlorococcus/genomes/MED4/protein.faa"
    with tempfile.TemporaryDirectory(dir=data_root) as tmp:
        tmp = Path(tmp)
        one = tmp / "one.faa"
        _write_limited_fasta(faa, one, 1)
        rc = subprocess.run(
            nextflow_cmd(nf, one, data_dir, tmp / "out", "one", tmp / "work",
                         cpus=2, applications=None, xrefs=False),
            cwd=tmp,
        ).returncode
    if rc != 0 or not data_ready(data_dir):
        print("ERROR: data download failed — re-run (completed databases are kept).",
              file=sys.stderr)
        return rc or 1
    print(f"\nData ready at {data_dir}")
    print("Databases:", ", ".join(sorted(p.name for p in data_dir.iterdir() if p.is_dir())))
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


def run_strain(strain: str, genome_dir: Path, data_dir: Path, nf: str,
               args) -> tuple[str, str]:
    """Run InterProScan on one strain. Returns (status, message)."""
    faa = genome_dir / "protein.faa"
    if not faa.exists():
        return "MISSING_INPUT", "protein.faa not found"

    out_dir = genome_dir / TOOL
    calls_p = tcio.calls_path(genome_dir, TOOL, strain, limited=args.limit)
    if calls_p.exists() and not args.force:
        return "SKIPPED", f"{calls_p.name} exists"

    out_dir.mkdir(parents=True, exist_ok=True)
    # Nextflow work dir + .nextflow* run logs. Removed after a successful run
    # (several GB of per-task intermediates); kept on failure for debugging.
    temp = out_dir / "temp"
    temp.mkdir(exist_ok=True)

    # Choose input FASTA (full or truncated for --limit).
    infix = f".limited_{args.limit}" if args.limit is not None else ""
    if args.limit is not None:
        in_faa = out_dir / f"{strain}{infix}.faa"
        input_proteins = _write_limited_fasta(faa, in_faa, args.limit)
    else:
        in_faa = faa
        input_proteins = _count_seqs(faa)

    out_base = f"{strain}.{TOOL}{infix}.raw"
    raw_json = out_dir / f"{out_base}.json"

    # Never let a re-run parse the previous run's output: InterProScan 5
    # refused to overwrite and silently wrote `<base>_1.json`, so a `--force`
    # run re-parsed stale data and reported success. Clear the canonical file
    # and any suffixed leftovers first.
    for stale in [*_raw_candidates(out_dir, out_base),
                  *out_dir.glob(f"{out_base}_[0-9]*.json")]:
        stale.unlink(missing_ok=True)

    cmd = nextflow_cmd(nf, in_faa.resolve(), data_dir.resolve(), out_dir.resolve(),
                       out_base, (temp / "work").resolve(), cpus=args.threads,
                       applications=args.applications, xrefs=not args.no_xrefs)

    log_dir = REPO_ROOT / "logs" / TOOL
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{strain}{infix}.log"

    print(f"\n{'='*60}\nInterProScan {IPS_VERSION} / InterPro {INTERPRO_RELEASE} — {strain} "
          f"({input_proteins} proteins)\n  log: {log_file}\n{'='*60}")
    t0 = time.time()
    with log_file.open("w") as lf:
        rc = subprocess.run(cmd, cwd=str(temp), stdout=lf, stderr=subprocess.STDOUT).returncode
    wallclock = time.time() - t0

    if rc != 0:
        return "FAILED", f"nextflow exit {rc}; see {log_file}"
    if not raw_json.exists():
        return "FAILED", f"no raw JSON at {raw_json.name}; see {log_file}"

    raw = _load_raw(raw_json)
    calls = parse_interproscan_json(raw)

    summary = summarize(
        calls, strain=strain, input_proteins=input_proteins,
        tool_version=_raw_tool_version(raw),
        applications=(args.applications or "ALL_DEFAULT"),
        image_digest=f"{PIPELINE}@{IPS_VERSION}", wallclock_s=wallclock,
        xrefs_requested=not args.no_xrefs,
    )
    tcio.save_calls(genome_dir, TOOL, strain, calls, limited=args.limit)
    tcio.save_skill_summary(genome_dir, TOOL, strain, summary, limited=args.limit)
    # Kept for --normalize (re-parse without re-scanning), never committed.
    compress_raw(raw_json)
    shutil.rmtree(temp, ignore_errors=True)
    if args.limit is not None:
        in_faa.unlink(missing_ok=True)

    return "OK", (f"{summary['calls_made']} calls, {summary['total_matches']} matches, "
                  f"{summary['proteins_no_match']} no-match, {wallclock:.0f}s")


def _raw_candidates(out_dir: Path, base: str) -> list[Path]:
    """The raw IPS JSON for `base`, compressed first: `<base>.json.gz`, `<base>.json`."""
    return [out_dir / f"{base}.json.gz", out_dir / f"{base}.json"]


def _raw_tool_version(raw: dict) -> str:
    """`<interproscan>-<interpro>` as recorded in the raw JSON (IPS6 carries
    both; IPS5's version string already embeds the InterPro release)."""
    ips = raw.get("interproscan-version", "unknown")
    ipr = raw.get("interpro-version")
    return f"{ips}-{ipr}" if ipr else ips


def _load_raw(path: Path) -> dict:
    import json
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return json.load(f)


def compress_raw(raw_json: Path) -> Path:
    """gzip `<base>.json` → `<base>.json.gz` (~9x smaller; `--normalize` reads
    either) and remove the plain file. Written to a temp name first, so an
    interrupted run never leaves a truncated `.gz` in place of the original."""
    gz = raw_json.with_name(raw_json.name + ".gz")
    tmp = gz.with_name(gz.name + ".tmp")
    with raw_json.open("rb") as src, gzip.open(tmp, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    tmp.replace(gz)
    raw_json.unlink()
    return gz


def normalize_strain(strain: str, data_dir: Path) -> tuple[str, str]:
    """Re-parse a strain's cached raw.json into the new faceted calls.json —
    no Docker, no scan. Drops the stale entry_xrefs.json sidecar if present.

    Never raises — a corrupt/truncated raw.json (or any other parse/write
    failure) yields a FAILED status row instead of aborting the whole batch.
    """
    out_dir = data_dir / TOOL
    raw_json = next((p for p in _raw_candidates(out_dir, f"{strain}.{TOOL}.raw")
                     if p.exists()), None)
    if raw_json is None:
        return "NO_RAW", "raw.json(.gz) missing — re-run the scan or use committed calls.json"

    try:
        t0 = time.time()
        raw = _load_raw(raw_json)
        calls = parse_interproscan_json(raw)
        # Match scan mode's input_proteins semantics (FASTA header count) so
        # sentinel_rate has the same denominator in both modes; fall back to
        # the fanned-out call count when protein.faa is not on disk.
        faa = data_dir / "protein.faa"
        n_prot = _count_seqs(faa) if faa.exists() else len(calls)
        summary = summarize(
            calls, strain=strain, input_proteins=n_prot,
            tool_version=_raw_tool_version(raw), applications="ALL_DEFAULT",
            wallclock_s=time.time() - t0,
        )
        tcio.save_calls(data_dir, TOOL, strain, calls)
        tcio.save_skill_summary(data_dir, TOOL, strain, summary)
        # Sidecar removed only on successful normalize — a FAILED strain's
        # artifacts (including any stale entry_xrefs.json) are left untouched.
        stale = out_dir / f"{strain}.{TOOL}.entry_xrefs.json"
        if stale.exists():
            stale.unlink()
        # A plain raw.json (scan interrupted before compression) is gzipped
        # here too, so the disk-only cache ends in one form either way.
        if raw_json.suffix != ".gz":
            compress_raw(raw_json)
    except Exception as exc:  # noqa: BLE001 — one bad strain must not kill the batch
        return "FAILED", f"{type(exc).__name__}: {exc}"

    return "NORMALIZED", f"{n_prot} proteins"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description="Run InterProScan 6 per strain (Phase 1).")
    p.add_argument("--strains", nargs="+", metavar="STRAIN", help="Only these strains")
    p.add_argument("--force", action="store_true", help="Re-run even if calls.json exists")
    p.add_argument("--limit", type=int, metavar="N", help="Only first N proteins (smoke test)")
    p.add_argument("--threads", type=int, default=8, help="CPUs per IPS6 task (--cpus); default 8")
    p.add_argument("--applications", metavar="APP,APP",
                   help="Restrict to these member DBs (default: all activated apps)")
    p.add_argument("--no-xrefs", action="store_true",
                   help="Omit --goterms/--pathways (GO + pathway xrefs are on by default)")
    p.add_argument("--data-dir", metavar="PATH",
                   help="InterProScan install root (default $INTERPROSCAN_DATA_DIR or ~/tools/InterProScan)")
    p.add_argument("--prepare-image", action="store_true",
                   help="nextflow pull the pinned pipeline release, then exit")
    p.add_argument("--refresh-data", action="store_true",
                   help=f"download the InterPro {INTERPRO_RELEASE} member-DB data, then exit")
    p.add_argument("--normalize", action="store_true",
                   help="Re-parse cached raw.json(.gz) into calls.json (new faceted format); no Docker, no scan")
    args = p.parse_args()

    data_root = Path(args.data_dir).expanduser() if args.data_dir else default_data_root()

    nf = nextflow_bin()
    if (args.prepare_image or args.refresh_data or not args.normalize) and nf is None:
        print("ERROR: nextflow not found (set $NEXTFLOW or install to ~/tools/nextflow).",
              file=sys.stderr)
        sys.exit(1)
    if args.prepare_image:
        sys.exit(prepare_pipeline(nf))
    if args.refresh_data:
        sys.exit(refresh_data(nf, data_root))

    if args.normalize:
        rows = load_genome_rows(strains=args.strains)
        results: list[tuple[str, str, str]] = []
        for row in rows:
            strain = row["strain_name"]
            genome_dir = REPO_ROOT / row["data_dir"].rstrip("/")
            status, msg = normalize_strain(strain, genome_dir)
            results.append((strain, status, msg))

        print(f"\n{'='*74}\n{'Strain':<12} {'Status':<14} {'Info'}\n{'-'*12} {'-'*14} {'-'*44}")
        for strain, status, msg in results:
            print(f"{strain:<12} {status:<14} {msg}")
        print(f"{'='*74}")

        failed = [s for s, st, _ in results if st == "FAILED"]
        if failed:
            print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
            sys.exit(1)
        sys.exit(0)

    # Normal run — refuse to start if the data is missing (the pipeline would
    # otherwise start a ~6.6 GB download inside the first strain's run).
    data_dir = v6_data_dir(data_root)
    if not data_ready(data_dir):
        print(f"ERROR: InterPro {INTERPRO_RELEASE} data not found under {data_dir}.\n"
              f"Run: uv run python {Path(__file__).relative_to(REPO_ROOT)} --refresh-data",
              file=sys.stderr)
        sys.exit(1)

    rows = load_genome_rows(strains=args.strains)
    results: list[tuple[str, str, str]] = []
    for row in rows:
        strain = row["strain_name"]
        genome_dir = REPO_ROOT / row["data_dir"].rstrip("/")
        status, msg = run_strain(strain, genome_dir, data_dir, nf, args)
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
