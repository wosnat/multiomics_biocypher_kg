# .claude/skills/tcdb-diamond/run_tcdb_diamond.py
"""Run diamond blastp vs. TCDB FASTA per strain. See spec
docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Ensure project root is on sys.path so `multiomics_kg` is importable when
# this script is run directly (uv run python .claude/skills/tcdb-diamond/...).
# Python adds the *script's* directory to sys.path, not the project root.
_REPO_ROOT_EARLY = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT_EARLY) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT_EARLY))

import dotenv

from multiomics_kg.utils.tcdb_diamond import build_strain_calls

REPO_ROOT = Path(__file__).resolve().parents[3]
GENOMES_CSV = REPO_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
LOGS_DIR = REPO_ROOT / "logs"
DEFAULT_TCDB_DATA_DIR = Path.home() / "tools" / "TCDB"

dotenv.load_dotenv(REPO_ROOT / ".env")


def resolve_tcdb_data_dir() -> Path:
    """Resolve the TCDB data dir (env override or default)."""
    env = os.environ.get("TCDB_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_TCDB_DATA_DIR


def build_tcdb_db(data_dir: Path, force: bool) -> Path:
    """Ensure ~/tools/TCDB/DB/tcdb.dmnd exists. Returns the path.

    Calls Saier Lab's extractTCDB.pl to download the FASTA and build the
    diamond DB. Skipped when the DB already exists and `force` is False.
    Requires:
      - extractTCDB.pl present at <data_dir>/TCDBtools/bin/
      - perl, wget, diamond on PATH
    """
    db_dir = data_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)
    dmnd = db_dir / "tcdb.dmnd"

    if dmnd.exists() and not force:
        return dmnd

    extract = data_dir / "TCDBtools" / "bin" / "extractTCDB.pl"
    if not extract.exists():
        print(
            f"ERROR: {extract} not found. Run:\n"
            f"  git clone https://github.com/SaierLaboratory/TCDBtools.git "
            f"{data_dir}/TCDBtools",
            file=sys.stderr,
        )
        sys.exit(1)

    for tool in ("perl", "wget", "diamond"):
        if not shutil.which(tool):
            print(f"ERROR: {tool} not on PATH (required by extractTCDB.pl)", file=sys.stderr)
            sys.exit(1)

    print(f"Building TCDB diamond DB in {db_dir} (this takes ~30s)...")
    for fmt in ("fasta", "diamond"):
        cmd = ["perl", str(extract), "-i", "tcdb", "-o", str(db_dir), "-f", fmt]
        result = subprocess.run(cmd, cwd=str(REPO_ROOT))
        if result.returncode != 0:
            print(f"ERROR: extractTCDB.pl -f {fmt} returned {result.returncode}", file=sys.stderr)
            sys.exit(1)

    if not dmnd.exists():
        print(f"ERROR: extractTCDB.pl completed but {dmnd} was not created", file=sys.stderr)
        sys.exit(1)

    return dmnd


def run_diamond(faa: Path, dmnd: Path, out_tsv: Path, threads: int, log_path: Path) -> bool:
    """Run diamond blastp for one strain. Returns True on success.

    Floor-only filtering at the diamond step: --evalue 0.001 only.
    Identity / coverage tiering happens in build_strain_calls (Python).

    Diamond's stdout + stderr are captured to `log_path` (overwritten on
    each run). The terminal only sees one progress line per strain.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "diamond", "blastp",
        "-q", str(faa),
        "-d", str(dmnd),
        "-o", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident",
        "qcovhsp", "scovhsp", "length", "evalue", "bitscore",
        "--evalue", "0.001",
        "--max-target-seqs", "5",
        "--more-sensitive",
        "--threads", str(threads),
    ]
    print(f"\n>>> diamond blastp {faa.name} → {out_tsv} (log: {log_path.relative_to(REPO_ROOT)})")
    with open(log_path, "w") as logf:
        logf.write(f"$ {' '.join(cmd)}\n\n")
        logf.flush()
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
    return result.returncode == 0


def truncate_faa(faa: Path, n_proteins: int, dest: Path) -> Path:
    """Copy the first N sequences of a FASTA to `dest`. Returns dest.

    Used by the --limit flag for fast end-to-end smoke tests against a small
    subset of one strain's proteome (~10-30s instead of ~5min).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(faa) as src, open(dest, "w") as out:
        seen = 0
        for line in src:
            if line.startswith(">"):
                seen += 1
                if seen > n_proteins:
                    break
            out.write(line)
    return dest


def load_genomes(strain_filter: str | None) -> list[dict]:
    """Parse cyanobacteria_genomes.csv, return rows for genome strains only.

    Skips reference_proteome_match and treatment organisms — only
    organism_type='genome_strain' rows have a `protein.faa`.
    """
    genomes: list[dict] = []
    with open(GENOMES_CSV) as f:
        reader = csv.DictReader(row for row in f if not row.strip().startswith("#"))
        for row in reader:
            if strain_filter and row["strain_name"] != strain_filter:
                continue
            if (row.get("organism_type") or "").strip() != "genome_strain":
                continue
            genomes.append(row)
    return genomes


def process_strain(
    strain: str,
    data_dir_genome: Path,
    dmnd: Path,
    threads: int,
    force: bool,
    limit: int | None = None,
) -> tuple[str, str, dict | None]:
    """Run the per-strain pipeline. Returns (strain, status, summary_or_None).

    When `limit` is set, runs diamond against the first `limit` proteins from
    protein.faa (truncated copy in /tmp), and writes outputs to
    `<strain>.tcdb.limited_<N>.{tsv,calls.json,skill_summary.json}` to avoid
    overwriting full-run artifacts. Skip-if-exists logic also keys on the
    limit-suffixed name so re-running with a different limit is a fresh run.

    Status values:
      OK / SKIP_NO_FAA / SKIP_EXISTS / FAILED_DIAMOND / FAILED_NO_EGGNOG
    """
    faa = data_dir_genome / "protein.faa"
    if not faa.exists():
        return strain, "SKIP_NO_FAA", None

    out_dir = data_dir_genome / "tcdb"
    suffix = f".limited_{limit}" if limit else ""
    out_tsv = out_dir / f"{strain}.tcdb{suffix}.tsv"
    out_calls = out_dir / f"{strain}.tcdb{suffix}.calls.json"
    out_summary = out_dir / f"{strain}.tcdb{suffix}.skill_summary.json"

    if out_calls.exists() and not force:
        with open(out_summary) as f:
            return strain, "SKIP_EXISTS", json.load(f)

    if limit:
        truncated = Path("/tmp") / f"tcdb_diamond_{strain}_first{limit}.faa"
        faa = truncate_faa(faa, limit, truncated)

    log_path = LOGS_DIR / f"tcdb_diamond_{strain}{suffix}.log"
    if not run_diamond(faa, dmnd, out_tsv, threads, log_path):
        print(f"  see log: {log_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        return strain, "FAILED_DIAMOND", None

    eggnog_path = data_dir_genome / "eggnog" / f"{strain}.emapper.annotations"
    if not eggnog_path.exists():
        # eggNOG missing → still proceed but agreement column will all be "extends"
        print(f"  WARN: {eggnog_path} missing — egn_agreement will all be 'extends'")

    calls, summary = build_strain_calls(out_tsv, eggnog_path)

    with open(out_calls, "w") as f:
        json.dump(calls, f, indent=2, sort_keys=True)
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    return strain, "OK", summary


def main():
    parser = argparse.ArgumentParser(
        description="Run diamond vs. TCDB FASTA per strain (Phase 1)."
    )
    parser.add_argument("--strain", help="Run only this strain (e.g. MED4)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if calls.json exists")
    parser.add_argument("--refresh-tcdb", action="store_true",
                        help="Re-download TCDB FASTA + diamond DB even if cached")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Diamond threads (default: os.cpu_count() or 4)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: run on first N proteins of each strain "
                             "only. Outputs go to <strain>.tcdb.limited_<N>.* "
                             "alongside (not replacing) full-run artifacts.")
    args = parser.parse_args()

    data_dir = resolve_tcdb_data_dir()
    print(f"TCDB data dir: {data_dir}")
    dmnd = build_tcdb_db(data_dir, force=args.refresh_tcdb)
    print(f"diamond DB: {dmnd}")

    genomes = load_genomes(args.strain)
    if not genomes:
        print(f"No genome_strain genomes found"
              f"{f' for strain {args.strain}' if args.strain else ''}.")
        sys.exit(1)

    results: list[tuple[str, str, dict | None]] = []
    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        results.append(
            process_strain(strain, data_dir_genome, dmnd, args.threads, args.force, args.limit)
        )

    # Status table
    print(f"\n{'='*92}")
    cols = ("Strain", "Status", "Hits", "T1", "T2", "T3",
            "confirms", "refines", "extends", "conflicts")
    print("{:<12} {:<14} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10}".format(*cols))
    print("-" * 92)
    for strain, status, summary in results:
        if summary is None:
            print(f"{strain:<12} {status:<14}")
            continue
        td = summary["tier_distribution"]
        ad = summary["agreement_distribution"]
        print(
            "{:<12} {:<14} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10}".format(
                strain, status, summary["proteins_with_call"],
                td.get("1", 0), td.get("2", 0), td.get("3", 0),
                ad.get("confirms", 0), ad.get("refines", 0),
                ad.get("extends", 0), ad.get("conflicts", 0),
            )
        )
    print("=" * 92)

    failed = [s for s, st, _ in results if st.startswith("FAILED")]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
