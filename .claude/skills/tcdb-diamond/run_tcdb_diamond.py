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
