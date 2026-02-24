#!/usr/bin/env python3
"""Run eggNOG-mapper on genome protein FASTA files for all configured strains."""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

import dotenv

REPO_ROOT = Path(__file__).resolve().parents[3]
GENOMES_CSV = REPO_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"

dotenv.load_dotenv(REPO_ROOT / ".env")


def load_genomes(strain_filter: str | None) -> list[dict]:
    """Parse cyanobacteria_genomes.csv, return list of genome dicts."""
    genomes = []
    with open(GENOMES_CSV) as f:
        reader = csv.DictReader(row for row in f if not row.strip().startswith("#"))
        for row in reader:
            if strain_filter and row["strain_name"] != strain_filter:
                continue
            genomes.append(row)
    return genomes


def run_emapper(strain: str, faa: Path, out_dir: Path, data_dir: str, cpu: int) -> bool:
    """Run emapper.py for a single strain. Returns True on success."""
    out_dir.mkdir(parents=True, exist_ok=True)

    emapper_bin = Path(sys.executable).parent / "emapper.py"
    cmd = [str(emapper_bin)] if emapper_bin.exists() else ["emapper.py"]
    cmd += [
        "--data_dir", data_dir,
        "-i", str(faa),
        "--itype", "proteins",
        "-o", strain,
        "--output_dir", str(out_dir),
        "--cpu", str(cpu),
        "--override",
    ]

    print(f"\n{'='*60}")
    print(f"Running eggNOG-mapper for {strain}")
    print(f"  input : {faa}")
    print(f"  output: {out_dir}")
    print(f"{'='*60}")

    env = os.environ.copy()
    env["EGGNOG_DATA_DIR"] = data_dir  # ensure expanded path reaches subprocess
    result = subprocess.run(cmd, cwd=str(REPO_ROOT), env=env)
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Run eggNOG-mapper on all genome strains.")
    parser.add_argument("--strain", help="Run only this strain (e.g. MED4)")
    parser.add_argument("--force", action="store_true", help="Re-run even if output exists")
    parser.add_argument("--cpu", type=int, default=4, help="CPUs per emapper run (default: 4)")
    args = parser.parse_args()

    data_dir = os.environ.get("EGGNOG_DATA_DIR", "")
    if not data_dir:
        print("ERROR: EGGNOG_DATA_DIR not set in .env", file=sys.stderr)
        sys.exit(1)
    data_dir = str(Path(data_dir).expanduser())
    if not Path(data_dir).exists():
        print(f"ERROR: EGGNOG_DATA_DIR does not exist: {data_dir}", file=sys.stderr)
        sys.exit(1)

    genomes = load_genomes(args.strain)
    if not genomes:
        print(f"No genomes found{f' for strain {args.strain}' if args.strain else ''}.")
        sys.exit(1)

    results = []  # (strain, status, message)

    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        faa = data_dir_genome / "protein.faa"
        out_dir = data_dir_genome / "eggnog"
        annotation_file = out_dir / f"{strain}.emapper.annotations"

        if not faa.exists():
            results.append((strain, "SKIP_NO_FAA", f"protein.faa not found"))
            continue

        if annotation_file.exists() and not args.force:
            n = sum(1 for line in annotation_file.read_text().splitlines()
                    if line and not line.startswith("#"))
            results.append((strain, "SKIP_EXISTS", f"{n} annotations already present"))
            continue

        ok = run_emapper(strain, faa, out_dir, data_dir, args.cpu)
        if ok and annotation_file.exists():
            n = sum(1 for line in annotation_file.read_text().splitlines()
                    if line and not line.startswith("#"))
            results.append((strain, "OK", f"{n} proteins annotated"))
        else:
            results.append((strain, "FAILED", "emapper.py returned non-zero exit code"))

    # Summary table
    print(f"\n{'='*60}")
    print(f"{'Strain':<15} {'Status':<15} {'Info'}")
    print(f"{'-'*15} {'-'*15} {'-'*30}")
    for strain, status, msg in results:
        print(f"{strain:<15} {status:<15} {msg}")
    print(f"{'='*60}")

    failed = [s for s, st, _ in results if st == "FAILED"]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
