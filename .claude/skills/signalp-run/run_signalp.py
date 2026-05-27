#!/usr/bin/env python3
"""Run SignalP 6.0 (via Docker) on genome protein FASTA files for all configured strains.

Usage:
    uv run python .claude/skills/signalp-run/run_signalp.py
    uv run python .claude/skills/signalp-run/run_signalp.py --strain HOT1A3
    uv run python .claude/skills/signalp-run/run_signalp.py --force

    # Normalize already-run raw output into the standard Phase-1 calls.json
    # (no Docker; reads each strain's signalp/prediction_results.txt):
    uv run python .claude/skills/signalp-run/run_signalp.py --normalize
    uv run python .claude/skills/signalp-run/run_signalp.py --normalize --strain MED4
"""

import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
GENOMES_CSV = REPO_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"

DOCKER_IMAGE = "signalp6"

# Make the project package importable when invoked as a script under `uv run`.
sys.path.insert(0, str(REPO_ROOT))
from multiomics_kg.utils.signalp import (  # noqa: E402
    OTHER_SENTINEL,
    is_kept,
    parse_prediction_results,
)


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


def check_docker_image() -> bool:
    """Check that the signalp6 Docker image exists."""
    result = subprocess.run(
        ["docker", "image", "inspect", DOCKER_IMAGE],
        capture_output=True,
    )
    return result.returncode == 0


def run_signalp(strain: str, genome_dir: Path) -> bool:
    """Run SignalP 6.0 in Docker for a single strain. Returns True on success."""
    signalp_dir = genome_dir / "signalp"
    signalp_dir.mkdir(parents=True, exist_ok=True)

    # Mount the genome directory into the container at /data
    abs_genome_dir = genome_dir.resolve()

    cmd = [
        "docker", "run", "--rm",
        "--user", f"{__import__('os').getuid()}:{__import__('os').getgid()}",
        "-v", f"{abs_genome_dir}:/data",
        DOCKER_IMAGE,
        "--fastafile", "/data/protein.faa",
        "--organism", "other",
        "--output_dir", "/data/signalp",
        # 'none' = summary files only (prediction_results.txt, output.gff3,
        # output.json, processed_entries.fasta, region_output.gff3). Avoids
        # writing one plot.txt per protein (~5K files per strain on the big
        # proteomes — pure noise for our use case).
        "--format", "none",
        "--mode", "fast",
        # 10 PyTorch threads on this 12-core host (leaves 2 cores free for
        # OS / docker overhead). Default is 8.
        "--torch_num_threads", "10",
    ]

    print(f"\n{'='*60}")
    print(f"Running SignalP 6.0 for {strain}")
    print(f"  input : {genome_dir / 'protein.faa'}")
    print(f"  output: {signalp_dir}")
    print(f"{'='*60}")

    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    return result.returncode == 0


def count_predictions(results_file: Path) -> dict[str, int]:
    """Count predictions by type from prediction_results.txt."""
    counts: dict[str, int] = {}
    with open(results_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                pred_type = parts[1]
                counts[pred_type] = counts.get(pred_type, 0) + 1
    return counts


def normalize_strain(strain: str, genome_dir: Path, force: bool) -> tuple[str, str]:
    """Normalize a strain's raw prediction_results.txt → standard calls.json.

    No Docker — reads the already-present raw SignalP output and writes
    ``<strain>.signalp.calls.json`` (WP_-keyed) + ``<strain>.signalp.skill_summary.json``.
    This brings ``signalp-run`` (which predates the add-a-tool calls.json
    convention) into line so Phase-2 KG integration can consume it like every
    other per-strain tool. Returns (status, message).
    """
    signalp_dir = genome_dir / "signalp"
    results_file = signalp_dir / "prediction_results.txt"
    calls_path = signalp_dir / f"{strain}.signalp.calls.json"
    summary_path = signalp_dir / f"{strain}.signalp.skill_summary.json"

    if not results_file.exists():
        return "MISSING_INPUT", "prediction_results.txt not found (run SignalP first)"
    if calls_path.exists() and not force:
        return "SKIP_EXISTS", f"{calls_path.name} already present"

    recs = parse_prediction_results(results_file.read_text())
    calls_path.write_text(json.dumps(recs, indent=2, sort_keys=True) + "\n")

    dist = dict(sorted(collections.Counter(r["signalp_type"] for r in recs.values()).items()))
    total = len(recs)
    kept = sum(1 for r in recs.values() if is_kept(r["signalp_type"]))
    summary = {
        "strain": strain,
        "tool_version": "SignalP-6.0",
        "input_proteins": total,
        "calls_made": total,                 # every record carries a call (OTHER incl.)
        "signal_peptide_count": kept,         # kept (non-OTHER) types
        "parse_failures": 0,
        "distribution": dist,
        "sentinel_rate": round(dist.get(OTHER_SENTINEL, 0) / total, 4) if total else 0.0,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return "OK", f"{total} proteins, {kept} signal peptides"


def main():
    parser = argparse.ArgumentParser(
        description="Run SignalP 6.0 (Docker) on all genome strains.",
    )
    parser.add_argument("--strain", help="Run only this strain (e.g. HOT1A3)")
    parser.add_argument("--force", action="store_true", help="Re-run even if output exists")
    parser.add_argument(
        "--normalize", action="store_true",
        help="Normalize existing raw output into <strain>.signalp.calls.json (no Docker).",
    )
    args = parser.parse_args()

    # ── Normalize mode: read raw output → calls.json (no Docker) ──────────────
    if args.normalize:
        genomes = load_genomes(args.strain)
        if not genomes:
            print(f"No genomes found{f' for strain {args.strain}' if args.strain else ''}.")
            sys.exit(1)
        results: list[tuple[str, str, str]] = []
        for g in genomes:
            strain = g["strain_name"]
            genome_dir = REPO_ROOT / g["data_dir"].rstrip("/")
            status, msg = normalize_strain(strain, genome_dir, args.force)
            results.append((strain, status, msg))
        print(f"\n{'='*70}")
        print(f"{'Strain':<12} {'Status':<15} {'Info'}")
        print(f"{'-'*12} {'-'*15} {'-'*40}")
        for strain, status, msg in results:
            print(f"{strain:<12} {status:<15} {msg}")
        print(f"{'='*70}")
        missing = [s for s, st, _ in results if st == "MISSING_INPUT"]
        if missing:
            print(f"\nMISSING_INPUT strains (run SignalP first): {', '.join(missing)}",
                  file=sys.stderr)
        return

    # Check Docker image exists
    if not check_docker_image():
        print(
            f"ERROR: Docker image '{DOCKER_IMAGE}' not found.\n"
            f"Build it first: bash .claude/skills/signalp-run/build.sh\n"
            f"See .claude/skills/signalp-run/INSTALL.md for instructions.",
            file=sys.stderr,
        )
        sys.exit(1)

    genomes = load_genomes(args.strain)
    if not genomes:
        print(f"No genomes found{f' for strain {args.strain}' if args.strain else ''}.")
        sys.exit(1)

    results: list[tuple[str, str, str]] = []  # (strain, status, message)

    for g in genomes:
        strain = g["strain_name"]
        genome_dir = REPO_ROOT / g["data_dir"].rstrip("/")
        faa = genome_dir / "protein.faa"
        results_file = genome_dir / "signalp" / "prediction_results.txt"

        if not faa.exists():
            results.append((strain, "SKIP_NO_FAA", "protein.faa not found"))
            continue

        if results_file.exists() and not args.force:
            counts = count_predictions(results_file)
            total = sum(counts.values())
            secreted = sum(v for k, v in counts.items() if k != "OTHER")
            results.append((
                strain,
                "SKIP_EXISTS",
                f"{total} proteins, {secreted} with signal peptides",
            ))
            continue

        ok = run_signalp(strain, genome_dir)
        if ok and results_file.exists():
            counts = count_predictions(results_file)
            total = sum(counts.values())
            secreted = sum(v for k, v in counts.items() if k != "OTHER")
            breakdown = ", ".join(f"{k}:{v}" for k, v in sorted(counts.items()) if k != "OTHER")
            results.append((
                strain,
                "OK",
                f"{total} proteins, {secreted} with signal peptides ({breakdown})",
            ))
        else:
            results.append((strain, "FAILED", "SignalP returned non-zero exit code"))

    # Summary table
    print(f"\n{'='*70}")
    print(f"{'Strain':<12} {'Status':<15} {'Info'}")
    print(f"{'-'*12} {'-'*15} {'-'*40}")
    for strain, status, msg in results:
        print(f"{strain:<12} {status:<15} {msg}")
    print(f"{'='*70}")

    failed = [s for s, st, _ in results if st == "FAILED"]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
