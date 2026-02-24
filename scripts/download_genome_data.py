#!/usr/bin/env python3
"""
Download pipeline for genome annotation data.

Reads data/Prochlorococcus/genomes/cyanobacteria_genomes.csv and for each genome:

  Step 1: Download NCBI genome (GFF, protein FASTA, GBFF)
           → cache/data/<org>/genomes/<strain>/
  Step 2: Download Cyanorak GFF and GBK annotations
           (only for strains with cyanorak_organism set)
           → cache/data/<org>/genomes/<strain>/cyanorak/
  Step 3: Download UniProt data per unique (org_group, taxid)
           → cache/data/<org_group>/uniprot/<taxid>/
  Step 4: Run eggNOG-mapper on protein FASTA
           → cache/data/<org>/genomes/<strain>/eggnog/

All steps skip existing cache files by default.
Use --force to re-download/re-run for specified strains.

Usage:
  uv run python scripts/download_genome_data.py
  uv run python scripts/download_genome_data.py --steps 1 2 3
  uv run python scripts/download_genome_data.py --strains MED4 MIT9313
  uv run python scripts/download_genome_data.py --strains MED4 --force
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import subprocess
import sys
from contextlib import ExitStack
from pathlib import Path

import dotenv

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
GENOMES_CSV = PROJECT_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"

dotenv.load_dotenv(PROJECT_ROOT / ".env")

log = logging.getLogger("download_genome_data")
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


# ── helpers ─────────────────────────────────────────────────────────────────────

def _read_genomes_csv(path: Path) -> list[dict]:
    """Parse cyanobacteria_genomes.csv, skipping comment lines."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(row for row in f if not row.strip().startswith("#"))
        for row in reader:
            rows.append(row)
    return rows


def _get_org_group(data_dir: str) -> str:
    """Extract organism group name from a data_dir path.

    e.g. 'cache/data/Prochlorococcus/genomes/MED4/' -> 'Prochlorococcus'
    Expects paths with structure: cache/data/<org_group>/genomes/<strain>/
    """
    parts = Path(data_dir.rstrip("/")).parts
    # cache(0) / data(1) / <org_group>(2) / genomes(3) / <strain>(4)
    if len(parts) >= 3 and parts[0] == "cache" and parts[1] == "data":
        return parts[2]
    # fallback: first non-structural segment
    for p in parts:
        if p not in (".", "cache", "data", "genomes"):
            return p
    return "unknown"


def _header(step_num: int, desc: str) -> None:
    log.info("─" * 60)
    log.info(f"Step {step_num}: {desc}")
    log.info("─" * 60)


# ── step 1: NCBI genome ──────────────────────────────────────────────────────────

def step1_ncbi(genomes: list[dict], force: bool) -> None:
    """Download NCBI genome files (GFF, protein FASTA, GBFF) for each strain."""
    _header(1, "NCBI genome download")
    from pypath.share import curl
    from multiomics_kg.adapters.cyanorak_ncbi_adapter import CyanorakNcbi

    for g in genomes:
        strain = g["strain_name"]
        data_dir = str(PROJECT_ROOT / g["data_dir"].rstrip("/"))
        gff_path = os.path.join(data_dir, "genomic.gff")
        prot_path = os.path.join(data_dir, "protein.faa")
        gbff_path = os.path.join(data_dir, "genomic.gbff")

        all_exist = all(os.path.exists(p) for p in [gff_path, prot_path, gbff_path])
        if all_exist and not force:
            log.info(f"  [SKIP] {strain} — NCBI files already cached")
            continue

        if force:
            for p in [gff_path, prot_path, gbff_path]:
                if os.path.exists(p):
                    os.remove(p)
                    log.debug(f"  Removed {p}")

        log.info(f"  Downloading NCBI genome for {strain} ({g['ncbi_accession']})...")
        os.makedirs(data_dir, exist_ok=True)
        adapter = CyanorakNcbi(
            ncbi_accession=g["ncbi_accession"],
            data_dir=data_dir,
        )
        with ExitStack() as stack:
            if force:
                stack.enter_context(curl.cache_off())
            adapter._download_ncbi_genome()
        log.info(f"  OK: {strain}")


# ── step 2: Cyanorak ─────────────────────────────────────────────────────────────

def step2_cyanorak(genomes: list[dict], force: bool) -> None:
    """Download Cyanorak GFF and GBK annotation files."""
    _header(2, "Cyanorak GFF/GBK download")
    from pypath.share import curl
    from multiomics_kg.adapters.cyanorak_ncbi_adapter import CyanorakNcbi

    for g in genomes:
        strain = g["strain_name"]
        cyan_org = g.get("cyanorak_organism", "").strip()
        if not cyan_org:
            log.info(f"  [SKIP] {strain} — no cyanorak_organism in config")
            continue

        data_dir = str(PROJECT_ROOT / g["data_dir"].rstrip("/"))
        cyan_dir = os.path.join(data_dir, "cyanorak")
        gff_path = os.path.join(cyan_dir, f"{cyan_org}.gff")
        gbk_path = os.path.join(cyan_dir, f"{cyan_org}.gbk")

        if os.path.exists(gff_path) and os.path.exists(gbk_path) and not force:
            log.info(f"  [SKIP] {strain} — Cyanorak files already cached")
            continue

        if force:
            for p in [gff_path, gbk_path]:
                if os.path.exists(p):
                    os.remove(p)
                    log.debug(f"  Removed {p}")

        log.info(f"  Downloading Cyanorak data for {strain} ({cyan_org})...")
        os.makedirs(cyan_dir, exist_ok=True)
        adapter = CyanorakNcbi(
            cyanorak_organism=cyan_org,
            data_dir=data_dir,
        )
        with ExitStack() as stack:
            if force:
                stack.enter_context(curl.cache_off())
            adapter._download_cyanorak_gff()
            adapter._download_cyanorak_gbk()
        log.info(f"  OK: {strain}")


# ── step 3: UniProt ──────────────────────────────────────────────────────────────

def step3_uniprot(genomes: list[dict], force: bool) -> None:
    """Download UniProt data per unique (org_group, taxid) pair.

    Uses rev=False to fetch all UniProt entries (TrEMBL + SwissProt) for maximum
    annotation coverage in the gene annotation merge step.

    Output files per taxid:
      cache/data/<org_group>/uniprot/<taxid>/uniprot_raw_data.json
      cache/data/<org_group>/uniprot/<taxid>/uniprot_preprocess_data.json
    """
    _header(3, "UniProt data download (per taxid)")
    from pypath.share import curl
    from multiomics_kg.adapters.uniprot_adapter import Uniprot

    # Deduplicate: one download per (org_group, taxid)
    # Use OrderedDict to preserve first-occurrence order
    seen: dict[tuple[str, int], dict] = {}
    for g in genomes:
        taxid = int(g["ncbi_taxon_id"])
        org_group = _get_org_group(g["data_dir"])
        key = (org_group, taxid)
        if key not in seen:
            seen[key] = g

    for (org_group, taxid), g in seen.items():
        strain = g["strain_name"]
        cache_dir = PROJECT_ROOT / "cache" / "data" / org_group / "uniprot" / str(taxid)
        raw_path = cache_dir / "uniprot_raw_data.json"
        pre_path = cache_dir / "uniprot_preprocess_data.json"

        if raw_path.exists() and pre_path.exists() and not force:
            log.info(
                f"  [SKIP] UniProt taxid={taxid} ({org_group}) — cache exists"
            )
            continue

        log.info(
            f"  Downloading UniProt for taxid={taxid} ({org_group}); "
            f"representative strain: {strain}"
        )
        cache_dir.mkdir(parents=True, exist_ok=True)

        # rev=False: fetch all entries (TrEMBL + SwissProt) for maximum annotation
        # coverage. The KG build uses rev=True (SwissProt only); here we want
        # every available functional annotation for the gene annotation merge.
        adapter = Uniprot(organism=taxid, rev=False)

        with ExitStack() as stack:
            stack.enter_context(curl.cache_off())
            adapter._download_uniprot_data()

        log.info(f"  Saving raw data ({len(adapter.uniprot_ids)} proteins)...")
        with open(raw_path, "w") as f:
            json.dump(adapter.data, f, indent=4, default=str)

        adapter._preprocess_uniprot_data()
        adapter._preprocess_organisms()

        log.info(f"  Saving preprocessed data...")
        with open(pre_path, "w") as f:
            json.dump(adapter.data, f, indent=4, default=str)

        log.info(f"  OK: UniProt taxid={taxid} → {cache_dir}")


# ── step 4: eggNOG-mapper ────────────────────────────────────────────────────────

def step4_eggnog(genomes: list[dict], force: bool, cpu: int) -> None:
    """Run eggNOG-mapper on protein FASTA for each strain.

    Requires EGGNOG_DATA_DIR set in .env.
    Skips strains where <strain>.emapper.annotations already exists (unless --force).
    """
    _header(4, "eggNOG-mapper annotation")

    data_dir_env = os.environ.get("EGGNOG_DATA_DIR", "")
    if not data_dir_env:
        log.error("EGGNOG_DATA_DIR not set in .env — skipping eggNOG step")
        return
    eggnog_db = str(Path(data_dir_env).expanduser())
    if not Path(eggnog_db).exists():
        log.error(f"EGGNOG_DATA_DIR does not exist: {eggnog_db} — skipping")
        return

    emapper_bin = Path(sys.executable).parent / "emapper.py"
    emapper_cmd = [str(emapper_bin)] if emapper_bin.exists() else ["emapper.py"]

    results: list[tuple[str, str, str]] = []

    for g in genomes:
        strain = g["strain_name"]
        genome_dir = PROJECT_ROOT / g["data_dir"].rstrip("/")
        faa = genome_dir / "protein.faa"
        out_dir = genome_dir / "eggnog"
        annotation_file = out_dir / f"{strain}.emapper.annotations"

        if not faa.exists():
            results.append((strain, "SKIP_NO_FAA", "protein.faa not found — run step 1 first"))
            continue

        if annotation_file.exists() and not force:
            n = sum(
                1 for line in annotation_file.read_text().splitlines()
                if line and not line.startswith("#")
            )
            results.append((strain, "SKIP_EXISTS", f"{n} annotations already present"))
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = emapper_cmd + [
            "--data_dir", eggnog_db,
            "-i", str(faa),
            "--itype", "proteins",
            "-o", strain,
            "--output_dir", str(out_dir),
            "--cpu", str(cpu),
            "--override",
        ]
        log.info(f"  Running eggNOG-mapper for {strain} (cpu={cpu})...")
        env = os.environ.copy()
        env["EGGNOG_DATA_DIR"] = eggnog_db
        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT), env=env)

        if result.returncode == 0 and annotation_file.exists():
            n = sum(
                1 for line in annotation_file.read_text().splitlines()
                if line and not line.startswith("#")
            )
            results.append((strain, "OK", f"{n} proteins annotated"))
        else:
            results.append((strain, "FAILED", "emapper.py returned non-zero exit"))

    # Summary table
    if results:
        print(f"\n{'─'*60}")
        print(f"{'Strain':<15} {'Status':<15} {'Info'}")
        print(f"{'─'*15} {'─'*15} {'─'*30}")
        for strain, status, msg in results:
            print(f"{strain:<15} {status:<15} {msg}")
        print(f"{'─'*60}")


# ── main ─────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download genome annotation data (NCBI, Cyanorak, UniProt, eggNOG).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Steps:
  1  NCBI genome (GFF + protein FASTA + GBFF)
  2  Cyanorak GFF/GBK (strains with cyanorak_organism only)
  3  UniProt (one download per unique taxid)
  4  eggNOG-mapper (requires EGGNOG_DATA_DIR in .env)

Examples:
  uv run python scripts/download_genome_data.py
  uv run python scripts/download_genome_data.py --steps 1 2 3
  uv run python scripts/download_genome_data.py --strains MED4 MIT9313
  uv run python scripts/download_genome_data.py --strains MED4 --force
        """,
    )
    parser.add_argument(
        "--steps", nargs="+", type=int, choices=[1, 2, 3, 4],
        default=[1, 2, 3, 4],
        help="Steps to run (default: all). 1=NCBI 2=Cyanorak 3=UniProt 4=eggNOG",
    )
    parser.add_argument(
        "--strains", nargs="+",
        help="Restrict to specific strain names (e.g. --strains MED4 MIT9313)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if cache files already exist",
    )
    parser.add_argument(
        "--cpu", type=int, default=4,
        help="CPU count for eggNOG-mapper step 4 (default: 4)",
    )
    args = parser.parse_args()

    genomes = _read_genomes_csv(GENOMES_CSV)
    if not genomes:
        log.error(f"No genomes found in {GENOMES_CSV}")
        sys.exit(1)

    if args.strains:
        known = {g["strain_name"] for g in genomes}
        unknown = set(args.strains) - known
        if unknown:
            log.warning(f"Unknown strain(s) ignored: {sorted(unknown)}")
        genomes = [g for g in genomes if g["strain_name"] in args.strains]
        if not genomes:
            log.error(f"No genomes matched strain filter: {args.strains}")
            sys.exit(1)
        log.info(f"Filtered to {len(genomes)} genome(s): {[g['strain_name'] for g in genomes]}")
    else:
        log.info(f"Processing all {len(genomes)} genome(s)")

    steps = sorted(set(args.steps))
    log.info(f"Steps to run: {steps}{' (force)' if args.force else ''}")

    if 1 in steps:
        step1_ncbi(genomes, force=args.force)
    if 2 in steps:
        step2_cyanorak(genomes, force=args.force)
    if 3 in steps:
        step3_uniprot(genomes, force=args.force)
    if 4 in steps:
        step4_eggnog(genomes, force=args.force, cpu=args.cpu)

    log.info("All steps complete.")


if __name__ == "__main__":
    main()
