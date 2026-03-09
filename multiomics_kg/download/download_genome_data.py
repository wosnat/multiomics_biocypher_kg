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
  Step 5: Build gene_mapping.csv from GFF/GBK files (requires steps 1+2)
           → cache/data/<org>/genomes/<strain>/gene_mapping.csv

All steps skip existing cache files by default.
Use --force to re-download/re-run for specified strains.

Usage:
  uv run python multiomics_kg/download/download_genome_data.py
  uv run python multiomics_kg/download/download_genome_data.py --steps 1 2 3
  uv run python multiomics_kg/download/download_genome_data.py --strains MED4 MIT9313
  uv run python multiomics_kg/download/download_genome_data.py --strains MED4 --force
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
import sys
from contextlib import ExitStack
from pathlib import Path

import dotenv

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # multiomics_kg/download/ -> multiomics_kg/ -> project root
GENOMES_CSV = PROJECT_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"

dotenv.load_dotenv(PROJECT_ROOT / ".env")

log = logging.getLogger("download_genome_data")
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


# ── helpers ──────────────────────────────────────────────────────────────────

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


# ── standalone download functions ────────────────────────────────────────────

def _ncbi_download_genome(ncbi_accession: str, data_dir: str, force: bool) -> bool:
    """Download NCBI genome zip and extract GFF, protein FASTA, GBFF, and CDS FASTA.

    Returns True if files were downloaded, False if already cached.
    Raises ConnectionError if the download fails.
    """
    from pypath.share import curl

    gff_path = os.path.join(data_dir, "genomic.gff")
    prot_path = os.path.join(data_dir, "protein.faa")
    gbff_path = os.path.join(data_dir, "genomic.gbff")
    cds_fna_path = os.path.join(data_dir, "cds_from_genomic.fna")

    if not force and all(os.path.exists(p) for p in [gff_path, prot_path, gbff_path]):
        return False

    if force:
        for p in [gff_path, prot_path, gbff_path, cds_fna_path]:
            if os.path.exists(p):
                os.remove(p)
                log.debug(f"  Removed {p}")

    os.makedirs(data_dir, exist_ok=True)

    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{ncbi_accession}/download"
        f"?include_annotation_type=GENOME_GFF"
        f"&include_annotation_type=GENOME_FASTA"
        f"&include_annotation_type=PROT_FASTA"
        f"&include_annotation_type=GENOME_GBFF"
        f"&include_annotation_type=SEQUENCE_REPORT"
        f"&include_annotation_type=CDS_FASTA"
        f"&hydrated=FULLY_HYDRATED"
    )

    with ExitStack() as stack:
        if force:
            stack.enter_context(curl.cache_off())
        c = curl.Curl(url, silent=False, compr='zip', large=False)

    if c.result is None:
        raise ConnectionError(
            f"Failed to download NCBI genome for {ncbi_accession} from {url}"
        )

    for name, content in c.result.items():
        if name.endswith('genomic.gff'):
            with open(gff_path, 'w') as f:
                f.write(content)
            log.debug(f"  NCBI GFF saved to {gff_path}")
        elif name.endswith('protein.faa'):
            with open(prot_path, 'w') as f:
                f.write(content)
            log.debug(f"  NCBI protein FASTA saved to {prot_path}")
        elif name.endswith('genomic.gbff'):
            with open(gbff_path, 'w') as f:
                f.write(content)
            log.debug(f"  NCBI GBFF saved to {gbff_path}")
        elif name.endswith('cds_from_genomic.fna'):
            with open(cds_fna_path, 'w') as f:
                f.write(content)
            log.debug(f"  NCBI CDS FASTA saved to {cds_fna_path}")

    if not os.path.exists(cds_fna_path):
        log.warning(f"  cds_from_genomic.fna not found in NCBI zip for {ncbi_accession}")

    # Also download GCA (GenBank) GFF — has original protein accessions (ABB*.1, etc.)
    # and old-style locus tags, useful for resolving IDs from older publications.
    if ncbi_accession.startswith("GCF_"):
        gca_accession = "GCA_" + ncbi_accession[4:]
        gca_gff_path = os.path.join(data_dir, "genomic_gca.gff")
        if force and os.path.exists(gca_gff_path):
            os.remove(gca_gff_path)
        if not os.path.exists(gca_gff_path):
            gca_url = (
                f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
                f"{gca_accession}/download"
                f"?include_annotation_type=GENOME_GFF"
                f"&hydrated=FULLY_HYDRATED"
            )
            try:
                with ExitStack() as stack:
                    if force:
                        stack.enter_context(curl.cache_off())
                    c_gca = curl.Curl(gca_url, silent=False, compr='zip', large=False)
            except Exception as e:
                log.warning(f"  Failed to download GCA GFF for {gca_accession}: {e}")
                c_gca = None
            if c_gca is not None and c_gca.result is not None:
                for name, content in c_gca.result.items():
                    if name.endswith('genomic.gff'):
                        with open(gca_gff_path, 'w') as f:
                            f.write(content)
                        log.debug(f"  NCBI GCA GFF saved to {gca_gff_path}")
                        break
                if not os.path.exists(gca_gff_path):
                    log.warning(f"  genomic_gca.gff not found for {gca_accession}")
            elif c_gca is not None:
                log.warning(f"  Failed to download GCA GFF for {gca_accession}")

    return True


def _cyanorak_download_file(
    cyanorak_organism: str, data_dir: str, ext: str, force: bool
) -> bool:
    """Download a single Cyanorak annotation file (gff or gbk).

    Returns True if the file was downloaded, False if already cached.
    Raises ConnectionError if the download fails.
    """
    import requests

    cyan_dir = os.path.join(data_dir, "cyanorak")
    os.makedirs(cyan_dir, exist_ok=True)
    out_path = os.path.join(cyan_dir, f"{cyanorak_organism}.{ext}")

    if not force and os.path.exists(out_path):
        return False

    if force and os.path.exists(out_path):
        os.remove(out_path)
        log.debug(f"  Removed {out_path}")

    url = f"https://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/{ext}/{cyanorak_organism}"

    response = requests.get(url, timeout=120)
    if not response.ok:
        raise ConnectionError(
            f"Failed to download Cyanorak {ext} for {cyanorak_organism}: "
            f"HTTP {response.status_code} from {url}"
        )

    with open(out_path, 'w') as f:
        f.write(response.text)
    log.debug(f"  Cyanorak {ext.upper()} saved to {out_path}")

    return True




# ── step 1: NCBI genome ──────────────────────────────────────────────────────

def step1_ncbi(genomes: list[dict], force: bool) -> None:
    """Download NCBI genome files (GFF, protein FASTA, GBFF) for each strain."""
    _header(1, "NCBI genome download")
    for g in genomes:
        strain = g["strain_name"]
        data_dir = str(PROJECT_ROOT / g["data_dir"].rstrip("/"))
        log.info(f"  Processing {strain} ({g['ncbi_accession']})...")
        if _ncbi_download_genome(g["ncbi_accession"], data_dir, force):
            log.info(f"  OK: {strain}")
        else:
            log.info(f"  [SKIP] {strain} — NCBI files already cached")


# ── step 2: Cyanorak ─────────────────────────────────────────────────────────

def step2_cyanorak(genomes: list[dict], force: bool) -> None:
    """Download Cyanorak GFF and GBK annotation files."""
    _header(2, "Cyanorak GFF/GBK download")
    for g in genomes:
        strain = g["strain_name"]
        cyan_org = g.get("cyanorak_organism", "").strip()
        if not cyan_org:
            log.info(f"  [SKIP] {strain} — no cyanorak_organism in config")
            continue
        data_dir = str(PROJECT_ROOT / g["data_dir"].rstrip("/"))
        gff_new = _cyanorak_download_file(cyan_org, data_dir, "gff", force)
        gbk_new = _cyanorak_download_file(cyan_org, data_dir, "gbk", force)
        if gff_new or gbk_new:
            log.info(f"  OK: {strain}")
        else:
            log.info(f"  [SKIP] {strain} — Cyanorak files already cached")


# ── step 3: UniProt ──────────────────────────────────────────────────────────

def step3_uniprot(genomes: list[dict], force: bool) -> None:
    """Download UniProt data per unique (org_group, taxid) pair."""
    _header(3, "UniProt data download (per taxid)")

    # Deduplicate: one download per (org_group, taxid)
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
        log.info(
            f"  Processing UniProt taxid={taxid} ({org_group}); "
            f"representative strain: {strain}"
        )
        from multiomics_kg.download.download_uniprot import download_uniprot
        if download_uniprot(taxid, cache_dir, force):
            log.info(f"  OK: UniProt taxid={taxid} → {cache_dir}")
        else:
            log.info(f"  [SKIP] UniProt taxid={taxid} ({org_group}) — cache exists")


# ── step 4: eggNOG-mapper ────────────────────────────────────────────────────

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


# ── step 5: gene_mapping.csv ─────────────────────────────────────────────────

def step5_gene_mapping(genomes: list[dict], force: bool) -> None:
    """Build gene_mapping.csv from downloaded GFF/GBK files.

    Requires steps 1 (NCBI) and 2 (Cyanorak) to have run first.
    Skips strains where gene_mapping.csv already exists (unless --force).
    """
    from multiomics_kg.download.build_gene_mapping import build_gene_mapping

    _header(5, "gene_mapping.csv generation")
    results: list[tuple[str, str, str]] = []

    for g in genomes:
        strain = g["strain_name"]
        data_dir = str(PROJECT_ROOT / g["data_dir"].rstrip("/"))
        mapping_path = os.path.join(data_dir, "gene_mapping.csv")

        if not force and os.path.exists(mapping_path):
            results.append((strain, "SKIP_EXISTS", "gene_mapping.csv already present"))
            continue

        ncbi_gff = os.path.join(data_dir, "genomic.gff")
        if not os.path.exists(ncbi_gff):
            results.append((strain, "SKIP_NO_GFF", "genomic.gff not found — run step 1 first"))
            continue

        # Resolve optional files (pass None if absent so build_gene_mapping skips them)
        ncbi_gbff = os.path.join(data_dir, "genomic.gbff")
        ncbi_gbff = ncbi_gbff if os.path.exists(ncbi_gbff) else None

        cyan_org = (g.get("cyanorak_organism") or "").strip() or None
        cyan_gff = os.path.join(data_dir, "cyanorak", f"{cyan_org}.gff") if cyan_org else None
        cyan_gbk = os.path.join(data_dir, "cyanorak", f"{cyan_org}.gbk") if cyan_org else None
        if cyan_gff and not os.path.exists(cyan_gff):
            cyan_gff = None
        if cyan_gbk and not os.path.exists(cyan_gbk):
            cyan_gbk = None

        log.info(f"  Building gene_mapping.csv for {strain}...")
        try:
            df = build_gene_mapping(
                ncbi_gff_file=ncbi_gff,
                ncbi_gbff_file=ncbi_gbff,
                cyan_gff_file=cyan_gff,
                cyan_gbk_file=cyan_gbk,
            )
            df.to_csv(mapping_path, index=False)
            log.info(f"  Gene mapping table saved to {mapping_path}")
            results.append((strain, "OK", f"{len(df)} genes"))
        except Exception as e:
            # Remove partially written file to avoid treating it as valid cache
            if os.path.exists(mapping_path):
                os.remove(mapping_path)
            log.debug(f"  Traceback for {strain}:", exc_info=True)
            results.append((strain, "FAILED", str(e)[:80]))

    if results:
        print(f"\n{'─'*60}")
        print(f"{'Strain':<15} {'Status':<15} {'Info'}")
        print(f"{'─'*15} {'─'*15} {'─'*30}")
        for strain, status, msg in results:
            print(f"{strain:<15} {status:<15} {msg}")
        print(f"{'─'*60}")


# ── main ─────────────────────────────────────────────────────────────────────

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
  5  gene_mapping.csv (requires steps 1+2)

Examples:
  uv run python multiomics_kg/download/download_genome_data.py
  uv run python multiomics_kg/download/download_genome_data.py --steps 1 2 3
  uv run python multiomics_kg/download/download_genome_data.py --strains MED4 MIT9313
  uv run python multiomics_kg/download/download_genome_data.py --strains MED4 --force
        """,
    )
    parser.add_argument(
        "--steps", nargs="+", type=int, choices=[1, 2, 3, 4, 5],
        default=[1, 2, 3, 4, 5],
        help="Steps to run (default: all). 1=NCBI 2=Cyanorak 3=UniProt 4=eggNOG 5=gene_mapping",
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
    if 5 in steps:
        step5_gene_mapping(genomes, force=args.force)

    log.info("All steps complete.")


if __name__ == "__main__":
    main()
