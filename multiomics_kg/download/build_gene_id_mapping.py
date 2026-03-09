#!/usr/bin/env python3
"""
Build per-strain gene_id_mapping.json (v2) from:
  1. gene_annotations_merged.json (base reference alt-IDs)
  2. paperconfig.yaml files (id_translation, annotation_gff, csv id_columns)

Uses iterative convergence via GeneIdGraph for true transitive closure.
No ordering dependency between sources — all are processed together.

Outputs per strain:
  cache/data/<org>/genomes/<strain>/gene_id_mapping.json    (v2)
  cache/data/<org>/genomes/<strain>/gene_id_mapping_report.json

Usage:
  uv run python -m multiomics_kg.download.build_gene_id_mapping [--strains STRAIN ...] [--force]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any
from urllib.parse import unquote

import pandas as pd
import yaml

from multiomics_kg.download.gene_id_graph import GeneIdGraph
from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.utils.gene_id_utils import get_genome_dir, load_gene_annotations

PAPERCONFIG_FILES_TXT = (
    PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
)

# ─── Paperconfig loading ──────────────────────────────────────────────────────


def load_all_paperconfigs() -> list[tuple[Path, dict]]:
    """Load all paperconfigs listed in paperconfig_files.txt."""
    results = []
    with open(PAPERCONFIG_FILES_TXT) as f:
        for line in f:
            path_str = line.strip()
            if not path_str or path_str.startswith('#'):
                continue
            path = PROJECT_ROOT / path_str
            if not path.exists():
                print(f"  [warn] paperconfig not found: {path}", file=sys.stderr)
                continue
            with open(path) as pf:
                config = yaml.safe_load(pf)
            if config:
                results.append((path, config))
    return results


def get_organism_for_entry(entry: dict) -> str | None:
    """Extract the organism string for a supplementary table entry.

    Priority: entry-level 'organism' > first statistical_analysis 'organism'.
    """
    org = entry.get("organism")
    if org:
        return str(org).strip().strip('"')
    for a in entry.get("statistical_analyses") or []:
        org = a.get("organism")
        if org:
            return str(org).strip().strip('"')
    return None


def collect_entries_for_genome_dir(
    paperconfigs: list[tuple[Path, dict]],
    target_genome_dir: Path,
) -> list[tuple[str, str, dict, Path]]:
    """Collect all supplementary table entries that match target_genome_dir.

    Returns list of (paper_name, table_key, entry_config, paperconfig_path).
    """
    results = []
    for pc_path, config in paperconfigs:
        pub = config.get("publication") or {}
        paper_name = pub.get("papername") or pc_path.parent.name
        supp = pub.get("supplementary_materials") or config.get("supplementary_materials") or {}
        for table_key, entry in supp.items():
            org = get_organism_for_entry(entry)
            if not org:
                continue
            resolved = get_genome_dir(org, str(PROJECT_ROOT))
            if resolved and Path(resolved).resolve() == target_genome_dir.resolve():
                results.append((paper_name, table_key, entry, pc_path))
    return results


# ─── Row extraction: source → list of (id_val, id_type) ──────────────────────


def _safe_str(val: Any) -> str:
    s = str(val).strip()
    if s.lower() == "nan":
        return ""
    return s


def extract_rows_from_id_translation(
    entry: dict,
    paper_name: str,
    table_key: str,
) -> list[tuple[list[tuple[str, str]], str]]:
    """Extract (id_pairs, source_label) rows from an id_translation entry."""
    filename = entry.get("filename", "")
    sep = entry.get("sep", ",")
    skip_rows = entry.get("skip_rows", 0)
    id_columns: list[dict] = entry.get("id_columns") or []

    if not id_columns:
        return []

    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] file not found: {path}", file=sys.stderr)
        return []

    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip_rows, dtype=str,
                         encoding="utf-8-sig")
    except Exception as e:
        print(f"    [warn] could not read {path}: {e}", file=sys.stderr)
        return []

    source_label = f"{paper_name}/{table_key}"
    result = []

    for _, row in df.iterrows():
        row_pairs: list[tuple[str, str]] = []
        for col_spec in id_columns:
            col = col_spec.get("column", "")
            id_type = col_spec.get("id_type", "other")
            if col not in df.columns:
                continue
            val = _safe_str(row.get(col, ""))
            if val:
                row_pairs.append((val, id_type))
        if row_pairs:
            result.append((row_pairs, source_label))

    return result


def _parse_gff_attrs(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def _parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def extract_rows_from_annotation_gff(
    entry: dict,
    paper_name: str,
    table_key: str,
) -> list[tuple[list[tuple[str, str]], str]]:
    """Extract rows from an annotation_gff entry.

    Each GFF feature row that has at least one recognized attribute becomes
    a source row linking (locus_tag, old_locus_tag, protein_id, Name, gene).
    """
    filename = entry.get("filename", "")
    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] GFF file not found: {path}", file=sys.stderr)
        return []

    is_gtf = path.suffix.lower() == ".gtf"
    source_label = f"{paper_name}/{table_key}"
    result: list[tuple[list[tuple[str, str]], str]] = []

    # GFF attribute → id_type mapping
    GFF_ATTR_TYPES = {
        "locus_tag": "locus_tag_ncbi",
        "old_locus_tag": "old_locus_tag",
        "protein_id": "protein_id_refseq",
        "Name": "locus_tag_ncbi",  # GFF Name attr is often the locus tag
        "gene": "gene_name",
    }
    VALID_FEATURES = {"gene", "CDS", "exon", "transcript"}

    try:
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] not in VALID_FEATURES:
                    continue
                attrs = _parse_gtf_attrs(parts[8]) if is_gtf else _parse_gff_attrs(parts[8])

                row_pairs: list[tuple[str, str]] = []
                for attr_name, id_type in GFF_ATTR_TYPES.items():
                    val = attrs.get(attr_name, "").strip()
                    if val:
                        row_pairs.append((val, id_type))

                # Extract "Alternative locus ID" from Note field (GCA GFF)
                # Format: Note=Alternative locus ID:P9313_01731
                # or Note=Alternative locus ID:P9313_00051%3B~Signal...
                note = attrs.get("Note", "")
                if "Alternative locus ID" in note:
                    note_decoded = unquote(note)
                    alt_match = re.search(
                        r"Alternative locus ID[: ]+([A-Za-z0-9_]+)",
                        note_decoded,
                    )
                    if alt_match:
                        row_pairs.append((alt_match.group(1), "alternative_locus_tag"))

                if len(row_pairs) >= 2:  # Need at least 2 IDs to be useful for linking
                    result.append((row_pairs, source_label))
    except Exception as e:
        print(f"    [warn] could not read GFF {path}: {e}", file=sys.stderr)

    return result


def extract_rows_from_cds_fna(genome_dir: Path) -> list[tuple[list[tuple[str, str]], str]]:
    """Extract rows from cds_from_genomic.fna FASTA headers.

    Each FASTA header encodes the CDS sequence ID plus bracketed attributes:
      >lcl|NC_007577.1_cds_WP_011375566.1_1 [gene=dnaN] [locus_tag=PMT9312_RS00005] [protein_id=WP_011375566.1]

    Yields one row per header containing:
      - cds_fna_id        : the full lcl| sequence identifier (Tier 1)
      - locus_tag_ncbi    : value of [locus_tag=...] attribute (Tier 1)
      - old_locus_tag     : value of [old_locus_tag=...] attribute (Tier 1, if present)
      - protein_id_refseq : value of [protein_id=...] attribute (Tier 2)
      - gene_name         : value of [gene=...] attribute (Tier 3)
    """
    fna_path = genome_dir / "cds_from_genomic.fna"
    if not fna_path.exists():
        return []

    source_label = "cds_from_genomic.fna"
    result: list[tuple[list[tuple[str, str]], str]] = []

    ATTR_TYPES = {
        "locus_tag": "locus_tag_ncbi",
        "old_locus_tag": "old_locus_tag",
        "protein_id": "protein_id_refseq",
        "gene": "gene_name",
    }

    try:
        with open(fna_path) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                header = line[1:].rstrip()
                parts = header.split(None, 1)
                seq_id = parts[0]
                rest = parts[1] if len(parts) > 1 else ""

                row_pairs: list[tuple[str, str]] = [(seq_id, "cds_fna_id")]

                for m in re.finditer(r"\[(\w+)=([^\]]+)\]", rest):
                    attr_name, attr_val = m.group(1), m.group(2).strip()
                    if attr_name in ATTR_TYPES and attr_val:
                        row_pairs.append((attr_val, ATTR_TYPES[attr_name]))

                if len(row_pairs) >= 2:
                    result.append((row_pairs, source_label))
    except Exception as e:
        print(f"    [warn] could not read {fna_path}: {e}", file=sys.stderr)

    return result



def extract_rows_from_csv_table(
    entry: dict,
    paper_name: str,
    table_key: str,
) -> list[tuple[list[tuple[str, str]], str]]:
    """Extract rows from a csv supplementary table (id_columns only, no DE data)."""
    id_columns: list[dict] = entry.get("id_columns") or []
    if not id_columns:
        return []

    filename = entry.get("original_filename") or entry.get("filename", "")
    sep = entry.get("sep", ",")
    skip_rows = entry.get("skip_rows", 0)
    analyses: list[dict] = entry.get("statistical_analyses") or []
    name_cols = list({a.get("name_col") for a in analyses if a.get("name_col")})

    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] csv file not found: {path}", file=sys.stderr)
        return []

    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip_rows, dtype=str,
                         encoding="utf-8-sig")
    except Exception as e:
        print(f"    [warn] could not read {path}: {e}", file=sys.stderr)
        return []

    source_label = f"{paper_name}/{table_key}"
    result: list[tuple[list[tuple[str, str]], str]] = []

    # Pre-filter id_columns to those actually present, and track which columns
    # are already covered by name_cols to avoid duplicate pairs
    name_col_set = set(name_cols)
    valid_id_columns = [
        cs for cs in id_columns
        if cs.get("column", "") in df.columns and cs.get("column", "") not in name_col_set
    ]

    for _, row in df.iterrows():
        row_pairs: list[tuple[str, str]] = []

        # Include name_col values with their declared or inferred id_type
        for nc in name_cols:
            if nc not in df.columns:
                continue
            val = _safe_str(row.get(nc, ""))
            if not val:
                continue
            nc_id_type = next(
                (c.get("id_type", "other") for c in id_columns if c.get("column") == nc),
                "locus_tag",  # name_col without explicit id_type is assumed locus_tag
            )
            row_pairs.append((val, nc_id_type))

        # Include declared id_columns not already covered by name_cols
        for col_spec in valid_id_columns:
            col = col_spec["column"]
            id_type = col_spec.get("id_type", "other")
            val = _safe_str(row.get(col, ""))
            if val:
                row_pairs.append((val, id_type))

        if row_pairs:
            result.append((row_pairs, source_label))

    return result


# ─── Seeding from annotations ─────────────────────────────────────────────────

# Fields in gene_annotations_merged.json → id_type mapping
_ANNOTATION_SCALAR_FIELDS: dict[str, str] = {
    "locus_tag_ncbi": "locus_tag_ncbi",
    "locus_tag_cyanorak": "locus_tag_cyanorak",
    "protein_id": "protein_id_refseq",
    "uniprot_accession": "uniprot_accession",
    "gene_name": "gene_name",
}

_ANNOTATION_LIST_FIELDS: dict[str, str] = {
    "old_locus_tags": "old_locus_tag",
    "alternative_locus_tags": "alternative_locus_tag",
    "gene_synonyms": "gene_synonym",
    "gene_name_synonyms": "gene_name",
}


def seed_graph_from_annotations(graph: GeneIdGraph, annotations: dict) -> int:
    """Seed the graph from gene_annotations_merged.json.

    Returns the number of genes seeded.
    """
    for locus_tag, entry in annotations.items():
        graph.add_anchor(locus_tag)

        for field, id_type in _ANNOTATION_SCALAR_FIELDS.items():
            val = entry.get(field)
            if val and isinstance(val, str):
                val = val.strip()
                if val:
                    graph.add_id_for_gene(locus_tag, val, id_type, "annotation")

        for field, id_type in _ANNOTATION_LIST_FIELDS.items():
            vals = entry.get(field)
            if isinstance(vals, list):
                for v in vals:
                    if isinstance(v, str) and v.strip():
                        graph.add_id_for_gene(locus_tag, v.strip(), id_type, "annotation")
            elif isinstance(vals, str) and vals.strip():
                for v in re.split(r",\s*", vals):
                    v = v.strip()
                    if v:
                        graph.add_id_for_gene(locus_tag, v, id_type, "annotation")

    return len(annotations)


# ─── Output writers ───────────────────────────────────────────────────────────


def write_gene_id_mapping(graph: GeneIdGraph, genome_dir: Path, organism: str, strain: str) -> None:
    """Write gene_id_mapping.json (v2) to genome_dir."""
    data = graph.to_json_structure(organism, strain)
    out_path = genome_dir / "gene_id_mapping.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    stats = data["stats"]
    print(
        f"  Wrote {out_path} "
        f"({stats['n_genes']} genes, {stats['n_specific']} specific, "
        f"{stats['n_multi']} multi, {stats['n_conflicts']} conflicts, "
        f"{stats['passes']} passes)"
    )


def write_diagnostic_report(graph: GeneIdGraph, genome_dir: Path, strain: str) -> None:
    """Write gene_id_mapping_report.json to genome_dir."""
    report = graph.build_diagnostic_report()
    report["strain"] = strain

    out_path = genome_dir / "gene_id_mapping_report.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    # Print warnings to stderr
    if report["warnings"]:
        print(f"  Diagnostic warnings for {strain}:", file=sys.stderr)
        for w in report["warnings"]:
            print(f"    {w}", file=sys.stderr)

    # Print unresolved summary
    unresolved = report.get("unresolved_rows_per_source", {})
    if unresolved:
        total_unresolved = sum(unresolved.values())
        print(f"  Unresolved rows (no anchor found): {total_unresolved} total")
        for src, n in sorted(unresolved.items(), key=lambda x: -x[1]):
            if n > 0:
                print(f"    {n} rows: {src}")

    print(f"  Wrote {out_path}")


# ─── Diamond protein match generation ────────────────────────────────────────


def generate_diamond_translations(
    entries: list[tuple[str, str, dict, Path]],
    genome_dir: Path,
    force: bool,
) -> None:
    """Pre-generate id_translation files that have a 'generate' block.

    Scans id_translation entries for a ``generate`` block with
    ``method: diamond_protein_match``.  When found, runs
    ``scripts/map_img_to_ncbi_proteins.py`` to produce the output CSV
    before ``extract_rows_from_id_translation`` tries to read it.
    """
    for paper_name, table_key, entry_config, _ in entries:
        if entry_config.get("type") != "id_translation":
            continue
        generate = entry_config.get("generate")
        if not generate:
            continue
        method = generate.get("method")
        if method != "diamond_protein_match":
            print(f"    [warn] unknown generate method: {method}")
            continue

        output_path = PROJECT_ROOT / entry_config["filename"]
        if output_path.exists() and not force:
            print(f"  [skip] {paper_name}/{table_key}: generated file exists")
            continue

        source_fasta = PROJECT_ROOT / generate["source_fasta"]
        source_id_col = generate.get("source_id_col", "source_id")
        ncbi_faa = genome_dir / "protein.faa"
        gene_mapping = genome_dir / "gene_mapping.csv"

        cmd = [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "map_img_to_ncbi_proteins.py"),
            "--img-faa", str(source_fasta),
            "--ncbi-faa", str(ncbi_faa),
            "--gene-mapping", str(gene_mapping),
            "--output", str(output_path),
            "--source-id-col", source_id_col,
        ]
        img_gff = generate.get("img_gff")
        if img_gff:
            cmd.extend(["--img-gff", str(PROJECT_ROOT / img_gff)])

        print(f"  [generate] {paper_name}/{table_key}: diamond protein match")
        print(f"    source: {source_fasta.name}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    [ERROR] diamond matching failed:\n{result.stderr}")
        else:
            for line in result.stdout.strip().split("\n")[-3:]:
                print(f"    {line}")


# ─── Per-strain orchestration ─────────────────────────────────────────────────


def process_strain(row: dict, paperconfigs: list, force: bool) -> None:
    strain = row["strain_name"]
    organism = row.get("organism_name", strain)
    genome_dir = (PROJECT_ROOT / row["data_dir"]).resolve()

    out_path = genome_dir / "gene_id_mapping.json"
    if out_path.exists() and not force:
        print(f"  [skip] {strain}: gene_id_mapping.json exists (use --force to overwrite)")
        return

    print(f"\n=== {strain} ===")

    annotations = load_gene_annotations(str(genome_dir))
    if annotations is None:
        print(f"  [warn] No gene_annotations_merged.json for {strain}, skipping")
        return

    graph = GeneIdGraph()
    n_seeded = seed_graph_from_annotations(graph, annotations)
    print(f"  Seeded {n_seeded} genes from annotations")
    print(f"  specific_lookup: {len(graph.specific_lookup)} entries, "
          f"multi_lookup: {len(graph.multi_lookup)} entries after seeding")

    # Auto-include genome-derived sources if present
    all_rows: list[tuple[list[tuple[str, str]], str]] = []

    # GCF cds_from_genomic.fna — maps current lcl| CDS IDs → current locus tags
    cds_rows = extract_rows_from_cds_fna(genome_dir)
    if cds_rows:
        all_rows.extend(cds_rows)
        print(f"  Loaded {len(cds_rows)} rows from cds_from_genomic.fna")
    else:
        print(f"  No cds_from_genomic.fna found for {strain} (run step 0 to download)")

    # GCA genomic_gca.gff — original GenBank annotation with old protein accessions
    # and old-style locus tags; reuses extract_rows_from_annotation_gff directly
    gca_gff = genome_dir / "genomic_gca.gff"
    if gca_gff.exists():
        gca_entry = {"filename": str(gca_gff.relative_to(PROJECT_ROOT))}
        gca_rows = extract_rows_from_annotation_gff(gca_entry, "genomic_gca", strain)
        if gca_rows:
            all_rows.extend(gca_rows)
            print(f"  Loaded {len(gca_rows)} rows from genomic_gca.gff")
    else:
        print(f"  No genomic_gca.gff found for {strain} (run step 0 to download)")

    # Collect all source rows from paperconfigs
    entries = collect_entries_for_genome_dir(paperconfigs, genome_dir)

    # Pre-generate any id_translation files that have a 'generate' block
    generate_diamond_translations(entries, genome_dir, force)

    if not entries:
        print(f"  No paperconfig entries for {strain}")
    else:
        for paper_name, table_key, entry_config, _ in entries:
            entry_type = entry_config.get("type", "csv")
            print(f"  [{entry_type}] {paper_name} / {table_key}")

            if entry_type == "id_translation":
                rows = extract_rows_from_id_translation(entry_config, paper_name, table_key)
                all_rows.extend(rows)
                print(f"    collected {len(rows)} rows")

            elif entry_type == "annotation_gff":
                rows = extract_rows_from_annotation_gff(entry_config, paper_name, table_key)
                all_rows.extend(rows)
                print(f"    collected {len(rows)} GFF rows")

            elif entry_type == "csv":
                rows = extract_rows_from_csv_table(entry_config, paper_name, table_key)
                all_rows.extend(rows)
                print(f"    collected {len(rows)} rows")

    if all_rows:
        print(f"  Processing {len(all_rows)} total rows (iterative convergence)...")
        passes = graph.process_all_rows(all_rows)
        print(f"  Converged after {passes} pass(es)")
        print(f"  specific_lookup: {len(graph.specific_lookup)} entries, "
              f"multi_lookup: {len(graph.multi_lookup)} entries after convergence")

    write_gene_id_mapping(graph, genome_dir, organism, strain)
    write_diagnostic_report(graph, genome_dir, strain)


# ─── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--strains", nargs="+", metavar="STRAIN", help="Only process these strains (default: all)"
    )
    parser.add_argument(
        "--force", action="store_true", help="Overwrite existing output files"
    )
    args = parser.parse_args()

    rows = load_genome_rows(args.strains)
    paperconfigs = load_all_paperconfigs()
    print(f"Loaded {len(paperconfigs)} paperconfigs")

    for row in rows:
        process_strain(row, paperconfigs, args.force)

    print("\nDone.")


if __name__ == "__main__":
    main()
