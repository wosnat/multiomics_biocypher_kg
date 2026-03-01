#!/usr/bin/env python3
"""
Build per-strain gene_id_mapping.json from:
  1. gene_annotations_merged.json (base reference alt-IDs)
  2. paperconfig.yaml files (id_translation, annotation_gff, csv id_columns)

Also writes backward-compatible gene_mapping_supp.csv.

Outputs per strain:
  cache/data/<org>/genomes/<strain>/gene_id_mapping.json
  cache/data/<org>/genomes/<strain>/gene_mapping_supp.csv

Usage:
  uv run python multiomics_kg/download/build_gene_id_mapping.py [--strains STRAIN ...] [--force]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.utils.gene_id_utils import get_genome_dir, load_gene_annotations

PAPERCONFIG_FILES_TXT = (
    PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
)

# ─── Scope inference ──────────────────────────────────────────────────────────

GENERIC_ID_TYPES = {"gene_name", "gene_synonym"}
ANNOTATION_SPECIFIC_ID_TYPES = {"rast_id", "annotation_specific"}


def infer_scope(id_type: str) -> str:
    if id_type in GENERIC_ID_TYPES:
        return "generic"
    if id_type in ANNOTATION_SPECIFIC_ID_TYPES:
        return "annotation_specific"
    return "organism_specific"


# ─── Paperconfig loading ──────────────────────────────────────────────────────


def load_all_paperconfigs() -> list[tuple[Path, dict]]:
    """Load all paperconfigs listed in paperconfig_files.txt."""
    results = []
    with open(PAPERCONFIG_FILES_TXT) as f:
        for line in f:
            path_str = line.strip()
            if not path_str:
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
        # Support both paper-level and strain-resource (no publication block) configs
        supp = pub.get("supplementary_materials") or config.get("supplementary_materials") or {}
        for table_key, entry in supp.items():
            org = get_organism_for_entry(entry)
            if not org:
                continue
            resolved = get_genome_dir(org, str(PROJECT_ROOT))
            if resolved and Path(resolved).resolve() == target_genome_dir.resolve():
                results.append((paper_name, table_key, entry, pc_path))
    return results


# ─── Mapping data structure ───────────────────────────────────────────────────


def _make_alt_id(id_val: str, id_type: str, **kwargs) -> dict:
    entry = {"id": id_val, "id_type": id_type, "scope": infer_scope(id_type)}
    entry.update(kwargs)
    return entry


def seed_from_annotations(annotations: dict) -> dict:
    """Build initial gene_id_mapping from gene_annotations_merged.json."""
    SCALAR_FIELDS = {
        "locus_tag_ncbi": "locus_tag_ncbi",
        "locus_tag_cyanorak": "locus_tag_cyanorak",
        "protein_id": "protein_id_refseq",
        "uniprot_accession": "uniprot_accession",
        "gene_name": "gene_name",
    }
    LIST_FIELDS = {
        "old_locus_tags": "old_locus_tag",
        "alternative_locus_tags": "alternative_locus_tag",
        "gene_synonyms": "gene_synonym",
        "gene_name_synonyms": "gene_name",
    }

    mapping = {}
    for locus_tag, entry in annotations.items():
        ref_ids = []

        for field, id_type in SCALAR_FIELDS.items():
            val = entry.get(field)
            if val and isinstance(val, str):
                val = val.strip()
                if val:
                    ref_ids.append(_make_alt_id(val, id_type))

        for field, id_type in LIST_FIELDS.items():
            vals = entry.get(field)
            if isinstance(vals, list):
                for v in vals:
                    if isinstance(v, str) and v.strip():
                        ref_ids.append(_make_alt_id(v.strip(), id_type))
            elif isinstance(vals, str) and vals.strip():
                for v in re.split(r",\s*", vals):
                    v = v.strip()
                    if v:
                        ref_ids.append(_make_alt_id(v, id_type))

        mapping[locus_tag] = {
            "locus_tag": locus_tag,
            "alt_ids": {"reference": ref_ids, "paper_ids": []},
            "product_synonyms": [],
        }
    return mapping


def build_reverse_lookup(mapping: dict) -> dict[str, str]:
    """Build alt_id → locus_tag reverse index from current mapping."""
    lookup: dict[str, str] = {}
    for locus_tag, entry in mapping.items():
        lookup[locus_tag] = locus_tag
        for alt in entry["alt_ids"]["reference"]:
            lookup[alt["id"]] = locus_tag
        for alt in entry["alt_ids"]["paper_ids"]:
            lookup[alt["id"]] = locus_tag
    return lookup


# ─── ID resolution ────────────────────────────────────────────────────────────


def _strip_uniprot_entry_suffix(value: str) -> str:
    """Strip trailing _ORGANISM suffix from UniProt entry name.

    "DNAA_PROM0" → "DNAA"
    "A3PBU0_PROM0" → "A3PBU0"
    """
    idx = value.rfind("_")
    if idx > 0:
        return value[:idx]
    return value


def resolve_value(value: str, id_type: str, lookup: dict) -> tuple[str | None, str]:
    """Try to resolve a value to a locus_tag.

    Returns (locus_tag, method) or (None, '').

    Strategies:
    - Direct lookup
    - uniprot_entry_name: strip _ORGANISM suffix
    - Whitespace-split fallback (handles "gene_name locus_tag" compound values)
    """
    value = str(value).strip()
    if not value or value in ("nan", ""):
        return None, ""

    if value in lookup:
        return lookup[value], "direct"

    if id_type == "uniprot_entry_name":
        stripped = _strip_uniprot_entry_suffix(value)
        if stripped and stripped != value and stripped in lookup:
            return lookup[stripped], "uniprot_entry_stripped"

    if " " in value:
        for token in value.split():
            if token in lookup:
                return lookup[token], "whitespace_split"

    return None, ""


# ─── GFF3 parsing ─────────────────────────────────────────────────────────────


def _parse_gff_attrs(attr_str: str) -> dict[str, str]:
    attrs = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def parse_annotation_gff(gff_path: Path, lookup: dict) -> list[tuple[str, str, str, str]]:
    """Parse a GFF3/GTF file and return novel ID bridges.

    Returns list of (locus_tag, id_val, id_type, attr_name) for IDs not already
    in the lookup (i.e., would be novel additions as paper_ids).
    """
    results: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    is_gtf = gff_path.suffix.lower() == ".gtf"

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ("gene", "CDS", "exon", "transcript"):
                continue

            if is_gtf:
                attrs = _parse_gtf_attrs(parts[8])
            else:
                attrs = _parse_gff_attrs(parts[8])

            locus_tag_attr = attrs.get("locus_tag", "")
            old_locus_tag_attr = attrs.get("old_locus_tag", "")
            name_attr = attrs.get("Name", attrs.get("gene_id", ""))
            gene_name_attr = attrs.get("gene", "")
            protein_id_attr = attrs.get("protein_id", "")

            # Find canonical locus_tag anchor
            anchor = None
            for val in (locus_tag_attr, old_locus_tag_attr, name_attr, gene_name_attr):
                if val and val in lookup:
                    anchor = lookup[val]
                    break
            if not anchor:
                continue

            # Record novel (not yet in lookup) ID bridges
            for val, id_type, attr_name in [
                (locus_tag_attr, "locus_tag_ncbi", "locus_tag"),
                (old_locus_tag_attr, "old_locus_tag", "old_locus_tag"),
                (protein_id_attr, "protein_id_refseq", "protein_id"),
                (name_attr, "gene_name", "Name"),
                (gene_name_attr, "gene_name", "gene"),
            ]:
                if val and val not in lookup:
                    key = (anchor, val, id_type)
                    if key not in seen:
                        seen.add(key)
                        results.append((anchor, val, id_type, attr_name))

    return results


def _parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    """Parse GTF attribute string (semicolon-separated key "value" pairs)."""
    attrs = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


# ─── Paper ID helpers ─────────────────────────────────────────────────────────


def _add_paper_id(
    mapping: dict,
    locus_tag: str,
    id_val: str,
    id_type: str,
    source_col: str,
    source_csv: str,
    paper: str,
) -> None:
    """Add a paper_id entry to mapping[locus_tag] if not already present."""
    if locus_tag not in mapping:
        return
    paper_ids = mapping[locus_tag]["alt_ids"]["paper_ids"]
    key = (id_val, id_type, source_col, source_csv, paper)
    existing = {
        (e["id"], e["id_type"], e.get("source_col"), e.get("source_csv"), e.get("paper"))
        for e in paper_ids
    }
    if key not in existing:
        paper_ids.append(
            _make_alt_id(id_val, id_type, source_col=source_col, source_csv=source_csv, paper=paper)
        )


def _add_product_synonym(
    mapping: dict, locus_tag: str, text: str, source_col: str, paper: str
) -> None:
    if locus_tag not in mapping or not text or text == "nan":
        return
    synonyms = mapping[locus_tag]["product_synonyms"]
    key = (text, source_col, paper)
    existing = {(s["text"], s.get("source_col"), s.get("paper")) for s in synonyms}
    if key not in existing:
        synonyms.append({"text": text, "source_col": source_col, "paper": paper})


# ─── Table processors ─────────────────────────────────────────────────────────


def process_id_translation(
    entry: dict,
    mapping: dict,
    lookup: dict,
    paper_name: str,
    source_csv_name: str,
) -> int:
    """Process an id_translation table. Returns the number of resolved rows."""
    filename = entry.get("filename", "")
    sep = entry.get("sep", ",")
    skip_rows = entry.get("skip_rows", 0)
    id_columns: list[dict] = entry.get("id_columns") or []
    product_columns: list[dict] = entry.get("product_columns") or []

    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] file not found: {path}", file=sys.stderr)
        return 0

    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip_rows, dtype=str)
    except Exception as e:
        print(f"    [warn] could not read {path}: {e}", file=sys.stderr)
        return 0

    source_csv = path.name
    resolved_count = 0

    for _, row in df.iterrows():
        anchor: str | None = None

        # Try each id_column in order to find a resolving anchor
        for col_spec in id_columns:
            col = col_spec.get("column", "")
            id_type = col_spec.get("id_type", "other")
            if col not in df.columns:
                continue
            val = str(row.get(col, "")).strip()
            lt, _ = resolve_value(val, id_type, lookup)
            if lt:
                anchor = lt
                break

        if not anchor:
            continue
        resolved_count += 1

        # Record all id_column values as paper_ids for this anchor
        for col_spec in id_columns:
            col = col_spec.get("column", "")
            id_type = col_spec.get("id_type", "other")
            if col not in df.columns:
                continue
            val = str(row.get(col, "")).strip()
            if not val or val == "nan":
                continue
            # For compound "gene_name locus_tag" values, record each token separately
            if " " in val:
                for token in val.split():
                    _add_paper_id(mapping, anchor, token, id_type, col, source_csv, paper_name)
            else:
                _add_paper_id(mapping, anchor, val, id_type, col, source_csv, paper_name)

        # Record product synonyms
        for col_spec in product_columns:
            col = col_spec.get("column", "")
            if col not in df.columns:
                continue
            val = str(row.get(col, "")).strip()
            _add_product_synonym(mapping, anchor, val, col, paper_name)

    return resolved_count


def process_annotation_gff_entry(
    entry: dict,
    mapping: dict,
    lookup: dict,
    paper_name: str,
) -> int:
    """Process an annotation_gff entry. Returns number of novel bridges added."""
    filename = entry.get("filename", "")
    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] GFF file not found: {path}", file=sys.stderr)
        return 0

    novel_bridges = parse_annotation_gff(path, lookup)
    source_csv = path.name
    for locus_tag, id_val, id_type, attr_name in novel_bridges:
        _add_paper_id(mapping, locus_tag, id_val, id_type, attr_name, source_csv, paper_name)
    return len(novel_bridges)


def process_csv_table(
    entry: dict,
    mapping: dict,
    lookup: dict,
    paper_name: str,
    source_csv_name: str,
) -> int:
    """Process a csv table for id_columns/product_columns metadata only.

    Does not emit expression edges — that is omics_adapter's job.
    Returns number of rows where an anchor locus_tag was found.
    """
    id_columns: list[dict] = entry.get("id_columns") or []
    product_columns: list[dict] = entry.get("product_columns") or []
    if not id_columns and not product_columns:
        return 0

    # Prefer original_filename (pre-fix-gene-ids) if declared
    filename = entry.get("original_filename") or entry.get("filename", "")
    sep = entry.get("sep", ",")
    skip_rows = entry.get("skip_rows", 0)
    analyses: list[dict] = entry.get("statistical_analyses") or []
    name_cols = list({a.get("name_col") for a in analyses if a.get("name_col")})

    path = PROJECT_ROOT / filename
    if not path.exists():
        print(f"    [warn] csv file not found: {path}", file=sys.stderr)
        return 0

    try:
        df = pd.read_csv(path, sep=sep, skiprows=skip_rows, dtype=str)
    except Exception as e:
        print(f"    [warn] could not read {path}: {e}", file=sys.stderr)
        return 0

    source_csv = path.name
    resolved_count = 0

    for _, row in df.iterrows():
        anchor: str | None = None

        # Try name_col first
        for nc in name_cols:
            if nc not in df.columns:
                continue
            val = str(row.get(nc, "")).strip()
            if not val or val == "nan":
                continue
            if nc == "locus_tag":
                if val in mapping:
                    anchor = val
                    break
            else:
                id_type = next(
                    (c.get("id_type", "other") for c in id_columns if c.get("column") == nc),
                    "other",
                )
                lt, _ = resolve_value(val, id_type, lookup)
                if lt:
                    anchor = lt
                    break

        # Fallback: try id_columns
        if not anchor:
            for col_spec in id_columns:
                col = col_spec.get("column", "")
                id_type = col_spec.get("id_type", "other")
                if col not in df.columns:
                    continue
                val = str(row.get(col, "")).strip()
                lt, _ = resolve_value(val, id_type, lookup)
                if lt:
                    anchor = lt
                    break

        if not anchor:
            continue
        resolved_count += 1

        for col_spec in id_columns:
            col = col_spec.get("column", "")
            id_type = col_spec.get("id_type", "other")
            if col not in df.columns:
                continue
            val = str(row.get(col, "")).strip()
            if not val or val == "nan":
                continue
            if " " in val:
                for token in val.split():
                    _add_paper_id(mapping, anchor, token, id_type, col, source_csv, paper_name)
            else:
                _add_paper_id(mapping, anchor, val, id_type, col, source_csv, paper_name)

        for col_spec in product_columns:
            col = col_spec.get("column", "")
            if col not in df.columns:
                continue
            val = str(row.get(col, "")).strip()
            _add_product_synonym(mapping, anchor, val, col, paper_name)

    return resolved_count


# ─── Output writers ───────────────────────────────────────────────────────────


def write_gene_id_mapping(mapping: dict, genome_dir: Path) -> None:
    out_path = genome_dir / "gene_id_mapping.json"
    with open(out_path, "w") as f:
        json.dump(mapping, f, indent=2)
    n_with_paper_ids = sum(1 for e in mapping.values() if e["alt_ids"]["paper_ids"])
    print(f"  Wrote {out_path} ({len(mapping)} genes, {n_with_paper_ids} with paper_ids)")


def write_gene_mapping_supp_csv(mapping: dict, genome_dir: Path) -> None:
    """Write backward-compatible gene_mapping_supp.csv."""
    out_path = genome_dir / "gene_mapping_supp.csv"
    rows: list[dict] = []
    for locus_tag, entry in mapping.items():
        for pid in entry["alt_ids"]["paper_ids"]:
            rows.append({
                "alt_id": pid["id"],
                "locus_tag": locus_tag,
                "source_col": pid.get("source_col", ""),
                "source_csv": pid.get("source_csv", ""),
                "paper": pid.get("paper", ""),
            })
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["alt_id", "locus_tag", "source_col", "source_csv", "paper"]
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {out_path} ({len(rows)} alt-ID rows)")


# ─── Per-strain orchestration ─────────────────────────────────────────────────


def process_strain(row: dict, paperconfigs: list, force: bool) -> None:
    strain = row["strain_name"]
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

    mapping = seed_from_annotations(annotations)
    lookup = build_reverse_lookup(mapping)
    print(f"  Seeded {len(mapping)} genes; lookup has {len(lookup)} entries")

    entries = collect_entries_for_genome_dir(paperconfigs, genome_dir)
    if not entries:
        print(f"  No paperconfig entries for {strain}")
    else:
        # Process id_translation first, then annotation_gff, then csv
        type_order = {"id_translation": 0, "annotation_gff": 1, "csv": 2}
        entries.sort(key=lambda x: type_order.get(x[2].get("type", "csv"), 2))

        for paper_name, table_key, entry_config, pc_path in entries:
            entry_type = entry_config.get("type", "csv")
            source_csv_name = Path(entry_config.get("filename", table_key)).name
            print(f"  [{entry_type}] {paper_name} / {table_key}")

            if entry_type == "id_translation":
                n = process_id_translation(
                    entry_config, mapping, lookup, paper_name, source_csv_name
                )
                print(f"    resolved {n} rows")
                # Rebuild lookup so subsequent tables benefit from new paper_ids
                lookup = build_reverse_lookup(mapping)
                print(f"    lookup now has {len(lookup)} entries")

            elif entry_type == "annotation_gff":
                n = process_annotation_gff_entry(entry_config, mapping, lookup, paper_name)
                print(f"    added {n} novel GFF bridges")

            elif entry_type == "csv":
                n = process_csv_table(
                    entry_config, mapping, lookup, paper_name, source_csv_name
                )
                print(f"    harvested alt-IDs from {n} rows")

    write_gene_id_mapping(mapping, genome_dir)
    write_gene_mapping_supp_csv(mapping, genome_dir)


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
