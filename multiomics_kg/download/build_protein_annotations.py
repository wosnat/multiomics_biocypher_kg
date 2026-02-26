#!/usr/bin/env python3
"""
Build per-taxid protein annotation tables from UniProt raw data.

Similar to build_gene_annotations.py but:
  - Single source: only uniprot_raw_data.json
  - Keyed by UniProt accession (not locus_tag)
  - Output per unique taxid (not per strain)

Input:  cache/data/<org_group>/uniprot/<taxid>/uniprot_raw_data.json
Output: cache/data/<org_group>/uniprot/<taxid>/protein_annotations.json

Usage:
  uv run python multiomics_kg/download/build_protein_annotations.py [--force]
  uv run python multiomics_kg/download/build_protein_annotations.py --strains MED4 EZ55
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import Any

from multiomics_kg.download.utils.annotation_helpers import (
    _coerce_to_tokens,
    _nonempty,
    _split,
)
from multiomics_kg.download.utils.annotation_transforms import _TRANSFORMS
from multiomics_kg.download.utils.cli import add_common_args, load_config
from multiomics_kg.download.utils.paths import GENOMES_CSV, PROJECT_ROOT, infer_organism_group

# ─── paths ────────────────────────────────────────────────────────────────────

DEFAULT_CONFIG = PROJECT_ROOT / "config/protein_annotations_config.yaml"


# ─── data loader ──────────────────────────────────────────────────────────────

def load_uniprot_columnar(path: str) -> dict[str, dict]:
    """Load column-oriented uniprot_raw_data.json → {uniprot_acc: {field: value}}.

    Stays keyed by UniProt accession (unlike build_gene_annotations.load_uniprot
    which re-indexes by RefSeq accession for gene-level joining).
    """
    with open(path) as f:
        col_data: dict[str, dict] = json.load(f)

    # Gather all UniProt IDs across all fields
    all_ids: set[str] = set()
    for col_values in col_data.values():
        if isinstance(col_values, dict):
            all_ids.update(col_values.keys())

    # Pivot to row-oriented dict keyed by UniProt accession
    rows: dict[str, dict] = {uid: {} for uid in all_ids}
    for field, col_values in col_data.items():
        if isinstance(col_values, dict):
            for uid, val in col_values.items():
                rows[uid][field] = val

    return rows


# ─── annotation builder ───────────────────────────────────────────────────────

class ProteinAnnotationBuilder:
    """Single-source builder: transforms one UniProt row into canonical protein fields."""

    def __init__(self, config: dict):
        self.field_configs: dict[str, dict] = config.get("fields", {})

    def _get_raw(self, fconf: dict, up: dict) -> Any:
        """Fetch raw value from the UniProt row."""
        field = fconf.get("field", "")
        raw = up.get(field)
        if not _nonempty(raw):
            return None
        return raw

    def _apply_transform(self, transform: str | None, value: Any) -> Any:
        if not transform or transform not in _TRANSFORMS:
            return value
        fn = _TRANSFORMS[transform]
        if isinstance(value, list):
            return [fn(v) for v in value if _nonempty(v)]
        return fn(value)

    def _resolve_passthrough(self, fconf: dict, up: dict) -> Any:
        raw = self._get_raw(fconf, up)
        if not _nonempty(raw):
            return None
        if isinstance(raw, str):
            # Sanitize characters that break Neo4j CSV import
            raw = raw.replace("|", ",").replace("'", "^")
        transform = fconf.get("transform")
        if transform:
            raw = self._apply_transform(transform, raw)
        return raw if _nonempty(raw) else None

    def _resolve_passthrough_list(self, fconf: dict, up: dict) -> list[str] | None:
        raw = self._get_raw(fconf, up)
        if not _nonempty(raw):
            return None
        delimiter = fconf.get("delimiter", ";")
        tokens = _coerce_to_tokens(raw, delimiter)
        transform = fconf.get("transform")
        if transform and transform in _TRANSFORMS:
            fn = _TRANSFORMS[transform]
            tokens = [fn(t) for t in tokens if _nonempty(t)]
        tokens = [t for t in tokens if _nonempty(t)]
        return tokens if tokens else None

    def _resolve_integer(self, fconf: dict, up: dict) -> int | None:
        raw = self._get_raw(fconf, up)
        if raw is None:
            return None
        try:
            return int(float(str(raw).replace(",", "").strip()))
        except (ValueError, TypeError):
            return None

    def _resolve_float(self, fconf: dict, up: dict) -> float | None:
        raw = self._get_raw(fconf, up)
        if raw is None:
            return None
        try:
            return float(str(raw).strip())
        except (ValueError, TypeError):
            return None

    def _resolve_bool(self, fconf: dict, up: dict) -> bool | None:
        raw = self._get_raw(fconf, up)
        if raw is None:
            return None
        s = str(raw).strip().lower()
        if s in ("reviewed", "true", "yes", "1"):
            return True
        if s in ("unreviewed", "false", "no", "0"):
            return False
        return None

    def build_merged(self, uid: str, up: dict) -> dict:
        """Apply field rules from config → canonical field set (sparse output)."""
        result: dict[str, Any] = {}

        for canonical_field, fconf in self.field_configs.items():
            ftype = fconf.get("type", "passthrough")

            if ftype == "passthrough":
                val = self._resolve_passthrough(fconf, up)
            elif ftype == "passthrough_list":
                val = self._resolve_passthrough_list(fconf, up)
            elif ftype == "integer":
                val = self._resolve_integer(fconf, up)
            elif ftype == "float":
                val = self._resolve_float(fconf, up)
            elif ftype == "bool":
                val = self._resolve_bool(fconf, up)
            else:
                continue

            # Keep False (bool) even though it's falsy; skip None and empty
            if val is not None and (_nonempty(val) or val is False):
                result[canonical_field] = val

        return result


# ─── per-taxid pipeline ───────────────────────────────────────────────────────

def process_taxid(
    org_group: str,
    ncbi_taxon_id: int,
    config: dict,
    force: bool = False,
) -> None:
    """Process one taxid's uniprot_raw_data.json → protein_annotations.json."""
    uniprot_dir = (
        PROJECT_ROOT / "cache" / "data" / org_group
        / "uniprot" / str(ncbi_taxon_id)
    )
    input_path = uniprot_dir / "uniprot_raw_data.json"
    output_path = uniprot_dir / "protein_annotations.json"

    if not input_path.exists():
        print(f"  [{org_group}/{ncbi_taxon_id}] No uniprot_raw_data.json — skipping")
        return

    if not force and output_path.exists():
        print(
            f"  [{org_group}/{ncbi_taxon_id}] Skipping (already exists). "
            "Use --force to rebuild."
        )
        return

    print(f"\n[{org_group}/{ncbi_taxon_id}] Loading {input_path} ...")
    rows = load_uniprot_columnar(str(input_path))
    print(f"  Loaded {len(rows)} UniProt entries")

    builder = ProteinAnnotationBuilder(config)
    output: dict[str, dict] = {}
    stats = dict(total=0, has_go=0, has_ec=0, has_function=0, reviewed=0)

    for uid, up_row in rows.items():
        merged = builder.build_merged(uid, up_row)
        if not merged:
            continue
        output[uid] = merged
        stats["total"] += 1
        if (merged.get("go_cellular_components")
                or merged.get("go_biological_processes")
                or merged.get("go_molecular_functions")):
            stats["has_go"] += 1
        if merged.get("ec_numbers"):
            stats["has_ec"] += 1
        if merged.get("function_description"):
            stats["has_function"] += 1
        if merged.get("is_reviewed"):
            stats["reviewed"] += 1

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  → {output_path}")

    n = stats["total"] or 1
    pct = lambda k: f"{100 * stats[k] // n}%"
    print(f"\n  === {org_group}/{ncbi_taxon_id} protein annotation coverage ===")
    print(f"  Proteins:       {stats['total']}")
    print(f"  Reviewed:       {stats['reviewed']} ({pct('reviewed')})")
    print(f"  Has GO terms:   {stats['has_go']} ({pct('has_go')})")
    print(f"  Has EC numbers: {stats['has_ec']} ({pct('has_ec')})")
    print(f"  Has function:   {stats['has_function']} ({pct('has_function')})")


# ─── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build protein annotation tables from UniProt preprocessed data."
    )
    add_common_args(parser, DEFAULT_CONFIG)
    args = parser.parse_args()

    config = load_config(args.config)

    # Discover unique (organism_group, taxid) pairs from the genomes CSV
    # (grouped so that e.g. all 3 Alteromonas strains sharing taxid 28108 are one entry)
    seen: dict[int, str] = {}  # taxid → org_group
    with open(GENOMES_CSV, newline="") as f:
        reader = csv.DictReader(
            (line for line in f if not line.strip().startswith("#"))
        )
        for row in reader:
            strain_name = row.get("strain_name", "")
            if args.strains and strain_name not in args.strains:
                continue
            taxon_id_str = (row.get("ncbi_taxon_id") or "").strip()
            if not taxon_id_str:
                continue
            ncbi_taxon_id = int(taxon_id_str)
            if ncbi_taxon_id in seen:
                continue  # already registered this taxid
            data_dir = (row.get("data_dir") or "").strip()
            seen[ncbi_taxon_id] = infer_organism_group(data_dir)

    if not seen:
        msg = f"no strains matched {args.strains}" if args.strains else "no strains found"
        print(f"Error: {msg}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(seen)} unique taxid(s) with config: {args.config}")
    for ncbi_taxon_id, org_group in seen.items():
        process_taxid(org_group, ncbi_taxon_id, config, force=args.force)


if __name__ == "__main__":
    main()
