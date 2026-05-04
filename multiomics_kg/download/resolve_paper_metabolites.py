#!/usr/bin/env python3
"""Step 7 — resolve metabolite names in paper CSV files.

For each metabolite_assays_table entry across all paperconfigs, opens the
source CSV, looks up each row's metabolite (by id_col first, then name_col
via metabolite_id_mapping.json), and writes <stem>_resolved.csv with two
added columns: metabolite_id + resolution_method.

Step 7 does NOT load the MNX resolver — it reads only the lightweight
metabolite_id_mapping.json that step 6 produced. This mirrors gene-side
step 4 (resolve_paper_ids.py).

Usage:
    uv run python -m multiomics_kg.download.resolve_paper_metabolites [--force] [--papers NAME ...]
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.utils.paperconfig_utils import (
    iter_metabolite_assays_tables,
    load_all_paperconfigs,
    get_paper_name,
)

log = logging.getLogger(__name__)

MAPPING_FILE = PROJECT_ROOT / "cache" / "data" / "metabolomics" / "metabolite_id_mapping.json"


def _resolve_row(
    row: pd.Series,
    name_col: str,
    id_col: str,
    id_type: str,
    name_lookup: dict,
) -> tuple[str | None, str]:
    """Return (primary_id_or_None, resolution_method)."""
    # 1. id_col direct
    if id_col and id_type:
        raw = (row.get(id_col) or "").strip()
        if raw:
            primary = raw if ":" in raw else f"{id_type}:{raw}"
            method = (
                "kegg_direct" if id_type == "kegg.compound"
                else "chebi_direct" if id_type == "chebi"
                else "mnx_direct" if id_type == "mnx"
                else "id_direct"
            )
            return primary, method

    # 2. name_col → mapping
    name = (row.get(name_col) or "").strip()
    if not name:
        return None, "unresolved"
    hits = name_lookup.get(name) or []
    if len(hits) == 1:
        return hits[0], "name_match"
    if len(hits) > 1:
        # ambiguous; v1 picks first deterministically and tags
        return sorted(hits)[0], "ambiguous_multi_id"
    return None, "unresolved"


def resolve_paper_metabolites_for_entry(
    entry: dict, mapping: dict, paperconfig_path: Path | None = None
) -> tuple[Path, dict]:
    """Resolve one metabolite_assays_table entry's CSV.

    Returns (resolved_csv_path, report_dict).
    """
    src_str = entry["filename"]
    src = Path(src_str) if Path(src_str).is_absolute() else (PROJECT_ROOT / src_str)
    if not src.exists():
        raise FileNotFoundError(f"Source CSV missing: {src}")

    name_col = entry["name_col"]
    id_col = entry.get("id_col") or ""
    id_type = entry.get("id_type") or ""
    name_lookup = mapping.get("name_lookup") or {}
    skip_rows = int(entry.get("skip_rows", 0) or 0)

    df = pd.read_csv(src, dtype=str, keep_default_na=False, skiprows=skip_rows or None)

    resolved_ids: list[str] = []
    methods: list[str] = []
    n_resolved = 0
    n_unresolved = 0

    for _, row in df.iterrows():
        primary, method = _resolve_row(row, name_col, id_col, id_type, name_lookup)
        resolved_ids.append(primary or "")
        methods.append(method)
        if primary:
            n_resolved += 1
        else:
            n_unresolved += 1

    df["metabolite_id"] = resolved_ids
    df["resolution_method"] = methods

    out = src.with_name(src.stem + "_resolved.csv")
    df.to_csv(out, index=False)

    method_counts: dict[str, int] = {}
    for m in methods:
        method_counts[m] = method_counts.get(m, 0) + 1

    report = {
        "source_csv": str(src),
        "resolved_csv": str(out),
        "total_rows": int(len(df)),
        "resolved": int(n_resolved),
        "unresolved": int(n_unresolved),
        "match_rate": (n_resolved / len(df)) if len(df) else 0.0,
        "method_counts": method_counts,
    }
    report_path = src.with_name(src.stem + "_resolution_report.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    return out, report


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Re-resolve even if <stem>_resolved.csv exists.")
    parser.add_argument("--papers", nargs="*", default=None,
                        help="Restrict to specific paper names (papername).")
    args = parser.parse_args(argv)

    if not MAPPING_FILE.exists():
        log.error(f"{MAPPING_FILE} not found — run step 6 first.")
        return 1
    mapping = json.loads(MAPPING_FILE.read_text())

    paperconfigs = load_all_paperconfigs()
    total_entries = 0
    for pc_path, cfg in paperconfigs:
        paper = get_paper_name(cfg) or "<unknown>"
        if args.papers and paper not in args.papers:
            continue
        for entry_key, entry in iter_metabolite_assays_tables(cfg):
            src_str = entry.get("filename") or ""
            if not src_str:
                log.warning(f"[step7] {paper}/{entry_key}: missing filename")
                continue
            src = Path(src_str) if Path(src_str).is_absolute() else (PROJECT_ROOT / src_str)
            if not src.exists():
                log.warning(f"[step7] {paper}/{entry_key}: source not found {src}")
                continue
            out = src.with_name(src.stem + "_resolved.csv")
            if out.exists() and not args.force:
                log.info(f"[step7] {paper}/{entry_key}: skipping (exists). Use --force.")
                continue
            try:
                _, report = resolve_paper_metabolites_for_entry(entry, mapping, pc_path)
            except Exception as e:
                log.error(f"[step7] {paper}/{entry_key}: failed: {e}")
                continue
            log.info(
                f"[step7] {paper}/{entry_key}: "
                f"{report['resolved']}/{report['total_rows']} resolved "
                f"({report['match_rate']:.1%})"
            )
            total_entries += 1

    log.info(f"[step7] Done. Processed {total_entries} entries.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
