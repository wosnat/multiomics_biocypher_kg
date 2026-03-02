#!/usr/bin/env python3
"""Resolve gene IDs in paper CSV files to canonical locus tags.

For each paperconfig.yaml, for each 'csv' supplementary table whose
name_col is not already 'locus_tag':
  - Reads the source CSV
  - Resolves name_col values to locus tags via gene_id_mapping.json
  - Writes <stem>_resolved.csv alongside the original, adding
    'locus_tag' and 'resolution_method' columns
  - Skips rows with empty name_col (they become NaN in locus_tag)
  - Reports resolution stats grouped by publication / organism

The omics_adapter probes for <stem>_resolved.csv and uses the
locus_tag column when present, eliminating dangling edges.

Usage:
  uv run python -m multiomics_kg.download.resolve_paper_ids [--force] [--papers NAME ...]
  uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Anjur 2025"
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from multiomics_kg.download.build_gene_id_mapping import (
    load_all_paperconfigs,
    get_organism_for_entry,
)
from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.utils.gene_id_utils import build_id_lookup, get_genome_dir, map_gene_id

MAX_UNRESOLVED_SHOWN = 20


# ─── Path helper ──────────────────────────────────────────────────────────────


def get_resolved_path(filename: str | Path) -> Path:
    """Return the pre-resolved CSV path: <stem>_resolved<ext>."""
    p = Path(filename)
    return p.parent / f"{p.stem}_resolved{p.suffix}"


# ─── Table resolution ─────────────────────────────────────────────────────────


def resolve_table(
    pub_name: str,
    table_key: str,
    table_config: dict,
    force: bool = False,
) -> dict | None:
    """Resolve gene IDs in one CSV supplementary table.

    Returns a result dict, or None if the table should be silently skipped.
    Result dict keys:
      skipped (bool), reason (str, on skip), pub_name, table_key, organism,
      name_col, src, resolved_path, total, resolved, unresolved_rows
    """
    # Only process csv tables with statistical_analyses
    if table_config.get("type", "csv") != "csv":
        return None

    stat_analyses = table_config.get("statistical_analyses") or []
    if not stat_analyses:
        return None

    filename = table_config.get("filename")
    if not filename:
        return None

    src = PROJECT_ROOT / filename
    if not src.exists():
        print(f"  [warn] CSV not found: {src}", file=sys.stderr)
        return None

    sep = table_config.get("sep", ",")
    skip_rows = int(table_config.get("skip_rows", 0) or 0)

    # Determine organism: table-level wins, else first analysis
    organism = (table_config.get("organism") or "").strip().strip('"')
    if not organism:
        for a in stat_analyses:
            org = (a.get("organism") or "").strip().strip('"')
            if org:
                organism = org
                break
    if not organism:
        print(f"  [warn] No organism declared for {pub_name}/{table_key} — skipping", file=sys.stderr)
        return None

    # Determine name_col: first analysis whose name_col is not 'locus_tag'
    name_col = None
    for a in stat_analyses:
        nc = a.get("name_col")
        if nc and nc != "locus_tag":
            name_col = nc
            break
    if not name_col:
        # All analyses already use locus_tag — nothing to do
        return {"skipped": True, "reason": "name_col is already locus_tag"}

    resolved_path = get_resolved_path(src)

    # Skip if _resolved.csv is already up to date
    if not force and resolved_path.exists():
        if resolved_path.stat().st_mtime >= src.stat().st_mtime:
            return {
                "skipped": True,
                "reason": "up to date",
                "resolved_path": str(resolved_path.relative_to(PROJECT_ROOT)),
            }

    # Load gene ID lookup for this organism
    genome_dir = get_genome_dir(organism, str(PROJECT_ROOT))
    if not genome_dir:
        print(
            f"  [warn] No genome dir for organism '{organism}' in {pub_name}/{table_key} — skipping",
            file=sys.stderr,
        )
        return None

    lookup, locus_tags, supp_keys = build_id_lookup(genome_dir)
    if lookup is None:
        print(
            f"  [warn] No gene annotation data for '{organism}' in {pub_name}/{table_key} — skipping",
            file=sys.stderr,
        )
        return None

    # Read source CSV
    try:
        df = pd.read_csv(src, sep=sep, skiprows=skip_rows if skip_rows else None)
    except Exception as e:
        print(f"  [error] Could not read {src}: {e}", file=sys.stderr)
        return None

    if name_col not in df.columns:
        print(
            f"  [warn] name_col '{name_col}' not in {src.name} columns — skipping {pub_name}/{table_key}",
            file=sys.stderr,
        )
        return None

    # Resolve each name_col value
    locus_tag_col: list[str | None] = []
    method_col: list[str] = []
    unresolved_rows: list[dict] = []
    total = 0
    n_resolved = 0

    for row_idx, val in enumerate(df[name_col]):
        if pd.isna(val) or str(val).strip() == "":
            locus_tag_col.append(None)
            method_col.append("empty")
            continue

        raw = str(val).strip().strip("*").strip()
        if not raw:
            locus_tag_col.append(None)
            method_col.append("empty")
            continue

        total += 1
        lt, method = map_gene_id(raw, lookup, locus_tags, supp_keys)
        if lt:
            locus_tag_col.append(lt)
            method_col.append(method or "lookup")
            n_resolved += 1
        else:
            locus_tag_col.append(None)
            method_col.append("unresolved")
            unresolved_rows.append({"row": row_idx, "raw_id": raw})

    # Write _resolved.csv
    df_out = df.copy()
    df_out["locus_tag"] = locus_tag_col
    df_out["resolution_method"] = method_col
    try:
        df_out.to_csv(resolved_path, index=False)
    except Exception as e:
        print(f"  [error] Could not write {resolved_path}: {e}", file=sys.stderr)
        return None

    return {
        "skipped": False,
        "pub_name": pub_name,
        "table_key": table_key,
        "organism": organism,
        "name_col": name_col,
        "src": str(src.relative_to(PROJECT_ROOT)),
        "resolved_path": str(resolved_path.relative_to(PROJECT_ROOT)),
        "total": total,
        "resolved": n_resolved,
        "unresolved_rows": unresolved_rows,
    }


# ─── Report ───────────────────────────────────────────────────────────────────


def print_report(results: list[dict]) -> None:
    """Print resolution stats grouped by publication / organism."""
    # Group by pub_name → organism
    by_pub: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for r in results:
        by_pub[r["pub_name"]][r["organism"]].append(r)

    total_resolved = sum(r["resolved"] for r in results)
    total_total = sum(r["total"] for r in results)
    pct_overall = 100 * total_resolved / max(total_total, 1)

    print()
    print("=" * 70)
    print("Gene ID Resolution Report")
    print(f"Overall: {total_resolved}/{total_total} resolved ({pct_overall:.1f}%)")
    print("=" * 70)

    for pub_name, by_org in sorted(by_pub.items()):
        print(f"\n  Publication: {pub_name}")
        for organism, tables in sorted(by_org.items()):
            print(f"    Organism: {organism}")
            for r in tables:
                pct = 100 * r["resolved"] / max(r["total"], 1)
                print(f"      [{r['table_key']}]  name_col={r['name_col']}")
                print(f"        Resolved:  {r['resolved']}/{r['total']} ({pct:.1f}%)")
                print(f"        Output:    {r['resolved_path']}")
                unresolved = r["unresolved_rows"]
                if unresolved:
                    n = len(unresolved)
                    shown = unresolved[:MAX_UNRESOLVED_SHOWN]
                    print(f"        Unresolved ({n}):")
                    for item in shown:
                        print(f"          row {item['row']}: {item['raw_id']}")
                    if n > MAX_UNRESOLVED_SHOWN:
                        print(f"          ... ({n - MAX_UNRESOLVED_SHOWN} more)")
    print()


# ─── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-resolve even if _resolved.csv is already up to date",
    )
    parser.add_argument(
        "--papers",
        nargs="+",
        metavar="NAME",
        help="Only process papers whose name contains one of these strings (case-insensitive)",
    )
    args = parser.parse_args()

    paperconfigs = load_all_paperconfigs()
    if not paperconfigs:
        print("No paperconfigs found.", file=sys.stderr)
        sys.exit(1)

    results: list[dict] = []

    for pc_path, config in paperconfigs:
        pub = config.get("publication") or {}
        pub_name = pub.get("papername") or pc_path.parent.name
        dir_name = pc_path.parent.name

        # --papers filter: match against papername OR directory name
        if args.papers:
            if not any(
                p.lower() in pub_name.lower() or p.lower() in dir_name.lower()
                for p in args.papers
            ):
                continue

        supp = (
            pub.get("supplementary_materials")
            or config.get("supplementary_materials")
            or {}
        )

        for table_key, table_config in supp.items():
            if not isinstance(table_config, dict):
                continue

            result = resolve_table(pub_name, table_key, table_config, force=args.force)
            if result is None:
                continue  # silently skipped (wrong type / missing file / no organism)
            if result.get("skipped"):
                reason = result.get("reason", "")
                # Only mention non-trivial skips
                if reason not in ("name_col is already locus_tag", "up to date"):
                    print(f"  Skipped {pub_name}/{table_key}: {reason}")
                continue

            results.append(result)

    if not results:
        print("No CSV tables resolved (all up to date or nothing to do).")
        return

    print_report(results)


if __name__ == "__main__":
    main()
