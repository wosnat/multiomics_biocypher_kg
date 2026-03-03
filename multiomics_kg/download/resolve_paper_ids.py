#!/usr/bin/env python3
"""Resolve gene IDs in paper CSV files to canonical locus tags.

For each paperconfig.yaml, for each 'csv' supplementary table whose
name_col is not already 'locus_tag':
  - Reads the source CSV
  - Resolves each row to a locus_tag via gene_id_mapping.json (v2 or v1)
    using a three-tier multi-column strategy:
    1. specific_lookup (Tier 1 gene-unique IDs), with list expansion
    2. Heuristics (zero-padding, asterisk stripping) → specific_lookup
    3. multi_lookup (Tier 2+3), only when singleton for this organism
  - Falls back to other columns in the row (id_columns from paperconfig)
    before giving up
  - Writes <stem>_resolved.csv alongside the original, adding
    'locus_tag' and 'resolution_method' columns
  - Writes <stem>_resolved_report.txt alongside for diagnostics
  - Never silently skips: every row gets an explicit resolution_method

The omics_adapter probes for <stem>_resolved.csv and uses the
locus_tag column when present, eliminating dangling edges.

Usage:
  uv run python -m multiomics_kg.download.resolve_paper_ids [--force] [--papers NAME ...]
  uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Anjur 2025"
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from multiomics_kg.download.build_gene_id_mapping import (
    load_all_paperconfigs,
    get_organism_for_entry,
)
from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.utils.gene_id_utils import (
    get_genome_dir,
    load_mapping_v2,
    resolve_row,
    MappingData,
)

MAX_UNRESOLVED_SHOWN = 20


# ─── Path helpers ──────────────────────────────────────────────────────────────


def get_resolved_path(filename: str | Path) -> Path:
    """Return the pre-resolved CSV path: <stem>_resolved<ext>."""
    p = Path(filename)
    return p.parent / f"{p.stem}_resolved{p.suffix}"


def get_report_path(resolved_path: Path) -> Path:
    """Return the diagnostic report path alongside a resolved CSV."""
    return resolved_path.parent / f"{resolved_path.stem}_report.txt"


# ─── Table resolution ─────────────────────────────────────────────────────────


def resolve_table(
    pub_name: str,
    table_key: str,
    table_config: dict,
    force: bool = False,
) -> dict | None:
    """Resolve gene IDs in one CSV supplementary table.

    Returns a result dict, or None if the table should be silently skipped.
    """
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

    # Determine organism
    organism = (table_config.get("organism") or "").strip().strip('"')
    if not organism:
        for a in stat_analyses:
            org = (a.get("organism") or "").strip().strip('"')
            if org:
                organism = org
                break
    if not organism:
        print(f"  [warn] No organism for {pub_name}/{table_key} — skipping", file=sys.stderr)
        return None

    # Determine name_col (first analysis whose name_col is not 'locus_tag')
    name_col = None
    for a in stat_analyses:
        nc = a.get("name_col")
        if nc and nc != "locus_tag":
            name_col = nc
            break
    if not name_col:
        return {"skipped": True, "reason": "name_col is already locus_tag"}

    resolved_path = get_resolved_path(src)

    if not force and resolved_path.exists():
        if resolved_path.stat().st_mtime >= src.stat().st_mtime:
            return {
                "skipped": True,
                "reason": "up to date",
                "resolved_path": str(resolved_path.relative_to(PROJECT_ROOT)),
            }

    # Load gene ID mapping for this organism
    genome_dir = get_genome_dir(organism, str(PROJECT_ROOT))
    if not genome_dir:
        print(
            f"  [warn] No genome dir for organism '{organism}' in {pub_name}/{table_key} — skipping",
            file=sys.stderr,
        )
        return None

    mapping_data = load_mapping_v2(genome_dir)
    if mapping_data is None:
        # Fallback: try old build_id_lookup for organisms without v2 mapping
        from multiomics_kg.utils.gene_id_utils import build_id_lookup, map_gene_id
        lookup, locus_tags, supp_keys = build_id_lookup(genome_dir)
        if lookup is None:
            print(
                f"  [warn] No gene annotation data for '{organism}' in {pub_name}/{table_key} — skipping",
                file=sys.stderr,
            )
            return None
        # Wrap legacy lookup in a MappingData for uniform handling
        mapping_data = MappingData(
            specific_lookup=lookup,
            multi_lookup={},
            conflicts={},
            locus_tags=locus_tags if locus_tags else set(),
            version=0,
        )

    # Collect id_columns from the table config
    id_columns: list[dict] = table_config.get("id_columns") or []

    # Read source CSV
    try:
        df = pd.read_csv(src, sep=sep, skiprows=skip_rows if skip_rows else None)
    except Exception as e:
        print(f"  [error] Could not read {src}: {e}", file=sys.stderr)
        return None

    if name_col not in df.columns:
        print(
            f"  [warn] name_col '{name_col}' not in {src.name} — skipping {pub_name}/{table_key}",
            file=sys.stderr,
        )
        return None

    # Resolve each row
    locus_tag_col: list[str | None] = []
    method_col: list[str] = []
    method_counts: Counter = Counter()
    unresolved_rows: list[dict] = []
    total = 0
    n_resolved = 0

    for row_idx, row in df.iterrows():
        name_val = row.get(name_col)
        if pd.isna(name_val) or str(name_val).strip() == "":
            locus_tag_col.append(None)
            method_col.append("empty")
            method_counts["empty"] += 1
            continue

        total += 1
        lt, method = resolve_row(row.to_dict(), name_col, id_columns, mapping_data)

        locus_tag_col.append(lt)
        method_col.append(method)
        method_counts[method] += 1

        if lt:
            n_resolved += 1
        else:
            unresolved_rows.append({
                "row": row_idx,
                "raw_id": str(name_val).strip(),
                "method": method,
                "id_col_vals": {
                    c.get("column", ""): str(row.get(c.get("column", ""), ""))
                    for c in id_columns
                    if c.get("column", "") in df.columns
                },
            })

    # Write _resolved.csv
    df_out = df.copy()
    df_out["locus_tag"] = locus_tag_col
    df_out["resolution_method"] = method_col
    try:
        df_out.to_csv(resolved_path, index=False)
    except Exception as e:
        print(f"  [error] Could not write {resolved_path}: {e}", file=sys.stderr)
        return None

    # Write diagnostic report
    _write_report(
        pub_name, table_key, name_col, id_columns, src, resolved_path,
        total, n_resolved, method_counts, unresolved_rows,
    )

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
        "method_counts": dict(method_counts),
        "unresolved_rows": unresolved_rows,
    }


def _write_report(
    pub_name: str,
    table_key: str,
    name_col: str,
    id_columns: list[dict],
    src: Path,
    resolved_path: Path,
    total: int,
    n_resolved: int,
    method_counts: Counter,
    unresolved_rows: list[dict],
) -> None:
    """Write per-table diagnostic report to <resolved_path>_report.txt."""
    report_path = get_report_path(resolved_path)
    pct = 100 * n_resolved / max(total, 1)
    lines = [
        f"Resolution report: {pub_name} / {table_key}",
        f"Source: {src.name}",
        f"name_col: {name_col}",
        f"id_columns: {[c.get('column') for c in id_columns]}",
        f"",
        f"Total rows with IDs: {total}",
        f"Resolved: {n_resolved} ({pct:.1f}%)",
        f"Unresolved: {total - n_resolved}",
        f"",
        f"Resolution method breakdown:",
    ]
    for method, count in sorted(method_counts.items()):
        lines.append(f"  {method}: {count}")

    if unresolved_rows:
        lines.append(f"")
        lines.append(f"Unresolved rows ({len(unresolved_rows)}):")
        shown = unresolved_rows[:MAX_UNRESOLVED_SHOWN]
        for item in shown:
            line = f"  row {item['row']}: raw_id={item['raw_id']!r}  reason={item['method']}"
            if item["id_col_vals"]:
                col_strs = ", ".join(f"{k}={v!r}" for k, v in item["id_col_vals"].items() if v and v != "nan")
                if col_strs:
                    line += f"  [{col_strs}]"
            lines.append(line)
        if len(unresolved_rows) > MAX_UNRESOLVED_SHOWN:
            lines.append(f"  ... ({len(unresolved_rows) - MAX_UNRESOLVED_SHOWN} more)")

    # Check if any column behaved unexpectedly
    multi_cols = {
        m.split(":", 1)[1]: count
        for m, count in method_counts.items()
        if m.startswith("ambiguous") or m.startswith("multi:")
    }
    if multi_cols:
        ambig_total = sum(
            c for m, c in method_counts.items() if m.startswith("ambiguous")
        )
        if ambig_total > 0:
            lines.append(f"")
            lines.append(
                f"NOTE: {ambig_total} rows were ambiguous (multi_lookup had >1 match). "
                f"Consider adding id_columns with more specific IDs, or check if the "
                f"name_col should be Tier 2 (protein-level) instead of Tier 1."
            )

    try:
        with open(report_path, "w") as f:
            f.write("\n".join(lines) + "\n")
    except Exception:
        pass  # Non-critical


# ─── Console report ────────────────────────────────────────────────────────────


def print_report(results: list[dict]) -> None:
    """Print resolution stats grouped by publication / organism."""
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

                # Show method breakdown for non-trivial cases
                methods = r.get("method_counts", {})
                non_trivial = {
                    m: c for m, c in methods.items()
                    if m not in ("direct", "empty") and c > 0
                }
                if non_trivial:
                    print(f"        Methods:   " + ", ".join(f"{m}={c}" for m, c in sorted(non_trivial.items())))

                unresolved = r["unresolved_rows"]
                if unresolved:
                    n = len(unresolved)
                    shown = unresolved[:MAX_UNRESOLVED_SHOWN]
                    print(f"        Unresolved ({n}) [{unresolved[0]['method'] if unresolved else '?'}]:")
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
                continue
            if result.get("skipped"):
                reason = result.get("reason", "")
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
