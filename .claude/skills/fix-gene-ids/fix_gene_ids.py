#!/usr/bin/env python
"""Map gene IDs in paper CSVs to locus tags using gene_annotations_merged.json.

For each analysis in a paperconfig, reads the paper CSV and the organism's
gene_annotations_merged.json (plus gene_mapping_supp.csv if available), maps gene IDs
to locus tags, and writes a new CSV with a locus_tag column.

New features:
- Alternative column scanning: for unmapped/mismatched rows, tries other CSV
  columns to find mappable gene IDs
- Import report integration: uses Neo4j import report to identify mismatched
  IDs that create dangling edges
- Patch mode: updates existing _with_locus_tag.csv files in place
- Report generation: writes fix_gene_ids_report.md per paper directory

Usage:
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml"
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --dry-run
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --import-report output/import.report
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --patch --import-report output/import.report
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

# Import shared utilities from the main package
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_PROJECT_ROOT))
from multiomics_kg.utils.gene_id_utils import (
    SKIP_PATTERNS,
    DESCRIPTION_COL_KEYWORDS,
    ID_COL_KEYWORDS,
    MAX_UNIQUE_FOR_ID_COL,
    get_genome_dir,
    build_id_lookup,
    load_supp_lookup,
    map_gene_id,
    is_id_like_column,
    load_import_report,
)


def get_candidate_columns(df, analysis):
    """Find columns in the CSV that might contain alternative gene IDs.

    Excludes name_col, logfc_col, adjusted_p_value_col, and non-ID columns.
    Returns list of column names.
    """
    name_col = analysis.get("name_col", "")
    logfc_col = analysis.get("logfc_col", "")
    adj_pval_col = analysis.get("adjusted_p_value_col", "")

    exclude_cols = {name_col}
    if logfc_col:
        exclude_cols.add(logfc_col)
    if adj_pval_col:
        exclude_cols.add(adj_pval_col)
    # Also exclude locus_tag column if it exists (from prior fix runs)
    exclude_cols.add("locus_tag")

    candidates = []
    for col in df.columns:
        if col in exclude_cols:
            continue
        if is_id_like_column(df[col], col, exclude_cols):
            candidates.append(col)

    return candidates


def process_analysis(paperconfig_path, analysis, supp_material, project_root,
                     dry_run=False, missing_gene_ids=None, patch_mode=False):
    """Process a single analysis: map gene IDs and write new CSV.

    Args:
        missing_gene_ids: set of gene IDs (without prefix) known to be missing
            from Neo4j import report. Used to identify rows that need alt-column recovery.
        patch_mode: if True, read existing _with_locus_tag.csv and only fix
            rows with empty/mismatched locus_tag.

    Returns dict with mapping stats.
    """
    analysis_id = analysis.get("id", "unknown")
    organism = analysis.get("organism", "")
    name_col = analysis.get("name_col", "")
    skip_rows = supp_material.get("skip_rows", analysis.get("skip_rows", 0))
    csv_path = supp_material.get("filename", "")

    result = {
        "analysis_id": analysis_id,
        "organism": organism,
        "csv_path": csv_path,
        "total": 0,
        "mapped": 0,
        "direct": 0,
        "lookup": 0,
        "composite": 0,
        "repadded": 0,
        "supp": 0,
        "skipped": 0,
        "unmapped": [],
        "output_file": None,
        # Alt-column recovery stats
        "alt_recovered": 0,
        "alt_by_column": {},  # col_name -> count
        "alt_details": [],  # list of (row_idx, original_id, alt_col, alt_value, locus_tag)
        "mismatch_recovered": 0,  # subset of alt_recovered that were import report mismatches
        "still_unmapped": [],  # IDs that remain unmapped after alt-column scan
    }

    if not csv_path or not os.path.exists(csv_path):
        print(f"  CSV not found: {csv_path}")
        return result

    if not organism:
        print(f"  No organism specified for analysis {analysis_id}")
        return result

    # Find genome dir
    genome_dir = get_genome_dir(organism, project_root)
    if not genome_dir:
        print(f"  No genome directory found for organism: {organism}")
        return result

    # Build lookup from gene_annotations_merged.json (includes supp fallback)
    lookup, locus_tags, supp_keys = build_id_lookup(genome_dir)
    if lookup is None:
        print(f"  Could not load gene_annotations_merged.json from {genome_dir}")
        return result

    # In patch mode, read the existing _with_locus_tag.csv
    if patch_mode:
        base = Path(csv_path)
        # If csv_path already ends with _with_locus_tag, use it directly
        if "_with_locus_tag" in base.stem:
            patch_file = base
        else:
            patch_file = base.parent / (base.stem + "_with_locus_tag" + base.suffix)
        if not patch_file.exists():
            print(f"  Patch target not found: {patch_file}")
            return result
        try:
            df = pd.read_csv(str(patch_file))
        except Exception as e:
            print(f"  Error reading {patch_file}: {e}")
            return result
        csv_path_for_output = str(patch_file)
    else:
        # Load paper CSV normally
        try:
            if skip_rows:
                df = pd.read_csv(csv_path, skiprows=skip_rows)
            else:
                df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"  Error reading {csv_path}: {e}")
            return result
        csv_path_for_output = csv_path

    if name_col not in df.columns:
        print(f"  Column '{name_col}' not found in {csv_path_for_output}")
        print(f"  Available columns: {list(df.columns)}")
        return result

    # ---- Phase 1: Primary name_col mapping ----
    mapped_locus_tags = []
    result["total"] = len(df)
    needs_alt_scan = []  # indices of rows needing alt-column recovery

    for idx, (_, row) in enumerate(df.iterrows()):
        raw_id = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""

        if not raw_id or raw_id == "nan":
            if patch_mode:
                # In patch mode, empty locus_tag rows need alt-column recovery
                mapped_locus_tags.append("")
                needs_alt_scan.append(idx)
                continue
            mapped_locus_tags.append("")
            continue

        # Skip tRNA/ncRNA/rRNA
        if SKIP_PATTERNS.match(raw_id):
            mapped_locus_tags.append("")
            result["skipped"] += 1
            continue

        # In patch mode, check existing locus_tag first
        if patch_mode and "locus_tag" in df.columns:
            existing_lt = str(row.get("locus_tag", "")).strip()
            if existing_lt and existing_lt != "nan":
                # Check if this locus_tag is in the mismatch set
                if missing_gene_ids and existing_lt in missing_gene_ids:
                    # Existing locus_tag creates a dangling edge — needs recovery
                    mapped_locus_tags.append(existing_lt)  # placeholder, will be overwritten
                    needs_alt_scan.append(idx)
                    continue
                else:
                    # Existing locus_tag is good
                    mapped_locus_tags.append(existing_lt)
                    result["mapped"] += 1
                    result["direct"] += 1
                    continue
            else:
                # Empty locus_tag — needs recovery
                mapped_locus_tags.append("")
                needs_alt_scan.append(idx)
                continue

        locus_tag, method = map_gene_id(raw_id, lookup, locus_tags, supp_keys)

        if locus_tag:
            # Check if this locus_tag is known to be missing from Neo4j
            if missing_gene_ids and locus_tag in missing_gene_ids:
                # Primary mapping produced a dangling ID — mark for alt-column recovery
                mapped_locus_tags.append(locus_tag)  # placeholder
                needs_alt_scan.append(idx)
            else:
                mapped_locus_tags.append(locus_tag)
                result["mapped"] += 1
                if method == "direct":
                    result["direct"] += 1
                elif method == "repadded":
                    result["repadded"] += 1
                elif method == "supp":
                    result["supp"] += 1
                elif method in ("lookup", "composite_lookup", "composite_direct"):
                    if "composite" in method:
                        result["composite"] += 1
                    else:
                        result["lookup"] += 1
        else:
            mapped_locus_tags.append("")
            needs_alt_scan.append(idx)

    # ---- Phase 2: Alternative column scanning for unmapped/mismatched rows ----
    candidate_cols = get_candidate_columns(df, analysis) if needs_alt_scan else []

    if needs_alt_scan and candidate_cols:
        print(f"    Alt-column candidates: {candidate_cols}")
        df_rows = df.reset_index(drop=True)

        for idx in needs_alt_scan:
            row = df_rows.iloc[idx]
            original_id = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""
            was_mismatch = (
                mapped_locus_tags[idx] != "" and
                missing_gene_ids and
                mapped_locus_tags[idx] in missing_gene_ids
            )

            recovered = False
            for col in candidate_cols:
                alt_val = str(row[col]).strip().strip('*').strip() if pd.notna(row[col]) else ""
                if not alt_val or alt_val == "nan":
                    continue

                alt_lt, alt_method = map_gene_id(alt_val, lookup, locus_tags, supp_keys)
                if alt_lt:
                    # Check that the recovered locus_tag is NOT itself missing
                    if missing_gene_ids and alt_lt in missing_gene_ids:
                        continue  # This alt also creates a dangling edge, skip

                    mapped_locus_tags[idx] = alt_lt
                    result["alt_recovered"] += 1
                    result["mapped"] += 1
                    result["alt_by_column"][col] = result["alt_by_column"].get(col, 0) + 1
                    result["alt_details"].append((idx, original_id, col, alt_val, alt_lt))
                    if was_mismatch:
                        result["mismatch_recovered"] += 1
                    recovered = True
                    break

            if not recovered:
                display_id = original_id if original_id else f"(row {idx}, empty name_col)"
                result["unmapped"].append(display_id)
                result["still_unmapped"].append({
                    "original_id": display_id,
                    "alt_values": {
                        col: str(row[col]).strip() if pd.notna(row[col]) else ""
                        for col in candidate_cols
                    },
                })
    elif needs_alt_scan:
        # No candidate columns — all go to unmapped
        df_rows = df.reset_index(drop=True)
        for idx in needs_alt_scan:
            row = df_rows.iloc[idx]
            original_id = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""
            display_id = original_id if original_id else f"(row {idx}, empty name_col)"
            result["unmapped"].append(display_id)
            result["still_unmapped"].append({
                "original_id": display_id,
                "alt_values": {},
            })

    # ---- Phase 3: Write output ----
    df["locus_tag"] = mapped_locus_tags

    if not dry_run:
        if patch_mode:
            output_file = csv_path_for_output
            df.to_csv(output_file, index=False)
            result["output_file"] = output_file
        else:
            base = Path(csv_path)
            output_file = base.parent / (base.stem + "_with_locus_tag" + base.suffix)
            df.to_csv(output_file, index=False)
            result["output_file"] = str(output_file)

    return result


def process_paperconfig(paperconfig_path, project_root, dry_run=False,
                        missing_gene_ids=None, patch_mode=False):
    """Process all analyses in a paperconfig file."""
    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    pub = config.get("publication", {})
    papername = pub.get("papername", "Unknown")
    supp_materials = pub.get("supplementary_materials", {})

    print(f"Processing: {papername}")
    print(f"  Config: {paperconfig_path}")
    print()

    results = []

    # Track which CSV files we've already processed (multiple analyses may share the same CSV)
    processed_csvs = {}

    for supp_key, supp_material in supp_materials.items():
        if supp_material.get("type") != "csv":
            continue

        csv_path = supp_material.get("filename", "")
        analyses = supp_material.get("statistical_analyses", [])

        for analysis in analyses:
            analysis_id = analysis.get("id", "unknown")

            # If we already processed this CSV for a different analysis with the same
            # organism and name_col, skip (the locus_tag column is already added)
            cache_key = (csv_path, analysis.get("organism", ""), analysis.get("name_col", ""))
            if cache_key in processed_csvs:
                print(f"  [{analysis_id}] Same CSV already processed, reusing result")
                results.append(processed_csvs[cache_key])
                continue

            print(f"  [{analysis_id}] {os.path.basename(csv_path)}")
            print(f"    Organism: {analysis.get('organism', 'N/A')}")
            print(f"    name_col: {analysis.get('name_col', 'N/A')}")

            result = process_analysis(
                paperconfig_path, analysis, supp_material, project_root,
                dry_run, missing_gene_ids, patch_mode
            )
            results.append(result)
            processed_csvs[cache_key] = result

            # Print stats
            total_valid = result["total"] - result["skipped"]
            if total_valid > 0:
                print(f"    Mapped: {result['mapped']}/{total_valid} "
                      f"({result['direct']} direct, {result['lookup']} via gene_names, "
                      f"{result['repadded']} via repadding, "
                      f"{result['composite']} via composite, "
                      f"{result['supp']} via supp mapping)")
                if result["alt_recovered"]:
                    print(f"    Alt-column recovered: {result['alt_recovered']} "
                          f"({result['mismatch_recovered']} were import-report mismatches)")
                    for col, cnt in sorted(result["alt_by_column"].items()):
                        print(f"      via '{col}': {cnt}")
                if result["skipped"]:
                    print(f"    Skipped: {result['skipped']} (tRNA/ncRNA/rRNA)")
                if result["unmapped"]:
                    unmapped_display = result["unmapped"][:20]
                    print(f"    Unmapped: {len(result['unmapped'])} — {', '.join(unmapped_display)}"
                          + (f" ... +{len(result['unmapped']) - 20} more" if len(result['unmapped']) > 20 else ""))
                if result["output_file"]:
                    print(f"    Wrote: {result['output_file']}")
                elif dry_run:
                    print(f"    (dry run — no file written)")
            print()

    return results


def write_report(paperconfig_path, papername, results, project_root, dry_run=False):
    """Write fix_gene_ids_report.md in the paper's directory."""
    paper_dir = Path(paperconfig_path).parent
    report_path = paper_dir / "fix_gene_ids_report.md"

    lines = []
    lines.append(f"# Fix Gene IDs Report: {papername}")
    lines.append(f"")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"Paperconfig: `{paperconfig_path}`")
    lines.append(f"")

    total_mapped = 0
    total_unmapped = 0
    total_alt_recovered = 0
    total_mismatch_recovered = 0

    for r in results:
        lines.append(f"## Analysis: `{r['analysis_id']}`")
        lines.append(f"")
        lines.append(f"- **CSV**: `{os.path.basename(r['csv_path'])}`")
        lines.append(f"- **Organism**: {r['organism']}")
        lines.append(f"- **Total rows**: {r['total']}")
        lines.append(f"- **Skipped (RNA)**: {r['skipped']}")
        lines.append(f"")

        total_valid = r["total"] - r["skipped"]
        lines.append(f"### Primary mapping (via name_col)")
        lines.append(f"")
        primary_mapped = r["mapped"] - r["alt_recovered"]
        lines.append(f"- Mapped: {primary_mapped}/{total_valid}")
        lines.append(f"  - Direct: {r['direct']}")
        lines.append(f"  - Via gene_names: {r['lookup']}")
        lines.append(f"  - Via repadding: {r['repadded']}")
        lines.append(f"  - Via composite: {r['composite']}")
        lines.append(f"  - Via supp mapping: {r['supp']}")
        lines.append(f"")

        if r["alt_recovered"] > 0 or r["alt_by_column"]:
            lines.append(f"### Alternative column recovery")
            lines.append(f"")
            lines.append(f"- Recovered: **{r['alt_recovered']}** rows")
            if r["mismatch_recovered"]:
                lines.append(f"- Of which **{r['mismatch_recovered']}** were import-report mismatches (dangling edges)")
            for col, cnt in sorted(r["alt_by_column"].items()):
                lines.append(f"  - Via column `{col}`: {cnt}")
            lines.append(f"")

            if r["alt_details"]:
                lines.append(f"#### Recovery details")
                lines.append(f"")
                lines.append(f"| Row | Original ID | Alt Column | Alt Value | Mapped Locus Tag |")
                lines.append(f"|-----|-------------|------------|-----------|-----------------|")
                for idx, orig, col, alt_val, lt in r["alt_details"][:50]:
                    lines.append(f"| {idx} | {orig} | {col} | {alt_val} | {lt} |")
                if len(r["alt_details"]) > 50:
                    lines.append(f"| ... | +{len(r['alt_details']) - 50} more | | | |")
                lines.append(f"")

        if r["still_unmapped"]:
            lines.append(f"### Still unmapped ({len(r['still_unmapped'])} rows)")
            lines.append(f"")
            lines.append(f"| Original ID | Alternative column values |")
            lines.append(f"|-------------|--------------------------|")
            for info in r["still_unmapped"][:30]:
                alt_str = "; ".join(f"{k}={v}" for k, v in info["alt_values"].items() if v)
                if not alt_str:
                    alt_str = "(no alt columns)"
                lines.append(f"| {info['original_id']} | {alt_str} |")
            if len(r["still_unmapped"]) > 30:
                lines.append(f"| ... | +{len(r['still_unmapped']) - 30} more |")
            lines.append(f"")

        if r["output_file"]:
            lines.append(f"**Output**: `{os.path.basename(r['output_file'])}`")
            lines.append(f"")

        total_mapped += r["mapped"]
        total_unmapped += len(r["unmapped"])
        total_alt_recovered += r["alt_recovered"]
        total_mismatch_recovered += r["mismatch_recovered"]

    # Summary
    lines.append(f"## Summary")
    lines.append(f"")
    lines.append(f"- Total mapped: {total_mapped}")
    lines.append(f"- Total unmapped: {total_unmapped}")
    lines.append(f"- Alt-column recovered: {total_alt_recovered}")
    if total_mismatch_recovered:
        lines.append(f"- Import-report mismatches recovered: {total_mismatch_recovered}")

    report_text = "\n".join(lines) + "\n"

    if not dry_run:
        with open(report_path, "w") as f:
            f.write(report_text)
        print(f"Report written to: {report_path}")
    else:
        print(f"(dry run — would write report to: {report_path})")

    return report_text


def main():
    parser = argparse.ArgumentParser(
        description="Map gene IDs in paper CSVs to locus tags using gene_annotations_merged.json"
    )
    parser.add_argument(
        "--paperconfig", required=True,
        help="Path to paperconfig.yaml file"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show mapping stats without writing files"
    )
    parser.add_argument(
        "--import-report",
        help="Path to import.report file (default: tries output/import.report, then Docker)",
    )
    parser.add_argument(
        "--no-import-report",
        action="store_true",
        help="Skip loading import report (disable mismatch-based alt-column recovery)",
    )
    parser.add_argument(
        "--patch",
        action="store_true",
        help="Patch existing _with_locus_tag.csv files (only fix empty/mismatched rows)",
    )
    args = parser.parse_args()

    # Determine project root (directory containing this script's repo)
    project_root = Path(__file__).resolve().parents[3]

    if not os.path.exists(args.paperconfig):
        print(f"Error: paperconfig not found: {args.paperconfig}")
        return 1

    # Load import report for mismatch detection
    missing_gene_ids = None
    if not args.no_import_report:
        missing_gene_ids = load_import_report(args.import_report, return_by_source=False)
        if missing_gene_ids is None:
            print("No import report available (use --import-report or run Docker)", file=sys.stderr)
        else:
            print(f"Import report: {len(missing_gene_ids)} unique missing gene IDs", file=sys.stderr)

    # Read papername for report
    with open(args.paperconfig) as f:
        config = yaml.safe_load(f)
    papername = config.get("publication", {}).get("papername", "Unknown")

    results = process_paperconfig(
        args.paperconfig, project_root, args.dry_run,
        missing_gene_ids, args.patch
    )

    # Write per-paper report
    write_report(args.paperconfig, papername, results, project_root, args.dry_run)

    # Summary
    total_mapped = sum(r["mapped"] for r in results)
    total_unmapped = sum(len(r["unmapped"]) for r in results)
    total_skipped = sum(r["skipped"] for r in results)
    total = sum(r["total"] for r in results)
    total_alt = sum(r["alt_recovered"] for r in results)

    print("=" * 50)
    print(f"Summary: {total_mapped} mapped, {total_unmapped} unmapped, "
          f"{total_skipped} skipped out of {total} total rows")
    if total_alt:
        print(f"  Alt-column recovered: {total_alt}")

    return 0


if __name__ == "__main__":
    exit(main())
