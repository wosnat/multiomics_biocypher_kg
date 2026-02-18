#!/usr/bin/env python
"""Map gene IDs in paper CSVs to locus tags using gene_mapping.csv.

For each analysis in a paperconfig, reads the paper CSV and the organism's
gene_mapping.csv, maps gene IDs to locus tags, and writes a new CSV with
a locus_tag column.

Usage:
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml"
    uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --dry-run
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd
import yaml


# Organism name -> genome directory mapping (same as check-gene-ids)
ORGANISM_TO_GENOME_DIR = {
    "prochlorococcus med4": "cache/data/Prochlorococcus/genomes/MED4",
    "prochlorococcus as9601": "cache/data/Prochlorococcus/genomes/AS9601",
    "prochlorococcus mit9301": "cache/data/Prochlorococcus/genomes/MIT9301",
    "prochlorococcus mit9312": "cache/data/Prochlorococcus/genomes/MIT9312",
    "prochlorococcus mit9313": "cache/data/Prochlorococcus/genomes/MIT9313",
    "prochlorococcus natl1a": "cache/data/Prochlorococcus/genomes/NATL1A",
    "prochlorococcus natl2a": "cache/data/Prochlorococcus/genomes/NATL2A",
    "prochlorococcus rsp50": "cache/data/Prochlorococcus/genomes/RSP50",
    "synechococcus cc9311": "cache/data/Synechococcus/genomes/CC9311",
    "synechococcus wh8102": "cache/data/Synechococcus/genomes/WH8102",
    "alteromonas macleodii mit1002": "cache/data/Alteromonas/genomes/MIT1002",
    "alteromonas mit1002": "cache/data/Alteromonas/genomes/MIT1002",
    "alteromonas macleodii ez55": "cache/data/Alteromonas/genomes/EZ55",
    "alteromonas ez55": "cache/data/Alteromonas/genomes/EZ55",
}

# Columns in gene_mapping.csv to build lookups from
MAPPING_COLS = [
    "locus_tag", "gene_names", "gene", "locus_tag_ncbi",
    "gene_names_cyanorak", "locus_tag_cyanoak", "old_locus_tags",
    "protein_id",
]

# Patterns to skip (tRNA, ncRNA, rRNA — intentionally not in gene nodes)
SKIP_PATTERNS = re.compile(r'^(tRNA|ncRNA|rRNA|Yfr\d|tmRNA)', re.IGNORECASE)


def get_genome_dir(organism, project_root):
    """Get the genome directory path for an organism."""
    key = organism.strip().lower()
    rel = ORGANISM_TO_GENOME_DIR.get(key)
    if rel:
        full = Path(project_root) / rel
        if full.exists():
            return str(full)
    return None


def build_id_lookup(genome_dir):
    """Build a unified lookup dict: raw_id -> locus_tag from gene_mapping.csv.

    Handles space-separated gene_names and comma-separated old_locus_tags.
    Returns (lookup_dict, locus_tag_set) or (None, None) if file not found.
    """
    mapping_file = Path(genome_dir) / "gene_mapping.csv"
    if not mapping_file.exists():
        return None, None

    try:
        df = pd.read_csv(mapping_file)
    except Exception as e:
        print(f"  Error reading {mapping_file}: {e}")
        return None, None

    lookup = {}
    locus_tags = set()

    for _, row in df.iterrows():
        locus = str(row["locus_tag"]).strip() if pd.notna(row.get("locus_tag")) else ""
        if not locus or locus == "nan":
            continue
        locus_tags.add(locus)

        for col in MAPPING_COLS:
            if col not in df.columns:
                continue
            val = str(row[col]).strip() if pd.notna(row[col]) else ""
            if not val or val == "nan":
                continue

            if col in ("gene_names", "gene_names_cyanorak"):
                # Space-separated list of gene names
                for name in val.split():
                    name = name.strip()
                    if name and name != "nan":
                        lookup[name] = locus
            elif col == "old_locus_tags":
                # Comma-separated list of old locus tags
                for tag in val.split(","):
                    tag = tag.strip()
                    if tag and tag != "nan":
                        lookup[tag] = locus
            else:
                lookup[val] = locus

    return lookup, locus_tags


def map_gene_id(raw_id, lookup, locus_tags):
    """Map a single gene ID to a locus_tag.

    Returns (locus_tag, method) or (None, None) if unmapped.
    method is one of: 'direct', 'lookup', 'composite_lookup'
    """
    raw_id = raw_id.strip()

    # Already a valid locus_tag
    if raw_id in locus_tags:
        return raw_id, "direct"

    # Direct lookup
    if raw_id in lookup:
        return lookup[raw_id], "lookup"

    # Try zero-padding normalization (e.g., MIT1002_0001 -> MIT1002_00001)
    m = re.match(r'^(.+)_(\d+)$', raw_id)
    if m:
        prefix, num = m.group(1), m.group(2)
        for pad in range(len(num) + 1, len(num) + 3):
            padded = f"{prefix}_{num.zfill(pad)}"
            if padded in locus_tags:
                return padded, "repadded"
            if padded in lookup:
                return lookup[padded], "repadded"

    # Try splitting by comma (composite gene names like "rps13,rpsM")
    if "," in raw_id:
        parts = [p.strip() for p in raw_id.split(",")]
        for part in parts:
            if part in locus_tags:
                return part, "composite_direct"
            if part in lookup:
                return lookup[part], "composite_lookup"

    return None, None


def process_analysis(paperconfig_path, analysis, supp_material, project_root, dry_run=False):
    """Process a single analysis: map gene IDs and write new CSV.

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
        "skipped": 0,
        "unmapped": [],
        "output_file": None,
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

    # Build lookup
    lookup, locus_tags = build_id_lookup(genome_dir)
    if lookup is None:
        print(f"  Could not load gene_mapping.csv from {genome_dir}")
        return result

    # Load paper CSV
    try:
        if skip_rows:
            df = pd.read_csv(csv_path, skiprows=skip_rows)
        else:
            df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"  Error reading {csv_path}: {e}")
        return result

    if name_col not in df.columns:
        print(f"  Column '{name_col}' not found in {csv_path}")
        print(f"  Available columns: {list(df.columns)}")
        return result

    # Map each gene ID
    mapped_locus_tags = []
    result["total"] = len(df)

    for _, row in df.iterrows():
        raw_id = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""

        if not raw_id or raw_id == "nan":
            mapped_locus_tags.append("")
            continue

        # Skip tRNA/ncRNA/rRNA
        if SKIP_PATTERNS.match(raw_id):
            mapped_locus_tags.append("")
            result["skipped"] += 1
            continue

        locus_tag, method = map_gene_id(raw_id, lookup, locus_tags)

        if locus_tag:
            mapped_locus_tags.append(locus_tag)
            result["mapped"] += 1
            if method == "direct":
                result["direct"] += 1
            elif method == "repadded":
                result["repadded"] += 1
            elif method in ("lookup", "composite_lookup", "composite_direct"):
                if "composite" in method:
                    result["composite"] += 1
                else:
                    result["lookup"] += 1
        else:
            mapped_locus_tags.append("")
            result["unmapped"].append(raw_id)

    # Add locus_tag column
    df["locus_tag"] = mapped_locus_tags

    # Write output
    if not dry_run:
        base = Path(csv_path)
        output_file = base.parent / (base.stem + "_with_locus_tag" + base.suffix)
        df.to_csv(output_file, index=False)
        result["output_file"] = str(output_file)

    return result


def process_paperconfig(paperconfig_path, project_root, dry_run=False):
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
                paperconfig_path, analysis, supp_material, project_root, dry_run
            )
            results.append(result)
            processed_csvs[cache_key] = result

            # Print stats
            total_valid = result["total"] - result["skipped"]
            if total_valid > 0:
                print(f"    Mapped: {result['mapped']}/{total_valid} "
                      f"({result['direct']} direct, {result['lookup']} via gene_names, "
                      f"{result['repadded']} via repadding, "
                      f"{result['composite']} via composite)")
                if result["skipped"]:
                    print(f"    Skipped: {result['skipped']} (tRNA/ncRNA/rRNA)")
                if result["unmapped"]:
                    print(f"    Unmapped: {len(result['unmapped'])} — {', '.join(result['unmapped'])}")
                if result["output_file"]:
                    print(f"    Wrote: {result['output_file']}")
                elif dry_run:
                    print(f"    (dry run — no file written)")
            print()

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Map gene IDs in paper CSVs to locus tags using gene_mapping.csv"
    )
    parser.add_argument(
        "--paperconfig", required=True,
        help="Path to paperconfig.yaml file"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show mapping stats without writing files"
    )
    args = parser.parse_args()

    # Determine project root (directory containing this script's repo)
    project_root = Path(__file__).resolve().parents[3]

    if not os.path.exists(args.paperconfig):
        print(f"Error: paperconfig not found: {args.paperconfig}")
        return 1

    results = process_paperconfig(args.paperconfig, project_root, args.dry_run)

    # Summary
    total_mapped = sum(r["mapped"] for r in results)
    total_unmapped = sum(len(r["unmapped"]) for r in results)
    total_skipped = sum(r["skipped"] for r in results)
    total = sum(r["total"] for r in results)

    print("=" * 50)
    print(f"Summary: {total_mapped} mapped, {total_unmapped} unmapped, "
          f"{total_skipped} skipped out of {total} total rows")

    return 0


if __name__ == "__main__":
    exit(main())
