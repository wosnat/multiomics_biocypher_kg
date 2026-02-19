#!/usr/bin/env python
"""Build gene_mapping_supp.csv for each organism from paper supplementary CSVs.

Scans all paperconfig.yaml files, finds supplementary CSVs for each organism,
and collects alternative gene ID columns (columns other than name_col that
look like gene identifiers). Writes gene_mapping_supp.csv to each organism's
cache directory.

Output columns: locus_tag, alt_id, source_col, source_csv, paper

Usage:
    uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py
    uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py --organism "Prochlorococcus MED4"
    uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py --dry-run
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml


# Patterns to skip (tRNA, ncRNA, rRNA — intentionally not in gene nodes)
SKIP_PATTERNS = re.compile(r'^(tRNA|ncRNA|rRNA|Yfr\d|tmRNA)', re.IGNORECASE)

# Columns in gene_mapping.csv to build lookups from (same as fix-gene-ids)
MAPPING_COLS = [
    "locus_tag", "gene_names", "gene", "locus_tag_ncbi",
    "gene_names_cyanorak", "locus_tag_cyanoak", "old_locus_tags",
    "protein_id",
]

# Column name keywords that strongly suggest gene ID content
ID_COL_KEYWORDS = [
    'locus', 'tag', 'gene', 'protein', 'orf', 'id', 'accession', 'homolog',
]

# Column name keywords that indicate description/annotation (not gene IDs)
DESCRIPTION_COL_KEYWORDS = [
    'product', 'description', 'function', 'role', 'cog', 'note', 'go',
    'ontology', 'kegg', 'pathway', 'domain', 'tigrfam', 'pfam', 'eggnog',
    'cluster', 'inference', 'species', 'organism', 'strain', 'taxon',
    'annotation', 'regulation', 'temperature'
]

# Columns with this many or fewer unique non-null values are categorical (e.g. Up/Down)
MAX_UNIQUE_FOR_ID_COL = 10


def load_genome_registry(project_root):
    """Load cyanobacteria_genomes.csv and build strain_name -> (data_dir, organism_type) mapping.

    Returns dict: strain_name_lower -> data_dir (absolute path)
    Also returns a list of all (strain_name, data_dir) pairs.
    """
    csv_path = Path(project_root) / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
    if not csv_path.exists():
        print(f"ERROR: genome registry not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    # Read CSV, skipping comment lines starting with #
    rows = []
    with open(csv_path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = line.split(',')
            else:
                parts = line.split(',')
                if len(parts) >= len(header):
                    rows.append(dict(zip(header, parts)))

    entries = []
    for row in rows:
        strain_name = row.get('strain_name', '').strip()
        data_dir = row.get('data_dir', '').strip()
        if strain_name and data_dir:
            full_dir = Path(project_root) / data_dir
            entries.append((strain_name, str(full_dir)))

    return entries


def build_organism_to_genome_dir(genome_entries):
    """Build a fuzzy-match dict: normalized organism name substring -> data_dir.

    Keys are normalized (lowercase) strain name variants so we can match
    organism names from paperconfigs like "Prochlorococcus MED4" or
    "Prochlorococcus marinus MED4".
    """
    mapping = {}
    for strain_name, data_dir in genome_entries:
        key = strain_name.lower()
        mapping[key] = data_dir
    return mapping


def get_genome_dir(organism_name, organism_to_dir):
    """Return genome data_dir for organism_name, or None if not found."""
    if not organism_name:
        return None
    norm = organism_name.strip().lower()
    # Try exact match on strain name
    if norm in organism_to_dir:
        return organism_to_dir[norm]
    # Try substring match
    for key, data_dir in organism_to_dir.items():
        if key in norm or norm in key:
            return data_dir
    return None


def build_id_lookup(genome_dir):
    """Build a unified lookup dict: raw_id -> locus_tag from gene_mapping.csv.

    Same logic as fix-gene-ids. Handles space-separated gene_names and
    comma-separated old_locus_tags.

    Returns (lookup_dict, locus_tag_set) or (None, None) if file not found.
    """
    mapping_file = Path(genome_dir) / "gene_mapping.csv"
    if not mapping_file.exists():
        return None, None

    try:
        df = pd.read_csv(mapping_file)
    except Exception as e:
        print(f"  Error reading {mapping_file}: {e}", file=sys.stderr)
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
                # Space-separated list
                for name in val.split():
                    name = name.strip()
                    if name and name != "nan":
                        lookup[name] = locus
            elif col == "old_locus_tags":
                # Comma-separated list
                for tag in val.split(","):
                    tag = tag.strip()
                    if tag and tag != "nan":
                        lookup[tag] = locus
            else:
                lookup[val] = locus

    return lookup, locus_tags


def map_gene_id(raw_id, lookup, locus_tags):
    """Map a single gene ID to a locus_tag. Same logic as fix-gene-ids.

    Returns locus_tag string or None.
    """
    raw_id = raw_id.strip().strip('*').strip()
    if not raw_id or raw_id == "nan":
        return None

    if raw_id in locus_tags:
        return raw_id
    if raw_id in lookup:
        return lookup[raw_id]

    # Try zero-padding normalization
    m = re.match(r'^(.+)_(\d+)$', raw_id)
    if m:
        prefix, num = m.group(1), m.group(2)
        for pad in range(len(num) + 1, len(num) + 3):
            padded = f"{prefix}_{num.zfill(pad)}"
            if padded in locus_tags:
                return padded
            if padded in lookup:
                return lookup[padded]

    # Try comma-split composite IDs
    if "," in raw_id:
        for part in raw_id.split(","):
            part = part.strip()
            if part in locus_tags:
                return part
            if part in lookup:
                return lookup[part]

    return None


def is_id_like_column(series, col_name, exclude_cols):
    """Return True if a column looks like it contains gene identifiers.

    Args:
        series: pandas Series of column values
        col_name: column name string
        exclude_cols: set of column names to skip (name_col, logfc_col, etc.)
    """
    if col_name in exclude_cols:
        return False

    col_lower = col_name.lower()

    # Skip columns that are clearly description/annotation (not gene IDs)
    if any(kw in col_lower for kw in DESCRIPTION_COL_KEYWORDS):
        return False

    # Skip columns that are obviously coordinate/position columns
    if any(kw in col_lower for kw in ['start', 'end', 'length', 'position', 'coord']):
        return False

    # Check values
    vals = series.dropna().astype(str)
    vals = vals[vals.str.strip() != '']
    vals = vals[vals != 'nan']
    if len(vals) == 0:
        return False

    # Skip low-cardinality columns (categorical values like "Up"/"Down"):
    # few unique values AND each one repeats many times (avg >= 5 occurrences per value)
    n_unique = vals.nunique()
    if n_unique <= MAX_UNIQUE_FOR_ID_COL and len(vals) / n_unique >= 5:
        return False

    # Skip purely numeric columns
    try:
        pd.to_numeric(vals)
        return False
    except (ValueError, TypeError):
        pass

    # Skip columns with very long values (descriptions, GO terms)
    avg_len = vals.str.len().mean()
    if avg_len > 50:
        return False

    # Skip if most values look like floating point fold-changes
    # (values like "-0.54", "1.8*" etc.)
    numeric_like = vals.str.replace(r'[*\s]', '', regex=True)
    try:
        pd.to_numeric(numeric_like)
        return False
    except (ValueError, TypeError):
        pass

    # Column name keyword hint — after passing description/position exclusions
    if any(kw in col_lower for kw in ID_COL_KEYWORDS):
        return True

    # No keyword hint: require that most values look like gene IDs
    # (no spaces, reasonable format)
    space_frac = vals.str.contains(' ').mean()
    if space_frac > 0.3:
        # Most values contain spaces — likely a description column
        return False

    return True


def parse_paperconfigs(paperconfig_list_path, project_root):
    """Parse all paperconfig files and return analysis info.

    Returns list of dicts with keys:
        papername, organism, name_col, logfc_col, adjusted_p_value_col,
        csv_path, skip_rows
    """
    analyses = []

    with open(paperconfig_list_path) as f:
        paths = [line.strip() for line in f if line.strip()]

    for rel_path in paths:
        full_path = Path(project_root) / rel_path
        if not full_path.exists():
            print(f"WARNING: paperconfig not found: {full_path}", file=sys.stderr)
            continue

        try:
            with open(full_path) as f:
                config = yaml.safe_load(f)
        except Exception as e:
            print(f"WARNING: Cannot parse {full_path}: {e}", file=sys.stderr)
            continue

        pub = config.get("publication", {})
        papername = pub.get("papername", full_path.parent.name)
        supp_materials = pub.get("supplementary_materials", {})

        for table_key, table_data in supp_materials.items():
            if table_data.get("type") != "csv":
                continue

            filename = table_data.get("filename", "")
            if not filename:
                continue
            csv_path = Path(project_root) / filename

            table_skip = table_data.get("skip_rows", 0)

            for analysis in table_data.get("statistical_analyses", []):
                organism = analysis.get("organism", "")
                name_col = analysis.get("name_col", "")
                logfc_col = analysis.get("logfc_col", "")
                adj_pval_col = analysis.get("adjusted_p_value_col", "")
                analysis_skip = analysis.get("skip_rows", table_skip)

                if not organism or not name_col:
                    continue

                analyses.append({
                    "papername": papername,
                    "organism": organism,
                    "name_col": name_col,
                    "logfc_col": logfc_col,
                    "adjusted_p_value_col": adj_pval_col,
                    "csv_path": str(csv_path),
                    "skip_rows": analysis_skip or 0,
                })

    return analyses


def process_organism(genome_dir, organism_analyses, dry_run=False):
    """Process all analyses for one organism and return records.

    Args:
        genome_dir: path to the organism's cache directory
        organism_analyses: list of analysis dicts for this organism
        dry_run: if True, don't write files

    Returns list of dicts: {locus_tag, alt_id, source_col, source_csv, paper}
    """
    # Load gene_mapping.csv lookup
    lookup, locus_tags = build_id_lookup(genome_dir)
    if lookup is None:
        print(f"  No gene_mapping.csv found in {genome_dir}, skipping", file=sys.stderr)
        return []

    records = []
    seen = set()  # Dedup key: (locus_tag, alt_id, source_col, source_csv, paper)

    # Track per-CSV processing (avoid duplicate work for same CSV+organism+name_col)
    processed_csvs = set()

    for analysis in organism_analyses:
        papername = analysis["papername"]
        name_col = analysis["name_col"]
        logfc_col = analysis.get("logfc_col", "")
        adj_pval_col = analysis.get("adjusted_p_value_col", "")
        csv_path = Path(analysis["csv_path"])
        skip_rows = analysis["skip_rows"]

        cache_key = (str(csv_path), name_col)
        if cache_key in processed_csvs:
            continue
        processed_csvs.add(cache_key)

        if not csv_path.exists():
            print(f"  WARNING: CSV not found: {csv_path}", file=sys.stderr)
            continue

        try:
            df = pd.read_csv(str(csv_path), skiprows=skip_rows, low_memory=False)
        except Exception as e:
            print(f"  WARNING: Cannot read {csv_path}: {e}", file=sys.stderr)
            continue

        if name_col not in df.columns:
            print(
                f"  WARNING: name_col '{name_col}' not in {csv_path.name}. "
                f"Columns: {list(df.columns)[:5]}...",
                file=sys.stderr,
            )
            continue

        # Determine excluded columns
        exclude_cols = {name_col}
        if logfc_col:
            exclude_cols.add(logfc_col)
        if adj_pval_col:
            exclude_cols.add(adj_pval_col)

        # Find ID-like columns
        id_cols = [
            col for col in df.columns
            if col != name_col and is_id_like_column(df[col], col, exclude_cols)
        ]

        if not id_cols:
            continue

        csv_basename = csv_path.name
        mapped_count = 0
        skipped_count = 0

        for _, row in df.iterrows():
            raw_id = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""
            if not raw_id or raw_id == "nan":
                continue

            # Skip RNA features
            if SKIP_PATTERNS.match(raw_id):
                skipped_count += 1
                continue

            locus_tag = map_gene_id(raw_id, lookup, locus_tags)
            if not locus_tag:
                continue

            mapped_count += 1

            for col in id_cols:
                alt_val = str(row[col]).strip().strip('*').strip() if pd.notna(row[col]) else ""
                if not alt_val or alt_val == "nan":
                    continue

                dedup_key = (locus_tag, alt_val, col, csv_basename, papername)
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)

                records.append({
                    "locus_tag": locus_tag,
                    "alt_id": alt_val,
                    "source_col": col,
                    "source_csv": csv_basename,
                    "paper": papername,
                })

        print(
            f"    {csv_path.name}: {mapped_count} rows resolved, "
            f"{skipped_count} skipped (RNA), "
            f"{len([r for r in records if r['source_csv'] == csv_basename and r['paper'] == papername])} new records",
            file=sys.stderr,
        )

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Build gene_mapping_supp.csv for each organism from paper supplementary CSVs."
    )
    parser.add_argument(
        "--organism",
        help="Only process this organism (fuzzy match on strain name, e.g. 'MED4', 'Prochlorococcus MED4')",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show stats without writing files",
    )
    parser.add_argument(
        "--paperconfig-list",
        default="data/Prochlorococcus/papers_and_supp/paperconfig_files.txt",
        help="Path to paperconfig_files.txt",
    )
    args = parser.parse_args()

    # Project root is 3 levels up from this script
    project_root = Path(__file__).resolve().parents[3]

    # Load genome registry
    genome_entries = load_genome_registry(project_root)
    organism_to_dir = build_organism_to_genome_dir(genome_entries)

    # Parse all paperconfigs
    list_path = Path(project_root) / args.paperconfig_list
    if not list_path.exists():
        print(f"ERROR: paperconfig list not found: {list_path}", file=sys.stderr)
        sys.exit(1)

    all_analyses = parse_paperconfigs(list_path, project_root)
    print(
        f"Parsed {len(all_analyses)} analyses from paperconfigs",
        file=sys.stderr,
    )

    # Group analyses by genome_dir
    dir_to_analyses = defaultdict(list)
    unmatched_organisms = set()

    for analysis in all_analyses:
        genome_dir = get_genome_dir(analysis["organism"], organism_to_dir)
        if genome_dir:
            dir_to_analyses[genome_dir].append(analysis)
        else:
            unmatched_organisms.add(analysis["organism"])

    if unmatched_organisms:
        print("Organisms with no genome directory (skipped):", file=sys.stderr)
        for org in sorted(unmatched_organisms):
            print(f"  {org}", file=sys.stderr)

    # Filter to requested organism if specified
    if args.organism:
        target_dir = get_genome_dir(args.organism, organism_to_dir)
        if not target_dir:
            print(
                f"ERROR: No genome directory found for organism '{args.organism}'",
                file=sys.stderr,
            )
            sys.exit(1)
        dirs_to_process = {target_dir: dir_to_analyses[target_dir]}
        print(f"Processing single organism: {args.organism} -> {target_dir}", file=sys.stderr)
    else:
        dirs_to_process = dict(dir_to_analyses)

    # Process each organism
    total_records = 0

    for genome_dir, analyses in sorted(dirs_to_process.items()):
        # Derive a display name from the path
        dir_path = Path(genome_dir)
        display_name = f"{dir_path.parent.parent.name}/{dir_path.name}"

        papers = sorted({a["papername"] for a in analyses})
        print(
            f"\nOrganism: {display_name} "
            f"({len(analyses)} analyses from {len(papers)} paper(s))",
            file=sys.stderr,
        )

        records = process_organism(genome_dir, analyses, dry_run=args.dry_run)

        if not records:
            print(f"  No records produced", file=sys.stderr)
            continue

        total_records += len(records)

        # Sort by locus_tag, then source
        df_out = pd.DataFrame(records)
        df_out = df_out.sort_values(["locus_tag", "paper", "source_col"])
        df_out = df_out.drop_duplicates()

        # Stats
        unique_loci = df_out["locus_tag"].nunique()
        by_paper = df_out.groupby("paper").size()
        print(f"  Records: {len(df_out)} total, {unique_loci} unique locus_tags", file=sys.stderr)
        for paper, count in by_paper.items():
            print(f"    {paper}: {count}", file=sys.stderr)

        if not args.dry_run:
            out_path = Path(genome_dir) / "gene_mapping_supp.csv"
            df_out.to_csv(out_path, index=False)
            print(f"  Written: {out_path}", file=sys.stderr)
            print(f"Written: {out_path}")
        else:
            print(f"  (dry run — would write {Path(genome_dir) / 'gene_mapping_supp.csv'})")

    print(f"\nTotal records: {total_records}", file=sys.stderr)


if __name__ == "__main__":
    main()
