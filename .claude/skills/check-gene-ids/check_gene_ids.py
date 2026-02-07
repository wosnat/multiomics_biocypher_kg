#!/usr/bin/env python
"""Validate gene ID matching between paper CSV data and gene nodes in the knowledge graph.

For each publication and supplementary table, checks whether gene IDs in the
CSV name_col would match existing Gene nodes. Reports per-table status and
suggests a fix strategy when IDs don't match.

Usage:
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/.../paperconfig.yaml"
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py --biocypher-dir biocypher-out/20260205200505
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml


# Number of sample IDs to check per analysis (no need to check every gene)
SAMPLE_SIZE = 30

# Organism name -> genome directory mapping
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

# Patterns to skip (tRNA, ncRNA, rRNA — intentionally not in gene nodes)
SKIP_PATTERNS = re.compile(r'^(tRNA|ncRNA|rRNA|Yfr\d|tmRNA)', re.IGNORECASE)

# Columns in gene_mapping.csv that could contain gene identifiers
GENE_MAPPING_ID_COLS = [
    "gene_names", "gene", "locus_tag_ncbi", "locus_tag",
    "protein_id", "locus_tag_cyanoak", "gene_names_cyanorak",
    "old_locus_tags",
]


def find_latest_biocypher_dir(base_dir="biocypher-out"):
    """Find the most recent biocypher output directory by timestamp name."""
    base = Path(base_dir)
    if not base.exists():
        return None
    dirs = sorted(
        [d for d in base.iterdir() if d.is_dir() and d.name.isdigit()],
        key=lambda d: d.name,
        reverse=True,
    )
    return str(dirs[0]) if dirs else None


def strip_quotes(val):
    """Strip BioCypher single-quote wrappers from a string value."""
    if isinstance(val, str) and len(val) >= 2 and val.startswith("'") and val.endswith("'"):
        return val[1:-1]
    return val


def load_gene_node_index(biocypher_dir):
    """Load gene node IDs and alternative identifier columns from Gene-part000.csv.

    Returns:
        dict with keys:
        - primary_ids: set of :ID values (e.g., "ncbigene:PMM0001")
        - locus_tag: dict of locus_tag -> primary_id
        - locus_tag_ncbi: dict of locus_tag_ncbi -> primary_id
        - locus_tag_cyanorak: dict of locus_tag_cyanorak -> primary_id
        - gene_names: dict of individual gene name -> primary_id
    """
    gene_file = Path(biocypher_dir) / "Gene-part000.csv"
    if not gene_file.exists():
        print(f"ERROR: Gene file not found: {gene_file}", file=sys.stderr)
        sys.exit(1)

    index = {
        "primary_ids": set(),
        "locus_tag": {},
        "locus_tag_ncbi": {},
        "locus_tag_cyanorak": {},
        "gene_names": {},
    }

    # Read header to find column indices
    header_file = Path(biocypher_dir) / "Gene-header.csv"
    with open(header_file) as f:
        headers = f.readline().strip().split("\t")

    col_idx = {h.split(":")[0]: i for i, h in enumerate(headers)}

    with open(gene_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if not fields:
                continue

            primary_id = fields[0]  # :ID column, not quoted
            index["primary_ids"].add(primary_id)

            # locus_tag
            if "locus_tag" in col_idx:
                val = strip_quotes(fields[col_idx["locus_tag"]])
                if val:
                    index["locus_tag"][val] = primary_id

            # locus_tag_ncbi
            if "locus_tag_ncbi" in col_idx:
                val = strip_quotes(fields[col_idx["locus_tag_ncbi"]])
                if val:
                    index["locus_tag_ncbi"][val] = primary_id

            # locus_tag_cyanorak
            if "locus_tag_cyanorak" in col_idx:
                val = strip_quotes(fields[col_idx["locus_tag_cyanorak"]])
                if val:
                    index["locus_tag_cyanorak"][val] = primary_id

            # gene_names (string[] field, pipe-delimited within quotes)
            if "gene_names" in col_idx:
                raw = strip_quotes(fields[col_idx["gene_names"]])
                if raw:
                    for name in raw.split("|"):
                        name = name.strip()
                        if name:
                            index["gene_names"][name] = primary_id

    return index


def is_organism_loaded(organism_name):
    """Check if an organism has a genome loaded (has a known genome directory)."""
    if not organism_name:
        return False
    norm = organism_name.strip().strip('"').lower()
    # Direct match
    if norm in ORGANISM_TO_GENOME_DIR:
        return True
    # Check if any known organism name is a substring
    for known in ORGANISM_TO_GENOME_DIR:
        if known in norm or norm in known:
            return True
    return False


def get_genome_dir(organism_name, project_root):
    """Get the genome directory for an organism, or None."""
    if not organism_name:
        return None
    norm = organism_name.strip().strip('"').lower()
    for key, rel_path in ORGANISM_TO_GENOME_DIR.items():
        if key in norm or norm in key:
            full = Path(project_root) / rel_path
            if full.exists():
                return str(full)
    return None


def scan_gff_for_ids(gff_path, sample_ids):
    """Scan a GFF file for sample IDs in any attribute field.

    Returns dict mapping attribute_name -> count of sample_ids found there.
    """
    matches = defaultdict(int)
    sample_set = set(sample_ids)
    if not sample_set or not Path(gff_path).exists():
        return matches

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs_str = parts[8]
            for attr in attrs_str.split(";"):
                if "=" not in attr:
                    continue
                key, val = attr.split("=", 1)
                for sid in sample_set:
                    if sid in val:
                        matches[key] += 1
    return matches


def check_other_csv_columns(csv_path, skip_rows, gene_index, exclude_col):
    """Check if any other column in the CSV has IDs matching gene nodes.

    Returns list of (column_name, sample_value, match_count) for columns that match.
    """
    try:
        df = pd.read_csv(csv_path, skiprows=skip_rows or 0, nrows=SAMPLE_SIZE)
    except Exception:
        return []

    results = []
    for col in df.columns:
        if col == exclude_col:
            continue
        vals = df[col].dropna().astype(str).head(SAMPLE_SIZE).tolist()
        if not vals:
            continue
        match_count = sum(
            1 for v in vals
            if f"ncbigene:{v}" in gene_index["primary_ids"]
        )
        if match_count > len(vals) * 0.5:  # >50% match
            sample_val = vals[0] if vals else ""
            results.append((col, sample_val, match_count, len(vals)))
    return results


def check_alt_gene_fields(sample_ids, gene_index):
    """Check if sample IDs match alternative gene node fields.

    Returns dict: field_name -> (match_count, sample_match)
    """
    results = {}
    for field in ["locus_tag", "locus_tag_ncbi", "locus_tag_cyanorak", "gene_names"]:
        lookup = gene_index[field]
        matched = [(sid, lookup[sid]) for sid in sample_ids if sid in lookup]
        if matched:
            results[field] = (len(matched), matched[0])
    return results


def load_gene_mapping(genome_dir):
    """Load gene_mapping.csv and build lookup dicts for all ID-like columns.

    Returns dict: col_name -> {id_value -> locus_tag}, or None if file not found.
    """
    mapping_file = Path(genome_dir) / "gene_mapping.csv"
    if not mapping_file.exists():
        return None

    try:
        df = pd.read_csv(mapping_file)
    except Exception:
        return None

    lookups = {}
    for col in GENE_MAPPING_ID_COLS:
        if col not in df.columns:
            continue
        lookup = {}
        for _, row in df.iterrows():
            val = str(row[col]).strip() if pd.notna(row[col]) else ""
            locus = str(row["locus_tag"]).strip() if pd.notna(row.get("locus_tag")) else ""
            if val and locus and val != "nan" and locus != "nan":
                lookup[val] = locus
        if lookup:
            lookups[col] = lookup
    return lookups


def check_gene_mapping(sample_ids, genome_dir):
    """Check if sample IDs match any column in gene_mapping.csv.

    Returns dict: col_name -> (match_count, sample_match_tuple), or None if not found.
    """
    lookups = load_gene_mapping(genome_dir)
    if not lookups:
        return None

    results = {}
    for col, lookup in lookups.items():
        matched = [(sid, lookup[sid]) for sid in sample_ids if sid in lookup]
        if matched:
            results[col] = (len(matched), matched[0])
    return results


def analyze_paperconfig(paperconfig_path, gene_index, project_root):
    """Analyze a single paperconfig.yaml and return per-analysis results.

    Returns:
        dict with papername and list of analysis results
    """
    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    pub = config.get("publication", {})
    papername = pub.get("papername", Path(paperconfig_path).parent.name)
    supp_materials = pub.get("supplementary_materials", {})

    analyses_results = []

    for table_key, table_data in supp_materials.items():
        filename = table_data.get("filename", "")
        csv_path = Path(project_root) / filename if filename else None
        stat_analyses = table_data.get("statistical_analyses", [])

        for analysis in stat_analyses:
            name_col = analysis.get("name_col", "")
            organism = analysis.get("organism", "")
            analysis_id = analysis.get("id", table_key)
            skip_rows = analysis.get("skip_rows") or table_data.get("skip_rows")
            csv_basename = Path(filename).name if filename else table_key

            result = {
                "table": csv_basename,
                "analysis_id": analysis_id,
                "organism": organism,
                "name_col": name_col,
                "status": "UNKNOWN",
                "fix_strategy": None,
                "details": "",
                "sample_ids": [],
            }

            # Read CSV and extract sample IDs
            if not csv_path or not csv_path.exists():
                result["status"] = "ERROR"
                result["details"] = f"CSV file not found: {filename}"
                analyses_results.append(result)
                continue

            try:
                df = pd.read_csv(str(csv_path), skiprows=skip_rows or 0, nrows=SAMPLE_SIZE)
            except Exception as e:
                result["status"] = "ERROR"
                result["details"] = f"Cannot read CSV: {e}"
                analyses_results.append(result)
                continue

            if name_col not in df.columns:
                result["status"] = "ERROR"
                result["details"] = f"Column '{name_col}' not found. Available: {list(df.columns)}"
                analyses_results.append(result)
                continue

            # Get sample IDs, filter out blanks and RNA features
            raw_ids = df[name_col].dropna().astype(str).tolist()
            raw_ids = [rid.strip() for rid in raw_ids if rid.strip()]
            # Filter out tRNA/ncRNA/rRNA
            raw_ids = [rid for rid in raw_ids if not SKIP_PATTERNS.match(rid)]

            if not raw_ids:
                result["status"] = "ERROR"
                result["details"] = "No valid gene IDs found in name_col after filtering"
                analyses_results.append(result)
                continue

            sample_ids = raw_ids[:SAMPLE_SIZE]
            result["sample_ids"] = sample_ids[:5]  # Store a few for display

            # Check primary ID match
            match_count = sum(
                1 for sid in sample_ids
                if f"ncbigene:{sid}" in gene_index["primary_ids"]
            )
            match_rate = match_count / len(sample_ids)

            if match_rate >= 0.8:
                result["status"] = "MATCH"
                result["details"] = f"{match_count}/{len(sample_ids)} IDs match gene nodes"
                analyses_results.append(result)
                continue

            # --- IDs don't match. Determine why and suggest fix. ---

            # Check if organism genome is loaded
            if not is_organism_loaded(organism):
                result["status"] = "NO MATCH"
                result["fix_strategy"] = "LOAD_ORGANISM"
                result["details"] = (
                    f"Organism '{organism}' genome not loaded. "
                    f"Add to cyanobacteria_genomes.csv to create gene nodes."
                )
                analyses_results.append(result)
                continue

            # Organism IS loaded but IDs don't match.
            # Check if another column in the CSV has matching IDs
            alt_cols = check_other_csv_columns(
                str(csv_path), skip_rows, gene_index, name_col
            )
            if alt_cols:
                best_col, sample_val, cnt, total = alt_cols[0]
                result["status"] = "NO MATCH"
                result["fix_strategy"] = "CHANGE_NAME_COL"
                result["details"] = (
                    f"Column '{best_col}' has IDs like '{sample_val}' "
                    f"that match gene nodes ({cnt}/{total}). "
                    f"Suggest: name_col: \"{best_col}\""
                )
                analyses_results.append(result)
                continue

            # Check gene_mapping.csv for the organism
            genome_dir = get_genome_dir(organism, project_root)
            if genome_dir:
                mapping_matches = check_gene_mapping(sample_ids, genome_dir)
                if mapping_matches:
                    best_col = max(mapping_matches, key=lambda k: mapping_matches[k][0])
                    cnt, (sample_id, mapped_to) = mapping_matches[best_col]
                    mapping_file = Path(genome_dir) / "gene_mapping.csv"
                    rel_mapping = os.path.relpath(mapping_file, project_root)
                    result["status"] = "NO MATCH"
                    result["fix_strategy"] = "CREATE_MAPPING_CSV"
                    result["details"] = (
                        f"IDs match column '{best_col}' in gene_mapping.csv "
                        f"({cnt}/{len(sample_ids)} match). "
                        f"E.g., '{sample_id}' -> locus_tag '{mapped_to}'. "
                        f"Run /fix-gene-ids on this paperconfig to create a "
                        f"_with_locus_tag.csv, then update name_col to 'locus_tag'."
                    )
                    analyses_results.append(result)
                    continue

                # gene_mapping.csv exists but IDs not found there — check GFF/GBK
                gff_files = list(Path(genome_dir).rglob("*.gff"))
                for gff_path in gff_files:
                    gff_matches = scan_gff_for_ids(str(gff_path), sample_ids[:10])
                    if gff_matches:
                        best_attr = max(gff_matches, key=gff_matches.get)
                        cnt = gff_matches[best_attr]
                        rel_mapping = os.path.relpath(
                            Path(genome_dir) / "gene_mapping.csv", project_root
                        )
                        result["status"] = "NO MATCH"
                        result["fix_strategy"] = "CREATE_MAPPING_GFF"
                        result["details"] = (
                            f"IDs found in GFF attribute '{best_attr}' "
                            f"in {gff_path.name} ({cnt} matches) "
                            f"but NOT in gene_mapping.csv. "
                            f"Add '{best_attr}' column to {rel_mapping}."
                        )
                        analyses_results.append(result)
                        break
                else:
                    result["status"] = "NO MATCH"
                    result["fix_strategy"] = "UNRELATED_IDS"
                    result["details"] = (
                        f"IDs not found in gene_mapping.csv or GFF/GBK files for {organism}. "
                        f"Manual investigation needed."
                    )
                    analyses_results.append(result)
            else:
                # No genome dir — fall back to gene node field check
                alt_fields = check_alt_gene_fields(sample_ids, gene_index)
                if alt_fields:
                    best_field = max(alt_fields, key=lambda k: alt_fields[k][0])
                    cnt, (sample_id, mapped_to) = alt_fields[best_field]
                    result["status"] = "NO MATCH"
                    result["fix_strategy"] = "CREATE_MAPPING_CSV"
                    result["details"] = (
                        f"IDs match gene node field '{best_field}' "
                        f"({cnt}/{len(sample_ids)} match). "
                        f"E.g., '{sample_id}' -> {mapped_to}. "
                        f"No gene_mapping.csv found — create one for this organism."
                    )
                else:
                    result["status"] = "NO MATCH"
                    result["fix_strategy"] = "UNRELATED_IDS"
                    result["details"] = (
                        f"IDs don't match gene nodes and no genome directory found for {organism}. "
                        f"Manual investigation needed."
                    )
                analyses_results.append(result)

            # Handle partial match (some IDs match, some don't)
            if match_rate > 0 and result["status"] == "NO MATCH":
                result["status"] = "PARTIAL MATCH"

    return {"papername": papername, "paperconfig_path": paperconfig_path, "analyses": analyses_results}


def format_report(all_results):
    """Format the full report as human-readable text."""
    lines = []
    lines.append("=" * 60)
    lines.append("Gene ID Validation Report")
    lines.append("=" * 60)
    lines.append("")

    summary_rows = []

    for paper in all_results:
        lines.append(f"## {paper['papername']}")
        lines.append(f"   Config: {paper['paperconfig_path']}")
        lines.append("")

        for a in paper["analyses"]:
            lines.append(f"  Table: {a['table']} (analysis: {a['analysis_id']})")
            lines.append(f"    Organism: {a['organism']} | name_col: {a['name_col']}")
            if a["sample_ids"]:
                lines.append(f"    Sample IDs: {', '.join(a['sample_ids'][:5])}")
            lines.append(f"    Status: {a['status']}")
            if a["fix_strategy"]:
                lines.append(f"    Fix strategy: {a['fix_strategy']} — {a['details']}")
            elif a["details"]:
                lines.append(f"    {a['details']}")
            lines.append("")

            summary_rows.append({
                "publication": paper["papername"],
                "table": a["table"],
                "status": a["status"],
                "fix_strategy": a["fix_strategy"] or "—",
            })

    # Summary table
    lines.append("")
    lines.append("=" * 60)
    lines.append("Summary")
    lines.append("=" * 60)
    lines.append("")

    # Column widths
    pub_w = max(len("Publication"), max((len(r["publication"]) for r in summary_rows), default=0))
    tbl_w = max(len("Table"), max((min(len(r["table"]), 40) for r in summary_rows), default=0))
    stat_w = max(len("Status"), max((len(r["status"]) for r in summary_rows), default=0))
    fix_w = max(len("Fix Strategy"), max((len(r["fix_strategy"]) for r in summary_rows), default=0))

    header = f"{'Publication':<{pub_w}} | {'Table':<{tbl_w}} | {'Status':<{stat_w}} | {'Fix Strategy':<{fix_w}}"
    lines.append(header)
    lines.append("-" * len(header))

    for r in summary_rows:
        table_short = r["table"][:40]
        lines.append(
            f"{r['publication']:<{pub_w}} | {table_short:<{tbl_w}} | {r['status']:<{stat_w}} | {r['fix_strategy']:<{fix_w}}"
        )

    # Stats
    total = len(summary_rows)
    matched = sum(1 for r in summary_rows if r["status"] == "MATCH")
    no_match = sum(1 for r in summary_rows if r["status"] in ("NO MATCH", "PARTIAL MATCH"))
    errors = sum(1 for r in summary_rows if r["status"] == "ERROR")

    lines.append("")
    lines.append(f"Total analyses: {total}")
    lines.append(f"  MATCH: {matched}")
    lines.append(f"  NO MATCH / PARTIAL: {no_match}")
    if errors:
        lines.append(f"  ERROR: {errors}")

    # Strategy breakdown
    strategies = defaultdict(int)
    for r in summary_rows:
        if r["fix_strategy"] != "—":
            strategies[r["fix_strategy"]] += 1
    if strategies:
        lines.append("")
        lines.append("Fix strategies needed:")
        for strat, count in sorted(strategies.items(), key=lambda x: -x[1]):
            lines.append(f"  {strat}: {count}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Validate gene ID matching between paper CSVs and gene nodes."
    )
    parser.add_argument(
        "--biocypher-dir",
        help="Path to biocypher output directory (default: latest in biocypher-out/)",
    )
    parser.add_argument(
        "--paperconfig",
        help="Check a single paperconfig.yaml (default: all from paperconfig_files.txt)",
    )
    parser.add_argument(
        "--paperconfig-list",
        default="data/Prochlorococcus/papers_and_supp/paperconfig_files.txt",
        help="Path to paperconfig_files.txt",
    )
    args = parser.parse_args()

    # Determine project root (directory containing biocypher-out/)
    project_root = Path(__file__).resolve().parents[3]  # .claude/skills/check-gene-ids/ -> project root
    os.chdir(project_root)

    # Find biocypher output dir
    bc_dir = args.biocypher_dir or find_latest_biocypher_dir("biocypher-out")
    if not bc_dir:
        print("ERROR: No biocypher output directory found.", file=sys.stderr)
        sys.exit(1)
    print(f"Using biocypher output: {bc_dir}", file=sys.stderr)

    # Load gene node index
    gene_index = load_gene_node_index(bc_dir)
    print(f"Loaded {len(gene_index['primary_ids'])} gene nodes", file=sys.stderr)

    # Determine which paperconfigs to check
    if args.paperconfig:
        paperconfig_paths = [args.paperconfig]
    else:
        list_file = Path(args.paperconfig_list)
        if not list_file.exists():
            print(f"ERROR: paperconfig list not found: {list_file}", file=sys.stderr)
            sys.exit(1)
        with open(list_file) as f:
            paperconfig_paths = [line.strip() for line in f if line.strip()]

    # Analyze each paper
    all_results = []
    for pc_path in paperconfig_paths:
        if not Path(pc_path).exists():
            print(f"WARNING: Paperconfig not found: {pc_path}", file=sys.stderr)
            continue
        result = analyze_paperconfig(pc_path, gene_index, project_root)
        all_results.append(result)

    # Print report to stdout and write to file
    report = format_report(all_results)
    print(report)

    report_file = project_root / "gene_id_validation_report.txt"
    with open(report_file, "w") as f:
        f.write(report)
        f.write("\n")
    print(f"\nReport written to {report_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
