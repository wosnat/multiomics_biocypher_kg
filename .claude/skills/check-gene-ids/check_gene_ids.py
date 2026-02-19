#!/usr/bin/env python
"""Validate gene ID matching between paper CSV data and gene nodes in the knowledge graph.

For each publication and supplementary table, checks ALL gene IDs in the
CSV name_col against existing Gene nodes. Reports per-table status with
detailed breakdowns (matched, RNA/tRNA/ncRNA skipped, mismatched) and
cross-references against the Docker import report if available.

Usage:
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/.../paperconfig.yaml"
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py --biocypher-dir biocypher-out/20260205200505
    uv run python .claude/skills/check-gene-ids/check_gene_ids.py --import-report output/import.report
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml


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
    "alteromonas macleodii hot1a3": "cache/data/Alteromonas/genomes/HOT1A3",
    "alteromonas hot1a3": "cache/data/Alteromonas/genomes/HOT1A3",
}

# Patterns for RNA features (tRNA, ncRNA, rRNA — intentionally not in gene nodes)
RNA_PATTERN = re.compile(
    r'^(tRNA|ncRNA|rRNA|Yfr\d|tmRNA|RNA_\d|PMT_ncRNA_|'
    r'\w+_tRNA\w+VIMSS|'  # e.g., A9601_tRNAAlaVIMSS1309073
    r'\w+_rr[ls]VIMSS)',  # e.g., A9601_rrlVIMSS1365720
    re.IGNORECASE
)

# Columns in gene_mapping.csv that could contain gene identifiers
GENE_MAPPING_ID_COLS = [
    "gene_names", "gene", "locus_tag_ncbi", "locus_tag",
    "protein_id", "locus_tag_cyanoak", "gene_names_cyanorak",
    "old_locus_tags",
]

# Sample size for expensive lookups (GFF scanning, alt column checks)
LOOKUP_SAMPLE_SIZE = 50


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


def load_import_report(import_report_path=None):
    """Load missing gene IDs from the Docker import report.

    Tries multiple sources:
    1. Explicit --import-report path
    2. Local file output/import.report
    3. Docker: docker compose exec deploy cat /data/build2neo/import.report

    Returns:
        dict: analysis_id_suffix -> set of missing ncbigene IDs (without prefix)
        The analysis_id_suffix is the part after the DOI, e.g., "fang_2019_vdom_addition"
    """
    report_lines = []

    # Try explicit path first
    if import_report_path and Path(import_report_path).exists():
        with open(import_report_path) as f:
            report_lines = f.readlines()
        print(f"Loaded import report from {import_report_path}", file=sys.stderr)
    else:
        # Try local file
        local_report = Path("output/import.report")
        if local_report.exists():
            with open(local_report) as f:
                report_lines = f.readlines()
            print(f"Loaded import report from {local_report}", file=sys.stderr)
        else:
            # Try Docker
            try:
                result = subprocess.run(
                    ["docker", "compose", "exec", "deploy", "cat", "/data/build2neo/import.report"],
                    capture_output=True, text=True, timeout=15
                )
                if result.returncode == 0 and result.stdout.strip():
                    report_lines = result.stdout.strip().split("\n")
                    print(f"Loaded import report from Docker ({len(report_lines)} lines)", file=sys.stderr)
            except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
                print(f"Could not load import report from Docker: {e}", file=sys.stderr)

    if not report_lines:
        return None

    # Parse: extract Affects_expression_of edges
    # Format: <source> (global id space)-[Affects_expression_of]-><target> (global id space) referring to missing node <target>
    missing_by_source = defaultdict(set)
    affects_pattern = re.compile(
        r'^(.+?) \(global id space\)-\[Affects_expression_of\]->ncbigene:(.+?) \(global id space\)'
    )

    for line in report_lines:
        m = affects_pattern.match(line.strip())
        if m:
            source_id = m.group(1)
            missing_gene = m.group(2)
            missing_by_source[source_id].add(missing_gene)

    # Also build a flat set of all missing gene IDs
    all_missing = set()
    for ids in missing_by_source.values():
        all_missing.update(ids)

    return {
        "by_source": missing_by_source,
        "all_missing": all_missing,
    }


def is_organism_loaded(organism_name):
    """Check if an organism has a genome loaded (has a known genome directory)."""
    if not organism_name:
        return False
    norm = organism_name.strip().strip('"').lower()
    if norm in ORGANISM_TO_GENOME_DIR:
        return True
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

    Returns list of (column_name, sample_value, match_count, total) for columns that match.
    """
    try:
        df = pd.read_csv(csv_path, skiprows=skip_rows or 0, nrows=LOOKUP_SAMPLE_SIZE)
    except Exception:
        return []

    results = []
    for col in df.columns:
        if col == exclude_col:
            continue
        vals = df[col].dropna().astype(str).head(LOOKUP_SAMPLE_SIZE).tolist()
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


def classify_ids(all_ids, gene_index, secondary_ids=None):
    """Classify all gene IDs into categories.

    Args:
        all_ids: list of gene ID strings from the CSV
        gene_index: gene node index dict
        secondary_ids: optional list of secondary gene IDs (same length as all_ids)

    Returns dict with:
        - matched: list of IDs that match gene nodes (via primary ncbigene:<id>)
        - matched_secondary: list of IDs matched via secondary column
        - rna_skipped: list of IDs matching RNA/tRNA/ncRNA patterns
        - mismatched: list of IDs that don't match anything
        - empty: count of blank/NA IDs filtered out
    """
    result = {
        "matched": [],
        "matched_secondary": [],
        "rna_skipped": [],
        "mismatched": [],
        "empty": 0,
    }

    primary_ids = gene_index["primary_ids"]

    for i, raw_id in enumerate(all_ids):
        sid = str(raw_id).strip() if pd.notna(raw_id) else ""
        if not sid or sid == "nan":
            result["empty"] += 1
            continue

        # Check RNA pattern
        if RNA_PATTERN.match(sid):
            result["rna_skipped"].append(sid)
            continue

        # Check primary match
        if f"ncbigene:{sid}" in primary_ids:
            result["matched"].append(sid)
            continue

        # Check secondary column match
        if secondary_ids is not None and i < len(secondary_ids):
            sec_id = str(secondary_ids[i]).strip() if pd.notna(secondary_ids[i]) else ""
            if sec_id and sec_id != "nan" and f"ncbigene:{sec_id}" in primary_ids:
                result["matched_secondary"].append(sid)
                continue

        result["mismatched"].append(sid)

    return result


def find_import_report_mismatches(analysis, import_report):
    """Find mismatches for this analysis in the import report.

    The import report source IDs are built by the omics adapter as:
    - "{doi}_{env_condition_id}" for environmental treatment edges
    - "insdc.gcf:{accession}" for organism-sourced edges
    - "ncbitaxon:{taxid}" for taxid-based organism edges

    Args:
        analysis: the analysis dict from paperconfig
        import_report: parsed import report dict or None

    Returns:
        set of missing gene IDs from import report for this analysis, or None
    """
    if not import_report:
        return None

    by_source = import_report["by_source"]
    matched_missing = set()

    # Build possible source keys this analysis could produce
    env_condition_id = analysis.get("environmental_treatment_condition_id", "")
    treatment_accession = analysis.get("treatment_assembly_accession", "")
    treatment_taxid = analysis.get("treatment_taxid", "")
    analysis_id = analysis.get("id", "")

    for source_key, missing_ids in by_source.items():
        # Match by env_condition_id (most common)
        if env_condition_id and env_condition_id in source_key:
            matched_missing.update(missing_ids)
        # Match by treatment accession
        elif treatment_accession and treatment_accession in source_key:
            matched_missing.update(missing_ids)
        # Match by treatment taxid
        elif treatment_taxid and str(treatment_taxid) in source_key:
            matched_missing.update(missing_ids)
        # Fallback: match by analysis_id
        elif analysis_id and analysis_id in source_key:
            matched_missing.update(missing_ids)
        # Note: not matching by DOI alone — too broad for papers with multiple analyses

    return matched_missing if matched_missing else None


def analyze_paperconfig(paperconfig_path, gene_index, project_root, import_report=None):
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
            secondary_name_col = analysis.get("secondary_name_col", "")
            organism = analysis.get("organism", "")
            analysis_id = analysis.get("id", table_key)
            skip_rows = analysis.get("skip_rows") or table_data.get("skip_rows")
            csv_basename = Path(filename).name if filename else table_key

            result = {
                "table": csv_basename,
                "analysis_id": analysis_id,
                "organism": organism,
                "name_col": name_col,
                "secondary_name_col": secondary_name_col,
                "status": "UNKNOWN",
                "fix_strategy": None,
                "details": "",
                "sample_ids": [],
                "total_ids": 0,
                "matched_count": 0,
                "matched_secondary_count": 0,
                "rna_count": 0,
                "mismatched_count": 0,
                "empty_count": 0,
                "mismatched_ids": [],
                "import_report_missing": None,
                "import_report_count": 0,
            }

            # Read CSV — read ALL rows
            if not csv_path or not csv_path.exists():
                result["status"] = "ERROR"
                result["details"] = f"CSV file not found: {filename}"
                analyses_results.append(result)
                continue

            try:
                df = pd.read_csv(str(csv_path), skiprows=skip_rows or 0)
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

            # Get ALL IDs from the CSV
            all_ids = df[name_col].tolist()
            result["total_ids"] = len(all_ids)

            # Get secondary IDs if configured
            secondary_ids = None
            if secondary_name_col and secondary_name_col in df.columns:
                secondary_ids = df[secondary_name_col].tolist()
            elif secondary_name_col and secondary_name_col not in df.columns:
                result["details"] = f"Warning: secondary_name_col '{secondary_name_col}' not found in CSV. "

            # Classify all IDs
            classification = classify_ids(all_ids, gene_index, secondary_ids)

            result["matched_count"] = len(classification["matched"])
            result["matched_secondary_count"] = len(classification["matched_secondary"])
            result["rna_count"] = len(classification["rna_skipped"])
            result["mismatched_count"] = len(classification["mismatched"])
            result["empty_count"] = classification["empty"]
            result["mismatched_ids"] = classification["mismatched"]
            result["sample_ids"] = classification["matched"][:3] + classification["mismatched"][:3]

            # Cross-reference with import report
            ir_missing = find_import_report_mismatches(analysis, import_report)
            if ir_missing is not None:
                result["import_report_missing"] = ir_missing
                result["import_report_count"] = len(ir_missing)

            total_gene_ids = result["matched_count"] + result["matched_secondary_count"] + result["mismatched_count"]
            if total_gene_ids == 0:
                if result["rna_count"] > 0:
                    result["status"] = "RNA_ONLY"
                    result["details"] = f"All {result['rna_count']} IDs are RNA/tRNA/ncRNA (expected, not loaded as gene nodes)"
                else:
                    result["status"] = "ERROR"
                    result["details"] = "No valid gene IDs found in name_col"
                analyses_results.append(result)
                continue

            total_matched = result["matched_count"] + result["matched_secondary_count"]
            match_rate = total_matched / total_gene_ids if total_gene_ids > 0 else 0

            if match_rate >= 0.95:
                result["status"] = "MATCH"
                parts = [f"{result['matched_count']}/{total_gene_ids} IDs match gene nodes"]
                if result["matched_secondary_count"] > 0:
                    parts.append(f"+{result['matched_secondary_count']} via secondary column '{secondary_name_col}'")
                if result["rna_count"] > 0:
                    parts.append(f"{result['rna_count']} RNA/tRNA/ncRNA skipped (OK)")
                if result["mismatched_count"] > 0:
                    parts.append(f"{result['mismatched_count']} unmatched")
                result["details"] = "; ".join(parts)
                analyses_results.append(result)
                continue

            if match_rate >= 0.5:
                result["status"] = "PARTIAL MATCH"
            else:
                result["status"] = "NO MATCH"

            # --- IDs don't match well. Determine why and suggest fix. ---

            mismatched = classification["mismatched"]
            sample_mismatched = mismatched[:LOOKUP_SAMPLE_SIZE]

            # Check if organism genome is loaded
            if not is_organism_loaded(organism):
                result["fix_strategy"] = "LOAD_ORGANISM"
                result["details"] = (
                    f"Organism '{organism}' genome not loaded. "
                    f"{result['matched_count']}/{total_gene_ids} match, "
                    f"{result['mismatched_count']} don't. "
                    f"Add to cyanobacteria_genomes.csv to create gene nodes."
                )
                analyses_results.append(result)
                continue

            # Check if another column in the CSV has matching IDs
            alt_cols = check_other_csv_columns(
                str(csv_path), skip_rows, gene_index, name_col
            )
            if alt_cols:
                best_col, sample_val, cnt, total = alt_cols[0]
                result["fix_strategy"] = "CHANGE_NAME_COL"
                result["details"] = (
                    f"{result['matched_count']}/{total_gene_ids} match, "
                    f"{result['mismatched_count']} don't. "
                    f"Column '{best_col}' has IDs like '{sample_val}' "
                    f"that match gene nodes ({cnt}/{total}). "
                    f"Suggest: name_col: \"{best_col}\""
                )
                analyses_results.append(result)
                continue

            # Check gene_mapping.csv for the organism
            genome_dir = get_genome_dir(organism, project_root)
            if genome_dir:
                mapping_matches = check_gene_mapping(sample_mismatched, genome_dir)
                if mapping_matches:
                    best_col = max(mapping_matches, key=lambda k: mapping_matches[k][0])
                    cnt, (sample_id, mapped_to) = mapping_matches[best_col]
                    mapping_file = Path(genome_dir) / "gene_mapping.csv"
                    rel_mapping = os.path.relpath(mapping_file, project_root)
                    result["fix_strategy"] = "CREATE_MAPPING_CSV"
                    result["details"] = (
                        f"{result['matched_count']}/{total_gene_ids} match, "
                        f"{result['mismatched_count']} don't. "
                        f"Mismatched IDs match column '{best_col}' in gene_mapping.csv "
                        f"({cnt}/{len(sample_mismatched)} of sample match). "
                        f"E.g., '{sample_id}' -> locus_tag '{mapped_to}'. "
                        f"Run /fix-gene-ids on this paperconfig to create a "
                        f"_with_locus_tag.csv, then update name_col to 'locus_tag'."
                    )
                    analyses_results.append(result)
                    continue

                # gene_mapping.csv exists but IDs not found there — check GFF/GBK
                gff_files = list(Path(genome_dir).rglob("*.gff"))
                found_gff = False
                for gff_path in gff_files:
                    gff_matches = scan_gff_for_ids(str(gff_path), sample_mismatched[:10])
                    if gff_matches:
                        best_attr = max(gff_matches, key=gff_matches.get)
                        cnt = gff_matches[best_attr]
                        rel_mapping = os.path.relpath(
                            Path(genome_dir) / "gene_mapping.csv", project_root
                        )
                        result["fix_strategy"] = "CREATE_MAPPING_GFF"
                        result["details"] = (
                            f"{result['matched_count']}/{total_gene_ids} match, "
                            f"{result['mismatched_count']} don't. "
                            f"IDs found in GFF attribute '{best_attr}' "
                            f"in {gff_path.name} ({cnt} matches) "
                            f"but NOT in gene_mapping.csv. "
                            f"Add '{best_attr}' column to {rel_mapping}."
                        )
                        analyses_results.append(result)
                        found_gff = True
                        break

                if not found_gff:
                    result["fix_strategy"] = "UNRELATED_IDS"
                    result["details"] = (
                        f"{result['matched_count']}/{total_gene_ids} match, "
                        f"{result['mismatched_count']} don't. "
                        f"IDs not found in gene_mapping.csv or GFF/GBK files for {organism}. "
                        f"Manual investigation needed."
                    )
                    analyses_results.append(result)
            else:
                # No genome dir — fall back to gene node field check
                alt_fields = check_alt_gene_fields(sample_mismatched, gene_index)
                if alt_fields:
                    best_field = max(alt_fields, key=lambda k: alt_fields[k][0])
                    cnt, (sample_id, mapped_to) = alt_fields[best_field]
                    result["fix_strategy"] = "CREATE_MAPPING_CSV"
                    result["details"] = (
                        f"{result['matched_count']}/{total_gene_ids} match, "
                        f"{result['mismatched_count']} don't. "
                        f"IDs match gene node field '{best_field}' "
                        f"({cnt}/{len(sample_mismatched)} match). "
                        f"E.g., '{sample_id}' -> {mapped_to}. "
                        f"No gene_mapping.csv found — create one for this organism."
                    )
                else:
                    result["fix_strategy"] = "UNRELATED_IDS"
                    result["details"] = (
                        f"{result['matched_count']}/{total_gene_ids} match, "
                        f"{result['mismatched_count']} don't. "
                        f"IDs don't match gene nodes and no genome directory found for {organism}. "
                        f"Manual investigation needed."
                    )
                analyses_results.append(result)

    return {"papername": papername, "paperconfig_path": paperconfig_path, "analyses": analyses_results}


def format_report(all_results, import_report=None):
    """Format the full report as human-readable text."""
    lines = []
    lines.append("=" * 80)
    lines.append("Gene ID Validation Report")
    lines.append("=" * 80)
    if import_report:
        total_ir = len(import_report["all_missing"])
        lines.append(f"Import report loaded: {total_ir} unique missing gene IDs")
    lines.append("")

    summary_rows = []

    for paper in all_results:
        lines.append(f"## {paper['papername']}")
        lines.append(f"   Config: {paper['paperconfig_path']}")
        lines.append("")

        for a in paper["analyses"]:
            lines.append(f"  Table: {a['table']} (analysis: {a['analysis_id']})")
            lines.append(f"    Organism: {a['organism']} | name_col: {a['name_col']}")
            if a.get("secondary_name_col"):
                lines.append(f"    Secondary name_col: {a['secondary_name_col']}")
            lines.append(f"    Status: {a['status']}")

            # Detailed breakdown
            if a["total_ids"] > 0:
                total_gene = a["matched_count"] + a.get("matched_secondary_count", 0) + a["mismatched_count"]
                parts = []
                parts.append(f"Total rows: {a['total_ids']}")
                if a["matched_count"]:
                    parts.append(f"Matched: {a['matched_count']}")
                if a.get("matched_secondary_count", 0):
                    parts.append(f"Matched (secondary): {a['matched_secondary_count']}")
                if a["rna_count"]:
                    parts.append(f"RNA/tRNA/ncRNA: {a['rna_count']} (OK)")
                if a["mismatched_count"]:
                    parts.append(f"Mismatched: {a['mismatched_count']}")
                if a["empty_count"]:
                    parts.append(f"Empty/NA: {a['empty_count']}")
                lines.append(f"    Breakdown: {' | '.join(parts)}")

            # Import report cross-reference
            if a.get("import_report_count", 0) > 0:
                lines.append(f"    Import report: {a['import_report_count']} missing gene IDs confirmed in Neo4j import")

            if a["fix_strategy"]:
                lines.append(f"    Fix strategy: {a['fix_strategy']} — {a['details']}")
            elif a["details"]:
                lines.append(f"    {a['details']}")

            # Show mismatched IDs (all of them, grouped)
            mismatched = a.get("mismatched_ids", [])
            if mismatched:
                # Categorize mismatched IDs by pattern
                categories = defaultdict(list)
                for mid in mismatched:
                    if mid.endswith("*"):
                        categories["trailing_asterisk"].append(mid)
                    elif mid.endswith(" ") or mid != mid.strip():
                        categories["trailing_whitespace"].append(mid)
                    elif "_pseudo" in mid:
                        categories["pseudo_suffix"].append(mid)
                    elif mid == "(unannotated)":
                        categories["unannotated"].append(mid)
                    elif "," in mid:
                        categories["composite_ids"].append(mid)
                    elif mid.startswith("contig"):
                        categories["contig_ids"].append(mid)
                    else:
                        categories["other"].append(mid)

                lines.append(f"    Mismatched IDs ({len(mismatched)} total):")
                for cat, ids in sorted(categories.items()):
                    if len(ids) <= 10:
                        lines.append(f"      {cat} ({len(ids)}): {', '.join(ids)}")
                    else:
                        shown = ', '.join(ids[:10])
                        lines.append(f"      {cat} ({len(ids)}): {shown} ... +{len(ids)-10} more")

            lines.append("")

            summary_rows.append({
                "publication": paper["papername"],
                "table": a["table"],
                "analysis_id": a["analysis_id"],
                "status": a["status"],
                "matched": a["matched_count"],
                "matched_sec": a.get("matched_secondary_count", 0),
                "rna": a["rna_count"],
                "mismatched": a["mismatched_count"],
                "import_report": a.get("import_report_count", 0),
                "fix_strategy": a["fix_strategy"] or "—",
            })

    # Summary table
    lines.append("")
    lines.append("=" * 80)
    lines.append("Summary")
    lines.append("=" * 80)
    lines.append("")

    if not summary_rows:
        lines.append("No analyses found.")
        return "\n".join(lines)

    # Column widths
    pub_w = max(len("Publication"), max(len(r["publication"]) for r in summary_rows))
    aid_w = max(len("Analysis"), max(min(len(r["analysis_id"]), 35) for r in summary_rows))
    stat_w = max(len("Status"), max(len(r["status"]) for r in summary_rows))

    header = (
        f"{'Publication':<{pub_w}} | {'Analysis':<{aid_w}} | {'Status':<{stat_w}} | "
        f"{'Match':>5} | {'2nd':>3} | {'RNA':>3} | {'Miss':>4} | {'IR':>3} | Fix Strategy"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for r in summary_rows:
        aid_short = r["analysis_id"][:35]
        ir_str = str(r["import_report"]) if r["import_report"] > 0 else ""
        sec_str = str(r["matched_sec"]) if r["matched_sec"] > 0 else ""
        lines.append(
            f"{r['publication']:<{pub_w}} | {aid_short:<{aid_w}} | {r['status']:<{stat_w}} | "
            f"{r['matched']:>5} | {sec_str:>3} | {r['rna']:>3} | {r['mismatched']:>4} | {ir_str:>3} | {r['fix_strategy']}"
        )

    # Stats
    total = len(summary_rows)
    matched = sum(1 for r in summary_rows if r["status"] == "MATCH")
    partial = sum(1 for r in summary_rows if r["status"] == "PARTIAL MATCH")
    no_match = sum(1 for r in summary_rows if r["status"] == "NO MATCH")
    rna_only = sum(1 for r in summary_rows if r["status"] == "RNA_ONLY")
    errors = sum(1 for r in summary_rows if r["status"] == "ERROR")

    lines.append("")
    lines.append(f"Total analyses: {total}")
    lines.append(f"  MATCH: {matched}")
    if partial:
        lines.append(f"  PARTIAL MATCH: {partial}")
    if no_match:
        lines.append(f"  NO MATCH: {no_match}")
    if rna_only:
        lines.append(f"  RNA_ONLY: {rna_only}")
    if errors:
        lines.append(f"  ERROR: {errors}")

    # Total counts
    total_matched = sum(r["matched"] + r["matched_sec"] for r in summary_rows)
    total_rna = sum(r["rna"] for r in summary_rows)
    total_mismatched = sum(r["mismatched"] for r in summary_rows)
    total_ir = sum(r["import_report"] for r in summary_rows)
    lines.append("")
    lines.append(f"Total gene IDs across all analyses: {total_matched + total_mismatched}")
    lines.append(f"  Matched: {total_matched}")
    lines.append(f"  RNA/tRNA/ncRNA (OK): {total_rna}")
    lines.append(f"  Mismatched: {total_mismatched}")
    if total_ir:
        lines.append(f"  Confirmed in import report: {total_ir}")

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
    parser.add_argument(
        "--import-report",
        help="Path to import.report file (default: tries Docker, then output/import.report)",
    )
    parser.add_argument(
        "--no-import-report",
        action="store_true",
        help="Skip loading import report",
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

    # Load import report
    import_report = None
    if not args.no_import_report:
        import_report = load_import_report(args.import_report)
        if import_report:
            print(f"Import report: {len(import_report['all_missing'])} unique missing gene IDs", file=sys.stderr)
        else:
            print("No import report available (use --import-report or run Docker)", file=sys.stderr)

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
        result = analyze_paperconfig(pc_path, gene_index, project_root, import_report)
        all_results.append(result)

    # Print report to stdout and write to file
    report = format_report(all_results, import_report)
    print(report)

    report_file = project_root / "gene_id_validation_report.txt"
    with open(report_file, "w") as f:
        f.write(report)
        f.write("\n")
    print(f"\nReport written to {report_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
