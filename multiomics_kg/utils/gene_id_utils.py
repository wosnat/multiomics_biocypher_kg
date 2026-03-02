"""Shared utilities for gene ID mapping.

Used by the check-gene-ids, fix-gene-ids, and build-gene-mapping-supp skills.
Primary data source: gene_annotations_merged.json per strain (falls back to _wide.json).
"""

import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd


# ─── Organism → genome cache directory ───────────────────────────────────────

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

# ─── Patterns & keywords ──────────────────────────────────────────────────────

# RNA features (tRNA, ncRNA, rRNA — intentionally not in gene nodes)
SKIP_PATTERNS = re.compile(
    r'^(tRNA|ncRNA|rRNA|Yfr\d|tmRNA|RNA_\d|PMT_ncRNA_|'
    r'\w+_tRNA\w+VIMSS|'
    r'\w+_rr[ls]VIMSS)',
    re.IGNORECASE,
)

# Fields in gene_annotations_merged.json that hold alternative gene identifiers
ANNOTATION_ID_FIELDS = [
    "locus_tag",               # primary key (identity mapping)
    "locus_tag_ncbi",          # RS-format NCBI locus tag (e.g. A9601_RS09110)
    "locus_tag_cyanorak",      # Cyanorak ID (e.g. CK_Pro_AS9601_00001)
    "protein_id",              # RefSeq WP_ accession
    "gene_name",               # primary gene symbol (single string)
    "gene_synonyms",           # list of alternative names (mixed: names + locus tags)
    "gene_name_synonyms",      # list of gene-name-like tokens only
    "alternative_locus_tags",  # list of locus-tag-like tokens only
    "old_locus_tags",          # list of old locus tag strings
]

# Column-detection heuristics (shared by fix-gene-ids and build-gene-mapping-supp)
ID_COL_KEYWORDS = [
    'locus', 'tag', 'gene', 'protein', 'orf', 'id', 'accession', 'homolog',
]
DESCRIPTION_COL_KEYWORDS = [
    'product', 'description', 'function', 'role', 'cog', 'note', 'go',
    'ontology', 'kegg', 'pathway', 'domain', 'tigrfam', 'pfam', 'eggnog',
    'cluster', 'inference', 'species', 'organism', 'strain', 'taxon',
    'annotation', 'regulation', 'temperature',
]
MAX_UNIQUE_FOR_ID_COL = 10


# ─── Genome directory resolution ──────────────────────────────────────────────

def get_genome_dir(organism_name, project_root):
    """Return the absolute genome cache directory for organism_name, or None."""
    if not organism_name:
        return None
    norm = organism_name.strip().strip('"').lower()
    for key, rel_path in ORGANISM_TO_GENOME_DIR.items():
        if key in norm or norm in key:
            full = Path(project_root) / rel_path
            if full.exists():
                return str(full)
    return None


def is_organism_loaded(organism_name):
    """Return True if the organism has a genome directory configured."""
    if not organism_name:
        return False
    norm = organism_name.strip().strip('"').lower()
    if norm in ORGANISM_TO_GENOME_DIR:
        return True
    for known in ORGANISM_TO_GENOME_DIR:
        if known in norm or norm in known:
            return True
    return False


# ─── JSON annotation loading ──────────────────────────────────────────────────

def load_gene_annotations(genome_dir):
    """Load gene_annotations_merged.json (falls back to _wide.json).

    Returns dict keyed by locus_tag, or None if neither file exists.
    """
    for filename in ("gene_annotations_merged.json", "gene_annotations_wide.json"):
        path = Path(genome_dir) / filename
        if path.exists():
            try:
                with open(path) as f:
                    return json.load(f)
            except Exception as e:
                print(f"  Error reading {path}: {e}")
    return None


# ─── gene_id_mapping.json loading ────────────────────────────────────────────


def load_gene_id_mapping(genome_dir):
    """Load gene_id_mapping.json. Returns None if the file does not exist.

    Falls back gracefully: callers should use load_gene_annotations() when None
    is returned.
    """
    path = Path(genome_dir) / "gene_id_mapping.json"
    if path.exists():
        try:
            with open(path) as f:
                return json.load(f)
        except Exception as e:
            print(f"  Error reading {path}: {e}")
    return None


def build_id_lookup_from_mapping(gene_id_mapping):
    """Build alt_id → locus_tag reverse index from gene_id_mapping.json.

    Only includes organism_specific and annotation_specific alt-IDs so the
    lookup remains strictly 1:1 (one ID → one gene).  Generic-scoped IDs
    (gene_name, gene_synonym) are intentionally excluded — they are shared
    across organisms and are accessed via search_genes_by_name() instead.

    Paper-derived IDs that are not already in the reference set are tracked in
    paper_id_keys (the returned third element, analogous to supp_keys from the
    original build_id_lookup).

    Returns:
        (lookup_dict, locus_tag_set, paper_id_keys)
        lookup_dict: raw_id -> locus_tag  (organism_specific / annotation_specific only)
        locus_tag_set: set of canonical locus tags
        paper_id_keys: set of IDs that came from paper_ids (not reference)
    """
    lookup: dict[str, str] = {}
    locus_tags: set[str] = set()
    paper_id_keys: set[str] = set()

    for locus_tag, entry in gene_id_mapping.items():
        locus_tags.add(locus_tag)
        lookup[locus_tag] = locus_tag

        for alt in (entry.get("alt_ids") or {}).get("reference") or []:
            if alt.get("scope") == "generic":
                continue
            aid = (alt.get("id") or "").strip()
            if aid:
                if aid in lookup and lookup[aid] != locus_tag:
                    print(f"Warning: ID '{aid}' maps to multiple locus tags: {lookup[aid]}, {locus_tag}")
                lookup[aid] = locus_tag

        for alt in (entry.get("alt_ids") or {}).get("paper_ids") or []:
            if alt.get("scope") == "generic":
                continue
            aid = (alt.get("id") or "").strip()
            if aid:
                if aid in lookup and lookup[aid] != locus_tag:
                    print(f"Warning: paper ID '{aid}' maps to multiple locus tags: {lookup[aid]}, {locus_tag}")
                if aid not in lookup:
                    lookup[aid] = locus_tag
                    paper_id_keys.add(aid)

    return lookup, locus_tags, paper_id_keys


# ─── ID lookup building ───────────────────────────────────────────────────────

def build_id_lookup(genome_dir):
    """Build a unified raw_id → locus_tag reverse index.

    Prefers gene_id_mapping.json (built by build_gene_id_mapping.py) when
    available — it includes all reference alt-IDs plus paper-derived IDs.
    Falls back to gene_annotations_merged.json + gene_mapping_supp.csv when
    the mapping file has not been built yet.

    Returns:
        (lookup_dict, locus_tag_set, supp_keys)
        lookup_dict: raw_id -> locus_tag
        locus_tag_set: set of canonical locus tags
        supp_keys: set of IDs from supplementary/paper sources (not canonical annotations)
        Returns (None, None, set()) if no annotation data found.
    """
    gene_id_mapping = load_gene_id_mapping(genome_dir)
    if gene_id_mapping is not None:
        return build_id_lookup_from_mapping(gene_id_mapping)

    annotations = load_gene_annotations(genome_dir)
    if annotations is None:
        return None, None, set()

    lookup = {}
    locus_tags = set()

    for locus_tag, entry in annotations.items():
        locus_tags.add(locus_tag)
        lookup[locus_tag] = locus_tag  # identity mapping

        # Single-value string fields
        for field in ("locus_tag_ncbi", "locus_tag_cyanorak", "protein_id", "gene_name"):
            val = entry.get(field)
            if val and isinstance(val, str):
                val = val.strip()
                if val:
                    lookup[val] = locus_tag

        # List fields (already Python lists in the JSON)
        for field in ("gene_synonyms", "gene_name_synonyms", "alternative_locus_tags", "old_locus_tags"):
            vals = entry.get(field)
            if isinstance(vals, list):
                for v in vals:
                    if isinstance(v, str):
                        v = v.strip()
                        if v:
                            lookup[v] = locus_tag
            elif isinstance(vals, str) and vals.strip():
                # Fallback: comma-delimited string (shouldn't happen with merged JSON)
                for v in re.split(r',\s*', vals):
                    v = v.strip()
                    if v:
                        lookup[v] = locus_tag

    # Merge supplementary mappings as fallback
    supp_lookup = load_supp_lookup(genome_dir)
    supp_keys = set()
    for alt_id, lt in supp_lookup.items():
        if alt_id not in lookup:
            lookup[alt_id] = lt
            supp_keys.add(alt_id)

    return lookup, locus_tags, supp_keys


def build_field_lookups(genome_dir):
    """Build per-field lookups for fix-strategy detection in check-gene-ids.

    Returns dict: field_name -> {id_value -> locus_tag}, or None if JSON not found.
    Used to determine *which* field the mismatched IDs match, so the right
    fix strategy can be reported.
    """
    annotations = load_gene_annotations(genome_dir)
    if annotations is None:
        return None

    lookups = {}

    for locus_tag, entry in annotations.items():
        # Single-value string fields
        for field in ("locus_tag_ncbi", "locus_tag_cyanorak", "protein_id", "gene_name"):
            val = entry.get(field)
            if val and isinstance(val, str):
                val = val.strip()
                if val:
                    lookups.setdefault(field, {})[val] = locus_tag

        # List fields
        for field in ("gene_synonyms", "gene_name_synonyms", "alternative_locus_tags", "old_locus_tags"):
            vals = entry.get(field)
            if isinstance(vals, list):
                for v in vals:
                    if isinstance(v, str):
                        v = v.strip()
                        if v:
                            lookups.setdefault(field, {})[v] = locus_tag

    return lookups if lookups else None


def load_supp_lookup(genome_dir):
    """Build alt_id → locus_tag dict from gene_mapping_supp.csv.

    Returns empty dict if file not found.
    """
    supp_file = Path(genome_dir) / "gene_mapping_supp.csv"
    if not supp_file.exists():
        return {}
    try:
        df = pd.read_csv(supp_file)
    except Exception:
        return {}

    lookup = {}
    for _, row in df.iterrows():
        locus = str(row["locus_tag"]).strip() if pd.notna(row.get("locus_tag")) else ""
        alt_id = str(row["alt_id"]).strip() if pd.notna(row.get("alt_id")) else ""
        if locus and locus != "nan" and alt_id and alt_id != "nan":
            if alt_id not in lookup:
                lookup[alt_id] = locus
    return lookup


# ─── Gene ID mapping ──────────────────────────────────────────────────────────

def map_gene_id(raw_id, lookup, locus_tags, supp_keys=None):
    """Map a single gene ID string to a canonical locus_tag.

    Tries in order: direct (already a locus_tag), lookup, zero-padding,
    comma-split composite.

    Returns:
        (locus_tag, method) where method is one of:
        'direct', 'lookup', 'supp', 'repadded', 'composite_direct', 'composite_lookup'
        Returns (None, None) if no mapping found.
    """
    if supp_keys is None:
        supp_keys = set()
    raw_id = raw_id.strip()

    if raw_id in locus_tags:
        return raw_id, "direct"

    if raw_id in lookup:
        method = "supp" if raw_id in supp_keys else "lookup"
        return lookup[raw_id], method

    # Zero-padding normalization (e.g. MIT1002_0001 → MIT1002_00001)
    m = re.match(r'^(.+)_(\d+)$', raw_id)
    if m:
        prefix, num = m.group(1), m.group(2)
        for pad in range(len(num) + 1, len(num) + 3):
            padded = f"{prefix}_{num.zfill(pad)}"
            if padded in locus_tags:
                return padded, "repadded"
            if padded in lookup:
                return lookup[padded], "repadded"

    # Comma-split composite IDs (e.g. "rps13,rpsM")
    if "," in raw_id:
        for part in [p.strip() for p in raw_id.split(",")]:
            if part in locus_tags:
                return part, "composite_direct"
            if part in lookup:
                return lookup[part], "composite_lookup"

    return None, None


# ─── Column detection ─────────────────────────────────────────────────────────

def is_id_like_column(series, col_name, exclude_cols):
    """Return True if a DataFrame column looks like it contains gene identifiers.

    Applies heuristics:
    - Excludes columns by name keywords (descriptions, coordinates)
    - Excludes purely numeric, low-cardinality categorical, or long-text columns
    - Includes columns with ID-keyword names or short non-space values
    """
    if col_name in exclude_cols:
        return False

    col_lower = col_name.lower()

    if any(kw in col_lower for kw in DESCRIPTION_COL_KEYWORDS):
        return False

    if any(kw in col_lower for kw in ('start', 'end', 'length', 'position', 'coord')):
        return False

    vals = series.dropna().astype(str)
    vals = vals[vals.str.strip() != '']
    vals = vals[vals != 'nan']
    if len(vals) == 0:
        return False

    # Low-cardinality categorical (e.g. Up/Down)
    n_unique = vals.nunique()
    if n_unique <= MAX_UNIQUE_FOR_ID_COL and len(vals) / max(n_unique, 1) >= 5:
        return False

    # Purely numeric
    try:
        pd.to_numeric(vals)
        return False
    except (ValueError, TypeError):
        pass

    # Very long values (descriptions)
    if vals.str.len().mean() > 50:
        return False

    # Float-like fold-changes (e.g. "-0.54", "1.8*")
    numeric_like = vals.str.replace(r'[*\s]', '', regex=True)
    try:
        pd.to_numeric(numeric_like)
        return False
    except (ValueError, TypeError):
        pass

    if any(kw in col_lower for kw in ID_COL_KEYWORDS):
        return True

    # No keyword: require most values to lack spaces
    if vals.str.contains(' ').mean() > 0.3:
        return False

    return True


# ─── Full-text gene search ────────────────────────────────────────────────────


def search_genes_by_name(genome_dir, query, fields=None):
    """Full-text search across gene names, synonyms, and product descriptions.

    Args:
        genome_dir: path to genome cache directory (absolute or relative to cwd)
        query: search string (case-insensitive substring match)
        fields: optional list of field types to restrict search; defaults to all.
                Recognised values: "gene_name", "gene_synonym", "product"

    Returns list of {"locus_tag", "match_field", "match_value"} dicts, deduplicated.
    Draws from gene_id_mapping.json (when present) for gene names / paper synonyms
    and from gene_annotations_merged.json for canonical gene_name and product fields.
    """
    query_lower = query.strip().lower()
    if not query_lower:
        return []

    search_names = fields is None or "gene_name" in fields
    search_synonyms = fields is None or "gene_synonym" in fields
    search_products = fields is None or "product" in fields

    results: list[dict] = []
    seen: set[tuple] = set()

    def _add(locus_tag: str, field: str, value: str) -> None:
        key = (locus_tag, field, value)
        if key not in seen:
            seen.add(key)
            results.append({"locus_tag": locus_tag, "match_field": field, "match_value": value})

    # Search gene_id_mapping.json (gene names from all sources + paper product synonyms)
    gene_id_mapping = load_gene_id_mapping(genome_dir)
    if gene_id_mapping:
        for locus_tag, entry in gene_id_mapping.items():
            for section in ("reference", "paper_ids"):
                for alt in (entry.get("alt_ids") or {}).get(section) or []:
                    id_type = alt.get("id_type", "")
                    val = (alt.get("id") or "").strip()
                    if not val:
                        continue
                    if search_names and id_type == "gene_name" and query_lower in val.lower():
                        _add(locus_tag, "gene_name", val)
                    elif search_synonyms and id_type == "gene_synonym" and query_lower in val.lower():
                        _add(locus_tag, "gene_synonym", val)

            if search_products:
                for syn in entry.get("product_synonyms") or []:
                    text = (syn.get("text") or "").strip()
                    if text and query_lower in text.lower():
                        _add(locus_tag, "product_synonym", text)

    # Also search gene_annotations_merged.json for canonical gene_name and product fields
    annotations = load_gene_annotations(genome_dir)
    if annotations:
        for locus_tag, entry in annotations.items():
            if search_names:
                gn = (entry.get("gene_name") or "").strip()
                if gn and query_lower in gn.lower():
                    _add(locus_tag, "gene_name", gn)
                for syn in entry.get("gene_synonyms") or []:
                    if isinstance(syn, str) and syn.strip() and query_lower in syn.lower():
                        _add(locus_tag, "gene_synonym", syn)

            if search_products:
                for field in ("product", "product_cyanorak"):
                    val = (entry.get(field) or "").strip()
                    if val and query_lower in val.lower():
                        _add(locus_tag, field, val)

    return results


# ─── Import report loading ────────────────────────────────────────────────────

def load_import_report(import_report_path=None, return_by_source=False):
    """Load missing gene IDs from the Neo4j import report.

    Tries (in order): explicit path, output/import.report, Docker exec.

    Args:
        import_report_path: explicit path to import.report, or None for auto-detect
        return_by_source: if True, return {"by_source": dict, "all_missing": set}
                          if False, return set of missing IDs (or None if unavailable)

    Returns None if no report could be loaded.
    """
    report_lines = []

    if import_report_path and Path(import_report_path).exists():
        with open(import_report_path) as f:
            report_lines = f.readlines()
    else:
        local_report = Path("output/import.report")
        if local_report.exists():
            with open(local_report) as f:
                report_lines = f.readlines()
        else:
            try:
                result = subprocess.run(
                    ["docker", "compose", "exec", "deploy", "cat", "/data/build2neo/import.report"],
                    capture_output=True, text=True, timeout=15,
                )
                if result.returncode == 0 and result.stdout.strip():
                    report_lines = result.stdout.strip().split("\n")
            except Exception:
                pass

    if not report_lines:
        return None

    affects_pattern = re.compile(
        r'^(.+?) \(global id space\)-\[Affects_expression_of\]->ncbigene:(.+?) \(global id space\)'
    )

    missing_by_source = defaultdict(set)
    all_missing = set()

    for line in report_lines:
        m = affects_pattern.match(line.strip())
        if m:
            source_id, missing_gene = m.group(1), m.group(2)
            missing_by_source[source_id].add(missing_gene)
            all_missing.add(missing_gene)

    if not all_missing:
        return None

    if return_by_source:
        return {"by_source": missing_by_source, "all_missing": all_missing}
    return all_missing
