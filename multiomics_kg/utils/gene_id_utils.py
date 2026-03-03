"""Shared utilities for gene ID mapping.

Used by the check-gene-ids, fix-gene-ids, and build-gene-mapping-supp skills.
Primary data source: gene_id_mapping.json v2 per strain (falls back to
gene_annotations_merged.json when the mapping has not been built yet).

v2 API
------
load_mapping_v2(genome_dir) -> MappingData
resolve_row(row, name_col, id_columns, md) -> (locus_tag | None, method_str)
expand_list(raw_val) -> list[str]

Legacy API (still supported for check-gene-ids / fix-gene-ids)
--------------------------------------------------------------
build_id_lookup(genome_dir) -> (lookup, locus_tags, supp_keys)
map_gene_id(raw_id, lookup, locus_tags, supp_keys) -> (locus_tag, method)
"""

import json
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

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


# ─── v2 API ───────────────────────────────────────────────────────────────────


@dataclass
class MappingData:
    """Gene ID mapping for one organism/strain (gene_id_mapping.json v2).

    Attributes:
        specific_lookup: Tier 1 alt-ID → locus_tag (1:1)
        multi_lookup: Tier 2+3 alt-ID → [locus_tag, ...] (1:many)
        conflicts: Tier 1 alt-ID → [lt1, lt2] (data errors)
        locus_tags: set of canonical locus_tags
        version: schema version (2 for new format, 1 for legacy)
    """
    specific_lookup: dict[str, str] = field(default_factory=dict)
    multi_lookup: dict[str, list[str]] = field(default_factory=dict)
    conflicts: dict[str, list[str]] = field(default_factory=dict)
    locus_tags: set[str] = field(default_factory=set)
    version: int = 2


def load_mapping_v2(genome_dir: str | Path) -> Optional["MappingData"]:
    """Load gene_id_mapping.json and return a MappingData instance.

    Handles both v2 (new graph-based) and v1 (legacy) formats.
    Returns None if no mapping file or annotation data exists.
    """
    raw = load_gene_id_mapping(genome_dir)
    if raw is None:
        return None

    version = raw.get("version", 1)

    if version == 2:
        # v2 format: has specific_lookup and multi_lookup pre-built
        md = MappingData(
            specific_lookup=raw.get("specific_lookup", {}),
            multi_lookup=raw.get("multi_lookup", {}),
            conflicts=raw.get("conflicts", {}),
            version=2,
        )
        # locus_tags: all keys in genes dict + keys that self-map in specific_lookup
        genes_dict = raw.get("genes", {})
        md.locus_tags = set(genes_dict.keys())
        # Also include anything in specific_lookup that maps to itself (canonical forms)
        for k, v in md.specific_lookup.items():
            md.locus_tags.add(v)
        return md

    # v1 format: legacy structure with alt_ids.reference / alt_ids.paper_ids
    # Convert to MappingData for uniform access
    md = MappingData(version=1)
    for locus_tag, entry in raw.items():
        md.locus_tags.add(locus_tag)
        md.specific_lookup[locus_tag] = locus_tag

        for alt in (entry.get("alt_ids") or {}).get("reference") or []:
            if alt.get("scope") == "generic":
                continue
            aid = (alt.get("id") or "").strip()
            if aid and aid not in md.specific_lookup:
                md.specific_lookup[aid] = locus_tag

        for alt in (entry.get("alt_ids") or {}).get("paper_ids") or []:
            if alt.get("scope") == "generic":
                continue
            aid = (alt.get("id") or "").strip()
            if aid and aid not in md.specific_lookup:
                md.specific_lookup[aid] = locus_tag

    return md


def expand_list(raw_val: str) -> list[str]:
    """Split a potentially list-valued cell into candidate ID strings.

    Always includes the full raw value first (in case the separator is part
    of the ID), then individual tokens split on ',' and ';'.

    Examples:
        "PMM0001" → ["PMM0001"]
        "PMM0001, PMM0002" → ["PMM0001, PMM0002", "PMM0001", "PMM0002"]
        "dnaA; dnaN" → ["dnaA; dnaN", "dnaA", "dnaN"]
    """
    raw_val = str(raw_val).strip()
    if not raw_val or raw_val.lower() in ("nan", ""):
        return []
    candidates: list[str] = [raw_val]
    # Split on comma and/or semicolon
    if "," in raw_val or ";" in raw_val:
        for part in re.split(r"[,;]", raw_val):
            part = part.strip()
            if part and part not in candidates:
                candidates.append(part)
    return candidates


def _heuristic_candidates(raw_val: str) -> list[str]:
    """Return heuristic normalized forms of a raw ID (in addition to the raw form).

    Heuristics:
    - Strip trailing '*' (footnote artifact)
    - Strip trailing/leading whitespace (already done by caller)
    - Zero-pad numeric suffix (e.g. MIT1002_0001 → MIT1002_00001)
    """
    candidates: list[str] = []
    # Strip trailing asterisk
    stripped = raw_val.rstrip("*").strip()
    if stripped and stripped != raw_val:
        candidates.append(stripped)
    # Zero-pad numeric suffix
    m = re.match(r'^(.+)_(\d+)$', raw_val)
    if m:
        prefix, num = m.group(1), m.group(2)
        for pad in range(len(num) + 1, len(num) + 3):
            padded = f"{prefix}_{num.zfill(pad)}"
            if padded not in candidates:
                candidates.append(padded)
    return candidates


def resolve_row(
    row: dict,
    name_col: str,
    id_columns: list[dict],
    mapping_data: "MappingData",
) -> tuple[Optional[str], str]:
    """Resolve a CSV row to a canonical locus_tag using three-tier strategy.

    Tries columns in order: name_col first, then id_columns.
    Within each column, tries in order:
      Pass 1: specific_lookup (Tier 1) with list expansion
      Pass 2: heuristics (zero-pad, strip asterisk) → specific_lookup
      Pass 3: multi_lookup (Tier 2+3), accept only singletons

    Never duplicates rows. When a column is ambiguous (multi_lookup > 1 match),
    tries other columns before giving up.

    Args:
        row: dict of column_name → cell_value
        name_col: primary ID column name
        id_columns: list of {column, id_type} dicts from paperconfig
        mapping_data: MappingData instance for the organism

    Returns:
        (locus_tag, method_str) where method_str describes how it was resolved.
        locus_tag is None when resolution fails; method_str gives the reason.
    """
    sl = mapping_data.specific_lookup
    ml = mapping_data.multi_lookup
    conflicts = mapping_data.conflicts

    # Build ordered column list: name_col first, then id_columns
    all_cols: list[str] = [name_col]
    for c in id_columns:
        col = c.get("column", "")
        if col and col not in all_cols:
            all_cols.append(col)

    diagnostic: dict[str, str] = {}  # col → failure reason

    # ── Pass 1: specific_lookup with list expansion (+ direct locus_tag check) ──
    for col in all_cols:
        raw = row.get(col)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        for val in expand_list(str(raw)):
            if not val:
                continue
            if val in sl:
                return sl[val], f"tier1:{col}"
            if val in mapping_data.locus_tags:
                return val, f"locus_tag:{col}"
            if val in conflicts:
                diagnostic[col] = "tier1_conflict"

    # ── Pass 2: heuristics → specific_lookup ──────────────────────────────────
    for col in all_cols:
        raw = row.get(col)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        for val in expand_list(str(raw)):
            for h_val in _heuristic_candidates(val):
                if h_val in sl:
                    return sl[h_val], f"heuristic:{col}"

    # ── Pass 3: multi_lookup, singletons only ─────────────────────────────────
    for col in all_cols:
        raw = row.get(col)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        for val in expand_list(str(raw)):
            if not val:
                continue
            matches = ml.get(val)
            if matches:
                if len(matches) == 1:
                    return matches[0], f"multi:{col}"
                else:
                    diagnostic[col] = f"ambiguous:{len(matches)}"

    # ── All failed — report most informative reason ────────────────────────────
    if any(v == "tier1_conflict" for v in diagnostic.values()):
        reason = "tier1_conflict"
    elif any(v.startswith("ambiguous") for v in diagnostic.values()):
        reason = "ambiguous"
    else:
        # Check if any list-expanded values from name_col had multiple parts resolving differently
        raw_name = row.get(name_col)
        if raw_name is not None and not (isinstance(raw_name, float) and pd.isna(raw_name)):
            parts = expand_list(str(raw_name))
            if len(parts) > 1:
                resolved_lts = {sl.get(p) or (ml.get(p, [None]) if p in ml and len(ml.get(p, [])) == 1 else [None])[0]
                                for p in parts[1:] if p}
                resolved_lts.discard(None)
                if len(resolved_lts) > 1:
                    return None, "multi_value_ambiguous"
        reason = "unresolved"

    return None, reason


def build_id_lookup_from_mapping(gene_id_mapping):
    """Build alt_id → locus_tag reverse index from gene_id_mapping.json.

    Supports both v1 and v2 formats. Returns (lookup, locus_tags, paper_id_keys)
    for backward compatibility with check-gene-ids / fix-gene-ids code.

    For v2: returns specific_lookup directly; paper_id_keys is empty (deprecated distinction).
    For v1: original logic (organism_specific only, excludes generic IDs).

    Returns:
        (lookup_dict, locus_tag_set, paper_id_keys)
    """
    version = gene_id_mapping.get("version", 1) if isinstance(gene_id_mapping.get("version"), int) else 1

    if version == 2:
        # v2: specific_lookup already built; locus_tags from genes dict
        specific = gene_id_mapping.get("specific_lookup", {})
        locus_tags: set[str] = set(gene_id_mapping.get("genes", {}).keys())
        # Add self-mapped entries (locus_tag → locus_tag)
        for v in specific.values():
            locus_tags.add(v)
        # Reconstruct full lookup with self-mappings for locus_tags
        lookup = dict(specific)
        for lt in locus_tags:
            lookup[lt] = lt
        return lookup, locus_tags, set()

    # v1 legacy format
    lookup: dict[str, str] = {}
    locus_tags_v1: set[str] = set()
    paper_id_keys: set[str] = set()

    for locus_tag, entry in gene_id_mapping.items():
        if not isinstance(entry, dict):
            continue
        locus_tags_v1.add(locus_tag)
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

    return lookup, locus_tags_v1, paper_id_keys


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
