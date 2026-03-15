#!/usr/bin/env python3
"""
Build per-strain gene annotation tables by merging three sources:
  1. gene_mapping.csv      (NCBI + Cyanorak merged)
  2. eggnog annotations    (.emapper.annotations)
  3. uniprot JSON          (uniprot_preprocess_data.json)

Merge rules are defined in config/gene_annotations_config.yaml.

Outputs per strain:
  cache/data/<org>/genomes/<strain>/gene_annotations_wide.json
  cache/data/<org>/genomes/<strain>/gene_annotations_merged.json

Usage:
  uv run python multiomics_kg/download/build_gene_annotations.py [--strains STRAIN ...] [--force]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from typing import Any

# Identifier-style gene names that should NOT appear in gene_summary.
# Matches locus-tag patterns like ALTBGP6_RS00025, MIT1002_00123, PMM0001, SYNW1033.
_IDENTIFIER_RE = re.compile(
    r'^[A-Za-z]+\d*_(?:RS)?\d+$'   # PREFIX_[RS]DIGITS  (e.g. TK37_RS12345, A9601_12345)
    r'|^[A-Z]{3,5}\d{4,}$'         # 3-5 uppercase + 4+ digits  (e.g. SYNW1033, PMM0001)
)
from urllib.parse import unquote

from multiomics_kg.download.utils.annotation_helpers import (
    _coerce_to_tokens,
    _nonempty,
    _split,
    extract_first_match_in_sources,
)
from multiomics_kg.download.utils.annotation_transforms import (
    _TRANSFORMS,
    _tx_add_go_prefix,
    _tx_extract_go_from_pipe,
    _tx_extract_pfam_ids,
    _tx_extract_pfam_names,
    _tx_first_token_space,
    _tx_strip_function_prefix,
    _tx_strip_prefix_ko,
)
from multiomics_kg.download.utils.cli import add_common_args, load_config, load_genome_rows
from multiomics_kg.download.utils.paths import PROJECT_ROOT, infer_organism_group

# ─── gene_category mapping tables ─────────────────────────────────────────────
# Map functional role annotations to ~26 controlled gene_category values.
# Three independent classification systems are used in priority order:
#   1. Cyanorak Role (cyanobacteria only, highest specificity)
#   2. TIGR Role (cyanobacteria only)
#   3. COG category (universal, from eggNOG)
#
# WARNING: COG and Cyanorak use the SAME single-letter codes for DIFFERENT
# functions. E.g., COG "E" = Amino acid metabolism, Cyanorak "E" = Central
# intermediary metabolism (N/P/S). This is intentional — the two classification
# systems are independent.

COG_TO_CATEGORY = {
    "A": "Transcription",
    "B": "Replication and repair",
    "C": "Energy production",
    "D": "Cell cycle and division",
    "E": "Amino acid metabolism",       # ≠ Cyanorak E
    "F": "Nucleotide metabolism",
    "G": "Carbohydrate metabolism",
    "H": "Coenzyme metabolism",
    "I": "Lipid metabolism",
    "J": "Translation",
    "K": "Transcription",
    "L": "Replication and repair",
    "M": "Cell wall and membrane",
    "N": "Cell motility",
    "O": "Post-translational modification",
    "P": "Inorganic ion transport",
    "Q": "Secondary metabolites",
    "R": "Unknown",
    "S": "Unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
    "W": "Cell wall and membrane",
    "X": "Mobile elements",
    "Y": "Unknown",
    "Z": "Cell cycle and division",
}

# Cyanorak top-level letter → gene_category.
# "D" is handled separately via CYANORAK_D_SUBCODES (it's a catch-all).
CYANORAK_TO_CATEGORY = {
    "A": "Amino acid metabolism",
    "B": "Coenzyme metabolism",
    "C": "Cell wall and membrane",
    # "D" handled by CYANORAK_D_SUBCODES below
    "E": "Central intermediary metabolism",  # ≠ COG E! N, P, S metabolism
    "F": "Replication and repair",
    "G": "Carbohydrate metabolism",          # incl. glycolysis, TCA, CO2 fixation
    "H": "Lipid metabolism",
    "I": "Mobile elements",
    "J": "Photosynthesis",
    "K": "Translation",
    "L": "Post-translational modification",
    "M": "Nucleotide metabolism",
    "N": "Regulatory functions",
    "O": "Signal transduction",
    "P": "Transcription",
    "Q": "Transport",
    "R": "Unknown",
}

CYANORAK_D_SUBCODES = {
    "D.1": "Stress response and adaptation",
    "D.2": "Cell cycle and division",
    "D.3": "Cellular processes",
    "D.4": "Post-translational modification",
    "D.5": "Cell motility",
    "D.6": "Cellular processes",
    "D.7": "Cellular processes",
}

TIGR_TO_CATEGORY = {
    "Amino acid biosynthesis": "Amino acid metabolism",
    "Biosynthesis of cofactors, prosthetic groups, and carriers": "Coenzyme metabolism",
    "Cell envelope": "Cell wall and membrane",
    "Cellular processes": "Cellular processes",
    "Central intermediary metabolism": "Central intermediary metabolism",
    "DNA metabolism": "Replication and repair",
    "Disrupted reading frame": "Unknown",
    "Energy metabolism": "Energy production",
    "Fatty acid and phospholipid metabolism": "Lipid metabolism",
    "Hypothetical proteins": "Unknown",
    "Mobile and extrachromosomal element functions": "Mobile elements",
    "Not Found": "Unknown",
    "Protein fate": "Post-translational modification",
    "Protein synthesis": "Translation",
    "Purines, pyrimidines, nucleosides, and nucleotides": "Nucleotide metabolism",
    "Regulatory functions": "Regulatory functions",
    "Signal transduction": "Signal transduction",
    "Transcription": "Transcription",
    "Transport and binding proteins": "Transport",
    "Unclassified": "Unknown",
    "Unknown function": "Unknown",
}

# All valid output values — used for build-time assertion
VALID_CATEGORIES = frozenset(
    set(COG_TO_CATEGORY.values())
    | set(CYANORAK_TO_CATEGORY.values())
    | set(CYANORAK_D_SUBCODES.values())
    | set(TIGR_TO_CATEGORY.values())
    | {"Unknown"}
)


def _compute_gene_category(result: dict) -> str:
    """Compute gene_category from Cyanorak Role → TIGR Role → COG category."""
    gene_category = None

    # Priority 1: Cyanorak Role (cyanobacteria only)
    cyanorak_roles = result.get("cyanorak_Role", [])
    if cyanorak_roles:
        code = cyanorak_roles[0]
        top_letter = code.split(".")[0]
        if top_letter == "D":
            parts = code.split(".")
            sub_key = ".".join(parts[:2]) if len(parts) >= 2 else "D"
            gene_category = CYANORAK_D_SUBCODES.get(sub_key, "Cellular processes")
        else:
            gene_category = CYANORAK_TO_CATEGORY.get(top_letter)

    # Priority 2: TIGR Role (skip if we already have a real category)
    if not gene_category or gene_category == "Unknown":
        tigr_descs = result.get("tIGR_Role_description", [])
        if tigr_descs:
            main_role = tigr_descs[0].split(" / ")[0].strip()
            cat = TIGR_TO_CATEGORY.get(main_role)
            if cat and cat != "Unknown":
                gene_category = cat

    # Priority 3: COG category (universal, from eggNOG)
    if not gene_category or gene_category == "Unknown":
        cog_cats = result.get("cog_category", [])
        if cog_cats:
            cat = COG_TO_CATEGORY.get(cog_cats[0])
            if cat:
                gene_category = cat

    return gene_category or "Unknown"


# ─── paths ────────────────────────────────────────────────────────────────────

DEFAULT_CONFIG = PROJECT_ROOT / "config/gene_annotations_config.yaml"


# ─── data loaders ─────────────────────────────────────────────────────────────

def load_gene_mapping(data_dir: str) -> dict[str, dict]:
    """Load gene_mapping.csv → {locus_tag: {col: value}}."""
    path = os.path.join(data_dir, "gene_mapping.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"gene_mapping.csv not found: {path}")
    result: dict[str, dict] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            lt = (row.get("locus_tag") or "").strip()
            if lt:
                result[lt] = dict(row)
    return result


def load_eggnog(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load .emapper.annotations → {query_wp_id: {col: value}}.

    Skips '##' comment lines; strips leading '#' from column names.
    """
    path = os.path.join(data_dir, "eggnog", f"{strain_name}.emapper.annotations")
    if not os.path.exists(path):
        return {}
    result: dict[str, dict] = {}
    with open(path, newline="") as f:
        lines = (line for line in f if not line.startswith("##"))
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            # Strip '#' prefix from '#query' column name
            clean_row = {k.lstrip("#"): v for k, v in row.items()}
            query = clean_row.get("query", "").strip()
            if query and query != "-":
                result[query] = clean_row
    return result


def load_uniprot(
    data_dir: str,
    ncbi_taxon_id: int | None,
    organism_group: str,
) -> dict[str, dict]:
    """Load protein_annotations.json → {refseq_wp_id: {field: value}}.

    The JSON is row-oriented: data[uniprot_id] = {field: value}.
    Re-indexes by refseq_ids (WP_ accessions) for joining with gene_mapping.
    """
    if not ncbi_taxon_id:
        print(f"  [uniprot] No taxon ID provided — skipping UniProt data")
        return {}

    path = str(
        PROJECT_ROOT / "cache" / "data" / organism_group
        / "uniprot" / str(ncbi_taxon_id) / "protein_annotations.json"
    )
    if not os.path.exists(path):
        print(f"  [uniprot] No data found at {path}")
        return {}

    with open(path) as f:
        rows: dict[str, dict] = json.load(f)

    # Re-index by RefSeq accession (refseq_ids → protein_id in gene_mapping)
    result: dict[str, dict] = {}
    for uid, row in rows.items():
        refseq_ids = row.get("refseq_ids", [])
        if not refseq_ids:
            continue
        entry = dict(row)
        entry["uniprot_accession"] = uid
        for rs_id in refseq_ids:
            rs_id = rs_id.strip()
            if rs_id and rs_id != "-":
                result[rs_id] = entry

    return result


# ─── annotation builder ───────────────────────────────────────────────────────

class AnnotationBuilder:
    """Applies merge rules from config YAML to produce wide + merged dicts."""

    def __init__(self, config: dict):
        self.field_configs: dict[str, dict] = config.get("fields", {})

    # ── raw value fetcher ──────────────────────────────────────────────────────

    def _get_raw(
        self,
        src_cfg: dict,
        gm: dict,
        eg: dict,
        up: dict,
    ) -> Any:
        """Fetch raw value from source row according to src_cfg spec."""
        source = src_cfg.get("source", "")
        field = src_cfg.get("field", "")

        if source == "gene_mapping":
            raw = gm.get(field)
        elif source == "eggnog":
            raw = eg.get(field)
        elif source == "uniprot":
            raw = up.get(field)
        else:
            return None

        # URL-decode string values (gene_mapping uses URL encoding)
        if isinstance(raw, str):
            raw = unquote(raw.strip())

        if not _nonempty(raw):
            return None

        return raw

    # ── apply transform to a single value ─────────────────────────────────────

    def _apply_transform(self, transform: str | None, value: Any) -> Any:
        """Apply a named transform to a value; returns empty string on failure."""
        if not transform or transform not in _TRANSFORMS:
            return value
        fn = _TRANSFORMS[transform]
        if isinstance(value, list):
            return [fn(v) for v in value if _nonempty(v)]
        return fn(value)

    # ── resolver: passthrough ──────────────────────────────────────────────────

    def _resolve_passthrough(
        self, fconf: dict, gm: dict, eg: dict, up: dict
    ) -> Any:
        raw = self._get_raw(fconf, gm, eg, up)
        if not _nonempty(raw):
            return None
        transform = fconf.get("transform")
        if transform:
            raw = self._apply_transform(transform, raw)
        return raw if _nonempty(raw) else None

    # ── resolver: passthrough_list ─────────────────────────────────────────────

    def _resolve_passthrough_list(
        self, fconf: dict, gm: dict, eg: dict, up: dict
    ) -> list[str] | None:
        raw = self._get_raw(fconf, gm, eg, up)
        if not _nonempty(raw):
            return None
        delimiter = fconf.get("delimiter", ",")
        tokens = _coerce_to_tokens(raw, delimiter)
        return tokens if tokens else None

    # ── resolver: single ──────────────────────────────────────────────────────

    def _resolve_single(
        self,
        fconf: dict,
        gm: dict,
        eg: dict,
        up: dict,
        source_tracking: dict,
        locus_tag: str = "",
    ) -> Any:
        """First non-empty candidate wins; record source if track_source set.

        Candidates may have an optional ``source_label`` key to override the
        recorded provenance string (e.g. 'cyanorak' instead of 'gene_mapping').

        When ``reject_identifiers`` is set on *fconf*, identifier-style values
        (matching ``_IDENTIFIER_RE`` or equal to the gene's locus_tag) are
        skipped so the next candidate can provide a real biological name.
        """
        track_key = fconf.get("track_source")
        reject_ids = fconf.get("reject_identifiers", False)
        for cand in fconf.get("candidates", []):
            raw = self._get_raw(cand, gm, eg, up)
            if not _nonempty(raw):
                continue
            transform = cand.get("transform")
            if transform == "first_token_space":
                # Handle list: take first element, then first token
                if isinstance(raw, list):
                    raw = raw[0] if raw else ""
                val = _tx_first_token_space(str(raw))
            elif transform == "strip_function_prefix":
                val = _tx_strip_function_prefix(
                    raw[0] if isinstance(raw, list) else str(raw)
                )
            elif transform:
                val = self._apply_transform(transform, raw)
            else:
                # Lists: join with space for display, or take first
                if isinstance(raw, list):
                    val = raw[0] if raw else ""
                else:
                    val = raw
            if _nonempty(val):
                # Skip identifier-style values when reject_identifiers is set
                if reject_ids and isinstance(val, str) and (
                    val == locus_tag or _IDENTIFIER_RE.match(val)
                ):
                    continue
                if track_key:
                    # Use source_label if set, otherwise fall back to source name
                    source_tracking[track_key] = cand.get("source_label", cand["source"])
                return val
        return None

    # ── resolver: union ────────────────────────────────────────────────────────

    def _resolve_union(
        self, fconf: dict, gm: dict, eg: dict, up: dict
    ) -> list[str] | None:
        """Merge tokens from all sources, deduplicate, apply global filter."""
        global_filter = fconf.get("filter")
        global_filter_not = fconf.get("filter_not")
        seen: dict[str, None] = {}  # ordered set

        for src_cfg in fconf.get("sources", []):
            raw = self._get_raw(src_cfg, gm, eg, up)
            if not _nonempty(raw):
                continue

            delimiter = src_cfg.get("delimiter", ",")
            transform = src_cfg.get("transform")

            # Special-case transforms that produce lists from a string
            if transform == "extract_pfam_ids":
                if isinstance(raw, list):
                    tokens = [t for t in raw if str(t).startswith("PF")]
                else:
                    tokens = _tx_extract_pfam_ids(str(raw))
            elif transform == "extract_pfam_names":
                if isinstance(raw, list):
                    tokens = [t for t in raw if t and not str(t).startswith("PF")]
                else:
                    tokens = _tx_extract_pfam_names(str(raw))
            elif transform == "extract_go_from_pipe":
                base_tokens = _coerce_to_tokens(raw, delimiter)
                tokens = [_tx_extract_go_from_pipe(t) for t in base_tokens]
            elif transform:
                base_tokens = _coerce_to_tokens(raw, delimiter)
                fn = _TRANSFORMS.get(transform)
                if fn:
                    raw_tokens = [fn(t) for t in base_tokens]
                    # Flatten: transforms may return lists (e.g. normalize_ec
                    # with multiple successors)
                    tokens = []
                    for t in raw_tokens:
                        if isinstance(t, list):
                            tokens.extend(t)
                        else:
                            tokens.append(t)
                else:
                    tokens = base_tokens
            else:
                tokens = _coerce_to_tokens(raw, delimiter)

            for tok in tokens:
                tok = str(tok).strip()
                if not tok or tok == "-":
                    continue
                if global_filter and not re.match(global_filter, tok):
                    continue
                if global_filter_not and re.match(global_filter_not, tok):
                    continue
                seen[tok] = None

        result = list(seen.keys())
        return result if result else None

    # ── resolver: integer / float ──────────────────────────────────────────────

    def _resolve_integer(
        self, fconf: dict, gm: dict, eg: dict, up: dict
    ) -> int | None:
        raw = self._get_raw(fconf, gm, eg, up)
        if raw is None:
            return None
        try:
            return int(float(str(raw).strip()))
        except (ValueError, TypeError):
            return None

    def _resolve_float(
        self, fconf: dict, gm: dict, eg: dict, up: dict
    ) -> float | None:
        raw = self._get_raw(fconf, gm, eg, up)
        if raw is None:
            return None
        try:
            return float(str(raw).strip())
        except (ValueError, TypeError):
            return None

    # ── build wide ────────────────────────────────────────────────────────────

    def build_wide(
        self,
        gm: dict,
        eg: dict,
        up: dict,
    ) -> dict:
        """All source fields, source-prefixed — full audit trail."""
        wide: dict[str, Any] = {}
        for k, v in gm.items():
            if _nonempty(v):
                wide[f"gene_mapping_{k}"] = v
        for k, v in eg.items():
            if _nonempty(v):
                wide[f"eggnog_{k}"] = v
        for k, v in up.items():
            if _nonempty(v):
                wide[f"uniprot_{k}"] = v
        return wide

    # ── build merged ──────────────────────────────────────────────────────────

    def build_merged(
        self,
        gm: dict,
        eg: dict,
        up: dict,
        organism_name: str | None = None,
    ) -> dict:
        """Apply merge rules → canonical field set."""
        result: dict[str, Any] = {}
        source_tracking: dict[str, str] = {}
        locus_tag = gm.get("locus_tag", "")

        for canonical_field, fconf in self.field_configs.items():
            ftype = fconf.get("type", "passthrough")

            if ftype == "single":
                val = self._resolve_single(fconf, gm, eg, up, source_tracking, locus_tag)
            elif ftype == "union":
                val = self._resolve_union(fconf, gm, eg, up)
            elif ftype == "passthrough":
                val = self._resolve_passthrough(fconf, gm, eg, up)
            elif ftype == "passthrough_list":
                val = self._resolve_passthrough_list(fconf, gm, eg, up)
            elif ftype == "integer":
                val = self._resolve_integer(fconf, gm, eg, up)
            elif ftype == "float":
                val = self._resolve_float(fconf, gm, eg, up)
            elif ftype == "extract_first_match":
                val = extract_first_match_in_sources(
                    fconf.get("sources", []), gm, eg, up,
                    fconf.get("pattern", ""),
                    fconf.get("extract_group", 0),
                )
            else:
                continue

            if _nonempty(val):
                result[canonical_field] = val

        # Remove canonical gene_name from synonym lists to avoid duplication
        gene_name = result.get("gene_name", "")
        if gene_name:
            for field in ("gene_synonyms", "gene_name_synonyms"):
                if field in result:
                    filtered = [s for s in result[field] if s != gene_name]
                    if filtered:
                        result[field] = filtered
                    else:
                        del result[field]

        # Add source-tracking fields collected during 'single' resolution
        result.update(source_tracking)

        # Compute annotation_quality (0–3) based on product content + structured annotations
        product = result.get("product", "")
        func_desc = result.get("function_description", "")
        is_hypothetical = not product or bool(re.match(
            r'^(hypothetical|conserved hypothetical|uncharacterized)\b', product, re.IGNORECASE
        ))

        if not is_hypothetical:
            # Count structured annotations: go_terms, kegg_ko, ec_numbers, pfam_ids
            structured_count = sum(
                1 for f in ("go_terms", "kegg_ko", "ec_numbers", "pfam_ids")
                if result.get(f)
            )
            quality = 3 if structured_count >= 2 else 2
        elif func_desc and func_desc != "-":
            quality = 1
        else:
            quality = 0
        result["annotation_quality"] = quality

        # Compute gene_category — high-level functional classification
        category = _compute_gene_category(result)
        assert category in VALID_CATEGORIES, (
            f"Invalid gene_category {category!r} for {result.get('locus_tag')}"
        )
        result["gene_category"] = category

        # Collect all source descriptions for LLM summaries
        alt_descriptions: list[str] = []
        alt_descriptions_set: set[str] = set()

        def _add_desc(label: str, text: str) -> None:
            s = f"[{label}] {text.strip()}" if text and text.strip() else None
            if s and s not in alt_descriptions_set:
                alt_descriptions.append(s)
                alt_descriptions_set.add(s)

        cyanorak_prod = unquote((gm.get("product_cyanorak") or "").strip())
        ncbi_prod = unquote((gm.get("product") or "").strip())
        eg_desc = (eg.get("Description") or "").strip()
        up_prod = (up.get("product") or "").strip()
        up_func = (up.get("function_description") or "").strip()
        if up_func.startswith("FUNCTION: "):
            up_func = up_func[len("FUNCTION: "):]
        up_family = (up.get("protein_family") or "").strip()

        _add_desc("cyanorak", cyanorak_prod)
        _add_desc("ncbi", ncbi_prod)
        _add_desc("eggnog", eg_desc)
        _add_desc("uniprot_product", up_prod)
        _add_desc("uniprot", up_func)
        _add_desc("protein_family", up_family)
        for desc in result.get("cyanorak_Role_description", []):
            _add_desc("cyanorak_role", desc)
        for desc in result.get("tIGR_Role_description", []):
            _add_desc("tigr_role", desc)
        for desc in result.get("eggnog_og_descriptions", []):
            _add_desc("cog", desc)
        for desc in result.get("kegg_ko_descriptions", []):
            _add_desc("kegg", desc)
        pfam_names_list = result.get("pfam_names", [])
        pfam_descs_list = result.get("pfam_descriptions", [])
        for i, name in enumerate(pfam_names_list):
            pfam_text = f"{name}: {pfam_descs_list[i]}" if i < len(pfam_descs_list) and pfam_descs_list[i] else name
            _add_desc("pfam", pfam_text)

        if alt_descriptions:
            result["alternate_functional_descriptions"] = alt_descriptions

        # ── Computed fields for MCP gene lookup ──────────────────────────────

        # organism_strain — preferred organism name (e.g., "Prochlorococcus MED4")
        if organism_name:
            result["organism_strain"] = organism_name

        # gene_summary — primary display field: "gene_name :: product :: description"
        gene_name = result.get("gene_name", "")
        locus_tag = result.get("locus_tag", "")
        # Clear gene_name when it's just an identifier, not a biological name
        if gene_name and (gene_name == locus_tag or _IDENTIFIER_RE.match(gene_name)):
            result["gene_name"] = None
            gene_name = ""
        product = result.get("product", "")
        best_desc = up_func or eg_desc or cyanorak_prod or ncbi_prod
        if best_desc == product:
            best_desc = ""
        # Skip uninformative "domain/protein of unknown function" descriptions
        if best_desc and re.match(r'^(Protein |Domain )of unknown function', best_desc):
            best_desc = ""
        summary_parts = [p for p in [gene_name, product, best_desc] if p]
        if summary_parts:
            result["gene_summary"] = " :: ".join(summary_parts)

        # all_identifiers — union of all alternative ID fields for get_gene lookup
        # Excludes locus_tag and gene_name (they have their own scalar indexes)
        scalar_indexed = {result.get("locus_tag"), gene_name} - {None, ""}
        all_ids = set(filter(None, [
            result.get("locus_tag_ncbi"),
            result.get("locus_tag_cyanorak"),
            result.get("protein_id"),
        ]))
        all_ids.update(result.get("old_locus_tags") or [])
        all_ids.update(result.get("alternative_locus_tags") or [])
        all_ids.update(result.get("gene_name_synonyms") or [])
        all_ids -= scalar_indexed
        if all_ids:
            result["all_identifiers"] = sorted(all_ids)

        return result


# ─── per-strain pipeline ──────────────────────────────────────────────────────

def process_strain(
    row: dict,
    config: dict,
    force: bool = False,
) -> None:
    strain_name = row["strain_name"]
    preferred_name = (row.get("preferred_name") or "").strip() or strain_name
    data_dir = row["data_dir"].rstrip("/")
    taxon_id_str = (row.get("ncbi_taxon_id") or "").strip()
    ncbi_taxon_id = int(taxon_id_str) if taxon_id_str else None
    organism_group = infer_organism_group(data_dir)

    wide_path = os.path.join(data_dir, "gene_annotations_wide.json")
    merged_path = os.path.join(data_dir, "gene_annotations_merged.json")

    if not force and os.path.exists(merged_path):
        print(f"[{strain_name}] Skipping (already exists). Use --force to rebuild.")
        return

    print(f"\n[{strain_name}] Loading sources...")

    gm_data = load_gene_mapping(data_dir)
    eg_data = load_eggnog(data_dir, strain_name)
    up_data: dict[str, dict] = {}
    if ncbi_taxon_id:
        up_data = load_uniprot(data_dir, ncbi_taxon_id, organism_group)

    print(f"  gene_mapping : {len(gm_data):>5} genes")
    print(f"  eggnog       : {len(eg_data):>5} entries")
    print(f"  uniprot      : {len(up_data):>5} entries (keyed by RefSeq)")

    builder = AnnotationBuilder(config)

    wide_out: dict[str, dict] = {}
    merged_out: dict[str, dict] = {}

    stats = dict(total=0, eggnog_hit=0, uniprot_hit=0,
                 has_product=0, has_go=0, has_cog=0, has_kegg_ko=0,
                 quality_0=0, quality_1=0, quality_2=0, quality_3=0)

    for locus_tag, gm_row in gm_data.items():
        protein_id = (gm_row.get("protein_id") or "").strip()
        eg_row = eg_data.get(protein_id, {})
        up_row = up_data.get(protein_id, {})

        stats["total"] += 1
        if eg_row:
            stats["eggnog_hit"] += 1
        if up_row:
            stats["uniprot_hit"] += 1

        wide_out[locus_tag] = builder.build_wide(gm_row, eg_row, up_row)
        merged = builder.build_merged(gm_row, eg_row, up_row, organism_name=preferred_name)
        merged_out[locus_tag] = merged

        if merged.get("product"):
            stats["has_product"] += 1
        if merged.get("go_terms"):
            stats["has_go"] += 1
        if merged.get("cog_category"):
            stats["has_cog"] += 1
        if merged.get("kegg_ko"):
            stats["has_kegg_ko"] += 1
        q = merged.get("annotation_quality", 0)
        stats[f"quality_{q}"] += 1
        cat = merged.get("gene_category", "Unknown")
        stats.setdefault(f"cat_{cat}", 0)
        stats[f"cat_{cat}"] += 1

    with open(wide_path, "w") as f:
        json.dump(wide_out, f, indent=2, sort_keys=True)
    print(f"  → {wide_path}")

    with open(merged_path, "w") as f:
        json.dump(merged_out, f, indent=2, sort_keys=True)
    print(f"  → {merged_path}")

    # Coverage report
    n = stats["total"] or 1
    pct = lambda k: f"{100 * stats[k] // n}%"
    print(f"\n  === {strain_name} coverage ===")
    print(f"  Genes:          {stats['total']}")
    print(f"  EggNOG matched: {stats['eggnog_hit']} ({pct('eggnog_hit')})")
    print(f"  UniProt matched:{stats['uniprot_hit']} ({pct('uniprot_hit')})")
    print(f"  Has product:    {stats['has_product']} ({pct('has_product')})")
    print(f"  Has GO terms:   {stats['has_go']} ({pct('has_go')})")
    print(f"  Has COG:        {stats['has_cog']} ({pct('has_cog')})")
    print(f"  Has KEGG KO:    {stats['has_kegg_ko']} ({pct('has_kegg_ko')})")
    print(f"  Quality 0/1/2/3: "
          f"{stats['quality_0']}/{stats['quality_1']}/"
          f"{stats['quality_2']}/{stats['quality_3']}")
    # gene_category distribution (sorted by count descending)
    cat_stats = {k[4:]: v for k, v in stats.items() if k.startswith("cat_")}
    if cat_stats:
        print(f"\n  === {strain_name} gene_category ===")
        for cat, count in sorted(cat_stats.items(), key=lambda x: -x[1]):
            print(f"    {cat:<40s} {count:>5} ({100 * count // n}%)")


# ─── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build gene annotation tables by merging NCBI/Cyanorak, EggNOG, and UniProt."
    )
    add_common_args(parser, DEFAULT_CONFIG)
    parser.add_argument(
        "--llm-summary", action="store_true",
        help="Generate LLM summaries per gene (Step 1C — not yet implemented)",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    rows = load_genome_rows(args.strains)

    print(f"Processing {len(rows)} strain(s) with config: {args.config}")
    for row in rows:
        process_strain(row, config, force=args.force)

    if args.llm_summary:
        print("\nLLM summary generation (Step 1C) not yet implemented.")


if __name__ == "__main__":
    main()
