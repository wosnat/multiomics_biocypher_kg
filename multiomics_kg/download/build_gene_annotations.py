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
import logging
import os
import re
from pathlib import Path
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
    _tx_first_token_space,
    _tx_strip_function_prefix,
    _tx_strip_prefix_ko,
    is_eggnog_description_stub,
)
from multiomics_kg.download.utils.cli import add_common_args, load_config, load_genome_rows
from multiomics_kg.download.utils.ortholog_group_utils import (
    extract_ortholog_groups,
    organism_group_from_path,
)
from multiomics_kg.utils.pfam_utils import PfamData, load_pfam_data
from multiomics_kg.download import build_interpro_reference
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


logger = logging.getLogger(__name__)


# ─── Pfam post-merge enrichment ──────────────────────────────────────────────

def enrich_pfam_fields(gene: dict, pfam_data: PfamData) -> list[str]:
    """Resolve raw pfam_ids tokens to clean PF* accessions and update related fields.

    pfam_ids at this point is a raw union of all tokens from Cyanorak protein_domains,
    eggNOG PFAMs, and UniProt pfam_ids -- a mix of PF* IDs, shortnames, TIGR*, IPR*.

    For each token:
    - PF* -> keep as-is (already a valid Pfam accession)
    - Found in pfam_data.by_shortname -> resolve to PF* ID
    - Otherwise -> drop (TIGR*, IPR*, or unresolved shortname)

    Side effects on gene dict:
    - pfam_ids: overwritten with clean PF* list (or deleted if empty)
    - alternate_functional_descriptions: appends "[pfam] shortname: description"
      entries from Pfam reference data for each resolved PF* ID

    Returns list of unresolved non-PF* tokens (for logging).
    """
    raw_tokens = gene.get("pfam_ids")
    if not raw_tokens:
        return []

    clean_ids: dict[str, None] = {}  # ordered set for dedup
    unresolved: list[str] = []
    # Per-token provenance map keyed by the RAW token (from _resolve_union); it
    # must be re-keyed to the resolved PF* accession so it stays aligned with the
    # rewritten pfam_ids list (and usable for edge `sources` downstream).
    raw_source_map: dict[str, list[str]] = gene.get("pfam_ids_source") or {}
    new_source_map: dict[str, set[str]] = {}

    for token in raw_tokens:
        token = str(token).strip()
        if not token:
            continue
        if token.startswith("PF"):
            resolved = token
        elif token in pfam_data.by_shortname:
            resolved = pfam_data.by_shortname[token]
        else:
            # TIGR*, IPR*, or unresolved shortname -- drop
            unresolved.append(token)
            continue
        clean_ids[resolved] = None
        if raw_source_map.get(token):
            new_source_map.setdefault(resolved, set()).update(raw_source_map[token])

    result = list(clean_ids.keys())
    if result:
        gene["pfam_ids"] = result
        if new_source_map:
            gene["pfam_ids_source"] = {k: sorted(v) for k, v in new_source_map.items()}
    elif "pfam_ids" in gene:
        del gene["pfam_ids"]
        gene.pop("pfam_ids_source", None)
        gene.pop("pfam_ids_evidence", None)

    # Recompute contributing_sources after enrichment (in case sources changed)
    gene["contributing_sources"] = _compute_contributing_sources(gene)

    # Add Pfam descriptions to alternate_functional_descriptions
    if result:
        alt_descs = gene.get("alternate_functional_descriptions", [])
        alt_descs_set = set(alt_descs)
        for pf_id in result:
            entry = pfam_data.by_accession.get(pf_id)
            if entry:
                text = f"{entry.shortname}: {entry.description}" if entry.description else entry.shortname
                desc = f"[pfam] {text}"
                if desc not in alt_descs_set:
                    alt_descs.append(desc)
                    alt_descs_set.add(desc)
        gene["alternate_functional_descriptions"] = alt_descs

    return unresolved


# ─── InterPro entry-xref propagation (Layer B) ────────────────────────────────

# Curated / direct sources; when any of these also asserts a token, the token is
# curated regardless of the InterPro inference that corroborates it.
_CURATED_SOURCES = {"ncbi", "cyanorak", "uniprot", "eggnog"}
# GO / CAZy propagate from FAMILY + DOMAIN entries; fold-level (HOMOLOGOUS_SUPERFAMILY)
# is shape-only and excluded. EC propagates from FAMILY only, single-EC (see below).
_INTERPRO_PROPAGATE_TYPES = {"FAMILY", "DOMAIN"}
_BARE_EC_RE = re.compile(r"^\d+\.\d+\.\d+$")            # 3.4.21  (no 4th field)
_VALID_EC_RE = re.compile(r"^\d+\.[\d\-]+\.[\d\-]+[\.\-]")  # mirrors the ec_numbers filter


def _normalize_interpro_ec(raw: str) -> list[str]:
    """Normalise one raw InterPro EC token → validated EC list (may be empty).

    Bare 3-level ECs (``3.4.21``) get ``.-`` appended so they satisfy the 4-field
    ``ec_numbers`` filter, then the shared ``normalize_ec`` transform remaps
    obsolete/transferred numbers. Anything still failing the EC shape is dropped.
    """
    raw = (raw or "").strip()
    if _BARE_EC_RE.match(raw):
        raw = raw + ".-"
    fn = _TRANSFORMS.get("normalize_ec")
    out = fn(raw) if fn else raw
    toks = out if isinstance(out, list) else [out]
    return [t for t in (str(x).strip() for x in toks) if t and _VALID_EC_RE.match(t)]


def _fold_interpro_field(gene: dict, field: str, new_tokens: dict[str, str]) -> None:
    """Merge InterPro-contributed *new_tokens* ({token: strength}) into *field*.

    strength ∈ {family, domain, signature}. For each token: append to the field
    list if new; add ``interpro`` to the per-token ``<field>_source`` map; and set
    ``<field>_evidence[token]`` to the strongest applicable evidence — ``curated``
    when any curated source also asserts it, else ``signature`` (direct Pfam HMM),
    ``family_inferred`` or ``domain_inferred``. The evidence map is sparse: only
    InterPro-touched tokens get an entry; consumers default the rest to ``curated``.
    """
    if not new_tokens:
        return
    lst = gene.get(field) or []
    seen = set(lst)
    src_map = dict(gene.get(f"{field}_source") or {})
    ev_map = dict(gene.get(f"{field}_evidence") or {})
    for tok, strength in new_tokens.items():
        if tok not in seen:
            lst.append(tok)
            seen.add(tok)
        srcs = set(src_map.get(tok, []))
        srcs.add("interpro")
        src_map[tok] = sorted(srcs)
        if srcs & _CURATED_SOURCES:
            ev = "curated"
        elif strength == "signature":
            ev = "signature"
        elif strength == "family":
            ev = "family_inferred"
        else:
            ev = "domain_inferred"
        ev_map[tok] = ev
    gene[field] = lst
    gene[f"{field}_source"] = src_map
    gene[f"{field}_evidence"] = ev_map


def enrich_interpro_fields(gene: dict, ipr_row: dict, interpro_ref: dict,
                            ncbifam_ref: dict | None = None) -> None:
    """Promote InterPro/NCBIfam entry-level xrefs into gene fields (Layer B).

    Noise-gated, type-aware propagation (design 2026-08-10 §2/§5.1, extended
    2026-08-17 for donor attribution + naming recovery):
    - ``go_terms``: donor-attributed — a GO is added iff at least one of its
      donor entries (``ipr_row["go_term_donors"]``, ``{GO: [IPR, ...]}``) is
      FAMILY or DOMAIN typed; evidence is ``family`` if any donor is FAMILY,
      else ``domain``. When *ipr_row* carries no ``go_term_donors`` key (old
      callers / rows predating Task 8), donors are derived from this gene's
      ``interpro_entries`` × *interpro_ref* exactly as the pre-donor code did.
    - ``cazy_ids``: FAMILY + DOMAIN entries (fold excluded) — unchanged.
    - ``ec_numbers``: FAMILY entries carrying **exactly one** EC (a multi-EC family
      is a candidate set, not a claim — those live in Layer A, Phase 3).
    - ``pfam_ids``: direct PFAM signature hits from calls.json (no inference).
    - ``alternate_functional_descriptions``: ``[interpro] <name>`` (FAMILY/DOMAIN,
      as before) plus naming recovery — ``[hamap] <desc>`` (from
      ``ipr_row["hamap_descriptions"]``) and ``[ncbifam] <product name>`` (from
      *ncbifam_ref* via the gene's ``ncbifam_ids``). hamap/ncbifam tokens are
      skipped when their text case-insensitively equals the gene's ``product``
      (in addition to the existing exact-string afd dedup, which still applies
      to all tokens).
    - ``gene_name``: fill-if-empty from the first (sorted) ``gene_symbol`` among
      the gene's ``ncbifam_ids`` reference entries; sets
      ``gene_name_source = "ncbifam"``. Never overwrites an existing gene_name.

    Each ontology token is tagged with ``interpro`` in ``<field>_source`` and an
    ``<field>_evidence`` strength via ``_fold_interpro_field``. Idempotent-ish:
    safe to run once per gene after ``build_merged``.
    """
    entries = gene.get("interpro_entries") or []

    # ── GO: donor-attributed gate ───────────────────────────────────────────
    donors = ipr_row.get("go_term_donors")
    if donors is None:
        # Fallback for rows without donor attribution: derive donors from this
        # gene's entries × the reference (matches the pre-donor behavior).
        donors = {}
        for ipr_id in entries:
            meta = interpro_ref.get(ipr_id)
            if not meta:
                continue
            for go in meta.get("go_terms", []):
                donors.setdefault(go, []).append(ipr_id)

    go_new: dict[str, str] = {}
    for go, donor_ids in donors.items():
        strength = None
        for ipr_id in donor_ids:
            meta = interpro_ref.get(ipr_id) or {}
            etype = (meta.get("type") or "").upper()
            if etype == "FAMILY":
                strength = "family"
                break
            if etype == "DOMAIN" and strength is None:
                strength = "domain"
        if strength is not None:
            go_new[go] = strength

    # ── EC / CAZy / interpro-name descriptions (unchanged logic) ───────────
    ec_new: dict[str, str] = {}
    cazy_new: dict[str, str] = {}
    desc_entries: list[str] = []
    for ipr_id in entries:
        meta = interpro_ref.get(ipr_id)
        if not meta:
            continue
        etype = (meta.get("type") or "").upper()
        strength = "family" if etype == "FAMILY" else "domain"
        if etype in _INTERPRO_PROPAGATE_TYPES:
            for cz in meta.get("cazy_ids", []):
                cazy_new.setdefault(cz, strength)
            name = meta.get("name")
            if name:
                desc_entries.append(f"[interpro] {name}")
        if etype == "FAMILY":
            ecs = meta.get("ec_numbers", [])
            if len(ecs) == 1:  # single-EC gate — noise > data otherwise
                for norm in _normalize_interpro_ec(ecs[0]):
                    ec_new.setdefault(norm, "family")

    # Direct PFAM signature hits (no inference), surfaced by load_interproscan().
    pfam_new: dict[str, str] = {}
    for sig in (ipr_row.get("pfam_signatures") or []):
        if sig.startswith("PF"):
            pfam_new.setdefault(sig, "signature")

    # ── Naming recovery: hamap + ncbifam, deduped against `product` ────────
    product_lc = (gene.get("product") or "").strip().lower()

    def _maybe_desc(tag: str, text: str) -> None:
        text = (text or "").strip()
        if not text or text.lower() == product_lc:
            return
        desc_entries.append(f"[{tag}] {text}")

    for h in (ipr_row.get("hamap_descriptions") or []):
        _maybe_desc("hamap", h)
    symbols: list[str] = []
    for acc in (gene.get("ncbifam_ids") or []):
        meta = (ncbifam_ref or {}).get(acc) or {}
        _maybe_desc("ncbifam", meta.get("name") or "")
        if meta.get("gene_symbol"):
            symbols.append(meta["gene_symbol"])
    if symbols and not gene.get("gene_name"):
        gene["gene_name"] = sorted(symbols)[0]
        gene["gene_name_source"] = "ncbifam"

    _fold_interpro_field(gene, "go_terms", go_new)
    _fold_interpro_field(gene, "ec_numbers", ec_new)
    _fold_interpro_field(gene, "cazy_ids", cazy_new)
    _fold_interpro_field(gene, "pfam_ids", pfam_new)

    if desc_entries:
        afd = gene.get("alternate_functional_descriptions") or []
        afd_set = set(afd)
        for d in desc_entries:
            if d not in afd_set:
                afd.append(d)
                afd_set.add(d)
        gene["alternate_functional_descriptions"] = afd


# ─── contributing_sources ────────────────────────────────────────────────────

def _has_source_label(gene: dict, label: str) -> bool:
    """Check whether a gene has any field provenance-tagged with `label`.

    Walks `gene["*_source"]` track fields plus `[label]` prefixes inside
    `alternate_functional_descriptions`. Two `*_source` shapes coexist:
    - **scalar** (single-resolver fields): `product_source == "cyanorak"`.
    - **per-token map** (union fields, e.g. `go_terms_source`):
      `{"GO:0003677": ["uniprot", "interpro"]}` — the label matches if it appears
      in any token's source list.
    """
    for k, v in gene.items():
        if not k.endswith("_source"):
            continue
        if isinstance(v, str):
            if v == label:
                return True
        elif isinstance(v, dict):
            if any(label in srcs for srcs in v.values()):
                return True
    # [label] prefix in alternate functional descriptions
    afd = gene.get("alternate_functional_descriptions") or []
    for entry in afd:
        if isinstance(entry, str) and entry.startswith(f"[{label}]"):
            return True
    return False


def _compute_contributing_sources(gene: dict) -> list[str]:
    """Return sorted list of data sources that contributed at least one field.

    Source presence rules:
    - 'ncbi': always present (every Gene comes from an NCBI GFF row).
    - 'cyanorak': locus_tag_cyanorak non-null OR any cyanorak-tagged field.
    - 'uniprot': uniprot_accession non-null OR any uniprot-tagged field.
    - 'eggnog': seed_ortholog or eggnog_ogs non-null OR any eggnog-tagged field.
    """
    sources = {"ncbi"}
    if gene.get("locus_tag_cyanorak") or _has_source_label(gene, "cyanorak"):
        sources.add("cyanorak")
    if gene.get("uniprot_accession") or _has_source_label(gene, "uniprot"):
        sources.add("uniprot")
    if (gene.get("seed_ortholog")
            or gene.get("eggnog_ogs")
            or _has_source_label(gene, "eggnog")):
        sources.add("eggnog")
    if (gene.get("psortb_localization")
            or gene.get("psortb_score") is not None
            or _has_source_label(gene, "psortb")):
        sources.add("psortb")
    if (gene.get("signalp_type")
            or gene.get("signalp_probability") is not None
            or _has_source_label(gene, "signalp")):
        sources.add("signalp")
    if (gene.get("interpro_entries")
            or _has_source_label(gene, "interproscan")):
        sources.add("interproscan")
    if (gene.get("tcdb_diamond_ids")
            or _has_source_label(gene, "tcdb_diamond")):
        sources.add("tcdb_diamond")
    if (gene.get("merops_ids")
            or _has_source_label(gene, "merops_diamond")):
        sources.add("merops_diamond")
    return sorted(sources)


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

    Skips '##' comment lines; strips leading '#' from column names. Replaces
    eggNOG-internal stub strings (e.g. "Alternative locus ID") in the
    Description column with the empty string so downstream consumers treat them
    as absent rather than as content.
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
                if is_eggnog_description_stub(clean_row.get("Description")):
                    clean_row["Description"] = ""
                result[query] = clean_row
    return result


def load_psortb(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load PSORTb Phase-1 calls.json → {protein_id_wp: {field: value}}.

    The artifact is a dict keyed by RefSeq WP_ accession (== gene_mapping.protein_id),
    each value carrying {localization, score, secondary_localization, secondary_score,
    is_multi_localized, is_unknown}. Returned verbatim so the field resolvers can pull
    `localization` / `score` by name. Missing file → {} (strain not yet PSORTb-run).
    """
    path = os.path.join(data_dir, "psortb", f"{strain_name}.psortb.calls.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    # Keys are already WP_ accessions; values are flat dicts. Pass through as-is.
    return {str(k).strip(): v for k, v in data.items() if isinstance(v, dict)}


def load_signalp(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load SignalP Phase-1 calls.json → {protein_id_wp: {field: value}}.

    The artifact is a dict keyed by RefSeq WP_ accession (== gene_mapping.protein_id),
    each value carrying {signalp_type, probability, cleavage_site, cleavage_probability}.
    Returned verbatim so the field resolvers can pull fields by name. "OTHER" calls are
    kept verbatim; the adapter skips them at edge-build time. Missing file → {} (strain
    not yet SignalP-normalized — run `signalp-run --normalize`).
    """
    path = os.path.join(data_dir, "signalp", f"{strain_name}.signalp.calls.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    # Keys are already WP_ accessions; values are flat dicts. Pass through as-is.
    return {str(k).strip(): v for k, v in data.items() if isinstance(v, dict)}


def load_interproscan(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load InterProScan calls.json (faceted format) → {protein_id_wp: {field: value}}.

    The artifact is a dict keyed by RefSeq WP_ accession (== gene_mapping.protein_id),
    with a per-protein shape of
    `{md5, match_count, libraries: {LIB: [{accession, name, ipr, start, end,
    evalue, score}, ...]}, interpro_entries: {IPR: {...}}, go_terms: {GO: [IPR, ...]}}`
    (accessions already version-stripped by the parser). We surface only the light
    per-gene summary the merge needs:
    - `interpro_entries` — distinct InterPro entry ids (drives contributing_sources
      + the DataSource node + Gene routing)
    - `pfam_signatures` — direct PFAM library HMM hits (a direct hit, no inference)
    - `ncbifam_ids` — direct NCBIFAM library hits
    - `hamap_descriptions` — HAMAP library match names (HAMAP has no stable
      accession→description reference, so the name is carried here)
    - `go_term_donors` — GO term → contributing InterPro entry id(s), as reported
      by InterProScan's own GO attribution

    The rich per-match evidence (coordinates / e-value / score / libraries) is NOT
    merged here; the interpro_adapter reads this same calls.json directly at
    KG-build time (like tcdb_adapter reads tcdb_pruned.json). All fields are sparse
    (omitted when empty); a protein with nothing to surface is dropped entirely.
    Missing file → {} (strain not yet InterProScan-run).
    See docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md.
    """
    path = os.path.join(data_dir, "interproscan", f"{strain_name}.interproscan.calls.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    result: dict[str, dict] = {}
    for wp, call in data.items():
        if not isinstance(call, dict):
            continue
        entries = sorted(call.get("interpro_entries") or {})
        libs = call.get("libraries") or {}
        pfam_sigs = sorted({
            r["accession"]
            for r in libs.get("PFAM", [])
            if (r.get("accession") or "").startswith("PF")
        })
        ncbifam = sorted({
            r["accession"]
            for r in libs.get("NCBIFAM", [])
            if r.get("accession")
        })
        hamap = sorted({
            r["name"]
            for r in libs.get("HAMAP", [])
            if r.get("name")
        })
        donors = call.get("go_terms") or {}
        if not (entries or pfam_sigs or ncbifam or hamap or donors):
            continue
        row: dict = {}
        if entries:
            row["interpro_entries"] = entries
        if pfam_sigs:
            row["pfam_signatures"] = pfam_sigs
        if ncbifam:
            row["ncbifam_ids"] = ncbifam
        if hamap:
            row["hamap_descriptions"] = hamap
        if donors:
            row["go_term_donors"] = donors
        result[str(wp).strip()] = row
    return result


def load_tcdb_diamond(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load tcdb-diamond Phase-1 calls.json → {protein_id_wp: {tcdb_ids: [...]}}.

    The artifact is a dict keyed by RefSeq WP_ accession (== gene_mapping.protein_id).
    We surface only the distinct TC ids this protein was called at — the per-call
    evidence (tier / confidence_score / identity / evalue / consensus_n) is NOT
    merged; `tcdb_adapter` reads this same calls.json directly at KG-build time
    for edge properties, exactly as `interpro_adapter` does.

    ALL candidates are surfaced. The tier policy in `classify_hit` is already the
    quality gate (e-value <= 0.001, HSP >= 50 aa, coverage floors) and is
    sibling-independent; the old post-hoc `filter_action` chain was removed in
    07357ac0 because its verdicts depended on unrelated sibling candidates.
    Weak calls are gated downstream instead — post-import folds 'tcdb' into
    annotation_types only for eggNOG-sourced or tier<=2 edges.

    Proteins with no call are dropped. Missing file → {} (strain not yet
    tcdb-diamond-run). See
    docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md.
    """
    path = os.path.join(data_dir, "tcdb", f"{strain_name}.tcdb.calls.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    result: dict[str, dict] = {}
    for wp, rec in data.items():
        if not isinstance(rec, dict):
            continue
        tcids = sorted({
            c["tcid"] for c in rec.get("calls", [])
            if isinstance(c, dict) and c.get("tcid")
        })
        if tcids:
            result[str(wp).strip()] = {"tcdb_ids": tcids}
    return result


def load_merops(data_dir: str, strain_name: str) -> dict[str, dict]:
    """Load merops-diamond Phase-1 calls.json → {protein_id_wp: {merops_ids: [...]}}.

    The artifact is a dict keyed by RefSeq WP_ accession (== gene_mapping.protein_id).
    We surface only the distinct called codes (family / subfamily / identifier,
    already tier-truncated by the runner) — the per-call evidence (tier /
    confidence_score / identity / evalue / call_class inputs) is NOT merged;
    `merops_adapter` reads this same calls.json directly at KG-build time for
    edge properties, exactly as `tcdb_adapter` / `interpro_adapter` do.

    ALL candidates are surfaced — the tier policy in
    `multiomics_kg/utils/merops_diamond.py` is already the quality gate; weak
    (tier-3) calls are gated downstream (post-import folds 'merops' into
    annotation_types only for tier<=2 edges). Proteins with no call are
    dropped. Missing file → {} (strain not yet merops-diamond-run). See
    docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md.
    """
    path = os.path.join(data_dir, "merops", f"{strain_name}.merops.calls.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    result: dict[str, dict] = {}
    for wp, rec in data.items():
        if not isinstance(rec, dict):
            continue
        codes = sorted({
            c["code"] for c in rec.get("calls", [])
            if isinstance(c, dict) and c.get("code")
        })
        if codes:
            result[str(wp).strip()] = {"merops_ids": codes}
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
        ps: dict | None = None,
        sp: dict | None = None,
        ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None,
    ) -> Any:
        """Fetch raw value from source row according to src_cfg spec."""
        ps = ps or {}
        sp = sp or {}
        ipr = ipr or {}
        tcd = tcd or {}
        mer = mer or {}
        source = src_cfg.get("source", "")
        field = src_cfg.get("field", "")

        if source == "gene_mapping":
            raw = gm.get(field)
        elif source == "eggnog":
            raw = eg.get(field)
        elif source == "uniprot":
            raw = up.get(field)
        elif source == "psortb":
            raw = ps.get(field)
        elif source == "signalp":
            raw = sp.get(field)
        elif source == "interproscan":
            raw = ipr.get(field)
        elif source == "tcdb_diamond":
            raw = tcd.get(field)
        elif source == "merops_diamond":
            raw = mer.get(field)
        else:
            return None

        # URL-decode string values (gene_mapping uses URL encoding)
        if isinstance(raw, str):
            raw = unquote(raw.strip())

        # allow_dash: keep a lone '-' (e.g. minus strand) instead of treating
        # it as the eggNOG "no value" sentinel.
        if not _nonempty(raw, keep_dash=bool(src_cfg.get("allow_dash", False))):
            return None

        return raw

    # ── apply transform to a single value ─────────────────────────────────────

    def _apply_transform(self, transform: str | None, value: Any) -> Any:
        """Apply a named transform to a value; returns empty string on failure."""
        if not transform or transform not in _TRANSFORMS:
            return value
        fn = _TRANSFORMS[transform]
        if isinstance(value, list):
            return [r for r in (fn(v) for v in value if _nonempty(v)) if r is not None]
        return fn(value)

    # ── resolver: passthrough ──────────────────────────────────────────────────

    def _resolve_passthrough(
        self, fconf: dict, gm: dict, eg: dict, up: dict, ps: dict | None = None,
        sp: dict | None = None, ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None
    ) -> Any:
        keep_dash = bool(fconf.get("allow_dash", False))
        raw = self._get_raw(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
        if not _nonempty(raw, keep_dash=keep_dash):
            return None
        transform = fconf.get("transform")
        if transform:
            raw = self._apply_transform(transform, raw)
        return raw if _nonempty(raw, keep_dash=keep_dash) else None

    # ── resolver: passthrough_list ─────────────────────────────────────────────

    def _resolve_passthrough_list(
        self, fconf: dict, gm: dict, eg: dict, up: dict, ps: dict | None = None,
        sp: dict | None = None, ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None
    ) -> list[str] | None:
        raw = self._get_raw(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
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
        ps: dict,
        source_tracking: dict,
        locus_tag: str = "",
        sp: dict | None = None,
        ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None,
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
            raw = self._get_raw(cand, gm, eg, up, ps, sp, ipr, tcd, mer)
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
        self, fconf: dict, gm: dict, eg: dict, up: dict, ps: dict | None = None,
        sp: dict | None = None, ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None,
        source_tracking: dict | None = None,
    ) -> list[str] | None:
        """Merge tokens from all sources, deduplicate, apply global filter.

        When *fconf* declares ``track_source`` and a *source_tracking* dict is
        passed, records a **per-token provenance map** ``{token: [source_label,
        …]}`` under that key — the honest shape for multi-source union fields
        (a GO term contributed by both UniProt and InterPro carries both). The
        label is each source's ``source_label`` (default: its ``source``); the
        map keys are exactly the surviving (post-filter) tokens.
        """
        global_filter = fconf.get("filter")
        global_filter_not = fconf.get("filter_not")
        track_key = fconf.get("track_source")
        seen: dict[str, None] = {}  # ordered set
        token_sources: dict[str, set[str]] = {}

        for src_cfg in fconf.get("sources", []):
            raw = self._get_raw(src_cfg, gm, eg, up, ps, sp, ipr, tcd, mer)
            if not _nonempty(raw):
                continue

            src_label = src_cfg.get("source_label", src_cfg.get("source", ""))
            delimiter = src_cfg.get("delimiter", ",")
            transform = src_cfg.get("transform")

            # Special-case transforms that produce lists from a string
            if transform == "extract_go_from_pipe":
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
                if tok is None:
                    continue
                tok = str(tok).strip()
                if not tok or tok == "-":
                    continue
                if global_filter and not re.match(global_filter, tok):
                    continue
                if global_filter_not and re.match(global_filter_not, tok):
                    continue
                seen[tok] = None
                if track_key and src_label:
                    token_sources.setdefault(tok, set()).add(src_label)

        result = list(seen.keys())
        if track_key and source_tracking is not None and token_sources:
            source_tracking[track_key] = {
                tok: sorted(token_sources[tok]) for tok in result if tok in token_sources
            }
        return result if result else None

    # ── resolver: integer / float ──────────────────────────────────────────────

    def _resolve_integer(
        self, fconf: dict, gm: dict, eg: dict, up: dict, ps: dict | None = None,
        sp: dict | None = None, ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None
    ) -> int | None:
        raw = self._get_raw(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
        if raw is None:
            return None
        try:
            return int(float(str(raw).strip()))
        except (ValueError, TypeError):
            return None

    def _resolve_float(
        self, fconf: dict, gm: dict, eg: dict, up: dict, ps: dict | None = None,
        sp: dict | None = None, ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None
    ) -> float | None:
        raw = self._get_raw(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
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
        ps: dict | None = None,
        sp: dict | None = None,
        ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None,
    ) -> dict:
        """All source fields, source-prefixed — full audit trail."""
        ps = ps or {}
        sp = sp or {}
        ipr = ipr or {}
        tcd = tcd or {}
        mer = mer or {}
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
        for k, v in ps.items():
            if _nonempty(v):
                wide[f"psortb_{k}"] = v
        for k, v in sp.items():
            if _nonempty(v):
                wide[f"signalp_{k}"] = v
        for k, v in ipr.items():
            if _nonempty(v):
                wide[f"interproscan_{k}"] = v
        for k, v in tcd.items():
            if _nonempty(v):
                wide[f"tcdb_diamond_{k}"] = v
        for k, v in mer.items():
            if _nonempty(v):
                wide[f"merops_diamond_{k}"] = v
        return wide

    # ── build merged ──────────────────────────────────────────────────────────

    def build_merged(
        self,
        gm: dict,
        eg: dict,
        up: dict,
        ps: dict | None = None,
        organism_name: str | None = None,
        sp: dict | None = None,
        ipr: dict | None = None, tcd: dict | None = None,
        mer: dict | None = None,
    ) -> dict:
        """Apply merge rules → canonical field set."""
        ps = ps or {}
        sp = sp or {}
        ipr = ipr or {}
        tcd = tcd or {}
        mer = mer or {}
        result: dict[str, Any] = {}
        source_tracking: dict[str, str] = {}
        locus_tag = gm.get("locus_tag", "")

        for canonical_field, fconf in self.field_configs.items():
            ftype = fconf.get("type", "passthrough")

            if ftype == "single":
                val = self._resolve_single(fconf, gm, eg, up, ps, source_tracking, locus_tag, sp=sp, ipr=ipr, tcd=tcd, mer=mer)
            elif ftype == "union":
                val = self._resolve_union(fconf, gm, eg, up, ps, sp, ipr, tcd, mer, source_tracking=source_tracking)
            elif ftype == "passthrough":
                val = self._resolve_passthrough(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
            elif ftype == "passthrough_list":
                val = self._resolve_passthrough_list(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
            elif ftype == "integer":
                val = self._resolve_integer(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
            elif ftype == "float":
                val = self._resolve_float(fconf, gm, eg, up, ps, sp, ipr, tcd, mer)
            elif ftype == "extract_first_match":
                val = extract_first_match_in_sources(
                    fconf.get("sources", []), gm, eg, up,
                    fconf.get("pattern", ""),
                    fconf.get("extract_group", 0),
                )
            else:
                continue

            if _nonempty(val, keep_dash=bool(fconf.get("allow_dash", False))):
                result[canonical_field] = val

        # Split comma-joined tokens in synonym lists (e.g. UniProt "pds,crtD"
        # arrives as a single token because the source uses space delimiter)
        for field in ("gene_synonyms", "gene_name_synonyms"):
            if field in result:
                expanded = []
                for tok in result[field]:
                    if "," in tok:
                        expanded.extend(t.strip() for t in tok.split(",") if t.strip())
                    else:
                        expanded.append(tok)
                result[field] = list(dict.fromkeys(expanded))  # dedupe, preserve order

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

        # Compute gene_category — high-level functional classification
        category = _compute_gene_category(result)
        assert category in VALID_CATEGORIES, (
            f"Invalid gene_category {category!r} for {result.get('locus_tag')}"
        )
        result["gene_category"] = category

        # Compute contributing_sources (F2 ship 2.4)
        result["contributing_sources"] = _compute_contributing_sources(result)

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
        # Pfam descriptions are added post-merge by enrich_pfam_fields()
        # (needs pfam_data reference + enriched pfam_ids)

        if alt_descriptions:
            result["alternate_functional_descriptions"] = alt_descriptions

        # ── Computed fields for MCP gene lookup ──────────────────────────────

        # organism_name — preferred organism name (e.g., "Prochlorococcus MED4")
        if organism_name:
            result["organism_name"] = organism_name

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
            result.get("uniprot_accession"),
        ]))
        all_ids.update(result.get("old_locus_tags") or [])
        all_ids.update(result.get("alternative_locus_tags") or [])
        all_ids.update(result.get("gene_name_synonyms") or [])
        all_ids -= scalar_indexed
        if all_ids:
            result["all_identifiers"] = sorted(all_ids)

        return result


# ─── per-strain pipeline ──────────────────────────────────────────────────────

def _write_metabolism_report(strain_name: str, data_dir: Path,
                              raw_eggnog_path: Path,
                              merged: dict) -> None:
    """Generate step2_metabolism_report.json for one strain.

    Counts raw values (from eggNOG TSV cols 15/18/19) vs kept/validated
    values (from the merged JSON's kegg_reactions / transporter_classification /
    cazy_ids fields).

    Note (Spec 1.2): kegg_reactions now stores raw KEGG R-numbers directly
    (no MNX resolution at build time).
    """
    raw_kegg: list[str] = []
    raw_tc: list[str] = []
    raw_cazy: list[str] = []

    if raw_eggnog_path.exists():
        with open(raw_eggnog_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= 18:
                    continue
                for col_idx, target in ((14, raw_kegg), (17, raw_tc), (18, raw_cazy)):
                    cell = cols[col_idx] if col_idx < len(cols) else ""
                    if cell in ("", "-"):
                        continue
                    for tok in cell.split(","):
                        tok = tok.strip()
                        if tok:
                            target.append(tok)

    kept_r_numbers: set[str] = set()
    validated_tc: set[str] = set()
    validated_cazy: set[str] = set()
    for gene in merged.values():
        kept_r_numbers.update(gene.get("kegg_reactions", []))
        validated_tc.update(gene.get("transporter_classification", []))
        validated_cazy.update(gene.get("cazy_ids", []))

    invalid_tc = sorted(set(raw_tc) - validated_tc)[:10]
    invalid_cazy = sorted(set(raw_cazy) - validated_cazy)[:10]

    report = {
        "strain": strain_name,
        "gene_count": len(merged),
        "kegg_reactions": {
            "raw_total":             len(raw_kegg),
            "raw_unique":            len(set(raw_kegg)),
            "kept_total":            sum(len(g.get("kegg_reactions", [])) for g in merged.values()),
            "kept_unique_r_numbers": len(kept_r_numbers),
        },
        "transporter_classification": {
            "raw_total":         len(raw_tc),
            "raw_unique":        len(set(raw_tc)),
            "validated_total":   sum(len(g.get("transporter_classification", [])) for g in merged.values()),
            "validated_unique":  len(validated_tc),
            "invalid_examples":  invalid_tc,
        },
        "cazy_ids": {
            "raw_total":         len(raw_cazy),
            "raw_unique":        len(set(raw_cazy)),
            "validated_total":   sum(len(g.get("cazy_ids", [])) for g in merged.values()),
            "validated_unique":  len(validated_cazy),
            "invalid_examples":  invalid_cazy,
        },
    }

    out_path = data_dir / "step2_metabolism_report.json"
    out_path.write_text(json.dumps(report, indent=2))


def process_strain(
    row: dict,
    config: dict,
    force: bool = False,
    pfam_data: PfamData | None = None,
    interpro_ref: dict | None = None,
    ncbifam_ref: dict | None = None,
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

    # Defensive check: seqid column was added by Task 2 (db065e0c); older caches lack it.
    # If missing, Gene.contig will be null until the cache is regenerated.
    first_row = next(iter(gm_data.values()), {})
    if gm_data and "seqid" not in first_row:
        logger.warning(
            f"gene_mapping.csv for {strain_name} lacks 'seqid' column — "
            f"cache is stale; re-run `bash scripts/prepare_data.sh --steps 0 --force`. "
            f"Gene.contig will be null for this strain until regenerated."
        )

    eg_data = load_eggnog(data_dir, strain_name)
    up_data: dict[str, dict] = {}
    if ncbi_taxon_id:
        up_data = load_uniprot(data_dir, ncbi_taxon_id, organism_group)
    ps_data = load_psortb(data_dir, strain_name)
    sp_data = load_signalp(data_dir, strain_name)
    ipr_data = load_interproscan(data_dir, strain_name)
    tcd_data = load_tcdb_diamond(data_dir, strain_name)
    mer_data = load_merops(data_dir, strain_name)

    print(f"  gene_mapping : {len(gm_data):>5} genes")
    print(f"  eggnog       : {len(eg_data):>5} entries")
    print(f"  uniprot      : {len(up_data):>5} entries (keyed by RefSeq)")
    print(f"  psortb       : {len(ps_data):>5} entries (keyed by RefSeq)")
    print(f"  signalp      : {len(sp_data):>5} entries (keyed by RefSeq)")
    print(f"  interproscan : {len(ipr_data):>5} entries (keyed by RefSeq)")
    print(f"  tcdb_diamond : {len(tcd_data):>5} entries (keyed by RefSeq)")
    print(f"  merops_diamond: {len(mer_data):>4} entries (keyed by RefSeq)")

    builder = AnnotationBuilder(config)

    wide_out: dict[str, dict] = {}
    merged_out: dict[str, dict] = {}

    stats = dict(total=0, eggnog_hit=0, uniprot_hit=0,
                 has_product=0, has_go=0, has_cog=0, has_kegg_ko=0)

    for locus_tag, gm_row in gm_data.items():
        protein_id = (gm_row.get("protein_id") or "").strip()
        eg_row = eg_data.get(protein_id, {})
        up_row = up_data.get(protein_id, {})
        ps_row = ps_data.get(protein_id, {})
        sp_row = sp_data.get(protein_id, {})
        ipr_row = ipr_data.get(protein_id, {})
        tcd_row = tcd_data.get(protein_id, {})
        mer_row = mer_data.get(protein_id, {})

        stats["total"] += 1
        if eg_row:
            stats["eggnog_hit"] += 1
        if up_row:
            stats["uniprot_hit"] += 1

        wide_out[locus_tag] = builder.build_wide(gm_row, eg_row, up_row, ps_row, sp_row, ipr_row, tcd_row, mer_row)
        merged = builder.build_merged(gm_row, eg_row, up_row, ps_row,
                                      organism_name=preferred_name, sp=sp_row, ipr=ipr_row,
                                      tcd=tcd_row, mer=mer_row)
        # Layer B: promote InterPro/NCBIfam entry xrefs into go/ec/cazy/pfam +
        # naming recovery (before the pfam enrichment loop, so InterPro's direct
        # PF* hits get re-keyed too).
        if interpro_ref is not None:
            enrich_interpro_fields(merged, ipr_row, interpro_ref, ncbifam_ref)
        merged_out[locus_tag] = merged

        if merged.get("product"):
            stats["has_product"] += 1
        if merged.get("go_terms"):
            stats["has_go"] += 1
        if merged.get("cog_category"):
            stats["has_cog"] += 1
        if merged.get("kegg_ko"):
            stats["has_kegg_ko"] += 1
        cat = merged.get("gene_category", "Unknown")
        stats.setdefault(f"cat_{cat}", 0)
        stats[f"cat_{cat}"] += 1

    # Compute ortholog_groups for each gene (post-merge enrichment)
    og_organism_group = organism_group_from_path(data_dir)
    og_stats = {"has_og": 0, "cyanorak": 0, "eggnog_bacteria": 0, "eggnog_intermediate": 0, "eggnog_family": 0}
    for locus_tag, gene in merged_out.items():
        ogs = extract_ortholog_groups(gene, og_organism_group)
        if ogs:
            gene["ortholog_groups"] = ogs
            og_stats["has_og"] += 1
            for og in ogs:
                if og["source"] == "cyanorak":
                    og_stats["cyanorak"] += 1
                elif og["taxon_id"] == 2:
                    og_stats["eggnog_bacteria"] += 1
                elif og["specificity_rank"] == 2:
                    og_stats["eggnog_intermediate"] += 1
                else:
                    og_stats["eggnog_family"] += 1
    n = stats["total"] or 1
    print(f"\n  === {strain_name} ortholog groups ===")
    print(f"  Genes with OG:        {og_stats['has_og']} ({100 * og_stats['has_og'] // n}%)")
    print(f"  Cyanorak memberships: {og_stats['cyanorak']}")
    print(f"  EggNOG bacteria:     {og_stats['eggnog_bacteria']}")
    print(f"  EggNOG intermediate: {og_stats['eggnog_intermediate']}")
    print(f"  EggNOG family:       {og_stats['eggnog_family']}")

    # Enrich pfam_ids: resolve shortnames to PF* accessions, drop non-Pfam tokens
    if pfam_data is not None:
        pfam_stats = {"has_pfam": 0, "resolved_shortnames": 0, "genes_before": 0}
        all_unresolved: dict[str, int] = {}  # token -> gene count
        for locus_tag, gene in merged_out.items():
            had_raw = bool(gene.get("pfam_ids"))
            if had_raw:
                pfam_stats["genes_before"] += 1
            unresolved = enrich_pfam_fields(gene, pfam_data)
            if gene.get("pfam_ids"):
                pfam_stats["has_pfam"] += 1
            # Count shortnames that were resolved (genes_after > genes with only PF* before)
            for tok in unresolved:
                all_unresolved[tok] = all_unresolved.get(tok, 0) + 1

        print(f"\n  === {strain_name} Pfam enrichment ===")
        print(f"  Genes with raw pfam tokens: {pfam_stats['genes_before']}")
        print(f"  Genes with clean PF* IDs:   {pfam_stats['has_pfam']}")
        if all_unresolved:
            top_unresolved = sorted(all_unresolved.items(), key=lambda x: -x[1])[:10]
            total_unresolved = sum(all_unresolved.values())
            print(f"  Unresolved tokens: {len(all_unresolved)} unique, {total_unresolved} total occurrences")
            for tok, count in top_unresolved:
                print(f"    {tok}: {count} genes")
        else:
            print(f"  Unresolved tokens: 0")

    with open(wide_path, "w") as f:
        json.dump(wide_out, f, indent=2, sort_keys=True)
    print(f"  → {wide_path}")

    with open(merged_path, "w") as f:
        json.dump(merged_out, f, indent=2, sort_keys=True)
    print(f"  → {merged_path}")

    _write_metabolism_report(
        strain_name=strain_name,
        data_dir=Path(data_dir),
        raw_eggnog_path=Path(data_dir) / "eggnog" / f"{strain_name}.emapper.annotations",
        merged=merged_out,
    )

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

    # Load Pfam reference data once for all strains
    cache_root = PROJECT_ROOT / "cache" / "data"
    pfam_data = load_pfam_data(cache_root)
    print(f"Pfam reference: {len(pfam_data.by_accession)} entries, "
          f"{len(pfam_data.by_shortname)} shortnames, {len(pfam_data.clans)} clans")

    # Load the InterPro reference once (lazy — returns the committed cache unless
    # --force/--refetch-raw; mirrors the Pfam precedent, no prepare_data renumber).
    interpro_ref = build_interpro_reference.build()
    n_ec = sum("ec_numbers" in m for m in interpro_ref.values())
    n_cazy = sum("cazy_ids" in m for m in interpro_ref.values())
    print(f"InterPro reference: {len(interpro_ref)} entries "
          f"({n_ec} with EC, {n_cazy} with CAZy)")

    # Load the NCBIfam reference once (lazy — reads the committed cache directly;
    # `{}` when missing so merge never blocks on a network fetch here).
    ncbifam_ref_path = cache_root / "ncbifam" / "ncbifam_reference.json"
    if ncbifam_ref_path.exists():
        with open(ncbifam_ref_path, encoding="utf-8") as fh:
            ncbifam_ref = json.load(fh)
    else:
        ncbifam_ref = {}
    print(f"NCBIfam reference: {len(ncbifam_ref)} entries")

    print(f"Processing {len(rows)} strain(s) with config: {args.config}")
    for row in rows:
        process_strain(row, config, force=args.force, pfam_data=pfam_data,
                       interpro_ref=interpro_ref, ncbifam_ref=ncbifam_ref)

    if args.llm_summary:
        print("\nLLM summary generation (Step 1C) not yet implemented.")


if __name__ == "__main__":
    main()
