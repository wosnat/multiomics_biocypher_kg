"""
KEGG REST API utilities: download and cache the data needed to build the KEGG
4-level hierarchy (KO → Pathway → Subcategory → Category).

All network calls use the public KEGG REST API (no auth required):
  https://rest.kegg.jp/

Two-level cache:
  <cache_root>/kegg/raw/          — raw REST responses (one file per endpoint)
  <cache_root>/kegg/kegg_data.json — parsed/aggregated dict (fast path)

The raw cache avoids re-downloading from KEGG when only the parser changes.
Pass ``force=True`` to re-download raw files and rebuild kegg_data.json.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

_KEGG_BASE = "https://rest.kegg.jp"
_BRITE_KO_URL = f"{_KEGG_BASE}/get/br:ko00001/json"
_KO_LIST_URL = f"{_KEGG_BASE}/list/ko"
_KO_PATHWAY_LINK_URL = f"{_KEGG_BASE}/link/pathway/ko"
_PATHWAY_KO_LIST_URL = f"{_KEGG_BASE}/list/pathway/ko"
_REACTION_LIST_URL = f"{_KEGG_BASE}/list/reaction"
_COMPOUND_LIST_URL = f"{_KEGG_BASE}/list/compound"
_LINK_COMPOUND_REACTION_URL = f"{_KEGG_BASE}/link/compound/reaction"
_LINK_PATHWAY_REACTION_URL = f"{_KEGG_BASE}/link/pathway/reaction"
_LINK_PATHWAY_COMPOUND_URL = f"{_KEGG_BASE}/link/pathway/compound"

_TIMEOUT = 120  # seconds per request

# Maps each REST URL to its raw-cache filename under <cache_root>/kegg/raw/
_RAW_FILES = {
    _KO_LIST_URL:                "list_ko.txt",
    _KO_PATHWAY_LINK_URL:        "link_pathway_ko.txt",
    _PATHWAY_KO_LIST_URL:        "list_pathway_ko.txt",
    _REACTION_LIST_URL:          "list_reaction.txt",
    _COMPOUND_LIST_URL:          "list_compound.txt",
    _LINK_COMPOUND_REACTION_URL: "link_compound_reaction.txt",
    _LINK_PATHWAY_REACTION_URL:  "link_pathway_reaction.txt",
    _LINK_PATHWAY_COMPOUND_URL:  "link_pathway_compound.txt",
    _BRITE_KO_URL:               "br_ko00001.json",
}


def _fetch_text(url: str) -> str:
    """Fetch plain text from a KEGG REST endpoint."""
    logger.info(f"Fetching {url}")
    resp = requests.get(url, timeout=_TIMEOUT)
    resp.raise_for_status()
    return resp.text


def _fetch_json(url: str) -> dict:
    """Fetch JSON from a KEGG REST endpoint."""
    logger.info(f"Fetching {url}")
    resp = requests.get(url, timeout=_TIMEOUT)
    resp.raise_for_status()
    return resp.json()


def _fetch_text_cached(url: str, raw_dir: Path, force: bool = False) -> str:
    """Fetch plain text from a KEGG REST endpoint, caching the response.

    If the cached file exists and ``force=False``, return its contents without
    hitting the network.  Otherwise download, save, and return.
    """
    filename = _RAW_FILES[url]
    cache_path = raw_dir / filename
    if cache_path.exists() and not force:
        logger.info(f"Reading cached {filename}")
        return cache_path.read_text(encoding="utf-8")
    text = _fetch_text(url)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(text, encoding="utf-8")
    return text


def _fetch_json_cached(url: str, raw_dir: Path, force: bool = False) -> dict:
    """Fetch JSON from a KEGG REST endpoint, caching the response.

    If the cached file exists and ``force=False``, return its contents without
    hitting the network.  Otherwise download, save, and return.
    """
    filename = _RAW_FILES[url]
    cache_path = raw_dir / filename
    if cache_path.exists() and not force:
        logger.info(f"Reading cached {filename}")
        return json.loads(cache_path.read_text(encoding="utf-8"))
    data = _fetch_json(url)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    return data


def _parse_ko_names(text: str) -> dict[str, str]:
    """Parse `list/ko` response into {K#####: name_str}."""
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.strip().split("\t", 1)
        if len(parts) != 2:
            continue
        # key like "ko:K02338"
        raw_id, name = parts
        ko_id = raw_id.removeprefix("ko:")
        if ko_id.startswith("K"):
            result[ko_id] = name.strip()
    logger.info(f"Parsed {len(result)} KO names")
    return result



def _parse_pathway_ko_names(text: str) -> dict[str, str]:
    """Parse `list/pathway/ko` response into {ko#####: name_str}.

    Unlike `list/pathway` (which returns map-prefixed IDs), this endpoint
    returns ko-prefixed reference pathway IDs and covers global/overview maps
    (ko01100–ko01320) that are absent from the BRITE hierarchy.
    """
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.strip().split("\t", 1)
        if len(parts) != 2:
            continue
        # key like "path:ko00230"
        raw_id, name = parts
        pw_id = raw_id.removeprefix("path:")
        if pw_id.startswith("ko"):
            result[pw_id] = name.strip()
    logger.info(f"Parsed {len(result)} pathway names from list/pathway/ko")
    return result


def _parse_ko_to_pathways(text: str) -> dict[str, list[str]]:
    """Parse `link/pathway/ko` response into {K#####: [ko#####, ...]}."""
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.strip().split("\t", 1)
        if len(parts) != 2:
            continue
        # "ko:K02338\tpath:ko03030"
        raw_ko, raw_pw = parts
        ko_id = raw_ko.removeprefix("ko:")
        pw_id = raw_pw.removeprefix("path:")
        if not ko_id.startswith("K") or not pw_id.startswith("ko"):
            continue
        result.setdefault(ko_id, []).append(pw_id)
    logger.info(f"Parsed KO→Pathway links for {len(result)} KOs")
    return result


def _parse_reaction_names(text: str) -> dict[str, str]:
    """Parse `/list/reaction` response into {R#####: name_str}."""
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_id, name = parts
        rxn_id = raw_id.removeprefix("rn:")
        if rxn_id.startswith("R") and rxn_id[1:].isdigit():
            result[rxn_id] = name.strip()
    logger.info(f"Parsed {len(result)} reaction names")
    return result


def _parse_compound_names(text: str) -> dict[str, str]:
    """Parse `/list/compound` into {C#####: first_synonym}.

    KEGG returns semicolon-separated synonyms; we keep the first as canonical.
    """
    result: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_id, names = parts
        cpd_id = raw_id.removeprefix("cpd:")
        if cpd_id.startswith("C") and cpd_id[1:].isdigit():
            first = names.split(";", 1)[0].strip()
            result[cpd_id] = first
    logger.info(f"Parsed {len(result)} compound names")
    return result


def _parse_reaction_to_compounds(text: str) -> dict[str, list[str]]:
    """Parse `/link/compound/reaction` into {R#####: [C#####, ...]}.

    KEGG returns lines in the form: `rn:R00200\tcpd:C00074`.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_rxn, raw_cpd = parts
        rxn_id = raw_rxn.removeprefix("rn:")
        cpd_id = raw_cpd.removeprefix("cpd:")
        if not (rxn_id.startswith("R") and cpd_id.startswith("C")):
            continue
        result.setdefault(rxn_id, []).append(cpd_id)
    logger.info(f"Parsed compound-reaction links for {len(result)} reactions")
    return result


def _parse_reaction_to_pathways(text: str) -> dict[str, list[str]]:
    """Parse `/link/pathway/reaction` into {R#####: [ko#####, ...]}.

    KEGG returns rn-prefixed pathway IDs (e.g. `path:rn00010`); we rewrite them
    to ko-prefixed form so they match existing KeggTerm pathway node IDs.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_rxn, raw_pw = parts
        rxn_id = raw_rxn.removeprefix("rn:")
        pw_id = raw_pw.removeprefix("path:")
        if not rxn_id.startswith("R"):
            continue
        # Normalize rn00010 → ko00010 (the form used by existing KeggTerm nodes)
        if pw_id.startswith("rn"):
            pw_id = "ko" + pw_id[2:]
        if not pw_id.startswith("ko"):
            continue
        result.setdefault(rxn_id, []).append(pw_id)
    logger.info(f"Parsed reaction-pathway links for {len(result)} reactions")
    return result


def _parse_compound_to_pathways(text: str) -> dict[str, list[str]]:
    """Parse `/link/pathway/compound` into {C#####: [ko#####, ...]}.

    KEGG returns map-prefixed pathway IDs; we rewrite to ko-prefixed for
    consistency with the rest of the KG.
    """
    result: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        raw_cpd, raw_pw = parts
        cpd_id = raw_cpd.removeprefix("cpd:")
        pw_id = raw_pw.removeprefix("path:")
        if not cpd_id.startswith("C"):
            continue
        if pw_id.startswith("map"):
            pw_id = "ko" + pw_id[3:]
        if not pw_id.startswith("ko"):
            continue
        result.setdefault(cpd_id, []).append(pw_id)
    logger.info(f"Parsed compound-pathway links for {len(result)} compounds")
    return result


def _parse_brite_hierarchy(brite_json: dict) -> tuple[
    dict[str, str],  # pathway_to_subcategory: {ko#####: subcat_code}
    dict[str, str],  # subcategory_names: {subcat_code: name}
    dict[str, str],  # subcategory_to_category: {subcat_code: cat_code}
    dict[str, str],  # category_names: {cat_code: name}
    dict[str, str],  # pathway_names: {ko#####: name_str}
]:
    """
    Parse the `br:ko00001` BRITE JSON hierarchy.

    Structure (4 levels):
      A (children of root): "09100 Metabolism"           → KeggCategory
      B (children of A):    "09101 Carbohydrate ..."     → KeggSubcategory
      C (children of B):    "00010 Glycolysis [PATH:ko00010]"  → extract ko#####
      D (children of C):    individual KOs (K#####) — ignored here

    Returns five dicts for building the Pathway→Subcategory→Category hierarchy,
    plus pathway names extracted from C-level node labels.
    """
    pathway_to_subcat: dict[str, str] = {}
    subcat_names: dict[str, str] = {}
    subcat_to_cat: dict[str, str] = {}
    cat_names: dict[str, str] = {}
    pw_names: dict[str, str] = {}

    # The top-level object has a "children" key with A-level entries
    a_entries = brite_json.get("children", [])

    for a_node in a_entries:
        a_name = a_node.get("name", "")
        # e.g. "09100 Metabolism"
        a_match = re.match(r"^(\d{5})\s+(.+)$", a_name)
        if not a_match:
            continue
        cat_code = a_match.group(1)
        cat_label = a_match.group(2).strip()
        cat_names[cat_code] = cat_label

        for b_node in a_node.get("children", []):
            b_name = b_node.get("name", "")
            b_match = re.match(r"^(\d{5})\s+(.+)$", b_name)
            if not b_match:
                continue
            subcat_code = b_match.group(1)
            subcat_label = b_match.group(2).strip()
            subcat_names[subcat_code] = subcat_label
            subcat_to_cat[subcat_code] = cat_code

            for c_node in b_node.get("children", []):
                c_name = c_node.get("name", "")
                # e.g. "00010 Glycolysis / Gluconeogenesis [PATH:ko00010]"
                pw_match = re.search(r"\[PATH:(ko\d+)\]", c_name)
                if not pw_match:
                    continue
                pw_id = pw_match.group(1)
                pathway_to_subcat[pw_id] = subcat_code
                # Extract pathway name from C-level label, e.g.
                # "00010 Glycolysis / Gluconeogenesis [PATH:ko00010]"
                name_match = re.match(r"^\d{5}\s+(.+?)\s*\[PATH:", c_name)
                if name_match:
                    pw_names[pw_id] = name_match.group(1).strip()

    logger.info(
        f"BRITE hierarchy parsed: {len(cat_names)} categories, "
        f"{len(subcat_names)} subcategories, "
        f"{len(pathway_to_subcat)} pathway→subcategory links, "
        f"{len(pw_names)} pathway names"
    )
    return pathway_to_subcat, subcat_names, subcat_to_cat, cat_names, pw_names


def download_kegg_data(cache_root: Path, force: bool = False) -> dict:
    """
    Download and cache KEGG hierarchy data needed for the 4-level KEGG graph,
    plus metabolism endpoints (reactions, compounds, and their associations).

    Returns a dict with keys:
      - ``ko_names``: ``{K#####: description_str}``
      - ``pathway_names``: ``{ko#####: name_str}``
      - ``ko_to_pathways``: ``{K#####: [ko#####, ...]}``
      - ``pathway_to_subcategory``: ``{ko#####: subcat_code_str}``
      - ``subcategory_names``: ``{subcat_code: name_str}``
      - ``subcategory_to_category``: ``{subcat_code: cat_code_str}``
      - ``category_names``: ``{cat_code: name_str}``
      - ``reaction_names``: ``{R#####: description_str}`` (Spec 1.2)
      - ``compound_names``: ``{C#####: name_str}`` (Spec 1.2)
      - ``reaction_to_compounds``: ``{R#####: [C#####, ...]}`` (Spec 1.2)
      - ``reaction_to_pathways``: ``{R#####: [ko#####, ...]}`` (Spec 1.2)
      - ``compound_to_pathways``: ``{C#####: [ko#####, ...]}`` (Spec 1.2)

    Two-level cache:
      - Fast path: if ``<cache_root>/kegg/kegg_data.json`` exists and
        ``force=False``, load and return it immediately.
      - Raw cache: each endpoint's response is stored under
        ``<cache_root>/kegg/raw/``.  When ``kegg_data.json`` is missing but
        raw files exist, the parser reads from disk without hitting the network.
      - ``force=True`` re-downloads all raw files and rebuilds ``kegg_data.json``.
    """
    cache_dir = Path(cache_root) / "kegg"
    cache_file = cache_dir / "kegg_data.json"
    raw_dir = cache_dir / "raw"

    if cache_file.exists() and not force:
        logger.info(f"Loading KEGG data from cache: {cache_file}")
        with open(cache_file, encoding="utf-8") as fh:
            return json.load(fh)

    logger.info("Building KEGG data (downloads from REST if raw cache is missing) ...")
    cache_dir.mkdir(parents=True, exist_ok=True)

    ko_names = _parse_ko_names(_fetch_text_cached(_KO_LIST_URL, raw_dir, force))
    ko_to_pathways = _parse_ko_to_pathways(_fetch_text_cached(_KO_PATHWAY_LINK_URL, raw_dir, force))
    brite_json = _fetch_json_cached(_BRITE_KO_URL, raw_dir, force)
    pathway_to_subcat, subcat_names, subcat_to_cat, cat_names, brite_pw_names = (
        _parse_brite_hierarchy(brite_json)
    )
    # Supplement BRITE pathway names with list/pathway/ko API names.
    # Global/overview maps (ko01100–ko01320) are absent from BRITE but present
    # in the pathway list API.
    api_pw_names = _parse_pathway_ko_names(_fetch_text_cached(_PATHWAY_KO_LIST_URL, raw_dir, force))
    pathway_names = {**api_pw_names, **brite_pw_names}  # BRITE wins on overlap

    # ── Spec 1.2 metabolism endpoints ─────────────────────────────────
    reaction_names = _parse_reaction_names(_fetch_text_cached(_REACTION_LIST_URL, raw_dir, force))
    compound_names = _parse_compound_names(_fetch_text_cached(_COMPOUND_LIST_URL, raw_dir, force))
    reaction_to_compounds = _parse_reaction_to_compounds(_fetch_text_cached(_LINK_COMPOUND_REACTION_URL, raw_dir, force))
    reaction_to_pathways = _parse_reaction_to_pathways(_fetch_text_cached(_LINK_PATHWAY_REACTION_URL, raw_dir, force))
    compound_to_pathways = _parse_compound_to_pathways(_fetch_text_cached(_LINK_PATHWAY_COMPOUND_URL, raw_dir, force))

    data = {
        "ko_names": ko_names,
        "pathway_names": pathway_names,
        "ko_to_pathways": ko_to_pathways,
        "pathway_to_subcategory": pathway_to_subcat,
        "subcategory_names": subcat_names,
        "subcategory_to_category": subcat_to_cat,
        "category_names": cat_names,
        # Spec 1.2 — metabolism
        "reaction_names": reaction_names,
        "compound_names": compound_names,
        "reaction_to_compounds": reaction_to_compounds,
        "reaction_to_pathways": reaction_to_pathways,
        "compound_to_pathways": compound_to_pathways,
    }

    with open(cache_file, "w", encoding="utf-8") as fh:
        json.dump(data, fh, sort_keys=True)
    logger.info(f"KEGG data cached to {cache_file}")
    return data
