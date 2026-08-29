"""Build the committed NCBIfam/PGAP reference cache.

Downloads NCBI's ``hmm_PGAP.tsv`` (the PGAP HMM reference — one row per
NCBIfam/TIGRFAM family) and writes
``cache/data/ncbifam/ncbifam_reference.json``
(``{unversioned_acc: {"name", "family_type"[, "gene_symbol"][, "description"]
[, "gene_synonyms"][, "pmids"][, "ec_numbers"][, "go_terms"]}}``), consumed by
``ncbifam_adapter`` at KG-build time for NCBIfam node names/metadata. Committed
so the Docker build runs offline (InterPro/Pfam precedent).

One source file:

- ``hmm_PGAP.tsv`` (NCBI FTP, ``https://ftp.ncbi.nlm.nih.gov/hmm/current/``) —
  tab-separated, one header line, ~38.4K rows (one per PGAP HMM family).
- ``TIGR_ROLE_NAMES`` + ``TIGRFAMS_ROLE_LINK`` (NCBI FTP,
  ``https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/``) — the frozen
  2018 JCVI role archive; written to ``cache/data/ncbifam/tigr_roles.json``
  (separate file so ``ncbifam_reference.json`` keeps its flat ``{acc: …}`` shape).

Pure parsing lives in ``multiomics_kg/utils/ncbifam.py``; this module only
handles download, caching, and I/O.

Usage (module):
    uv run python -m multiomics_kg.download.build_ncbifam_reference [--force] [--refetch-raw]

- default: build the JSON from the cached raw TSV (or download it if absent).
- ``--force``: rebuild the JSON even if it already exists (from cached raw).
- ``--refetch-raw``: re-download the raw TSV first (only on an NCBIfam
  release).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path

import requests

from multiomics_kg.utils.ncbifam import (
    HYPOTH_FAMILY_TYPES,
    parse_hmm_pgap_rows,
    parse_tigr_role_link,
    parse_tigr_role_names,
)

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

HMM_PGAP_URL = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv"

CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "ncbifam"
RAW_DIR = CACHE_DIR / "raw"
HMM_PGAP_RAW = RAW_DIR / "hmm_PGAP.tsv"
REFERENCE_JSON = CACHE_DIR / "ncbifam_reference.json"

# TIGRFAMs 15.0 role archive (frozen 2018) — the only surviving family→role map.
TIGRFAMS_15_BASE = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0"
TIGR_ROLE_NAMES_URL = f"{TIGRFAMS_15_BASE}/TIGR_ROLE_NAMES"
TIGRFAMS_ROLE_LINK_URL = f"{TIGRFAMS_15_BASE}/TIGRFAMS_ROLE_LINK"
TIGR_ROLE_NAMES_RAW = RAW_DIR / "TIGR_ROLE_NAMES"
TIGRFAMS_ROLE_LINK_RAW = RAW_DIR / "TIGRFAMS_ROLE_LINK"
TIGR_ROLES_JSON = CACHE_DIR / "tigr_roles.json"
TIGR_ROLES_RELEASE = "TIGRFAMs 15.0 (frozen 2018)"


def _download(url: str, dest: Path) -> None:
    logger.info("Downloading %s → %s", url, dest)
    resp = requests.get(url, stream=True, timeout=180)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    logger.info("Wrote %s (%d KB)", dest.name, dest.stat().st_size // 1024)


def _ensure_raw(refetch: bool) -> None:
    if refetch or not HMM_PGAP_RAW.exists():
        _download(HMM_PGAP_URL, HMM_PGAP_RAW)


def _load_tigr_roles_json() -> dict | None:
    if TIGR_ROLES_JSON.exists():
        with open(TIGR_ROLES_JSON, encoding="utf-8") as fh:
            return json.load(fh)
    return None


def build_tigr_roles(force: bool = False, refetch_raw: bool = False) -> dict:
    """Build (and cache) ``tigr_roles.json`` from the frozen TIGRFAMs 15.0 archive.

    ``{"release", "roles": {role_id: {mainrole, sub1role}}, "family_role":
    {TIGR_acc: role_id}}``. Unnamed roles (``719``) are excluded from both maps.
    On a download failure the committed file is reused with a warning (TCDB
    outage precedent); only a missing file is fatal.
    """
    existing = _load_tigr_roles_json()
    if existing is not None and not force and not refetch_raw:
        logger.info("TIGR roles cache exists: %s (use --force to rebuild)", TIGR_ROLES_JSON)
        return existing

    try:
        for url, dest in ((TIGR_ROLE_NAMES_URL, TIGR_ROLE_NAMES_RAW),
                          (TIGRFAMS_ROLE_LINK_URL, TIGRFAMS_ROLE_LINK_RAW)):
            if refetch_raw or not dest.exists():
                _download(url, dest)
    except Exception as exc:  # network / FTP outage
        if existing is not None:
            logger.warning("TIGRFAMs archive download failed (%s); reusing committed %s",
                           exc, TIGR_ROLES_JSON)
            return existing
        raise

    with open(TIGR_ROLE_NAMES_RAW, encoding="utf-8") as fh:
        roles = parse_tigr_role_names(fh)
    with open(TIGRFAMS_ROLE_LINK_RAW, encoding="utf-8") as fh:
        family_role = parse_tigr_role_link(fh, roles)

    out = {"release": TIGR_ROLES_RELEASE, "roles": roles, "family_role": family_role}
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(TIGR_ROLES_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    logger.info("Wrote %s: %d named roles, %d family→role links",
                TIGR_ROLES_JSON, len(roles), len(family_role))
    return out


def build(force: bool = False, refetch_raw: bool = False) -> dict[str, dict]:
    """Build (and cache) the reference dict. Returns it."""
    if REFERENCE_JSON.exists() and not force and not refetch_raw:
        logger.info("NCBIfam reference cache exists: %s (use --force to rebuild)", REFERENCE_JSON)
        build_tigr_roles(force=False, refetch_raw=False)
        with open(REFERENCE_JSON, encoding="utf-8") as fh:
            return json.load(fh)

    _ensure_raw(refetch_raw)
    with open(HMM_PGAP_RAW, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        ref = parse_hmm_pgap_rows(reader)

    # QC: report family_type distribution + xref coverage; warn on collapse —
    # a zero count almost certainly means the TSV column layout changed, not
    # that NCBIfam genuinely stopped annotating.
    type_counts: dict[str, int] = {}
    with_gene_symbol = with_desc = with_synonyms = with_pmids = with_ec = with_go = 0
    hypoth_count = 0
    for meta in ref.values():
        ft = meta.get("family_type", "")
        type_counts[ft] = type_counts.get(ft, 0) + 1
        if ft in HYPOTH_FAMILY_TYPES:
            hypoth_count += 1
        if meta.get("gene_symbol"):
            with_gene_symbol += 1
        if meta.get("description"):
            with_desc += 1
        if meta.get("gene_synonyms"):
            with_synonyms += 1
        if meta.get("pmids"):
            with_pmids += 1
        if meta.get("ec_numbers"):
            with_ec += 1
        if meta.get("go_terms"):
            with_go += 1
    if not ref:
        logger.warning("Parsed 0 entries — did hmm_PGAP.tsv column layout change?")
    if not with_gene_symbol:
        logger.warning("No entry carries gene_symbol — did the TSV column name change?")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(REFERENCE_JSON, "w", encoding="utf-8") as fh:
        json.dump(ref, fh, indent=1, sort_keys=True)

    logger.info(
        "Wrote %s: %d entries (%d hypoth, %d gene_symbol, %d description, "
        "%d gene_synonyms, %d pmids, %d ec_numbers, %d go_terms); "
        "family_types=%s",
        REFERENCE_JSON, len(ref), hypoth_count, with_gene_symbol, with_desc,
        with_synonyms, with_pmids, with_ec, with_go,
        dict(sorted(type_counts.items())),
    )
    build_tigr_roles(force=force, refetch_raw=refetch_raw)
    return ref


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build the NCBIfam reference cache")
    ap.add_argument("--force", action="store_true", help="rebuild the JSON from cached raw TSV")
    ap.add_argument("--refetch-raw", action="store_true", help="re-download the raw TSV first")
    args = ap.parse_args()
    build(force=args.force, refetch_raw=args.refetch_raw)


if __name__ == "__main__":
    main()
