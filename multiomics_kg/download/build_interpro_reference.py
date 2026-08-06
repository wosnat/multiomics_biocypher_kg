"""prepare_data step 9 — build the committed InterPro reference cache.

Downloads the two InterPro ``current_release`` text files and writes
``cache/data/interpro/interpro_reference.json``
(``{IPRxxxxxx: {"name", "type", "parent", "level"}}``), consumed by
``interpro_adapter`` at KG-build time for InterProEntry node names/types and the
``Interpro_entry_is_a_interpro_entry`` hierarchy. Committed so the Docker build
runs offline (Pfam/KEGG precedent).

Pure parsing lives in ``multiomics_kg/utils/interpro_reference.py``; this module
only handles download, caching, and I/O.

Usage (module):
    uv run python -m multiomics_kg.download.build_interpro_reference [--force] [--refetch-raw]

- default: build the JSON from cached raw files (or download them if absent).
- ``--force``: rebuild the JSON even if it already exists (from cached raw).
- ``--refetch-raw``: re-download the raw FTP files first (only on an InterPro
  release).

See ``docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md``.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import requests

from multiomics_kg.utils.interpro_reference import KNOWN_INTERPRO_TYPES, build_reference

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

INTERPRO_BASE = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release"
ENTRY_LIST_URL = f"{INTERPRO_BASE}/entry.list"
PARENT_CHILD_URL = f"{INTERPRO_BASE}/ParentChildTreeFile.txt"

CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "interpro"
RAW_DIR = CACHE_DIR / "raw"
ENTRY_LIST_RAW = RAW_DIR / "entry.list"
PARENT_CHILD_RAW = RAW_DIR / "ParentChildTreeFile.txt"
REFERENCE_JSON = CACHE_DIR / "interpro_reference.json"


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
    if refetch or not ENTRY_LIST_RAW.exists():
        _download(ENTRY_LIST_URL, ENTRY_LIST_RAW)
    if refetch or not PARENT_CHILD_RAW.exists():
        _download(PARENT_CHILD_URL, PARENT_CHILD_RAW)


def build(force: bool = False, refetch_raw: bool = False) -> dict[str, dict]:
    """Build (and cache) the reference dict. Returns it."""
    if REFERENCE_JSON.exists() and not force and not refetch_raw:
        logger.info("InterPro reference cache exists: %s (use --force to rebuild)", REFERENCE_JSON)
        with open(REFERENCE_JSON, encoding="utf-8") as fh:
            return json.load(fh)

    _ensure_raw(refetch_raw)
    entry_text = ENTRY_LIST_RAW.read_text(encoding="utf-8")
    tree_text = PARENT_CHILD_RAW.read_text(encoding="utf-8")

    ref = build_reference(entry_text, tree_text)

    # QC: report type distribution + hierarchy coverage; warn on unexpected types.
    type_counts: dict[str, int] = {}
    in_tree = 0
    for meta in ref.values():
        type_counts[meta["type"]] = type_counts.get(meta["type"], 0) + 1
        if meta["parent"] is not None:
            in_tree += 1
    unexpected = {t for t in type_counts if t and t not in KNOWN_INTERPRO_TYPES}
    if unexpected:
        logger.warning("Unexpected InterPro entry types (new release?): %s", sorted(unexpected))

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(REFERENCE_JSON, "w", encoding="utf-8") as fh:
        json.dump(ref, fh, indent=1, sort_keys=True)

    logger.info(
        "Wrote %s: %d entries (%d with a parent); types=%s",
        REFERENCE_JSON, len(ref), in_tree, dict(sorted(type_counts.items())),
    )
    return ref


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build the InterPro reference cache (prepare_data step 9)")
    ap.add_argument("--force", action="store_true", help="rebuild the JSON from cached raw files")
    ap.add_argument("--refetch-raw", action="store_true", help="re-download the raw FTP files first")
    args = ap.parse_args()
    build(force=args.force, refetch_raw=args.refetch_raw)


if __name__ == "__main__":
    main()
