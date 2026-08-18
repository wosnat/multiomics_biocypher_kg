"""Build the committed MEROPS reference cache.

Turns MEROPS's ``family.txt`` + ``clan.txt`` (release 12.5) into
``cache/data/merops/merops_reference.json``::

    {"families": {code: {"name", "clan", "family_type"}},
     "clans":    {code: {"description", "family_type"}}}

consumed by ``merops_adapter`` at KG-build time for MeropsFamily node names /
clan descriptions / peptidase-vs-inhibitor typing. Committed so the Docker
build runs offline (InterPro/NCBIfam precedent).

Input resolution order:

1. ``$MEROPS_DATA_DIR/DB/{family,clan}.txt`` (default ``~/tools/MEROPS``) —
   already present on any machine that ran the ``/merops-diamond`` Phase-1
   skill (its self-install downloads them).
2. Otherwise self-download (~300 KB total) from the EBI FTP
   ``current_release`` into ``cache/data/merops/raw/`` (gitignored).

Pure parsing lives in ``multiomics_kg/utils/merops_diamond.py``
(``parse_family_txt`` / ``parse_clan_txt``); this module only handles
locating/downloading the raw files and I/O. Effectively run-once: MEROPS has
been in maintenance mode since release 12.5 (Sept 2023).

Usage (module):
    uv run python -m multiomics_kg.download.build_merops_reference [--force] [--refetch-raw]

- default: build the JSON from local raw files (or download if absent).
- ``--force``: rebuild the JSON even if it already exists.
- ``--refetch-raw``: re-download the raw txt files first (only on a MEROPS
  release).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
from pathlib import Path

import requests

from multiomics_kg.utils.merops_diamond import (
    aggregate_cleavages,
    cleavage_properties,
    parse_clan_txt,
    parse_family_txt,
    parse_interpro_txt_stream,
    type_example_names,
)

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

MEROPS_RELEASE_URL = "https://ftp.ebi.ac.uk/pub/databases/merops/current_release"

CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "merops"
RAW_DIR = CACHE_DIR / "raw"
REFERENCE_JSON = CACHE_DIR / "merops_reference.json"

_SUBFAMILY_RE = re.compile(r"^[A-Z]\d+[A-Z]$")

# name → (URL, filename under $MEROPS_DATA_DIR/DB/ from the Phase-1 install)
RAW_FILES = {
    "family.txt": (f"{MEROPS_RELEASE_URL}/database_files/family.txt", "family.txt"),
    "clan.txt": (f"{MEROPS_RELEASE_URL}/database_files/clan.txt", "clan.txt"),
    # the Phase-1 runner keeps the pristine download as merops_scan.raw.lib
    "merops_scan.lib": (f"{MEROPS_RELEASE_URL}/merops_scan.lib", "merops_scan.raw.lib"),
    # Pfam bridge + cleavage specificity (2026-08-18 follow-up). Neither ships
    # with the Phase-1 tool install — they download straight into raw/.
    "interpro.txt": (f"{MEROPS_RELEASE_URL}/interpro.txt", "interpro.txt"),
    "Substrate_search.txt": (f"{MEROPS_RELEASE_URL}/Substrate_search.txt", "Substrate_search.txt"),
}


def _merops_data_dir() -> Path:
    return Path(os.environ.get("MEROPS_DATA_DIR", "~/tools/MEROPS")).expanduser()


def _download(url: str, dest: Path) -> None:
    logger.info("Downloading %s → %s", url, dest)
    resp = requests.get(url, stream=True, timeout=180)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    logger.info("Wrote %s (%d KB)", dest.name, max(1, dest.stat().st_size // 1024))


def _locate_raw(name: str, refetch: bool) -> Path:
    """Prefer the Phase-1 tool install; else the gitignored raw cache (download)."""
    url, tool_name = RAW_FILES[name]
    tool_copy = _merops_data_dir() / "DB" / tool_name
    if not refetch and tool_copy.exists():
        return tool_copy
    raw_copy = RAW_DIR / name
    if refetch or not raw_copy.exists():
        _download(url, raw_copy)
    return raw_copy


def build(force: bool = False, refetch_raw: bool = False) -> dict[str, dict]:
    """Build (and cache) the reference dict. Returns it."""
    if REFERENCE_JSON.exists() and not force and not refetch_raw:
        logger.info("MEROPS reference cache exists: %s (use --force to rebuild)", REFERENCE_JSON)
        with open(REFERENCE_JSON, encoding="utf-8") as fh:
            return json.load(fh)

    # latin-1: MEROPS files carry 0x96 en-dashes in names (scan-lib precedent)
    family_path = _locate_raw("family.txt", refetch_raw)
    clan_path = _locate_raw("clan.txt", refetch_raw)
    scan_lib_path = _locate_raw("merops_scan.lib", refetch_raw)
    families = parse_family_txt(family_path.read_text(encoding="latin-1"))
    clans = parse_clan_txt(clan_path.read_text(encoding="latin-1"))
    # family.txt names only ~27% of families; fall back to each family's /
    # subfamily's type-example name derived from the scan library (MEROPS's
    # own naming convention — the .001 holotype, e.g. M01 → "aminopeptidase N").
    example_names = type_example_names(scan_lib_path.read_text(encoding="latin-1"))

    interpro_path = _locate_raw("interpro.txt", refetch_raw)
    with open(interpro_path, encoding="latin-1", newline="") as fh:
        pfam_bridge = parse_interpro_txt_stream(fh)

    substrate_path = _locate_raw("Substrate_search.txt", refetch_raw)
    with open(substrate_path, encoding="latin-1") as fh:
        cleavage_agg = aggregate_cleavages(fh)
    cleavage = {
        fam: props
        for fam, agg in sorted(cleavage_agg.items())
        if (props := cleavage_properties(agg))
    }

    ref = {
        "families": {
            code: {
                "name": meta.get("name") or example_names.get(code),
                "clan": meta.get("clan"),
                "family_type": meta.get("entry_type"),
            }
            for code, meta in sorted(families.items())
        },
        # subfamily codes only exist in the scan library, not family.txt;
        # clan + family_type are derivable from the parent family at read time.
        "subfamily_names": {
            code: name for code, name in sorted(example_names.items())
            if _SUBFAMILY_RE.match(code)
        },
        "clans": clans,
        "pfam_bridge": pfam_bridge,
        "cleavage": cleavage,
    }

    # QC — zero counts almost certainly mean the txt column layout changed.
    named = sum(1 for m in ref["families"].values() if m["name"])
    with_clan = sum(1 for m in ref["families"].values() if m["clan"])
    inhibitors = sum(1 for m in ref["families"].values() if m["family_type"] == "inhibitor")
    if not ref["families"]:
        logger.warning("Parsed 0 families — did family.txt column layout change?")
    if not ref["clans"]:
        logger.warning("Parsed 0 clans — did clan.txt column layout change?")
    if not inhibitors:
        logger.warning("No inhibitor families parsed — did the type column move?")
    if named < len(ref["families"]) // 2:
        logger.warning(
            "Only %d/%d families named — did the scan-lib fallback break?",
            named, len(ref["families"]),
        )
    if not pfam_bridge:
        logger.warning("Parsed 0 pfam_bridge pairs — did interpro.txt column layout change?")
    if not cleavage:
        logger.warning("Parsed 0 cleavage profiles — did Substrate_search.txt column layout change?")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(REFERENCE_JSON, "w", encoding="utf-8") as fh:
        json.dump(ref, fh, indent=1, sort_keys=True)

    logger.info(
        "Wrote %s: %d families (%d named, %d clan-assigned, %d inhibitor), "
        "%d subfamily names, %d clans, %d pfam-bridge families, %d cleavage families",
        REFERENCE_JSON, len(ref["families"]), named, with_clan, inhibitors,
        len(ref["subfamily_names"]), len(ref["clans"]),
        len(pfam_bridge), len(cleavage),
    )
    return ref


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build the MEROPS reference cache")
    ap.add_argument("--force", action="store_true", help="rebuild the JSON from local raw files")
    ap.add_argument("--refetch-raw", action="store_true", help="re-download the raw txt files first")
    args = ap.parse_args()
    build(force=args.force, refetch_raw=args.refetch_raw)


if __name__ == "__main__":
    main()
