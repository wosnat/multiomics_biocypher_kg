"""prepare_data step 9 — build the committed InterPro reference cache.

Downloads the InterPro ``current_release`` reference files and writes
``cache/data/interpro/interpro_reference.json``
(``{IPRxxxxxx: {"name", "type", "parent", "level"[, "go_terms"][, "pathways"]}}``),
consumed by ``interpro_adapter`` at KG-build time for InterProEntry node
names/types and the ``Interpro_entry_is_a_interpro_entry`` hierarchy. Committed
so the Docker build runs offline (Pfam/KEGG precedent).

Four source files:

- ``entry.list`` (names + types) and ``ParentChildTreeFile.txt`` (is-a hierarchy)
- ``interpro2go`` (~3 MB) — entry→GO
- ``interpro.xml.gz`` (~42 MB) — entry→pathway; the only source of these xrefs.
  Streamed gzip-wise, never fully decompressed.

**No KEGG pathways exist in InterPro** (verified against the full release,
2026-08-06: 0 ``db="KEGG"`` xrefs; MetaCyc 79,683; Reactome 506,724). The KG's
KEGG pathway layer remains eggNOG-KO-derived. Reactome is excluded by default —
its xrefs are species-expanded and irrelevant to marine bacteria.

GO/pathway xrefs are entry-level in InterProScan's own output too (they come from
``signature.entry``), so building them from the reference release is equivalent
to re-running InterProScan with ``--goterms --pathways`` over all strains, minus
~27 wallclock hours of Docker.

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
import gzip
import json
import logging
from pathlib import Path

import requests

from multiomics_kg.utils.interpro_reference import (
    DEFAULT_PATHWAY_DBS,
    KNOWN_INTERPRO_TYPES,
    build_reference,
    parse_interpro2go,
    parse_pathway_xrefs,
)

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

INTERPRO_BASE = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release"
ENTRY_LIST_URL = f"{INTERPRO_BASE}/entry.list"
PARENT_CHILD_URL = f"{INTERPRO_BASE}/ParentChildTreeFile.txt"
INTERPRO2GO_URL = f"{INTERPRO_BASE}/interpro2go"
INTERPRO_XML_URL = f"{INTERPRO_BASE}/interpro.xml.gz"

CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "interpro"
RAW_DIR = CACHE_DIR / "raw"
ENTRY_LIST_RAW = RAW_DIR / "entry.list"
PARENT_CHILD_RAW = RAW_DIR / "ParentChildTreeFile.txt"
INTERPRO2GO_RAW = RAW_DIR / "interpro2go"
INTERPRO_XML_RAW = RAW_DIR / "interpro.xml.gz"
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
    for url, dest in (
        (ENTRY_LIST_URL, ENTRY_LIST_RAW),
        (PARENT_CHILD_URL, PARENT_CHILD_RAW),
        (INTERPRO2GO_URL, INTERPRO2GO_RAW),
        (INTERPRO_XML_URL, INTERPRO_XML_RAW),
    ):
        if refetch or not dest.exists():
            _download(url, dest)


def build(
    force: bool = False,
    refetch_raw: bool = False,
    pathway_dbs: tuple[str, ...] = DEFAULT_PATHWAY_DBS,
) -> dict[str, dict]:
    """Build (and cache) the reference dict. Returns it."""
    if REFERENCE_JSON.exists() and not force and not refetch_raw:
        logger.info("InterPro reference cache exists: %s (use --force to rebuild)", REFERENCE_JSON)
        with open(REFERENCE_JSON, encoding="utf-8") as fh:
            return json.load(fh)

    _ensure_raw(refetch_raw)
    entry_text = ENTRY_LIST_RAW.read_text(encoding="utf-8")
    tree_text = PARENT_CHILD_RAW.read_text(encoding="utf-8")

    go_map = parse_interpro2go(INTERPRO2GO_RAW.read_text(encoding="utf-8"))
    logger.info(
        "interpro2go: %d entries → %d GO mappings",
        len(go_map), sum(len(v) for v in go_map.values()),
    )

    # Stream the 42 MB gzip rather than decompressing it (uncompressed ~1 GB).
    with gzip.open(INTERPRO_XML_RAW, "rt", encoding="utf-8", errors="replace") as fh:
        pathway_map = parse_pathway_xrefs(fh, include_dbs=pathway_dbs)
    logger.info(
        "interpro.xml (%s): %d entries → %d pathway xrefs",
        ",".join(pathway_dbs), len(pathway_map), sum(len(v) for v in pathway_map.values()),
    )

    ref = build_reference(entry_text, tree_text, go_map=go_map, pathway_map=pathway_map)

    # QC: report type distribution + hierarchy/xref coverage; warn on unexpected types.
    type_counts: dict[str, int] = {}
    in_tree = with_go = with_pw = 0
    for meta in ref.values():
        type_counts[meta["type"]] = type_counts.get(meta["type"], 0) + 1
        if meta["parent"] is not None:
            in_tree += 1
        if meta.get("go_terms"):
            with_go += 1
        if meta.get("pathways"):
            with_pw += 1
    unexpected = {t for t in type_counts if t and t not in KNOWN_INTERPRO_TYPES}
    if unexpected:
        logger.warning("Unexpected InterPro entry types (new release?): %s", sorted(unexpected))
    # A collapse here means the release changed shape (renamed file / new schema),
    # not that InterPro genuinely stopped annotating — fail loudly rather than
    # silently committing a reference cache with no xrefs.
    if not with_go:
        logger.warning("No entry carries go_terms — did interpro2go change format?")
    if pathway_dbs and not with_pw:
        logger.warning("No entry carries pathways — did interpro.xml change format?")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(REFERENCE_JSON, "w", encoding="utf-8") as fh:
        json.dump(ref, fh, indent=1, sort_keys=True)

    logger.info(
        "Wrote %s: %d entries (%d with a parent, %d with GO, %d with pathways); types=%s",
        REFERENCE_JSON, len(ref), in_tree, with_go, with_pw, dict(sorted(type_counts.items())),
    )
    return ref


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build the InterPro reference cache (prepare_data step 9)")
    ap.add_argument("--force", action="store_true", help="rebuild the JSON from cached raw files")
    ap.add_argument("--refetch-raw", action="store_true", help="re-download the raw FTP files first")
    ap.add_argument(
        "--pathway-dbs", default=",".join(DEFAULT_PATHWAY_DBS),
        help=("comma-separated InterPro pathway db tokens to keep (default: %(default)s). "
              "REACTOME is available but species-expanded (~507K xrefs) and noisy for "
              "bacteria. KEGG does not exist in InterPro. Empty string = no pathways."),
    )
    args = ap.parse_args()
    pathway_dbs = tuple(d.strip().upper() for d in args.pathway_dbs.split(",") if d.strip())
    build(force=args.force, refetch_raw=args.refetch_raw, pathway_dbs=pathway_dbs)


if __name__ == "__main__":
    main()
