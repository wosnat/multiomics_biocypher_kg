"""prepare_data step 9 — build the committed InterPro reference cache.

Downloads the InterPro ``current_release`` reference files and writes
``cache/data/interpro/interpro_reference.json``
(``{IPRxxxxxx: {"name", "type", "parent", "level"[, "go_terms"][, "pathways"]
[, "description"]}}``), consumed by ``interpro_adapter`` at KG-build time for
InterProEntry node names/types/descriptions and the
``Interpro_entry_is_a_interpro_entry`` hierarchy. Committed so the Docker
build runs offline (Pfam/KEGG precedent).

Four source files:

- ``entry.list`` (names + types) and ``ParentChildTreeFile.txt`` (is-a hierarchy)
- ``interpro2go`` (~3 MB) — entry→GO
- ``interpro.xml.gz`` (~42 MB) — entry→pathway, and (Task 6) entry→description
  (first ``<abstract>`` paragraph). Streamed gzip-wise in two passes, never
  fully decompressed: db_xrefs (pathways/EC/CAZy) in one, descriptions in the
  other — see ``parse_entry_db_xrefs`` / ``parse_entry_descriptions``.

**Size gate**: descriptions on all ~54K entries can push the committed JSON
past ~25 MB. When that happens, ``build()`` falls back to keeping descriptions
only for InterPro accessions observed in a committed per-strain
``<strain>.interproscan.calls.json`` plus their is-a ancestors (see
``_observed_interpro_ids`` / ``_with_ancestors``) and logs which mode
(``full`` vs ``observed-only``) shipped.

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
    PATHWAY_DB_NAMES,
    build_reference,
    parse_entry_db_xrefs,
    parse_entry_descriptions,
    parse_interpro2go,
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

# Task 6 size gate (spec): the committed reference must stay ~25 MB or smaller.
# Full first-paragraph descriptions on all ~54K entries push it past that; the
# documented fallback keeps descriptions only for entries actually observed in
# a committed per-strain calls.json plus their is-a ancestors (see build()).
MAX_REFERENCE_JSON_BYTES = 25 * 1024 * 1024
INTERPROSCAN_CALLS_GLOB = "*/genomes/*/interproscan/*.interproscan.calls.json"


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


def _observed_interpro_ids(cache_root: Path = PROJECT_ROOT / "cache" / "data") -> set[str]:
    """Union of InterPro accessions observed in any committed per-strain
    ``<strain>.interproscan.calls.json`` (Task 6 size-gate fallback only —
    the reference itself is pruned separately, by ``interpro_adapter``, at
    KG-build time)."""
    ids: set[str] = set()
    for calls_path in cache_root.glob(INTERPROSCAN_CALLS_GLOB):
        with open(calls_path, encoding="utf-8") as fh:
            calls = json.load(fh)
        for protein in calls.values():
            ids.update(protein.get("interpro_entries", {}).keys())
    return ids


def _with_ancestors(ids: set[str], ref: dict[str, dict]) -> set[str]:
    """Extend *ids* with every is-a ancestor reachable via ``ref[acc]["parent"]``."""
    out = set(ids)
    for acc in ids:
        cur = ref.get(acc, {}).get("parent")
        while cur and cur not in out:
            out.add(cur)
            cur = ref.get(cur, {}).get("parent")
    return out


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
    # One pass extracts pathways + EC + CAZy; EC/CAZy are InterPro-XML-only.
    with gzip.open(INTERPRO_XML_RAW, "rt", encoding="utf-8", errors="replace") as fh:
        raw_xrefs = parse_entry_db_xrefs(fh, include_dbs=(*pathway_dbs, "EC", "CAZY"))
    pathway_wanted = {db.upper() for db in pathway_dbs}
    pathway_map: dict[str, list[str]] = {}
    ec_map: dict[str, list[str]] = {}
    cazy_map: dict[str, list[str]] = {}
    for acc, dbs in raw_xrefs.items():
        pws: set[str] = set()
        for db, keys in dbs.items():
            if db in pathway_wanted:
                label = PATHWAY_DB_NAMES.get(db, db)
                pws.update(f"{label}:{k}" for k in keys)
            elif db == "EC":
                ec_map[acc] = keys
            elif db == "CAZY":
                cazy_map[acc] = keys
        if pws:
            pathway_map[acc] = sorted(pws)
    logger.info(
        "interpro.xml: pathways(%s) %d entries/%d xrefs; EC %d entries/%d; CAZy %d entries/%d",
        ",".join(pathway_dbs),
        len(pathway_map), sum(len(v) for v in pathway_map.values()),
        len(ec_map), sum(len(v) for v in ec_map.values()),
        len(cazy_map), sum(len(v) for v in cazy_map.values()),
    )

    # Second streaming pass for entry abstracts (first-paragraph descriptions).
    # Kept separate from the pathway/EC/CAZy pass above because it needs its
    # own multi-line accumulation state per entry — see parse_entry_descriptions().
    with gzip.open(INTERPRO_XML_RAW, "rt", encoding="utf-8", errors="replace") as fh:
        description_map = parse_entry_descriptions(fh)
    logger.info("interpro.xml: descriptions for %d entries", len(description_map))

    ref = build_reference(
        entry_text, tree_text,
        go_map=go_map, pathway_map=pathway_map, ec_map=ec_map, cazy_map=cazy_map,
        description_map=description_map,
    )

    # Task 6 size gate: drop to observed-only descriptions if the full set
    # would push the committed JSON over ~25 MB. "Observed" = any InterPro
    # accession appearing in a committed per-strain calls.json, plus its is-a
    # ancestors (interpro_adapter's hierarchy walk can surface an ancestor
    # node even though genes never attach to it directly).
    description_mode = "full"
    # Must match the indent=1/sort_keys=True kwargs used for the actual write
    # below — a compact (no-indent) estimate here undercounts by several MB
    # and can silently miss the gate.
    approx_size = len(json.dumps(ref, indent=1, sort_keys=True).encode("utf-8"))
    if approx_size > MAX_REFERENCE_JSON_BYTES:
        keep = _with_ancestors(_observed_interpro_ids(), ref)
        dropped = 0
        for acc, meta in ref.items():
            if "description" in meta and acc not in keep:
                del meta["description"]
                dropped += 1
        description_mode = "observed-only"
        logger.warning(
            "Full descriptions (~%d bytes) exceed the %d-byte size gate — "
            "restricted descriptions to %d observed+ancestor entries (dropped %d)",
            approx_size, MAX_REFERENCE_JSON_BYTES, len(keep), dropped,
        )

    # QC: report type distribution + hierarchy/xref coverage; warn on unexpected types.
    type_counts: dict[str, int] = {}
    in_tree = with_go = with_pw = with_ec = with_cazy = with_description = 0
    for meta in ref.values():
        type_counts[meta["type"]] = type_counts.get(meta["type"], 0) + 1
        if meta["parent"] is not None:
            in_tree += 1
        if meta.get("go_terms"):
            with_go += 1
        if meta.get("pathways"):
            with_pw += 1
        if meta.get("ec_numbers"):
            with_ec += 1
        if meta.get("cazy_ids"):
            with_cazy += 1
        if meta.get("description"):
            with_description += 1
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
    # EC/CAZy are new (2026-08-10). A zero count almost certainly means the db
    # token differs from the assumed "EC"/"CAZY" — surface it, don't guess.
    if not with_ec:
        logger.warning("No entry carries ec_numbers — is the interpro.xml db token still 'EC'?")
    if not with_cazy:
        logger.warning("No entry carries cazy_ids — is the interpro.xml db token still 'CAZY'?")
    if not with_description:
        logger.warning("No entry carries description — did the interpro.xml <abstract> shape change?")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(REFERENCE_JSON, "w", encoding="utf-8") as fh:
        json.dump(ref, fh, indent=1, sort_keys=True)
    written_bytes = REFERENCE_JSON.stat().st_size

    logger.info(
        "Wrote %s (%.1f MB): %d entries (%d parent, %d GO, %d pathways, %d EC, %d CAZy, "
        "%d description [%s]); types=%s",
        REFERENCE_JSON, written_bytes / (1024 * 1024), len(ref), in_tree, with_go, with_pw,
        with_ec, with_cazy, with_description, description_mode,
        dict(sorted(type_counts.items())),
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
