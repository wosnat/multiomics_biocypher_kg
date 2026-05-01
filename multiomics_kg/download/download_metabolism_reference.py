"""Download MNX and/or TCDB reference data.

Default behaviour downloads both groups (4 MNX TSVs + 3 TCDB TSVs). Use
``--sources mnx`` or ``--sources tcdb`` to download a single group. MNX dominates
total size (~1.5 GB unzipped); TCDB tables are <1 MB each. Files are cached
under cache/data/{mnx,tcdb}/ and skipped on re-run unless --force.

Used by:
- prepare_data.sh step 0 sub-step 6 (TCDB only — MNX moved to scripts/refresh_mnx.sh)
- scripts/refresh_mnx.sh (MNX only)

CAZy family hierarchy is NOT downloaded here. It is bootstrapped in Phase 1.1B
from observed eggNOG `CAZy` columns plus mechanical ID parsing — the format
itself encodes the hierarchy (`GH13_1` → family `GH13` → class `GH`), and we
only care about families our genes actually reference.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import requests

log = logging.getLogger(__name__)

# (URL, cache-relative path)
MNX_SOURCES: dict[str, tuple[str, str]] = {
    "mnx_chem_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",       "mnx/chem_prop.tsv"),
    "mnx_chem_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",       "mnx/chem_xref.tsv"),
    "mnx_reac_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",       "mnx/reac_prop.tsv"),
    "mnx_reac_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",       "mnx/reac_xref.tsv"),
}

TCDB_SOURCES: dict[str, tuple[str, str]] = {
    "tcdb_families":      ("https://www.tcdb.org/cgi-bin/projectv/public/families.py",           "tcdb/families.tsv"),
    "tcdb_substrates":    ("https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py",           "tcdb/substrates.tsv"),
    "tcdb_superfamilies": ("https://www.tcdb.org/cgi-bin/substrates/listSuperfamilies.py",       "tcdb/superfamilies.tsv"),
}

SOURCES_BY_GROUP: dict[str, dict[str, tuple[str, str]]] = {
    "mnx": MNX_SOURCES,
    "tcdb": TCDB_SOURCES,
}

# Backward-compat flat view of all sources (union of both groups)
SOURCES: dict[str, tuple[str, str]] = {**MNX_SOURCES, **TCDB_SOURCES}

DEFAULT_CACHE_ROOT = Path("cache/data")


def download_one(url: str, dest: Path, force: bool = False) -> bool:
    """Download URL → dest. Returns True if downloaded, False if skipped."""
    if dest.exists() and not force:
        log.info(f"  skip (cached): {dest}")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"  GET {url} → {dest}")
    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=64 * 1024):
                if chunk:
                    f.write(chunk)
    return True


def download_all(
    cache_root: Path = DEFAULT_CACHE_ROOT,
    force: bool = False,
    sources: list[str] | None = None,
) -> None:
    """Download requested source groups into cache_root.

    Args:
        cache_root: where to write files (default cache/data)
        force: re-download even if cached
        sources: list of source-group keys (e.g. ['mnx', 'tcdb']); None = all groups
    """
    selected_groups = sources or list(SOURCES_BY_GROUP.keys())
    invalid = [s for s in selected_groups if s not in SOURCES_BY_GROUP]
    if invalid:
        raise ValueError(
            f"Unknown source group(s): {invalid}. "
            f"Valid: {sorted(SOURCES_BY_GROUP.keys())}"
        )

    log.info(
        f"download_metabolism_reference: cache_root={cache_root} "
        f"force={force} groups={selected_groups}"
    )
    n_downloaded = 0
    n_total = 0
    for group in selected_groups:
        for key, (url, rel) in SOURCES_BY_GROUP[group].items():
            n_total += 1
            if download_one(url, cache_root / rel, force=force):
                n_downloaded += 1
    log.info(f"  done — {n_downloaded}/{n_total} downloaded, "
             f"{n_total - n_downloaded} cached.")


def main(force: bool = False, sources: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    download_all(force=force, sources=sources)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if cache files exist.")
    parser.add_argument(
        "--sources", nargs="+", choices=sorted(SOURCES_BY_GROUP.keys()),
        default=None,
        help=f"Source groups to download (default: all = {sorted(SOURCES_BY_GROUP.keys())})",
    )
    args = parser.parse_args()
    main(force=args.force, sources=args.sources)
