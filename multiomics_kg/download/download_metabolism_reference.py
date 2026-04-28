"""Step 0 sub-step 6 — Download MNX and TCDB reference data.

Seven downloads (4 MNX TSVs + 3 TCDB TSVs). MNX dominates total size (~1.5 GB
unzipped); TCDB tables are <1 MB each. Files are cached under
cache/data/{mnx,tcdb}/ and skipped on re-run unless --force.

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
SOURCES: dict[str, tuple[str, str]] = {
    "mnx_chem_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv",       "mnx/chem_prop.tsv"),
    "mnx_chem_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",       "mnx/chem_xref.tsv"),
    "mnx_reac_prop":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv",       "mnx/reac_prop.tsv"),
    "mnx_reac_xref":      ("https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",       "mnx/reac_xref.tsv"),
    "tcdb_families":      ("https://www.tcdb.org/cgi-bin/projectv/public/families.py",           "tcdb/families.tsv"),
    "tcdb_substrates":    ("https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py",           "tcdb/substrates.tsv"),
    "tcdb_superfamilies": ("https://www.tcdb.org/cgi-bin/substrates/listSuperfamilies.py",       "tcdb/superfamilies.tsv"),
}

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


def download_all(cache_root: Path = DEFAULT_CACHE_ROOT, force: bool = False) -> None:
    """Download every source in SOURCES into cache_root."""
    log.info(f"download_metabolism_reference: cache_root={cache_root} force={force}")
    n_downloaded = 0
    for key, (url, rel) in SOURCES.items():
        if download_one(url, cache_root / rel, force=force):
            n_downloaded += 1
    log.info(f"  done — {n_downloaded}/{len(SOURCES)} downloaded, "
             f"{len(SOURCES) - n_downloaded} cached.")


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    download_all(force=force)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if cache files exist.")
    args = parser.parse_args()
    main(force=args.force)
