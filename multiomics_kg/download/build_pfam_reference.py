"""Build the committed Pfam reference cache (prepare_data step 9).

Writes ``cache/data/pfam/pfam_reference.json``
(``{release, by_accession, by_shortname, clans, dead}``), consumed by the step-2
gene-annotation merge (``enrich_pfam_fields``: short name -> accession, killed
families forwarded to their successor) and by ``MultiPfamAnnotationAdapter`` at
KG-build time for Pfam / PfamClan node names. Committed so the Docker build runs
offline; the KG build never downloads Pfam.

Raw files (EBI FTP ``Pfam/current_release`` + ``releases/Pfam33.1``), kept in the
gitignored ``cache/data/pfam/raw/``:

- ``Pfam-A.clans.tsv.gz`` — every family: accession, clan, short name, description.
- ``Pfam-A.dead.gz`` — killed families: short name, the family each was merged
  into (``FW``, often none) and why.
- ``Pfam.version.gz`` — the release number, stored as ``release``.
- Pfam 33.1 ``Pfam-A.clans.tsv.gz`` — eggNOG-mapper 2.1.13 annotates with Pfam 33.1
  short names; renamed ones are added to ``by_shortname``.

Usage (module):
    uv run python -m multiomics_kg.download.build_pfam_reference [--force] [--refetch-raw]

- default: build the JSON from the cached raw files (downloading any that are
  missing) unless it already exists.
- ``--force``: rebuild the JSON even if it already exists (from cached raw).
- ``--refetch-raw``: re-download the raw files first (only on a Pfam release;
  keep it in step with the InterProScan data release).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from multiomics_kg.utils.pfam_utils import build_pfam_reference, reference_path

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CACHE_ROOT = PROJECT_ROOT / "cache" / "data"


def build(force: bool = False, refetch_raw: bool = False, cache_root: Path = CACHE_ROOT) -> None:
    out = reference_path(cache_root)
    if out.exists() and not (force or refetch_raw):
        logger.info("%s exists; skipping (use --force to rebuild)", out)
        return
    build_pfam_reference(cache_root, refetch_raw=refetch_raw)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build the Pfam reference cache")
    ap.add_argument("--force", action="store_true", help="rebuild the JSON from cached raw files")
    ap.add_argument("--refetch-raw", action="store_true", help="re-download the raw files first")
    args = ap.parse_args()
    build(force=args.force, refetch_raw=args.refetch_raw)


if __name__ == "__main__":
    main()
