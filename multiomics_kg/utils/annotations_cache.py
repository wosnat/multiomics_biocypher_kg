"""Per-strain shared cache for ``gene_annotations_merged.json``.

Eight per-strain adapters (functional GO/EC/KEGG/COG/Pfam, plus Cazy, TCDB, and
OrthologGroup) each parse the same JSON file. With ~27 strains × ~6 adapters
that ran for a given build, that's ~160+ redundant parses of multi-MB JSONs
(~30 s in ``json.decoder.raw_decode`` under cProfile).

This module exposes a cached loader. The returned dict is shared across callers
— **do not mutate it**.
"""

from __future__ import annotations

import functools
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


@functools.lru_cache(maxsize=None)
def _load_cached(path_str: str) -> dict:
    path = Path(path_str)
    if not path.exists():
        logger.warning(f"gene_annotations_merged.json not found at {path}")
        return {}
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load_merged_annotations(genome_dir) -> dict:
    """Load (cached) ``gene_annotations_merged.json`` from ``genome_dir``.

    Args:
        genome_dir: directory containing the JSON (str or Path).

    Returns:
        Parsed dict, or ``{}`` if the file is missing.

    The returned dict is shared across all callers for the same path.
    Treat it as read-only.
    """
    return _load_cached(str(Path(genome_dir) / "gene_annotations_merged.json"))
