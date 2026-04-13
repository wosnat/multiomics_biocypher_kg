"""
KEGG BRITE functional hierarchy utilities.

Downloads 12 BRITE trees via KEGG REST JSON endpoint and caches them at
<cache_root>/kegg/brite_<tree_code>.json (same kegg/ subdirectory as
kegg_data.json).

Reuses _fetch_json from kegg_utils to avoid duplicating HTTP logic.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

from multiomics_kg.utils.kegg_utils import _fetch_json

logger = logging.getLogger(__name__)

# 12 configured BRITE trees: tree_code → canonical tree name used as node property
BRITE_TREES: dict[str, str] = {
    "ko01000": "enzymes",
    "ko02000": "transporters",
    "ko01002": "peptidases",
    "ko03000": "transcription_factors",
    "ko02044": "secretion",
    "ko02022": "two_component",
    "ko02048": "defense",
    "ko03110": "chaperones",
    "ko03011": "ribosome",
    "ko03012": "translation_factors",
    "ko03016": "trna_biogenesis",
    "ko03032": "dna_replication",
}

_LEVEL_KINDS: dict[int, str] = {
    0: "brite_class",
    1: "brite_subclass",
    2: "brite_family",
    3: "brite_subfamily",
}

_KEGG_BRITE_URL = "https://rest.kegg.jp/get/br:{tree_code}/json"
_RATE_LIMIT_SLEEP = 0.2  # seconds between sequential tree fetches


def compute_level_kind(depth: int) -> str:
    """Map nesting depth (0–3) to BriteCategory level_kind label.

    depth 0 → 'brite_class'    (A-level, broadest, no parent edge)
    depth 1 → 'brite_subclass' (B-level)
    depth 2 → 'brite_family'   (C-level)
    depth 3 → 'brite_subfamily' (D-level, non-KO only)

    Raises ValueError for depth >= 4 (not expected for the 12 configured trees).
    """
    if depth not in _LEVEL_KINDS:
        raise ValueError(
            f"Unexpected BRITE nesting depth {depth}; maximum supported is 3. "
            f"Check whether the tree has unexpected extra levels."
        )
    return _LEVEL_KINDS[depth]


def download_brite_tree(
    tree_code: str,
    cache_root: Path,
    cache: bool = True,
) -> dict:
    """Fetch and cache one KEGG BRITE tree in JSON format.

    Args:
        tree_code: KEGG BRITE tree ID, e.g. ``"ko02000"``
        cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
        cache: if False, re-download even if cache file exists

    Returns:
        Parsed JSON dict from KEGG REST (keys: ``"name"``, ``"children"``).

    Raises:
        requests.HTTPError: on non-200 HTTP response
    """
    cache_dir = Path(cache_root) / "kegg"
    cache_file = cache_dir / f"brite_{tree_code}.json"

    if cache_file.exists() and cache:
        logger.info(f"Loading BRITE tree {tree_code} from cache: {cache_file}")
        with open(cache_file, encoding="utf-8") as fh:
            return json.load(fh)

    url = _KEGG_BRITE_URL.format(tree_code=tree_code)
    data = _fetch_json(url)

    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    logger.info(f"BRITE tree {tree_code} cached to {cache_file}")
    return data


def load_brite_trees(
    cache_root: Path,
    trees: list[str],
    cache: bool = True,
) -> dict[str, dict]:
    """Load all requested BRITE trees, sleeping between sequential fetches.

    Args:
        cache_root: project-level cache directory
        trees: list of tree codes, e.g. ``["ko02000", "ko01002"]``
        cache: passed through to ``download_brite_tree``

    Returns:
        ``{tree_code: parsed_json}`` for every requested tree.
    """
    result: dict[str, dict] = {}
    for i, tree_code in enumerate(trees):
        if i > 0:
            time.sleep(_RATE_LIMIT_SLEEP)
        result[tree_code] = download_brite_tree(tree_code, cache_root, cache=cache)
    return result
