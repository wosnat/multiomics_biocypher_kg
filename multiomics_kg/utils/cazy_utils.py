"""CAZy hierarchy accessor API.

Reads cache/data/cazy/cazy_hierarchy.json (built by multiomics_kg.download.build_metabolite_resolver
from observed eggNOG `CAZy` columns).
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT_PATH = Path("cache/data/cazy/cazy_hierarchy.json")
_CACHE: dict[str, dict] | None = None


def load_cazy() -> dict[str, dict]:
    """Load the CAZy hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_cazy(cazy_id: str) -> bool:
    if not cazy_id:
        return False
    return cazy_id in load_cazy()


def cazy_ancestors(cazy_id: str) -> list[str]:
    h = load_cazy()
    if cazy_id not in h:
        return []
    chain: list[str] = []
    parent = h[cazy_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))
