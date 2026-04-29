"""TCDB hierarchy accessor API.

Reads cache/data/tcdb/tcdb_hierarchy.json (built by sub-step 7).
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT_PATH = Path("cache/data/tcdb/tcdb_hierarchy.json")
_CACHE: dict[str, dict] | None = None


def load_tcdb() -> dict[str, dict]:
    """Load the TCDB hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_tcdb(tc_id: str) -> bool:
    """Return True if `tc_id` exists as a key in the hierarchy."""
    if not tc_id:
        return False
    return tc_id in load_tcdb()


def tcdb_ancestors(tc_id: str) -> list[str]:
    """Return root-to-parent ancestor chain for `tc_id`. Empty for unknown / root."""
    h = load_tcdb()
    if tc_id not in h:
        return []
    chain: list[str] = []
    parent = h[tc_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))
