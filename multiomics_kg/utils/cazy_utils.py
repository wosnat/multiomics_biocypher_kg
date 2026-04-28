"""CAZy hierarchy accessor API.

Phase 1.1A skeleton — every function raises NotImplementedError until 1.1B.

API:
    load_cazy()                  -> dict[str, dict]     # cached at module level
    cazy_ancestors(cazy_id: str) -> list[str]
    is_valid_cazy(cazy_id: str)  -> bool
"""
from __future__ import annotations


def load_cazy() -> dict[str, dict]:
    raise NotImplementedError("Phase 1.1B")


def cazy_ancestors(cazy_id: str) -> list[str]:
    raise NotImplementedError("Phase 1.1B")


def is_valid_cazy(cazy_id: str) -> bool:
    raise NotImplementedError("Phase 1.1B")
