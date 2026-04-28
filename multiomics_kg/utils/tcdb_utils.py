"""TCDB hierarchy accessor API.

Phase 1.1A skeleton — every function raises NotImplementedError until 1.1B.

API:
    load_tcdb()                 -> dict[str, dict]      # cached at module level
    tcdb_ancestors(tc_id: str)  -> list[str]            # ["1", "1.A", "1.A.1", "1.A.1.1"] for "1.A.1.1.1"
    is_valid_tcdb(tc_id: str)   -> bool
"""
from __future__ import annotations


def load_tcdb() -> dict[str, dict]:
    raise NotImplementedError("Phase 1.1B")


def tcdb_ancestors(tc_id: str) -> list[str]:
    raise NotImplementedError("Phase 1.1B")


def is_valid_tcdb(tc_id: str) -> bool:
    raise NotImplementedError("Phase 1.1B")
