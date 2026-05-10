# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein post-steps.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.
"""
from __future__ import annotations


def truncate_tcid(tcid: str, n_parts: int) -> str:
    """Truncate a TCID to its first n_parts dot-separated segments.

    A TCID has up to 5 parts (class.subclass.family.subfamily.specificity).
    Returns the input unchanged if it already has <= n_parts segments.
    """
    if not tcid:
        return ""
    parts = tcid.split(".")
    return ".".join(parts[:n_parts])
