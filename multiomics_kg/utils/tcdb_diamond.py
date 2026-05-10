# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein post-steps.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.
"""
from __future__ import annotations


def classify_hit(hit: dict) -> int | None:
    """Assign a confidence tier (1/2/3) to a parsed diamond hit row.

    Returns None when the hit fails the tier-3 floor (drop it).

    Floor (gblast3-style): e-value <= 0.001, HSP length >= 50, AND
    (qcov >= 40 OR scov >= 40). Above the floor:
      - tier 1: identity >= 70 AND qcov >= 70
      - tier 2: identity >= 40 AND qcov >= 60
      - tier 3: floor only
    """
    if hit["evalue"] > 0.001:
        return None
    if hit["length"] < 50:
        return None
    if hit["qcov"] < 40.0 and hit["scov"] < 40.0:
        return None

    if hit["identity"] >= 70.0 and hit["qcov"] >= 70.0:
        return 1
    if hit["identity"] >= 40.0 and hit["qcov"] >= 60.0:
        return 2
    return 3


def truncate_tcid(tcid: str, n_parts: int) -> str:
    """Truncate a TCID to its first n_parts dot-separated segments.

    A TCID has up to 5 parts (class.subclass.family.subfamily.specificity).
    Returns the input unchanged if it already has <= n_parts segments.
    """
    if not tcid:
        return ""
    parts = tcid.split(".")
    return ".".join(parts[:n_parts])
