"""Fold-change convention checks shared by omics_adapter and validate_paperconfig."""

import math


def signed_linear_violations(values) -> list[float]:
    """Values a signed linear table cannot hold: |x| < 1 (no change is ±1).

    Converting one would invert its direction (0.5 -> -1, -0.5 -> +1), so a table
    declared ``signed_linear`` that contains any is mis-declared. Shared by
    ``omics_adapter`` (refuses the analysis) and ``validate_paperconfig.py`` (error).
    Accepts raw cells; stray ``*`` markers are stripped, non-numeric cells ignored.
    """
    bad = []
    for v in values:
        try:
            x = float(str(v).strip().strip('*').strip())
        except (ValueError, TypeError):
            continue
        if math.isfinite(x) and -1 < x < 1:
            bad.append(x)
    return bad
