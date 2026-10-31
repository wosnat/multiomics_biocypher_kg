"""Fold-change convention checks shared by omics_adapter and validate_paperconfig."""

import math

# A log2 table with no value in (-1, 1) and |x| above this is almost certainly signed linear.
SIGNED_LINEAR_SUSPECT_MAX = 20.0


def _numbers(values) -> list[float]:
    out = []
    for v in values:
        try:
            x = float(str(v).strip().strip('*').strip())
        except (ValueError, TypeError):
            continue
        if math.isfinite(x):
            out.append(x)
    return out


def signed_linear_violations(values) -> list[float]:
    """Values a signed linear table cannot hold: |x| < 1 (no change is ±1).

    Converting one would invert its direction (0.5 -> -1, -0.5 -> +1), so a table
    declared ``signed_linear`` that contains any is mis-declared. Shared by
    ``omics_adapter`` (refuses the analysis) and ``validate_paperconfig.py`` (error).
    Accepts raw cells; stray ``*`` markers are stripped, non-numeric cells ignored.
    """
    return [x for x in _numbers(values) if -1 < x < 1]


def signed_linear_suspect(values) -> float | None:
    """max |x| when a table declared ``log2`` looks signed linear, else None.

    A signed linear table never holds a value in (-1, 1); a log2 table almost always
    does (most genes barely change). No value in (-1, 1) plus a maximum |x| above
    ``SIGNED_LINEAR_SUSPECT_MAX`` (2^20-fold) is the signature. Small-effect signed
    tables (Lindell 2007, max 14) slip under it, so it is a warning, not a proof.
    """
    xs = _numbers(values)
    if not xs or any(-1 < x < 1 for x in xs):
        return None
    top = max(abs(x) for x in xs)
    return top if top > SIGNED_LINEAR_SUSPECT_MAX else None
