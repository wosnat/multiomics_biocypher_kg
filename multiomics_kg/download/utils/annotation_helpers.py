"""Shared value-testing and splitting helpers for annotation build scripts."""

from __future__ import annotations

import re
from typing import Any


def _nonempty(value: Any) -> bool:
    """Return True if value is non-None, non-empty, and not the eggnog sentinel '-'."""
    if value is None:
        return False
    if isinstance(value, str):
        s = value.strip()
        return bool(s) and s != "-"
    if isinstance(value, (list, tuple)):
        return any(_nonempty(v) for v in value)
    return True


def _split(value: str, delimiter: str) -> list[str]:
    """Split string and strip whitespace; skip empty tokens and '-' sentinels.

    For comma delimiter, uses smart splitting: does NOT split on ', ' (comma
    followed by space). This correctly handles URL-decoded Cyanorak GFF values
    where internal commas are encoded as %2C, so after decoding they appear as
    ', ' while list separators remain plain ','.
    """
    if not isinstance(value, str) or not value.strip():
        return []
    if delimiter == ",":
        # Smart comma: only split on ',' NOT followed by space
        parts = re.split(r",(?! )", value)
    else:
        parts = value.split(delimiter)
    return [v.strip() for v in parts if v.strip() and v.strip() != "-"]


def _coerce_to_tokens(raw: Any, delimiter: str) -> list[str]:
    """Convert a raw value (str or list) to a flat list of non-empty string tokens."""
    if isinstance(raw, list):
        tokens = []
        for item in raw:
            if _nonempty(item):
                tokens.append(str(item).strip())
        return tokens
    if isinstance(raw, str):
        return _split(raw, delimiter)
    return []
