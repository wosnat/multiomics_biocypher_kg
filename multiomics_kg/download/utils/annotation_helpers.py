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


def extract_first_match_in_sources(
    sources: list[dict],
    gm: dict,
    eg: dict,
    up: dict,
    pattern: str,
    extract_group: int = 0,
) -> str | None:
    """Return first token (across sources) whose full string matches ``pattern``.

    Tries each source in order; within a source iterates tokens left-to-right.
    Returns the captured group ``extract_group`` (0 = full match).

    Intended for config entries of type ``extract_first_match``.
    """
    compiled = re.compile(pattern)
    for src_cfg in sources:
        source = src_cfg.get("source", "")
        field = src_cfg.get("field", "")
        if source == "gene_mapping":
            raw = gm.get(field)
        elif source == "eggnog":
            raw = eg.get(field)
        elif source == "uniprot":
            raw = up.get(field)
        else:
            continue
        if isinstance(raw, str):
            raw = raw.strip()
        if not _nonempty(raw):
            continue
        delimiter = src_cfg.get("delimiter", ",")
        tokens = _coerce_to_tokens(raw, delimiter)
        for tok in tokens:
            m = compiled.match(tok.strip())
            if m:
                return m.group(extract_group)
    return None


def _coerce_to_tokens(raw: Any, delimiter: str) -> list[str]:
    """Convert a raw value (str or list) to a flat list of non-empty string tokens.

    For list inputs, each item is further split by the delimiter so that
    space-concatenated entries like ['dnaN BFV95_0002'] (common in UniProt
    gene_names) are expanded when delimiter=' '.
    """
    if isinstance(raw, list):
        tokens = []
        for item in raw:
            if _nonempty(item):
                tokens.extend(_split(str(item).strip(), delimiter))
        return tokens
    if isinstance(raw, str):
        return _split(raw, delimiter)
    return []
