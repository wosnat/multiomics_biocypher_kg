"""Cached wrappers around ``bioregistry.normalize_curie`` plus BioCypher-safe text cleanup.

``bioregistry.normalize_curie`` runs a pydantic ``Resource.model_dump`` on every
call. Profiling shows ~3 M calls during a KG build with ~6 M ``model_dump``
invocations dominating ~45% of wall-clock time. The KG only references a small
set of distinct CURIEs (bounded by node count), so an unbounded ``lru_cache``
collapses that cost to a one-time miss per unique input.

``clean_text`` uses ``str.translate`` (single pass with a precomputed table)
instead of chained ``.replace()`` calls.
"""

from __future__ import annotations

import functools

from bioregistry import normalize_curie as _normalize_curie


@functools.lru_cache(maxsize=None)
def normalize_curie(curie: str) -> str | None:
    """Cached pass-through to ``bioregistry.normalize_curie``."""
    return _normalize_curie(curie)


@functools.lru_cache(maxsize=None)
def get_id(prefix: str, identifier: str, sep: str = ":") -> str:
    """Build and normalize a CURIE; fall back to ``{prefix}_{identifier}`` if normalization fails.

    Centralizes the pattern used in every adapter:
    ``normalize_curie(f"{prefix}:{identifier}") or f"{prefix}_{identifier}"``.
    """
    raw = f"{prefix}{sep}{identifier}"
    return _normalize_curie(raw) or f"{prefix}_{identifier}"


# Pre-built translation table — faster than chained .replace().
_BIOCYPHER_SAFE_TRANSLATION = str.maketrans({
    "|": ",",  # pipe is the list-value separator in BioCypher CSV output
    "'": "^",  # single quote wraps string values in BioCypher CSV output
})


def clean_text(text):
    """Strip BioCypher-special characters from a string or list of strings."""
    if isinstance(text, str):
        return text.translate(_BIOCYPHER_SAFE_TRANSLATION)
    if isinstance(text, list):
        return [
            t.translate(_BIOCYPHER_SAFE_TRANSLATION) if isinstance(t, str) else t
            for t in text
        ]
    return text
