"""Loader for config/data_source_citations.yaml.

The yaml is the single source of per-DataSource metadata (name, description,
citations, dois, url, license, version) and of the build framework's
citations stamped on ``Schema_info.build_citations``.

Versions are either a literal ``version`` or a ``version_from`` reader over a
committed cache. A missing cache yields ``""`` with a warning (fresh checkout),
never an error.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[2]
CITATIONS_PATH = _REPO_ROOT / "config" / "data_source_citations.yaml"
CHEBI_CACHE_PATH = _REPO_ROOT / "cache" / "data" / "chebi" / "chebi_metabolite_xrefs.json"
GO_RELEASE_PATH = _REPO_ROOT / "cache" / "data" / "go_terms" / "go_release.txt"

# `'` → `^` and `|` deleted by the adapters' _clean_str; `"` is swapped for `'`
# when post-import injects build_citations as a Cypher list; `\` escapes there.
FORBIDDEN_CHARS = ("'", '"', "|", "\\")


def load_citations(path: str | Path | None = None) -> dict[str, Any]:
    """Return ``{"sources": {id: entry}, "build_citations": [str]}``."""
    with open(path or CITATIONS_PATH) as f:
        data = yaml.safe_load(f) or {}
    return {
        "sources": data.get("sources") or {},
        "build_citations": data.get("build_citations") or [],
    }


def _read_chebi_release() -> str:
    if not CHEBI_CACHE_PATH.exists():
        logger.warning(f"ChEBI cache {CHEBI_CACHE_PATH} missing; DataSource chebi gets no version")
        return ""
    return str(json.loads(CHEBI_CACHE_PATH.read_text()).get("release", ""))


def _read_go_release() -> str:
    if not GO_RELEASE_PATH.exists():
        logger.warning(f"{GO_RELEASE_PATH} missing; DataSource go gets no version")
        return ""
    return GO_RELEASE_PATH.read_text().strip()


_VERSION_READERS = {
    "chebi_cache": _read_chebi_release,
    "go_release": _read_go_release,
}


def resolve_version(entry: dict[str, Any]) -> str:
    """Literal ``version``, else the ``version_from`` reader's value, else ``""``."""
    reader = entry.get("version_from")
    if reader is not None:
        if reader not in _VERSION_READERS:
            raise ValueError(f"unknown version_from {reader!r}; known: {sorted(_VERSION_READERS)}")
        return _VERSION_READERS[reader]()
    return str(entry.get("version") or "")
