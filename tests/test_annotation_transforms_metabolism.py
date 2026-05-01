"""Tests for the Phase 1.1B metabolism transforms (validate_tcdb).

Note: resolve_kegg_reaction_to_mnxr was removed in Spec 1.2 pivot — raw KEGG
R-numbers are now kept as-is (no MNX resolution at build time).

Note: validate_cazy was removed when CAZy hierarchy was promoted to a
pure-Python in-process table; CAZy IDs are passthrough at this layer and
filtered at the cazy_adapter level.
"""
from __future__ import annotations

import json

import pytest

from multiomics_kg.download.utils import annotation_transforms as at


@pytest.fixture(autouse=True)
def patch_caches(monkeypatch, tmp_path):
    """Set up tiny TCDB hierarchy, redirect module caches."""
    # TCDB hierarchy
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1": {}, "3.A.1.1": {}, "3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    yield


def test_resolve_kegg_reaction_transform_removed():
    """Spec 1.2 pivot: KEGG reactions stay as raw R-numbers (no MNX resolution)."""
    assert "resolve_kegg_reaction_to_mnxr" not in at._TRANSFORMS


def test_validate_tcdb_keeps_known():
    assert at._tx_validate_tcdb("3.A.1.1.1") == "3.A.1.1.1"
    assert at._tx_validate_tcdb("3.A.1") == "3.A.1"


def test_validate_tcdb_drops_unknown():
    assert at._tx_validate_tcdb("99.X.99") is None


def test_validate_cazy_transform_removed():
    """CAZy hierarchy moved to pure-Python utils; transform no longer exists."""
    assert "validate_cazy" not in at._TRANSFORMS
