"""Tests for the Phase 1.1B metabolism transforms (validate_tcdb, validate_cazy).

Note: resolve_kegg_reaction_to_mnxr was removed in Spec 1.2 pivot — raw KEGG
R-numbers are now kept as-is (no MNX resolution at build time).
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download.utils import annotation_transforms as at


@pytest.fixture(autouse=True)
def patch_caches(monkeypatch, tmp_path):
    """Set up tiny TCDB/CAZy hierarchies, redirect module caches."""
    # TCDB hierarchy
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1": {}, "3.A.1.1": {}, "3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    # CAZy hierarchy
    from multiomics_kg.utils import cazy_utils as cu
    cazy_path = tmp_path / "cazy_hierarchy.json"
    cazy_path.write_text(json.dumps({"GH": {}, "GH13": {}, "GH13_1": {}, "GT": {}, "GT19": {}}))
    monkeypatch.setattr(cu, "DEFAULT_PATH", cazy_path)
    monkeypatch.setattr(cu, "_CACHE", None)

    yield


def test_resolve_kegg_reaction_transform_removed():
    """Spec 1.2 pivot: KEGG reactions stay as raw R-numbers (no MNX resolution)."""
    assert "resolve_kegg_reaction_to_mnxr" not in at._TRANSFORMS


def test_validate_tcdb_keeps_known():
    assert at._tx_validate_tcdb("3.A.1.1.1") == "3.A.1.1.1"
    assert at._tx_validate_tcdb("3.A.1") == "3.A.1"


def test_validate_tcdb_drops_unknown():
    assert at._tx_validate_tcdb("99.X.99") is None


def test_validate_cazy_keeps_known():
    assert at._tx_validate_cazy("GH13_1") == "GH13_1"
    assert at._tx_validate_cazy("GT19") == "GT19"


def test_validate_cazy_drops_unknown():
    assert at._tx_validate_cazy("XX99") is None
