"""Tests for utils/tcdb_utils.py — TCDB hierarchy accessor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.utils import tcdb_utils as tu


TINY_TCDB = {
    "1":         {"name": "Channels and Pores", "level": 0, "level_kind": "tc_class", "parent": None},
    "1.A":       {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
    "1.A.1":     {"name": "VIC Family", "level": 2, "level_kind": "tc_family", "parent": "1.A", "abbreviation": "VIC"},
    "1.A.1.1":   {"name": "", "level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.1.1": {"name": "", "level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.1",
                  "substrate_classes": ["potassium(1+)"]},
}


@pytest.fixture(autouse=True)
def patch_tcdb(monkeypatch, tmp_path):
    """Point load_tcdb at a tiny synthetic JSON; reset the module cache."""
    p = tmp_path / "tcdb_hierarchy.json"
    p.write_text(json.dumps(TINY_TCDB))
    monkeypatch.setattr(tu, "DEFAULT_PATH", p)
    monkeypatch.setattr(tu, "_CACHE", None)
    yield


def test_load_tcdb_returns_dict():
    h = tu.load_tcdb()
    assert isinstance(h, dict)
    assert h["1"]["level"] == 0


def test_load_tcdb_caches():
    """Second call returns the same object without re-reading the file."""
    a = tu.load_tcdb()
    b = tu.load_tcdb()
    assert a is b


def test_is_valid_tcdb_present():
    assert tu.is_valid_tcdb("1.A.1.1.1") is True
    assert tu.is_valid_tcdb("1.A") is True
    assert tu.is_valid_tcdb("1") is True


def test_is_valid_tcdb_absent():
    assert tu.is_valid_tcdb("99.X.99") is False
    assert tu.is_valid_tcdb("") is False


def test_tcdb_ancestors_full_chain():
    """Specificity-level TCID → list of all ancestors (root → parent)."""
    assert tu.tcdb_ancestors("1.A.1.1.1") == ["1", "1.A", "1.A.1", "1.A.1.1"]


def test_tcdb_ancestors_partial():
    assert tu.tcdb_ancestors("1.A.1") == ["1", "1.A"]
    assert tu.tcdb_ancestors("1") == []


def test_tcdb_ancestors_unknown_returns_empty():
    assert tu.tcdb_ancestors("99.X.99") == []


# ============================================================================
# Cross-ontology bridges (Pfam / GO -> TcdbFamily)
# ============================================================================

from multiomics_kg.utils.tcdb_utils import build_tcdb_bridges

# 1.A.1.5.2 is kept (gene-annotated); 1.A.1.9.1 is not, and rolls up to 1.A.1.
_H = {
    "1": {"level": 0, "level_kind": "tc_class", "parent": None},
    "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
    "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
    "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.5.2": {"level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.5"},
    "1.A.1.9": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.9.1": {"level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.9"},
    "2": {"level": 0, "level_kind": "tc_class", "parent": None},
    "2.A": {"level": 1, "level_kind": "tc_subclass", "parent": "2"},
    "2.A.9.1.1": {"level": 4, "level_kind": "tc_specificity", "parent": "2.A"},
}
_KEPT = {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2", "2", "2.A"}


def test_bridge_attaches_directly_when_tcid_survives():
    b = build_tcdb_bridges({"PF00001": {"1.A.1.5.2"}}, _H, _KEPT)
    assert b == {"1.A.1.5.2": {"PF00001": ["1.A.1.5.2"]}}


def test_bridge_rolls_up_to_nearest_surviving_ancestor():
    """1.A.1.9.1 has no node; its nearest kept ancestor is the family 1.A.1."""
    b = build_tcdb_bridges({"PF00002": {"1.A.1.9.1"}}, _H, _KEPT)
    assert b == {"1.A.1": {"PF00002": ["1.A.1.9.1"]}}


def test_bridge_drops_targets_shallower_than_family():
    """2.A.9.1.1 rolls up only as far as the subclass 2.A -> uninformative, dropped.

    An edge saying "this domain relates to Electrochemical Potential-driven
    Transporters" is true but carries no usable signal.
    """
    assert build_tcdb_bridges({"PF00003": {"2.A.9.1.1"}}, _H, _KEPT) == {}


def test_bridge_records_all_curated_tcids_after_rollup():
    """Roll-up must not lose provenance: every original TCID is retained."""
    b = build_tcdb_bridges({"PF00004": {"1.A.1.9.1", "1.A.1.5.2"}}, _H, _KEPT)
    assert b["1.A.1"]["PF00004"] == ["1.A.1.9.1"]
    assert b["1.A.1.5.2"]["PF00004"] == ["1.A.1.5.2"]


def test_bridge_curated_tcids_are_sorted_and_deduped():
    b = build_tcdb_bridges({"PF00005": {"1.A.1.9.1"}}, _H, _KEPT | {"1.A.1.9"})
    assert b["1.A.1.9"]["PF00005"] == ["1.A.1.9.1"]


def test_bridge_skips_tcids_absent_from_the_hierarchy():
    assert build_tcdb_bridges({"PF00006": {"9.Z.99.9.9"}}, _H, _KEPT) == {}


def test_bridge_handles_go_ids_identically():
    """The builder is xref-agnostic — GO ids roll up the same way as Pfam."""
    b = build_tcdb_bridges({"GO:0005215": {"1.A.1.9.1"}}, _H, _KEPT)
    assert b == {"1.A.1": {"GO:0005215": ["1.A.1.9.1"]}}


def test_bridge_empty_input_returns_empty():
    assert build_tcdb_bridges({}, _H, _KEPT) == {}
