"""Pure-function CAZy utility tests (post-cazy_hierarchy.json removal)."""
from __future__ import annotations

from multiomics_kg.utils.cazy_utils import (
    CAZY_CLASSES,
    cazy_ancestors,
    is_valid_cazy,
    parse_cazy_id,
)


def test_cazy_classes_table_has_six_classes():
    assert set(CAZY_CLASSES.keys()) == {"GH", "GT", "PL", "CE", "AA", "CBM"}
    assert CAZY_CLASSES["GH"] == "Glycoside Hydrolases"


def test_parse_cazy_id_family_no_subfamily():
    assert parse_cazy_id("GH13") == ("GH13", None)


def test_parse_cazy_id_with_subfamily():
    assert parse_cazy_id("GH13_5") == ("GH13", "GH13_5")


def test_parse_cazy_id_cbm():
    assert parse_cazy_id("CBM48") == ("CBM48", None)


def test_parse_cazy_id_malformed_returns_none():
    assert parse_cazy_id("GH") is None
    assert parse_cazy_id("XYZ12") is None
    assert parse_cazy_id("") is None
    assert parse_cazy_id("GH13_") is None


def test_parse_cazy_id_strips_whitespace():
    assert parse_cazy_id("  GH13  ") == ("GH13", None)


def test_is_valid_cazy_recognizes_class_family_subfamily():
    assert is_valid_cazy("GH")
    assert is_valid_cazy("GH13")
    assert is_valid_cazy("GH13_5")


def test_is_valid_cazy_rejects_garbage():
    assert not is_valid_cazy("not-a-cazy")
    assert not is_valid_cazy("")
    assert not is_valid_cazy("GH-13")


def test_cazy_ancestors_subfamily():
    assert cazy_ancestors("GH13_5") == ["GH", "GH13"]


def test_cazy_ancestors_family():
    assert cazy_ancestors("GH13") == ["GH"]


def test_cazy_ancestors_class():
    assert cazy_ancestors("GH") == []


def test_cazy_ancestors_unknown_returns_empty():
    assert cazy_ancestors("not-a-cazy") == []
