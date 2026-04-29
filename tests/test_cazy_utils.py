"""Tests for utils/cazy_utils.py — CAZy hierarchy accessor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.utils import cazy_utils as cu


TINY_CAZY = {
    "GH":     {"name": "Glycoside Hydrolases", "level": 0, "level_kind": "cazy_class",     "parent": None,    "class": "GH"},
    "GT":     {"name": "GlycosylTransferases", "level": 0, "level_kind": "cazy_class",     "parent": None,    "class": "GT"},
    "GH13":   {"name": "",                     "level": 1, "level_kind": "cazy_family",    "parent": "GH",    "class": "GH"},
    "GH13_1": {"name": "",                     "level": 2, "level_kind": "cazy_subfamily", "parent": "GH13",  "class": "GH"},
    "GT19":   {"name": "",                     "level": 1, "level_kind": "cazy_family",    "parent": "GT",    "class": "GT"},
}


@pytest.fixture(autouse=True)
def patch_cazy(monkeypatch, tmp_path):
    p = tmp_path / "cazy_hierarchy.json"
    p.write_text(json.dumps(TINY_CAZY))
    monkeypatch.setattr(cu, "DEFAULT_PATH", p)
    monkeypatch.setattr(cu, "_CACHE", None)
    yield


def test_load_cazy_returns_dict():
    h = cu.load_cazy()
    assert isinstance(h, dict)
    assert h["GH"]["class"] == "GH"


def test_is_valid_cazy_present():
    assert cu.is_valid_cazy("GH") is True
    assert cu.is_valid_cazy("GH13") is True
    assert cu.is_valid_cazy("GH13_1") is True


def test_is_valid_cazy_absent():
    assert cu.is_valid_cazy("XX99") is False
    assert cu.is_valid_cazy("GH99999") is False  # plausible format but not in hierarchy
    assert cu.is_valid_cazy("") is False


def test_cazy_ancestors_subfamily():
    assert cu.cazy_ancestors("GH13_1") == ["GH", "GH13"]


def test_cazy_ancestors_family():
    assert cu.cazy_ancestors("GH13") == ["GH"]


def test_cazy_ancestors_class():
    assert cu.cazy_ancestors("GH") == []


def test_cazy_ancestors_unknown():
    assert cu.cazy_ancestors("XX99") == []
