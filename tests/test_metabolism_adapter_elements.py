"""Unit tests for the _parse_elements helper in metabolism_adapter."""
from __future__ import annotations

import pytest

from multiomics_kg.adapters.metabolism_adapter import _parse_elements


def test_parse_elements_water():
    assert _parse_elements("H2O") == ["H", "O"]


def test_parse_elements_complex_phosphate():
    assert _parse_elements("C10H12N5O13P3") == ["C", "H", "N", "O", "P"]


def test_parse_elements_two_letter_no_substring_clash():
    """Hill parsing must not split 'Na' into ['N','a'] or 'Cl' into ['C','l']."""
    assert _parse_elements("NaCl") == ["Cl", "Na"]


def test_parse_elements_kegg_charge_suffix():
    """KEGG uses '*N' suffix for charge/radical state; chemparse must tolerate it
    (parse the elements correctly, ignore the trailing *N)."""
    elts = _parse_elements("C42H44FeN8O8S2*4")
    assert "Fe" in elts
    assert "C" in elts
    assert "S" in elts
    assert "*" not in elts
    # Sorted, unique
    assert elts == sorted(set(elts))


def test_parse_elements_iron_sulfur_cluster():
    """Another KEGG-style charge suffix case."""
    elts = _parse_elements("Fe2S2*8")
    assert "Fe" in elts
    assert "S" in elts


def test_parse_elements_null_returns_empty():
    assert _parse_elements(None) == []


def test_parse_elements_empty_string_returns_empty():
    assert _parse_elements("") == []


def test_parse_elements_malformed_returns_empty():
    """Garbage input must not crash the adapter."""
    assert _parse_elements("???") == []


def test_parse_elements_returns_sorted_list():
    """Result is always sorted alphabetically by element symbol."""
    result = _parse_elements("ZnO")
    assert result == sorted(result)


def test_metabolite_yielder_emits_elements_when_formula_present(tmp_path):
    """Smoke test: a metabolite with formula gets an 'elements' key in props."""
    import json
    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter

    # Minimal kegg_data.json fixture
    fake_data = {
        "reactions": {},
        "compounds": {
            "C00001": {"name": "H2O", "formula": "H2O"},
            "C99999": {"name": "Generic", "formula": None},
        },
    }
    fake_path = tmp_path / "kegg_data.json"
    fake_path.write_text(json.dumps(fake_data))

    adapter = MetabolismAdapter(kegg_data_path=fake_path, test_mode=False)
    nodes = list(adapter.get_nodes())

    # Find the H2O node — it should have elements=["H", "O"]
    h2o = next((p for nid, lbl, p in nodes if "C00001" in nid), None)
    assert h2o is not None
    assert h2o.get("elements") == ["H", "O"]

    # Find the no-formula node — elements should be absent (sparse, dropped by _drop_nulls)
    generic = next((p for nid, lbl, p in nodes if "C99999" in nid), None)
    assert generic is not None
    assert "elements" not in generic
