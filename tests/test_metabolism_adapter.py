"""Unit tests for MetabolismAdapter and MultiMetabolismAdapter."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.adapters import metabolism_adapter as ma


XREFS_FIXTURE = {
    "reactions": {
        "R00200": {
            "name": "pyruvate kinase reaction",
            "compound_ids": ["C00074", "C00008"],
            "kegg_pathway_ids": ["ko00010", "ko00710"],
            "ec_numbers": ["2.7.1.40"],
            "mnxr_id": "MNXR101234",
            "rhea_ids": ["10828"],
            "mass_balance": "balanced",
            "reaction_class": "chemical",
        },
        "R12345": {
            "name": "test reaction",
            "compound_ids": ["C99999"],
            "kegg_pathway_ids": [],
            "ec_numbers": [],
            "mnxr_id": None,
            "rhea_ids": [],
            "mass_balance": "unbalanced",
            "reaction_class": "chemical",
        },
    },
    "compounds": {
        "C00074": {
            "name": "phosphoenolpyruvate",
            "formula": "C3H5O6P",
            "mass": 168.04,
            "inchikey": "DTBNBXWJWCWCIK-UHFFFAOYSA-N",
            "smiles": "OC(=O)C(=C)OP(O)(O)=O",
            "mnxm_id": "MNXM73",
            "chebi_id": "44897",
            "hmdb_id": "HMDB0000263",
        },
        "C00008": {
            "name": "ADP",
            "formula": "C10H15N5O10P2",
            "mass": 427.20,
            "inchikey": "XTWYTFMLZFPYCI-KQYNXXCUSA-N",
            "smiles": "Nc1ncnc2c1ncn2[CH]1O[CH]...",
            "mnxm_id": "MNXM7",
            "chebi_id": "16761",
            "hmdb_id": "HMDB0001341",
        },
        "C99999": {
            "name": "obscure compound",
            "formula": None,
            "mass": None,
            "inchikey": None,
            "smiles": None,
            "mnxm_id": None,
            "chebi_id": None,
            "hmdb_id": None,
        },
    },
}


@pytest.fixture
def xrefs_file(tmp_path):
    p = tmp_path / "kegg_metabolism_xrefs.json"
    p.write_text(json.dumps(XREFS_FIXTURE))
    return p


def test_reaction_nodes_have_kegg_primary_id(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    nodes = list(adapter.get_nodes())

    rxn_ids = [n[0] for n in nodes if n[1] == "reaction"]
    assert "kegg.reaction:R00200" in rxn_ids
    assert "kegg.reaction:R12345" in rxn_ids


def test_reaction_node_properties(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    rxn = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.reaction:R00200" and n[1] == "reaction")
    _, _, props = rxn
    assert props["kegg_reaction_id"] == "R00200"
    assert props["name"] == "pyruvate kinase reaction"
    assert props["ec_numbers"] == ["2.7.1.40"]
    assert props["kegg_pathway_ids"] == ["ko00010", "ko00710"]
    assert props["mnxr_id"] == "MNXR101234"
    assert props["rhea_ids"] == ["10828"]
    assert props["mass_balance"] == "balanced"
    assert props["reaction_class"] == "chemical"


def test_metabolite_nodes_have_kegg_primary_id(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    nodes = list(adapter.get_nodes())
    cpd_ids = [n[0] for n in nodes if n[1] == "metabolite"]
    assert "kegg.compound:C00074" in cpd_ids
    assert "kegg.compound:C99999" in cpd_ids


def test_metabolite_node_properties_drop_null_fields(xrefs_file):
    """Sparse properties: nulls should be omitted from the property dict."""
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    cpd = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.compound:C99999" and n[1] == "metabolite")
    _, _, props = cpd
    assert props["kegg_compound_id"] == "C99999"
    assert props["name"] == "obscure compound"
    # All other fields are null in the fixture; should be absent from the dict
    assert "formula" not in props
    assert "mnxm_id" not in props


def test_metabolite_full_properties_preserved(xrefs_file):
    adapter = ma.MetabolismAdapter(xrefs_path=xrefs_file)
    cpd = next(n for n in adapter.get_nodes()
               if n[0] == "kegg.compound:C00074" and n[1] == "metabolite")
    _, _, props = cpd
    assert props["formula"] == "C3H5O6P"
    assert props["mass"] == 168.04
    assert props["mnxm_id"] == "MNXM73"
    assert props["chebi_id"] == "44897"
    assert props["hmdb_id"] == "HMDB0000263"
