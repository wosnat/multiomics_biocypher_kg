import pytest

from multiomics_kg.utils.metabolite_identity.classes import parse_compound_brite, prune_compound_brite

RAW = {"name": "br08001", "children": [
    {"name": "Organic acids", "children": [
        {"name": "Carboxylic acids [Fig]", "children": [
            {"name": "Monocarboxylic acids", "children": [
                {"name": "C00058  Formate; Methanoate"}, {"name": "C00033  Acetate"}]}]}]},
    {"name": "Peptides", "children": [
        {"name": "Amino acids", "children": [{"name": "C00041  L-Alanine"}, {"name": "C00033  Acetate"}]}]}]}


def test_parse_ids_names_levels():
    nodes, members = parse_compound_brite("br08001", RAW)
    assert nodes["kegg.brite:br08001.A1"] == {"parent": None, "name": "Organic acids", "level": 0}
    assert nodes["kegg.brite:br08001.A1.B1"]["name"] == "Carboxylic acids"
    assert ("C00041", "kegg.brite:br08001.A2.B1") in members
    assert [m for m in members if m[0] == "C00033"] == [("C00033", "kegg.brite:br08001.A1.B1.C1"),
                                                        ("C00033", "kegg.brite:br08001.A2.B1")]


def test_parse_fails_loud_on_empty():
    with pytest.raises(ValueError):
        parse_compound_brite("br08001", {"name": "br08001", "children": [{"name": "X", "children": []}]})


def test_prune_ancestors_only_and_alias_lands_on_canonical():
    nodes, members = parse_compound_brite("br08001", RAW)
    out = prune_compound_brite("br08001", nodes, members,
                               canonical_of={"kegg.compound:C00058": "kegg.compound:C00033"},
                               written={"kegg.compound:C00033"})
    assert {n["id"] for n in out["nodes"]} == {"kegg.brite:br08001.A1", "kegg.brite:br08001.A1.B1",
                                                "kegg.brite:br08001.A1.B1.C1", "kegg.brite:br08001.A2",
                                                "kegg.brite:br08001.A2.B1"}
    assert {"metabolite_id": "kegg.compound:C00033", "class_id": "kegg.brite:br08001.A1.B1.C1"} in out["members"]
    assert len(out["members"]) == len({(m["metabolite_id"], m["class_id"]) for m in out["members"]})
    assert out["tree"] == "compounds_with_biological_roles"
