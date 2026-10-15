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


@pytest.mark.parametrize("raw,clean", [
    ("Enkephalin [KO:K18832]", "Enkephalin"), ("Bradykinin [KO:K03898]", "Bradykinin"),
    ("Carboxylic acids [Fig]", "Carboxylic acids"), ("Lipids [BR:br08002]", "Lipids"),
    ("Glycolysis [PATH:map00010]", "Glycolysis"), ("[2Fe-2S] cluster", "[2Fe-2S] cluster")])
def test_clean_strips_cross_reference_brackets(raw, clean):
    """Final review item 6: any ``[XX:…]`` cross-reference bracket is stripped from class names."""
    from multiomics_kg.utils.metabolite_identity.classes import _clean
    assert _clean(raw) == clean


def test_download_compound_brite_drops_only_the_failed_tree(tmp_path, monkeypatch):
    """One tree's download failure leaves that tree out (committed fallback per tree); the
    trees that downloaded are kept and cached."""
    from multiomics_kg.utils import kegg_utils
    from multiomics_kg.utils.metabolite_identity.classes import COMPOUND_BRITE_TREES, download_compound_brite

    def fetch(url):
        if "br08021" in url:
            raise RuntimeError("HTTP 500")
        return {"name": url}
    monkeypatch.setattr(kegg_utils, "_fetch_json", fetch)
    out = download_compound_brite(tmp_path, force=True)
    assert set(out) == set(COMPOUND_BRITE_TREES) - {"br08021"}
    assert not (tmp_path / "br_br08021.json").exists() and (tmp_path / "br_br08001.json").exists()
