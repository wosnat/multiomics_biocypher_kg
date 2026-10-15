import json

from multiomics_kg.adapters.metabolite_class_adapter import MetaboliteClassAdapter


def _write(tmp_path):
    kegg = tmp_path / "kegg"
    kegg.mkdir()
    (kegg / "brite_br08001.json").write_text(json.dumps({
        "tree_code": "br08001", "tree": "compounds_with_biological_roles",
        "nodes": [{"id": "kegg.brite:br08001.A1", "parent": None, "name": "Peptides", "level": 0},
                  {"id": "kegg.brite:br08001.A1.B1", "parent": "kegg.brite:br08001.A1",
                   "name": "Amino acids <i>x</i>", "level": 1}],
        "members": [{"metabolite_id": "kegg.compound:C00041", "class_id": "kegg.brite:br08001.A1.B1"},
                    {"metabolite_id": "kegg.compound:C99999", "class_id": "kegg.brite:br08001.A1.B1"}]}))
    return tmp_path


def test_nodes_edges_and_dangling_proofing(tmp_path):
    a = MetaboliteClassAdapter(cache_root=_write(tmp_path), metabolite_node_ids={"kegg.compound:C00041"})
    a.download_data()
    nodes = {n[0]: n for n in a.get_nodes()}
    assert nodes["kegg.brite:br08001.A1.B1"][1] == "metabolite class"
    assert nodes["kegg.brite:br08001.A1.B1"][2]["level_kind"] == "kegg_compound_subclass"
    edges = list(a.get_edges())
    member = [e for e in edges if e[3] == "metabolite_in_metabolite_class"]
    assert [(e[1], e[4]) for e in member] == [("kegg.compound:C00041", {"sources": ["kegg"], "evidence": "curated"})]
    assert any(e[3] == "metabolite_class_is_a_metabolite_class" and e[2] == "kegg.brite:br08001.A1" for e in edges)
