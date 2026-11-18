"""MultiKeggModuleAdapter: KeggModule nodes, hierarchy / KO / pathway edges, organism completeness."""
import json

import pytest

from multiomics_kg.adapters.kegg_module_adapter import MultiKeggModuleAdapter

NITRO = ["Pathway modules", "Energy metabolism", "Nitrogen metabolism"]
DRUG = ["Signature modules", "Drug resistance", "Beta-Lactam resistance"]
KEGG_DATA = {
    "modules": {
        "M00901": {"name": "Uses M00902", "module_type": "pathway", "definition": "M00902 K00002",
                   "class_path": NITRO, "pathways": ["ko00910", "ko99999"]},
        "M00902": {"name": "One step", "module_type": "pathway", "definition": "K00001",
                   "class_path": NITRO, "pathways": []},
        "M00903": {"name": "Marker's set", "module_type": "signature",
                   "definition": "K00003+K00004 K00005", "class_path": DRUG, "pathways": []},
        "M00904": {"name": "Deleted KO only", "module_type": "pathway", "definition": "K09999",
                   "class_path": NITRO, "pathways": []},
    },
    "module_ec_to_kos": {"1.1.1.1": ["K00001"]},
    "kegg_release": "pathway 2026-11-18",
}
GENES = {
    "LT_1": {"ec_numbers": ["1.1.1.1"]},                        # EC-only gene
    "LT_2": {"kegg_ko": ["K00002", "K00003"], "ec_numbers": ["1.1.1.1"]},   # has a KO: never EC
    "LT_3": {"kegg_ko": ["K09999"]},                            # deleted KO
    "LT_4": {"kegg_ko": ["K00005"], "ec_numbers": ["1.1.1.-"]},
}
KO_IDS = {"K00001", "K00002", "K00003", "K00004", "K00005"}


@pytest.fixture
def adapter(tmp_path):
    cache = tmp_path / "cache"
    (cache / "kegg").mkdir(parents=True)
    (cache / "kegg" / "kegg_data.json").write_text(json.dumps(KEGG_DATA))
    a_dir, b_dir = tmp_path / "sA", tmp_path / "sB"
    a_dir.mkdir()
    b_dir.mkdir()                                               # no gene_annotations_merged.json
    (a_dir / "gene_annotations_merged.json").write_text(json.dumps(GENES))
    reg = tmp_path / "genomes.csv"
    reg.write_text("ncbi_accession,strain_name,data_dir\n"
                   f"GCF_000000001.1,A,{a_dir}\n"
                   f"#GCF_000000003.1,C,{a_dir}\n"
                   f"GCF_000000002.1,B,{b_dir}\n")
    return MultiKeggModuleAdapter(str(reg), cache, ko_node_ids=KO_IDS, pathway_node_ids={"ko00910"})


def _nodes(a):
    return {n[0]: n for n in a.get_nodes()}


def _edges(a, label):
    return [e for e in a.get_edges() if e[3] == label]


def test_module_and_class_nodes(adapter):
    nodes = _nodes(adapter)
    assert {i for i in nodes if i.startswith("kegg.module:M")} == {
        "kegg.module:M00901", "kegg.module:M00902", "kegg.module:M00903"}   # M00904: deleted KO only
    m = nodes["kegg.module:M00903"]
    assert m[1] == "kegg module"
    assert m[2] == {"name": "Marker^s set", "level": 3, "level_kind": "kegg_module",
                    "module_type": "signature", "module_id": "M00903",
                    "definition": "K00003+K00004 K00005", "step_count": 2}
    sub = nodes["kegg.module:signature_modules.drug_resistance.beta_lactam_resistance"][2]
    assert (sub["level"], sub["level_kind"], sub["module_type"], sub["name"]) == (
        2, "kegg_module_subclass", "signature", "Beta-Lactam resistance")
    top = nodes["kegg.module:pathway_modules"][2]
    assert (top["level"], top["level_kind"], top["module_type"]) == (0, "kegg_module_type", "pathway")
    assert len(nodes) == 3 + 6
    assert adapter.module_node_ids() == set(nodes)


def test_is_a_edges(adapter):
    pairs = {(e[1], e[2]) for e in _edges(adapter, "kegg_module_is_a_kegg_module")}
    assert ("kegg.module:M00902", "kegg.module:pathway_modules.energy_metabolism.nitrogen_metabolism") in pairs
    assert ("kegg.module:pathway_modules.energy_metabolism", "kegg.module:pathway_modules") in pairs
    assert len(pairs) == 3 + 4


def test_ko_edges_dangling_proof(adapter):
    ko = {(e[1], e[2]): e[4] for e in _edges(adapter, "kegg_term_in_kegg_module")}
    assert ko[("kegg.orthology:K00002", "kegg.module:M00901")] == {"steps": [2], "component": "required"}
    assert ko[("kegg.orthology:K00004", "kegg.module:M00903")] == {"steps": [1], "component": "required"}
    assert not any("K09999" in src for src, _ in ko)
    assert len(ko) == 5
    e = _edges(adapter, "kegg_term_in_kegg_module")[0]
    assert e[0] == f"{e[1].split(':')[1]}-in_module-{e[2].split(':')[1]}"


def test_pathway_edges_pruned(adapter):
    assert [(e[1], e[2]) for e in _edges(adapter, "kegg_module_in_kegg_pathway")] == [
        ("kegg.module:M00901", "kegg.pathway:ko00910")]


def test_organism_edges(adapter):
    org = {e[2]: e for e in _edges(adapter, "organism_has_kegg_module")}
    assert set(org) == {"kegg.module:M00901", "kegg.module:M00902", "kegg.module:M00903"}
    assert all(e[1] == "insdc.gcf:GCF_000000001.1" for e in org.values())   # B has no annotations
    p = org["kegg.module:M00902"][4]
    assert (p["status"], p["steps_present_via_ec"], p["ec_genes"], p["ec_filled_step_indices"]) == (
        "complete", 1, ["LT_1"], [1])
    p = org["kegg.module:M00901"][4]
    assert (p["status"], p["ec_genes"]) == ("complete", ["M00902"])
    p = org["kegg.module:M00903"][4]
    assert p == {"status": "one_step_missing", "steps_total": 2, "steps_present": 1,
                 "steps_present_via_ec": 0, "missing_step_indices": [1],
                 "missing_step_definitions": ["K00003+K00004"], "missing_kos": ["K00003", "K00004"]}
    assert org["kegg.module:M00903"][0] == "insdc.gcf:GCF_000000001.1-has_module-M00903"


def test_no_quote_or_pipe_in_strings(adapter):
    def strings(v):
        if isinstance(v, str):
            yield v
        elif isinstance(v, list):
            for x in v:
                yield from strings(x)
    items = [n[2] for n in adapter.get_nodes()] + [e[4] for e in adapter.get_edges()]
    assert all("'" not in s and "|" not in s for props in items for v in props.values()
               for s in strings(v))
