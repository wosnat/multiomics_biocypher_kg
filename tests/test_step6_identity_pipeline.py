import copy
import io
import json

import pandas as pd
import pytest

from multiomics_kg.utils.metabolite_identity.pipeline import (
    Step6Paths, Step6ValidationError, render_outputs, run_identity_stages)
from tests.step6_fixture import PAPER, TABLE, make_inputs

K = "kegg.compound:"


def _paths(tmp_path):
    return Step6Paths(tmp_path / "kegg_data.json", tmp_path / "tcdb_pruned.json", tmp_path,
                      tmp_path / "chebi.json", tmp_path / "map.json", tmp_path / "lock.json",
                      tmp_path / "report.json")


def _render(tmp_path, inp):
    out = run_identity_stages(inp)
    return out, render_outputs(out, _paths(tmp_path))


def test_canonical_graph_shape(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    kd = out.kegg_data
    assert "C00244" in kd["compounds"] and "chebi:17632" not in kd["additional_compounds"]
    assert kd["compounds"]["C00244"]["identity_basis"] == ["chebi_conjugate"]
    assert "metabolomics" in kd["compounds"]["C00041"]["evidence_sources"]
    assert "C00525" not in kd["compounds"]                                    # measurement-only, dropped
    assert kd["compounds"]["C00078"]["recurated_from_ids"] == ["kegg.compound:C00525"]
    assert "kegg.compound:C00525" not in kd["compounds"]["C00078"]["all_identifiers"]
    assert kd["compounds"]["C00041"]["stereo_form_of"] == ["kegg.compound:C01401"]
    assert out.tcdb_pruned["lumping_pin"]["lumping_tcdb_ids"]                 # seeded
    assert out.class_trees["br08001"]["members"] == [
        {"metabolite_id": "kegg.compound:C00041", "class_id": "kegg.brite:br08001.A1.B1"}]


def test_metabolomics_tag_only_on_final_targets(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    tagged = {k for k, v in out.kegg_data["compounds"].items() if "metabolomics" in v["evidence_sources"]}
    assert tagged == {"C00041", "C00078"}


def test_rerun_with_committed_lock_and_pin_is_byte_identical(tmp_path):
    inp = make_inputs(tmp_path)
    out1, files1 = _render(tmp_path, copy.deepcopy(inp))
    inp2 = make_inputs(tmp_path)
    inp2.lock, inp2.lumping_pin, inp2.today = out1.lock, out1.tcdb_pruned["lumping_pin"], "2099-01-01"
    _, files2 = _render(tmp_path, inp2)
    assert files1 == files2                                                   # frozen_on unchanged


def test_ambiguous_name_fails_and_renders_nothing(tmp_path):
    inp = make_inputs(tmp_path, extra_rows=[("Mystery", None)],
                      extra_candidates={"mystery": {"kegg.compound:C00244", "kegg.compound:C00041"}})
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("ambiguous" in m for m in e.value.errors)


def test_resolved_csv_targets_are_written_compounds(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    for path, text in out.resolved_files.items():
        if path.name.endswith("_resolved.csv"):
            df = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
            for mid in df["metabolite_id"]:
                assert not mid or mid.split(":", 1)[1] in out.kegg_data["compounds"] \
                    or mid in out.kegg_data["additional_compounds"]


# ---- beyond the brief -------------------------------------------------------------------------

def test_inputs_are_not_mutated(tmp_path):
    inp = make_inputs(tmp_path)
    before = copy.deepcopy((inp.table, inp.kegg_rest, inp.subtree_substrates, inp.lock))
    run_identity_stages(inp)
    assert (inp.table, inp.kegg_rest, inp.subtree_substrates, inp.lock) == before


def test_rewrite_and_stereo_on_survivors(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    kd = out.kegg_data
    assert out.tcdb_pruned["subtree_substrates"]["2.A.1"] == [f"{K}C00244", f"{K}C01401"]
    assert kd["compounds"]["C00041"]["enantiomer_of"] == [f"{K}C00133"]
    assert "enantiomer_of" not in kd["compounds"]["C00078"]          # D-Trp dropped (stage L)
    assert kd["compounds"]["C00244"]["evidence_sources"] == ["metabolism", "transport"]
    assert all("mnx_chebi_aliases" not in e for e in kd["compounds"].values())
    assert out.tcdb_pruned["chemical_roles"] and out.tcdb_pruned["subtree_roles"]
    assert out.tcdb_pruned["kept_tcdb_ids"] == ["2.A.1", "3.A.3"]
    assert all(not t["members"] for c, t in out.class_trees.items() if c != "br08001")


def test_lock_mapping_and_report(tmp_path):
    out, files = _render(tmp_path, make_inputs(tmp_path))
    assert out.lock["names"]["alanine"]["id"] == f"{K}C00041"
    assert out.lock["names"]["alanine"]["frozen_on"] == "2026-10-15"
    assert out.lock["names"]["tryptophan"]["id"] == f"{K}C00078"
    dec = out.mapping["decisions"][PAPER][TABLE]
    assert dec["alanine"] == {"metabolite_id": f"{K}C00041", "resolution_method": "dominant_stereo_form",
                              "stereo_assignment": "dominant_form_assumed"}
    assert out.mapping["id_aliases"]["chebi:17632"] == f"{K}C00244"
    assert "kegg.compound:C00525" not in out.mapping["id_aliases"]
    rep = out.report
    remap = {r["old_id"]: r for r in rep["remap"]}
    assert remap["chebi:17632"]["new_id"] == f"{K}C00244" and remap["chebi:17632"]["basis"] == "chebi_conjugate"
    assert remap[f"{K}C00525"] == {"old_id": f"{K}C00525", "new_id": f"{K}C00078",
                                   "basis": "remeasured_elsewhere", "role_id": ""}
    assert remap["chebi:10004"]["basis"] == "removed_role"
    assert remap["chebi:16977"]["basis"] == "rekey_chebi_kegg_xref"
    assert rep["unaccounted_old_ids"] == []
    assert rep["dropped"] == {f"{K}C00525": f"{K}C00078"}
    assert {r["paper_name"]: r["reason"] for r in rep["repoints"]} == {
        "Alanine": "dominant_stereo_form", "Tryptophan": "dominant_stereo_form"}
    assert rep["merge_groups"][f"{K}C00244"]["members"] == ["chebi:17632", f"{K}C00244"]
    # every rendered JSON is indent=2 / sorted / newline-terminated
    for path, text in files.items():
        if path.suffix == ".json" and "resolution_report" not in path.name:
            assert text == json.dumps(json.loads(text), indent=2, sort_keys=True) + "\n"
    assert _paths(tmp_path).chebi_cache in files
    assert tmp_path / "brite_br08001.json" in files


def test_chebi_cache_only_when_built_from_raw(tmp_path):
    inp = make_inputs(tmp_path)
    inp.chebi_from_raw = False
    out, files = _render(tmp_path, inp)
    assert out.chebi_cache is None and _paths(tmp_path).chebi_cache not in files


def test_brite_falls_back_to_committed_tree(tmp_path):
    inp = make_inputs(tmp_path)
    del inp.brite_raw["br08002"]
    inp.brite_committed = {"br08002": {
        "tree_code": "br08002", "tree": "lipids",
        "nodes": [{"id": "kegg.brite:br08002.A1", "parent": None, "name": "Fatty acids", "level": 0},
                  {"id": "kegg.brite:br08002.A2", "parent": None, "name": "Gone", "level": 0}],
        "members": [{"metabolite_id": "kegg.compound:C00244", "class_id": "kegg.brite:br08002.A1"},
                    {"metabolite_id": "kegg.compound:C00525", "class_id": "kegg.brite:br08002.A2"}]}}
    out, _ = _render(tmp_path, inp)
    tree = out.class_trees["br08002"]
    assert tree["members"] == [{"metabolite_id": f"{K}C00244", "class_id": "kegg.brite:br08002.A1"}]
    assert [n["id"] for n in tree["nodes"]] == ["kegg.brite:br08002.A1"]


def test_brite_missing_everywhere_fails(tmp_path):
    inp = make_inputs(tmp_path)
    del inp.brite_raw["br08002"]
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("br08002" in m for m in e.value.errors)


def test_identifier_collision_becomes_validation_error(tmp_path):
    inp = make_inputs(tmp_path)
    # Two unrelated nodes sharing one HMDB id: stage H refuses, surfaced as a validation error.
    inp.table[f"{K}C00078"]["hmdb_id"] = "HMDB0000001"
    inp.table[f"{K}C00133"]["hmdb_id"] = "HMDB0000001"
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("HMDB0000001" in m for m in e.value.errors)


def test_old_chebi_id_rekeyed_onto_dropped_loser_is_linked(tmp_path):
    # Previous build resolved row 1 to chebi:16296, which stage D now re-keys onto D-Trp C00525.
    inp = make_inputs(tmp_path)
    inp.old_resolved[(PAPER, TABLE, 1)] = "chebi:16296"
    inp.rekeys["chebi:16296"] = f"{K}C00525"
    inp.previous_ids.add("chebi:16296")
    out = run_identity_stages(inp)
    rep = out.report
    assert rep["dropped"] == {f"{K}C00525": f"{K}C00078"} and rep["dropped_without_rows"] == []
    assert out.kegg_data["compounds"]["C00078"]["recurated_from_ids"] == [f"{K}C00525"]
    remap = {r["old_id"]: r for r in rep["remap"]}
    assert remap["chebi:16296"] == {"old_id": "chebi:16296", "new_id": f"{K}C00078",
                                    "basis": "remeasured_elsewhere", "role_id": ""}
    assert rep["unaccounted_old_ids"] == []


def test_remap_onto_unwritten_node_fails(tmp_path):
    # No old row moved off D-Trp, so it is dropped without a target; an old id re-keyed onto it
    # would remap to a node step 6 does not write.
    inp = make_inputs(tmp_path)
    del inp.old_resolved[(PAPER, TABLE, 1)]
    inp.rekeys["chebi:16296"] = f"{K}C00525"
    inp.previous_ids.add("chebi:16296")
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("chebi:16296" in m and "C00525" in m for m in e.value.errors)


def test_measured_set_target_not_written_is_validation_error():
    from multiomics_kg.utils.metabolite_identity.paper import Decision, ResolutionResult
    from multiomics_kg.utils.metabolite_identity.pipeline import _measured_set
    table = {f"{K}C1": {"evidence_sources": []}}
    res = ResolutionResult(decisions={("p", "t", 0): Decision(f"{K}C9", "name_match", "as_reported")})
    with pytest.raises(Step6ValidationError) as e:
        _measured_set(table, res, {("p", "t", 0): f"{K}C1"}, lambda x: x, {})
    assert any("C9" in m for m in e.value.errors)


def test_pathway_reprune(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    kd = out.kegg_data
    assert set(kd["pathways"]) == {"ko00250", "ko00380"}
    assert set(kd["subcategories"]) == {"09105"} and set(kd["categories"]) == {"09100"}
    assert kd["kos"]["K00001"]["pathways"] == ["ko00250"]
    assert kd["reactions"]["R3"]["pathways"] == ["ko00380"]
    assert kd["compounds"]["C00078"]["pathways"] == ["ko00380"]


def test_report_sections_and_sparse_filter(tmp_path):
    inp = make_inputs(tmp_path)
    inp.table[f"{K}C00133"]["smiles"] = ""
    inp.table[f"{K}C00133"]["extra"] = {}
    out = run_identity_stages(inp)
    e = out.kegg_data["compounds"]["C00133"]
    assert "smiles" not in e and "extra" not in e
    rep = out.report
    assert rep["structure"]["counts"] == {"defined": 5}               # survivors only
    assert rep["id_aliases"] == out.mapping["id_aliases"]
    assert rep["merge_groups_by_basis"] == {"chebi_conjugate": [f"{K}C00244"]}
