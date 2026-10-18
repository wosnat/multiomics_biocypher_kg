"""GO transport targets through the pure identity pipeline (stages E-P)."""
import pytest

from multiomics_kg.utils.metabolite_identity.pipeline import Step6ValidationError, run_identity_stages
from tests.metabolite_identity_helpers import cpd
from tests.step6_fixture import make_inputs

K = "kegg.compound:"
ROLE = "chebi.role:33229"


def _go_inputs(tmp_path, **kw):
    inp = make_inputs(tmp_path)
    # chebi:17632 (nitrate) merges onto C00244 at stage E.
    inp.go_substrates = {"GO:0015112": {"chebi:17632": ["CHEBI:17632"], f"{K}C00244": ["CHEBI:48107"]},
                         "GO:0015000": {f"{K}C01401": ["CHEBI:16449"]}}
    inp.go_roles = {"GO:0015343": {ROLE: ["CHEBI:33229"]}}
    inp.go_release = "2026-07-26"
    for k, v in kw.items():
        setattr(inp, k, v)
    return inp


def test_pipeline_rewrites_and_validates_go_targets(tmp_path):
    out = run_identity_stages(_go_inputs(tmp_path))
    kd = out.kegg_data
    assert kd["go_transport_substrates"] == {
        "GO:0015000": {f"{K}C01401": ["CHEBI:16449"]},
        "GO:0015112": {f"{K}C00244": ["CHEBI:17632", "CHEBI:48107"]},
        "GO:0015343": {ROLE: ["CHEBI:33229"]}}
    assert kd["go_transport_release"] == "2026-07-26"
    assert kd["go_transport_excluded_targets"] == []

    bad = _go_inputs(tmp_path)
    bad.go_substrates["GO:0015112"][f"{K}C09999"] = ["CHEBI:9"]
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(bad)
    assert any("go_transport_substrates" in m and "C09999" in m for m in e.value.errors)


def test_go_role_without_chemical_role_fails(tmp_path):
    inp = _go_inputs(tmp_path, go_roles={"GO:1": {"chebi.role:99": ["CHEBI:99"]}})
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("go_transport_substrates" in m and "chebi.role:99" in m for m in e.value.errors)


def test_target_merged_onto_filtered_node_is_dropped(tmp_path):
    # The filter names chebi:17632; it canonicalises onto C00244, so both of GO:0015112's
    # targets are dropped after canonicalisation; an unresolved filter entry keeps its ChEBI id.
    inp = _go_inputs(tmp_path, go_excluded_primaries={"chebi:17632", "chebi:15378"},
                     go_filter_primaries={"chebi:17632": "chebi:17632", "chebi:15378": None})
    out = run_identity_stages(inp)
    kd = out.kegg_data
    assert "GO:0015112" not in kd["go_transport_substrates"]
    assert kd["go_transport_excluded_targets"] == ["chebi:15378", f"{K}C00244"]
    assert out.report["go_transport"]["drops"]["filtered_after_canonical"] == 2


def test_report_block_and_distribution(tmp_path):
    inp = _go_inputs(tmp_path, go_distribution_input={
        "s1": [["GO:0015112"], ["GO:0015112", "GO:0015343"], []], "s2": [["GO:0015000"]]},
        go_drops={"filtered": 3, "unresolved": 1}, go_axioms_in_corpus=4)
    rep = run_identity_stages(inp).report["go_transport"]
    assert rep["release"] == "2026-07-26"
    assert rep["axioms_in_corpus"] == 4
    assert rep["edges"] == 3 and rep["targets"] == 3 and rep["role_targets"] == 1
    assert rep["drops"]["filtered"] == 3 and rep["drops"]["unresolved"] == 1
    rows = {r["target"]: r for r in rep["distribution"]}
    assert rows[f"{K}C00244"]["genes"] == 2 and rows[f"{K}C00244"]["strains"] == 1
    assert rows[ROLE]["genes"] == 1
    assert rows[f"{K}C01401"]["strains"] == 1
    assert rep["new_metabolites"] == 0


def test_pipeline_no_pathway_pull_in(tmp_path):
    inp = _go_inputs(tmp_path)
    inp.table[f"{K}C00800"] = cpd("GO-only compound", ev=("go_transport",), pathways=("ko00999",))
    inp.compound_to_pathways["C00800"] = ["ko00999"]
    inp.go_substrates["GO:0015999"] = {f"{K}C00800": ["CHEBI:800"]}
    out = run_identity_stages(inp)
    kd = out.kegg_data
    assert "C00800" in kd["compounds"]
    assert kd["compounds"]["C00800"]["evidence_sources"] == ["go_transport"]
    assert not kd["compounds"]["C00800"].get("pathways")
    assert "ko00999" not in kd["pathways"]
    assert out.report["go_transport"]["new_metabolites"] == 1


def test_no_go_inputs_still_runs(tmp_path):
    out = run_identity_stages(make_inputs(tmp_path))
    assert out.kegg_data["go_transport_substrates"] == {}
