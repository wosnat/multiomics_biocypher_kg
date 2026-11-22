"""Adapter emits DE-quality Experiment props and the low_signal_reason edge prop."""
import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.omics_adapter import OMICSAdapter

CONDITION = "condition_at_column_minimum"


def _adapter(tmp_path, analyses, exp_extra=None, df=None, with_exp2=False):
    data_file = tmp_path / "de.csv"
    if df is None:
        df = pd.DataFrame({"gene": ["PMM0001", "PMM0002"], "fc": [1.5, -2.0],
                           "padj": [0.001, 0.001]})
    df.to_csv(data_file, index=False)
    exp = {"name": "T", "organism": "Prochlorococcus MED4", "treatment_condition": "T",
           "control_condition": "C", "omics_type": "RNASEQ", "test_type": "DESeq2",
           "treatment_type": ["nitrogen"], "background_factors": ["axenic"]}
    exp.update(exp_extra or {})
    experiments = {"e1": exp}
    if with_exp2:
        experiments["e2"] = dict(exp, name="no analyses")
    base = {"name_col": "gene", "logfc_col": "fc", "adjusted_p_value_col": "padj",
            "pvalue_threshold": 0.05, "logfc_threshold": 1.0, "experiment": "e1"}
    analyses = [dict(base, id=f"a{i}", **a) for i, a in enumerate(analyses)]
    config = {"publication": {
        "papername": "DEQ", "doi": "10.1234/deq", "control_timing": "matched",
        "experiments": experiments,
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(data_file), "statistical_analyses": analyses}},
    }}
    cfg = tmp_path / "paperconfig.yaml"
    cfg.write_text(yaml.dump(config))
    adapter = OMICSAdapter(config_file=str(cfg))
    adapter.extracted_data = {"title": "T", "authors": [], "journal": "", "year": 2024,
                              "abstract": "", "description": ""}
    return adapter


def _exp_props(adapter):
    return {n[0].split("_", 1)[1]: n[2] for n in adapter.get_nodes() if n[1] == "experiment"}


def test_experiment_carries_caveats_and_fold_change_type(tmp_path):
    ad = _adapter(tmp_path, [{"fold_change_type": "linear"}], exp_extra={
        "data_caveats": ["fold_change_inflated_at_low_counts"],
        "fold_change_provenance": ["computed_by_kg"],
        "data_notes": "x", "data_caveats_reviewed_on": "2026-11-19"})
    p = _exp_props(ad)["e1"]
    assert p["data_caveats"] == ["fold_change_inflated_at_low_counts"]
    assert p["fold_change_provenance"] == ["computed_by_kg"]
    assert p["data_notes"] == "x"
    assert p["data_caveats_reviewed_on"] == "2026-11-19"
    assert p["fold_change_type"] == "linear"


def test_de_experiment_defaults(tmp_path):
    p = _exp_props(_adapter(tmp_path, [{}]))["e1"]
    assert p["data_caveats"] == []
    assert p["fold_change_provenance"] == []
    assert p["fold_change_type"] == "log2"
    assert "data_caveats_reviewed_on" not in p
    assert "low_signal_basis" not in p


def test_non_de_experiment_has_no_de_quality_props(tmp_path):
    p = _exp_props(_adapter(tmp_path, [{}], with_exp2=True))["e2"]
    assert "data_caveats" not in p
    assert "fold_change_provenance" not in p
    assert "fold_change_type" not in p


def test_mixed_fold_change_type_raises(tmp_path):
    ad = _adapter(tmp_path, [{"fold_change_type": "log2"}, {"fold_change_type": "linear"}])
    with pytest.raises(ValueError, match="fold_change_type"):
        ad.get_nodes()


def test_conflicting_low_signal_basis_raises(tmp_path):
    ad = _adapter(tmp_path, [
        {"treatment_level_cols": ["t"], "control_level_cols": ["c"]},
        {"mean_level_cols": ["t"]}],
        df=pd.DataFrame({"gene": ["PMM0001"], "fc": [1.0], "padj": [0.001],
                         "t": [1.0], "c": [1.0]}))
    with pytest.raises(ValueError, match="low_signal_basis|level_basis|basis"):
        ad.get_nodes()


LEVEL_DF = pd.DataFrame({"gene": ["PMM0001", "PMM0002", "PMM0003"],
                         "t": [0.01, 5, 6], "c": [4, 5, 0.01],
                         "fc": [-8.6, 0.0, 9.1], "padj": [0.001] * 3})
LEVELS = {"treatment_level_cols": ["t"], "control_level_cols": ["c"]}


def _reasons(ad):
    return {e[2].split(":")[-1]: e[4].get("low_signal_reason")
            for e in ad.get_edges() if "log2_fold_change" in e[4]}


def test_low_signal_reason_on_edges(tmp_path):
    ad = _adapter(tmp_path, [LEVELS], df=LEVEL_DF)
    r = _reasons(ad)
    assert r["PMM0001"] == CONDITION
    assert r["PMM0002"] is None
    assert r["PMM0003"] == CONDITION
    assert _exp_props(ad)["e1"]["low_signal_basis"] == CONDITION


def test_low_signal_survives_unresolved_rows(tmp_path):
    df = LEVEL_DF.copy()
    df["gene"] = ["", "PMM0002", "PMM0003"]
    r = _reasons(_adapter(tmp_path, [LEVELS], df=df))
    assert "PMM0001" not in r
    assert r["PMM0002"] is None
    assert r["PMM0003"] == CONDITION


def test_missing_level_column_raises(tmp_path):
    ad = _adapter(tmp_path, [{"treatment_level_cols": ["t_typo"], "control_level_cols": ["c"]}],
                  df=LEVEL_DF)
    with pytest.raises(ValueError, match="t_typo"):
        ad.get_edges()


def test_significant_only_keeps_low_signal_reason(tmp_path):
    ad = _adapter(tmp_path, [LEVELS], df=LEVEL_DF)
    ad.significance_mode = "significant_only"
    r = _reasons(ad)
    assert r["PMM0001"] == CONDITION and r["PMM0003"] == CONDITION


def test_one_arm_level_declaration_names_the_analysis(tmp_path):
    ad = _adapter(tmp_path, [{"treatment_level_cols": ["t"]}], df=LEVEL_DF)
    with pytest.raises(ValueError, match=r"a0.*level|level.*a0"):
        ad.get_nodes()


def test_null_data_caveats_is_empty_list(tmp_path):
    ad = _adapter(tmp_path, [{}], exp_extra={"data_caveats": None})
    assert _exp_props(ad)["e1"]["data_caveats"] == []


def test_date_valued_reviewed_on_is_emitted_as_iso_string(tmp_path):
    import datetime
    ad = _adapter(tmp_path, [{}], exp_extra={
        "data_caveats_reviewed_on": datetime.date(2026, 11, 19)})
    assert _exp_props(ad)["e1"]["data_caveats_reviewed_on"] == "2026-11-19"
