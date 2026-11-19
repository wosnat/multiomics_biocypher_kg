"""validate_de_quality: data_caveats, the unreviewed allowlist, level columns, scale."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from validate_paperconfig import validate_de_quality, load_caveat_allowlist  # noqa: E402

REVIEWED = {"data_caveats": [], "data_caveats_reviewed_on": "2026-11-18"}


def _run(exp_extra=None, analyses=None, allowlist=None, filename="missing.csv", de=True):
    exp = {"name": "T", **REVIEWED, **(exp_extra or {})}
    exp = {k: v for k, v in exp.items() if v is not Ellipsis}
    pub = {"experiments": {"e1": exp}}
    if analyses is None:
        analyses = [{"id": "a1", "experiment": "e1", "logfc_col": "lfc"}]
    supp = {"s1": {"type": "csv", "filename": filename,
                   "statistical_analyses": analyses if de else []}}
    return validate_de_quality(pub, supp, "x.yaml", allowlist or {})


def _only(msgs, sub):
    assert any(sub in m for m in msgs), msgs
    assert len(msgs) == 1, msgs


def test_clean_reviewed_experiment():
    assert _run() == ([], [])


def test_unknown_caveat_value():
    _only(_run({"data_caveats": ["nope"], "data_caveat_note": "n"})[0], "data_caveats")


def test_probe_caveat_rejected():
    errs, _ = _run({"data_caveats": ["paralogs_share_one_array_probe"], "data_caveat_note": "n"})
    _only(errs, "paralogs_share_one_array_probe")


def test_caveats_need_note():
    _only(_run({"data_caveats": ["no_biological_replicates"]})[0], "data_caveat_note")


def test_note_forbids_pipe_and_quote():
    for bad in ("a|b", "it's"):
        errs, _ = _run({"data_caveats": ["no_biological_replicates"], "data_caveat_note": bad})
        _only(errs, "must not contain")


def test_reviewed_on_must_be_iso_date():
    _only(_run({"data_caveats_reviewed_on": "2026-13-01"})[0], "data_caveats_reviewed_on")


def test_new_de_experiment_requires_review_fields():
    errs, _ = _run({"data_caveats": Ellipsis, "data_caveats_reviewed_on": Ellipsis})
    assert any("missing data_caveats " in e for e in errs)
    assert any("missing data_caveats_reviewed_on" in e for e in errs)
    assert len(errs) == 2


def test_allowlisted_experiment_may_omit_review_fields():
    errs, warns = _run({"data_caveats": Ellipsis, "data_caveats_reviewed_on": Ellipsis},
                       allowlist={"x.yaml": {"e1"}})
    assert errs == [] and warns == []


def test_stale_allowlist_entry_errors():
    _only(_run(allowlist={"x.yaml": {"e1"}})[0], "remove it from config/data_caveats_unreviewed.yaml")


def test_stale_allowlist_key_warns():
    _, warns = _run({"data_caveats": Ellipsis, "data_caveats_reviewed_on": Ellipsis},
                    allowlist={"x.yaml": {"e1", "gone"}})
    _only(warns, "gone")


def test_one_arm_level_cols_errors():
    a = [{"id": "a1", "experiment": "e1", "treatment_level_cols": ["t1"]}]
    _only(_run(analyses=a)[0], "must be given together")


def test_missing_level_column_errors(tmp_path):
    f = tmp_path / "de.csv"
    pd.DataFrame({"g": [1], "t1": [1], "c1": [1]}).to_csv(f, index=False)
    a = [{"id": "a1", "experiment": "e1", "treatment_level_cols": ["t1"],
          "control_level_cols": ["c_missing"]}]
    _only(_run(analyses=a, filename=str(f))[0], "c_missing")


def test_mixed_fold_change_type_errors():
    a = [{"id": "a1", "experiment": "e1"},
         {"id": "a2", "experiment": "e1", "fold_change_type": "linear"}]
    _only(_run(analyses=a)[0], "fold_change_type")


def test_mixed_level_basis_errors():
    a = [{"id": "a1", "experiment": "e1", "mean_level_cols": ["m"]},
         {"id": "a2", "experiment": "e1", "treatment_level_cols": ["t"], "control_level_cols": ["c"]}]
    _only(_run(analyses=a)[0], "level bases")


def test_caveats_on_non_de_experiment_warns():
    errs, warns = _run({"data_caveats": ["no_biological_replicates"], "data_caveat_note": "n"},
                       de=False)
    assert errs == []
    _only(warns, "no statistical_analyses")


def test_missing_allowlist_file_is_empty(tmp_path):
    assert load_caveat_allowlist(tmp_path / "none.yaml") == {}


# --- DE QA report -----------------------------------------------------------
from validate_paperconfig import de_quality_report, validate_paperconfig_content  # noqa: E402


def _report_config(tmp_path, scope="", detail=None, name="t.csv"):
    df = pd.DataFrame({
        "gene": [f"g{i}" for i in range(6)],
        "lfc": [-11.0, 1.2, 0.3, -0.4, 2.0, 0.1],
        "padj": [0.5, 0.01, 0.9, 0.8, 0.001, 0.7],
        "baseMean": [1, 50, 60, 70, 80, 90],
    })
    path = tmp_path / name
    df.to_csv(path, index=False)
    exp = {"name": "T", "table_scope": scope, **REVIEWED}
    if detail:
        exp["table_scope_detail"] = detail
    return {"publication": {
        "experiments": {"e1": exp},
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(path),
            "statistical_analyses": [{
                "id": "a1", "experiment": "e1", "name_col": "gene", "logfc_col": "lfc",
                "adjusted_p_value_col": "padj", "pvalue_threshold": 0.05,
                "logfc_threshold": 1.0, "timepoint": "24h"}]}}}}


def test_report_lists_statistics_and_prompts(tmp_path):
    out = de_quality_report(_report_config(tmp_path), str(tmp_path / "p.yaml"))
    assert "QA a1 (e1, tp 24h)" in out
    assert "extreme_fc_count: 1" in out
    assert "fold_change_inflated_at_low_counts" in out
    assert "candidate level columns: ['baseMean']" in out
    assert "never fails validation" in out


def test_report_prompts_for_modified_table_and_scope_count(tmp_path):
    cfg = _report_config(tmp_path, scope="all_detected_genes",
                         detail="40 significant genes", name="t_modified.csv")
    out = de_quality_report(cfg, str(tmp_path / "p.yaml"))
    assert "fold_change_recomputed_by_kg" in out
    assert "fold_change_sign_reversed_by_kg" in out
    assert "significance_calls_differ_from_paper" in out


def test_report_never_adds_errors(tmp_path):
    cfg = _report_config(tmp_path)
    cfg_path = str(tmp_path / "p.yaml")
    before = validate_paperconfig_content(cfg, cfg_path)
    de_quality_report(cfg, cfg_path)
    assert validate_paperconfig_content(cfg, cfg_path) == before


def test_report_uses_source_csv_without_resolved(tmp_path):
    cfg = _report_config(tmp_path)
    assert not list(tmp_path.glob("*_resolved.csv"))
    assert "rows: 6" in de_quality_report(cfg, str(tmp_path / "p.yaml"))
