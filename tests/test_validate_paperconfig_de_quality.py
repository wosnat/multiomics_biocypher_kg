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
