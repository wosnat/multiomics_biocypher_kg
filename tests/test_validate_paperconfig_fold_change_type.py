"""validate_paperconfig.py: fold_change_type values and the signed_linear full-column check."""
import sys
from pathlib import Path

import pandas as pd

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from validate_paperconfig import validate_paperconfig_content  # noqa: E402


def _config(tmp_path, ratios, fold_change_type):
    data_file = tmp_path / "de.csv"
    pd.DataFrame({
        "Synonym": [f"PMM{i:04d}" for i in range(len(ratios))],
        "ratio": ratios,
        "padj": [0.001] * len(ratios),
    }).to_csv(data_file, index=False)
    config = {"publication": {
        "papername": "Signed FC", "doi": "10.1234/signed.fc", "control_timing": "matched",
        "experiments": {"e1": {
            "name": "T", "organism": "Prochlorococcus MED4", "treatment_condition": "T",
            "control_condition": "C", "omics_type": "RNASEQ", "test_type": "DESeq2",
            "treatment_type": ["nitrogen"], "background_factors": ["axenic"]}},
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(data_file),
            "statistical_analyses": [{
                "id": "a1", "experiment": "e1", "name_col": "Synonym", "logfc_col": "ratio",
                "fold_change_type": fold_change_type, "adjusted_p_value_col": "padj",
                "pvalue_threshold": 0.05}]}},
    }}
    return config, str(tmp_path / "paperconfig.yaml")


def _fc_errors(errors):
    return [e for e in errors if "fold_change_type" in e]


def test_signed_linear_accepted(tmp_path):
    errors, _ = validate_paperconfig_content(*_config(tmp_path, [-3.0, 2.0, 5.0], "signed_linear"))
    assert _fc_errors(errors) == []


def test_signed_linear_value_in_unit_interval_is_error(tmp_path):
    # the offending value sits past the 5-row sample the validator reads for columns
    ratios = [-3.0, 2.0, 5.0, -4.0, 6.0, 2.5, 0.5]
    errors, _ = validate_paperconfig_content(*_config(tmp_path, ratios, "signed_linear"))
    assert any("lie in (-1, 1)" in e for e in _fc_errors(errors))


def test_unknown_fold_change_type_is_error(tmp_path):
    errors, _ = validate_paperconfig_content(*_config(tmp_path, [2.0], "ratio"))
    assert any("must be 'log2', 'linear' or 'signed_linear'" in e for e in errors)


def _analysis_errors(tmp_path, extra):
    config, path = _config(tmp_path, [-3.0, 2.0], "log2")
    a = config["publication"]["supplementary_materials"]["s1"]["statistical_analyses"][0]
    a.update(extra)
    errors, _ = validate_paperconfig_content(config, path)
    return errors


def test_logfc_threshold_required(tmp_path):
    assert any("'logfc_threshold' is required" in e for e in _analysis_errors(tmp_path, {}))
    assert not any("'logfc_threshold' is required" in e for e in _analysis_errors(tmp_path, {"logfc_threshold": 0}))
    assert not any("'logfc_threshold' is required" in e for e in _analysis_errors(tmp_path, {"prefiltered": True}))


def test_asterisk_min_level_validated(tmp_path):
    errs = _analysis_errors(tmp_path, {"logfc_threshold": 1.0, "asterisk_min_level": 2})
    assert any("needs 'pvalue_asterisk_in_logfc: true'" in e for e in errs)
    errs = _analysis_errors(tmp_path, {"pvalue_asterisk_in_logfc": True, "asterisk_min_level": 0})
    assert any("must be a positive integer" in e for e in errs)


def _suspect_warnings(warnings):
    return [w for w in warnings if "look like signed linear" in w]


def test_log2_table_that_looks_signed_linear_warns(tmp_path):
    # no value in (-1, 1) and max |x| 120: a signed linear table declared log2 (Kaur 2018 shape)
    ratios = [-3.0, 2.0, 5.0, -4.0, 6.0, 120.0]
    _, warnings = validate_paperconfig_content(*_config(tmp_path, ratios, None))
    assert len(_suspect_warnings(warnings)) == 1


def test_ordinary_log2_table_does_not_warn(tmp_path):
    for ratios in ([-3.0, 0.2, 25.0], [-3.0, 2.0, 14.0]):  # a near-zero value / small effects
        _, warnings = validate_paperconfig_content(*_config(tmp_path, ratios, "log2"))
        assert _suspect_warnings(warnings) == []
