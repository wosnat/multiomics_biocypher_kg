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
