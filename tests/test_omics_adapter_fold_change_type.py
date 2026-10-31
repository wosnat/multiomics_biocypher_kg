"""Tests for fold_change_type support in omics adapter."""
import math
import pytest


def test_linear_fc_converted_to_log2():
    """Linear fold-change 2.0 should become log2(2.0) = 1.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(2.0, "linear") == pytest.approx(1.0)


def test_linear_fc_less_than_one():
    """Linear fold-change 0.5 should become log2(0.5) = -1.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(0.5, "linear") == pytest.approx(-1.0)


def test_linear_fc_of_one():
    """Linear fold-change 1.0 (no change) should become 0.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(1.0, "linear") == pytest.approx(0.0)


def test_linear_fc_zero_returns_none():
    """Linear fold-change 0 is undefined in log2 -- return None."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(0.0, "linear") is None


def test_linear_fc_negative_returns_none():
    """Negative linear fold-change is invalid -- return None."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(-1.5, "linear") is None


def test_log2_fc_passthrough():
    """log2 fold-change should pass through unchanged."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(-2.5, "log2") == pytest.approx(-2.5)


def test_default_fc_type_is_log2():
    """When fold_change_type is None/omitted, treat as log2."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(1.5, None) == pytest.approx(1.5)


# --- _validate_fc_range tests ---

def test_validate_linear_with_negatives_warns(caplog):
    """Linear FC with negative values should warn."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([2.0, -0.5, 1.0], "linear", "test_analysis")
    assert "negative values" in caplog.text.lower()


def test_validate_log2_all_positive_above_one_warns(caplog):
    """log2 FC where all values > 1.0 should warn (likely linear)."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([2.5, 3.0, 1.5], "log2", "test_analysis")
    assert "may be linear" in caplog.text.lower()


def test_validate_log2_all_positive_above_one_no_warn_if_significant_only(caplog):
    """log2 FC all > 1.0 should NOT warn if table_scope is significant_only."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([2.5, 3.0, 1.5], "log2", "test_analysis", table_scope="significant_only")
    assert "may be linear" not in caplog.text.lower()


def test_validate_log2_mixed_sign_no_warn(caplog):
    """log2 FC with both positive and negative values should not warn."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([-1.5, 0.5, 2.0], "log2", "test_analysis")
    assert caplog.text == ""


def test_validate_linear_all_positive_no_warn(caplog):
    """Linear FC with all positive values should not warn."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([0.5, 2.0, 1.0], "linear", "test_analysis")
    assert caplog.text == ""


def test_validate_empty_values_no_warn(caplog):
    """Empty values list should not warn."""
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    import logging
    with caplog.at_level(logging.WARNING):
        _validate_fc_range([], "linear", "test_analysis")
    assert caplog.text == ""


def test_linear_downregulated_gene_can_be_significant(tmp_path):
    """Significance is judged on the log2 value, not the raw linear ratio.

    A ratio of 0.25 is log2 -2: past a |log2FC| >= 1 threshold, so significant.
    Before the fix the raw 0.25 was compared with the threshold and failed.
    """
    import pandas as pd
    import yaml
    from multiomics_kg.adapters.omics_adapter import OMICSAdapter

    data_file = tmp_path / "de.csv"
    pd.DataFrame({
        "Synonym": ["PMM0001", "PMM0002", "PMM0003"],
        "ratio": [4.0, 0.25, 1.2],
        "padj": [0.001, 0.001, 0.001],
    }).to_csv(data_file, index=False)
    config = {"publication": {
        "papername": "Linear FC", "doi": "10.1234/linear.fc", "control_timing": "matched",
        "experiments": {"e1": {
            "name": "T", "organism": "Prochlorococcus MED4", "treatment_condition": "T",
            "control_condition": "C", "omics_type": "RNASEQ", "test_type": "DESeq2",
            "treatment_type": ["nitrogen"], "background_factors": ["axenic"]}},
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(data_file),
            "statistical_analyses": [{
                "id": "a1", "experiment": "e1", "name_col": "Synonym", "logfc_col": "ratio",
                "fold_change_type": "linear", "adjusted_p_value_col": "padj",
                "pvalue_threshold": 0.05, "logfc_threshold": 1.0}]}},
    }}
    config_file = tmp_path / "paperconfig.yaml"
    config_file.write_text(yaml.dump(config))
    adapter = OMICSAdapter(config_file=str(config_file))
    adapter.extracted_data = {"title": "T", "authors": [], "journal": "", "year": 2024,
                              "abstract": "", "description": ""}
    edges = [e for e in adapter.get_edges() if "log2_fold_change" in e[4]]
    by_gene = {e[2].split(":")[-1]: e[4] for e in edges}
    assert by_gene["PMM0002"]["log2_fold_change"] == pytest.approx(-2.0)
    assert by_gene["PMM0002"]["significant"] == "significant"
    assert by_gene["PMM0001"]["significant"] == "significant"
    assert by_gene["PMM0003"]["significant"] == "not significant"



# --- signed_linear ("+3 = 3-fold up, -3 = 3-fold down") ---

@pytest.mark.parametrize("raw, expected", [
    (3.0, 1.584962500721156), (-3.0, -1.584962500721156),
    (2.0, 1.0), (-2.0, -1.0), (1.0, 0.0), (-1.0, 0.0),
])
def test_signed_linear_conversion(raw, expected):
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(raw, "signed_linear") == pytest.approx(expected)


def test_signed_linear_violations():
    from multiomics_kg.utils.fold_change import signed_linear_violations
    assert signed_linear_violations([3, -2, "1.0", "-1", "4.5*", "NA", None]) == []
    assert signed_linear_violations([3, 0.5, "-0.25 *", 0]) == [0.5, -0.25, 0.0]


def test_log2_heuristic_flags_signed_linear_shaped_table(caplog):
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    with caplog.at_level("WARNING"):
        _validate_fc_range([-8463.9, -3.0, 2.0, 353.9], "log2", "kaur_like")
    assert "signed_linear" in caplog.text


def test_log2_heuristic_quiet_on_normal_log2(caplog):
    from multiomics_kg.adapters.omics_adapter import _validate_fc_range
    with caplog.at_level("WARNING"):
        # no value in (-1, 1) but small magnitudes: a |log2FC| >= 1 filtered table
        _validate_fc_range([-6.6, -1.2, 1.1, 4.5], "log2", "filtered_log2")
        # large magnitudes but values in (-1, 1): ordinary log2
        _validate_fc_range([-25.0, 0.3, 2.0], "log2", "wide_log2")
    assert "signed_linear" not in caplog.text


def _run_adapter(tmp_path, ratios, analysis_extra):
    import pandas as pd
    import yaml
    from multiomics_kg.adapters.omics_adapter import OMICSAdapter

    data_file = tmp_path / "de.csv"
    pd.DataFrame({
        "Synonym": [f"PMM000{i + 1}" for i in range(len(ratios))],
        "ratio": ratios,
        "padj": [0.001] * len(ratios),
    }).to_csv(data_file, index=False)
    analysis = {"id": "a1", "experiment": "e1", "name_col": "Synonym", "logfc_col": "ratio",
                "adjusted_p_value_col": "padj", "pvalue_threshold": 0.05, "logfc_threshold": 1.0}
    analysis.update(analysis_extra)
    config = {"publication": {
        "papername": "Signed FC", "doi": "10.1234/signed.fc", "control_timing": "matched",
        "experiments": {"e1": {
            "name": "T", "organism": "Prochlorococcus MED4", "treatment_condition": "T",
            "control_condition": "C", "omics_type": "RNASEQ", "test_type": "DESeq2",
            "treatment_type": ["nitrogen"], "background_factors": ["axenic"]}},
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(data_file),
            "statistical_analyses": [analysis]}},
    }}
    config_file = tmp_path / "paperconfig.yaml"
    config_file.write_text(yaml.dump(config))
    adapter = OMICSAdapter(config_file=str(config_file))
    adapter.extracted_data = {"title": "T", "authors": [], "journal": "", "year": 2024,
                              "abstract": "", "description": ""}
    edges = [e for e in adapter.get_edges() if "log2_fold_change" in e[4]]
    return {e[2].split(":")[-1]: e[4] for e in edges}


def test_signed_linear_end_to_end(tmp_path):
    by_gene = _run_adapter(tmp_path, [-8.0, 3.0, 1.5],
                           {"fold_change_type": "signed_linear", "logfc_threshold": 1.0})
    assert by_gene["PMM0001"]["log2_fold_change"] == pytest.approx(-3.0)
    assert by_gene["PMM0001"]["expression_direction"] == "down"
    assert by_gene["PMM0002"]["expression_direction"] == "up"
    # 1.5-fold: log2 0.585, below the 2-fold threshold
    assert by_gene["PMM0003"]["significant"] == "not significant"


def test_signed_linear_with_unit_interval_value_refuses_analysis(tmp_path, caplog):
    with caplog.at_level("ERROR"):
        by_gene = _run_adapter(tmp_path, [-8.0, 3.0, 0.5],
                               {"fold_change_type": "signed_linear"})
    assert by_gene == {}
    assert "mis-declared" in caplog.text


def test_explicit_zero_logfc_threshold_is_honoured(tmp_path):
    """logfc_threshold: 0 means "q/p only", not the 1.0 default."""
    by_gene = _run_adapter(tmp_path, [-1.5, 1.2],
                           {"fold_change_type": "signed_linear", "logfc_threshold": 0})
    assert by_gene["PMM0001"]["significant"] == "significant"
    assert by_gene["PMM0002"]["significant"] == "significant"


# --- significance rules (2026-10-31 audit) ---

def _run_raw(tmp_path, rows, analysis_extra):
    """rows: list of (fc cell, p cell) strings; cells may be '' (empty)."""
    import pandas as pd
    import yaml
    from multiomics_kg.adapters.omics_adapter import OMICSAdapter

    data_file = tmp_path / "de.csv"
    pd.DataFrame({
        "Synonym": [f"PMM000{i + 1}" for i in range(len(rows))],
        "fc": [r[0] for r in rows],
        "padj": [r[1] for r in rows],
    }).to_csv(data_file, index=False)
    analysis = {"id": "a1", "experiment": "e1", "name_col": "Synonym", "logfc_col": "fc"}
    analysis.update(analysis_extra)
    config = {"publication": {
        "papername": "Sig rules", "doi": "10.1234/sig.rules", "control_timing": "matched",
        "experiments": {"e1": {
            "name": "T", "organism": "Prochlorococcus MED4", "treatment_condition": "T",
            "control_condition": "C", "omics_type": "RNASEQ", "test_type": "DESeq2",
            "treatment_type": ["nitrogen"], "background_factors": ["axenic"]}},
        "supplementary_materials": {"s1": {
            "type": "csv", "filename": str(data_file), "statistical_analyses": [analysis]}},
    }}
    config_file = tmp_path / "paperconfig.yaml"
    config_file.write_text(yaml.dump(config))
    adapter = OMICSAdapter(config_file=str(config_file))
    adapter.extracted_data = {"title": "T", "authors": [], "journal": "", "year": 2024,
                              "abstract": "", "description": ""}
    edges = [e for e in adapter.get_edges() if "log2_fold_change" in e[4]]
    return {e[2].split(":")[-1]: e[4] for e in edges}


def test_empty_p_value_cell_is_not_significant(tmp_path):
    """A declared p column with an empty cell (DESeq2 NA padj) is untested, not significant."""
    by_gene = _run_raw(tmp_path, [("3.0", "0.001"), ("3.0", "")],
                       {"adjusted_p_value_col": "padj", "pvalue_threshold": 0.05, "logfc_threshold": 0})
    assert by_gene["PMM0001"]["significant"] == "significant"
    assert by_gene["PMM0002"]["significant"] == "not significant"
    assert "adjusted_p_value" not in by_gene["PMM0002"]


def test_no_p_column_keeps_fold_only_rule(tmp_path):
    """A table with no p-value column at all is still judged on fold change alone."""
    by_gene = _run_raw(tmp_path, [("3.0", ""), ("0.5", "")], {"logfc_threshold": 1.0})
    assert by_gene["PMM0001"]["significant"] == "significant"
    assert by_gene["PMM0002"]["significant"] == "not significant"


def test_missing_logfc_threshold_refuses_analysis(tmp_path, caplog):
    with caplog.at_level("ERROR"):
        by_gene = _run_raw(tmp_path, [("3.0", "0.001")],
                           {"adjusted_p_value_col": "padj", "pvalue_threshold": 0.05})
    assert by_gene == {}
    assert "no logfc_threshold" in caplog.text


def test_prefiltered_and_asterisk_need_no_logfc_threshold(tmp_path):
    assert _run_raw(tmp_path, [("3.0", "")], {"prefiltered": True})["PMM0001"]["significant"] == "significant"
    by_gene = _run_raw(tmp_path, [("3.0*", ""), ("0.2", "")], {"pvalue_asterisk_in_logfc": True})
    assert by_gene["PMM0001"]["significant"] == "significant"
    assert by_gene["PMM0002"]["significant"] == "not significant"


def test_asterisk_min_level(tmp_path):
    rows = [("*2.0", ""), ("**2.0", ""), ("1.5 *", ""), ("1.5 **", ""), ("1.5", "")]
    default = _run_raw(tmp_path, rows, {"pvalue_asterisk_in_logfc": True})
    assert [default[f"PMM000{i}"]["significant"] for i in range(1, 6)] == \
        ["significant", "significant", "significant", "significant", "not significant"]
    strict = _run_raw(tmp_path, rows, {"pvalue_asterisk_in_logfc": True, "asterisk_min_level": 2})
    assert [strict[f"PMM000{i}"]["significant"] for i in range(1, 6)] == \
        ["not significant", "significant", "not significant", "significant", "not significant"]
    assert strict["PMM0001"]["log2_fold_change"] == pytest.approx(2.0)


def test_zero_fold_change_has_no_direction(tmp_path):
    by_gene = _run_raw(tmp_path, [("0", "0.5"), ("0.4", "0.5"), ("-0.4", "0.5")],
                       {"adjusted_p_value_col": "padj", "pvalue_threshold": 0.05, "logfc_threshold": 0})
    assert "expression_direction" not in by_gene["PMM0001"]
    assert by_gene["PMM0002"]["expression_direction"] == "up"
    assert by_gene["PMM0003"]["expression_direction"] == "down"


def test_signed_linear_suspect():
    from multiomics_kg.utils.fold_change import signed_linear_suspect
    assert signed_linear_suspect([-3, 2, "120*", "NA"]) == 120.0
    assert signed_linear_suspect([-3, 2, 14]) is None       # small effects: under the cutoff
    assert signed_linear_suspect([-3, 0.4, 120]) is None    # a value in (-1, 1): ordinary log2
    assert signed_linear_suspect([]) is None
