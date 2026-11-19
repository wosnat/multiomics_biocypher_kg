import math

import pytest

from multiomics_kg.utils.de_calls import (
    RowCall, call_row, check_significance, expression_status,
)

AST = {"logfc_col": "fc", "pvalue_asterisk_in_logfc": True,
       "asterisk_min_level": 2, "pvalue_threshold": 0.05}


def test_call_row_asterisk_levels():
    r = call_row({"fc": "**1.5"}, AST, {"fc"})
    assert r.significant == "significant" and r.adjusted_p_value == 0.05
    r = call_row({"fc": "1.5 *"}, AST, {"fc"})
    assert r.significant == "not significant" and r.adjusted_p_value == 1.0


def test_call_row_empty_pvalue_is_not_significant():
    a = {"logfc_col": "fc", "adjusted_p_value_col": "padj", "logfc_threshold": 1}
    r = call_row({"fc": 3.0, "padj": float("nan")}, a, {"fc", "padj"})
    assert r.significant == "not significant"


@pytest.mark.parametrize("fc", ["NA", "", "abc", float("inf"), float("nan")])
def test_call_row_unparseable_fc_returns_none(fc):
    assert call_row({"fc": fc}, {"logfc_col": "fc"}, {"fc"}) is None


def test_call_row_linear_conversion():
    a = {"logfc_col": "fc", "fold_change_type": "linear", "logfc_threshold": 0}
    assert call_row({"fc": 0.25}, a, {"fc"}).log2_fold_change == pytest.approx(-2.0)
    r = call_row({"fc": 0}, a, {"fc"})
    assert isinstance(r, RowCall) and r.raw_fc == 0.0 and r.log2_fold_change is None


def test_check_significance_explicit_zero_threshold():
    a = {"logfc_threshold": 0, "pvalue_threshold": 0.05}
    assert check_significance(0.01, 0.01, None, a) == "significant"


def test_expression_status_zero_fc():
    assert expression_status("significant", 0.0) == "not_significant"
    assert expression_status("significant", 1.0) == "significant_up"
    assert expression_status("significant", -1.0) == "significant_down"


# ---------------------------------------------------------------------------
# Low-signal rule + per-time-point statistics
# ---------------------------------------------------------------------------
import pandas as pd

from multiomics_kg.utils.de_calls import (
    LOW_SIGNAL_CONDITION, LOW_SIGNAL_MEAN, level_basis, low_signal_mask,
    low_signal_reason, time_point_stats,
)

ARMS = {"treatment_level_cols": ["t"], "control_level_cols": ["c"]}


def test_level_basis_kinds_and_errors():
    assert level_basis(ARMS) == LOW_SIGNAL_CONDITION
    assert level_basis({"mean_level_cols": ["m"]}) == LOW_SIGNAL_MEAN
    assert level_basis({}) is None
    with pytest.raises(ValueError):
        level_basis({"treatment_level_cols": ["t"]})
    with pytest.raises(ValueError):
        level_basis({**ARMS, "mean_level_cols": ["m"]})
    assert level_basis({"mean_level_cols": "baseMean"}) == LOW_SIGNAL_MEAN


def test_mask_rule1_minimum_and_empty():
    df = pd.DataFrame({"t": [0.01, 5, 7, None], "c": [3, 0.01, 4, 2]})
    assert low_signal_mask(df, ARMS) == [True, True, False, True]


def test_mask_rule1_replicates_need_all_at_minimum():
    df = pd.DataFrame({"r1": [0, 0, 2], "r2": [0, 3, 0], "c": [5, 6, 7]})
    a = {"treatment_level_cols": ["r1", "r2"], "control_level_cols": ["c"]}
    assert low_signal_mask(df, a) == [True, False, False]


def test_mask_constant_column_is_all_minimum():
    df = pd.DataFrame({"t": [1, 2, 3], "c": [5, 5, 5]})
    assert low_signal_mask(df, ARMS) == [True, True, True]


def test_mask_non_numeric_text_is_empty():
    df = pd.DataFrame({"t": ["n.d.", 4, 5], "c": [3, 4, 5]})
    assert low_signal_mask(df, ARMS)[0] is True


def test_mask_all_empty_column():
    df = pd.DataFrame({"t": [None, None], "c": [1, 2]})
    assert low_signal_mask(df, ARMS) == [True, True]


def test_mask_rule2_bottom_decile():
    df = pd.DataFrame({"m": list(range(1, 21))})
    mask = low_signal_mask(df, {"mean_level_cols": ["m"]})
    assert mask == [v <= 2 for v in range(1, 21)]
    assert isinstance(mask, list) and all(isinstance(x, bool) for x in mask)


def test_mask_none_without_level_cols():
    assert low_signal_mask(pd.DataFrame({"a": [1]}), {}) is None


def test_low_signal_reason_rule2_needs_large_fc():
    assert low_signal_reason(LOW_SIGNAL_MEAN, True, 5.0) is None
    assert low_signal_reason(LOW_SIGNAL_MEAN, True, -5.01) == LOW_SIGNAL_MEAN
    assert low_signal_reason(LOW_SIGNAL_CONDITION, True, 0.1) == LOW_SIGNAL_CONDITION
    assert low_signal_reason(LOW_SIGNAL_CONDITION, False, 9) is None
    assert low_signal_reason(None, True, 9) is None


def test_time_point_stats_values():
    s = time_point_stats([
        (-1, "not_significant", None),
        (2, "significant_up", None),
        (-6, "not_significant", "condition_at_column_minimum"),
        (11, "significant_up", None),
        (0, "not_significant", None),
    ])
    assert s.total == 5
    assert s.negative_fraction == pytest.approx(0.4)
    assert s.median_log2fc == 0.0
    assert s.max_abs_log2fc == 11.0
    assert s.large_fc_not_significant == 1
    assert s.extreme_fc_count == 1
    assert s.low_signal_count == 1


def test_time_point_stats_all_zero():
    s = time_point_stats([(0.0, "not_significant", None)] * 3)
    assert s.negative_fraction == 0.0 and s.median_log2fc == 0.0


def test_time_point_stats_median_interpolates_and_empty():
    s = time_point_stats([(1, "x", None), (2, "x", None)])
    assert s.median_log2fc == 1.5
    e = time_point_stats([])
    assert e.total == 0 and e.median_log2fc == 0.0 and e.low_signal_count == 0
