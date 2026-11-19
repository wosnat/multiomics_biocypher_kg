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
