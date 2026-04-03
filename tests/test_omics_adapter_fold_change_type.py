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
