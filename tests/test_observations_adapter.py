"""Tests for ObservationsAdapter (Plan 2)."""
import os
from pathlib import Path

import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.observations_adapter import (
    ObservationsAdapter,
    MultiObservationsAdapter,
    _clean_str,
    _make_derived_metric_id,
    _resolve_csv_path,
    _parse_boolean_cell,
)


def test_clean_str_replaces_single_quote_and_pipe():
    assert _clean_str("it's | great") == "it^s , great"


def test_clean_str_handles_none():
    assert _clean_str(None) == ""


def test_clean_str_passes_non_string_through():
    # ints stringify but don't go through replacement
    assert _clean_str(42) == "42"


def test_make_derived_metric_id_uses_doi_short():
    dm_id = _make_derived_metric_id(
        doi="10.1128/mSystems.00040-18",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:mSystems.00040-18:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_make_derived_metric_id_falls_back_to_paper_slug():
    dm_id = _make_derived_metric_id(
        doi="",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:biller_2018:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_resolve_csv_path_prefers_resolved(tmp_path):
    src = tmp_path / "table.csv"
    resolved = tmp_path / "table_resolved.csv"
    src.write_text("a,b\n1,2\n")
    resolved.write_text("a,b,resolved_locus_tag\n1,2,PMM0001\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == resolved
    assert used_resolved is True


def test_resolve_csv_path_falls_back_to_original(tmp_path):
    src = tmp_path / "table.csv"
    src.write_text("a,b\n1,2\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == src
    assert used_resolved is False


def _bp(value, blank_policy="skip"):
    return _parse_boolean_cell(
        value,
        true_tokens=["Y", "yes"],
        false_tokens=["N", "no"],
        skip_tokens=["NA", "n/a"],
        blank_policy=blank_policy,
    )


def test_boolean_true_token():
    assert _bp("Y") == "true"
    assert _bp("yes") == "true"


def test_boolean_false_token():
    assert _bp("N") == "false"
    assert _bp("no") == "false"


def test_boolean_skip_token_returns_none():
    assert _bp("NA") is None
    assert _bp("n/a") is None


def test_boolean_blank_with_skip_policy():
    assert _bp("", blank_policy="skip") is None
    assert _bp(None, blank_policy="skip") is None
    assert _bp(float("nan"), blank_policy="skip") is None


def test_boolean_blank_with_true_policy():
    assert _bp("", blank_policy="true") == "true"
    assert _bp(None, blank_policy="true") == "true"


def test_boolean_blank_with_false_policy():
    assert _bp("", blank_policy="false") == "false"
    assert _bp(float("nan"), blank_policy="false") == "false"


def test_boolean_unknown_token_raises():
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp("maybe")


def test_boolean_invalid_blank_policy_raises():
    with pytest.raises(ValueError, match="Invalid blank_policy"):
        _parse_boolean_cell("", true_tokens=["Y"], false_tokens=[], skip_tokens=[], blank_policy="whatever")


def test_boolean_whitespace_treated_as_blank():
    # Empty-after-strip → blank → skip_policy default returns None
    assert _bp("   ") is None
    # With blank_policy="false", whitespace is blank → "false"
    assert _bp("   ", blank_policy="false") == "false"
