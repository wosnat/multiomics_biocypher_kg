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
