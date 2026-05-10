# tests/test_tcdb_diamond.py
"""Unit tests for multiomics_kg.utils.tcdb_diamond."""
from multiomics_kg.utils import tcdb_diamond  # noqa: F401
from multiomics_kg.utils.tcdb_diamond import truncate_tcid


def test_truncate_tcid_keeps_first_n_parts():
    assert truncate_tcid("1.A.11.1.5", 3) == "1.A.11"
    assert truncate_tcid("1.A.11.1.5", 4) == "1.A.11.1"
    assert truncate_tcid("1.A.11.1.5", 5) == "1.A.11.1.5"


def test_truncate_tcid_passthrough_when_already_short():
    assert truncate_tcid("1.A.11", 5) == "1.A.11"
    assert truncate_tcid("1.A", 3) == "1.A"


def test_truncate_tcid_empty_input_returns_empty():
    assert truncate_tcid("", 3) == ""
