"""Smoke tests for Phase 1.1A skeleton modules.

These verify the modules import and expose the expected API surface.
Real behaviour ships in Phase 1.1B; until then every stub raises
NotImplementedError.
"""
import sqlite3
import pytest


def test_metabolite_utils_imports_and_raises():
    from multiomics_kg.utils import metabolite_utils

    # API surface is present
    assert hasattr(metabolite_utils, "open_resolver")
    assert hasattr(metabolite_utils, "resolve_metabolite")
    assert hasattr(metabolite_utils, "resolve_reaction")

    # Stubs raise NotImplementedError
    with pytest.raises(NotImplementedError):
        metabolite_utils.open_resolver()
    fake_conn = sqlite3.connect(":memory:")
    with pytest.raises(NotImplementedError):
        metabolite_utils.resolve_metabolite("C00031", fake_conn)
    with pytest.raises(NotImplementedError):
        metabolite_utils.resolve_reaction("R00200", fake_conn)


def test_tcdb_utils_imports_and_raises():
    from multiomics_kg.utils import tcdb_utils

    assert hasattr(tcdb_utils, "load_tcdb")
    assert hasattr(tcdb_utils, "tcdb_ancestors")
    assert hasattr(tcdb_utils, "is_valid_tcdb")

    with pytest.raises(NotImplementedError):
        tcdb_utils.load_tcdb()
    with pytest.raises(NotImplementedError):
        tcdb_utils.tcdb_ancestors("1.A.1.1.1")
    with pytest.raises(NotImplementedError):
        tcdb_utils.is_valid_tcdb("1.A.1.1.1")
