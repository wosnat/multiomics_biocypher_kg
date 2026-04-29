"""Smoke tests for Phase 1.1A skeleton modules.

These verify the modules import and expose the expected API surface.
Real behaviour ships in Phase 1.1B; until then every stub raises
NotImplementedError.
"""
import pytest


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


def test_cazy_utils_imports_and_raises():
    from multiomics_kg.utils import cazy_utils

    assert hasattr(cazy_utils, "load_cazy")
    assert hasattr(cazy_utils, "cazy_ancestors")
    assert hasattr(cazy_utils, "is_valid_cazy")

    with pytest.raises(NotImplementedError):
        cazy_utils.load_cazy()
    with pytest.raises(NotImplementedError):
        cazy_utils.cazy_ancestors("GH13_1")
    with pytest.raises(NotImplementedError):
        cazy_utils.is_valid_cazy("GH13_1")
