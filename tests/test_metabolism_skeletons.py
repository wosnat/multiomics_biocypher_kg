"""Smoke tests for Phase 1.1A skeleton modules.

These verify the modules import and expose the expected API surface.
Real behaviour ships in Phase 1.1B; until then every stub raises
NotImplementedError.
"""
import pytest


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
