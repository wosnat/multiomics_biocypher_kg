"""Tests for the Phase 1.1B metabolism transforms registry.

Note: resolve_kegg_reaction_to_mnxr was removed in Spec 1.2 pivot — raw KEGG
R-numbers are now kept as-is (no MNX resolution at build time).

Note: validate_cazy was removed when CAZy hierarchy was promoted to a
pure-Python in-process table; CAZy IDs are passthrough at this layer and
filtered at the cazy_adapter level.

Note: validate_tcdb was removed as part of the TCDB ontology promotion (commit
2 of the TCDB/CAZy ontologies feature); TCDB IDs are passthrough at this layer
and unknown IDs are filtered at adapter time.
"""
from __future__ import annotations

from multiomics_kg.download.utils import annotation_transforms as at


def test_resolve_kegg_reaction_transform_removed():
    """Spec 1.2 pivot: KEGG reactions stay as raw R-numbers (no MNX resolution)."""
    assert "resolve_kegg_reaction_to_mnxr" not in at._TRANSFORMS


def test_validate_tcdb_transform_removed():
    """TCDB hierarchy promoted to first-class ontology; transform no longer exists."""
    assert "validate_tcdb" not in at._TRANSFORMS
    assert not hasattr(at, "_tx_validate_tcdb")


def test_validate_cazy_transform_removed():
    """CAZy hierarchy moved to pure-Python utils; transform no longer exists."""
    assert "validate_cazy" not in at._TRANSFORMS
