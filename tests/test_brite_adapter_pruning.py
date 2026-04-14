"""Unit tests for MultiBriteAdapter subhierarchy pruning.

Uses synthetic in-memory tree JSON to exercise the pruning algorithm
without hitting the 12-tree cache on disk.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from multiomics_kg.adapters.brite_adapter import MultiBriteAdapter


# Minimal synthetic BRITE tree with 3 A-level branches, mixed KO leaves.
# Structure:
#   A1 Metabolism
#     B1 Glycolysis
#       C1 oxidoreductase
#         D1 K00001 alcohol dehydrogenase   (leaf KO)
#         D2 K00002 other dehydrogenase     (leaf KO)
#     B2 TCA
#       C1 K00004 tricarboxylic enzyme      (leaf KO)
#   A2 Genetic Information
#     B1 Translation
#       C1 K00100 ribosomal protein          (leaf KO)
#   A3 EmptyBranch
#     B1 NothingUseful
#       C1 K99999 unreferenced              (leaf KO)
SYNTHETIC_TREE = {
    "name": "synth01000",
    "children": [
        {
            "name": "A1 Metabolism",
            "children": [
                {
                    "name": "B1 Glycolysis",
                    "children": [
                        {
                            "name": "C1 oxidoreductase",
                            "children": [
                                {"name": "K00001 alcohol dehydrogenase"},
                                {"name": "K00002 other dehydrogenase"},
                            ],
                        }
                    ],
                },
                {
                    "name": "B2 TCA",
                    "children": [{"name": "K00004 tricarboxylic enzyme"}],
                },
            ],
        },
        {
            "name": "A2 Genetic Information",
            "children": [
                {
                    "name": "B1 Translation",
                    "children": [{"name": "K00100 ribosomal protein"}],
                }
            ],
        },
        {
            "name": "A3 EmptyBranch",
            "children": [
                {
                    "name": "B1 NothingUseful",
                    "children": [{"name": "K99999 unreferenced"}],
                }
            ],
        },
    ],
}


@pytest.fixture
def synth_cache(tmp_path: Path) -> Path:
    """Write synthetic tree to a fake cache dir and return the cache root."""
    cache_root = tmp_path / "cache"
    kegg_dir = cache_root / "kegg"
    kegg_dir.mkdir(parents=True)
    (kegg_dir / "brite_synth01000.json").write_text(json.dumps(SYNTHETIC_TREE))
    return cache_root


@pytest.fixture(autouse=True)
def patch_brite_trees(monkeypatch):
    """Override the BRITE_TREES registry to point at our synthetic tree only."""
    from multiomics_kg.utils import brite_utils

    monkeypatch.setattr(
        brite_utils,
        "BRITE_TREES",
        {"synth01000": "synthetic"},
    )
    # Also patch the re-exported reference in brite_adapter.
    from multiomics_kg.adapters import brite_adapter

    monkeypatch.setattr(
        brite_adapter,
        "BRITE_TREES",
        {"synth01000": "synthetic"},
    )


def test_empty_known_kos_prunes_everything(synth_cache: Path, caplog):
    """With known_ko_ids = ∅, no nodes or edges should be emitted."""
    caplog.set_level(logging.INFO)
    adapter = MultiBriteAdapter(cache_root=synth_cache, known_ko_ids=set())
    adapter.download_data(cache=True)

    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    assert nodes == [], f"Expected 0 nodes, got {len(nodes)}"
    assert edges == [], f"Expected 0 edges, got {len(edges)}"
    assert any(
        "pruned entirely" in rec.message for rec in caplog.records
    ), "expected 'pruned entirely' log message for empty tree"
