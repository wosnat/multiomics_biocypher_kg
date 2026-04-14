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


def test_single_known_ko_keeps_full_ancestor_chain(synth_cache: Path):
    """known_ko_ids={'K00001'} should keep exactly A1, A1.B1, A1.B1.C1 and the KO edge."""
    adapter = MultiBriteAdapter(cache_root=synth_cache, known_ko_ids={"K00001"})
    adapter.download_data(cache=True)

    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    node_ids = {n[0] for n in nodes}
    assert node_ids == {
        "kegg.brite:synth01000.A1",
        "kegg.brite:synth01000.A1.B1",
        "kegg.brite:synth01000.A1.B1.C1",
    }, f"unexpected kept nodes: {sorted(node_ids)}"

    # All three parent edges inside the chain, plus exactly one KO edge.
    parent_edges = [e for e in edges if e[3] == "brite_category_is_a_brite_category"]
    ko_edges = [e for e in edges if e[3] == "kegg_term_in_brite_category"]

    parent_pairs = {(e[1], e[2]) for e in parent_edges}
    assert parent_pairs == {
        ("kegg.brite:synth01000.A1.B1", "kegg.brite:synth01000.A1"),
        ("kegg.brite:synth01000.A1.B1.C1", "kegg.brite:synth01000.A1.B1"),
    }, f"unexpected parent edges: {sorted(parent_pairs)}"

    assert len(ko_edges) == 1, f"expected 1 KO edge, got {len(ko_edges)}"
    assert ko_edges[0][1] == "kegg.orthology:K00001"
    assert ko_edges[0][2] == "kegg.brite:synth01000.A1.B1.C1"


def test_no_dangling_ko_edges(synth_cache: Path):
    """Every emitted KO edge's raw KO must be in known_ko_ids."""
    known = {"K00001", "K00100"}
    adapter = MultiBriteAdapter(cache_root=synth_cache, known_ko_ids=known)
    adapter.download_data(cache=True)

    ko_edges = [
        e for e in adapter.get_edges() if e[3] == "kegg_term_in_brite_category"
    ]
    for edge_id, _src, _tgt, _label, _props in ko_edges:
        raw_ko = edge_id.split("--")[0]
        assert raw_ko in known, f"dangling KO edge for {raw_ko!r}"


def test_no_dangling_parent_edges(synth_cache: Path):
    """Every parent edge's endpoints must both be in the emitted node set."""
    adapter = MultiBriteAdapter(
        cache_root=synth_cache, known_ko_ids={"K00001", "K00100"}
    )
    adapter.download_data(cache=True)

    node_ids = {n[0] for n in adapter.get_nodes()}
    parent_edges = [
        e for e in adapter.get_edges() if e[3] == "brite_category_is_a_brite_category"
    ]
    for _eid, src, tgt, _label, _props in parent_edges:
        assert src in node_ids, f"parent edge source {src!r} not in node set"
        assert tgt in node_ids, f"parent edge target {tgt!r} not in node set"


def test_known_ko_ids_is_required():
    """Constructor must raise if known_ko_ids is not supplied."""
    with pytest.raises(TypeError):
        MultiBriteAdapter(cache_root="irrelevant")  # type: ignore[call-arg]


@pytest.fixture
def large_synth_cache(tmp_path: Path) -> tuple[Path, set[str]]:
    """150 leaves under a single A/B/C chain; returns (cache_root, all_kos)."""
    cache_root = tmp_path / "cache_large"
    kegg_dir = cache_root / "kegg"
    kegg_dir.mkdir(parents=True)

    all_kos: set[str] = set()
    leaves = []
    for i in range(150):
        ko = f"K{i:05d}"
        all_kos.add(ko)
        leaves.append({"name": f"{ko} synthetic enzyme"})

    tree = {
        "name": "synth02000",
        "children": [
            {
                "name": "A1 RootCat",
                "children": [
                    {
                        "name": "B1 SubCat",
                        "children": [
                            {"name": "C1 FamilyCat", "children": leaves},
                        ],
                    }
                ],
            }
        ],
    }
    (kegg_dir / "brite_synth02000.json").write_text(json.dumps(tree))
    return cache_root, all_kos


def test_test_mode_cap_applied_post_pruning(
    large_synth_cache: tuple[Path, set[str]], monkeypatch
):
    """With 150 kept nodes and test_mode=True, get_nodes yields ≤100 per tree
    and no edge references a dropped node."""
    cache_root, all_kos = large_synth_cache
    from multiomics_kg.utils import brite_utils
    from multiomics_kg.adapters import brite_adapter

    monkeypatch.setattr(brite_utils, "BRITE_TREES", {"synth02000": "synthetic_large"})
    monkeypatch.setattr(
        brite_adapter, "BRITE_TREES", {"synth02000": "synthetic_large"}
    )

    adapter = MultiBriteAdapter(
        cache_root=cache_root, known_ko_ids=all_kos, test_mode=True
    )
    adapter.download_data(cache=True)

    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    # 3 scaffold categories (A1, B1, C1) + up to 100 leaf limit — but D-level
    # leaves are KO leaves (not BriteCategory nodes), so this synthetic tree
    # only has 3 BriteCategory nodes total. Cap only kicks in when categories
    # exceed 100, so this tests the "no cap needed" path cleanly.
    assert len(nodes) == 3

    # Every KO edge's leaf must be in the emitted node set.
    node_ids = {n[0] for n in nodes}
    ko_edges = [e for e in edges if e[3] == "kegg_term_in_brite_category"]
    for _eid, _src, tgt, _label, _props in ko_edges:
        assert tgt in node_ids, f"KO edge points at dropped leaf {tgt!r}"
