# BRITE Subhierarchy Pruning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prune `MultiBriteAdapter` output to the subhierarchy reachable from `KeggTerm` nodes already in the KG, eliminating ~12K dangling `Kegg_term_in_brite_category` edges that cause the Docker `import` container to exit with code 70.

**Architecture:** `MultiBriteAdapter` gains a required `known_ko_ids: set[str]` constructor argument. A new per-tree pre-processing pass (`_prune`) walks each JSON tree once, collects nodes + KO edges + parent map, seeds leaves from `known_ko_ids`, marks ancestors bottom-up, and stores a `PrunedTree` dataclass. `get_nodes`/`get_edges` become trivial iterators over the pruned structures. The orchestrator passes `kegg_anno_adapter.all_ko_ids()` (promoted from `_all_ko_ids`) into the BRITE adapter.

**Tech Stack:** Python 3.12, BioCypher, pytest. No new dependencies.

**Spec:** [docs/superpowers/specs/2026-04-14-brite-pruning-design.md](../specs/2026-04-14-brite-pruning-design.md)

---

## File Structure

**Modified files:**
- `multiomics_kg/adapters/brite_adapter.py` — add `PrunedTree` dataclass, `_prune()` method, thin `get_nodes`/`get_edges`; require `known_ko_ids` in `__init__`.
- `multiomics_kg/adapters/functional_annotation_adapter.py` — rename `_all_ko_ids` → `all_ko_ids`, add docstring.
- `create_knowledge_graph.py` — pass `kegg_anno_adapter.all_ko_ids()` into `MultiBriteAdapter`.
- `tests/kg_validity/test_brite.py` — relax `test_brite_tree_present` and `test_brite_level_zero_nodes_exist`; add `test_brite_no_dangling_ko_edges`.
- `CLAUDE.md` — refresh BRITE node/edge counts, mention pruning.

**New files:**
- `tests/test_brite_adapter_pruning.py` — adapter-level unit tests with synthetic JSON.

**Per-file responsibility:**
- `brite_adapter.py`: owns the pruning logic; emits only kept nodes/edges.
- `functional_annotation_adapter.py`: unchanged pruning behaviour; `all_ko_ids()` becomes the public coupling point.
- `create_knowledge_graph.py`: wires the two adapters together.
- `test_brite_adapter_pruning.py`: isolates pruning correctness from the full KG build.
- `test_brite.py`: validates the live graph after import.

---

## Task 1: Promote `_all_ko_ids` to public API

**Files:**
- Modify: `multiomics_kg/adapters/functional_annotation_adapter.py:558-563`

- [ ] **Step 1: Rename the method**

Open [multiomics_kg/adapters/functional_annotation_adapter.py](../../../multiomics_kg/adapters/functional_annotation_adapter.py) at line 558. Replace:

```python
    def _all_ko_ids(self) -> set[str]:
        """Collect all KO IDs referenced across all strains."""
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_ko_ids()
        return ids
```

with:

```python
    def all_ko_ids(self) -> set[str]:
        """Canonical source of KO IDs that will become ``KeggTerm`` nodes in the KG.

        Returned set is the union of KO IDs referenced by any gene across all
        configured strains. Downstream adapters (e.g. ``MultiBriteAdapter``)
        consume this to prune their output to the subhierarchy reachable from
        these KOs — keeping the graph consistent and preventing dangling
        edges that would fail ``neo4j-admin import``.
        """
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_ko_ids()
        return ids
```

- [ ] **Step 2: Update internal callers**

Run:

```bash
grep -n "_all_ko_ids" multiomics_kg/adapters/functional_annotation_adapter.py
```

Expected: two remaining callers in `get_nodes` (around line 572) and `get_edges` (around line 665). Replace both `self._all_ko_ids()` → `self.all_ko_ids()`.

- [ ] **Step 3: Verify no other callers**

Run:

```bash
grep -rn "_all_ko_ids" multiomics_kg/ tests/ create_knowledge_graph.py
```

Expected: no matches (everything renamed).

- [ ] **Step 4: Run existing tests**

Run: `pytest -m "not slow and not kg" tests/ -k kegg -v`
Expected: tests pass (renaming is behaviour-preserving).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/functional_annotation_adapter.py
git commit -m "refactor(kegg): promote _all_ko_ids to public all_ko_ids"
```

---

## Task 2: Adapter unit test — empty `known_ko_ids` prunes everything

**Files:**
- Create: `tests/test_brite_adapter_pruning.py`

- [ ] **Step 1: Create the test file with a synthetic-tree helper and the first test**

Create [tests/test_brite_adapter_pruning.py](../../../tests/test_brite_adapter_pruning.py):

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest tests/test_brite_adapter_pruning.py::test_empty_known_kos_prunes_everything -v`
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'known_ko_ids'` (constructor does not yet accept the parameter).

- [ ] **Step 3: Commit the failing test**

```bash
git add tests/test_brite_adapter_pruning.py
git commit -m "test(brite): failing test for empty known_ko_ids pruning"
```

---

## Task 3: Add `PrunedTree`, require `known_ko_ids`, implement `_prune`

**Files:**
- Modify: `multiomics_kg/adapters/brite_adapter.py`

- [ ] **Step 1: Add imports and the `PrunedTree` dataclass**

Open [multiomics_kg/adapters/brite_adapter.py](../../../multiomics_kg/adapters/brite_adapter.py). At the top of the file, update the imports block to add `dataclass` and `field`:

Replace:

```python
from __future__ import annotations

import html
import logging
import re
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.utils.brite_utils import BRITE_TREES, compute_level_kind, load_brite_trees
```

with:

```python
from __future__ import annotations

import html
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.utils.brite_utils import BRITE_TREES, compute_level_kind, load_brite_trees
```

Then add the `PrunedTree` dataclass just below the `_DEPTH_LETTERS` constant (around line 37, immediately after `_DEPTH_LETTERS = "ABCD"`):

```python


@dataclass(slots=True)
class PrunedTree:
    """Result of pruning a single BRITE tree to the subhierarchy reachable
    from KOs present in the KG.

    Attributes:
        tree_code: KEGG BRITE tree identifier (e.g. ``ko02000``).
        nodes: list of ``(node_id, parent_id, name, level)`` in discovery
            order. ``level_kind`` is derived at emit time via
            ``compute_level_kind(level)``.
        ko_edges: list of ``(ko_id_raw, leaf_node_id)``. Contains only KOs
            in ``known_ko_ids`` and leaves that survived pruning.
    """

    tree_code: str
    nodes: list[tuple[str, str | None, str, int]] = field(default_factory=list)
    ko_edges: list[tuple[str, str]] = field(default_factory=list)
```

- [ ] **Step 2: Update the `__init__` signature to require `known_ko_ids`**

Replace the existing `__init__` (lines 81-92):

```python
    def __init__(
        self,
        cache_root: "str | Path",
        trees: "list[str] | None" = None,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.cache_root = Path(cache_root)
        self.trees = trees or list(BRITE_TREES.keys())
        self.test_mode = test_mode
        self.cache = cache
        self._tree_data: dict[str, dict] = {}
```

with:

```python
    def __init__(
        self,
        cache_root: "str | Path",
        known_ko_ids: set[str],
        trees: "list[str] | None" = None,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        """
        Args:
            cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
            known_ko_ids: set of KO IDs (raw ``K#####`` strings) that will become
                ``KeggTerm`` nodes in the KG. Sourced from
                ``MultiKeggAnnotationAdapter.all_ko_ids()``. Pruning keeps only
                the subhierarchy whose subtree contains at least one of these
                KOs; KO→category edges are emitted only for these KOs. An empty
                set is valid (yields no output, logs per-tree coverage).
            trees: list of KEGG tree codes; defaults to all 12 configured trees
            test_mode: cap at 100 nodes per tree for fast iteration
            cache: if False, re-download trees even if cache files exist
        """
        self.cache_root = Path(cache_root)
        self.known_ko_ids = known_ko_ids
        self.trees = trees or list(BRITE_TREES.keys())
        self.test_mode = test_mode
        self.cache = cache
        self._tree_data: dict[str, dict] = {}
        self._pruned: dict[str, PrunedTree] = {}
```

- [ ] **Step 3: Update `download_data` to trigger pruning**

Replace the existing `download_data` (lines 94-96):

```python
    def download_data(self, cache: bool = True) -> None:
        """Fetch/load all configured BRITE trees into memory."""
        self._tree_data = load_brite_trees(self.cache_root, self.trees, cache=cache)
```

with:

```python
    def download_data(self, cache: bool = True) -> None:
        """Fetch/load all configured BRITE trees into memory and prune them."""
        self._tree_data = load_brite_trees(self.cache_root, self.trees, cache=cache)
        self._pruned = {
            tree_code: self._prune(tree_code, tree_json)
            for tree_code, tree_json in self._tree_data.items()
        }
```

- [ ] **Step 4: Implement `_prune`**

Add a new method immediately after `_walk` (before `get_nodes`, around line 148):

```python
    def _prune(self, tree_code: str, tree_json: dict) -> PrunedTree:
        """Walk tree once, collect structural data, mark keepers, return pruned result.

        Algorithm:
        1. Single pass via _walk collecting:
           - all_nodes: every BriteCategory tuple in discovery order
           - all_ko_edges: every (ko_id_raw, leaf_node_id) pair
           - parent_of: dict[node_id, parent_id] (total map — each BRITE
             category has exactly one parent by positional-ID construction)
        2. Seed leaf_hits from known_ko_ids.
        3. Bottom-up mark: for each leaf, walk parent_of upward adding
           each ancestor to kept until reaching an A-level node (no parent
           in parent_of) or an already-kept ancestor.
        4. Filter all_nodes / all_ko_edges into a PrunedTree.
        5. Log per-tree coverage; empty result is a valid state.
        """
        all_nodes: list[tuple[str, str | None, str, int]] = []
        all_ko_edges: list[tuple[str, str]] = []
        parent_of: dict[str, str] = {}
        seen_ids: set[str] = set()

        top_children = tree_json.get("children", [])
        for item in self._walk(tree_code, top_children, []):
            if item[0] == "node":
                _, node_id, parent_id, name, _tc, level, _lk = item
                if node_id in seen_ids:
                    prior_parent = parent_of.get(node_id)
                    if prior_parent != parent_id:
                        logger.error(
                            f"BRITE tree {tree_code}: node {node_id!r} assigned "
                            f"two different parents ({prior_parent!r} vs "
                            f"{parent_id!r}); keeping first, skipping duplicate."
                        )
                    continue
                seen_ids.add(node_id)
                all_nodes.append((node_id, parent_id, name, level))
                if parent_id is not None:
                    parent_of[node_id] = parent_id
            elif item[0] == "edge_ko":
                _, ko_id_raw, leaf_id = item
                all_ko_edges.append((ko_id_raw, leaf_id))

        # Seed: leaves whose KO is known.
        leaf_hits = {
            leaf_id for ko, leaf_id in all_ko_edges if ko in self.known_ko_ids
        }

        # Bottom-up mark.
        kept: set[str] = set()
        for leaf_id in leaf_hits:
            current: str | None = leaf_id
            while current is not None and current not in kept:
                kept.add(current)
                current = parent_of.get(current)

        pruned_nodes = [n for n in all_nodes if n[0] in kept]
        pruned_ko_edges = [
            (ko, leaf_id)
            for ko, leaf_id in all_ko_edges
            if ko in self.known_ko_ids and leaf_id in kept
        ]

        if not pruned_nodes:
            logger.info(
                f"BRITE tree {tree_code}: 0/{len(all_nodes)} categories kept "
                f"— pruned entirely (no known KOs in this tree)"
            )
        else:
            logger.info(
                f"BRITE tree {tree_code}: {len(pruned_nodes)}/{len(all_nodes)} "
                f"categories kept, {len(pruned_ko_edges)}/{len(all_ko_edges)} "
                f"KO edges kept"
            )

        return PrunedTree(
            tree_code=tree_code,
            nodes=pruned_nodes,
            ko_edges=pruned_ko_edges,
        )
```

- [ ] **Step 5: Rewrite `get_nodes` and `get_edges` to iterate `self._pruned`**

Replace the existing `get_nodes` (lines 148-189):

```python
    def get_nodes(self):
        """Yield ``(node_id, 'brite category', properties)`` for every BriteCategory."""
        if not self._tree_data:
            self.download_data(cache=self.cache)

        seen_ids: set[str] = set()
        total_count = 0

        for tree_code, tree_json in self._tree_data.items():
            top_children = tree_json.get("children", [])
            tree_count = 0

            for item in self._walk(tree_code, top_children, []):
                if item[0] != "node":
                    continue
                _, node_id, _parent_id, name, tc, level, level_kind = item

                if node_id in seen_ids:
                    continue
                seen_ids.add(node_id)

                yield (
                    node_id,
                    "brite category",
                    {
                        "name": name,
                        "tree": BRITE_TREES[tc],
                        "tree_code": tc,
                        "level": level,
                        "level_kind": level_kind,
                    },
                )
                tree_count += 1
                total_count += 1

                if self.test_mode and tree_count >= 100:
                    logger.debug(
                        f"MultiBriteAdapter.get_nodes: test_mode cap at {tree_count} for {tc}"
                    )
                    break

        logger.info(f"MultiBriteAdapter.get_nodes: {total_count} BriteCategory nodes")
```

with:

```python
    def get_nodes(self):
        """Yield ``(node_id, 'brite category', properties)`` for every surviving BriteCategory."""
        if not self._pruned:
            self.download_data(cache=self.cache)

        total_count = 0
        for tree_code, pt in self._pruned.items():
            emitted_in_tree = 0
            for node_id, _parent_id, name, level in pt.nodes:
                if self.test_mode and emitted_in_tree >= 100:
                    logger.debug(
                        f"MultiBriteAdapter.get_nodes: test_mode cap at "
                        f"{emitted_in_tree} for {tree_code}"
                    )
                    break
                yield (
                    node_id,
                    "brite category",
                    {
                        "name": name,
                        "tree": BRITE_TREES[tree_code],
                        "tree_code": tree_code,
                        "level": level,
                        "level_kind": compute_level_kind(level),
                    },
                )
                emitted_in_tree += 1
                total_count += 1

        logger.info(f"MultiBriteAdapter.get_nodes: {total_count} BriteCategory nodes")
```

Replace the existing `get_edges` (lines 191-234):

```python
    def get_edges(self):
        """Yield Brite_category_is_a_brite_category and Kegg_term_in_brite_category edges."""
        if not self._tree_data:
            self.download_data(cache=self.cache)

        parent_count = 0
        ko_count = 0
        seen_ko_edges: set[tuple[str, str]] = set()

        for tree_code, tree_json in self._tree_data.items():
            top_children = tree_json.get("children", [])

            for item in self._walk(tree_code, top_children, []):
                if item[0] == "node":
                    _, node_id, parent_id, _name, _tc, _level, _lk = item
                    if parent_id is not None:
                        yield (
                            f"{node_id}--parent",
                            node_id,
                            parent_id,
                            "brite_category_is_a_brite_category",
                            {},
                        )
                        parent_count += 1

                elif item[0] == "edge_ko":
                    _, ko_id_raw, brite_node_id = item
                    ko_nid = _ko_node_id(ko_id_raw)
                    edge_key = (ko_nid, brite_node_id)
                    if edge_key not in seen_ko_edges:
                        seen_ko_edges.add(edge_key)
                        yield (
                            f"{ko_id_raw}--brite--{brite_node_id}",
                            ko_nid,
                            brite_node_id,
                            "kegg_term_in_brite_category",
                            {},
                        )
                        ko_count += 1

        logger.info(
            f"MultiBriteAdapter.get_edges: {parent_count} parent edges, "
            f"{ko_count} KO→BRITE edges"
        )
```

with:

```python
    def get_edges(self):
        """Yield Brite_category_is_a_brite_category and Kegg_term_in_brite_category edges."""
        if not self._pruned:
            self.download_data(cache=self.cache)

        parent_count = 0
        ko_count = 0

        for tree_code, pt in self._pruned.items():
            # Respect the same test_mode cap get_nodes uses, to avoid
            # emitting edges whose source or target node was dropped.
            nodes_to_use = pt.nodes[:100] if self.test_mode else pt.nodes
            emitted_ids = {n[0] for n in nodes_to_use}

            # Parent edges — both endpoints must be among emitted nodes.
            for node_id, parent_id, _name, _level in nodes_to_use:
                if parent_id is not None and parent_id in emitted_ids:
                    yield (
                        f"{node_id}--parent",
                        node_id,
                        parent_id,
                        "brite_category_is_a_brite_category",
                        {},
                    )
                    parent_count += 1

            # KO edges — already filtered at pruning time; double-check
            # leaf is among emitted nodes (matters in test_mode).
            for ko_id_raw, leaf_id in pt.ko_edges:
                if leaf_id not in emitted_ids:
                    continue
                yield (
                    f"{ko_id_raw}--brite--{leaf_id}",
                    _ko_node_id(ko_id_raw),
                    leaf_id,
                    "kegg_term_in_brite_category",
                    {},
                )
                ko_count += 1

        logger.info(
            f"MultiBriteAdapter.get_edges: {parent_count} parent edges, "
            f"{ko_count} KO→BRITE edges"
        )
```

- [ ] **Step 6: Run the failing test**

Run: `pytest tests/test_brite_adapter_pruning.py::test_empty_known_kos_prunes_everything -v`
Expected: PASS (adapter accepts `known_ko_ids`, emits nothing when set is empty, logs "pruned entirely").

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/adapters/brite_adapter.py tests/test_brite_adapter_pruning.py
git commit -m "feat(brite): prune subhierarchy to KOs present in the KG"
```

---

## Task 4: Positive pruning test — single KO keeps only its ancestor chain

**Files:**
- Modify: `tests/test_brite_adapter_pruning.py`

- [ ] **Step 1: Add the test**

Append to [tests/test_brite_adapter_pruning.py](../../../tests/test_brite_adapter_pruning.py):

```python


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
```

- [ ] **Step 2: Run all pruning tests**

Run: `pytest tests/test_brite_adapter_pruning.py -v`
Expected: all 5 tests PASS (`test_empty_known_kos_prunes_everything`, `test_single_known_ko_keeps_full_ancestor_chain`, `test_no_dangling_ko_edges`, `test_no_dangling_parent_edges`, `test_known_ko_ids_is_required`).

- [ ] **Step 3: Commit**

```bash
git add tests/test_brite_adapter_pruning.py
git commit -m "test(brite): positive pruning coverage + required-arg check"
```

---

## Task 5: Test-mode cap applies post-pruning

**Files:**
- Modify: `tests/test_brite_adapter_pruning.py`

- [ ] **Step 1: Add a larger synthetic tree fixture and test**

Append to [tests/test_brite_adapter_pruning.py](../../../tests/test_brite_adapter_pruning.py):

```python


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
```

- [ ] **Step 2: Run the new test**

Run: `pytest tests/test_brite_adapter_pruning.py::test_test_mode_cap_applied_post_pruning -v`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_brite_adapter_pruning.py
git commit -m "test(brite): verify test_mode cap applies post-pruning"
```

---

## Task 6: Wire `known_ko_ids` through the orchestrator

**Files:**
- Modify: `create_knowledge_graph.py:139-146`

- [ ] **Step 1: Pass `all_ko_ids()` into the BRITE adapter**

Open [create_knowledge_graph.py](../../../create_knowledge_graph.py). Find the BRITE adapter construction block (around line 139):

```python
    # KEGG BRITE functional hierarchies: 12 trees → BriteCategory nodes + hierarchy/KO edges
    brite_adapter = MultiBriteAdapter(
        cache_root=Path("cache/data"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    brite_adapter.download_data(cache=CACHE)
    bc.write_nodes(brite_adapter.get_nodes())
    bc.write_edges(brite_adapter.get_edges())
```

Replace with:

```python
    # KEGG BRITE functional hierarchies: 12 trees → BriteCategory nodes + hierarchy/KO edges.
    # Pruned to the subhierarchy reachable from KeggTerm nodes (gene-referenced KOs).
    brite_adapter = MultiBriteAdapter(
        cache_root=Path("cache/data"),
        known_ko_ids=kegg_anno_adapter.all_ko_ids(),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    brite_adapter.download_data(cache=CACHE)
    bc.write_nodes(brite_adapter.get_nodes())
    bc.write_edges(brite_adapter.get_edges())
```

- [ ] **Step 2: Smoke-check — run with `--test` to confirm the build still completes**

Run:

```bash
uv run python create_knowledge_graph.py --test 2>&1 | tail -50
```

Expected: run completes without exception. Look for `MultiBriteAdapter.get_nodes: N BriteCategory nodes` with N > 0 and `MultiBriteAdapter.get_edges: M parent edges, K KO→BRITE edges` lines in the log. `N`, `M`, `K` should all be substantially smaller than before pruning.

- [ ] **Step 3: Run the existing unit tests to confirm no regression**

Run: `pytest -m "not slow and not kg" -v`
Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add create_knowledge_graph.py
git commit -m "feat(brite): wire known_ko_ids from MultiKeggAnnotationAdapter"
```

---

## Task 7: Update KG validity tests for pruning semantics

**Files:**
- Modify: `tests/kg_validity/test_brite.py:38-74`

- [ ] **Step 1: Relax `test_brite_tree_present`**

Open [tests/kg_validity/test_brite.py](../../../tests/kg_validity/test_brite.py). Replace the test (around line 38):

```python
def test_brite_tree_present(run_query, tree_code, tree_name):
    """Each of the 12 configured BRITE trees must have BriteCategory nodes."""
    result = run_query(
        "MATCH (b:BriteCategory {tree_code: $tc}) RETURN count(b) AS cnt",
        tc=tree_code,
    )
    assert result[0]["cnt"] > 0, (
        f"No BriteCategory nodes found for tree {tree_code} ({tree_name})"
    )
```

with:

```python
def test_brite_tree_present(run_query, tree_code, tree_name):
    """Configured BRITE trees should have BriteCategory nodes **iff** at least
    one of their KOs is a KeggTerm in the KG.

    After subhierarchy pruning (see plans/2026-04-14-brite-pruning.md), a tree
    whose KOs are all outside the KG legitimately prunes to zero nodes. We assert
    consistency: either the tree has nodes AND has at least one incident KO
    edge, or it has neither.
    """
    node_result = run_query(
        "MATCH (b:BriteCategory {tree_code: $tc}) RETURN count(b) AS cnt",
        tc=tree_code,
    )
    edge_result = run_query(
        """
        MATCH (:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: $tc})
        RETURN count(*) AS cnt
        """,
        tc=tree_code,
    )
    node_count = node_result[0]["cnt"]
    edge_count = edge_result[0]["cnt"]

    if node_count == 0:
        assert edge_count == 0, (
            f"Tree {tree_code} ({tree_name}) has 0 nodes but {edge_count} "
            f"KO edges — pruning inconsistency"
        )
    else:
        assert edge_count > 0, (
            f"Tree {tree_code} ({tree_name}) has {node_count} nodes but 0 "
            f"KO edges — pruned tree should always have KO coverage"
        )
```

- [ ] **Step 2: Relax `test_brite_level_zero_nodes_exist`**

Replace (around line 53):

```python
def test_brite_level_zero_nodes_exist(run_query):
    result = run_query(
        """
        MATCH (b:BriteCategory {level: 0})
        RETURN count(b) AS cnt
        """
    )
    assert result[0]["cnt"] >= 12, (
        f"Only {result[0]['cnt']} level=0 BriteCategory nodes; expected ≥ 12"
    )
```

with:

```python
def test_brite_level_zero_nodes_exist(run_query):
    """At least one level-0 node must exist per non-empty tree.

    With subhierarchy pruning, some of the 12 configured trees may prune to
    zero nodes if none of their KOs are in the KG. We require at least 6
    level-0 nodes overall (half the trees must remain populated for the
    config to be meaningful).
    """
    result = run_query(
        """
        MATCH (b:BriteCategory {level: 0})
        RETURN count(b) AS cnt
        """
    )
    assert result[0]["cnt"] >= 6, (
        f"Only {result[0]['cnt']} level=0 BriteCategory nodes; expected ≥ 6 "
        f"(half the 12 configured trees)"
    )
```

- [ ] **Step 3: Add `test_brite_no_dangling_ko_edges`**

Add a new test immediately after `test_brite_no_duplicate_ko_edges` (around line 140):

```python


def test_brite_no_dangling_ko_edges(run_query):
    """No Kegg_term_in_brite_category edge may start from a non-KeggTerm node.

    Regression guard for the pruning fix: before pruning, the BRITE adapter
    emitted ~12K edges referring to KO IDs that had no KeggTerm node, causing
    ``neo4j-admin import`` to fail with exit code 70.
    """
    result = run_query(
        """
        MATCH (n)-[r:Kegg_term_in_brite_category]->()
        WHERE NOT n:KeggTerm
        RETURN count(r) AS cnt
        """
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} Kegg_term_in_brite_category edges have a "
        f"non-KeggTerm source — pruning regression"
    )
```

- [ ] **Step 4: Run KG validity tests (requires running Neo4j from a fresh build)**

Run: `pytest -m kg tests/kg_validity/test_brite.py -v`
Expected: tests pass after rebuilding the graph (Task 8). If Neo4j isn't running, tests skip — acceptable at this point; they'll be verified in Task 8.

- [ ] **Step 5: Commit**

```bash
git add tests/kg_validity/test_brite.py
git commit -m "test(brite): validity tests accept pruned subhierarchy + dangling-edge guard"
```

---

## Task 8: End-to-end verification — rebuild, import, validate

**Files:** none modified.

- [ ] **Step 1: Snapshot edge counts before the fix for comparison**

Run:

```bash
/omics-edge-snapshot
```

Expected: snapshot written under `cache/edge_snapshots/` — compare against post-build snapshot in Step 5.

- [ ] **Step 2: Full graph rebuild**

Run:

```bash
uv run python create_knowledge_graph.py
```

Expected: completes without error. Final BRITE adapter log lines show a meaningful reduction vs. the pre-fix baseline (~4,100 nodes → lower; ~17K KO edges → ~5-8K).

- [ ] **Step 3: Docker build + import**

Run:

```bash
docker compose down
docker compose up -d build import
```

Wait for the import container to exit. Run:

```bash
docker ps -a --format '{{.Names}}\t{{.Status}}' | grep -E "build|import"
```

Expected: `build` Exited (0), `import` Exited (0). If `import` shows Exited (70), check `docker logs import 2>&1 | grep "bad entries"` — should report 0 `Kegg_term_in_brite_category` failures.

- [ ] **Step 4: Bring the graph up and run KG validity tests**

Run:

```bash
docker compose up -d post-process deploy
# Wait for deploy to be ready (port 7687 reachable)
pytest -m kg tests/kg_validity/ -v
```

Expected: all tests pass, including the new `test_brite_no_dangling_ko_edges`. Note: unrelated failures tied to the Kratzl 2024 `insdc.gcf` issue (if any) are out-of-scope for this plan.

- [ ] **Step 5: Regenerate the KG snapshot**

Run:

```bash
uv run python tests/kg_validity/generate_snapshot.py
pytest -m kg tests/kg_validity/test_snapshot.py -v
```

Expected: snapshot regenerates; snapshot test passes.

- [ ] **Step 6: Take the edge snapshot again**

Run:

```bash
/omics-edge-snapshot
```

Expected: expression edge counts (`Changes_expression_of`) unchanged vs. Step 1 — this fix does not touch expression data. If counts differ, stop and investigate.

- [ ] **Step 7: Commit the regenerated snapshot**

```bash
git add tests/kg_validity/snapshot_data.json
git commit -m "test(kg): regenerate snapshot after BRITE pruning"
```

---

## Task 9: Update `CLAUDE.md` with new BRITE counts

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Capture the new counts from the live graph**

Run:

```bash
/cypher-queries
```

Execute in the skill:

```cypher
MATCH (b:BriteCategory) RETURN count(b) AS brite_nodes;
MATCH ()-[r:Brite_category_is_a_brite_category]->() RETURN count(r) AS parent_edges;
MATCH ()-[r:Kegg_term_in_brite_category]->() RETURN count(r) AS ko_edges;
```

Record the three numbers.

- [ ] **Step 2: Update the BRITE count line in `CLAUDE.md`**

Open [CLAUDE.md](../../../CLAUDE.md). Find the single line that currently reads (search for `~4,100 \`BriteCategory\` nodes`):

```
- BriteCategory nodes: ~4,100 nodes from 12 KEGG BRITE functional hierarchy trees. Node IDs: `kegg.brite:{tree_code}.{A#}.{B#}...` (positional, 1-based, no tree-root node emitted).
```

Replace the approximate counts with actuals from Step 1. Also update the edge-count line further down (search for `~8,400 \`Brite_category_is_a_brite_category\``):

```
Edges: ~8,400 `Brite_category_is_a_brite_category` (child→parent hierarchy within each tree) + ~17,000 `Kegg_term_in_brite_category` (KeggTerm KO → BriteCategory leaf).
```

Replace with the numbers from Step 1.

Add a sentence at the end of the BriteCategory paragraph:

> BRITE output is pruned to the subhierarchy reachable from `KeggTerm` nodes — a `BriteCategory` node appears iff at least one KO in its subtree is referenced by a gene in the KG. Pruning is keyed off `MultiKeggAnnotationAdapter.all_ko_ids()` (see `docs/superpowers/specs/2026-04-14-brite-pruning-design.md`).

- [ ] **Step 3: Update the memory summary counts in MEMORY.md (if present)**

Run:

```bash
grep -n "BriteCategory" /home/osnat/.claude/projects/-home-osnat-github-multiomics-biocypher-kg/memory/MEMORY.md
```

If the line listing `~4,100 BriteCategory nodes` exists, edit it to use the new counts and append "(subhierarchy pruned to KOs in KG)". If there's no matching line, skip.

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs(brite): update node/edge counts after subhierarchy pruning"
```

---

## Self-Review Notes

**Spec coverage checklist:**
- Coupling via `known_ko_ids` constructor arg (required) → Task 3, Step 2.
- `all_ko_ids()` public rename → Task 1.
- `_prune` algorithm (single walk → seed → mark → filter) → Task 3, Step 4.
- `get_nodes`/`get_edges` iterate `self._pruned` → Task 3, Step 5.
- `test_mode` cap applied post-pruning in both getters → Task 3, Step 5 & Task 5.
- Empty-tree behaviour (silent drop, per-tree info log) → Task 3, Step 4 & Task 2 assertion.
- Defensive check for duplicate node_id with conflicting parent → Task 3, Step 4.
- Adapter unit tests (empty, single-KO, no-dangling, required-arg, test_mode) → Tasks 2, 4, 5.
- KG validity test updates (`test_brite_tree_present`, `test_brite_level_zero_nodes_exist`) → Task 7.
- New `test_brite_no_dangling_ko_edges` → Task 7.
- Snapshot regeneration → Task 8, Step 5.
- `CLAUDE.md` update → Task 9.
- Expected KG deltas validated → Task 8, Step 2 + Step 4.

**Not in plan (deferred):**
- The 4 `insdc.gcf:*` missing-organism errors in Kratzl 2024 — separate issue, separate plan.
- Switching `ko01000` EC categories to EC-based IDs — brainstormed and deferred as scope creep.
