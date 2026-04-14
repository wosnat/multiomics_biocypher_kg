"""
BRITE adapter: ingests 12 KEGG BRITE functional hierarchy trees.

Yields:
- BriteCategory nodes for every non-KO entry in all 12 trees
- Brite_category_is_a_brite_category edges (child → parent, within tree)
- Kegg_term_in_brite_category edges (KeggTerm KO → BRITE leaf category)

Level assignment — no synthetic tree-root node; A-level entries are broadest:
  A-level entry → level=0, level_kind='brite_class',    no parent edge
  B-level entry → level=1, level_kind='brite_subclass', parent = A-level
  C-level entry → level=2, level_kind='brite_family',   parent = B-level
  D-level non-KO→ level=3, level_kind='brite_subfamily',parent = C-level

Node ID scheme (positional, 1-based, stable):
  A-level entry 3       → kegg.brite:ko02000.A3
  B-level entry 2 of A3 → kegg.brite:ko02000.A3.B2
  C-level entry 1       → kegg.brite:ko02000.A3.B2.C1
"""

from __future__ import annotations

import html
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.utils.brite_utils import BRITE_TREES, compute_level_kind, load_brite_trees

logger = logging.getLogger(__name__)

# KO leaf detection: entries whose name starts with K followed by exactly 5 digits
_KO_RE = re.compile(r"^K\d{5}")
# Positional depth labels: index 0 → 'A', 1 → 'B', 2 → 'C', 3 → 'D'
_DEPTH_LETTERS = "ABCD"


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


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV: replace single quotes → ^, strip pipes."""
    return value.replace("'", "^").replace("|", "")


def _brite_node_id(tree_code: str, path: list[int]) -> str:
    """Build a BriteCategory node ID from tree code and 1-based positional path.

    path=[]       → kegg.brite:ko02000  (tree namespace; never emitted as a node)
    path=[3]      → kegg.brite:ko02000.A3
    path=[3, 2]   → kegg.brite:ko02000.A3.B2
    path=[3, 2, 1]→ kegg.brite:ko02000.A3.B2.C1
    """
    if path:
        path_str = ".".join(
            f"{_DEPTH_LETTERS[i]}{path[i]}" for i in range(len(path))
        )
        raw = f"kegg.brite:{tree_code}.{path_str}"
    else:
        raw = f"kegg.brite:{tree_code}"
    return normalize_curie(raw) or raw


def _ko_node_id(ko_id: str) -> str:
    """KeggTerm node ID for a KO accession — must match MultiKeggAnnotationAdapter."""
    raw = f"kegg.orthology:{ko_id}"
    return normalize_curie(raw) or raw


class MultiBriteAdapter:
    """
    Multi-tree adapter: loads 12 KEGG BRITE trees and yields BriteCategory
    nodes plus two edge streams.

    Args:
        cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
        trees: list of KEGG tree codes; defaults to all 12 configured trees
        test_mode: cap at 100 nodes per tree for fast iteration during development
        cache: if False, re-download trees even if cache files exist
    """

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
        self.known_ko_ids = frozenset(known_ko_ids)
        self.trees = trees or list(BRITE_TREES.keys())
        self.test_mode = test_mode
        self.cache = cache
        self._tree_data: dict[str, dict] = {}
        self._pruned: dict[str, PrunedTree] = {}

    def download_data(self, cache: bool = True) -> None:
        """Fetch/load all configured BRITE trees into memory and prune them."""
        self._tree_data = load_brite_trees(self.cache_root, self.trees, cache=cache)
        self._pruned = {
            tree_code: self._prune(tree_code, tree_json)
            for tree_code, tree_json in self._tree_data.items()
        }

    def _walk(self, tree_code: str, children: list, path: list[int]):
        """Recursively walk BRITE JSON children, yielding typed tuples.

        Yields:
            ("node",    node_id, parent_id_or_None, name, tree_code, level, level_kind)
            ("edge_ko", ko_id_raw, parent_brite_node_id)
        """
        for i, child in enumerate(children, start=1):
            raw_name = child.get("name", "")
            if not raw_name:
                logger.warning(
                    f"BRITE tree {tree_code}: entry at depth {len(path)}, "
                    f"index {i} has no name, skipping"
                )
                continue

            name = html.unescape(raw_name)
            current_path = path + [i]

            if _KO_RE.match(name):
                # KO leaf: emit a KO→BRITE edge; parent must exist
                if path:
                    ko_id_raw = name.split()[0]  # "K02036  pstB; ..." → "K02036"
                    yield ("edge_ko", ko_id_raw, _brite_node_id(tree_code, path))
                else:
                    logger.warning(
                        f"BRITE tree {tree_code}: KO leaf {name!r} at root level "
                        f"(no parent category), skipping edge"
                    )
            else:
                # BriteCategory node
                node_id = _brite_node_id(tree_code, current_path)
                parent_id = _brite_node_id(tree_code, path) if path else None
                depth = len(current_path) - 1  # A-level=0, B-level=1, C-level=2, D-level=3
                level = depth
                level_kind = compute_level_kind(depth)
                yield (
                    "node",
                    node_id,
                    parent_id,
                    _clean_str(name),
                    tree_code,
                    level,
                    level_kind,
                )

                sub_children = child.get("children", [])
                if sub_children:
                    yield from self._walk(tree_code, sub_children, current_path)

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

    def get_nodes(self):
        """Yield ``(node_id, 'brite category', properties)`` for every surviving BriteCategory."""
        if not self._pruned:
            self.download_data(cache=self.cache)

        total_count = 0
        for tree_code, pt in self._pruned.items():
            nodes_to_emit = pt.nodes[:100] if self.test_mode else pt.nodes
            if self.test_mode and len(pt.nodes) > 100:
                logger.debug(
                    f"MultiBriteAdapter.get_nodes: test_mode cap at 100 for {tree_code}"
                )
            for node_id, _parent_id, name, level in nodes_to_emit:
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
                total_count += 1

        logger.info(f"MultiBriteAdapter.get_nodes: {total_count} BriteCategory nodes")

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
