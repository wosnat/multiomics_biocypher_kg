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
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.utils.brite_utils import BRITE_TREES, compute_level_kind, load_brite_trees

logger = logging.getLogger(__name__)

# KO leaf detection: entries whose name starts with K followed by exactly 5 digits
_KO_RE = re.compile(r"^K\d{5}")
# Positional depth labels: index 0 → 'A', 1 → 'B', 2 → 'C', 3 → 'D'
_DEPTH_LETTERS = "ABCD"


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
        trees: "list[str] | None" = None,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.cache_root = Path(cache_root)
        self.trees = trees or list(BRITE_TREES.keys())
        self.test_mode = test_mode
        self.cache = cache
        self._tree_data: dict[str, dict] = {}

    def download_data(self, cache: bool = True) -> None:
        """Fetch/load all configured BRITE trees into memory."""
        self._tree_data = load_brite_trees(self.cache_root, self.trees, cache=cache)

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
