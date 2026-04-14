# BRITE Subhierarchy Pruning — Design

**Date:** 2026-04-14
**Status:** design, pre-implementation
**Scope:** `MultiBriteAdapter`, `MultiKeggAnnotationAdapter` (minor), orchestrator, tests, docs.
**Non-scope:** the 4 `insdc.gcf` missing-organism errors from Kratzl 2024 (separate issue, tracked independently).

## Problem

The Docker `import` container exits with code 70. `import.report` contains 1001 bad relationship entries (fault-tolerance cap), of which **996** are `Kegg_term_in_brite_category` edges referring to missing `kegg.orthology:K*` nodes. The BRITE adapter emits ~17,000 KO→category edges covering all KEGG KOs in the 12 configured BRITE trees, but the KG only materialises ~4,700 `KeggTerm` nodes (those referenced by gene annotations). The ~12K delta dangles.

Beyond the dangling edges, the current output also includes thousands of `BriteCategory` nodes whose entire subtree contains no KO that any gene in our KG uses. These add noise to queries and bloat the graph.

## Goal

`MultiBriteAdapter` should emit only the **subhierarchy relevant to our KG**: a `BriteCategory` node appears iff at least one KO in its subtree is a `KeggTerm` in the KG; KO→category edges appear iff the KO is a `KeggTerm`.

## Design

### Coupling

`MultiBriteAdapter` requires a `known_ko_ids: set[str]` constructor argument (no default). `create_knowledge_graph.py` sources it from `MultiKeggAnnotationAdapter.all_ko_ids()` (renamed from the private `_all_ko_ids`).

The KEGG adapter already prunes its own hierarchy to ancestors of known KOs (`pw_ids` / `subcat_ids` / `cat_ids` derived from `ko_ids`); BRITE is the symmetric pattern on a nested-tree data shape.

### Data shape

- `BriteCategory` internal hierarchy: **strict tree** by construction. Node IDs are positional (`kegg.brite:ko02000.A3.B2.C1`), so a category that appears under multiple parents in the raw JSON becomes multiple distinct nodes — correctly, since same-name duplicates at different positions are semantically distinct categories (e.g. "Small subunit" under bacterial/archaeal/eukaryotic ribosome A-levels are non-homologous).
- `KeggTerm ↔ BriteCategory`: **many-to-many** — a KO can sit under several BRITE leaves; KeggTerm nodes are already deduplicated by `kegg.orthology:K#####`.
- Overall subgraph in the KG: **DAG** — fine, algorithm below handles it.

### Pruning algorithm

One pre-processing pass per tree, performed at the end of `MultiBriteAdapter.download_data()`. Result stored on `self._pruned: dict[str, PrunedTree]`; `get_nodes()` and `get_edges()` iterate that structure.

```python
@dataclass(slots=True)
class PrunedTree:
    tree_code: str
    nodes: list[tuple[str, str | None, str, int, str]]
        # (node_id, parent_id, name, level, level_kind) — discovery order
    ko_edges: list[tuple[str, str]]
        # (ko_id_raw, leaf_node_id) — only for KOs in known_ko_ids
```

Per-tree pruning:

1. **Walk once.** Reuse the existing `_walk` generator to collect:
   - `all_nodes`: every BriteCategory discovered (with parent_id, name, level, level_kind).
   - `all_ko_edges`: every `(ko_id_raw, leaf_node_id)` pair.
   - `parent_of: dict[str, str]`: node_id → parent_id (total map — each BRITE category has exactly one parent by positional-ID construction).
2. **Seed leaves.** `leaf_hits = {leaf_id for (ko, leaf_id) in all_ko_edges if ko in known_ko_ids}`.
3. **Mark keepers bottom-up.** For each `leaf_id` in `leaf_hits`, add it to `kept: set[str]`, then walk `parent_of` upward, adding each ancestor to `kept`, until either `parent_of[current]` is undefined (reached an A-level node with no parent) or the current node is already in `kept` (another seed already marked the rest of the chain).
4. **Materialise.** Build `PrunedTree(nodes=[n for n in all_nodes if n.id in kept], ko_edges=[(ko, lid) for (ko, lid) in all_ko_edges if ko in known_ko_ids and lid in kept])`.
5. **Empty tree.** If `len(nodes) == 0`, log `f"BRITE tree {tree_code}: 0/{len(all_nodes)} categories kept — pruned entirely"` and store the empty `PrunedTree`. Downstream iteration yields nothing. No failure.

Defensive check inside the walker: if a `node_id` would be assigned twice with different `parent_id`s, log `logger.error` and skip the duplicate. Guards against any hypothetical DAG introduction on the category side.

### Getters

```python
def get_nodes(self):
    for tree_code, pt in self._pruned.items():
        for i, (nid, _pid, name, level, lk) in enumerate(pt.nodes):
            if self.test_mode and i >= 100:
                break
            yield (
                nid,
                "brite category",
                {
                    "name": name,
                    "tree": BRITE_TREES[tree_code],
                    "tree_code": tree_code,
                    "level": level,
                    "level_kind": lk,
                },
            )

def get_edges(self):
    for tree_code, pt in self._pruned.items():
        nodes_to_emit = pt.nodes[:100] if self.test_mode else pt.nodes
        emitted_ids = {n[0] for n in nodes_to_emit}
        # parent edges between kept nodes
        for nid, pid, *_ in nodes_to_emit:
            if pid is not None and pid in emitted_ids:
                yield (f"{nid}--parent", nid, pid, "brite_category_is_a_brite_category", {})
        # KO edges (already filtered at pruning time)
        for ko_raw, leaf_id in pt.ko_edges:
            if leaf_id in emitted_ids:
                yield (
                    f"{ko_raw}--brite--{leaf_id}",
                    _ko_node_id(ko_raw),
                    leaf_id,
                    "kegg_term_in_brite_category",
                    {},
                )
```

The test-mode cap is applied **after** pruning. It's a safety net — pruned output is already small by construction.

## Impact on the KEGG adapter

The KEGG adapter already performs the same kind of pruning (upward propagation from `ko_ids` through flat `ko_to_pathways` / `pathway_to_subcategory` / `subcategory_to_category` maps). No structural change needed. Two cosmetic changes:

1. Rename `MultiKeggAnnotationAdapter._all_ko_ids()` → `all_ko_ids()` (public), with a one-line docstring noting it's the canonical source of KO IDs for downstream pruning (BRITE).
2. Add a module-level note to both `functional_annotation_adapter.py` (around `MultiKeggAnnotationAdapter`) and `brite_adapter.py` (around `MultiBriteAdapter`) stating they follow the shared "keep only ancestors of known KOs" convention.

## Tests

### Adapter unit tests (new file `tests/test_brite_adapter_pruning.py`)

Using synthetic 2-tree JSON fixtures:

- `test_empty_known_kos` — every tree prunes to empty; `get_nodes()`/`get_edges()` yield nothing; info log emitted per tree.
- `test_single_ko_seed` — `known_ko_ids = {"K00001"}`; verify the full root-to-leaf ancestor chain for every BRITE leaf containing K00001 is kept, and nothing else.
- `test_no_dangling_ko_edges` — every emitted `Kegg_term_in_brite_category` edge has its raw KO in `known_ko_ids`.
- `test_no_dangling_parent_edges` — every parent edge has both endpoints in the emitted node set.
- `test_empty_tree_logged` — tree with zero used KOs logs the coverage summary.
- `test_mode_cap_post_pruning` — with large synthetic tree, 150 kept nodes, `test_mode=True` emits ≤100 node and parent edges stay consistent.
- `test_missing_known_ko_ids_is_type_error` — constructor without `known_ko_ids` raises at construction (required arg).

### KG validity tests (`tests/kg_validity/test_brite.py`)

- Update `test_brite_tree_present` — allow zero nodes for trees whose KOs are all outside the KG; assert that every tree with ≥1 used KO has ≥1 node.
- Update `test_brite_level_zero_nodes_exist` — lower-bound becomes ≥1 per non-empty tree; total ≥6 (allow up to 6 trees to be entirely pruned).
- New `test_brite_no_dangling_ko_edges` —
  ```cypher
  MATCH ()-[r:Kegg_term_in_brite_category]->()
  WHERE NOT (startNode(r):KeggTerm) RETURN count(r)
  ```
  must be 0.
- Re-verify `test_brite_transporters_ko_count`, `test_brite_peptidases_ko_count`, `test_brite_ftsh_in_peptidases` — expected still to pass since pruning keeps exactly the leaves that were previously valid.

### Snapshot

After rebuild, regenerate `tests/kg_validity/snapshot_data.json` via `uv run python tests/kg_validity/generate_snapshot.py` and commit. Existing snapshot entries referencing pruned BriteCategory nodes would otherwise fail `test_snapshot.py`.

## Post-import interactions

Post-import Cypher ([scripts/post-import.sh](../../../scripts/post-import.sh), [scripts/post-import.cypher](../../../scripts/post-import.cypher)) computes `member_ko_count`, `gene_count`, `organism_count` on `BriteCategory`. It operates on whatever nodes survive import — no changes needed. Recursive KO counts sum only kept descendants by construction. Full-text and scalar indexes on BriteCategory continue to apply.

## Expected KG deltas

- `BriteCategory` nodes: ~4,100 → expected notably lower (exact number depends on per-tree KO coverage; likely in the range 1,500–2,500).
- `Brite_category_is_a_brite_category` edges: ~8,400 → proportional drop.
- `Kegg_term_in_brite_category` edges: ~17,000 → ~5,000–8,000 (edges preserved only for the ~4,700 known KOs, with the same many-to-many factor).
- Import report: zero `Kegg_term_in_brite_category` bad entries. (The 4 `insdc.gcf:*` entries from Kratzl 2024 remain — separate fix.)

## Docs updates

- [CLAUDE.md](../../../CLAUDE.md) BriteCategory section: refresh node/edge counts, add one line that BRITE content is pruned to the subhierarchy reachable from `KeggTerm` nodes, reference this spec.
- No changes to `docs/kg-changes/brite-categories.md` (it describes the broader feature; a small addendum on pruning may be useful but is optional).

## Open questions — none

All design decisions resolved during brainstorming:
- Scope: full subtree pruning.
- Source: `MultiKeggAnnotationAdapter.all_ko_ids()` via orchestrator.
- Argument: required constructor param, `known_ko_ids: set[str]`.
- Empty tree: drop silently, log coverage.
- Node IDs: positional (unchanged).
- `test_mode`: cap applied post-pruning (safety net).
