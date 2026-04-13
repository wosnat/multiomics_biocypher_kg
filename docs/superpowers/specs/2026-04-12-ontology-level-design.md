# Unified `level` on Ontology Terms — KG Implementation Design

Date: 2026-04-12

## Context

The multiomics_explorer's enrichment surface (see
`multiomics_explorer/docs/kg-specs/2026-04-12-ontology-hierarchy-level.md`)
needs to roll genes up to a chosen hierarchy level across ten ontology labels.
Today the explorer dispatches four different strategies per ontology (string
lookup, dot-count, id-segment parse, BFS), and query-time BFS over GO has
already hit Neo4j's 1.4 GiB transaction memory cap on one workload.

This design implements the KG-side half of that plan: precompute a canonical
`level: int` property on every ontology-term node at KG-build time.

## Scope

**In scope.** Add two properties to ten node labels:

- `level: int` on every term.
- `level_is_best_effort: "true"` (str, sparse — only emitted on GO DAG terms
  where `min_depth != max_depth` across root paths).

Labels affected: `BiologicalProcess`, `MolecularFunction`, `CellularComponent`,
`EcNumber`, `KeggTerm`, `CyanorakRole`, `TigrRole`, `CogFunctionalCategory`,
`Pfam`, `PfamClan`.

One additional rename on `KeggTerm`: the existing string `level`
(`{ko, pathway, subcategory, category}`) becomes `level_kind`; the new integer
`level` takes its place with mapping `category=0, subcategory=1, pathway=2,
ko=3`.

**Out of scope.**

- Ancestor closure (`HAS_ANCESTOR` / `ANNOTATED_AT`) — parent spec's phases 2/3,
  deferred pending profiling of `level`-pruned BFS.
- Per-term gene counts.
- Explorer-side consumer changes — coordinated via KG-rebuild checkpoint.

## Approach

Pure-adapter implementation. No post-import Cypher changes. Every `level` value
is written by the adapter that owns the node's source data; each adapter
already has the structural information needed.

| Label | Source | `level` derivation |
|---|---|---|
| `BiologicalProcess`, `MolecularFunction`, `CellularComponent` | `go_utils.load_go_data()` | BFS from canonical roots (`GO:0008150`, `GO:0003674`, `GO:0005575`) through `is_a` + `part_of` parents; min-depth |
| `EcNumber` | `ec_adapter.get_nodes()` four-level nested loop | Nesting depth: level_1→0, level_2→1, level_3→2, level_4→3 |
| `KeggTerm` | `functional_annotation_adapter` KEGG emitters | Fixed by kind: category=0, subcategory=1, pathway=2, ko=3 |
| `CyanorakRole` | `parse_cyanorak_role_tree()` — `{code: {parent}}` | Walk `parent` pointers to root |
| `Pfam` | `functional_annotation_adapter` Pfam emitter | Fixed = 1 (including 1,935 clan-less Pfams) |
| `PfamClan` | `functional_annotation_adapter` PfamClan emitter | Fixed = 0 |
| `TigrRole` | flat | 0 |
| `CogFunctionalCategory` | flat | 0 |

`level_is_best_effort` is set (`"true"`) only by the GO adapter, only on
terms whose BFS encountered both a shorter and a longer root path. All other
labels omit the property entirely; explorer queries use
`t.level_is_best_effort IS NOT NULL` rather than coalesce.

We emit only when true to sidestep past biocypher bugs around bool-typed
properties and to keep the CSV export tight.

### GO BFS algorithm

One pass per namespace:

1. Seed frontier with the canonical root (level 0).
2. For each term popped, scan children (terms that list it in `parents`);
   update each child's `min_depth = min(current, parent_depth + 1)` and
   `max_depth = max(current, parent_depth + 1)`.
3. After convergence, emit `level = min_depth`, and
   `level_is_best_effort = "true"` iff `min_depth != max_depth`.
4. Any term not reached from its namespace root is an **orphan** — logged and
   excluded from node emission. Expected count: 0.

Traversed relations: `is_a` and `part_of` only — not `regulates*`, which are
not hierarchy-carrying under the explorer's semantics.

## Components

### Schema (`config/schema_config.yaml`)

- Add `level: int` to all ten node types.
- Add `level_is_best_effort: str` to all ten node types (sparse; explorer uses
  `IS NOT NULL`).
- On `kegg term`: rename `level: str` → `level_kind: str`, then add `level: int`.

### `multiomics_kg/utils/go_utils.py`

New public helper:

```python
def compute_go_levels(
    go_data: dict[str, dict],
) -> tuple[dict[str, tuple[int, bool]], list[str]]:
    """
    For each GO term reachable from its namespace root, return
    (min_depth, is_best_effort). Also return a list of orphan term IDs.
    """
```

Unit-tested separately from adapter integration. Canonical root IDs
(`GO:0008150`, `GO:0003674`, `GO:0005575`) are module-level constants.

### `multiomics_kg/adapters/go_adapter.py`

Call `compute_go_levels(go_data)` once during node materialisation. On each
`(node_id, label, props)` emit, attach `props["level"] = min_depth`, and
`props["level_is_best_effort"] = "true"` only when the flag is set.

### `multiomics_kg/adapters/ec_adapter.py`

In `get_nodes()` ([ec_adapter.py:218-286](/home/osnat/github/second_multiomics/multiomics_kg/adapters/ec_adapter.py#L218-L286)), set `props["level"] = 0/1/2/3` inside each
of the four nested loops.

### `multiomics_kg/adapters/functional_annotation_adapter.py`

- KEGG node emitters (`_ko_node_id`, `_pathway_node_id`, `_subcat_node_id`,
  `_cat_node_id` at [functional_annotation_adapter.py:601-630](/home/osnat/github/second_multiomics/multiomics_kg/adapters/functional_annotation_adapter.py#L601-L630)):
  rename the existing `"level"` key to `"level_kind"`; add
  `"level": 3/2/1/0` respectively.
- `MultiCogRoleAnnotationAdapter.get_nodes()`: add a local
  `_role_depth(code, role_tree)` helper that walks `parent` pointers; set
  `props["level"]` for every CyanorakRole node; set `level = 0` on TigrRole
  and CogFunctionalCategory nodes.
- Pfam emitter: set `props["level"] = 1`.
- PfamClan emitter: set `props["level"] = 0`.

### Tests

**Unit test** — new `tests/test_go_levels.py`:

- Three-root linear chain → monotonic levels, no flag set.
- Diamond with equal-length arms → no flag.
- Diamond with unequal arms → flag set on confluence node, `level = min`.
- `is_a` and `part_of` both traversed; `regulates` ignored.
- Orphan term (no parents, not a canonical root) excluded + surfaced in the
  orphan list returned by `compute_go_levels`.

**KG-validity tests** — new `tests/kg_validity/test_ontology_level.py`:

1. Coverage: every node of the ten labels has non-null `level`.
2. Flat ontologies: `TigrRole` (n=114), `CogFunctionalCategory` (n=26), all
   `level = 0`.
3. `Pfam`: all `level = 1` (count=5,471 including clan-less); `PfamClan`: all
   `level = 0` (count=509).
4. Roots: `CyanorakRole {level: 0}` count = 19; `EcNumber {level: 0}` count = 7;
   `KeggTerm {level: 0}` count = 6.
5. GO exactly one level-0 term per namespace.
6. GO min-depth invariant: for every non-root GO term `t`, `t.level = 1 +
   min(p.level)` over its `is_a`/`part_of` parent edges.
7. GO BP level distribution matches spec table (16/98/274/565/721/650/462/211/
   42/10/2 at depths 1–11). Warn on ≤5% drift; fail on >5% drift.
8. GO `level_is_best_effort = "true"` exists on some but not all terms; at
   least one known multi-depth term (hard-coded at test-write time) is flagged.
9. `KeggTerm`: every node has both `level_kind` and `level`; mapping
   `{category:0, subcategory:1, pathway:2, ko:3}` holds for every node.

**Existing test update** — `tests/test_kegg_annotation_adapter.py`:

- Any assertion on the string `"level"` property becomes `"level_kind"`.
- Add assertions for the new integer `"level"`.

### Rollout

Single KG rebuild lands the change. The parent spec on the explorer side
already commits to coordinating via this rebuild checkpoint. No migration
script, no post-import script changes.

## Change note for downstream consumers

Summarised in `docs/kg-changes/ontology-level.md` per the repo convention.

## Verification

Per the parent spec's Verification section, plus the distribution and
invariant checks listed under Tests above. All run as pytest against the live
deployed Neo4j instance; auto-skip when Neo4j is unreachable.

## Decision log

- **Hybrid vs pure-adapter:** pure-adapter chosen. GO BFS is cheap (~6k terms)
  in Python because `go_utils` already holds the full `parents` dict in memory,
  and keeping all `level` emissions on one side of the build boundary removes
  a whole class of drift bug.
- **`level_source` vs `level_is_best_effort`:** dropped `level_source` in favour
  of the per-term flag. More precise (a GO term with two equal-length parent
  paths is not ambiguous) and self-describing in the graph.
- **Bool vs sparse string for the flag:** sparse string. Past biocypher bugs
  with bool-typed properties; omitting the key when false also keeps CSV output
  tight.
- **GO root handling:** canonical three GO roots only. Orphans are a data
  anomaly and surface as a pytest failure rather than silently getting
  `level=0`.
- **CyanorakRole / EcNumber derivation:** structural (parent chain, nested
  loop depth) rather than string parsing. The adapters already expose the
  structure.
