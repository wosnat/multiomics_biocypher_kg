# KG Change: Unified `level` on Ontology Terms

Date: 2026-04-12

## Summary

Add a canonical `level: int` property to every ontology-term node. Add an
optional sparse marker `level_is_best_effort: "true"` on GO DAG terms whose
min-path differs from their max-path to the namespace root. Rename the existing
`KeggTerm.level` (string) to `KeggTerm.level_kind`.

## Motivation

The multiomics_explorer enrichment surface needs to roll genes up to a chosen
hierarchy level across ten ontology labels. Without a uniform `level` property
the explorer dispatches four strategies (string lookup, dot-count, id-segment
parse, BFS) and has already hit Neo4j's 1.4 GiB transaction memory cap on a
GO query. One int property + one sparse marker replaces all four strategies.

## New Node Properties

| Label | `level: int` | `level_is_best_effort: "true"` |
|---|---|---|
| `BiologicalProcess`, `MolecularFunction`, `CellularComponent` | min-depth from namespace root | emitted iff min-path ≠ max-path |
| `EcNumber` | 0/1/2/3 by segment depth | never emitted |
| `KeggTerm` | 0=category, 1=subcategory, 2=pathway, 3=ko | never emitted |
| `CyanorakRole` | 0/1/2 by parent-chain depth | never emitted |
| `Pfam` | always 1 | never emitted |
| `PfamClan` | always 0 | never emitted |
| `TigrRole`, `CogFunctionalCategory` | always 0 | never emitted |

`level = 0` is always the broadest (root) level.

`level_is_best_effort` is only emitted when its value would be `"true"`; the
property is absent otherwise. Downstream consumers check with
`t.level_is_best_effort IS NOT NULL` (no coalesce needed). This avoids
past biocypher bugs with bool-typed properties.

## Renames on KeggTerm

| Before | After | Notes |
|---|---|---|
| `level: str` ∈ `{category, subcategory, pathway, ko}` | `level_kind: str` (same values) | semantic label |
| — | `level: int` ∈ `{0, 1, 2, 3}` | new hierarchy depth |

The `ONTOLOGY_CONFIG['kegg']['gene_connects_to_level'] = 'ko'` convention in
the explorer remains valid — gene→KEGG edges still terminate on the
leaf-kind node, which now also has `level = 3`.

## Derivation

Pure-adapter — no post-import Cypher changes.

- GO: Python BFS from canonical roots (`GO:0008150`, `GO:0003674`,
  `GO:0005575`) through `is_a` + `part_of` parents, done once in
  `go_utils.compute_go_levels`. `regulates*` is not traversed.
- EcNumber: structural — the adapter's four-level nested iteration assigns
  `level` directly per loop.
- CyanorakRole: structural — walk `parent` pointers in the parsed role tree.
- KEGG, Pfam, PfamClan, TigrRole, COG: fixed by kind.

## No Changes To

- Node labels or edge types.
- Hierarchy edges (all `*_is_a_*` and `*_part_of_*` relations retained).
- Post-import scripts.
- Docker compose / prepare_data.sh.

## KG Repo Files Changed

| File | Change |
|---|---|
| `config/schema_config.yaml` | `level` + `level_is_best_effort` added to 10 node types; `KeggTerm.level` renamed to `level_kind` |
| `multiomics_kg/utils/go_utils.py` | New `compute_go_levels` helper |
| `multiomics_kg/adapters/go_adapter.py` | Emit `level` (+ sparse `level_is_best_effort`) on GO nodes |
| `multiomics_kg/adapters/ec_adapter.py` | Emit `level` inside each nested loop |
| `multiomics_kg/adapters/functional_annotation_adapter.py` | KEGG/Cyanorak/Tigr/COG/Pfam/PfamClan emitters |
| `tests/test_go_levels.py` | New unit test for `compute_go_levels` |
| `tests/kg_validity/test_ontology_level.py` | New live-graph tests |
| `tests/test_kegg_annotation_adapter.py` | Updated for KEGG rename |
| `CLAUDE.md` | One bullet under "Actual Neo4j labels" documenting the convention |

## Consumer Impact

- Explorer's hierarchy helper: `WHERE t.level = $target_level` works uniformly
  across all ten labels. The four existing dispatch strategies can be
  collapsed into one query template.
- Anyone reading `KeggTerm.level` as a string must switch to `level_kind` (or
  start using the new integer `level`). Coordinated via KG-rebuild checkpoint
  as documented in the parent spec.

## Rollout

Single KG rebuild. Pure addition except for the KEGG rename; no migration
script.
