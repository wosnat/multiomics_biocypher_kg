# <Tool> ontology (<NodeLabel> nodes) — TEMPLATE

<!-- Copy to docs/kg-changes/<tool>-extension.md. Replace every <…>. Modeled on
     docs/kg-changes/tcdb-cazy-ontologies.md. This is the release-notes / what-changed doc —
     Step 6, the durable hand-off artifact (incl. to the out-of-scope MCP follow-up). -->

**Date:** <YYYY-MM-DD>
**Spec:** [`docs/superpowers/specs/<date>-<tool>-kg-integration-design.md`](../superpowers/specs/<date>-<tool>-kg-integration-design.md)
**Track:** <3A ontology | 3B Gene property> · <flat | hierarchical> · <scored | bare> · <functional | structural>

## What's changing

One paragraph: what the tool predicts, what graph elements it adds, and the merged-annotation field
it now flows through (`gene_annotations_merged.json` ← `<tool>` source).

## New node label: `<NodeLabel>`

- Count: <N> nodes. Levels: <flat: 1 / hierarchical: list>. ID prefix `<prefix>:` (e.g. `<prefix>:OuterMembrane`).
- Properties (adapter-emitted): `name`, `<tool>_id`, `level` (0=root)<, `level_kind` if hierarchical>.
- Computed (post-import): `gene_count`, `organism_count`<, `member_count` if hierarchical>.
- Vocabulary (skipped sentinels: <Unknown / Other>):
  | id | name |
  |---|---|
  | <RAW> | <display> |

## New edge types

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_<x>` | Gene → <NodeLabel> | <N> | <`score: float` / none `{}`> |
| `<X>_is_a_<x>` (hierarchical only) | <NodeLabel> → <NodeLabel> | <N> | none |

## Cross-node rollups

- **<functional>**: folded into `Gene.annotation_types` + `Gene.informative_annotation_types` as
  `<source>`<; added to the `annotation_quality` 8-bucket count — shifts `genes_by_function` `min_quality`>.
- **<structural>**: kept OUT of `annotation_types`/`annotation_quality`. Optional `Gene.<tool>_count` routing signal.
- `OrganismTaxon` / `Publication`: no per-ontology rollup (unchanged).

## Indexes

Scalar: `<node>_level_idx`, `<node>_id_idx`<, `<node>_level_kind_idx`>. Full-text: `<node>FullText` on `name`, `<tool>_id`.

## Validation results

- `/omics-edge-snapshot` before/after: expression edges unchanged (<N>); new `Gene_has_<x>` edges appear (<N>).
- `pytest -m kg`: <pass> (new assertions: node count > 0, edge presence, <score ∈ [lo,hi]>, no-orphan).
- `snapshot_data.json` regenerated.
- `post-import-validate.sh` baseline/after diff: additions-only.

## What this enables

- <e.g. "all OuterMembrane genes in organism X" — a one-hop MATCH (g:Gene)-[:Gene_has_<x>]->(:<NodeLabel> {<tool>_id:'OuterMembrane'})>
- <score-ranked shortlists via rank_by_score, if scored>

## What does NOT change

- Phase-1 (`<tool>-run`) untouched. No existing node/edge/property modified (additions-only).
- **MCP / explorer integration is out of scope** — querying this ontology via `ontology_landscape` /
  `search_ontology` / `genes_by_ontology` requires adding it to the explorer's `ONTOLOGY_CONFIG`
  (`label` + `gene_rel` + `hierarchy_rels` + `fulltext_index`) + the 12→13 `ontology` enum, tracked in
  the `multiomics_explorer` repo. **Follow-up owner: <…>.**

## Graph-size impact (measured)

| | before | after |
|---|---|---|
| `<NodeLabel>` nodes | 0 | <N> |
| `Gene_has_<x>` edges | 0 | <N> |

## See also

- Spec: `docs/superpowers/specs/<date>-<tool>-kg-integration-design.md`
- Skill: `.claude/skills/integrate-a-tool/SKILL.md`
- Phase-1 runner: `.claude/skills/<tool>-run/SKILL.md`
