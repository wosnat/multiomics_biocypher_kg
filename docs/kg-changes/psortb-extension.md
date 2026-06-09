# PSORTb subcellular localization (`SubcellularLocalization` nodes)

**Date:** 2026-05-26
**Spec:** [`docs/superpowers/specs/2026-05-26-psortb-kg-integration-design.md`](../superpowers/specs/2026-05-26-psortb-kg-integration-design.md)
**Track:** 3A ontology · flat · scored · **structural** (where the gene *is*, not what it does)
**Skill:** produced with [`/integrate-a-tool`](../../.claude/skills/integrate-a-tool/SKILL.md) (Phase-2 partner to `/add-a-tool`)

## What's changing

PSORTb v3.0.3 predicts the subcellular localization of each protein (Gram-negative model).
Phase 1 (`/psortb-run`) wrote committed `<strain>.psortb.calls.json` artifacts; this change wires
them into the KG. The localization call + confidence score flow through the gene-annotation merge
(`psortb` source → `psortb_localization` / `psortb_score` fields in `gene_annotations_merged.json`,
joined on the RefSeq WP_ `protein_id`), then surface as a small **flat ontology**: the five real
PSORTb classes become `SubcellularLocalization` nodes, and each gene with a confident call gets one
`Gene_has_subcellular_localization` edge carrying its `score`.

## New node label: `SubcellularLocalization`

- **Count:** 5 nodes (all real classes observed across 40 organisms). Flat — single level.
  Post-import `gene_count`: Cytoplasmic 49,251 · CytoplasmicMembrane 25,118 · OuterMembrane 2,049 ·
  Periplasmic 1,871 · Extracellular 1,072.
- **ID prefix `psortb_`** (underscore, not colon: `psortb` is not a registered bioregistry prefix,
  so `normalize_curie('psortb:X')` returns `None` and the id falls back to `psortb_OuterMembrane`).
- **Properties (adapter-emitted):** `name`, `psortb_id`, `level` (always `0` — flat, no `level_kind`).
- **Computed (post-import):** `gene_count`, `organism_count`.
- **Vocabulary** (the `Unknown` sentinel is skipped — absence of an edge encodes "no confident call"):

  | psortb_id | name |
  |---|---|
  | `Cytoplasmic` | Cytoplasmic |
  | `CytoplasmicMembrane` | Cytoplasmic membrane |
  | `Periplasmic` | Periplasmic |
  | `OuterMembrane` | Outer membrane |
  | `Extracellular` | Extracellular |

## New edge type

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_subcellular_localization` | Gene → SubcellularLocalization | 79,361 | `score: float` (PSORTb confidence, observed range 7.5–10.0); post-import `rank_by_score: int` |

This is the KG's **first scored ontology edge** — the `score` rides on the edge (no existing
gene→ontology edge carries properties), modeled on the `Changes_expression_of` property block.
1:1 (≤1 edge per gene); no hierarchy (`*_is_a_*`) edges.

## Cross-node rollups

- **Structural, not functional:** deliberately **NOT** folded into `Gene.annotation_types` /
  `informative_annotation_types` / `annotation_quality` (localization is *where* a gene is, not
  *what it does* — folding it in would inflate functional-annotation coverage and skew
  `genes_by_function` `min_quality`).
- **`Gene.subcellular_localization`** (post-import, str): denormalized 1:1 routing string = the
  gene's single PSORTb class (absent when no confident call). Surfaces via `gene_details` for free.
- `OrganismTaxon` / `Publication`: no per-ontology rollup (unchanged).

## Indexes

Scalar: `subcellular_localization_level_idx`, `subcellular_localization_id_idx`.
Full-text: `subcellularLocalizationFullText` on `name`, `psortb_id`.

## Merge wiring + provenance

- `config/gene_annotations_config.yaml`: `psortb` source (join `protein_id`) + `psortb_localization`
  / `psortb_score` fields + a `psortb` `logical_source` (`scope: gene_level`, `provenance: tool_run`).
- `multiomics_kg/download/build_gene_annotations.py`: `load_psortb()` + the source threaded through
  `_get_raw` / resolvers / `build_wide` / `build_merged` / the `process_strain` row join /
  `_compute_contributing_sources`.
- A `DataSource` node `data_source:psortb` emits automatically (named via
  `data_source_adapter.py`'s `_name_for`/`_description_for`); `Gene.contributing_sources` lists `psortb`.

## Validation results

- `/omics-edge-snapshot` before/after: expression (`Changes_expression_of`) edge counts **unchanged**
  (232,758; no publication lost edges).
- `pytest -m kg`: **967 passed, 4 skipped** — incl. `tests/kg_validity/test_psortb.py` (node set ⊆ the
  5 classes, no `Unknown`, flat level 0, `score` ∈ [7.5, 10.0], `rank_by_score` present, no orphan
  edges, routing-string consistency, psortb kept out of `annotation_types`) and the DataSource
  node-count test updated 4 → 5.
- `pytest -m "not slow and not kg"`: green (incl. `tests/test_psortb.py`, 10 unit tests).
- `post-import-validate.sh` baseline/after diff: additions-only.
- `snapshot_data.json` regenerated.

## What this enables

- "All outer-membrane genes in organism X" — one hop:
  `MATCH (g:Gene)-[:Gene_has_subcellular_localization]->(:SubcellularLocalization {psortb_id:'OuterMembrane'})`.
- Score-ranked shortlists per localization via `rank_by_score`.
- A `Gene.subcellular_localization` routing string for gene-overview / detail surfaces.

## What does NOT change

- Phase 1 (`/psortb-run`) untouched. No existing node/edge/property modified (additions-only).
- **MCP / explorer integration is out of scope** — querying this ontology via `ontology_landscape` /
  `search_ontology` / `genes_by_ontology` requires adding it to the explorer's `ONTOLOGY_CONFIG`
  + the `ontology` enum, tracked in the `multiomics_explorer` repo. This doc is the hand-off artifact.

## See also

- Spec: `docs/superpowers/specs/2026-05-26-psortb-kg-integration-design.md`
- Skill: `.claude/skills/integrate-a-tool/SKILL.md`
- Phase-1 runner: `.claude/skills/psortb-run/SKILL.md`
- Phase-1 design: `docs/superpowers/specs/2026-05-10-psortb-localization-design.md` (its §6 deferred
  Phase-2 sketch proposed a Gene-property approach; this integration supersedes it with the ontology track)
