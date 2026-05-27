# SignalP signal-peptide type (`SignalPeptideType` nodes)

**Date:** 2026-05-27
**Spec:** [`docs/superpowers/specs/2026-05-27-signalp-kg-integration-design.md`](../superpowers/specs/2026-05-27-signalp-kg-integration-design.md)
**Track:** 3A ontology · flat · scored · **structural** (how the gene's product is handled, not what it does)
**Skill:** produced with [`/integrate-a-tool`](../../.claude/skills/integrate-a-tool/SKILL.md) (Phase-2 partner to `/add-a-tool`)

## What's changing

SignalP 6.0 predicts the signal-peptide type of each protein (Sec/Tat translocon × SPase I/II/III).
Phase 1 (`/signalp-run`) ran SignalP on 40 strains, writing raw `signalp/prediction_results.txt`.
Because `signalp-run` predates the `add-a-tool` calls.json convention, a new `--normalize` mode
parses that raw output (no Docker re-run) into committed `<strain>.signalp.calls.json` artifacts;
this change then wires them into the KG. The type call + winning-class probability + cleavage site
flow through the gene-annotation merge (`signalp` source → `signalp_type` / `signalp_probability` /
`signalp_cleavage_site` / `signalp_cleavage_probability` fields in `gene_annotations_merged.json`,
joined on the RefSeq WP_ `protein_id`), then surface as a small **flat ontology**: the five real
SignalP types become `SignalPeptideType` nodes, and each gene with a confident call gets one
`Gene_has_signal_peptide_type` edge carrying its probability + cleavage info.

This mirrors the PSORTb integration ([`psortb-extension.md`](psortb-extension.md)) end-to-end —
flat, scored, 1:1, structural.

## New node label: `SignalPeptideType`

- **Count:** 5 nodes (all real types observed across 40 organisms). Flat — single level.
  Post-import `gene_count`: SP 9,710 · LIPO 3,065 · TAT 441 · PILIN 314 · TATLIPO 83.
- **ID prefix `signalp_`** (underscore, not colon: `signalp` is not a registered bioregistry prefix,
  so `normalize_curie('signalp:X')` returns `None` and the id falls back to `signalp_SP`).
- **Properties (adapter-emitted):** `name`, `signalp_id`, `level` (always `0` — flat, no `level_kind`).
- **Computed (post-import):** `gene_count`, `organism_count`.
- **Vocabulary** (the `OTHER` sentinel is skipped — absence of an edge encodes "no signal peptide"):

  | signalp_id | name | translocon / SPase |
  |---|---|---|
  | `SP` | Signal peptide (Sec/SPI) | Sec, SPase I |
  | `LIPO` | Lipoprotein signal peptide (Sec/SPII) | Sec, SPase II |
  | `TAT` | TAT signal peptide (Tat/SPI) | Tat, SPase I |
  | `TATLIPO` | TAT lipoprotein signal peptide (Tat/SPII) | Tat, SPase II |
  | `PILIN` | Pilin-like signal peptide (Sec/SPIII) | Sec, SPase III |

## New edge type

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_signal_peptide_type` | Gene → SignalPeptideType | 13,613 | `probability: float` (winning-class likelihood, ∈[0,1]); `cleavage_site: int` (residue before the cut; absent when no CS reported); `cleavage_probability: float` (absent when no CS); post-import `rank_by_probability: int` |

The KG's **second scored ontology edge** (after PSORTb) — the score rides on the edge, modeled on
the `Changes_expression_of` property block. 1:1 (≤1 edge per gene); no hierarchy (`*_is_a_*`) edges.
The full 6-class likelihood vector is intentionally **not** stored (only the winning class's
probability). `cleavage_site` / `cleavage_probability` are omitted (absent) when SignalP reports no
cleavage site.

## Cross-node rollups

- **Structural, not functional:** deliberately **NOT** folded into `Gene.annotation_types` /
  `informative_annotation_types` / `annotation_quality` (signal-peptide presence is *how/where* a
  gene's product is handled, not *what it does* — folding it in would inflate functional-annotation
  coverage and skew `genes_by_function` `min_quality`).
- **`Gene.signal_peptide_type`** (post-import, str): denormalized 1:1 routing string = the gene's
  single SignalP type (absent when no confident call). Surfaces via `gene_details` for free.
- `OrganismTaxon` / `Publication`: no per-ontology rollup (unchanged).

## Indexes

Scalar: `signal_peptide_type_level_idx`, `signal_peptide_type_id_idx`.
Full-text: `signalPeptideTypeFullText` on `name`, `signalp_id`.

## Merge wiring + provenance

- `config/gene_annotations_config.yaml`: `signalp` source (join `protein_id`) + `signalp_type` /
  `signalp_probability` / `signalp_cleavage_site` / `signalp_cleavage_probability` fields + a
  `signalp` `logical_source` (`scope: gene_level`, `provenance: tool_run`).
- `multiomics_kg/download/build_gene_annotations.py`: `load_signalp()` + the source threaded through
  `_get_raw` / all six resolvers / `build_wide` / `build_merged` / the `process_strain` row join /
  `_compute_contributing_sources`. (No new `annotation_transforms.py` transform — normalizing to
  calls.json first makes the merge a clean PSORTb mirror; the decision-tree's `_tx_parse_cs_pos` is
  obviated because CS parsing happens in the normalizer's pure parser, `multiomics_kg/utils/signalp.py`.)
- A `DataSource` node `data_source:signalp` emits automatically (named via
  `data_source_adapter.py`'s `_name_for`/`_description_for`); `Gene.contributing_sources` lists `signalp`.

## Validation results (live rebuild, 2026-05-27)

- **Build:** create_knowledge_graph.py wrote 5 `SignalPeptideType` nodes + **13,613**
  `Gene_has_signal_peptide_type` edges; full `docker compose up -d` chain exited 0.
- **Live counts:** node `gene_count` SP 9,710 · LIPO 3,065 · TAT 441 · PILIN 314 · TATLIPO 83;
  `organism_count` 30–40; every edge carries `probability` (range 0.30–1.0), `cleavage_site`,
  `rank_by_probability`; every source gene carries the `signal_peptide_type` routing string.
- `/omics-edge-snapshot` before/after: expression (`Changes_expression_of`) edge counts, all
  DerivedMetric edges, and all metabolism nodes/edges **unchanged** (no regressions, exit 0).
- `pytest -m kg`: **983 passed, 4 skipped** — incl. `tests/kg_validity/test_signalp.py` (node set ⊆
  the 5 types, no `OTHER`, flat level 0, `probability` ∈ [0,1], `cleavage_site` a positive int when
  present, `rank_by_probability` present, no orphan edges, routing-string consistency, signalp kept
  out of `annotation_types`) and the DataSource node-count test updated 5 → 6.
- `pytest -m "not slow and not kg"`: **1966 passed** (incl. `tests/test_signalp.py`, 20 unit tests).
- `snapshot_data.json` regenerated; `tests/kg_validity/test_snapshot.py` green.

## What this enables

- "All secreted (signal-peptide-bearing) genes in organism X" — one hop:
  `MATCH (g:Gene)-[:Gene_has_signal_peptide_type]->(:SignalPeptideType)`.
- Type-specific shortlists (e.g. all lipoproteins): `{signalp_id:'LIPO'}`.
- Probability-ranked shortlists per type via `rank_by_probability`.
- A `Gene.signal_peptide_type` routing string for gene-overview / detail surfaces.

## What does NOT change

- Phase 1 (`/signalp-run`) Docker run untouched (only a `--normalize` mode added). No existing
  node/edge/property modified (additions-only).
- **MCP / explorer integration is out of scope** — querying this ontology via `ontology_landscape` /
  `search_ontology` / `genes_by_ontology` requires adding it to the explorer's `ONTOLOGY_CONFIG`
  + the `ontology` enum, tracked in the `multiomics_explorer` repo. This doc is the hand-off artifact.

## See also

- Spec: `docs/superpowers/specs/2026-05-27-signalp-kg-integration-design.md`
- Skill: `.claude/skills/integrate-a-tool/SKILL.md`
- Phase-1 runner: `.claude/skills/signalp-run/SKILL.md` (`--normalize` mode + calls.json schema)
- Sibling integration: `docs/kg-changes/psortb-extension.md` (the flat-scored-structural template)
