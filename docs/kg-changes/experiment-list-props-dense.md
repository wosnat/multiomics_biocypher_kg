# Dense `treatment_type` / `background_factors` on experiment-shaped nodes

**Date:** 2026-08-27 · **Trigger:** explorer bug report against KG-SYNC-006
(`gene_clusters_by_gene` → pydantic `treatment_type: Input should be a valid
list, input_value=None` on the Steglich 2010 decay clusters).

## Root cause

The adapters (`omics_adapter`, `cluster_adapter`, `observations_adapter`,
`metabolite_assay_adapter`) already emit `treatment_type = []` /
`background_factors = []` when a paperconfig declares an empty list. BioCypher
writes that as an empty `''` cell in a `string[]` column, and **`neo4j-admin
import` materializes an empty array cell as *no property***. So the three
characterization experiments (Steglich 2010 mRNA half-lives, Voigt 2014 TSS
maps for MED4 and MIT9313) imported with `treatment_type IS NULL` — plus the
denormalized copies on 12 `ClusteringAnalysis` (Steglich decay clusters, Hackl
2023 genomic islands ×11) and 14 `DerivedMetric` nodes, and one Experiment
(Bernstein 2017) with `background_factors IS NULL`.

`growth_phases` never had this problem only because post-import already sets
it unconditionally.

## Semantics (the contract) — revised the same day

| Property | Density | Empty list? |
|---|---|---|
| `Experiment.treatment_type` | **dense, non-empty** (`min_size: 1`) | never. A non-empty list is the explorer's indicator that the node is a real experiment. A study with no perturbation names *what was measured*: `rna_decay` (Steglich 2010), `tss_mapping` (Voigt 2014), `genomic_analysis` (Hackl 2023 islands). Mint a new short categorical value when a new paper fits nothing. |
| `Experiment.background_factors` | **dense, non-empty** (`min_size: 1`) | never — an experiment always has a held-constant context (`axenic`, continuous `light`, `coculture`, …) |
| `ClusteringAnalysis.treatment_type` | dense, non-empty (validator rule on `gene_clusters`) | never (`genomic_analysis` for sequence-only analyses) |
| `ClusteringAnalysis.background_factors` | dense | `[]` allowed — a `genomic_analysis` has no experimental context (Hackl 2023) |
| `DerivedMetric.*`, `MetaboliteAssay.*` | dense (copied from the parent Experiment) | inherits the Experiment rule |

Contrast with `table_scope` (KG-SYNC-006 ORG-003), which is *sparse* because a
metabolomics-only experiment genuinely has no DE table — the treatment axis is
applicable to every experiment.

The first cut of this fix (commit `33772b9b`) allowed `treatment_type = []` on
the three characterization experiments; it was revised the same day because
"non-empty ⇒ real experiment" is a more useful invariant for the explorer than
"`[]` ⇒ characterization".

### Vocabulary mechanics

- The denormalized copies are registered too: `ClusteringAnalysis` / `DerivedMetric` /
  `MetaboliteAssay` × `treatment_type` / `background_factors` (6 closed entries, same value
  lists; `min_size: 1` everywhere except `ClusteringAnalysis.background_factors`), so
  `genomic_analysis` — which only occurs on `ClusteringAnalysis` — is graph-verified.

- `config/controlled_vocabularies.yaml` gains a `min_size` key (string_array
  only; the loader rejects it on scalars). `test_controlled_vocabularies.py::
  test_min_size_lists_are_dense_and_long_enough` asserts it generically on every
  carrier, so no bespoke test is needed per property.
- `scripts/validate_paperconfig.py` loads `CANONICAL_TREATMENT_TYPES` and
  `CANONICAL_BACKGROUND_FACTORS` from the yaml. The two sets are disjoint except
  `light`, `darkness`, `diel`, `coculture`, `chemical`, `viral` (meaning depends
  on the field). `mutant` was validator-only and unused — dropped.

## What changed

- `scripts/post-import.sh` + `scripts/post-import.cypher`: new block before
  `Experiment growth_phases` — `SET x.treatment_type = coalesce(x.treatment_type,
  [])`, same for `background_factors`, on Experiment and the three denormalized
  labels.
- `scripts/validate_paperconfig.py`: **error** on empty/missing
  `background_factors` and on empty `treatment_type` (experiments and
  `gene_clusters`); per-field canonical sets loaded from the yaml. Fixed a
  `TypeError` on an explicit `background_factors: null`. New vocabulary
  values `oxygen`, `rna_decay`, `tss_mapping`, `genomic_analysis`.
- `config/controlled_vocabularies.yaml`: `Experiment.treatment_type` gains
  `oxygen`; both descriptions document the dense rule. (Hash drift until the
  next build — expected.)
- `tests/kg_validity/test_expression.py`: `test_experiment_list_props_dense`,
  `test_experiment_background_factors_non_empty`,
  `test_denormalized_experiment_list_props_dense[ClusteringAnalysis|DerivedMetric|MetaboliteAssay]`;
  `oxygen` added to the canonical set.
- `tests/test_paperconfig_validation.py`: five validator tests for the new
  rules; fixtures now declare `background_factors`.
- Bernstein 2017 paperconfig relabelled — see CHANGELOG `### Data`.

## Explorer side

No code change. `coalesce(e.treatment_type, [])` is not needed at any
projection site; the edge-case scenario becomes an assertion that
`treatment_type == ['rna_decay']` on the Steglich analysis. Add `oxygen`,
`rna_decay`, `tss_mapping`, `genomic_analysis` to any hard-coded
treatment-type enum (or read `ControlledVocabulary`). Verification:
`MATCH (e:Experiment) RETURN count(e) = count(e.treatment_type) AS dense` → `true`.
