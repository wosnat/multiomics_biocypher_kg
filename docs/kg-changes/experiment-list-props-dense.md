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

## Semantics (the contract)

| Property | Density | Meaning of `[]` |
|---|---|---|
| `Experiment.treatment_type` | **dense** | characterization experiment — no perturbation whose response is reported (rifampicin decay, TSS mapping). *Not* missing data. |
| `Experiment.background_factors` | **dense and non-empty** | never `[]` — an experiment always has a held-constant context (`axenic`, continuous `light`, `coculture`, …) |
| same two props on `ClusteringAnalysis` / `DerivedMetric` / `MetaboliteAssay` | **dense** | `[]` allowed (e.g. Hackl 2023 genomic-island "clusters" are predicted from sequence, no experiment) |

Contrast with `table_scope` (KG-SYNC-006 ORG-003), which is *sparse* because a
metabolomics-only experiment genuinely has no DE table — the treatment axis is
applicable to every experiment, so its empty case is `[]`, not absence.

## What changed

- `scripts/post-import.sh` + `scripts/post-import.cypher`: new block before
  `Experiment growth_phases` — `SET x.treatment_type = coalesce(x.treatment_type,
  [])`, same for `background_factors`, on Experiment and the three denormalized
  labels.
- `scripts/validate_paperconfig.py`: **error** on empty/missing
  `background_factors`; **error** on empty `treatment_type` when the
  experiment is referenced by any `csv` statistical analysis (DE ⇒
  perturbation); `[]` with no DE analyses prints a "characterization
  experiment" note. Fixed a `TypeError` on an explicit `background_factors:
  null`. New vocabulary value `oxygen`.
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
`treatment_type == []` on the Steglich analysis. Verification:
`MATCH (e:Experiment) RETURN count(e) = count(e.treatment_type) AS dense` → `true`.
