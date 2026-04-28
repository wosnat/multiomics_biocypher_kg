# Non-DE evidence — DerivedMetric extension (Plan 3)

**Merged:** 2026-04-21 (Plan 3 completion); Waldbauer 2012 first real numeric-DM integration landed 2026-04-23.
**Scope:** DerivedMetric nodes + 6 new edge types + `PAIRED_RNASEQ_PROTEOME` extended omics_type. AbundanceAnalysis deferred to a follow-up slice.

## What changed

The knowledge graph now represents non-DE evidence from differential expression publications as **DerivedMetric** nodes, parallel to the existing `ClusteringAnalysis` pattern.

A DerivedMetric captures a column-level scalar summary per gene (periodicity flag, lag coefficient, fold-change category, etc.) that the paper reports **alongside or instead of** raw fold-change.

### New node type

| Node | Count (Biller 2018 + Waldbauer 2012) | Node ID format | Key properties |
|---|---|---|---|
| `DerivedMetric` | 13 (7 Biller 2018 + 6 Waldbauer 2012) | `derived_metric:{doi_short}:{entry_key}:{metric_type}` | `name`, `metric_type`, `value_kind ∈ {numeric, boolean, categorical}`, `rankable`, `has_p_value`, `p_value_threshold`, `unit`, `allowed_categories`, `field_description`, `total_gene_count` (post-import), `growth_phases` (post-import), plus denormalized `experiment_id` / `organism_name` / `publication_doi` / `compartment` / `omics_type` / `treatment` / `light_condition` / `experimental_context` / `treatment_type` / `background_factors` |

### New edge types

Three **binding** edges (1:many parent → DM):

| Edge | Source | Target | Cardinality |
|---|---|---|---|
| `PublicationHasDerivedMetric` | Publication | DerivedMetric | many DMs per Publication |
| `ExperimentHasDerivedMetric` | Experiment | DerivedMetric | many DMs per Experiment |
| `DerivedMetricBelongsToOrganism` | DerivedMetric | OrganismTaxon | one organism per DM |

Three **measurement** edges (DM → Gene). Each DerivedMetric emits exactly ONE edge type, chosen by its `value_kind`:

| Edge | `value_kind` | Edge properties |
|---|---|---|
| `Derived_metric_quantifies_gene` | `numeric` | `metric_type`, `value` (float), `p_value`, `adjusted_p_value`, `rank_by_metric` (post-import, only if parent `rankable="true"`), `metric_percentile` (post-import), `metric_bucket` (post-import: `top_decile`/`top_quartile`/`mid`/`low`), `significant` (post-import, only if parent `has_p_value="true"` and edge has non-null `adjusted_p_value`) |
| `Derived_metric_flags_gene` | `boolean` | `metric_type`, `value_flag ∈ {"true","false"}` |
| `Derived_metric_classifies_gene` | `categorical` | `metric_type`, `value_text` (must be in parent `allowed_categories`) |

### New `Experiment.compartment` property

Every Experiment now has `compartment` (string). Default `"whole_cell"`. Vocab (5 values, defined in `multiomics_kg/vocab/non_de_evidence.py`): `whole_cell`, `vesicle`, `exoproteome`, `spent_medium`, `lysate`.

### New extended omics_type

`PAIRED_RNASEQ_PROTEOME` joins the existing 5 omics_types (`RNASEQ`, `MICROARRAY`, `PROTEOMICS`, `EXOPROTEOMICS`, `METABOLOMICS`) in `EXTENDED_OMICS_TYPES` (`multiomics_kg/vocab/non_de_evidence.py`). Use when a single experiment jointly reports transcript and protein measurements that the paper analyzed as a paired dataset (e.g., Waldbauer 2012 diel cycle, per-gene peak/lag/amplitude across both modalities). DerivedMetric rows attach to one PAIRED_RNASEQ_PROTEOME Experiment instead of two separate RNASEQ + PROTEOMICS Experiments so downstream queries can filter on the cross-modality relationship directly.

## New post-import-computed rollups

Plan 3 post-import Cypher (`scripts/post-import.sh` + `scripts/post-import.cypher`) computes:

### Per DerivedMetric
- `total_gene_count` (int) — count of outgoing measurement edges
- `growth_phases` (str[]) — union from parent Experiment

### Per Experiment
- `reports_fold_change` (str `"true"`/`"false"`) — `"true"` iff has outgoing `Changes_expression_of`
- `reports_derived_metric_types` (str[])
- `derived_metric_count` (int)
- `derived_metric_value_kinds` (str[])
- `derived_metric_gene_count` (int) — distinct genes reachable via any DM edge type

### Per Publication
- `derived_metric_count` (int)
- `derived_metric_gene_count` (int)
- `compartments` (str[]) — from child Experiments
- `derived_metric_types` (str[])
- `derived_metric_value_kinds` (str[])

### Per OrganismTaxon
Same 5 fields as Publication.

### Per Gene (routing signals for MCP dispatch)
- `numeric_metric_count` / `boolean_metric_count` / `categorical_metric_count` (int)
- `numeric_metric_types_observed` / `boolean_metric_types_observed` / `categorical_metric_types_observed` (str[])
- `compartments_observed` (str[])

### Empty-state defaults

All new int properties default to `0`; all new str[] properties default to `[]`. **Never null** — downstream queries can skip null-guards.

## New indexes

- Scalar (8): `derived_metric_metric_type_idx`, `derived_metric_value_kind_idx`, `derived_metric_compartment_idx`, `derived_metric_omics_type_idx`, `derived_metric_treatment_type_idx`, `derived_metric_organism_idx`, `derived_metric_experiment_idx`, `experiment_compartment_idx`
- Full-text (1): `derivedMetricFullText` on `DerivedMetric(name, field_description)`

## Paperconfig surface

New supplementary-materials entry type: `derived_metrics_table`. See `.claude/skills/paperconfig/SKILL.md` for schema.

Biller 2018 (`10.1128/mSystems.00040-18`) was the first real-paper integration: 7 DerivedMetric nodes (6 boolean periodicity flags from Tables S4A/S4B across NATL2A axenic/coculture and MIT1002 coculture; 1 categorical darkness-survival class from Table S5). This retrofits evidence that was previously encoded as `GeneCluster` nodes pre-2026-04-20.

Waldbauer 2012 (`10.1371/journal.pone.0043432`) is the first real numeric-DM integration (landed 2026-04-23): a single `PAIRED_RNASEQ_PROTEOME` Experiment with 6 numeric DerivedMetrics from Table S2 (transcript/protein peak times, lag, diel amplitudes, damping ratio) across 312 cycling genes. No DE comparison — the paper reports diel-cycle summary statistics only.

Edge counts (production graph):
- `Derived_metric_flags_gene`: 4,160 (Biller 2018, 6 boolean DMs)
- `Derived_metric_classifies_gene`: 258 (Biller 2018, 1 categorical DM)
- `Derived_metric_quantifies_gene`: 1,872 (Waldbauer 2012, 6 numeric DMs × 312 genes)

## Example Cypher

```cypher
// Find all genes flagged as periodic in NATL2A axenic L:D (boolean)
MATCH (dm:DerivedMetric {metric_type: 'periodic_in_axenic_LD'})
  -[r:Derived_metric_flags_gene]->(g:Gene)
WHERE r.value_flag = 'true'
  AND g.organism_name = 'Prochlorococcus NATL2A'
RETURN g.locus_tag, g.product;

// Top-decile diel-amplitude genes (numeric, rankable) — Waldbauer 2012
MATCH (dm:DerivedMetric {metric_type: 'diel_amplitude_protein_log2'})
  -[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.metric_bucket = 'top_decile'
RETURN g.organism_name, g.locus_tag, r.value, r.metric_percentile
ORDER BY r.value DESC;

// Classify genes by darkness-survival category (categorical)
MATCH (dm:DerivedMetric {metric_type: 'darkness_survival_class'})
  -[r:Derived_metric_classifies_gene]->(g:Gene)
RETURN r.value_text AS category, count(g) AS gene_count
ORDER BY gene_count DESC;
```

## Out of scope / follow-ups

- **AbundanceAnalysis** node (parent spec §1) — deferred to a later slice. Will cover per-sample abundance (spectral counts, copy-number-normalized transcript counts, community-proteomics fractions). Expected drivers: Biller 2022 vesicles, Oleza 2015/2017 exoproteome.
- **Metabolite layer** (parent spec Appendix) — own spec, own timeline.
- **zinser 2009 Fourier-metric retrofit** — pure paperconfig-authoring task now that the numeric-DM infrastructure is exercised by real data.

## MCP surface

Unchanged in this slice. MCP-team TODO: implement `gene_derived_metrics(locus_tag)` / `derived_metric_ranked_genes(metric_type, bucket)` / `publication_derived_metrics(doi)` tools reading the post-import-computed properties above.
