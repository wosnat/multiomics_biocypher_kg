# Non-DE Evidence Extension — Design

**Date:** 2026-04-17
**Status:** Draft (awaiting user review)
**Scope:** Schema extension for paper-grounded quantitative and qualitative gene-level data that doesn't fit the current `Changes_expression_of` (DE) shape.

## Problem

The KG currently represents paper evidence primarily through `Changes_expression_of` edges (Experiment → Gene, with `log2_fold_change` + `adjusted_p_value`). This works for papers that publish a differential-expression table comparing treatment vs control.

A survey of the 45 papers under `data/Prochlorococcus/papers_and_supp/` and `data/Synechococcus/papers_and_supp/` identified categories of published evidence that do not fit the DE shape:

- **Per-sample abundance without a DE comparison** (10 papers): NSAF / spectral counts / transcripts-per-cell — quantitative values that the paper did *not* subject to a treatment-vs-control statistical test (e.g., biller 2022 vesicle proteome, Oleza 2015 exoproteome).
- **Paper-computed per-gene derived metrics** (5+ papers): diel phase, Fourier periodicity score, protein-transcript lag, damping ratio — not abundance values at all, but statistics the paper calculates per gene from an underlying time course.

**Out of scope, deferred to separate spec:** Metabolite / lipid / non-gene-entity representation (Kujawinski 2023, biller 2022 metabolite tables, Capovilla 2023 metabolomics, etc.). See the **Appendix — Deferred: Metabolite layer** at the end of this document.

**Gene-only targets.** All new edges point at `Gene`. Proteomics resolves protein IDs to encoding genes at adapter time (same pattern as existing `Changes_expression_of` for proteomics DE). Protein nodes are often missing per `CLAUDE.md`'s orphan-proteins issue, so Gene is the reliable hub.

## Design principle

**The KG contains paper-grounded evidence; analysis-level interpretation lives outside the KG.**

- No cross-paper *interpretation* on Gene nodes (signature gene sets, stress classifications, "stress response consistency" properties). MCP / Python compute these on demand.
- Cross-paper *routing counts* on Gene nodes are allowed and already in use (`expression_edge_count`, `cluster_membership_count`, `closest_ortholog_genera`); this spec adds parallel non-DE counts so the routing tool isn't blind to half the new evidence (see §5).
- No pathway enrichment, matching-across-evidence, or interpretation pre-computed in the graph.
- Schema extensions exist to make per-paper evidence *interpretable enough* that an LLM/MCP can combine it; combining itself is query-time.

**Gene nodes are the aggregation hubs.** Every evidence edge points at a Gene. Cross-paper combining happens by traversing edges off a gene, not by reading precomputed summaries.

**Rank is the cross-paper comparability device.** Raw values (NSAF, spectral count, TPM, Fourier score) are not meaningfully comparable across studies that use different methods. Rank within the analysis scope is. Mirrors existing `rank_by_effect` / `rank_up` / `rank_down` on DE edges.

**Analysis nodes carry column-level metadata.** Following the cluster pattern (`Publication → ClusteringAnalysis → GeneCluster → Gene`), new non-DE evidence uses intermediate analysis nodes:
- `Publication → Experiment → AbundanceAnalysis → Gene` for per-sample abundance
- `Publication → Experiment → DerivedMetric → Gene` for paper-computed per-gene metrics

Column-level metadata (value_type, unit, rankable, timepoint) lives on the analysis node; per-gene measurement data (value, rank) lives on the edge. This prevents duplicating column description across ~N_genes edges.

## Schema changes

### §1. Abundance evidence

**New node: `AbundanceAnalysis`** — one per (Experiment × time_point). Carries the "column description" metadata shared across all edges from this analysis.

```
(:AbundanceAnalysis {
  id: str,                        // unique; "{doi}:{paperconfig_key}" for single-timepoint,
                                  //         "{doi}:{paperconfig_key}:{time_point_order}" for timeseries
  name: str,                      // human-readable; e.g., "MIT9312 HL vesicle proteome" (parallels ClusteringAnalysis.name)
  experiment_id: str,             // explicit parent Experiment ID; avoids string-parsing the analysis ID
  organism_name: str,             // denormalized from Experiment for MCP query parallelism
  // Filter symmetry — denormalized from parent Experiment so list_abundance_analyses(...) filters
  // collapse to 1-hop. Mirrors ClusteringAnalysis pattern. Validator + post-import invariant
  // assert these match the parent Experiment.
  publication_doi: str,
  compartment: str,               // controlled vocabulary (see §3); inherited from Experiment
  omics_type: str,
  treatment_type: str[],
  background_factors: str[],
  treatment: str | null,          // denormalized from Experiment
  light_condition: str | null,    // denormalized from Experiment
  experimental_context: str | null, // denormalized from Experiment
  // Column-level metadata
  value_type: str,                // controlled vocabulary (see below)
  time_point: str | null,         // label for timeseries; null for single-timepoint analyses
  time_point_order: int | null,
  time_point_hours: float | null, // no sentinel values — null when unknown
  growth_phase: str | null,       // physiological state at this timepoint; parallels Changes_expression_of.growth_phase
  field_description: str | null,  // describes what value represents (e.g., "biovolume-normalized spectral counts"); orthogonal to Experiment.experimental_context
  // Post-import computed (mirror ClusteringAnalysis):
  total_gene_count: int | null,   // count of outgoing measurement edges; parallels ClusteringAnalysis.total_gene_count
  growth_phases: str[] | null     // union of parent Experiment's growth_phases (see post-import block #3)
})
```

**New edges.** Three binding edges mirror the cluster pattern (Publication/Experiment/Organism → ClusteringAnalysis), plus one measurement edge to Gene:

```
(Publication)-[:publication_has_abundance_analysis]->(AbundanceAnalysis)   // mirrors Publication_has_clustering_analysis
(Experiment)-[:experiment_has_abundance_analysis]->(AbundanceAnalysis)      // mirrors Experiment_has_clustering_analysis
(AbundanceAnalysis)-[:abundance_analysis_belongs_to_organism]->(OrganismTaxon)  // mirrors Clusteringanalysis_belongs_to_organism

(AbundanceAnalysis)-[:abundance_analysis_measures_gene {
  value: float,                      // central tendency (mean by convention; median flagged in AbundanceAnalysis.field_description if paper uses it), in value_type units. Non-null — abundance edges always carry a measurement. Detection-only tables use derived_metric_flags_gene with metric_type=detected instead; see §2.
  value_sd: float | null,            // standard deviation across replicates; null when n_replicates < 2 or paper doesn't report; adapter converts SEM → SD via SEM·√n when n is known
  n_replicates: int | null,          // biological replicate count folded into value; null if unknown
  value_type: str,                   // denormalized from parent AbundanceAnalysis; collapses cross-paper value-type queries to 1-hop
  // Post-import (grouped by AbundanceAnalysis — every edge has non-null rank/percentile/bucket):
  rank_in_sample: int,               // 1 = most abundant
  abundance_percentile: float,       // 0-100
  abundance_bucket: str              // discrete enum from abundance_percentile. Vocabulary + thresholds: top_decile (≥90), top_quartile (≥75 and <90), mid (≥25 and <75), low (<25). Materialized for index-friendly cross-paper queries; mirrors expression_status on DE edges.
}]->(Gene)
```

**Cardinality:** one `AbundanceAnalysis` per (Experiment × time_point). Non-timeseries experiments have one AbundanceAnalysis. Timeseries experiments have N analyses, one per timepoint. All timepoints of a given series share the same parent Experiment; the Experiment's `is_time_course` flag and `time_point_*` arrays (existing post-import aggregations) remain the timeseries container. See §3 for how abundance-side timepoints are surfaced on the Experiment.

**Edge-direction convention.** `abundance_analysis_measures_gene` and `derived_metric_quantifies_gene` point *from* the analysis *to* the gene — the analysis node is the source. This mirrors the existing `Changes_expression_of` (Experiment → Gene) convention: **measurement edges** put the analytic source upstream and the measured target (Gene) downstream. It is the reverse of `Gene_in_gene_cluster` (gene-as-source) — but cluster edges express categorical *membership*, not measurement, so gene-as-source is natural there. Rule of thumb: edges that carry numeric per-gene properties (value, rank, fold change) go *analysis → gene*; edges that express "this gene belongs to X" go *gene → X*.

**Target is `Gene` only.** Proteomics papers resolve protein IDs to encoding genes at ingest.

**`abundance_bucket` is materialized** as a string enum on every `abundance_analysis_measures_gene` edge (see edge spec above). Mirrors `expression_status` on DE edges; gives queries an index-friendly handle without a CASE per call site. Thresholds are pinned in the spec to keep post-import iterations from drifting.

**Detection-only tables (papers listing proteins present without quantification)** are NOT modeled as `abundance_table` with a special flag. Use `derived_metrics_table` with `value_kind=boolean` and the single generic `metric_type=detected`. Compartment context comes from the parent Experiment (denormalized onto DerivedMetric, plus `Gene.compartments_observed` rolls up across all detection edges per gene). Compartment-specific variants (`detected_in_vesicle`, etc.) are NOT used — they would balloon the vocabulary table without adding query power. Positive-only lists use `blank_policy: skip` and no `false_tokens`. This keeps the abundance path clean — every abundance edge carries a real numeric value — and uses the boolean primitive we already have for Y/N per-gene evidence.

**Compartment:** lives on Experiment (see §3) and is denormalized onto every child AbundanceAnalysis. Every AbundanceAnalysis inherits from its parent Experiment.

**Use when:** paper reports abundance/detection but no treatment-vs-control fold change (vesicle proteome, exoproteome, per-condition detection lists). Do *not* use as a universal retrofit onto time-series papers that already have DE edges.

### §2. Derived-metrics evidence

**New node: `DerivedMetric`** — one per (Experiment × metric_type). Carries column-level metric semantics.

```
(:DerivedMetric {
  id: str,                        // unique; "{doi}:{paperconfig_key}:{metric_type}"
  name: str,                      // human-readable; e.g., "Fourier periodicity score" (parallels GeneCluster.name)
  experiment_id: str,             // explicit parent Experiment ID; avoids string-parsing the analysis ID
  organism_name: str,             // denormalized from Experiment for MCP query parallelism
  // Filter symmetry — denormalized from parent Experiment so list_derived_metrics(...) filters
  // collapse to 1-hop. Validator + post-import invariant assert these match the parent Experiment.
  publication_doi: str,
  compartment: str,
  omics_type: str,
  treatment_type: str[],
  background_factors: str[],
  treatment: str | null,          // denormalized from Experiment
  light_condition: str | null,    // denormalized from Experiment
  experimental_context: str | null, // denormalized from Experiment
  // Column-level metric semantics
  metric_type: str,               // controlled vocabulary (see below)
  value_kind: str,                // "numeric" | "boolean" | "categorical". Determines which edge type the analysis emits. Fixed per metric_type (not per-paper).
  unit: str | null,               // e.g., "h", "dimensionless". Required for value_kind=numeric; null for boolean/categorical.
  rankable: str,                  // "true" | "false" — string enum (not bool; BioCypher booleans unreliable). Governs post-import rank_by_metric/metric_percentile/metric_bucket. Always "false" when value_kind != "numeric".
  has_p_value: str,               // "true" | "false" — string enum. "true" iff outgoing edges carry p_value / adjusted_p_value. Always "false" when value_kind != "numeric".
  p_value_threshold: float | null,// significance cutoff used by the paper; null when no p-value or no threshold reported. Used by post-import to derive `significant` on edges.
  allowed_categories: str[] | null, // for value_kind=categorical: enumerated set of permitted value_text values (e.g., ["present_in_axenic","unique_to_coculture",...]). Validator + adapter reject edges with out-of-set values. Null for numeric/boolean.
  field_description: str | null,  // describes what value represents (e.g., "Fourier periodicity score; higher = more cyclic"); orthogonal to Experiment.experimental_context
  // Post-import computed (mirror ClusteringAnalysis):
  total_gene_count: int | null,   // count of outgoing measurement edges; parallels ClusteringAnalysis.total_gene_count
  growth_phases: str[] | null     // union of parent Experiment's growth_phases (see post-import block #3)
})
```

**New edges.** Three binding edges mirror the cluster pattern, plus three sibling measurement edges to Gene:

```
(Publication)-[:publication_has_derived_metric]->(DerivedMetric)   // mirrors Publication_has_clustering_analysis
(Experiment)-[:experiment_has_derived_metric]->(DerivedMetric)      // mirrors Experiment_has_clustering_analysis
(DerivedMetric)-[:derived_metric_belongs_to_organism]->(OrganismTaxon)  // mirrors Clusteringanalysis_belongs_to_organism
```

Three sibling edges from `DerivedMetric` to `Gene`, one per `value_kind`. Each DerivedMetric emits exactly one edge type for all its outgoing measurements (driven by its `value_kind`). No `value_kind` discriminator on the edges themselves — the edge type IS the discriminator.

```
// value_kind=numeric: paper-computed scalar measurement
(DerivedMetric)-[:derived_metric_quantifies_gene {
  metric_type: str,               // denormalized from parent DerivedMetric
  value: float,                   // paper-reported numeric value (non-null)
  p_value: float | null,          // raw significance (e.g., Fisher's g-test for fourier_score); null when DerivedMetric.has_p_value="false" OR per-row missing
  adjusted_p_value: float | null, // multiple-testing-corrected p-value when paper reports it; null otherwise. Adapter does NOT re-correct.
  // Post-import (only when DerivedMetric.rankable="true"):
  rank_by_metric: int | null,     // 1 = highest value
  metric_percentile: float | null,// 0-100; gene's percentile within this DerivedMetric by value. Parallels abundance_percentile. Null when rankable="false".
  metric_bucket: str | null,      // discrete enum from metric_percentile; null when rankable="false". Vocabulary + thresholds match abundance_bucket: top_decile (≥90), top_quartile (≥75 and <90), mid (≥25 and <75), low (<25). Materialized for index-friendly cross-paper queries like "top-decile periodic genes."
  // Post-import (only when DerivedMetric.has_p_value="true"):
  significant: str | null         // "true" | "false" — string enum (matches existing Changes_expression_of.significant convention). Computed as (adjusted_p_value < DerivedMetric.p_value_threshold). Null when has_p_value="false" OR threshold null OR adjusted_p_value null.
}]->(Gene)

// value_kind=boolean: per-gene Y/N flag from a per-gene test (e.g., RAIN periodicity FDR<0.05)
(DerivedMetric)-[:derived_metric_flags_gene {
  metric_type: str,               // denormalized from parent DerivedMetric
  value_flag: str                 // "true" | "false" — string enum (not bool; BioCypher boolean edge properties unreliable). Non-null.
                                  // Positive-only mode: when paperconfig blank_policy="skip" and false_tokens is empty,
                                  // every emitted edge carries value_flag="true" — the property is constant across
                                  // edges of that DerivedMetric. Kept for schema uniformity (queries always use
                                  // `WHERE r.value_flag = 'true'` regardless of mode), and so future mixed-mode
                                  // tables ingested under the same metric_type don't require a schema change.
}]->(Gene)

// value_kind=categorical: paper-defined mutually-exclusive class assignment
(DerivedMetric)-[:derived_metric_classifies_gene {
  metric_type: str,               // denormalized from parent DerivedMetric
  value_text: str                 // class label (non-null); must be in DerivedMetric.allowed_categories
}]->(Gene)
```

**Cardinality:** one `DerivedMetric` per (Experiment × metric_type). Experiments with multiple metrics (zinser 2009: Fourier + peak_time + peak_R²; Biller 2018 S4A: 4 boolean periodicity flags) have multiple DerivedMetric children.

**Edge-direction convention.** All three sibling measurement edges (`derived_metric_quantifies_gene`, `derived_metric_flags_gene`, `derived_metric_classifies_gene`) point *from* DerivedMetric *to* Gene — same convention as `abundance_analysis_measures_gene` and `Changes_expression_of`. See §1 "Edge-direction convention" for the full rule (measurement edges → analysis-as-source; membership edges → gene-as-source).

**Target is `Gene` only.**

**Rankability lives on the DerivedMetric node**, not a per-edge decision. Post-import populates `rank_by_metric` only when the parent DerivedMetric has `value_kind='numeric'` AND `rankable="true"`.

**Three value kinds supported, three sibling edge types.** `DerivedMetric.value_kind` selects which edge type the analysis emits:
- `numeric` → `derived_metric_quantifies_gene` (paper-computed scalar like Fourier score)
- `boolean` → `derived_metric_flags_gene` (per-gene Y/N flag from a per-gene test like RAIN periodicity FDR<0.05)
- `categorical` → `derived_metric_classifies_gene` (mutually-exclusive class assignment like darkness-survival category)

Each metric_type has a fixed value_kind in the vocabulary — a future numeric variant of an existing categorical metric gets a new metric_type rather than flipping value_kind. A given DerivedMetric emits exactly one edge type for all its outgoing measurements; mixed-kind metrics are not allowed (split into multiple DerivedMetric nodes).

**Boolean semantics.** Real paper tables use inconsistent tokens (`Y` / `N` / blank / `NA` / `n/a` / `#N/A` / `?`). Paperconfig declares explicit token lists; adapter hard-errors on unexpected tokens.

Paperconfig fields on a boolean metric entry:
- `true_tokens: str[]` (required) — literal cell values that map to `value_flag="true"`. E.g., `["Y"]`, or `["Y", "yes", "1"]` if the column mixes conventions. Case-sensitive literal match — no case folding.
- `false_tokens: str[]` (optional; default `[]`) — literal cell values that map to `value_flag="false"`. When empty, no cell produces `value_flag="false"` via explicit token match.
- `skip_tokens: str[]` (optional; default `["NA", "N/A", "n/a", "#N/A"]`) — literal cell values meaning "not tested" → no edge emitted.
- `blank_policy: str` (optional; default `"skip"`) — one of `"skip"` (blank → no edge; paper listed only positive cases means absent = not tested), `"false"` (blank → `value_flag="false"`), `"true"` (blank → `value_flag="true"`, rare).

Adapter processing per cell:
1. If cell value ∈ `true_tokens` → emit edge with `value_flag="true"`.
2. Else if cell value ∈ `false_tokens` → emit edge with `value_flag="false"`.
3. Else if cell value ∈ `skip_tokens` → no edge (gene = not tested).
4. Else if cell value is empty/whitespace → apply `blank_policy`.
5. Else → **hard error at ingest** (unexpected token; paperconfig author must explicitly classify it). Prevents silent miscategorization.

"Gene not in CSV at all" always means "not tested" — no edge is emitted, regardless of the token policy.

**Categorical semantics.**
- Source CSV column carries the category string per row.
- `DerivedMetric.allowed_categories` enumerates permitted values; adapter rejects rows whose value isn't in the set.
- One edge per gene per metric (mutually exclusive).

**Use when:** paper reports a per-gene statistic, flag, or classification as part of its own analysis — diel phase, Fourier score, protein-transcript lag, RAIN periodicity Y/N, darkness-survival category, half-life, correlation coefficient. Use `gene_clusters` only when the paper presents results as named gene groups with their own descriptions/labels (e.g., k-means clusters with functional descriptions per cluster), not when the "groups" are just composite labels of independent boolean flags.

### §3. Experiment node extensions

Eight new scalar/array properties (one adapter-emitted, seven post-import-computed) plus four parallel time-point arrays:

**Adapter-emitted:**
- `compartment: str` — controlled vocabulary. Default `"whole_cell"` for experiments that don't specify. Adapter-emitted (so Experiment-level filtering works pre-post-import); also denormalized onto child analysis nodes.

**Post-import-computed** (see post-import block #4):
- `reports_fold_change: str` — `"true"` | `"false"` (string enum, not bool — same convention as `Changes_expression_of.significant`). Set to `"true"` iff Experiment has outgoing `Changes_expression_of` edges.
- `reports_derived_metric_types: str[]` — distinct `metric_type` values reachable via `experiment_has_derived_metric` → DerivedMetric. Empty if none. (Renamed from `reports_calculated_metrics` for vocabulary alignment with `DerivedMetric.metric_type`.)
- `abundance_analysis_count: int` — count of child `AbundanceAnalysis` nodes (via `experiment_has_abundance_analysis`). Parallels existing `clustering_analysis_count`.
- `derived_metric_count: int` — count of child `DerivedMetric` nodes (via `experiment_has_derived_metric`).
- `derived_metric_value_kinds: str[]` — sorted distinct `value_kind` values across child DerivedMetrics (subset of `{"numeric","boolean","categorical"}`).
- `abundance_gene_count: int` — distinct genes reachable via any child `AbundanceAnalysis`. Parallel to existing `gene_count` (DE-only); ensures downstream summary queries don't silently report `gene_count: 0` for non-DE experiments (e.g., vesicle proteome) that actually have measurement coverage.
- `derived_metric_gene_count: int` — distinct genes reachable via any child `DerivedMetric` (across all three derived-metric edge types).

**Intra-experiment rollups**, not cross-paper rollups. An Experiment summarizes its own one-hop neighbourhood (its child analyses and DE edges) — not aggregating evidence across papers for a gene. The "no cross-paper interpretation" principle (see Design principle) applies to gene-hub interpretations (`Gene.stress_response_consistency`, etc.), not to intra-experiment filters or routing counts.

**Time-point arrays — DE-only today, extended with parallel abundance arrays.**

Existing `Experiment.time_point_labels / time_point_orders / time_point_hours / time_point_totals / time_point_significant_up / time_point_significant_down` aggregate over `Changes_expression_of` edges. Their semantics (per-tp significant counts) are DE-specific.

For timeseries `AbundanceAnalysis` children (e.g., a future Waldbauer per-tp transcript ingestion), a non-DE experiment would otherwise have empty time_point arrays — undermining the §3 honesty fix (`abundance_gene_count` non-zero, `time_point_labels=[]`).

Resolution:
- **Existing DE-centric arrays stay DE-scoped** — breaking their semantics (e.g., redefining `time_point_totals` as a DE+abundance union) would silently change MCP tool outputs that currently rely on DE-only counts. Non-breaking.
- **New parallel abundance arrays**, post-import-computed from `AbundanceAnalysis` children (existing arrays ordered by `time_point_order`, parallel semantics):
  - `abundance_time_point_labels: str[]`
  - `abundance_time_point_orders: int[]`
  - `abundance_time_point_hours: float[]` (uses `-1.0` sentinel per existing array convention)
  - `abundance_time_point_totals: int[]` — per-tp gene count (outgoing `abundance_analysis_measures_gene` edges per AbundanceAnalysis)
- `is_time_course` remains a single boolean that flips true if *either* DE or abundance children are multi-timepoint.

DerivedMetric children are scalar summaries (no time_point) and do not contribute to time_point arrays.

Downstream consequence: any tool or query that surfaces per-experiment timepoint summaries can present DE time-point arrays alongside abundance time-point arrays without conflation. Non-DE timeseries experiments report a non-empty `abundance_time_point_labels` instead of a misleading empty DE array. (How that surfaces in MCP envelopes is downstream — see "What's NOT in this spec".)

Validator invariant: for timeseries experiments, the union of DE timepoints and abundance timepoints must form a single ordered sequence (no two different timepoint labels sharing the same `time_point_order`). Enforces invariant #3.

### §4. `Changes_expression_of` — unchanged

No retrofit in this spec. Existing DE edges stay flat: `(Experiment)-[:Changes_expression_of]->(Gene)` with edge-level `time_point`, `log2_fold_change`, `adjusted_p_value`, `rank_by_effect`, etc. All existing MCP queries continue to work unchanged.

**Known asymmetry after this spec lands:**
- DE evidence: `Experiment → Gene` (flat; edge carries timepoint + stats)
- Non-DE evidence: `Experiment → Analysis → Gene` (two-hop; analysis carries timepoint/column-metadata)

The asymmetry is documented transition state. A follow-up spec may retrofit DE to the analysis-node pattern once the new shape has been validated in practice. Listed under "What's NOT in this spec."

### §5. Gene routing fields for non-DE evidence

The KG already carries cross-paper *routing counts* on Gene nodes (`expression_edge_count`, `significant_up_count`, `significant_down_count`, `cluster_membership_count`, `cluster_types`, `closest_ortholog_genera`, etc. — set by post-import). The MCP `gene_overview` routing tool reads these to decide which detail tools to call.

Without parallel non-DE counts, that routing tool is blind to half the new evidence — `gene_overview` would not know whether to recommend `abundance_by_gene` or `derived_metrics_by_gene` for any given gene.

**New post-import-derived properties on Gene** (split per derived-metric edge type so the routing tool can pick the right MCP tool):

- `abundance_analysis_count: int` — distinct `AbundanceAnalysis` nodes that quantify this gene (incoming `abundance_analysis_measures_gene` edges).
- `numeric_metric_count: int` — distinct `DerivedMetric` nodes (value_kind=numeric) reaching this gene via `derived_metric_quantifies_gene`.
- `classifier_flag_count: int` — distinct `DerivedMetric` nodes (value_kind=boolean) reaching this gene via `derived_metric_flags_gene`.
- `classifier_label_count: int` — distinct `DerivedMetric` nodes (value_kind=categorical) reaching this gene via `derived_metric_classifies_gene`.
- `compartments_observed: str[]` — sorted distinct `compartment` values across all incoming non-DE evidence (parallels `cluster_types`).
- `numeric_metric_types_observed: str[]` — sorted distinct `metric_type` values from incoming `derived_metric_quantifies_gene` edges.
- `classifier_flag_types_observed: str[]` — sorted distinct `metric_type` values from incoming `derived_metric_flags_gene` edges.
- `classifier_label_types_observed: str[]` — sorted distinct `metric_type` values from incoming `derived_metric_classifies_gene` edges.

These are routing counts, not interpretations. They mirror the existing pattern. Cross-paper *interpretation* (signature gene sets, stress classifications, stress-response consistency) remains out of scope and is computed on demand by MCP/Python.

## Controlled vocabularies (grounded in actual paper content)

### `compartment`

| Value | Meaning | Source papers |
|---|---|---|
| `whole_cell` | Intracellular (default) | Most DE papers |
| `vesicle` | Extracellular vesicle fraction | biller 2022 |
| `exoproteome` | Secreted proteins in medium | Oleza 2015, Oleza 2017, kaur 2018 |
| `spent_medium` | Supernatant from culture | Thompson 2016 |
| `lysate` | Cell lysate | Thompson 2016 |

New values added only when a concrete paper requires one.

### `value_type` (on `AbundanceAnalysis`)

| Value | Source papers |
|---|---|
| `tpm` | alonso 2023, Capovilla 2023 |
| `rpkm` | Thompson 2016 |
| `transcript_count_norm` | Capovilla 2023 |
| `transcripts_per_cell` | alonso 2023 |
| `transcripts_per_biovolume` | alonso 2023 |
| `rma_normalized` | zinser 2009 (microarray) |
| `spectral_count` | Oleza 2017, biller 2022 S2 |
| `spectral_count_biovolume_norm` | biller 2022 S3 |
| `nsaf` | Oleza 2015, Oleza 2017 |

Each `value_type` has documented semantics (what it measures, whether cross-sample comparison is meaningful, units). Documented in `docs/value_types.md` (to be written alongside adapter implementation).

### `metric_type` (on `DerivedMetric`)

The `value_kind` column selects which edge-value field is populated. `rankable` and `has_p_value` apply only to `value_kind=numeric`; both are forced to false when value_kind is boolean or categorical (significance for boolean/categorical metrics is encoded in the value itself — no separate p-value, no rank). `allowed_categories` is required when value_kind=categorical.

| Value | value_kind | Unit | Rankable? | Has p-value? | Allowed categories | Source papers |
|---|---|---|---|---|---|---|
| `fourier_score` | numeric | dimensionless | yes (higher = more periodic) | yes (Fisher's g-test) | — | zinser 2009 |
| `peak_time_h` | numeric | hours | **no** (phase, not magnitude) | no | — | zinser 2009 (Waldbauer 2012 on extraction) |
| `peak_fit_r_squared` | numeric | dimensionless (0–1) | yes (higher = better fit) | yes (F-test) | — | zinser 2009 |
| `protein_transcript_lag_h` | numeric | hours | yes (higher = more translationally regulated) | optionally | — | Waldbauer 2012 |
| `damping_ratio` | numeric | dimensionless | yes (higher = more damped) | no | — | Waldbauer 2012 |
| `diel_amplitude` | numeric | dimensionless | yes | no | — | Waldbauer 2012 |
| `periodic_in_axenic_LD` | boolean | — | no | no | — | Biller 2018 (S4A) |
| `periodic_in_coculture_LD` | boolean | — | no | no | — | Biller 2018 (S4A, S4B) |
| `periodic_in_axenic_extended_darkness` | boolean | — | no | no | — | Biller 2018 (S4A) |
| `periodic_in_coculture_extended_darkness` | boolean | — | no | no | — | Biller 2018 (S4A, S4B) |
| `darkness_survival_class` | categorical | — | no | no | `{present_in_axenic, present_in_coculture, unique_to_axenic, unique_to_coculture}` | Biller 2018 (S5) |
| `detected` | boolean | — | no | no | — | Papers reporting only a presence list with no quantification (e.g., future Oleza / Biller detection cases). Single generic metric_type — compartment-specific variants are explicitly NOT used; compartment lives on the parent Experiment. |

New values added only when a concrete paper requires one. Each new `metric_type` declares its `value_kind`, `rankable`, `has_p_value`, and (if categorical) `allowed_categories` in the vocabulary. These flags are metric-level (not per-paper) — a future paper that reports a numeric variant of an existing boolean metric (e.g., a continuous periodicity score) gets a new `metric_type` (`periodicity_score_continuous` or similar) rather than retroactively changing value_kind on the existing one. Same rule applies to `has_p_value` and `rankable`.

### `omics_type` additions (on Experiment)

| Value | Use |
|---|---|
| `PAIRED_RNASEQ_PROTEOME` (new) | Experiments whose only outgoing evidence is paired RNA-seq × proteomics derived metrics (e.g., Waldbauer 2012 transcript-vs-protein dynamics). Invoked per invariant #4. Future pairings of other modalities (e.g., microarray × metabolomics) will introduce their own paired-modality value; naming stays explicit about which two modalities are paired. |

Existing values (`RNASEQ`, `MICROARRAY`, `PROTEOMICS`, `EXOPROTEOMICS`, `METABOLOMICS`) unchanged.

## Paperconfig extension

Reuses existing paperconfig shape: top-level `experiments:` block stays authoritative; new data tables are declared as entries under the existing `supplementary_materials:` block. Three new `type:` values join the existing `csv` / `id_translation` / `annotation_gff` / `gene_clusters` vocabulary. One optional new field (`compartment`) on the `experiments:` block. No new top-level sections.

### Changes to `experiments:` block

**One new optional field:** `compartment` (defaults to `"whole_cell"` when omitted — preserves existing paperconfigs).

```yaml
experiments:
  biller_2022_mit9312_vesicle_proteome:
    name: "MIT9312 HL vesicle proteome"
    organism: "Prochlorococcus MIT9312"
    omics_type: PROTEOMICS
    compartment: vesicle              # NEW — optional, defaults to whole_cell
    treatment_condition: "HL vesicle fraction (no contrast)"
    # control_condition: omitted — no DE contrast
    # test_type: omitted — no statistical test
    treatment_type: []                # empty when no treatment variable (observation-only)
    experimental_context: "MIT9312 HL cells grown in Pro99, vesicles isolated by ultracentrifugation"
    medium: "Pro99"
    temperature: "24C"
    light_condition: "high_light"

  waldbauer_2012_med4_paired_diel:
    name: "MED4 diel paired transcript-protein analysis (Waldbauer 2012)"
    organism: "Prochlorococcus MED4"
    omics_type: PAIRED_RNASEQ_PROTEOME     # NEW value
    compartment: whole_cell
    treatment_condition: "Light-dark synchronized diel cycle (14:10)"
    treatment_type: [diel]
    experimental_context: "Derived per-gene metrics (lag, damping, phase) from paired mRNA + protein time courses in 312 genes"
```

**Validator relaxation** (mirrors existing cluster-only relaxation in `scripts/validate_paperconfig.py`):
- `control_condition` and `test_type` become optional for experiments whose only supplementary entries are `abundance_table` / `timeseries_abundance_table` / `derived_metrics_table`.
- `treatment_type` stays required but may be `[]`.
- `table_scope` stays recommended for DE only.

### New `supplementary_materials:` entry types

**Type `abundance_table`** — per-sample abundance, single sample. Adapter emits one `AbundanceAnalysis` node + N `abundance_analysis_measures_gene` edges.

```yaml
supplementary_materials:
  biller_2022_mit9312_vesicles:
    type: abundance_table
    filename: "data/.../emi15834-sup-0004-tables3.xlsx"
    sheet: "9312HL_CellsVesicles_COMETa"
    organism: "Prochlorococcus MIT9312"
    id_columns:
      - column: "Gene Number"
        id_type: locus_tag
    experiment: biller_2022_mit9312_vesicle_proteome
    name_col: "Gene Number"
    replicate_cols: ["vesicle_rep_A", "vesicle_rep_B", "vesicle_rep_C"]
    value_type: spectral_count_biovolume_norm
    field_description: "Biovolume-normalized spectral counts"  # optional; becomes AbundanceAnalysis.field_description
    # Optional single-timepoint label:
    # timepoint: "Day 2"
    # timepoint_hours: 48
    # Detection-only (paper reports presence without quantification): do NOT use abundance_table.
    # Use derived_metrics_table with value_kind=boolean and metric_type=detected (or a
    # compartment-specific variant). Compartment comes from the parent Experiment.
```

**Type `timeseries_abundance_table`** — per-sample abundance across timepoints. Adapter emits one `AbundanceAnalysis` node per timepoint, plus N `abundance_analysis_measures_gene` edges per analysis, all under the same Experiment.

```yaml
supplementary_materials:
  waldbauer_med4_diel_transcript_timeseries:
    type: timeseries_abundance_table
    filename: "..."
    organism: "Prochlorococcus MED4"
    id_columns:
      - column: "locus_tag"
        id_type: locus_tag
    experiment: waldbauer_med4_rnaseq_diel
    name_col: "locus_tag"
    value_type: tpm
    timepoints:
      - timepoint: "0h"
        timepoint_hours: 0.0
        replicate_cols: ["t0_rep1", "t0_rep2"]
      - timepoint: "2h"
        timepoint_hours: 2.0
        replicate_cols: ["t2_rep1", "t2_rep2"]
      # ...
```

Governed by the "no redundant per-timepoint abundance" invariant. Used only when a paper's per-timepoint signal isn't already captured by DE or derived metrics. Waldbauer 2012 raw time course is NOT ingested in v1 (only derived metrics).

**Type `derived_metrics_table`** — paper-computed per-gene metrics. Adapter emits one `DerivedMetric` node per `metric_type` column + N `derived_metric_quantifies_gene` edges per metric.

```yaml
supplementary_materials:
  zinser_2009_diel_metrics:
    type: derived_metrics_table
    filename: "data/.../Table_S1.csv"
    organism: "Prochlorococcus MED4"
    id_columns:
      - column: "Gene or region"
        id_type: locus_tag_cyanorak
    experiment: zinser_2009_med4_diel_microarray
    name_col: "Gene or region"
    metrics:
      - metric_type: fourier_score             # → DerivedMetric.metric_type
        value_kind: numeric                    # required; matches vocabulary entry → DerivedMetric.value_kind
        value_col: "Fourier"                   # CSV column carrying the value → edge.value
        unit: null                             # → DerivedMetric.unit
        field_description: null                # → DerivedMetric.field_description
        # Optional p-value columns when has_p_value="true" for this metric_type:
        p_value_col: "Fourier P"               # raw p-value column → edge.p_value
        adjusted_p_value_col: "Fourier FDR"    # multiple-testing-corrected column → edge.adjusted_p_value
        p_value_threshold: 0.05                # paper-declared significance cutoff → DerivedMetric.p_value_threshold
      - metric_type: peak_time_h
        value_kind: numeric
        value_col: "Peak"
        unit: "h"
        field_description: "Peak expression time within 48h diel cycle (paper's Fourier-fit)"
        # No p_value_col — has_p_value="false" for peak_time_h
      - metric_type: peak_fit_r_squared
        value_kind: numeric
        value_col: "Peak R value"
        unit: null
        field_description: null
        p_value_col: "Peak P"
        p_value_threshold: 0.05
```

**`derived_metrics_table` — boolean example (Biller 2018 S4B, MIT1002 periodicity in coculture):**

```yaml
supplementary_materials:
  biller_2018_mit1002_periodicity:
    type: derived_metrics_table
    filename: "data/.../table s4B sys003182233st4.csv"
    organism: "Alteromonas macleodii MIT1002"
    id_columns:
      - column: "Locus ID"
        id_type: locus_tag
    experiment: darkness_extended_darkness_mit1002_rnaseq
    name_col: "Locus ID"
    metrics:
      - metric_type: periodic_in_coculture_LD
        value_kind: boolean
        value_col: "Periodic in co-cultured, L:D cultures"   # → edge.value_flag (after token mapping)
        true_tokens: ["Y"]
        false_tokens: []                  # no explicit false cells in this table
        skip_tokens: ["NA", "N/A"]        # default, shown here for clarity
        blank_policy: skip                # CSV lists ONLY positive cases; absence = "not in periodic set"
        field_description: "RAIN periodicity FDR<0.05 in coculture L:D"
      - metric_type: periodic_in_coculture_extended_darkness
        value_kind: boolean
        value_col: "Periodic in co-cultured, extended darkness cultures"
        true_tokens: ["Y"]
        blank_policy: skip
        field_description: "RAIN periodicity FDR<0.05 in coculture extended darkness"
```

(S4A is the same shape with 4 metrics for axenic_LD / coculture_LD / axenic_dark / coculture_dark, organism `Prochlorococcus NATL2A`. The `periodicity_cluster` composite column is dropped — it was a synthetic collapse of these flags.)

**`derived_metrics_table` — categorical example (Biller 2018 S5, NATL2A darkness survival):**

```yaml
supplementary_materials:
  biller_2018_natl2a_darkness_survival:
    type: derived_metrics_table
    filename: "data/.../table S5 sys003182233st5.csv"
    organism: "Prochlorococcus NATL2A"
    id_columns:
      - column: "NCBI ID_2"
        id_type: locus_tag_ncbi
    experiment: darkness_extended_darkness_natl2a_rnaseq_axenic   # link to one experiment; second linkage handled at adapter level if needed
    name_col: "NCBI ID_2"
    metrics:
      - metric_type: darkness_survival_class
        value_kind: categorical
        value_col: "darkness_cluster"      # → edge.value_text (must be in DerivedMetric.allowed_categories)
        # allowed_categories declared in metric_type vocabulary; validator + adapter
        # reject rows whose value isn't in the set.
        field_description: "Transcript presence at 72-144h extended darkness, axenic vs coculture (Biller 2018 Table S5 categories)"
```

### Vocabulary alignments with existing paperconfig

- **`timepoint` / `timepoint_hours`** (singular, no underscore) in paperconfig — reused from existing `statistical_analyses` conventions. At ingest, the adapter stamps `AbundanceAnalysis.time_point` / `time_point_order` / `time_point_hours` (with underscore, matching existing Experiment-level arrays and DE edge properties).
- **`id_columns` with `column` + `id_type`**: reused. `id_type` uses existing vocabulary (`locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `old_locus_tag`, `uniprot_accession`, `protein_id_refseq`, `gene_name`, …).
- **`name_col`**: reused — the column whose values are resolved to Gene nodes via the existing v2 gene-id mapping.
- **`organism`**: reused at supplementary_materials level.
- **`experiment:` reference by key**: reused from `statistical_analyses` pattern.
- **Paperconfig `timepoint` sentinel convention** (`"unknown"` for unknown) applies to `timeseries_abundance_table` entries too.

### Validator + skill updates required (in scope)

- `scripts/validate_paperconfig.py`:
  - Register new types: `abundance_table`, `timeseries_abundance_table`, `derived_metrics_table`.
  - Add `compartment` to valid experiment fields with controlled-vocabulary check against `{whole_cell, vesicle, exoproteome, spent_medium, lysate}`.
  - Add `value_type` controlled-vocabulary check (on `value_type` field of abundance entries).
  - Add `metric_type` controlled-vocabulary check with `value_kind`, `rankable`, `has_p_value`, and `allowed_categories` metadata pulled from the single-source vocabulary table in §2. Validator behavior per value_kind (each value_kind drives a different edge type at adapter time — `derived_metric_quantifies_gene` / `derived_metric_flags_gene` / `derived_metric_classifies_gene`):
    - `value_kind=numeric` → emits `derived_metric_quantifies_gene`. Entry must declare `value_col`, may declare `unit`. When `has_p_value="true"` for the metric_type, validator allows `p_value_col` / `adjusted_p_value_col` / `p_value_threshold`; when `has_p_value="false"`, validator rejects those fields (mismatch is a paperconfig error).
    - `value_kind=boolean` → emits `derived_metric_flags_gene`. Entry must declare `value_col` + `true_tokens` (non-empty list). May declare `false_tokens` (default `[]`), `skip_tokens` (default `["NA", "N/A", "n/a", "#N/A"]`), `blank_policy` (default `"skip"`; one of `{"skip", "false", "true"}`). Rejects `unit`, `p_value_col`, `adjusted_p_value_col`, `p_value_threshold` (meaningless for boolean). Vocabulary entries for boolean metric_types must declare `rankable="false"` and `has_p_value="false"`. Validator dry-runs the source CSV: collects the set of distinct values in `value_col`, and errors if any value isn't in `true_tokens ∪ false_tokens ∪ skip_tokens ∪ {blank}`.
    - `value_kind=categorical` → emits `derived_metric_classifies_gene`. Entry must declare `value_col`. Vocabulary must declare non-null `allowed_categories`. Rejects `unit`, `true_tokens`, `false_tokens`, `skip_tokens`, `blank_policy`, p-value fields. Validator scans the source CSV's `value_col` and warns if any cell's value falls outside `allowed_categories` (adapter then rejects those rows at ingest).
  - When a paperconfig entry's declared `value_kind` disagrees with the vocabulary's `value_kind` for that metric_type, validator rejects — keeps metric_type → value_kind stable.
  - A single `derived_metrics_table` paperconfig entry MAY mix metrics with different `value_kind` values. Each `metrics:` item creates its own DerivedMetric node, and each DerivedMetric emits exactly one edge type matching its own value_kind. Validator does not need to enforce homogeneity per-entry.
  - `abundance_table` entries require `value_col` + `value_type`. `detection_only` mode is NOT supported — detection evidence uses `derived_metrics_table` with `value_kind=boolean` + `metric_type=detected` (validator guides authors to the right shape when a paperconfig tries to use `detection_only`).
  - Extend VALID_TYPES with `PAIRED_RNASEQ_PROTEOME`.
  - Extend DE-only vs observation/metric-only relaxation for `control_condition` / `test_type`.
  - Enforce one-omics_type-per-Experiment and one-timeseries-per-Experiment invariants.
  - **Controlled vocabularies must read from a single source of truth.** The three vocabulary tables in this spec (§3 `compartment`, §3 `value_type`, §3 `metric_type` + `rankable`) are authoritative. Implement them as module-level constants imported by both the validator and the adapter so they can't drift. New vocabulary values added in one place only.
  - **Cross-entry invariant: "no redundant per-timepoint abundance."** No `abundance_table` OR `timeseries_abundance_table` entry may duplicate a `csv` DE entry at the same `(Experiment, Gene, timepoint)` triple. The current validator is single-entry-scoped; enforcing this invariant requires cross-entry awareness (scan all entries in a paperconfig sharing the same `experiment:` key). Implementation options: (a) extend validator with a second pass that groups entries by experiment, or (b) lift the check to KG-validity post-build. Spec prefers (a) so errors are caught pre-ingest; implementation may fall back to (b) if (a) proves impractical.

#### Adapter-side consistency checks (separate from validator — require row-level inspection)

These checks live in the adapter, not the validator, because they need to inspect actual CSV row values during ingest:

- Every `abundance_analysis_measures_gene` edge must have a non-null numeric `value` (and, post-import, non-null `rank_in_sample` / `abundance_percentile` / `abundance_bucket`). Adapter rejects rows with missing or non-numeric values at ingest — paperconfig author either drops the row, provides a value, or switches to `derived_metrics_table` + `metric_type=detected` if the paper only reports presence.
- Every emitted `derived_metric_classifies_gene` edge must have `value_text` ∈ parent `DerivedMetric.allowed_categories`. Out-of-set values cause a hard error at ingest (validator may have warned earlier; adapter is the last line of defense).
- Boolean token parsing fails fast on unexpected tokens (per §2 boolean semantics) — adapter does not silently coerce.
- Denormalized fields (`compartment`, `omics_type`, `treatment_type`, `background_factors`, `treatment`, `light_condition`, `experimental_context`, `publication_doi`, `organism_name`) on every emitted analysis node must equal the parent Experiment's values. Adapter computes them from the parent rather than re-reading the paperconfig — drift impossible.
- `.claude/skills/paperconfig/SKILL.md`:
  - New "Type `abundance_table`" section under Supplementary Materials.
  - New "Type `timeseries_abundance_table`" section.
  - New "Type `derived_metrics_table`" section.
  - `compartment` field in Experiments Block section (optional; default `whole_cell`).
  - `PAIRED_RNASEQ_PROTEOME` omics_type documented with guidance on when to use.
  - `value_type` and `metric_type` vocabulary tables alongside existing Canonical Vocabulary.
  - Edge routing note: `abundance_table` → `AbundanceAnalysis` + `abundance_analysis_measures_gene`; `derived_metrics_table` → `DerivedMetric` + `derived_metric_quantifies_gene`.
  - Examples: biller 2022 vesicles (abundance_table), zinser 2009 Fourier metrics (derived_metrics_table), Waldbauer paired (derived_metrics_table on PAIRED_RNASEQ_PROTEOME experiment).

## Post-import computations

Add to `scripts/post-import.sh` (and `scripts/post-import.cypher` reference copy) the following new blocks:

1. **Rank `abundance_analysis_measures_gene` edges** grouped by `(AbundanceAnalysis)`. Every edge has a non-null `value` (detection-only evidence lives on `derived_metric_flags_gene` instead). Compute `rank_in_sample` (1 = most abundant), `abundance_percentile` (0–100), and `abundance_bucket` (`top_decile` if percentile ≥90; `top_quartile` if ≥75 and <90; `mid` if ≥25 and <75; `low` otherwise). Thresholds pinned in this spec — do not drift.

2. **Rank + percentile + significance for numeric derived metrics** (boolean/categorical edges have no rank or p-value passes — significance is encoded directly in `value_flag` / `value_text`).
   - Rank `derived_metric_quantifies_gene` edges grouped by `(DerivedMetric)` — only when `DerivedMetric.rankable="true"`. Compute `rank_by_metric` (1 = highest value), `metric_percentile` (0–100), and `metric_bucket` using the same pinned thresholds as `abundance_bucket`: `top_decile` (≥90), `top_quartile` (≥75 and <90), `mid` (≥25 and <75), `low` (<25).
   - Derive `significant` for `derived_metric_quantifies_gene` edges whose parent has `has_p_value="true"` AND `p_value_threshold IS NOT NULL`: `significant = "true"` iff `adjusted_p_value < p_value_threshold`, else `"false"`. Null when adjusted_p_value is null on the row (string enum, not bool; matches existing DE convention).

3. **Analysis-node rollups.** Compute identical computed-field sets on both `AbundanceAnalysis` and `DerivedMetric` (mirrors `ClusteringAnalysis`):
   - `total_gene_count`: count of outgoing measurement edges. For `AbundanceAnalysis`: `abundance_analysis_measures_gene` edges. For `DerivedMetric`: the single applicable edge type — `derived_metric_quantifies_gene` (if value_kind=numeric), `derived_metric_flags_gene` (if value_kind=boolean), or `derived_metric_classifies_gene` (if value_kind=categorical). Each DerivedMetric only emits one edge type, so the count is unambiguous.
   - `growth_phases: str[]`: union of the parent Experiment's `growth_phases` (parallels the existing `ClusteringAnalysis.growth_phases` post-import in `scripts/post-import.cypher`). Empty array when the parent has none.

4. **Experiment rollups** (parallels existing `Experiment.clustering_analysis_count` / `cluster_types` / `cluster_count` computed from child ClusteringAnalysis nodes).
   - `Experiment.reports_fold_change`: `"true"` iff outgoing `Changes_expression_of` edges exist, else `"false"` (string enum).
   - `Experiment.reports_derived_metric_types`: distinct `metric_type` values reachable via `experiment_has_derived_metric` → DerivedMetric. Parallels `cluster_types`.
   - `Experiment.abundance_analysis_count`: count of child `AbundanceAnalysis` nodes (via `experiment_has_abundance_analysis`). Parallels `clustering_analysis_count`.
   - `Experiment.derived_metric_count`: count of child `DerivedMetric` nodes (via `experiment_has_derived_metric`).
   - `Experiment.derived_metric_value_kinds: str[]`: sorted distinct `value_kind` values across child DerivedMetrics (`{"numeric","boolean","categorical"}` subset). Parallels the cluster_types rollup on Experiment.
   - `Experiment.abundance_gene_count`: distinct genes reachable via any child `AbundanceAnalysis`.
   - `Experiment.derived_metric_gene_count`: distinct genes reachable via any child `DerivedMetric` (across all three derived-metric edge types — quantifies/flags/classifies).
   - Abundance time-point arrays (parallel to existing DE time-point arrays, see §3):
     - `Experiment.abundance_time_point_labels: str[]` (ordered by `time_point_order`; `""` when unlabeled)
     - `Experiment.abundance_time_point_orders: int[]`
     - `Experiment.abundance_time_point_hours: float[]` (uses `-1.0` sentinel)
     - `Experiment.abundance_time_point_totals: int[]` — per-tp gene count
   - `Experiment.is_time_course` may need re-evaluation to include abundance multi-tp children (true iff either DE or abundance children span >1 timepoint). Existing adapter-level derivation from paperconfig is authoritative where present; post-import promotes to true if abundance children extend it.

5. **Publication rollups** (parallels existing `Publication.clustering_analysis_count` / `cluster_types` / `cluster_count`, using the new `publication_has_abundance_analysis` + `publication_has_derived_metric` binding edges).
   - `Publication.abundance_analysis_count: int`: count of child AbundanceAnalysis nodes.
   - `Publication.derived_metric_count: int`: count of child DerivedMetric nodes.
   - `Publication.abundance_gene_count: int`: distinct genes reachable via any child AbundanceAnalysis.
   - `Publication.derived_metric_gene_count: int`: distinct genes reachable via any child DerivedMetric (across all 3 edge types).
   - `Publication.compartments: str[]`: sorted distinct `compartment` values across child Experiments. Parallels existing `Publication.omics_types`.
   - `Publication.derived_metric_types: str[]`: sorted distinct `metric_type` values across child DerivedMetrics. Parallels `cluster_types`.
   - `Publication.derived_metric_value_kinds: str[]`: sorted distinct `value_kind` values across child DerivedMetrics.

6. **OrganismTaxon rollups** (parallels existing `OrganismTaxon.clustering_analysis_count` / `cluster_types` / `cluster_count`, using the new `abundance_analysis_belongs_to_organism` + `derived_metric_belongs_to_organism` binding edges).
   - `OrganismTaxon.abundance_analysis_count: int`
   - `OrganismTaxon.derived_metric_count: int`
   - `OrganismTaxon.abundance_gene_count: int`
   - `OrganismTaxon.derived_metric_gene_count: int`
   - `OrganismTaxon.compartments: str[]`
   - `OrganismTaxon.derived_metric_types: str[]`
   - `OrganismTaxon.derived_metric_value_kinds: str[]`

7. **Gene routing-count rollups** (parallels existing `expression_edge_count`, `cluster_membership_count`, etc.). Split per derived-metric edge type so the routing tool can dispatch to the right MCP tool:
   - `Gene.abundance_analysis_count`: distinct AbundanceAnalysis nodes that quantify this gene.
   - `Gene.numeric_metric_count`, `Gene.classifier_flag_count`, `Gene.classifier_label_count`: distinct DerivedMetric nodes reaching this gene via `derived_metric_quantifies_gene` / `derived_metric_flags_gene` / `derived_metric_classifies_gene` respectively.
   - `Gene.compartments_observed`: sorted distinct `compartment` values across all incoming non-DE edges (via parent analysis nodes).
   - `Gene.numeric_metric_types_observed`, `Gene.classifier_flag_types_observed`, `Gene.classifier_label_types_observed`: sorted distinct `metric_type` values per edge type.

8. **Indexes** (scalar):
   - `abundance_analysis_value_type_idx` on `AbundanceAnalysis(value_type)`
   - `abundance_analysis_compartment_idx` on `AbundanceAnalysis(compartment)`
   - `abundance_analysis_omics_type_idx` on `AbundanceAnalysis(omics_type)`
   - `abundance_analysis_treatment_type_idx` on `AbundanceAnalysis(treatment_type)`
   - `abundance_analysis_organism_idx` on `AbundanceAnalysis(organism_name)`
   - `abundance_analysis_experiment_idx` on `AbundanceAnalysis(experiment_id)`
   - `derived_metric_metric_type_idx` on `DerivedMetric(metric_type)`
   - `derived_metric_value_kind_idx` on `DerivedMetric(value_kind)`
   - `derived_metric_compartment_idx` on `DerivedMetric(compartment)`
   - `derived_metric_omics_type_idx` on `DerivedMetric(omics_type)`
   - `derived_metric_treatment_type_idx` on `DerivedMetric(treatment_type)`
   - `derived_metric_organism_idx` on `DerivedMetric(organism_name)`
   - `derived_metric_experiment_idx` on `DerivedMetric(experiment_id)`
   - `experiment_compartment_idx` on `Experiment(compartment)`

9. **Full-text indexes** (parallel `clusteringAnalysisFullText`):
   - `abundanceAnalysisFullText` on `AbundanceAnalysis(name, field_description)`
   - `derivedMetricFullText` on `DerivedMetric(name, field_description)`

Follow existing post-import style: group multiple statements into one `cypher-shell` invocation per concern (ranks / rollups / indexes); emit `[timing]` lines.

### Rank property conventions (cross-edge mapping)

| Edge | Property | Grouped by | Ranks by | 1 means |
|---|---|---|---|---|
| `Changes_expression_of` (existing, unchanged) | `rank_by_effect` | (Experiment, time_point_order) | `|log2_fold_change|` | strongest |FC| |
| `Changes_expression_of` (existing, unchanged) | `rank_up` | (Experiment, time_point_order), among `expression_direction='up'` | `log2_fold_change` | most upregulated |
| `Changes_expression_of` (existing, unchanged) | `rank_down` | (Experiment, time_point_order), among `expression_direction='down'` | `|log2_fold_change|` | most downregulated |
| `abundance_analysis_measures_gene` (new) | `rank_in_sample` | (AbundanceAnalysis), every edge | `value` | most abundant |
| `derived_metric_quantifies_gene` (new, numeric) | `rank_by_metric` | (DerivedMetric), only if `rankable="true"` | `value` | highest value |
| `derived_metric_quantifies_gene` (new, numeric) | `metric_percentile` (float 0-100) | (DerivedMetric), only if `rankable="true"` | `value` | higher percentile = higher value |
| `derived_metric_quantifies_gene` (new, numeric) | `metric_bucket` (enum) | (DerivedMetric), only if `rankable="true"` | derived from `metric_percentile` | top_decile/top_quartile/mid/low |
| `derived_metric_quantifies_gene` (new, numeric) | `significant` (str `"true"` \| `"false"`) | per-edge, only if `has_p_value="true"` | `adjusted_p_value < DerivedMetric.p_value_threshold` | `"true"` = paper-significant |
| `derived_metric_flags_gene` (new, boolean) | none | n/a | n/a — `value_flag` (str `"true"` \| `"false"`) is the answer | n/a |
| `derived_metric_classifies_gene` (new, categorical) | none | n/a | n/a — `value_text` is the answer | n/a |

## Ingest invariants (schema-enforced policies)

Invariants that KG validity tests check:

1. **Compartment uniqueness per Experiment.** Exactly one `compartment` value per Experiment. Different compartments from one paper → separate Experiments.
2. **One omics_type per Experiment.** Paperconfig convention (Weissberg 2025 already splits RNA-seq vs PROTEOMICS as separate Experiments). Validator enforces.
3. **One timeseries per Experiment.** If `Experiment.is_time_course=true`, its `experiment_has_abundance_analysis` / `experiment_has_derived_metric` / `Changes_expression_of` children describe a single ordered timepoint sequence. Multi-dimensional designs (e.g., 2 conditions × 4 timepoints) split into multiple Experiments.
4. **Multi-modality paired metrics → comparison Experiment.** Papers whose evidence is paired-modality derived metrics (Waldbauer 2012 transcript × protein lag/phase/damping) attach those metrics to an Experiment with `omics_type=PAIRED_RNASEQ_PROTEOME`, not to one of the source-modality Experiments. For Waldbauer v1 where source-modality raw time courses aren't ingested, this is the *only* Experiment for the paper.
5. **No redundant per-timepoint abundance.** An `AbundanceAnalysis` + `abundance_analysis_measures_gene` edge is NOT produced for an `(Experiment, Gene, time_point)` triple if that same triple is already covered by a `Changes_expression_of` edge. Prevents double-encoding of the same measurement. (Derived metrics are scalar summaries without a time_point property, so they cannot collide with per-timepoint abundance and are not checked by this invariant.)
6. **Gene-only targets.** `abundance_analysis_measures_gene` and `derived_metric_quantifies_gene` edges target `:Gene` only. No Protein, no Metabolite in this spec.
7. **Replicate SD requires n ≥ 2.** `value_sd` is null when `n_replicates < 2` or when the paper doesn't report spread. Adapter converts SEM → SD via SEM·√n when `n_replicates` is known.
8. **Time-point hours uses null, not sentinel values.** Scalar `time_point_hours` (on `AbundanceAnalysis` node and `Changes_expression_of` edge) is null when unknown — no `-1.0` sentinel. Array-valued `Experiment.time_point_hours: float[]` still uses `-1.0` per existing convention (Neo4j arrays cannot contain nulls).
9. **Denormalized fields match their parent.** All Experiment-denormalized fields on `AbundanceAnalysis` and `DerivedMetric` (`experiment_id`, `organism_name`, `publication_doi`, `compartment`, `omics_type`, `treatment_type`, `background_factors`, `treatment`, `light_condition`, `experimental_context`) match the parent Experiment's values exactly. Post-import `growth_phases` on both analysis nodes equals the parent Experiment's `growth_phases`. Edge-level denormalizations (`abundance_analysis_measures_gene.value_type`, `derived_metric_quantifies_gene.metric_type`, `derived_metric_flags_gene.metric_type`, `derived_metric_classifies_gene.metric_type`) match their parent analysis node. KG validity tests assert the equality.
10. **One edge type per DerivedMetric, matching value_kind.** Each `DerivedMetric` emits exactly one of the three derived-metric edge types, matching its own `value_kind`:
    - `value_kind=numeric` → only `derived_metric_quantifies_gene` outgoing edges.
    - `value_kind=boolean` → only `derived_metric_flags_gene` outgoing edges.
    - `value_kind=categorical` → only `derived_metric_classifies_gene` outgoing edges.
    No DerivedMetric emits two different edge types. KG validity tests assert this.
11. **Categorical values must be in allowed set.** For every `derived_metric_classifies_gene` edge, `value_text` is a member of the parent `DerivedMetric.allowed_categories`.

## What's in scope for this spec

The scope is **code changes to ingest the new schema** plus a **"what changed" doc for downstream consumers** (explorer / MCP teams own their tool-layer changes separately — see "Out of scope" below).

### Schema / graph layer

- **`config/schema_config.yaml`** — register two new node types (`AbundanceAnalysis`, `DerivedMetric`), ten new edge types (6 binding + 4 measurement), all properties listed below. Update existing `Experiment`, `Publication`, `OrganismTaxon`, `Gene` node definitions to declare new post-import-computed properties.
- **Node/edge additions** (formal summary):
  - Nodes: `AbundanceAnalysis`, `DerivedMetric`
  - Binding edges (mirror ClusteringAnalysis pattern): `publication_has_abundance_analysis`, `experiment_has_abundance_analysis`, `abundance_analysis_belongs_to_organism`, `publication_has_derived_metric`, `experiment_has_derived_metric`, `derived_metric_belongs_to_organism`
  - Measurement edges: `abundance_analysis_measures_gene`, `derived_metric_quantifies_gene` (numeric), `derived_metric_flags_gene` (boolean), `derived_metric_classifies_gene` (categorical)
- **New Experiment properties**: `compartment` (adapter-emitted) + post-import `reports_fold_change`, `reports_derived_metric_types`, `abundance_analysis_count`, `derived_metric_count`, `derived_metric_value_kinds`, `abundance_gene_count`, `derived_metric_gene_count`, and 4 abundance time-point arrays. DE time-point arrays stay DE-scoped.
- **New Publication + OrganismTaxon properties** (post-import): 7 fields each — see §3 and post-import blocks #5/#6.
- **New Gene routing properties** (post-import): 8 fields — see §5.
- **New `omics_type` value**: `PAIRED_RNASEQ_PROTEOME`.

### Adapter layer

- **`multiomics_kg/adapters/omics_adapter.py`** (or a sibling module — decide during implementation): emit the two new analysis nodes + 10 new edge types; resolve protein IDs → genes for proteomics abundance; aggregate replicates; stamp denormalized fields from parent Experiment.
- **Token-parsing logic** for boolean metrics: implement `true_tokens` / `false_tokens` / `skip_tokens` / `blank_policy` with hard-error on unexpected tokens.
- **Invariants enforced at emission**: exactly one edge type per DerivedMetric matching its `value_kind`; `allowed_categories` membership for categorical; every `abundance_analysis_measures_gene` edge has non-null numeric `value`; denormalized fields equal parent Experiment.

### Post-import layer

- **`scripts/post-import.sh`** and **`scripts/post-import.cypher`** (must stay byte-identical after diff test per CLAUDE.md): add 9 blocks — rank passes for abundance + numeric derived metrics, analysis-node rollups, Experiment rollups, Publication rollups, OrganismTaxon rollups, Gene routing rollups, scalar indexes, full-text indexes.

### Paperconfig touchpoints (everywhere paperconfig is read)

- **`scripts/validate_paperconfig.py`** — register new types (`abundance_table`, `timeseries_abundance_table`, `derived_metrics_table`), new `compartment` field on experiments block, controlled-vocabulary checks (`compartment`, `value_type`, `metric_type`, `value_kind`, `allowed_categories`), per-value_kind field gating, cross-entry "no redundant per-timepoint abundance" check.
- **`multiomics_kg/download/build_gene_id_mapping.py`** — harvest `id_columns` from the three new entry types (abundance + timeseries-abundance + derived-metrics) so paper-specific IDs get added to `gene_id_mapping.json` v2. Parallels existing handling of `csv` `id_columns`.
- **`multiomics_kg/download/resolve_paper_ids.py`** — extend pre-resolution to the new entry types: if `name_col` on an abundance/derived-metric table isn't already `locus_tag`, write a `<stem>_resolved.csv` with `locus_tag` + `resolution_method` columns. Adapter then probes for the resolved file (mirrors existing `csv`-type handling).
- Any other paperconfig-scanning code surfaces discovered during implementation — update to recognize the new entry types (at minimum: silently skip; at best: surface in reports).

### Controlled vocabularies (single source of truth)

Create a shared Python module at **`multiomics_kg/vocab/non_de_evidence.py`** exposing:
- `COMPARTMENTS: set[str]`
- `VALUE_TYPES: set[str]` (AbundanceAnalysis)
- `METRIC_TYPES: dict[str, MetricSpec]` where each entry declares `value_kind`, `unit`, `rankable`, `has_p_value`, `allowed_categories` (as a `dataclass` or `TypedDict`)
- `OMICS_TYPES: set[str]` (extend existing omics-type set with `PAIRED_RNASEQ_PROTEOME`)

Both the validator (`scripts/validate_paperconfig.py`) and the adapter (`multiomics_kg/adapters/omics_adapter.py` or its sibling) import from this module — no drift possible. If an existing config-loader pattern elsewhere in the repo is more idiomatic, the implementer may move the definitions there provided the "one source, two consumers" property holds.

### Test layer

- **Adapter unit tests** (`tests/test_omics_adapter_*.py` or new file): per-value_kind emission; token parsing (true/false/skip/blank/unexpected → error); denormalization correctness; `experiment_id` plumbing. One of the fixtures MUST exercise `metric_type=detected` in positive-only mode (`blank_policy=skip`, `false_tokens=[]`) end-to-end — this is the detection-only replacement pattern and is a new path the adapter must cover.
- **KG validity tests** (`tests/kg_validity/`): the success-criteria assertions in this spec — node counts, edge type homogeneity per DerivedMetric, `allowed_categories` membership, denormalized-fields-match-parent, rank / percentile / bucket / significant consistency, Publication + OrganismTaxon + Experiment rollup consistency, binding-edge cardinality, no redundant per-timepoint abundance.
- **Post-import regression** (`scripts/post-import-validate.sh`): existing deterministic-dump comparison must remain byte-identical for the DE-only subset; new dump includes the new properties and rollups.
- **`/omics-edge-snapshot` skill**: extend to include counts + per-paper breakdowns for the three new measurement edge types alongside existing `Changes_expression_of`.
- **Snapshot fixture**: regenerate `tests/kg_validity/snapshot_data.json` to cover at least one node/edge of each new type (or guard with a follow-up task note if no paperconfigs exist yet to exercise them at spec-close).

### Downstream communication — "what changed" doc

- **New artifact**: `docs/kg-changes/non-de-evidence-extension.md` (path parallels existing `docs/kg-changes/*.md` series — see `reference-proteome-match-organisms.md`, `ontology-level.md`, `brite-categories.md`). Single consolidated reference for the explorer + MCP teams: new nodes, new edges, new properties, new vocabularies, Cypher query patterns (subset of this spec's "Cross-paper query patterns"), the flat-vs-analysis-node asymmetry between DE and non-DE evidence, and how to dispatch MCP envelopes by edge type / `value_kind`. This is a **doc deliverable**, not implementation guidance — MCP team consumes it to plan their tool changes separately.

### Skill docs (internal)

- **`.claude/skills/paperconfig/SKILL.md`** — document the three new entry types + `compartment` field + `PAIRED_RNASEQ_PROTEOME` omics_type + controlled-vocabulary tables. Needed so a human (or a future `/paperconfig` skill invocation) can *author* new paperconfigs against the new schema. (Actually filling paperconfigs out against existing data is a follow-up task — out of scope here.)
- **`.claude/skills/cypher-queries/SKILL.md`** — add query templates for the new edge types (mirrors the cross-paper query section of this spec).
- **`.claude/skills/omics-edge-snapshot/SKILL.md`** — update to reflect the new edge types the snapshot should cover.

## What's explicitly NOT in this spec

- **MCP tool changes.** The explorer / MCP layer consumes the new schema but owns its own tool-surface changes. This spec describes the schema exhaustively (§MCP surface deltas kept for reference, but those are *downstream* tool-design sketches, not work for this spec). Implementation here stops at the graph + validator + adapter + post-import + "what changed" doc. The doc (`docs/kg-changes/non-de-evidence-extension.md`) is how the MCP team learns what's new; they file their own follow-up to update tools.
- **Creating / converting paperconfigs for real papers.** Biller 2018 S4A/S4B/S5 retrofit, zinser 2009 metric retrofit, Biller 2022 vesicle ingestion, Oleza 2015 / Oleza 2017 / Waldbauer 2012 ingestion — all these require actual paperconfig authoring + data-file prep + validation passes + KG rebuild. Those are **follow-up tasks**, not this spec's implementation. This spec ships the *infrastructure*; the backlog table (§Backlog impact) is a map for the follow-ups. Default behavior at spec close: the code + schema + validator support the new entry types; no paperconfig yet uses them; KG build still succeeds; tests pass (with the new-entry-type tests exercising synthetic fixtures, not real papers).
- **DE retrofit to analysis-node pattern.** `Changes_expression_of` stays flat. A follow-up spec may restructure existing DE edges to mirror the new `Experiment → Analysis → Gene` shape; deferred because existing DE works, MCP is built around it, and the migration risk is non-trivial. Asymmetry documented in §4 and invariants.
- **Metabolite / lipid / non-gene entity nodes** — see Appendix.
- **Protein nodes as edge targets** — Gene-only (orphan-proteins issue).
- **`Finding` / narrative-claim nodes** (for Zhang 2021-style narrative-only papers) — separate spec.
- **Retrofit of per-timepoint raw abundance** for existing papers (zinser 2009 25-timepoint RMA, alonso 2023 per-condition transcripts, Capovilla 2023 per-sample counts). zinser derived-metric retrofit *is* in scope; per-timepoint abundance retrofit is deferred.
- **Cross-paper *interpretation* on Gene nodes** (`stress_response_consistency`, signature gene set flags, "observation modalities" categorizations) — analysis-layer. Routing *counts* (`abundance_analysis_count`, `compartments_observed`, etc.) are in scope per §5 — they follow the existing `expression_edge_count` / `cluster_membership_count` pattern.
- **Pathway enrichment, signature gene sets, stress classifications** — MCP / Python.
- **Numeric-only restriction on derived metrics** — superseded; v1 covers numeric, boolean, and categorical (see §2 and `value_kind`).

## Backlog impact (follow-up tasks, NOT in this spec's implementation)

The tables below are a **map for follow-up paperconfig-authoring tasks**, not work for this spec. This spec ships the infrastructure; each row becomes its own small task that (a) creates or edits a paperconfig, (b) preps source data files if needed, (c) runs the validator, (d) rebuilds the KG, (e) runs `/omics-edge-snapshot` + KG validity tests, (f) updates relevant skill docs if new patterns surface.

| Paper | New Experiments | New analysis nodes | New edges |
|---|---|---|---|
| biller 2022 (proteomics S2/S3) | 2 (cells, vesicles) | 2 AbundanceAnalysis | ~3,500 abundance_analysis_measures_gene |
| Oleza 2015 | 1 (exoproteome) | 1 AbundanceAnalysis | ~500 abundance_analysis_measures_gene |
| Waldbauer 2012 | 1 (PAIRED_RNASEQ_PROTEOME) | 3 DerivedMetric numeric (lag, phase, damping) | ~900 derived_metric_quantifies_gene |
| Oleza 2017 (optional R. pomeroyi abundance) | 1 | 1 AbundanceAnalysis | ~1,000 abundance_analysis_measures_gene |
| zinser 2009 (retrofit on existing Experiment) | 0 new | 3 DerivedMetric numeric (Fourier, peak_time, peak_R²) | ~11,000 derived_metric_quantifies_gene |
| Biller 2018 S4A (NATL2A periodicity, retrofit) | 0 new | 4 DerivedMetric boolean (axenic/coculture × LD/dark) | ~4 × 1,800 = ~7,200 derived_metric_flags_gene |
| Biller 2018 S4B (MIT1002 periodicity, retrofit) | 0 new | 2 DerivedMetric boolean (coculture × LD/dark) | ~2 × 530 = ~1,060 derived_metric_flags_gene |
| Biller 2018 S5 (NATL2A darkness survival, retrofit) | 0 new | 1 DerivedMetric categorical (`darkness_survival_class`) | ~270 derived_metric_classifies_gene |

**Plumbing-only (no schema change — straight paperconfig work):**

| Paper | Edges |
|---|---|
| ziegler 2025 | already has paperconfig; add path to list |
| McDonagh 2012 | plain DE (XLS→CSV) |
| Pandhal 2007 | plain DE (XLS→CSV) |
| Thompson 2016 | 34 DE CSVs as multiple statistical_analyses |
| Szul 2019 | inspect DOCX; likely plain DE |

**Totals:** ~5 new Experiments, ~17 analysis nodes, ~25,400 new measurement edges across the three derived-metric edge types + `abundance_analysis_measures_gene`. Plus ~51 new binding edges (3 per analysis node: Publication + Experiment + Organism, mirroring cluster pattern). Existing ~227K `Changes_expression_of` unchanged. Biller 2018 retrofit additionally **drops** the 3 force-fit `gene_clusters` entries (`natl2a_periodicity`, `mit1002_periodicity`, `natl2a_darkness_survival`) and their composite `periodicity_cluster` / `darkness_cluster` synthetic columns from the existing paperconfig — net: cleaner schema for those tables.

**Deferred (own specs):**

- biller 2022 metabolite/lipid tables (S1, S5, S6, S7) — metabolite layer spec.
- Kujawinski 2023 — metabolite layer.
- Capovilla 2023 metabolite table — metabolite layer.
- Ma 2022 metabolomics (if any) — metabolite layer.
- Zhang 2021 — narrative-only, `Finding` nodes.
- Labban 2022 — blocked on custom gene ID annotation from authors.

## Cross-paper query patterns this enables

All queries traverse `(Experiment) → (Analysis) → (Gene)` for new evidence; existing DE queries traverse `(Experiment) → (Gene)`.

> **Cypher label convention.** The examples below use schema-level snake_case for readability (e.g., `experiment_has_abundance_analysis`). Actual Neo4j relationship labels are PascalCase-with-underscores per BioCypher conversion (e.g., `Experiment_has_abundance_analysis`, parallel to existing `Experiment_has_clustering_analysis`). When copy-pasting into cypher-shell, capitalize the first letter of each label. CLAUDE.md "Actual Neo4j labels" is the authoritative list.

1. **Compartment-scoped evidence.** "Which genes are detected in the vesicle fraction across any study":
   ```
   MATCH (e:Experiment {compartment: 'vesicle'})-[:experiment_has_abundance_analysis]->(a:AbundanceAnalysis)-[:abundance_analysis_measures_gene]->(g:Gene)
   RETURN g, a, e
   ```

2. **Cross-paper diel phase.** "Genes with peak expression between hour 4 and 8":
   ```
   MATCH (d:DerivedMetric {metric_type: 'peak_time_h'})-[r:derived_metric_quantifies_gene]->(g:Gene)
   WHERE r.value >= 4 AND r.value <= 8
   RETURN g, d, r
   ```
   Aggregated across zinser 2009 and any future paper that reports `peak_time_h`.

3. **Translational regulation signals.** "Genes with transcript-protein lag > 3h":
   ```
   MATCH (d:DerivedMetric {metric_type: 'protein_transcript_lag_h'})-[r:derived_metric_quantifies_gene]->(g:Gene)
   WHERE r.value > 3
   ```

3b. **Cross-paper top-decile by metric** (percentile + bucket normalize across studies with different N). "Genes in the top decile for Fourier periodicity in any study":
   ```
   MATCH (d:DerivedMetric {metric_type: 'fourier_score'})-[r:derived_metric_quantifies_gene {metric_bucket: 'top_decile'}]->(g:Gene)
   RETURN g, d, r.value, r.metric_percentile
   ```
   Same pattern as query #4 on abundance — denormalized `metric_bucket` means index-friendly lookup.

4. **Cross-paper abundance dominance.** "Genes consistently top-decile in the exoproteome" — using denormalized `compartment` and materialized `abundance_bucket` for index-friendly lookups:
   ```
   MATCH (a:AbundanceAnalysis {compartment: 'exoproteome'})-[r:abundance_analysis_measures_gene {abundance_bucket: 'top_decile'}]->(g:Gene)
   RETURN g, a, r
   ```

5. **Cross-modality paired evidence** (per gene). "Genes with DE edges under N-stress AND exoproteome detection":
   ```
   MATCH (g:Gene)
   WHERE EXISTS { (e1:Experiment)-[:Changes_expression_of]->(g) WHERE 'nitrogen' IN e1.treatment_type }
     AND EXISTS { (a:AbundanceAnalysis {compartment: 'exoproteome'})-[:abundance_analysis_measures_gene]->(g) }
   ```
   Combines the two shapes (flat DE, two-hop non-DE). The denormalized `compartment` on AbundanceAnalysis collapses the lookup to one hop. Queries that blend the shapes commonly will benefit from any future DE retrofit.

6. **Cross-paper boolean classifications** (Biller 2018 retrofit query). "Genes flagged as periodic in coculture L:D conditions across any paper":
   ```
   MATCH (d:DerivedMetric {metric_type: 'periodic_in_coculture_LD'})-[r:derived_metric_flags_gene]->(g:Gene)
   WHERE r.value_flag = 'true'
   RETURN g, d
   ```
   Single edge type, single property check — no cluster-string substring match.

7. **Cross-paper combined boolean filter.** "Genes periodic in coculture L:D AND coculture darkness":
   ```
   MATCH (g:Gene)
   WHERE EXISTS {
           (:DerivedMetric {metric_type: 'periodic_in_coculture_LD'})-[r1:derived_metric_flags_gene]->(g)
           WHERE r1.value_flag = 'true'
         }
     AND EXISTS {
           (:DerivedMetric {metric_type: 'periodic_in_coculture_extended_darkness'})-[r2:derived_metric_flags_gene]->(g)
           WHERE r2.value_flag = 'true'
         }
   ```
   Replaces the awkward `WHERE gc.id CONTAINS 'coculture_LD' AND gc.id CONTAINS 'coculture_darkness'` pattern that the current `gene_clusters` force-fit requires.

8. **Cross-paper categorical class lookup** (Biller 2018 S5). "NATL2A genes that survive in extended darkness uniquely under coculture":
   ```
   MATCH (d:DerivedMetric {metric_type: 'darkness_survival_class'})-[r:derived_metric_classifies_gene]->(g:Gene)
   WHERE r.value_text = 'unique_to_coculture'
   RETURN g, d
   ```

Routing counts (`Gene.abundance_analysis_count`, `Gene.numeric_metric_count`, `Gene.classifier_flag_count`, `Gene.classifier_label_count`, `Gene.compartments_observed`, plus the three `*_types_observed` arrays) are materialized on Gene per §5; all other gene-level interpretation is query-time.

## MCP surface deltas (downstream reference only — NOT in scope for this spec)

**This section is not implementation work for this spec.** It sketches how the explorer / MCP team would likely consume the new schema once this spec's code lands. Treat it as a starting point for their own follow-up spec. The source of truth for the schema itself is §1–§5 above; this section exists so the "what changed" doc (`docs/kg-changes/non-de-evidence-extension.md`) has a reasonable sketch to cite.

Concrete integration checklist.

**New tools** (split per edge type — clean envelopes, no `value_kind` dispatch in the consumer):

- `abundance_by_gene(locus_tag, organism=None, compartment=None, value_type=None)` — uses denormalized `value_type` on the edge and `compartment` on the analysis for 1-hop filtering. Returns per-AbundanceAnalysis rank/percentile/abundance_bucket + value + value_type + parent Experiment + paper. (Detection-only evidence is served separately by `classifier_flags_by_gene` on `metric_type=detected`.)
- `numeric_metrics_by_gene(locus_tag, metric_type=None)` — only `derived_metric_quantifies_gene` edges; returns `value` + `rank_by_metric` / `metric_percentile` / `metric_bucket` (may be null for non-rankable) + `p_value` / `adjusted_p_value` / `significant` (when applicable) + `unit` + DerivedMetric.field_description + parent Experiment + paper.
- `classifier_flags_by_gene(locus_tag, metric_type=None)` — only `derived_metric_flags_gene` edges; returns `value_flag` (string `"true"`/`"false"`) + DerivedMetric.field_description + parent Experiment + paper. No rank/p-value fields in envelope.
- `classifier_labels_by_gene(locus_tag, metric_type=None)` — only `derived_metric_classifies_gene` edges; returns `value_text` (with `allowed_categories` from DerivedMetric for context) + field_description + parent Experiment + paper.
- `list_abundance_analyses(compartment=None, value_type=None, organism=None, treatment_type=None, search=None)` — parallels `list_clustering_analyses`; uses denormalized fields + `abundanceAnalysisFullText` for `search`.
- `list_derived_metrics(metric_type=None, value_kind=None, compartment=None, organism=None, treatment_type=None, search=None)` — single list tool spans all three derived-metric kinds (DerivedMetric is one node type). Filter by `value_kind` to scope to numeric/boolean/categorical.

**Extended tools:**
- `list_experiments` gains filters: `compartment`, `reports_fold_change`, `reports_derived_metric_types`, `is_time_course`, `omics_type`. Summary mode reports `gene_count` (DE), `abundance_gene_count`, `derived_metric_gene_count` so non-DE experiments don't appear empty.
- `describe_experiment` response envelope includes `compartment`, `reports_fold_change`, `reports_derived_metric_types`, `abundance_gene_count`, `derived_metric_gene_count`, and lists child analysis nodes (AbundanceAnalysis / DerivedMetric, the latter grouped by value_kind).
- `gene_overview` (routing tool) reads new Gene rollups (`abundance_analysis_count`, `numeric_metric_count`, `classifier_flag_count`, `classifier_label_count`, `compartments_observed`, `numeric_metric_types_observed`, `classifier_flag_types_observed`, `classifier_label_types_observed`) so it knows which of the four detail tools to recommend without first hitting them.
- `describe_gene` response envelope surfaces counts: `n_abundance_analyses`, `n_numeric_metrics`, `n_classifier_flags`, `n_classifier_labels`, `n_de_edges`, partitioned by `compartment` / `metric_type` (read from Gene rollups; deeper breakdowns query-time).

No existing gene-scoped DE tools require changes (this spec does not modify `Changes_expression_of`).

## Risks and open questions

- **PDF extraction for Waldbauer 2012** required to harvest its derived metrics. Approach: manual transcription of Tables S1/S2 into CSV, or pdfplumber / tabula / OCR. Scope for a single paper; tooling shouldn't block the spec.
- **Protein ID → Gene resolution for proteomics abundance tables** (biller 2022 S3, Oleza 2015 polypeptide list). Existing v2 gene-id mapping handles `uniprot_accession`, `protein_id_refseq`, `locus_tag`, etc. Adapter reuses this path. Proteins that don't resolve to a Gene are dropped at ingest (mirrors current proteomics DE behavior).
- **Adapter location** — new edge shapes could live in the existing `omics_adapter.py` (grows large) or a sibling `observations_adapter.py` (mirrors the cluster_adapter separation). TBD during implementation; does not affect schema.
- **Validator boolean-token dry-run cost.** Boolean entry validation reads the source CSV to enumerate distinct cell tokens in `value_col` (rejects the paperconfig if any token isn't classified). For Biller-2018-sized tables (<1K rows), this is sub-second. For larger CSVs (>100K rows) the read could be noticeable; if it ever becomes a concern, the validator can sample (e.g., 10K rows) and warn instead of fail. Not addressed in v1 since current backlog tables are small.
- **DE retrofit asymmetry.** Until a follow-up spec retrofits `Changes_expression_of` to the analysis-node pattern, queries that blend DE + non-DE evidence (query pattern 5 above) must handle two shapes. Documented.
- **Compartment semantics for compartment-contrast DE.** Hypothetical "vesicle-vs-whole_cell enrichment" DE experiments would need a compartment choice (the single-compartment-per-Experiment invariant forces one). Sidestepped in v1 since biller S4 enrichment isn't in backlog. Revisit if a future paper requires it.
- **Future paired-modality values.** `PAIRED_RNASEQ_PROTEOME` names one specific pairing. Future papers pairing other modalities will add explicit values (e.g., `PAIRED_MICROARRAY_PROTEOME`, `PAIRED_RNASEQ_METABOLOMICS`). Keeping the pair explicit avoids an opaque generic bucket.

## Success criteria

**Spec-close requirements** (this spec's implementation is done when all of these are green):

- Schema changes merged in `config/schema_config.yaml`.
- Adapter emits new nodes + edges when fed a synthetic paperconfig fixture (not a real paper — real-paper integration is follow-up).
- Validator accepts all three new entry types and rejects each documented misuse (bad value_kind, unexpected token, out-of-allowed-categories value, mixed value_kind fields on the wrong entry, redundant per-timepoint abundance).
- Gene-id-mapping and paper-id-resolution code recognize the new entry types (harvests `id_columns`, writes `_resolved.csv` when applicable).
- `scripts/post-import.sh` + `scripts/post-import.cypher` are byte-identical on the DE-only subset against the pre-spec dump (no regression on existing computed properties); new blocks execute successfully on a fixture paperconfig.
- Full build pipeline (`uv run python create_knowledge_graph.py` + docker compose up) succeeds with no real paperconfig using the new entry types (baseline KG unchanged), and with a fixture paperconfig added (exercises all new nodes/edges/rollups).
- `/omics-edge-snapshot` reports zero `Changes_expression_of` regression; new edge type counts > 0 when fixture is active, 0 otherwise.
- KG validity tests updated with the assertions below; pass green.
- `docs/kg-changes/non-de-evidence-extension.md` written and self-contained for the explorer / MCP team.
- Skill docs (`.claude/skills/paperconfig/SKILL.md`, `/cypher-queries`, `/omics-edge-snapshot`) updated.

**Out of scope at spec close** (follow-up tasks):
- Biller 2022, Oleza 2015/2017, Waldbauer 2012, zinser 2009, Biller 2018 retrofit — all deferred to per-paper follow-up tasks.
- MCP tool surface changes — owned by the MCP team, unblocked by the "what changed" doc.

### KG validity assertions (part of spec-close requirements)

The KG validity test suite must assert all of the following. These are the concrete checks behind the "KG validity tests updated" spec-close bullet above.

- Every `abundance_analysis_measures_gene`, `derived_metric_quantifies_gene`, `derived_metric_flags_gene`, and `derived_metric_classifies_gene` edge targets a `:Gene` node.
- Every `abundance_analysis_measures_gene` edge has non-null `value`, `rank_in_sample`, `abundance_percentile`, and `abundance_bucket` (detection-only evidence is on `derived_metric_flags_gene`, not abundance edges).
- Every Experiment has exactly one `compartment` value.
- Every Experiment has exactly one `omics_type` value.
- Every Experiment with `is_time_course=true` has child analyses describing a single ordered timepoint sequence.
- `value_sd` is null whenever `n_replicates < 2` (or when paper doesn't report spread).
- No `abundance_analysis_measures_gene` edge duplicates a `Changes_expression_of` edge at the same `(Experiment, Gene, time_point)` triple. (DerivedMetric edges are scalar summaries without a time_point — not checked here.)
- `rank_by_metric`, `metric_percentile`, and `metric_bucket` are non-null iff the parent `DerivedMetric.rankable="true"` (only meaningful on `derived_metric_quantifies_gene`; the fields don't exist on flags/classifies edges).
- `metric_bucket` matches the pinned thresholds: `top_decile` iff `metric_percentile ≥ 90`; `top_quartile` iff `75 ≤ metric_percentile < 90`; `mid` iff `25 ≤ metric_percentile < 75`; `low` iff `metric_percentile < 25`. Same vocabulary + thresholds as `abundance_bucket`.
- `p_value` and `adjusted_p_value` on `derived_metric_quantifies_gene` edges are non-null only when parent `DerivedMetric.has_p_value="true"`. (Per-row null is allowed within a `has_p_value="true"` analysis when the paper's table has missing entries.) These properties do not exist on `derived_metric_flags_gene` or `derived_metric_classifies_gene`.
- `significant` on `derived_metric_quantifies_gene` is non-null iff `DerivedMetric.has_p_value="true"` AND `DerivedMetric.p_value_threshold IS NOT NULL` AND the edge's `adjusted_p_value` is non-null. When non-null, ∈ `{"true", "false"}` and equals `"true"` iff `adjusted_p_value < p_value_threshold`. Property does not exist on flags/classifies edges (significance is encoded directly in the value).
- Each `DerivedMetric` has outgoing edges of exactly one of the three derived-metric edge types, matching its `value_kind`.
- Every `derived_metric_classifies_gene` edge's `value_text` is a member of the parent `DerivedMetric.allowed_categories`.
- `value` is non-null on every `derived_metric_quantifies_gene` edge; `value_flag` is non-null and ∈ `{"true","false"}` on every `derived_metric_flags_gene` edge; `value_text` is non-null on every `derived_metric_classifies_gene` edge. (No null sentinel — the edge wouldn't exist if there were no measurement.)
- Each of the six new binding edges (Publication/Experiment/Organism × AbundanceAnalysis/DerivedMetric) has cardinality 1:many from parent to analysis (each analysis has exactly one Publication parent, one Experiment parent, one Organism parent — invariant #9 enforces the denormalized fields match).
- All boolean-semantic properties are stored as string enums `{"true", "false"}`, never as bool: `DerivedMetric.rankable`, `DerivedMetric.has_p_value`, `derived_metric_quantifies_gene.significant`, `derived_metric_flags_gene.value_flag`, `Experiment.reports_fold_change`. KG validity tests assert each property's value is always in its enum.
- `abundance_bucket` matches the pinned thresholds on every abundance edge: `top_decile` iff `abundance_percentile ≥ 90`; `top_quartile` iff `75 ≤ abundance_percentile < 90`; `mid` iff `25 ≤ abundance_percentile < 75`; `low` iff `abundance_percentile < 25`. Never null.
- Denormalized fields on `AbundanceAnalysis` and `DerivedMetric` (`compartment`, `omics_type`, `treatment_type`, `background_factors`, `publication_doi`, `organism_name`) equal their parent Experiment's values.
- Edge-denormalized fields (`abundance_analysis_measures_gene.value_type`, `derived_metric_quantifies_gene.metric_type`) equal their parent analysis node's value.
- `Experiment.abundance_gene_count` equals `count(DISTINCT g)` over `(exp)-[:experiment_has_abundance_analysis]->(:AbundanceAnalysis)-[:abundance_analysis_measures_gene]->(g:Gene)`; same for `derived_metric_gene_count`.
- `Experiment.abundance_analysis_count` / `derived_metric_count` / `derived_metric_value_kinds` are consistent with direct query-time counts of child analysis nodes.
- Parallel Publication rollups (`Publication.abundance_analysis_count`, `Publication.derived_metric_count`, `Publication.abundance_gene_count`, `Publication.derived_metric_gene_count`, `Publication.compartments`, `Publication.derived_metric_types`, `Publication.derived_metric_value_kinds`) are consistent with direct query-time counts over the corresponding 1-hop binding edges.
- Parallel OrganismTaxon rollups (same seven fields) are consistent with direct query-time counts over the corresponding 1-hop binding edges.
- Existing DE time-point arrays (`time_point_labels`, `time_point_totals`, `time_point_significant_up`, `time_point_significant_down`) reflect DE edges only — independent of any AbundanceAnalysis children.
- `abundance_time_point_*` arrays are parallel (same length, ordered by `time_point_order`); `abundance_time_point_totals[i]` equals the count of `abundance_analysis_measures_gene` edges from the AbundanceAnalysis at `abundance_time_point_orders[i]`.
- For timeseries experiments with both DE and AbundanceAnalysis children, no two distinct timepoint labels share the same `time_point_order` across the union of DE and abundance arrays.
- Gene routing counts (`abundance_analysis_count`, `numeric_metric_count`, `classifier_flag_count`, `classifier_label_count`, `compartments_observed`, `numeric_metric_types_observed`, `classifier_flag_types_observed`, `classifier_label_types_observed`) are consistent with a direct query-time count for a sampled set of genes.
- **Empty-state defaults**: all new int rollup properties (`abundance_analysis_count`, `derived_metric_count`, `*_gene_count` on Experiment/Publication/OrganismTaxon, and Gene routing counts) default to `0` — never null. All new str[] rollup properties (`derived_metric_value_kinds`, `compartments`, `derived_metric_types` on Publication/OrganismTaxon, `compartments_observed`, `*_types_observed` on Gene) default to `[]` — never null. Post-import MUST run a "set defaults on unvisited nodes" pass before any `list_*` tool could observe them.
- `Changes_expression_of` edge counts are unchanged by this spec (verified by `/omics-edge-snapshot` before/after).

---

## Appendix — Deferred: Metabolite layer (own spec)

These notes are harvested so the follow-up spec doesn't start from zero. **Nothing in this appendix is in scope for the current spec.**

### Why deferred

- **ID normalization non-trivial.** Papers reference compounds via KEGG compound, ChEBI, HMDB, trivial name, m/z + retention_time tuples (untargeted metabolomics). Structural deduplication (InChIKey) requires chemical reasoning.
- **Biochemistry linkage requires its own pipeline.** `(Gene)-[:catalyzes_reaction {role: substrate|product}]->(Metabolite)` via KEGG reactions + existing EC/KO on genes. Needs reaction download, substrate/product normalization, generic-vs-specific handling.
- **Untargeted metabolomics features without identities** (biller 2022 Table S6: ~4,064 mass features, most unidentified) need their own representation — feature-as-node vs. skip-until-identified.
- **Mass-spec-specific metadata** (analytical fraction, m/z, retention time, MS-DIAL parameters) adds schema surface orthogonal to the gene-evidence design here.

### Sketch of future Metabolite node

```
(:Metabolite {
  id: str,                           // bioregistry CURIE: "kegg.compound:C00001" or "chebi:15377"
  name: str,
  name_synonyms: str[],
  formula: str | null,
  monoisotopic_mass: float | null,   // Da; unambiguous (not average mass)
  inchi_key: str | null              // structural deduplication
})
```

### Sketch of future evidence edges

Under the analysis-node pattern established by this spec, metabolite evidence would likely use:
- `(Experiment) → (AbundanceAnalysis {..., entity_type: metabolite}) → (Metabolite)` for per-sample metabolite abundance, OR a parallel `MetaboliteAbundanceAnalysis` node if entity-type-specific metadata diverges.
- Metabolite DE: target-widening of `Changes_expression_of` to Metabolite, *OR* a new `Changes_abundance_of` edge, *OR* a metabolite-specific analysis node. The choice affects whether existing Gene-scoped MCP queries need a `target_kind` filter (the breaking-change issue the explorer flagged). Belongs in the metabolite spec, not here.
- `(Gene)-[:catalyzes_reaction {role: substrate|product, reaction_id: str}]->(Metabolite)` — biochemistry layer.

### Backlog papers awaiting this layer

- **Kujawinski 2023** — metabolite diversity across Prochlorococcus ecotypes.
- **biller 2022 Tables S1, S5, S6, S7** — lipids, targeted + untargeted metabolites.
- **Capovilla 2023 Table S3** — metabolomics under chitosan.
- **Ma 2022** — any metabolite content.

### Analytical value

"Nutrients uptaken / released / exchanged" — the motivating query #2 from brainstorming — most depends on metabolite nodes + biochemistry edges. Current KG answers it only via transporter gene expression; metabolite identity closes the loop.

### Controlled vocabularies likely needed

- Metabolite ID type: `kegg_compound`, `chebi`, `hmdb`, `inchi_key`, `mass_feature`.
- Analytical fraction: `rp_organic_positive`, `rp_aqueous_positive`, `hilic_positive`, `hilic_negative` (biller 2022 fractions).
- Additional `value_type`: `peak_area`, `peak_area_biovolume_norm`, `concentration_um` (+ unit handling).

### Open design questions

- Feature-as-node (store every unidentified mass feature)? Or skip until identified?
- Metabolite grouping (e.g., lipid class) as separate nodes or properties?
- Scope of biochemistry edges: all KEGG reactions, or only those where ≥1 gene in the KG catalyzes?
- Reuse `AbundanceAnalysis` for metabolites (target widening) or parallel node type?
