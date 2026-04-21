---
name: paperconfig
description: Create and validate paperconfig.yaml files for integrating omics publication data (RNA-seq, proteomics, metabolomics, microarray) into the knowledge graph. Use when adding a new paper's differential expression data.
argument-hint: [paper-directory-name]
user-invocable: true
allowed-tools: Read, Grep, Glob, Write, Edit, Bash(uv *), Bash(python *)
---

# Paperconfig YAML Instructions

How to create a `paperconfig.yaml` file for integrating omics publication data into the knowledge graph.

## Overview

Each publication with differential expression or omics data gets its own `paperconfig.yaml` file under `data/Prochlorococcus/papers_and_supp/<AuthorName Year>/`. This config tells the OMICSAdapter how to parse supplementary data tables and create knowledge graph nodes and edges.

## Directory Setup

1. Create a folder: `data/Prochlorococcus/papers_and_supp/<AuthorName Year>/`
2. Place the publication PDF in this folder
3. Place the supplementary data CSV file(s) in this folder (or a subfolder)
4. Create `paperconfig.yaml` in the same folder

## YAML Structure

A paperconfig has three top-level sections under `publication:`:

```yaml
publication:
  papername: "<Author Year>"
  doi: "10.NNNN/xxxxx"  # optional — overrides PDF-extracted DOI
  papermainpdf: "data/Prochlorococcus/papers_and_supp/<AuthorName Year>/<filename>.pdf"

  experiments:
    # ... (required — defines experiment-level metadata)

  supplementary_materials:
    # ... (required — defines data tables and per-timepoint analyses)
```

### 1. Publication Metadata (Required)

```yaml
publication:
  papername: "Author Year"
  doi: "10.NNNN/xxxxx"  # optional — overrides PDF-extracted DOI
  papermainpdf: "data/Prochlorococcus/papers_and_supp/<folder>/<pdf_filename>.pdf"
```

- `papername`: Short citation label (e.g., "Biller 2016"). Used for display only.
- `doi`: *(optional)* Explicit DOI override. When present, takes precedence over the PDF-extracted DOI for the Publication node ID and all Experiment node IDs. Must match pattern `10.NNNN/...`. Use this when the PDF extractor fails to find the correct DOI.
- `papermainpdf`: Relative path to the publication PDF. The adapter extracts metadata (title, DOI, authors, abstract) from this PDF automatically using an LLM-based extractor.

### 2. Experiments Block (Required)

The `experiments` block defines experiment-level metadata that is shared across all analyses (timepoints) within that experiment. Each experiment entry becomes the basis for `Experiment` nodes and `Changes_expression_of` edges in the knowledge graph.

```yaml
  experiments:
    <experiment_key>:
      name: "<human-readable experiment description>"
      organism: "Prochlorococcus MED4"
      omics_type: RNASEQ              # RNASEQ | PROTEOMICS | METABOLOMICS | MICROARRAY
      test_type: DESeq2               # statistical method
      treatment_type: [coculture]     # canonical treatment category list (see below)
      background_factors: []          # experimental context factors list (see below, optional)
      compartment: "whole_cell"       # optional; default "whole_cell". Values: whole_cell | vesicle | exoproteome | secretome
                                      # Different compartments from the same paper MUST be split into separate Experiments.
      treatment_condition: "Coculture with Alteromonas HOT1A3"
      control_condition: "Axenic"
      experimental_context: "in Pro99 medium under continuous light"
      table_scope: all_detected_genes   # all_detected_genes | significant_any_timepoint | significant_only | top_n | filtered_subset
      # table_scope_detail: ""          # free-text clarification (optional, use with filtered_subset)
      # Optional fields:
      medium: "Pro99 natural seawater medium"
      temperature: "24C"
      light_condition: "continuous light"
      light_intensity: "55 umol photons m-2 s-1"
      # Coculture-specific (required when treatment_type is 'coculture' or 'viral'):
      treatment_organism: "Alteromonas macleodii HOT1A3"
      treatment_taxid: 28108
      treatment_assembly_accession: GCF_901457835.2  # optional; use for genome-loaded organisms
```

#### Required Experiment Fields

| Field | Description | Examples |
|-------|-------------|----------|
| `name` | Human-readable description of the experiment | `"DE of Prochlorococcus MED4 coculture vs axenic"` |
| `organism` | Organism being profiled | `"Prochlorococcus MED4"`, `"Alteromonas MIT1002"` |
| `omics_type` | Omics data type | `RNASEQ`, `MICROARRAY`, `PROTEOMICS`, `METABOLOMICS` |
| `test_type` | Statistical method used | `DESeq2`, `edgeR`, `Rockhopper`, `DESeq`, `microarray`, `microarray_Cyber-T`, `microarray_LPE`, `microarray_Goldenspike` |
| `treatment_type` | Canonical treatment category (list) | See **Treatment Types** below |
| `treatment_condition` | Description of the experimental/test condition | `"Coculture with Alteromonas HOT1A3"`, `"Salt-acclimated (5% salt)"` |
| `control_condition` | Description of the baseline/reference condition | `"Axenic"`, `"Normal seawater (3.8% salt)"` |

#### Optional Experiment Fields (Recommended)

| Field | Description | When to Use |
|-------|-------------|-------------|
| `experimental_context` | Other factors held constant | Always helpful; describes medium, light, temperature, etc. |
| `medium` | Growth medium | When known (e.g., `"Pro99 natural seawater medium"`) |
| `temperature` | Temperature | When known (e.g., `"24C"`) |
| `light_condition` | Light regime | When known (e.g., `"continuous light"`, `"13:11 light:dark cycle"`, `"constant darkness"`) |
| `light_intensity` | Light intensity | When known (e.g., `"55 umol photons m-2 s-1"`) |
| `table_scope` | What genes the source DE table contains | Always recommended. Values: `all_detected_genes`, `significant_any_timepoint`, `significant_only`, `top_n`, `filtered_subset` |
| `table_scope_detail` | Free-text clarification for `table_scope` | When `table_scope` is `filtered_subset` or ambiguous (e.g., `"Top 50% of genes by expression level"`) |
| `background_factors` | Experimental context factors not being compared in DE (list) | When conditions like coculture/axenic status or light regime are relevant but not the treatment variable. Same vocabulary as `treatment_type` (see Canonical Vocabulary table below) |

#### Formatting conventions

- **`light_condition`**: Use spaces, not underscores: `"continuous light"` (not `"continuous_light"`)
- **`light_intensity`**: Use ASCII `umol` (not Unicode `µmol`): `"55 umol photons m-2 s-1"`
- **`medium`**: Include treatment-specific modifications when relevant: `"Pro99 with no added iron"` (not just `"Pro99"`). For the control medium, describe what was used: `"Pro99 with 1 uM Fe total"`
- **`temperature`**: Use format `"24C"` (no degree symbol)
| `treatment_organism` | Name of the organism causing the effect | Required for coculture/viral experiments |
| `treatment_taxid` | NCBI Taxonomy ID of the treatment organism | Required for coculture/viral experiments |
| `treatment_assembly_accession` | NCBI RefSeq assembly accession | Use when the treatment organism has a loaded genome in the KG |

#### Canonical Vocabulary (`treatment_type` and `background_factors`)

`treatment_type` and `background_factors` share the same canonical vocabulary. `treatment_type` = "what environmental variable is being manipulated." `background_factors` = "what conditions are held constant." The specific values (e.g., "55 umol photons" for `light`) live in the `treatment_condition`, `control_condition`, and `light_condition` fields.

| Value | As `treatment_type` | As `background_factors` |
|---|---|---|
| `nitrogen` | N-limitation / N-starvation | N-replete medium |
| `phosphorus` | P-limitation / P-starvation | P-replete medium |
| `iron` | Fe-limitation | Fe-replete medium |
| `carbon` | CO2 / carbon source manipulation / chitosan | Fixed carbon source |
| `light` | Light quality or intensity manipulation | Continuous light regime |
| `darkness` | Extended dark treatment | Dark regime |
| `diel` | Diel light-dark cycling (circadian) | Diel light-dark regime |
| `temperature` | Thermal shift / acclimation | Fixed temperature |
| `salt` | Salinity / osmotic changes | Fixed salinity |
| `coculture` | Co-cultivation with another organism (NOT phage) | Coculture context |
| `viral` | Phage infection or vDOM exposure | Phage present |
| `chemical` | Chemical treatment (e.g., DCMU) | Chemical inhibitor present |
| `plastic` | Plastic leachate exposure | — |
| `growth_phase` | Growth state / multi-condition comparison | — |
| `mutant` | Mutant or evolved strain comparison | Mutant background |
| `axenic` | — (background only) | Pure culture, no other organisms |

**Rules for assigning `background_factors`:**
- **Axenic/coculture status:** Always in `background_factors` unless `coculture` IS the DE comparison (treatment_type). Non-coculture, non-viral experiments: `axenic` if no partner organism.
- **Viral experiments:** Do NOT add `axenic` or `coculture` to background_factors. Viral experiments infect otherwise-axenic cells by default. Only add `coculture` if the experiment is genuinely in a multi-species community.
- **Light regime:** Add `light` when under continuous light and treatment_type is not light-related. Add `diel` when under a light:dark cycle. Add `darkness` when in dark regime. Applies to ALL experiment types including coculture and viral.
- **Coculture experiments:** Add `light` or `diel` for the light regime, but do NOT add `axenic` or `coculture`.

**Examples:**
- Darkness experiment in coculture → `treatment_type: [darkness]`, `background_factors: [coculture, diel]`
- Phage infection under continuous light → `treatment_type: [viral]`, `background_factors: [light]`
- N-starvation, axenic, continuous light → `treatment_type: [nitrogen]`, `background_factors: [axenic, light]`
- Coculture vs axenic under continuous light → `treatment_type: [coculture]`, `background_factors: [light]`
- DCMU + light vs light → `treatment_type: [chemical]`, `background_factors: [axenic, light]`
- Multi-temperature acclimation → `treatment_type: [temperature]`, `control_condition: "Multi-temperature comparison (17-30C) — pairwise contrasts"`

#### Extending the Canonical Vocabulary

The canonical vocabulary is intentionally minimal. When creating a paperconfig for a new paper:

1. **Check existing values first.** Can this condition map to an existing category?
   - "UV exposure" → `light`
   - "CO2 enrichment" → `carbon`
   - "DCMU treatment" → `chemical`
2. **Prefer general categories over specific ones.** `chemical` is better than `dcmu_inhibitor`. The specifics go in `treatment_condition` and `experimental_context`.
3. **If no existing value fits**, flag it to the user: "This paper studies [X], which doesn't map cleanly to any existing canonical value. Closest match: `Y`. Should I use `Y` or propose a new value?"
4. **User decides.** If a new value is approved, update all three locations in a single commit:
   - `CANONICAL_CONDITION_TYPES` in `scripts/validate_paperconfig.py`
   - Canonical Vocabulary table in this SKILL.md
   - `CLAUDE.md` key graph facts section
5. **Never invent new canonical values without user approval.**

#### Grouping Analyses into Experiments

Each experiment represents **one biological question**: same organism, same treatment
type, same omics, same experimental context. Multiple timepoints of the same
comparison share one experiment. Key rules:

1. **Same biological trajectory = same experiment.** If analyses compare different
   timepoints of the same starvation/infection/treatment, they belong in one
   experiment. The `timepoint` field on each analysis captures the time variation.
   Example: P-limited at 4h, 24h, 46h, 59h → one experiment with 4 analyses.

2. **Rescue/re-addition = same experiment as the starvation.** If P is re-added
   after P-starvation, those analyses belong in the same experiment. The `timepoint`
   label captures the intervention (e.g., `"50h (P added)"`).

3. **Different inoculum density = different experiment.** If two analyses compare
   the same organisms but at different inoculum concentrations, they are separate
   experiments — put the density in `experimental_context`.

4. **Phase-varying treatment text = same experiment.** Papers like Weissberg 2025
   use different treatment labels per timepoint ("Nutrient starvation", "Long-term
   starvation", "Decline"). If these all describe the same biological trajectory
   (same environmental condition), group them into one experiment. Use the condition
   name (not the per-timepoint label) as `treatment_condition`.

5. **Always preserve `timepoint`.** Every analysis from a time-course experiment
   MUST have `timepoint` (the original string label from the paper) and
   `timepoint_hours` (numeric conversion or null if unparseable).

6. **Axenic vs coculture context = separate experiments.** If the same organism
   and treatment are measured in axenic and coculture conditions, those are
   separate experiments (different `experimental_context`).

7. **Phage = `viral`, not `coculture`.** When `treatment_organism` is Phage,
   use `treatment_type: viral`. Reserve `coculture` for bacterial co-cultivation.
   Note: if phage infection is the *context* but not the *treatment variable*
   (e.g., Lin 2015 studies P-limitation in phage-infected cells), use the
   actual treatment type (`phosphorus`) and describe phage in
   `experimental_context`.

8. **Use the most specific `treatment_type`.** Prefer the specific environmental
   variable over generic labels. `growth_phase` is a last resort when no
   specific category fits. Examples:
   - N starvation time course → `nitrogen`
   - Alternative N sources (cyanate, urea) vs ammonium → `nitrogen`
     (still nitrogen metabolism, even without deprivation)
   - Iron rescue (return to replete after starvation) → `iron`
     (part of the same iron experiment)
   - Nutrient starvation growth phases (exponential → decline → death) →
     `nitrogen` (if N is the limiting nutrient)
   - Dark vs light → `light` or `darkness` depending on which is treatment
   - DCMU or other inhibitor → `chemical`
   - Use `growth_phase` only for growth phase comparisons without a specific stressor

### 3. Supplementary Materials (Required)

Each supplementary table is a keyed entry under `supplementary_materials`. Five entry types are supported:

| Type | Purpose |
|---|---|
| `csv` | Differential expression results table with `statistical_analyses` |
| `id_translation` | Pure ID mapping (no DE data); bridges non-standard IDs to locus tags |
| `annotation_gff` | GFF3/GTF file; adds protein_id/Name bridges for a strain |
| `gene_clusters` | Cluster assignment table; creates ClusteringAnalysis + GeneCluster nodes |
| `derived_metrics_table` | Column-level scalar summaries (periodicity flags, classifiers, numeric scores); creates DerivedMetric nodes + one of 3 measurement edge types per metric |


#### Type `csv` -- Expression data table (required for omics edges)

```yaml
  supplementary_materials:
    <table_key>:
      type: csv
      filename: "data/Prochlorococcus/papers_and_supp/<folder>/<data_file>.csv"
      sep: ","                  # optional; column delimiter, default "," (use "\t" for TSV)
      organism: "Prochlorococcus MED4"   # optional; inferred from experiment if absent
      original_filename: "data/.../original.csv"   # optional; pre-fix-gene-ids CSV for reference

      # optional; if absent, heuristic column detection is used (backward compatible)
      id_columns:
        - column: "NCBI ID"
          id_type: locus_tag_ncbi
        - column: "Gene Name"
          id_type: gene_name
      product_columns:           # optional; columns containing functional descriptions
        - column: "Gene description"

      statistical_analyses:
        - <analysis_1>
        - <analysis_2>
        # ... multiple analyses can reference the same file
```

- `<table_key>`: A descriptive key for the table (e.g., `supp_table_1`, `supp_table_2_timepoint_A`).
- `filename`: Relative path to the data CSV file.
- `sep`: Column delimiter; default `","`. Use `"\t"` for tab-separated files.
- `organism`: Organism for this table; if absent, inferred from the experiment referenced by the analysis.
- `original_filename`: Path to the original (pre-fix-gene-ids) CSV. Useful when `filename` was updated to a `_with_locus_tag.csv` by fix-gene-ids -- the original is still available for ID column inspection.
- `id_columns`: Documents which columns contain gene identifiers and their type (see **ID Types** below). If absent, heuristic detection is used.
- `product_columns`: Documents which columns contain functional descriptions to harvest as product synonyms.
- Multiple analyses can share the same file if different columns represent different comparisons.

#### Type `id_translation` -- Pure ID mapping resource (no expression data)

Use when a paper includes a supplementary table that maps gene IDs without DE results (e.g., a genome annotation CSV, a cross-reference table). The build script (`build_gene_id_mapping.py`) uses this to enrich the per-strain ID lookup before processing DE CSVs. **Not processed by omics_adapter -- no expression edges are emitted.**

```yaml
    <table_key>:
      type: id_translation
      filename: "data/Prochlorococcus/papers_and_supp/<folder>/<id_table>.csv"
      sep: ","           # optional; default ","
      organism: "Prochlorococcus MIT9301"   # required
      id_columns:        # required
        - column: "NCBI locus tag"
          id_type: locus_tag_ncbi
        - column: "Old locus tag"
          id_type: old_locus_tag
        - column: "Gene name"
          id_type: gene_name
        - column: "JGI ID"
          id_type: jgi_id
      product_columns:   # optional
        - column: "Gene product"
```

- `organism` and `id_columns` are **required** for this type.
- List `id_translation` entries **before** `csv` entries in the same paper so JGI IDs / alt-IDs are in the lookup when the DE table is processed.

#### Type `annotation_gff` -- Paper-specific GFF/GTF reannotation (ID bridging only)

Use when a paper includes a GFF3 or GTF file from a custom genome reannotation. The build script extracts locus tags, old locus tags, gene names, and protein IDs to bridge paper-specific IDs to canonical locus tags. **Not processed by omics_adapter and skipped by fix-gene-ids / check-gene-ids.**

```yaml
    <table_key>:
      type: annotation_gff
      filename: "data/Prochlorococcus/papers_and_supp/<folder>/<reannotation>.gff"
      organism: "Prochlorococcus MED4"   # required
```

- `organism` is **required** for this type.
- Both GFF3 (`.gff`, `.gff3`) and GTF (`.gtf`) formats are supported.

### `type: gene_clusters`

Co-expression cluster membership tables. Processed by `cluster_adapter.py` (NOT by `omics_adapter`).

**Entry key** must be short, meaningful, and unique within the paper (used in graph node IDs):
- Good: `med4_kmeans_nstarvation`, `mit9313_mfuzz_diel`
- Bad: `cluster_table_1`, `supp_table_clusters`

**Required fields (analysis level):**

| Field | Type | Description |
|---|---|---|
| `name` | str | Human-readable label for the clustering analysis (e.g., "MED4 K-means N-starvation clusters") |
| `filename` | str | Path to cluster membership CSV |
| `organism` | str | Target organism (canonical name) |
| `gene_id_col` | str | CSV column with gene identifiers |
| `cluster_col` | str | CSV column with cluster assignment |
| `cluster_type` | str | `diel_cycle` \| `time_series_dynamics` \| `response_pattern` |

**Optional fields (analysis level):**

| Field | Type | Description |
|---|---|---|
| `score_col` | str | CSV column with fuzzy membership score |
| `cluster_method` | str | Algorithm (e.g., "K-means (K=9)", "Mfuzz") |
| `omics_type` | str | MICROARRAY, RNASEQ, PROTEOMICS, METABOLOMICS |
| `light_condition` | str | Light regime (e.g., "14:10 L:D") |
| `treatment_type` | str[] | Array of canonical treatment types |
| `background_factors` | str[] | Array of background condition factors (same vocabulary as `treatment_type`) |
| `treatment` | str | Experiment description |
| `experimental_context` | str | Additional context |
| `experiments` | str[] | **Recommended.** List of experiment keys from the same paperconfig (links ClusteringAnalysis → Experiment). Every cluster analysis should link to at least one experiment. Create an experiment entry if one doesn't exist yet. A cluster analysis may link to multiple experiments if relevant. |

**Per-cluster data** comes from extraction JSON (`cluster_extraction_{entry_key}.json`), not from the paperconfig. The extraction pipeline reads cluster keys from the CSV `cluster_col` and produces: `id`, `name`, `functional_description`, `temporal_pattern`, `expression_dynamics`.

**Data flow:** paperconfig (analysis metadata) + CSV (cluster membership) → extraction pipeline → extraction JSON (per-cluster descriptions) → adapter reads all three.

**Example:**

```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"
  filename: "data/.../med4_kmeans_clusters.csv"
  organism: "Prochlorococcus MED4"
  gene_id_col: "gene_id"
  cluster_col: "cluster"
  cluster_type: "response_pattern"
  cluster_method: "K-means (K=9)"
  omics_type: MICROARRAY
  treatment_type: ["nitrogen"]
  background_factors: ["axenic", "light"]
  treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
  light_condition: "continuous light"
  experimental_context: "Custom Affymetrix microarray..."
  experiments: [n_starvation_med4]
```

**Notes:**
- Gene IDs go through the same step 4 resolution pipeline as DE tables
- For papers with separate clusters per organism, use separate `type: gene_clusters` entries (one per organism/analysis)
- `treatment_type` must be an array (same enum as experiments)

### `type: derived_metrics_table`

Column-level scalar summaries per gene (periodicity flags, classifiers, numeric scores). Processed by `observations_adapter.py`. Creates one DerivedMetric node per metric and emits one of 3 measurement edge types to Gene based on `value_kind`.

**Entry key** must be short, meaningful, unique within the paper (used in graph node IDs as `derived_metric:{doi_short}:{entry_key}:{metric_type}`).

**Required fields (entry level):**

| Field | Type | Description |
|---|---|---|
| `filename` | str | Path to metrics CSV |
| `organism` | str | Target organism (canonical name) |
| `experiment` | str | Key into the `experiments:` block — parent Experiment for all metrics in this entry |
| `name_col` | str | CSV column with gene identifiers (usually `locus_tag`) |
| `id_columns` | list | Maps `name_col` to `id_type` (see ID Types below) |
| `metrics` | list | One or more metric definitions |

**Per-metric required fields:**

| Field | Type | Description |
|---|---|---|
| `metric_type` | str | Short metric identifier (e.g. `periodic_in_axenic_LD`, `darkness_survival_class`, `fourier_score`) |
| `value_kind` | str | `numeric` \| `boolean` \| `categorical` — determines which edge type is emitted |
| `value_col` | str | CSV column holding the metric value |
| `field_description` | str | Free-text description of what the metric measures |

**Per-metric optional fields:**

| Field | Type | Description |
|---|---|---|
| `name` | str | Human-readable metric label (defaults to `metric_type` if absent) |
| `unit` | str | Unit of measurement (numeric only) |
| `blank_policy` | str | **REQUIRED when `value_kind: boolean`.** `skip` \| `false` \| `true` — how to treat blank/NaN cells |
| `allowed_categories` | str[] | **REQUIRED when `value_kind: categorical`.** Validator rejects rows with values outside this list |
| `rankable` | str | `"true"` \| `"false"` (**string**, never bool) — numeric only; governs post-import `rank_by_metric` / `metric_percentile` / `metric_bucket` |
| `has_p_value` | str | `"true"` \| `"false"` (**string**) — numeric only; enables `significant` derivation |
| `p_value_col` | str | CSV column for raw p-values (numeric + `has_p_value="true"`) |
| `adjusted_p_value_col` | str | CSV column for adjusted p-values (numeric + `has_p_value="true"`) |
| `p_value_threshold` | float | Cutoff for `significant` derivation (numeric + `has_p_value="true"`) |

**Emitted edges:**

| `value_kind` | Edge type | Edge properties |
|---|---|---|
| `numeric` | `Derived_metric_quantifies_gene` | `value`, `p_value`, `adjusted_p_value`, plus post-import `rank_by_metric` / `metric_percentile` / `metric_bucket` (if rankable) / `significant` (if has_p_value) |
| `boolean` | `Derived_metric_flags_gene` | `value_flag ∈ {"true","false"}` |
| `categorical` | `Derived_metric_classifies_gene` | `value_text` (must match `allowed_categories`) |

**Examples:**

```yaml
# BOOLEAN: periodicity flags per gene
s4a_natl2a_axenic:
  type: derived_metrics_table
  filename: "data/.../Biller 2018/table_s4a.csv"
  organism: "Prochlorococcus NATL2A"
  experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
  name_col: locus_tag
  id_columns:
    - column: locus_tag
      id_type: locus_tag
  metrics:
    - metric_type: periodic_in_axenic_LD
      name: "Periodic in NATL2A axenic L:D (Table S4A)"
      value_kind: boolean
      value_col: periodic_in_LD
      blank_policy: skip
      field_description: "Boolean flag: periodic under 12:12 L:D in axenic NATL2A"

# CATEGORICAL: survival classes
s5_natl2a_survival:
  type: derived_metrics_table
  filename: "data/.../Biller 2018/table_s5.csv"
  organism: "Prochlorococcus NATL2A"
  experiment: darkness_extended_darkness_natl2a_rnaseq_axenic
  name_col: locus_tag
  id_columns:
    - column: locus_tag
      id_type: locus_tag
  metrics:
    - metric_type: darkness_survival_class
      name: "NATL2A darkness survival class (Table S5)"
      value_kind: categorical
      value_col: survival_class
      allowed_categories:
        - "strongly_upregulated"
        - "moderately_upregulated"
        - "unchanged"
        - "moderately_downregulated"
        - "strongly_downregulated"
      field_description: "Five-level survival-class assignment"

# NUMERIC: rankable periodicity score with p-value
fourier_metrics:
  type: derived_metrics_table
  filename: "data/.../fourier_scores.csv"
  organism: "Prochlorococcus MED4"
  experiment: diel_rnaseq_med4
  name_col: locus_tag
  id_columns:
    - column: locus_tag
      id_type: locus_tag
  metrics:
    - metric_type: fourier_score
      value_kind: numeric
      value_col: fourier
      unit: ""
      rankable: "true"
      has_p_value: "true"
      p_value_col: p_value
      adjusted_p_value_col: adj_p_value
      p_value_threshold: 0.05
      field_description: "Fourier-transform periodicity score"
```

**Notes:**
- Each DerivedMetric emits exactly ONE edge type — determined by `value_kind`. Mixing edge types under one DM would be a bug.
- Gene IDs go through the same step 4 resolution pipeline as DE tables.
- `allowed_categories` is mandatory for categorical metrics. Rows with values outside the list are skipped (and logged).
- For papers with paired-modality metrics (e.g., transcript × protein lag coefficients), link all metrics to a single `PAIRED_RNASEQ_PROTEOME` Experiment rather than to individual source-modality experiments.

#### ID Types (`id_type` values)

| `id_type` | Description | Example |
|-----------|-------------|---------|
| `locus_tag` | Canonical locus tag | `PMM0001`, `P9301_RS09095` |
| `locus_tag_ncbi` | NCBI RefSeq RS-format locus tag | `TX50_RS00020` |
| `locus_tag_cyanorak` | Cyanorak cluster-specific locus tag | `CK_Pro_MED4_00001` |
| `old_locus_tag` | Legacy locus tag (deprecated NCBI format) | `PMM0001`, `P9301_05911` |
| `alternative_locus_tag` | Alt tag from a different annotation | `PMED4_00071` |
| `gene_name` | Gene symbol (generic; shared across organisms) | `dnaA`, `dnaN` |
| `gene_synonym` | Alternative gene symbol | `beta-clamp` |
| `protein_id_refseq` | RefSeq WP_ protein accession | `WP_011131639.1` |
| `uniprot_accession` | UniProt accession | `Q7V6L1` |
| `uniprot_entry_name` | UniProt entry name (build script strips `_ORGANISM` suffix) | `DNAA_PROM0` |
| `jgi_id` | JGI IMG gene catalog ID (integer string) | `2626311743` |
| `probeset` | Microarray probeset ID | `MED4_ARR_0008_x_at` |
| `rast_id` | RAST annotation FIG ID | `fig\|59919.17.peg.1` |
| `annotation_specific` | ID unique within a non-standard annotation system | |
| `other` | Any other identifier type | |

#### Anti-patterns — do NOT declare free-text columns as `locus_tag`

**Never** declare a free-text column (product descriptions, UniProt OS strings, narrative text) with `id_type: locus_tag`. The `build_gene_id_mapping` builder tokenizes every value in the column and registers each token as an alt_id in the strain's `specific_lookup`. A single paper can pollute the strain's `gene_id_mapping.json` with thousands of junk alt_ids — every subsequent row match collapses to the one anchor gene the builder happened to pick up first.

Symptoms seen in practice (Domínguez 2017 + Moreno 2023 S2 BL107 + Moreno 2023 S4 Marinobacter, 2026-04-15):
- Paperconfig ends up with a handful of `Changes_expression_of` edges instead of hundreds
- `specific_lookup.values()` has one anchor locus_tag with thousands of alt_ids (e.g. 2,148 keys all mapping to `Pro1395` in SS120)
- The junk keys include whole sentences, punctuation (`'/'`), and individual words (`'1'`, `'10'`)

Detection:

```bash
# any paperconfig that declares a likely-free-text column as locus_tag:
grep -A1 -E 'column: "(Description|description|product|Product|Name|name|Protein|Gene Name)"' \
    data/*/papers_and_supp/*/paperconfig.yaml | grep 'id_type: locus_tag'

# pollution in an existing gene_id_mapping:
python3 -c "
import json, collections
m = json.load(open('cache/data/<Organism>/genomes/<Strain>/gene_id_mapping.json'))
print(collections.Counter(m['specific_lookup'].values()).most_common(3))"
# healthy: top-3 anchors have <~20 alt_ids each; polluted: one anchor has thousands
```

The `scripts/validate_paperconfig.py` validator warns when any of the known free-text column names is declared `id_type: locus_tag` (but the warning is advisory — it won't prevent a broken config from validating).

Correct pattern — extract the locus tag into a clean column via the `_modified.csv` builder:

```python
# scripts/build_modified_csv/build_<paper>_modified_csv.py
import re
GN_RE = re.compile(r"GN=(\S+)")
df["extracted_gn"] = df["Description"].fillna("").map(
    lambda s: (m := GN_RE.search(s)) and m.group(1) or ""
)
df.to_csv(src.with_name(src.stem + "_modified.csv"), index=False)
```

```yaml
# paperconfig points at _modified.csv, uses the extracted column
supp_table_1:
  type: csv
  filename: "...table s3 Combined_modified.csv"
  organism: "..."
  id_columns:
    - column: "extracted_gn"
      id_type: gene_name       # Tier 3 — safe for mixed locus tag / gene symbol content
    - column: "Accession"
      id_type: uniprot_entry_name
  product_columns:
    - column: "Description"    # narrative-only; product_columns don't feed build_gene_id_mapping
  statistical_analyses:
    - id: "..."
      name_col: "Accession"    # NOT "Description"
      logfc_col: "log2_fold_change"
      ...
```

**When the extracted column has mixed content** (locus tags like `Pro_1040` plus gene symbols like `atpB`), declare it `id_type: gene_name` (Tier 3 multi_lookup, only used when a value singleton-resolves). Use `id_type: locus_tag` only when every value in the column is a real locus tag.

See `scripts/build_modified_csv/README.md` for the builder convention and `docs/community_proteomics_marref_saga.md` for the historical context that surfaced this anti-pattern.

### 4. Statistical Analysis Fields

Each analysis entry describes one **per-timepoint** differential expression comparison within an experiment. Experiment-level metadata (organism, omics_type, test_type, conditions, etc.) is defined in the `experiments` block; analyses reference an experiment and add only per-timepoint details.

#### Required Fields

| Field | Description | Examples |
|-------|-------------|----------|
| `id` | Unique identifier for this analysis | `DE_coculture_vs_axenic_med4_20h` |
| `experiment` | Reference to a key in the `experiments` block | `coculture_hot1a3` |
| `name_col` | Column name in CSV containing the gene/protein identifier | `"Gene"`, `"Locus tag2"`, `"Gene ID"`, `"Synonym"` |
| `logfc_col` | Column name in CSV containing the log2 fold change value | `"log2FoldChange"`, `"logFC"` |

#### Timepoint and Growth Phase Fields

Every analysis must have `timepoint`, `timepoint_hours`, and `growth_phase`. Look for these first in the paperconfig's own free-text fields (`treatment_condition`, `control_condition`, `experimental_context`, `treatment`, `light_condition`); fall back to the paper's methods section only if still unclear. **If neither source tells you, use the sentinel and move on — do not guess.**

- `timepoint` is **never null** — always a non-empty string. Use `"unknown"` as the sentinel.
- `timepoint_hours` may be `null` when the paper doesn't report a sampling time.
- `growth_phase` uses `unknown` (not null) as the sentinel.

| Field | Description | Examples |
|-------|-------------|----------|
| `timepoint` | Short label suitable as a figure-axis tag. Use the sampling time if a time course, else a concise physiological context tag. | `"24h"`, `"120 min"`, `"log phase"`, `"pre-infection"`, `"unknown"` |
| `timepoint_hours` | Numeric conversion to hours (or `null` if the paper doesn't report a time). | `20`, `0.5`, `48`, `null` |
| `growth_phase` | Physiological state at sampling — value from the enum below. | `exponential`, `nutrient_limited`, `infected`, `unknown` |

**`growth_phase` enum** (use an exact value or the `other:<slug>` escape):

| Value | When to use |
|-------|-------------|
| `exponential` | mid-log, cells dividing normally |
| `stationary` | post-log, no more division, no overt stress named |
| `nutrient_limited` | cells clearly starved of N, P, Fe, C, etc. |
| `acclimated_steady_state` | chronic / ≥5-generation exposure to a non-lethal condition |
| `infected` | post-infection timepoint in a phage/viral study |
| `recovery` | post-rescue / post-readdition |
| `diel` | sampled across a light:dark cycle |
| `darkness` | prolonged dark exposure (hours+, not the dark half of diel) |
| `death` | sampled past viable-cell peak |
| `acute_stress` | short (≤6h) stress at still-dividing cells, only when no more-specific phase fits |
| `unknown` | paper genuinely gives no info |
| `other:<slug>` | paper describes a phase the enum doesn't cover (e.g. `other:heat_acclimated`) |

**Default bias:** `unknown` > guessing. Prefer `other:<slug>` over `unknown` when the paper gives positional information the enum misses.

**Minimal-diff rule:** when backfilling existing paperconfigs, only add these three fields to each analysis entry (add `timepoint_hours` only if missing). Do not reflow, re-quote, or reorder other fields — the backfill is reviewed via `git diff`.

#### Optional Fields

| Field | Description | When to Use |
|-------|-------------|-------------|
| `adjusted_p_value_col` | Column name for adjusted p-value | When the CSV has a p-value column (e.g., `"padj"`, `"FDR"`, `"q Value"`) |
| `timepoint` | Human-readable time point label | For display (e.g., `"20h"`, `"48h"`, `"-12h"`) |
| `skip_rows` | Number of header rows to skip when reading CSV | When the CSV has extra header rows before the data columns (e.g., `3`) |
| `pvalue_asterisk_in_logfc` | Boolean, `true` if significance is marked by `*` appended to fold-change values | When there is no separate p-value column and the fold-change column has asterisks (e.g., `"2.5*"`) |
| `prefiltered` | Boolean, `true` if the table only contains significant results | When authors pre-filtered the table to significant genes only |
| `pvalue_threshold` | Adjusted p-value threshold used/stated in the publication (e.g., `0.05`) | When the paper states a specific significance cutoff |
| `logfc_threshold` | Absolute log2 fold-change threshold used/stated in the publication (e.g., `0.8`, `1.0`) | When the paper states a specific fold-change cutoff |

**Note on tables where all values are significant:** Some supplementary tables only list genes that are already significantly differentially expressed (pre-filtered by the authors). In these cases, set `prefiltered: true` on the analysis. This tells the adapter that all rows are significant. You may also omit `adjusted_p_value_col` and `pvalue_asterisk_in_logfc` if not applicable.

**Note on significance thresholds:** If the paper or table legend states specific significance criteria (e.g., "adjusted p-value < 0.05 and |log2FC| >= 1"), record these as `pvalue_threshold` and `logfc_threshold`. The adapter uses these to mark edges as significant or not. If the paper does not state thresholds, omit these fields -- the adapter can apply default thresholds from its constructor settings.

## Edge Label Routing

The OMICSAdapter creates `Changes_expression_of` edges in the knowledge graph. The edge source node is determined by the experiment's `treatment_type` and coculture fields:

| Experiment fields | Edge label | Source node type |
|---|---|---|
| `treatment_organism` + `treatment_taxid` (coculture/viral) | `Changes_expression_of` | `Experiment` (linked to `OrganismTaxon` via `Tests_coculture_with`) |
| Environmental treatment types (no treatment_organism) | `Changes_expression_of` | `Experiment` |

All expression edges now use the unified `Changes_expression_of` label, with the `Experiment` node as the source. The experiment's `treatment_type` property enables filtering by category (e.g., find all coculture experiments, all nitrogen stress experiments).

For coculture experiments, the `Experiment` node is linked to the treatment organism via a `Tests_coculture_with` edge. This replaces the old split between `Condition_changes_expression_of` and `Coculture_changes_expression_of`.

The adapter also emits a `Has_experiment` edge from the `Publication` node to the `Experiment` node for each experiment, enabling navigation from publications to the experiments they describe.

## Organism Reference Data

Organism names, taxids, and assembly accessions must be consistent with these two canonical CSV files:

- `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` -- strains with loaded genomes
- `data/Prochlorococcus/treatment_organisms.csv` -- genus-level or non-specific organisms (no genome loaded)

The `organism` field format is `"<Genus> <strain_name>"` where `strain_name` matches the `strain_name` column in `cyanobacteria_genomes.csv`, or just the genus name for treatment-only organisms.

### Genome-loaded strains (from `cyanobacteria_genomes.csv`)

Always include `treatment_assembly_accession` on the experiment when using one of these as `treatment_organism`:

| organism / treatment_organism value | ncbi_taxon_id | ncbi_accession (assembly) |
|---|---|---|
| Prochlorococcus MED4 | 59919 | GCF_000011465.1 |
| Prochlorococcus AS9601 | 146891 | GCF_000015645.1 |
| Prochlorococcus MIT9301 | 167546 | GCF_000015965.1 |
| Prochlorococcus MIT9312 | 74546 | GCF_000012645.1 |
| Prochlorococcus MIT9313 | 74547 | GCF_000011485.1 |
| Prochlorococcus NATL1A | 167555 | GCF_000015685.1 |
| Prochlorococcus NATL2A | 59920 | GCF_000012465.1 |
| Prochlorococcus RSP50 | 1924285 | GCF_001989415.1 |
| Synechococcus CC9311 | 64471 | GCF_000014585.1 |
| Synechococcus WH8102 | 84588 | GCF_000195975.1 |
| Alteromonas macleodii MIT1002 | 28108 | GCF_901457835.2 |
| Alteromonas macleodii EZ55 | 28108 | GCF_901457815.2 |
| Alteromonas macleodii HOT1A3 | 28108 | GCF_001578515.1 |

### Non-genome organisms (from `treatment_organisms.csv`)

Omit `treatment_assembly_accession` for these genus-level organisms:

| treatment_organism value | treatment_taxid | notes |
|---|---|---|
| Phage | 10239 | no assembly |
| Marinobacter | 413470 | genus-level |
| Thalassospira | 191411 | genus-level |
| Pseudohoeflea | 398581 | genus-level |
| Alteromonas | 28108 | genus-level (when no specific strain) |
| Synthetic heterotroph community | 413470 | placeholder taxid (Marinobacter) |

If an organism name and its taxid/assembly disagree, validate by checking the paper's supplementary legend files first, then use the central CSVs as authoritative for taxid and assembly accession.

## CSV Data File Requirements

- Must be a valid CSV file readable by pandas
- Must contain at minimum:
  - A column with gene/protein identifiers (referenced by `name_col`)
  - A column with log2 fold-change values (referenced by `logfc_col`)
- Optionally: a column with adjusted p-values (referenced by `adjusted_p_value_col`)
- Gene identifiers should ideally be locus tags or gene IDs that can be matched to existing gene nodes in the knowledge graph
- Rows with empty/NA gene identifiers or non-numeric fold-change values are automatically skipped

## Complete Examples

### Example 1: Coculture Experiment (Single Timepoint)

```yaml
publication:
  papername: "Aharonovich 2016"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paper.pdf"

  experiments:
    coculture_hot1a3:
      name: "DE of Prochlorococcus MED4 in Coculture with Alteromonas HOT1A3 vs Axenic"
      organism: "Prochlorococcus MED4"
      omics_type: RNASEQ
      test_type: Rockhopper
      treatment_type: coculture
      control_condition: Axenic
      treatment_condition: "Coculture with Alteromonas HOT1A3"
      experimental_context: "in Pro99 medium under continuous light"
      treatment_organism: "Alteromonas macleodii HOT1A3"
      treatment_taxid: 28108

  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/de_genes_med4_t20h.csv"
      statistical_analyses:
        - id: DE_coculture_vs_axenic_med4_20h
          experiment: coculture_hot1a3
          timepoint_hours: 20
          name_col: Synonym
          logfc_col: "log 2 fold change coculture/axenic"
          adjusted_p_value_col: "q Value"
```

### Example 2: Environmental Stress (Salt Stress)

```yaml
publication:
  papername: "Al-Hosani 2015"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Al-Hosani 2015/paper.pdf"

  experiments:
    salt_acclimation:
      name: "DE of Prochlorococcus AS9601 under salt acclimation vs normal seawater"
      organism: "Prochlorococcus AS9601"
      omics_type: RNASEQ
      test_type: DESeq
      treatment_type: salt
      control_condition: "Normal seawater (3.8% salt)"
      treatment_condition: "Salt-acclimated (5% salt)"
      experimental_context: "Axenic cells in PRO99 medium at 22C under continuous light"
      medium: "PRO99 natural seawater medium"
      temperature: "22C"
      light_condition: continuous light
      light_intensity: "46 umol photons m-2 s-1"

  supplementary_materials:
    supp_table_3:
      type: csv
      filename: "data/Prochlorococcus/papers_and_supp/Al-Hosani 2015/de_genes_salt.csv"
      statistical_analyses:
        - id: DE_salt_acclimated_vs_seawater_AS9601
          experiment: salt_acclimation
          timepoint_hours: null
          name_col: "Gene id"
          logfc_col: log2FoldChange
          adjusted_p_value_col: padj
```

### Example 3: Time-Course Experiment (Multiple Analyses Referencing Same Experiment)

When a single CSV contains columns for multiple timepoints, all analyses reference the same experiment:

```yaml
publication:
  papername: "Biller 2016"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Biller 2016/paper.pdf"

  experiments:
    coculture_mit1002_natl2a:
      name: "DE Prochlorococcus NATL2A coculture vs axenic"
      organism: "Prochlorococcus NATL2A"
      omics_type: RNASEQ
      test_type: DESeq2
      treatment_type: coculture
      control_condition: Axenic
      treatment_condition: "Co-culture with Alteromonas macleodii MIT1002"
      treatment_organism: "Alteromonas macleodii MIT1002"
      treatment_taxid: 28108
      treatment_assembly_accession: GCF_901457835.2

  supplementary_materials:
    supp_table_2:
      type: csv
      filename: "data/.../PRO_DE_genes.csv"
      statistical_analyses:
        - id: DE_coculture_vs_axenic_NATL2A_2h
          experiment: coculture_mit1002_natl2a
          timepoint_hours: 2
          name_col: Original_NCBI_ID
          logfc_col: "2 hours"
          skip_rows: 3
          pvalue_asterisk_in_logfc: true
        - id: DE_coculture_vs_axenic_NATL2A_12h
          experiment: coculture_mit1002_natl2a
          timepoint_hours: 12
          name_col: Original_NCBI_ID
          logfc_col: "12 hours"
          skip_rows: 3
          pvalue_asterisk_in_logfc: true
        - id: DE_coculture_vs_axenic_NATL2A_24h
          experiment: coculture_mit1002_natl2a
          timepoint_hours: 24
          name_col: Original_NCBI_ID
          logfc_col: "24 hours"
          skip_rows: 3
          pvalue_asterisk_in_logfc: true
```

### Example 4: Multiple Experiments from Same Paper

When a paper has separate experiments for different organisms or conditions:

```yaml
publication:
  papername: "Biller 2018"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Biller 2018/paper.pdf"

  experiments:
    extended_darkness_axenic:
      name: "DE Prochlorococcus NATL2A axenic in extended darkness vs diel"
      organism: "Prochlorococcus NATL2A"
      omics_type: RNASEQ
      test_type: DESeq2
      treatment_type: darkness
      control_condition: "13:11 diel light:dark cycle"
      treatment_condition: "Extended darkness"
      experimental_context: "Axenic Prochlorococcus NATL2A in Pro99 at 24C"
      medium: "Pro99 natural seawater medium"
      temperature: "24C"
      light_condition: constant darkness
      light_intensity: "0 umol photons m-2 s-1"

    extended_darkness_coculture:
      name: "DE Prochlorococcus NATL2A coculture in extended darkness vs diel"
      organism: "Prochlorococcus NATL2A"
      omics_type: RNASEQ
      test_type: DESeq2
      treatment_type: darkness
      control_condition: "13:11 diel light:dark cycle"
      treatment_condition: "Extended darkness"
      experimental_context: "Prochlorococcus NATL2A co-cultured with Alteromonas macleodii MIT1002 in Pro99 at 24C"
      medium: "Pro99 natural seawater medium"
      temperature: "24C"
      light_condition: constant darkness
      light_intensity: "0 umol photons m-2 s-1"

  supplementary_materials:
    supp_table_s3:
      type: csv
      filename: "data/.../de_genes.csv"
      statistical_analyses:
        - id: DE_extended_darkness_vs_diel_axenic_NATL2A_1h
          experiment: extended_darkness_axenic
          timepoint_hours: 1
          name_col: "NCBI ID"
          logfc_col: "Axenic, 36 hours"
          pvalue_asterisk_in_logfc: true
        - id: DE_extended_darkness_vs_diel_coculture_NATL2A_1h
          experiment: extended_darkness_coculture
          timepoint_hours: 1
          name_col: "NCBI ID"
          logfc_col: "Co-culture, 36 hours"
          pvalue_asterisk_in_logfc: true
        # ... repeat for additional timepoints
```

### Example 5: Multiple Tables from Same Paper

When a paper has separate CSV files for different organisms:

```yaml
  supplementary_materials:
    supp_table_pro:
      type: csv
      filename: "data/.../pro_de_genes.csv"
      statistical_analyses:
        - id: DE_pro_experiment_1
          experiment: pro_experiment
          # ...
    supp_table_alt:
      type: csv
      filename: "data/.../alt_de_genes.csv"
      statistical_analyses:
        - id: DE_alt_experiment_1
          experiment: alt_experiment
          # ...
```

## Strain-Level Resource Paperconfigs

Some ID translation files are not associated with any specific paper — they apply to all papers for a given strain (e.g., a strain-wide locus tag cross-reference table). These belong in a **strain resource paperconfig**: a `paperconfig.yaml` with no `publication` block, containing only `id_translation` and/or `annotation_gff` entries.

**Placement:** `data/Prochlorococcus/papers_and_supp/<Strain>_resources/paperconfig.yaml`
**Registration:** listed in `paperconfig_files.txt` like any other paperconfig.

```yaml
# Strain-level shared gene ID resources for Prochlorococcus MIT9313
# No publication block — this is not a paper
supplementary_materials:
  mit9313_id_translation:
    type: id_translation
    filename: "data/Prochlorococcus/papers_and_supp/MIT9313_resources/MIT9313_genbank.tsv"
    sep: "\t"
    organism: "Prochlorococcus MIT9313"
    id_columns:
      - column: "NCBI locus tag"
        id_type: locus_tag_ncbi
      - column: "Old locus tag"
        id_type: old_locus_tag
      - column: "Gene name"
        id_type: gene_name
```

`omics_adapter` ignores strain resource paperconfigs entirely (no `publication` block → no publication node emitted; no `statistical_analyses` → no edges emitted).

## Registration

After creating the paperconfig, register it in the pipeline:

1. Add the path to `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` (one path per line)
2. Or instantiate it directly in `create_knowledge_graph.py`:

```python
omics = OMICSAdapter(config_file="data/Prochlorococcus/papers_and_supp/<AuthorName Year>/paperconfig.yaml")
omics.download_data(cache=True)
bc.write_nodes(omics.get_nodes())
bc.write_edges(omics.get_edges())
```

## Validation Checklist

Before submitting a paperconfig:

- [ ] `papermainpdf` path points to an existing PDF file (omit for strain-level resource paperconfigs with no `publication` block)
- [ ] All `filename` paths point to existing CSV/GFF files
- [ ] `experiments` block exists with at least one experiment entry (for publication configs)
- [ ] Each experiment has all required fields: `name`, `organism`, `omics_type`, `test_type`, `treatment_type`, `treatment_condition`, `control_condition`
- [ ] `organism` values are canonical (match `cyanobacteria_genomes.csv` or `treatment_organisms.csv`)
- [ ] `treatment_type` is a canonical value (see Treatment Types table)
- [ ] `omics_type` is one of: RNASEQ, MICROARRAY, PROTEOMICS, METABOLOMICS
- [ ] `test_type` is a canonical value (DESeq2, DESeq, edgeR, Rockhopper, microarray, etc.)
- [ ] Each analysis has `id`, `experiment`, `name_col`, `logfc_col`
- [ ] Each analysis `experiment` reference points to a valid key in `experiments` block
- [ ] Each analysis has `timepoint_hours` (number or `null`)
- [ ] `name_col` matches an actual column name in the CSV
- [ ] `logfc_col` matches an actual column name in the CSV
- [ ] `adjusted_p_value_col` (if specified) matches an actual column name in the CSV
- [ ] Each analysis `id` is unique within the file
- [ ] `skip_rows` is set correctly if the CSV has extra header rows
- [ ] Gene/protein identifiers in `name_col` use locus tags or IDs that can match existing gene nodes
- [ ] `prefiltered` is set to `true` if the table only contains significant results
- [ ] `pvalue_threshold` and `logfc_threshold` are set if the paper states specific significance criteria
- [ ] For coculture/viral experiments: `treatment_organism` and `treatment_taxid` are set on the experiment
- [ ] `treatment_taxid` matches `ncbi_taxon_id` in `cyanobacteria_genomes.csv` or `treatment_organisms.csv`
- [ ] `treatment_assembly_accession` (if present) matches `ncbi_accession` in `cyanobacteria_genomes.csv` -- omit for genus-level organisms
- [ ] Any new organism not already in the central CSVs has been added to the appropriate file (see "Registering New Organisms" below)
- [ ] For `id_translation` entries: `organism` and `id_columns` are declared; all column names match actual CSV headers
- [ ] For `annotation_gff` entries: `organism` is declared and `filename` exists
- [ ] `id_translation` entries are listed **before** `csv` entries for the same organism so alt-IDs are in the lookup when DE data is processed
- [ ] `id_columns` column names match actual headers in the CSV (check with `head -1` or pandas read)

## Workflow

When the user invokes this skill with a paper directory name (e.g., `/paperconfig "Author Year"`):

1. Look for the paper directory under `data/Prochlorococcus/papers_and_supp/$ARGUMENTS/`
2. List all files in that directory (PDF, CSV, XLSX, TXT, DOCX, GFF, GTF, etc.)
3. Read any legend/description text files first -- these explain what the supplementary tables contain (column meanings, significance criteria, experimental details)
4. Read the PDF to understand the experiment (organisms, conditions, methods, statistical approach)
5. For each data file (CSV, TSV, XLSX):
   - Read the headers and a few sample rows
   - Check: are there extra header rows before the data? Set `skip_rows` if so
   - Check: are fold-change values appended with `*`? Set `pvalue_asterisk_in_logfc: true`
   - Determine the delimiter (`sep`: `","` or `"\t"`)
   - **Classify every column:**
     - Which columns contain **gene identifiers**? → declare in `id_columns` with `id_type`:
       - Locus tags (PMM0001, P9301_RS09095 style) → `locus_tag` or `locus_tag_ncbi`
       - Old/legacy locus tags (P9301_05911 style) → `old_locus_tag`
       - UniProt accessions (Q7V6L1) → `uniprot_accession`
       - UniProt entry names (DNAA_PROM0) → `uniprot_entry_name`
       - RefSeq WP_ accessions → `protein_id_refseq`
       - JGI IMG integer IDs → `jgi_id`
       - Gene symbols (dnaA, dnaN) → `gene_name`
       - Microarray probesets → `probeset`
       - RAST fig|...|peg.N IDs → `rast_id`
       - Paper-specific IDs without a standard type → `annotation_specific`
     - Which columns contain **functional descriptions** (product names, annotations)? → declare in `product_columns`
     - Which columns are the **fold-change** and **p-value** columns? → `logfc_col`, `adjusted_p_value_col`
   - Determine whether this file is a DE results table (→ `type: csv`) or a pure ID mapping table with no expression data (→ `type: id_translation`)
6. Check for GFF/GTF files in the directory. If present:
   - Add an `annotation_gff` entry with `type: annotation_gff`, `filename`, and `organism`
   - These are processed by `build_gene_id_mapping.py` to bridge paper-specific IDs to canonical locus tags; the omics adapter ignores them
7. Draft the `paperconfig.yaml` following the schema above:
   - Define the `experiments` block first with all experiment-level metadata (organism, omics_type, test_type, treatment_type, conditions, coculture fields)
   - Each analysis references an experiment key and adds only per-timepoint fields (timepoint_hours, name_col, logfc_col, etc.)
   - List `id_translation` and `annotation_gff` entries **before** any `csv` entries for the same organism so alt-IDs are in the lookup when DE data is processed
   - If the `name_col` IDs are non-standard (JGI IDs, probesets, RAST IDs), ensure an `id_translation` entry maps them to locus tags before the DE CSV is processed
   - Declare `id_columns` on every `csv` entry to document the gene ID provenance, even when `name_col` is already `locus_tag` (post-fix-gene-ids state)
8. **Register any new organisms** in the central CSVs (see "Registering New Organisms" below)
9. Run the validation script (see below) to check all paths, columns, references, and ID uniqueness
10. Register the config in `paperconfig_files.txt`

### Registering New Organisms

After drafting the paperconfig, read both central CSV files and check every `organism` and `treatment_organism` value used:

```bash
# Check current contents
cat data/Prochlorococcus/genomes/cyanobacteria_genomes.csv
cat data/Prochlorococcus/treatment_organisms.csv
```

For each organism name that is **not already present**:

**If the organism has a specific strain with a known genome assembly** → add to `cyanobacteria_genomes.csv`:
- `ncbi_accession`: RefSeq assembly accession (GCF_...) from NCBI — look this up from the paper's methods or NCBI
- `cyanorak_organism`: leave empty if not in CyanoRAK
- `ncbi_taxon_id`: NCBI Taxonomy ID from the paper or NCBI
- `strain_name`: strain identifier as used in the `organism` field (e.g., `MIT9313`, `HOT1A3`)
- `data_dir`: `cache/data/<genus>/genomes/<strain_name>/` — fill in the expected path even if data hasn't been downloaded yet

**If genus-level or no specific genome** (e.g., just "Marinobacter", "Phage") → add to `treatment_organisms.csv`:
- `ncbi_taxon_id`: NCBI Taxonomy ID
- `organism_name`: the name as used in the `treatment_organism` field

Add a comment line above the new row explaining the context (e.g., which paper first introduced it).

### Validation Script

After creating the paperconfig, always run the validation script:

```bash
uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/<Author Year>/paperconfig.yaml"
```

This script checks:
- YAML is parseable and has required top-level structure
- PDF file exists
- All CSV files exist and are readable
- `experiments` block exists and is well-formed
- Each experiment has required fields: `name`, `organism`, `omics_type`, `test_type`, `treatment_type`, `treatment_condition`, `control_condition`
- Canonical vocabulary: `organism`, `omics_type`, `test_type`, `treatment_type`, `treatment_organism`
- All CSV files exist and are readable
- Column names (`name_col`, `logfc_col`, `adjusted_p_value_col`) match actual CSV headers
- `skip_rows` is applied correctly when reading CSV headers
- Each analysis has required fields: `id`, `experiment`, `name_col`, `logfc_col`
- Each analysis `experiment` reference points to a valid key in `experiments` block
- Each analysis has `timepoint_hours` (number or null; warns if missing)
- All analysis IDs are unique
- `logfc_col` values look numeric (warns on non-numeric values)
- `id_columns` and `product_columns` column names match CSV headers
- `id_translation` entries have `organism` and `id_columns`
- `annotation_gff` entries have `organism` and valid file extension

The script exits with code 0 on success, 1 on validation failure.

### Gene ID Mapping

If the paper CSV uses standard gene names or locus-tag variants (e.g., `secE`, `rps13,rpsM`, old-format `PMM0814`), use the `/fix-gene-ids` skill:

1. Run `/check-gene-ids` to confirm the fix strategy is `CREATE_MAPPING_CSV`
2. Run `/fix-gene-ids` on the paperconfig to create `_with_locus_tag.csv` files
3. Update the paperconfig: change `filename` to the new CSV and `name_col` to `"locus_tag"`; add `original_filename` pointing to the original CSV

If the paper CSV uses **non-standard IDs** (JGI catalog IDs, probesets, paper-specific annotation IDs):

1. Declare the primary ID column in `id_columns` on the `csv` entry (e.g., `id_type: jgi_id`)
2. Add an `id_translation` entry (or `annotation_gff`) to the same paperconfig to bridge those IDs to locus tags
3. Run `uv run python multiomics_kg/download/build_gene_id_mapping.py` to rebuild the per-strain mapping
4. Then run `/fix-gene-ids` — it will now find matches via the enriched mapping
