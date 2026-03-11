# Methods: Differential Expression Integration via the Omics Adapter

## Overview

Integrating differential expression results from diverse publications into a unified knowledge graph requires a declarative system that can accommodate variation in experimental designs, statistical reporting conventions, gene identifier vocabularies, and edge semantics. We developed a two-component system: (1) a YAML-based configuration schema (`paperconfig.yaml`) that declaratively describes each publication's supplementary data tables and their experimental metadata, and (2) an omics adapter that reads these configurations, generates graph nodes for publications and experimental conditions, and emits `Condition_changes_expression_of` or `Coculture_changes_expression_of` edges linking causal factors to genes with quantitative expression properties. The system currently integrates 24 publications spanning RNA-seq, proteomics, metabolomics, and microarray experiments across 13 bacterial strains.

## Paperconfig YAML Schema

Each publication is described by a single `paperconfig.yaml` file stored alongside its supplementary data files. The schema has three levels: publication metadata, environmental condition definitions, and supplementary material declarations.

### Publication metadata

The top-level `publication` block declares the paper's identity and the path to its main PDF file:

```yaml
publication:
  papername: "Author Year"
  doi: "10.xxxx/..."
  papermainpdf: "data/.../paper.pdf"
```

The `papername` serves as a human-readable identifier; the DOI, when available, becomes the canonical node ID in the graph.

### Environmental condition nodes

Publications studying environmental stresses (light, nutrients, gas composition) define reusable environmental condition blocks:

```yaml
  environmental_conditions:
    condition_id:
      name: "Nitrogen starvation"
      condition_type: "nutrient_stress"
      nitrogen_level: "starved"
      light_condition: "continuous white light"
      description: "N-starved Pro99 at 22C"
```

Each condition becomes an `EnvironmentalCondition` node in the graph, with a unique ID formed by concatenating the publication DOI and the condition's local key. Statistical analyses reference these conditions by their local key, creating a many-to-one relationship: multiple time points or replicates from the same experimental condition all share a single source node.

### Supplementary material declarations

The `supplementary_materials` block maps each data file to one or more statistical analyses:

```yaml
  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "data/.../results.csv"
      organism: "Prochlorococcus MED4"
      sep: ","
      skip_rows: 0
      statistical_analyses:
        - id: "unique_analysis_id"
          type: RNASEQ
          name_col: "Gene"
          logfc_col: "log2FoldChange"
          adjusted_p_value_col: "padj"
          control_condition: "Axenic"
          treatment_condition: "Coculture"
          experimental_context: "Pro99 medium, continuous light"
```

Each analysis within a table must carry a unique `id` field that, combined with the publication DOI, forms a globally unique edge ID prefix. This design supports parallel edges: the same gene may have multiple expression edges from different time points, conditions, or publications, each with distinct properties.

File-level settings (`organism`, `sep`, `skip_rows`) propagate to individual analyses unless overridden, reducing repetition in multi-analysis tables.

## Edge Source Determination

The causal factor in a differential expression experiment is modeled as the source node of each expression edge. Two source types are supported, selected per analysis by the fields present in the paperconfig; they produce different edge labels to enable unambiguous queries by experiment type:

**Environmental condition source.** When the analysis declares an `environmental_treatment_condition_id`, the adapter uses the corresponding `EnvironmentalCondition` node as the edge source and emits a `Condition_changes_expression_of` edge. This applies to stress experiments (light shifts, nutrient limitation, gas perturbation) where the causal factor is an abiotic condition rather than a biological entity. Approximately 95,000 edges in the current graph use this source type.

**Organism taxon source.** When the analysis instead declares `treatment_organism` and `treatment_assembly_accession` (or `treatment_taxid`), the adapter uses the corresponding `OrganismTaxon` node as the edge source and emits a `Coculture_changes_expression_of` edge. This applies to coculture experiments where gene expression changes in one organism are attributed to the presence of another organism. The assembly accession is preferred for source identification (yielding an INSDC GCF-based CURIE), with NCBI taxonomy ID as a fallback. Approximately 12,000 edges in the current graph use this source type.

LLM agents and Cypher queries that need to span both experiment types can use `r:Condition_changes_expression_of|Coculture_changes_expression_of` in relationship patterns.

## Publication Metadata Extraction

Publication nodes carry structured metadata (title, authors, journal, DOI, abstract, study type, organisms) extracted from the paper's main PDF using LangChain-based LLM processing. The `PDFPublicationExtractor` reads the first five pages of each PDF (stopping at the References section), cleans OCR artifacts and formatting noise, and prompts an LLM to return structured JSON with publication fields. Results are cached in a persistent JSON file to avoid redundant LLM API calls across builds. When the PDF is unavailable or extraction fails, the adapter falls back to the `doi` and `papername` fields from the paperconfig for node identification.

## Edge Properties

Each `Condition_changes_expression_of` and `Coculture_changes_expression_of` edge carries the following properties:

| Property | Type | Source | Required |
|----------|------|--------|----------|
| `log2_fold_change` | float | Data CSV column (`logfc_col`) | Yes |
| `adjusted_p_value` | float | Data CSV column (`adjusted_p_value_col`) or inferred from asterisk convention | No |
| `expression_direction` | string | Derived: `"up"` if log2FC > 0, `"down"` otherwise | Yes |
| `significant` | boolean | Computed from thresholds or asterisk convention | When determinable |
| `control_condition` | string | Paperconfig analysis field | When available |
| `experimental_context` | string | Paperconfig analysis field | When available |
| `time_point` | string | Paperconfig analysis field | When available |
| `publications` | string[] | Publication DOI | Yes |

The `log2_fold_change` value is the only strictly required numeric property; rows with missing, non-numeric, or non-finite (NaN, Inf) fold-change values are skipped. The `adjusted_p_value` may be null when the original study did not report per-gene significance values (e.g., pre-filtered gene lists).

### Expression direction

The `expression_direction` property is derived deterministically from the sign of `log2_fold_change`: positive values indicate upregulation (`"up"`), negative values indicate downregulation (`"down"`). This derived field simplifies downstream queries by allowing agents to filter by direction without interpreting numeric signs.

## Significance Determination

Statistical significance is determined through a three-priority cascade, evaluated per row:

**Priority 1 — Prefiltered tables.** When the analysis declares `prefiltered: true`, all rows are considered significant. This applies to publications that report only genes passing their own significance thresholds.

**Priority 2 — Asterisk convention.** Some publications encode significance as asterisks appended to fold-change values (e.g., `"1.23 *"` for significant, `"0.45"` for non-significant). When `pvalue_asterisk_in_logfc: true` is set, the adapter strips asterisks from the fold-change string and records their presence as a boolean significance flag. A synthetic `adjusted_p_value` is assigned: the configured `pvalue_threshold` (default 0.05) when significant, or 1.0 when not. This convention is used by several time-course studies in the dataset (e.g., Biller 2016, Coe 2024) where per-gene p-values were computed but only the binary significance outcome was reported in supplementary tables.

**Priority 3 — Threshold-based classification.** When neither prefiltering nor asterisk conventions apply, significance is computed by comparing the row's `adjusted_p_value` against a threshold (default 0.05, overridable per analysis via `pvalue_threshold`) and the absolute `log2_fold_change` against a threshold (default 1.0, overridable via `logfc_threshold`). A row is significant only when both criteria are met. When neither p-value nor fold-change threshold can be evaluated (both values absent), significance is recorded as null.

The adapter supports two operating modes: in `"all"` mode (default), all rows produce edges regardless of significance, with a `significant` boolean property for downstream filtering; in `"significant_only"` mode, rows classified as non-significant are excluded from edge generation.

## Pre-resolved CSV Integration

When gene identifiers in a publication's CSV file are not canonical locus tags, the adapter integrates with the gene ID resolution pipeline (described in the gene ID mapping methods). The `resolve_paper_ids` module pre-computes a `_resolved.csv` file alongside each source CSV, containing two additional columns: `locus_tag` (the resolved canonical identifier, or null when unresolvable) and `resolution_method` (the tier and column that produced the match).

At edge generation time, the adapter probes for the existence of a `_resolved.csv` file corresponding to each source CSV. When present and containing a `locus_tag` column, the adapter reads the pre-resolved file and uses the `locus_tag` column as the gene identifier, ignoring the original `name_col`. Rows where `locus_tag` is null (unresolved genes) are skipped, producing no edge rather than a dangling reference to a non-existent gene node. The resolved file is always a clean comma-separated format, so the adapter bypasses the original file's `sep` and `skip_rows` settings.

When no pre-resolved file exists (e.g., for papers whose `name_col` is already `locus_tag`), the adapter reads the original CSV directly and uses the `name_col` values as gene identifiers.

## String Sanitization

All string property values are sanitized before yielding nodes and edges. Single quotes are replaced with carets (`'` → `^`) and pipe characters are removed (`|` → `,`), because the BioCypher CSV export wraps string fields in single quotes and uses pipe as the array delimiter. Unsanitized values containing these characters would corrupt the `neo4j-admin import` input.

## Multi-Paper Aggregation

The `MultiOMICSAdapter` wrapper reads a text file (`paperconfig_files.txt`) listing all paperconfig paths, instantiates one `OMICSAdapter` per publication, and aggregates their nodes and edges. This design allows each adapter instance to operate independently — loading its own config, extracting its own publication metadata, and generating its own edges — while the wrapper provides the unified interface expected by the knowledge graph build pipeline. The current deployment processes 24 paperconfig files, generating approximately ~188K expression edges (`Condition_changes_expression_of` + `Coculture_changes_expression_of`) across all publications.

## Edge Identity and Parallel Edges

Edge IDs follow the pattern `{publication_doi}_{analysis_id}_{gene_curie}`, ensuring global uniqueness. The `analysis_id` component (required and validated for uniqueness within each publication) distinguishes edges from different experimental conditions targeting the same gene. This allows the graph to represent, for example, a gene's response to nitrogen stress at six different time points as six distinct edges, each with its own `log2_fold_change`, `time_point`, and `significant` values.

The adapter validates analysis ID uniqueness at load time: duplicate IDs within a single publication raise an immediate error, preventing silent edge collisions during graph construction.

## Data Validation

The adapter performs multi-level validation during edge generation:

1. **Column validation**: Required columns (`name_col`, `logfc_col`) must exist in the CSV; missing columns cause the entire table to be skipped with a warning. Optional columns (`adjusted_p_value_col`) produce a warning when absent but do not prevent edge generation.

2. **Row-level filtering**: Rows are skipped when the gene identifier is empty or null, the fold-change value is missing, non-numeric, or non-finite (NaN, Inf), or (in `significant_only` mode) the row fails the significance test. Skipped row counts are logged per table for diagnostics.

3. **Gene identifier cleaning**: Whitespace and trailing asterisks are stripped from gene IDs before CURIE construction. Empty strings after cleaning are treated as missing.

4. **Type coercion**: Fold-change and p-value strings are parsed to float with explicit error handling; unparseable values cause the row to be skipped rather than propagating errors.
