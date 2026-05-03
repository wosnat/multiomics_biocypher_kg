# Metabolomics paper integration — Phase 2 of metabolite scaffold

**Date:** 2026-05-03
**Status:** Design approved; implementation plan TBD
**Replaces:** Phase 2 sketch in `2026-04-28-metabolite-scaffold-design.md` (lines 29–46)
**Validation papers:** Capovilla 2023, Kujawinski 2023

## Context

Phase 1 of the metabolite scaffold (Foundation, Reactions, Transport+CAZy) has landed: `Metabolite` and `Reaction` nodes exist; the MNX resolver and `kegg_data.json` cache drive node emission; `Organism_has_metabolite` is materialized post-import from the catalysis + transport gene paths.

This spec activates Phase 2: integration of paper-reported metabolite measurements (concentrations + presence flags) into the KG via a new `MetaboliteAssay` node type plus measurement edges anchored on `Metabolite`. It fully replaces the brief Phase-2 sketch with a concrete schema, paperconfig entry shape, resolution flow, post-import Cypher, and commit-gated test plan.

Phase 2 also extends Phase 1 in two places: `Metabolite.evidence_sources` gains the `"metabolomics"` value, and paper-discovered compounds become a third path that `Metabolite` nodes (and `Organism_has_metabolite` edges) materialize from.

**Out of this spec, deferred:** paired-comparison fold-change/p-value edges (no current paper has them); `MassFeature` node type for unidentified peaks; isotopologue tracing (Szul 2019); `Assay_classifies_metabolite` (no current paper needs categorical metabolite values).

## 1 — Architecture

Six moving parts; each independently testable.

### 1.1 — Step 6 extension (`build_kegg_metabolism_xrefs.py`)

Walks all paperconfigs for `metabolite_assays_table` entries; resolves each metabolite name via the MNX resolver. Three effects:

- **`kegg_data.json` extended:** paper-measured compounds added to the extended set. KEGG-primary compounds enter the existing compounds bucket; non-KEGG (`chebi:*`, `mnx:*`) compounds enter `additional_compounds`. `evidence_sources` is unioned per compound (e.g. `["metabolism", "metabolomics"]` when also gene-reachable; `["metabolomics"]` for measurement-only).
- **`cache/data/metabolomics/metabolite_id_mapping.json`** (new lightweight cache): name-alias → primary-id lookup harvested across all papers, plus per-paper resolution traces. Step 7 reads only this; the 4 GB MNX resolver is loaded once per build, in step 6 only.
- **Pathway extension:** the existing extended pathway set (gene-catalysis-reachable ∪ transport-reachable) is unioned with paper-measurement-reachable. Existing `Metabolite_in_pathway` materialization picks up the new edges automatically. Only KEGG-primary measured compounds gain pathway memberships in v1; non-KEGG measured compounds (`chebi:*`/`mnx:*`) get no pathway edges.

Resolution priority per row, for each `metabolite_assays_table` entry:
1. `id_col` cell (when present and non-empty), parsed per `id_type` → primary ID
2. `aliases_file` override (paper-local YAML) on `name_col`
3. `formula_col` disambiguation among multiple MNX hits (when `formula_col` is set and name resolution returns multiple candidates)
4. MNX resolver against `name_col`
5. Else → `unresolved`, logged in per-paper resolution report

### 1.2 — Step 7 (new — `multiomics_kg/download/resolve_paper_metabolites.py`)

Pure CSV rewriter. For each `metabolite_assays_table` source CSV, looks up each row's `(name_col, id_col?)` in `metabolite_id_mapping.json`, writes `<stem>_resolved.csv` with two added columns:

- `metabolite_id` — e.g. `kegg.compound:C00031`, NaN when unresolved
- `resolution_method` — see enum in §3.3

Also writes `<stem>_resolution_report.json` (per-paper match-rate stats; unresolved + fuzzy + ambiguous-multi-hit names listed).

Run order: `... 5 → 6 → 7`. Step 7 does **not** load the MNX resolver. Mirrors gene-side step 4. Docker builds skip the resolver entirely.

### 1.3 — `metabolism_adapter.py` — verify only, no code change required

The adapter already reads `kegg_data.json` and emits `Metabolite` nodes plus `Metabolite_in_pathway` / `Reaction_has_metabolite` edges. With paper-measured compounds now in the cache (carrying `evidence_sources` that includes `"metabolomics"`), they flow through automatically.

### 1.4 — New adapter `multiomics_kg/adapters/metabolite_assay_adapter.py`

Mirrors `observations_adapter.py` (the DerivedMetric adapter) module structure: single file with `MetaboliteAssayAdapter` (single-paper) + `MultiMetaboliteAssayAdapter` (registry-driven).

Mirroring details:

- New `iter_metabolite_assays_tables` helper added to `paperconfig_utils.py`
- `_make_metabolite_assay_id(doi, entry_key, assay_key)` → node ID `metabolite_assay:{doi_short}:{entry_key}:{assay_key}`
- Same denormalized parent-Experiment block on each `MetaboliteAssay` node: `experiment_id`, `organism_name`, `publication_doi`, `omics_type`, `treatment_type`, `background_factors`, `treatment`, `light_condition`, `experimental_context`, `compartment`
- Same `_clean_str` sanitizer + `_resolve_csv_path` helper
- Replicate aggregation helper `_aggregate_replicates(values: list[float], null_values: set[str]) → (mean, sd, n_replicates, n_non_zero, replicate_values)`
- Cell-format parser for `embedded_mean_sd_n` (handles `"0.00054 (8.8e-05), n=2"`, `"nd"`, `"NA"`); default `numeric` parser is a plain float coercion
- Same registry-walk filter in `MultiMetaboliteAssayAdapter`

Emits:
- `MetaboliteAssay` nodes (one per metabolite_assays_table entry × assay key)
- `Assay_quantifies_metabolite` edges (numeric, replicate-aggregated)
- `Assay_flags_metabolite` edges (boolean presence)
- `PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism` binding edges (BioCypher PascalCase form, mirrors `PublicationHasDerivedMetric`)

**Does not** emit `Publication`, `Experiment`, or `OrganismTaxon`. Those stay owned by `omics_adapter` (with a small fix so its experiments-block iteration is unconditional, covering metabolomics-only papers like Kujawinski 2023). `create_knowledge_graph.py` instantiates `MultiMetaboliteAssayAdapter` after `MultiObservationsAdapter`.

`Assay_classifies_metabolite` deferred for v1.

### 1.5 — Post-import additions

| Target | Property/edge | Behavior |
|---|---|---|
| `Assay_quantifies_metabolite` | `detection_status` (str) | derived: `detected` (n_non_zero == n_replicates) \| `sporadic` (0 < n_non_zero < n_replicates) \| `not_detected` (n_non_zero == 0) |
| `Assay_quantifies_metabolite` | `rank_by_metric`, `metric_percentile`, `metric_bucket` | parallel to DM numeric edges; only when parent assay `rankable="true"` |
| `MetaboliteAssay` | `total_metabolite_count` | distinct count of `Assay_*_metabolite` outgoing edges |
| `MetaboliteAssay` | `value_min`, `value_max`, `value_q1`, `value_median`, `value_q3` | numeric distribution stats (null for boolean assays) |
| `MetaboliteAssay` | `flag_true_count`, `flag_false_count` | boolean flag counts (null for numeric assays) |
| `MetaboliteAssay` | `growth_phases` (str[]) | derived from parent Experiment |
| `Metabolite` | `measured_assay_count`, `measured_organisms`, `measured_paper_count` | new, derived from incoming Assay edges |
| `Metabolite` | `organism_count`, `organism_names` | UNION extended to include measurement path |
| `Metabolite` | `gene_count`, `transporter_count`, `pathway_*` | unchanged |
| `Organism_has_metabolite` | edge materialization | extended to measurement-only pairs |
| `Organism_has_metabolite` | `evidence_sources` (str[]) | new — values `metabolism` \| `transport` \| `measured` |
| `Organism_has_metabolite` | `measured_assay_count`, `measured_compartments`, `measured_paper_count` | new augmentations (0 / [] for non-measured edges) |
| `Experiment` | `metabolite_assay_count`, `metabolite_compartments`, `metabolite_count` | new |
| `Publication` | `metabolite_assay_count`, `metabolite_compartments`, `metabolite_count` | new |
| `OrganismTaxon` | `measured_metabolite_count` | new |
| `Gene.*` | — | unchanged |

### 1.6 — Validator + skill updates

- `tools/validate_paperconfig.py` gains a `metabolite_assays_table` schema (see §3.4) and a `--report-backlog` mode (scans for `# BACKLOG:` markers across all paperconfigs)
- `.claude/skills/paperconfig/SKILL.md` gets a metabolomics section + extraction templates
- `CLAUDE.md` "paperconfig supplementary_materials entry types" table gains the new row

## 2 — Schema additions to `config/schema_config.yaml`

Closely mirrors the DerivedMetric block (lines 119–277 of current schema_config.yaml).

### 2.1 — New node type `metabolite assay`

```yaml
# MetaboliteAssay node: one per (Experiment × compartment × value_kind).
# Carries assay-level metabolomics semantics. Edge type emitted to Metabolite
# depends on value_kind. Adapter emits in metabolite_assay_adapter.py;
# post-import computes total_metabolite_count + distribution stats + growth_phases.
metabolite assay:
  is_a: information content entity
  represented_as: node
  preferred_id: metabolite_assay_id
  label_in_input: metabolite_assay
  properties:
    name: str
    experiment_id: str
    organism_name: str
    publication_doi: str
    compartment: str
    omics_type: str
    treatment_type: str[]
    background_factors: str[]
    treatment: str
    light_condition: str
    experimental_context: str
    metric_type: str
    value_kind: str                  # "numeric" | "boolean"
    unit: str
    rankable: str                    # "true" | "false"
    aggregation_method: str          # "mean_across_replicates" | "pre_aggregated" | "single_measurement"
    field_description: str
    total_metabolite_count: int      # post-import
    growth_phases: str[]             # post-import: from parent Experiment
    # Numeric distribution stats (post-import; null for boolean assays)
    value_min: float
    value_max: float
    value_q1: float
    value_median: float
    value_q3: float
    # Boolean flag counts (post-import; null for numeric assays)
    flag_true_count: int
    flag_false_count: int
```

Differences from `derived metric`:
- Added: `aggregation_method`
- Removed: `has_p_value`, `p_value_threshold`, `allowed_categories` (no significance variant in v1)
- `total_gene_count` → `total_metabolite_count`

### 2.2 — Binding edges

```yaml
publication has metabolite assay:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: publication_has_metabolite_assay
  source: publication
  target: metabolite assay

experiment has metabolite assay:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: experiment_has_metabolite_assay
  source: experiment
  target: metabolite assay

metabolite assay belongs to organism:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: metabolite_assay_belongs_to_organism
  source: metabolite assay
  target: organism taxon
```

### 2.3 — Measurement edges

```yaml
# MetaboliteAssay → Metabolite (value_kind="numeric"; replicate-aggregated)
assay quantifies metabolite:
  is_a: Association
  represented_as: edge
  label_as_edge: assay_quantifies_metabolite
  label_in_input: assay_quantifies_metabolite
  source: metabolite assay
  target: metabolite
  properties:
    metric_type: str
    condition_label: str
    time_point: str                  # nullable
    time_point_order: int            # nullable
    time_point_hours: float          # sentinel -1.0 when unknown
    value: float                     # mean across replicates
    value_sd: float                  # stdev across replicates
    n_replicates: int
    n_non_zero: int
    replicate_values: float[]        # raw per-replicate values, ordered
    detection_status: str            # post-import
    rank_by_metric: int              # post-import (rankable="true" only)
    metric_percentile: float         # post-import (rankable="true" only)
    metric_bucket: str               # post-import (rankable="true" only)

# MetaboliteAssay → Metabolite (value_kind="boolean"; presence/absence)
assay flags metabolite:
  is_a: Association
  represented_as: edge
  label_as_edge: assay_flags_metabolite
  label_in_input: assay_flags_metabolite
  source: metabolite assay
  target: metabolite
  properties:
    metric_type: str
    condition_label: str             # nullable for whole-paper flags
    flag_value: str                  # "true" | "false"
    n_replicates: int                # default 1 when source is whole-paper
    n_positive: int                  # default = n_replicates when flag_value="true"
```

### 2.4 — Property additions on existing entries

**`metabolite` node** — three new computed properties (post-import):
```yaml
    measured_assay_count: int
    measured_organisms: str[]
    measured_paper_count: int
```
Existing `evidence_sources: str[]` gains `"metabolomics"` as a possible value (no schema change). `organism_count` and `organism_names` semantics extend to include the measurement path (post-import logic only, no schema change).

**`organism has metabolite` edge** — four new properties:
```yaml
    evidence_sources: str[]          # values: metabolism | transport | measured
    measured_assay_count: int
    measured_compartments: str[]
    measured_paper_count: int
```

**`experiment` node** — three new properties:
```yaml
    metabolite_assay_count: int
    metabolite_compartments: str[]
    metabolite_count: int
```

**`publication` node** — same three:
```yaml
    metabolite_assay_count: int
    metabolite_compartments: str[]
    metabolite_count: int
```

**`organism taxon` node** — one:
```yaml
    measured_metabolite_count: int
```

### 2.5 — Compartment vocab additions

Existing values: `whole_cell` (default), `vesicle`, `exoproteome`, `secretome`. Phase 2 adds: `intracellular`, `extracellular`, `exometabolome`, `column_extract`. Validator's enum updated accordingly.

## 3 — paperconfig entry shape

### 3.1 — `metabolite_assays_table` entry type

```yaml
<entry_key>:
  type: metabolite_assays_table
  filename: "data/.../source.csv"

  # ID resolution columns
  name_col: "Compound Name"          # column holding metabolite name (always required)
  id_col: ""                         # column holding pre-mapped ID (optional; e.g. "KEGG_ID")
  id_type: kegg.compound             # "kegg.compound" | "chebi" | "mnx" (only meaningful if id_col set)
  formula_col: ""                    # optional; molecular formula column (used as disambiguator only)

  # MS-feature placeholders (v0 ignored; reserved for future MS-matching)
  mz_col: ""
  retention_time_col: ""
  adduct_col: ""

  # Header / cell parsing
  skip_rows: 0                       # int; rows to skip before header
  header_rows: []                    # list of row indices forming a multi-row header (empty = single header at skip_rows offset)
  null_values: ["nd", "NA", "n.d."]  # cells meaning "not detected" → contribute n_non_zero=0
  missing_values: [""]               # cells meaning "not measured" → exclude from aggregation
  cell_format: numeric               # "numeric" | "embedded_mean_sd_n"
  extra_columns_to_edge: []          # list of column names whose value gets copied to each emitted edge

  # Aggregation
  aggregation_method: mean_across_replicates    # node-level default; per-assay override allowed
  drop_undetected: false             # default false: keep all-zero condition edges

  # Compound-class hint (optional; routes resolution differently for lipids, etc.)
  compound_class: ""                 # "" | "lipid" | "peptide" | ...

  # Manual ID curation
  aliases_file: "metabolite_aliases.yaml"   # optional, relative to paper directory

  # Cross-paper ID matching (placeholder; v0 ignored, v1+ activated)
  id_columns: []                     # list[{column, id_type, tier}]; mirrors gene-side id_translation

  # The actual assays
  assays:
    - id: <assay_key>                       # unique within entry
      experiment: <experiment_key>          # references paperconfig.experiments
      organism: "Prochlorococcus MIT9303"
      compartment: intracellular
      metric_type: intracellular_concentration
      value_kind: numeric                   # "numeric" | "boolean"
      unit: "fg/cell"
      rankable: true
      field_description: "intracellular concentration in fg/cell, blank-corrected"
      aggregation_method: mean_across_replicates    # optional override
      sample_columns:
        # Numeric variant: aggregates a list of replicate columns into one edge per metabolite
        - condition_label: control
          time_point: "T=4"
          time_point_order: 1
          time_point_hours: -1.0
          replicate_columns:
            - "9303 Control T=4, replicate 1"
            - "9303 Control T=4, replicate 2"
            - "9303 Control T=4, replicate 3"

        # Boolean variant: one column → one edge per metabolite, flag_true_value matches
        - condition_label: ""
          flag_column: "intracellular"
          flag_true_value: "yes"
```

### 3.2 — Per-paper `metabolite_aliases.yaml`

Lives next to the paper's `paperconfig.yaml`. Free-text key → primary ID:

```yaml
# data/Prochlorococcus/papers_and_supp/Kujawinski 2023/metabolite_aliases.yaml
"γ-aminobutyric acid": "kegg.compound:C00334"
"GABA": "kegg.compound:C00334"
"AmMP": "kegg.compound:C20267"
"4-amino-5-aminomethyl-2-methylpyrimidine": "kegg.compound:C20267"
"2-3-dihydroxypropane1sulfonate": "kegg.compound:C19675"   # paper had a typo; missing hyphen
```

Step 6 reads aliases during the harvest phase and merges them into `metabolite_id_mapping.json`. Per-paper overrides win over MNX resolver hits when both produce a result. Comments in the YAML capture *why* the override exists.

### 3.3 — Resolution-method enum

Step 7 writes `resolution_method` per row. Initial values:

| Value | Trigger |
|---|---|
| `kegg_direct` | `id_col` cell parsed cleanly |
| `chebi_direct` / `mnx_direct` | id_col with non-KEGG `id_type` |
| `name_match` | MNX resolver hit on `name_col` (single candidate) |
| `name_match_formula_disambiguated` | multiple name hits, narrowed to one by `formula_col` |
| `alias_override` | matched a key in `metabolite_aliases.yaml` |
| `fuzzy_match` | resolver returned a best-effort match (similarity above `--fuzzy-similarity-min` threshold, default 0.85) |
| `unresolved` | none of the above |
| `formula_only_unresolved` | formula present but no name; resolution declined |

Future additions expected: `class_label_drop`, `manual_curation`, `isomer_collapsed`, `charge_state_normalized`, `ambiguous_multi_id`. Treat the enum as iterative.

### 3.4 — Validator schema for `metabolite_assays_table`

`validate_metabolite_assays_table(entry)` enforces:

- Required keys: `type`, `filename`, `name_col`, `assays`
- `assays` non-empty list; each entry must have `id`, `experiment`, `organism`, `compartment`, `metric_type`, `value_kind`, `sample_columns`
- `value_kind ∈ {"numeric", "boolean"}`
- `compartment ∈ {whole_cell, intracellular, extracellular, exometabolome, column_extract, vesicle, exoproteome, secretome}`
- `aggregation_method ∈ {mean_across_replicates, pre_aggregated, single_measurement}`
- `cell_format ∈ {numeric, embedded_mean_sd_n}`
- For `value_kind: numeric` — every `sample_columns` entry must have `replicate_columns: list[str]` (length ≥ 1)
- For `value_kind: boolean` — every `sample_columns` entry must have `flag_column: str` and `flag_true_value: str`
- `experiment` must reference a key in this paperconfig's `experiments:` block
- `organism` must match a known OrganismTaxon `preferred_name`
- `aliases_file` (if set) must exist and be valid YAML mapping `str → str`
- `formula_col` (if set) must reference a CSV column

### 3.5 — Resolution complications observed in biller 2022

biller 2022 surfaces complications hidden by Capovilla and Kujawinski. Catalogued:

| Complication | Example | v1 stance |
|---|---|---|
| **Mass features without compound IDs** (Table S6) | 4065 rows like `"47_CyanoAq", m/z=120.08111, RT=1.139` — no compound name | Out of scope for v1; would need a `MassFeature` node type. |
| **Embedded-uncertainty cells** (Table S7) | `"0.00054 (8.8e-05), n=2"` — one cell encodes mean + sd + n_replicates | In scope. `cell_format: embedded_mean_sd_n` parser. |
| **`nd` / `NA` / blank distinct semantics** | "nd" = explicitly not detected; "NA" = couldn't compute sd; blank = not applicable | In scope. `null_values` (treat as not-detected → n_non_zero=0) vs `missing_values` (exclude from aggregation). |
| **Lipid-class shorthand** (Table S1) | `"DGDG-30:1-1OH"` — LipidMaps notation, often not in KEGG/MNX | Accept low resolution rate; `compound_class: lipid` hint reserved; per-paper `metabolite_aliases.yaml` is the manual fallback. |
| **Multi-row headers** | Row 4 = strain; row 5 = compartment; data row 6+ | In scope. `header_rows: [4, 5]` lets the user spell it out. |
| **Per-row analytical fraction** | "Adenine: HILICPos", "Aconitic Acid: HILICNeg" | In scope. `extra_columns_to_edge: ["Analytical Fraction"]` carries it through to the edge. |
| **Variable per-cell replicate count** | `n=1`, `n=2`, `n=3` cells in same paper | In scope (the `embedded_mean_sd_n` parser extracts per-cell n). |
| **Formula degeneracy** | `C6H12O6` matches glucose, fructose, galactose, … | Formula-alone is unresolvable; formula-as-disambiguator is in scope (`formula_col`). |

### 3.6 — Cross-paper resolution: placeholder + iteration plan

Gene-side build_gene_id_mapping.py harvests ID columns from all paperconfigs, builds a transitive closure linking IDs across papers via shared genes. For metabolites:

- **Structural cross-references** (KEGG ↔ ChEBI ↔ MNX ↔ BiGG ↔ HMDB ↔ MetaCyc ↔ ModelSEED ↔ Reactome ↔ Rhea ↔ LipidMaps ↔ SwissLipids ↔ VMH ↔ enviPath) — already encoded in the MNX resolver. Free.
- **Cross-paper synonym pooling** — step 6 walks **all** paperconfigs and merges aliases into one `metabolite_id_mapping.json`. Curated overrides in paper A automatically benefit paper B.
- **Paper-specific internal IDs** (placeholder for later) — `id_columns: []` on the entry; v0 ignored; v1+ reserves the shape.

`metabolite_id_mapping.json` schema reserves a three-tier structure mirroring `gene_id_mapping.json`:

```jsonc
{
  "specific_lookup":  { ... },   // tier 1, 1:1, conflicts = data errors
  "multi_lookup":     { ... },   // tier 2, singletons usable
  "name_lookup":      { ... },   // tier 3, singletons usable, ambiguity expected
  "conflicts":        { ... },   // tier-1 collisions for human review
  "compounds":        { ... }    // per-compound alias catalogue
}
```

**v1 fills `name_lookup` only.** `specific_lookup` and `multi_lookup` activate when papers introduce ID columns. `conflicts` activates immediately if step 6 sees two papers map the same name to incompatible primaries.

## 4 — Post-import Cypher

Mirrors `scripts/post-import.cypher` pattern: ~3 grouped invocations (indexes / small aggregations / `CALL IN TRANSACTIONS` writes). All snippets must remain byte-identical between `scripts/post-import.cypher` and `scripts/post-import.sh` per the existing CLAUDE.md invariant.

### 4.1 — Indexes

```cypher
CREATE INDEX metabolite_assay_organism_idx     IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.organism_name);
CREATE INDEX metabolite_assay_compartment_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.compartment);
CREATE INDEX metabolite_assay_metric_type_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.metric_type);
CREATE INDEX metabolite_assay_value_kind_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.value_kind);
CREATE INDEX metabolite_assay_experiment_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.experiment_id);

CREATE FULLTEXT INDEX metaboliteAssayFullText  IF NOT EXISTS FOR (a:MetaboliteAssay)
  ON EACH [a.name, a.field_description, a.treatment, a.experimental_context];
```

### 4.2 — Edge-level computed properties

```cypher
// detection_status (every numeric edge)
MATCH ()-[r:Assay_quantifies_metabolite]->()
CALL { WITH r
  SET r.detection_status =
    CASE
      WHEN r.n_non_zero = 0                  THEN 'not_detected'
      WHEN r.n_non_zero = r.n_replicates     THEN 'detected'
      ELSE 'sporadic'
    END
} IN TRANSACTIONS OF 50000 ROWS;

// rank_by_metric / metric_percentile / metric_bucket — per assay, only when rankable=true
MATCH (a:MetaboliteAssay {rankable:'true'})-[r:Assay_quantifies_metabolite]->()
WITH a, r ORDER BY a.id, r.value DESC
WITH a, collect(r) AS edges, count(r) AS total
UNWIND range(0, size(edges)-1) AS i
WITH edges[i] AS r, i+1 AS rk, total
SET r.rank_by_metric = rk,
    r.metric_percentile = 100.0 * (1.0 - (rk - 1.0) / total),
    r.metric_bucket = CASE
      WHEN 100.0 * (1.0 - (rk - 1.0) / total) >= 90 THEN 'top_decile'
      WHEN 100.0 * (1.0 - (rk - 1.0) / total) >= 75 THEN 'top_quartile'
      WHEN 100.0 * (1.0 - (rk - 1.0) / total) <  25 THEN 'low'
      ELSE 'mid'
    END;
```

### 4.3 — `MetaboliteAssay` node rollups

```cypher
MATCH (a:MetaboliteAssay)
OPTIONAL MATCH (a)-[r:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH a, count(DISTINCT m) AS cnt
SET a.total_metabolite_count = cnt;

// Numeric distribution stats
MATCH (a:MetaboliteAssay {value_kind:'numeric'})-[r:Assay_quantifies_metabolite]->()
WITH a,
     min(r.value) AS vmin, max(r.value) AS vmax,
     percentileDisc(r.value, 0.25) AS q1,
     percentileDisc(r.value, 0.5)  AS med,
     percentileDisc(r.value, 0.75) AS q3
SET a.value_min=vmin, a.value_max=vmax, a.value_q1=q1, a.value_median=med, a.value_q3=q3;

// Boolean flag counts
MATCH (a:MetaboliteAssay {value_kind:'boolean'})-[r:Assay_flags_metabolite]->()
WITH a,
     sum(CASE WHEN r.flag_value='true'  THEN 1 ELSE 0 END) AS t,
     sum(CASE WHEN r.flag_value='false' THEN 1 ELSE 0 END) AS f
SET a.flag_true_count=t, a.flag_false_count=f;

// growth_phases from parent Experiment
MATCH (a:MetaboliteAssay)<-[:ExperimentHasMetaboliteAssay]-(e:Experiment)
SET a.growth_phases = coalesce(e.growth_phases, []);
```

### 4.4 — `Metabolite` rollups

```cypher
// New metabolomics-only properties
MATCH (m:Metabolite)
OPTIONAL MATCH (m)<-[:Assay_quantifies_metabolite|Assay_flags_metabolite]-(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
WITH m,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT o.preferred_name) AS orgs,
     count(DISTINCT p) AS pcnt
SET m.measured_assay_count = acnt,
    m.measured_organisms   = apoc.coll.sort([x IN orgs WHERE x IS NOT NULL]),
    m.measured_paper_count = pcnt;

// Extend existing organism_count / organism_names UNION (modify the existing block)
MATCH (m:Metabolite)
OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(:TcdbFamily {level_kind:'tc_specificity'})<-[:Gene_has_tcdb_family]-(:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
OPTIONAL MATCH (m)<-[:Assay_quantifies_metabolite|Assay_flags_metabolite]-(:MetaboliteAssay)-[:MetaboliteAssayBelongsToOrganism]->(o3:OrganismTaxon)
WITH m, apoc.coll.toSet([o IN collect(DISTINCT o1) + collect(DISTINCT o2) + collect(DISTINCT o3) WHERE o IS NOT NULL | o.preferred_name]) AS names
SET m.organism_count = size(names),
    m.organism_names = apoc.coll.sort(names);
```

### 4.5 — `Organism_has_metabolite` materialization extension + augmentation

```cypher
// NEW measurement-path materialization (alongside existing metabolism + transport blocks)
MATCH (a:MetaboliteAssay)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH DISTINCT o, m
CALL { WITH o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['measured']
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'measured' IN r.evidence_sources
         THEN r.evidence_sources
         ELSE r.evidence_sources + 'measured' END
} IN TRANSACTIONS OF 10000 ROWS;

// Existing metabolism + transport blocks gain ON CREATE/ON MATCH evidence_sources
// (refactor existing materialization to use the same ON CREATE/ON MATCH pattern with 'metabolism'/'transport' values)

// Augmentation properties on every edge (set 0/[] for non-measured)
MATCH (o:OrganismTaxon)-[r:Organism_has_metabolite]->(m:Metabolite)
OPTIONAL MATCH (o)<-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m)
OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
WITH r,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT p) AS pcnt
SET r.measured_assay_count   = acnt,
    r.measured_compartments  = [c IN comps WHERE c IS NOT NULL],
    r.measured_paper_count   = pcnt;
```

### 4.6 — Experiment / Publication / OrganismTaxon rollups

```cypher
MATCH (e:Experiment)-[:ExperimentHasMetaboliteAssay]->(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH e,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT m) AS mcnt
SET e.metabolite_assay_count   = acnt,
    e.metabolite_compartments  = [c IN comps WHERE c IS NOT NULL],
    e.metabolite_count         = mcnt;

MATCH (p:Publication)-[:PublicationHasMetaboliteAssay]->(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH p,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT m) AS mcnt
SET p.metabolite_assay_count   = acnt,
    p.metabolite_compartments  = [c IN comps WHERE c IS NOT NULL],
    p.metabolite_count         = mcnt;

MATCH (o:OrganismTaxon)<-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH o, count(DISTINCT m) AS mcnt
SET o.measured_metabolite_count = mcnt;
```

### 4.7 — Defaults for absent rollups

```cypher
MATCH (e:Experiment) WHERE e.metabolite_assay_count IS NULL
SET e.metabolite_assay_count = 0,
    e.metabolite_compartments = [],
    e.metabolite_count = 0;

MATCH (p:Publication) WHERE p.metabolite_assay_count IS NULL
SET p.metabolite_assay_count = 0,
    p.metabolite_compartments = [],
    p.metabolite_count = 0;

MATCH (o:OrganismTaxon) WHERE o.measured_metabolite_count IS NULL
SET o.measured_metabolite_count = 0;

MATCH (m:Metabolite) WHERE m.measured_assay_count IS NULL
SET m.measured_assay_count = 0,
    m.measured_organisms = [],
    m.measured_paper_count = 0;
```

### 4.8 — Validation harness extensions

`scripts/post-import-validate.sh` adds parallel queries for every property added in §1.5 + §4. Output stays deterministic so `diff baseline.txt after.txt` is byte-clean for refactor-only changes.

## 5 — Test plan, acceptance criteria, commit gates

### 5.1 — Acceptance criteria (binary pass/fail)

| ID | Criterion |
|---|---|
| A1 | Both validation papers integrated: Capovilla 2023 + Kujawinski 2023 emit MetaboliteAssay nodes + edges in the production KG. |
| A2 | Every property in §1.5 + §2 is present on the deployed graph; types match. |
| A3 | `omics-edge-snapshot` taken before and after shows zero net loss in `Changes_expression_of` and `Derived_metric_*` edges. Per-paper deltas all 0 or positive. |
| A4 | Each `<stem>_resolution_report.json` reports match-rate; both validation papers ≥ 70% resolved (Capovilla expected ~85–90%; Kujawinski expected ~95% via KEGG IDs). |
| A5 | All pre-existing tests pass: `pytest -m "not slow and not kg"` and `pytest -m kg` both green. |
| A6 | New tests pass (see §5.2). |
| A7 | `scripts/post-import-validate.sh` baseline diff is byte-clean for refactor commits. |
| A8 | No orphan MetaboliteAssay nodes: every assay has all 3 binding edges (Publication, Experiment, Organism); every measurement edge points to an existing Metabolite. |
| A9 | Every row in `<stem>_resolved.csv` has either `metabolite_id` populated or `resolution_method = unresolved`. |

### 5.2 — Test additions

**Unit tests** (`tests/`, `pytest -m "not slow and not kg"`):
- `test_paperconfig_utils.py` — `iter_metabolite_assays_tables` filtering + schema validation
- `test_metabolite_assay_adapter.py` — single-paper adapter emits expected nodes/edges; replicate aggregation correctness; `embedded_mean_sd_n` parser; boolean flag parser; `null_values` vs `missing_values` distinction
- `test_resolve_paper_metabolites.py` — step 7 CSV rewriter; `resolution_method` enum coverage
- `test_build_kegg_metabolism_xrefs.py` — paper metabolite harvesting; alias merging; `evidence_sources` union; pathway extension
- `test_validate_paperconfig.py` — new `metabolite_assays_table` validation rules

**KG validity tests** (`tests/kg_validity/`, `pytest -m kg`):
- `test_metabolomics.py` (NEW):
  - MetaboliteAssay nodes present (≥1 per validation paper × ≥1 per compartment-assay split)
  - `value_kind ∈ {"numeric", "boolean"}`
  - `aggregation_method` enum check
  - `compartment` ∈ approved vocab
  - All assays have all 3 binding edges
  - `Assay_quantifies_metabolite` edges: `value`, `value_sd`, `n_replicates`, `n_non_zero`, `replicate_values` all present and consistent (n_replicates == size(replicate_values), n_non_zero ≤ n_replicates, mean(replicate_values) ≈ value within tolerance)
  - `detection_status` ∈ {detected, sporadic, not_detected} and consistent
  - `Assay_flags_metabolite` edges: `flag_value ∈ {"true", "false"}`
  - `MetaboliteAssay.total_metabolite_count > 0`; numeric distribution stats present iff value_kind=numeric; flag counts present iff value_kind=boolean
  - Spot checks: a known-resolved compound measured in MIT9303 via Capovilla; presence flag on `2,3-dihydroxypropane-1-sulfonate` → extracellular in Kujawinski
- `test_structure.py` — extend orphan checks to MetaboliteAssay; verify Metabolite nodes for paper-introduced compounds with `evidence_sources` containing `"metabolomics"`
- `test_organism.py` — verify `OrganismTaxon.measured_metabolite_count` populated for organisms with assays; defaults to 0 elsewhere
- `test_snapshot.py` — regenerate snapshot fixture
- `test_post_import.py` — verify all new computed properties present; verify `Organism_has_metabolite.evidence_sources` array contents; verify augmentation defaults on non-measured edges

**Spot-check Cypher templates** added to `.claude/skills/cypher-queries/SKILL.md`:
- "List all metabolite assays per paper"
- "List metabolites measured in organism X with detection frequency"
- "Find compounds measured in both metabolomics and predicted by metabolism (gene-reachable)"

### 5.3 — Commit gate sequence

8 commits, each with its own validation. Mirrors the prior `non-de-evidence-extension` plan structure.

| # | Commit | Scope | Validation |
|---|--------|-------|------------|
| 1 | **Schema + paperconfig_utils + validator** | schema_config.yaml additions; `iter_metabolite_assays_tables`; validate_paperconfig.py rules. No data change. | unit tests pass; existing KG unchanged |
| 2 | **Step 6 extension + step 7** | `build_kegg_metabolism_xrefs.py` walks paperconfigs for metabolite_assays_tables; writes metabolite_id_mapping.json + extends kegg_data.json. New `resolve_paper_metabolites.py` writes `_resolved.csv` files. prepare_data.sh adds step 7. | unit tests; manual diff of kegg_data.json shows new entries; resolved CSVs exist |
| 3 | **metabolism_adapter verification** | Confirm new Metabolite nodes appear in KG; pathway materialization picks them up. Likely no code change. | KG rebuild; new Metabolite nodes for Kujawinski-only compounds visible |
| 4 | **MetaboliteAssayAdapter + create_knowledge_graph wiring** | New adapter file; small fix in omics_adapter for unconditional experiments iteration. No paperconfig populated yet. | unit tests; KG rebuild produces zero new MetaboliteAssay nodes; pre-existing tests green |
| 5 | **Post-import Cypher + validate.sh extensions** | scripts/post-import.cypher + post-import.sh additions; post-import-validate.sh extensions. No KG content yet. | post-import-validate.sh baseline diff shows ONLY new property enumerations |
| 6 | **Capovilla 2023 metabolomics integration** | Capovilla paperconfig.yaml gets `metabolite_assays_table` entries; experiments block reused; metabolite_aliases.yaml as needed. | Step 6/7 produce non-empty resolution; KG build emits MetaboliteAssay + edges; KG validity tests pass; omics-edge-snapshot shows no DE regression |
| 7 | **Kujawinski 2023 metabolomics integration** | Kujawinski paperconfig.yaml created from scratch; experiments block declared; numeric (KEGG export) + boolean (s2 presence flags) entries; aliases for paper-internal abbreviations; **MIT0801 entries pre-written + commented out** per §5.5 | Step 6/7 with ≥95% resolution; boolean flag edges visible; cross-strain MIT0801 row dropped (organism not in KG) with logged skip |
| 8 | **Docs + skills** | CLAUDE.md updates; `.claude/skills/paperconfig/SKILL.md` metabolomics section; `docs/kg-changes/metabolomics-extension.md` changelog | manual review; tests still green |

Each commit includes its own validation step; gate fails → previous commit reverts, not this one.

### 5.4 — Risks / open items

- **Resolution-rate uncertainty for Capovilla compound names** — table headers should resolve cleanly via MNX. Worst case → small `metabolite_aliases.yaml`. No design change.
- **Performance** — `Organism_has_metabolite` materialization extension adds a third sub-query path; verify with `EXPLAIN` that it stays under existing post-import budget. Likely fine; metabolites are small (~1.5K nodes).
- **Snapshot regeneration cadence** — `test_snapshot.py` will need re-baselining after commit 6 + commit 7. Standard procedure.

### 5.5 — MIT0801 backlog + commented-out paperconfig pattern

**Backlog file** (`plans/strain_deployment_backlog.md`, new — created in commit 7) tracks strains referenced by paperconfigs but not yet deployed. Initial entry:

```markdown
# Strain deployment backlog

Strains referenced by paperconfigs whose data remains uningested pending KG deployment.
Track via `/deploy-strain` workflow when picked up.

## MIT0801 (Prochlorococcus, LLI ecotype)
- **Source NCBI taxid:** 1110371 (or similar — verify via NCBI Assembly)
- **Referenced by:**
  - Kujawinski 2023 — 3 sample columns in `ChisholmPro_cellSpecific_KEGGexport.csv`
    (replete_extracellular_s0801ax_10, replete_filter_s0801ax_10) + glycine-betaine row in table s3
- **Blockers:** strain not in `cyanobacteria_genomes.csv`; no genome assembly downloaded.
- **paperconfig state:** Kujawinski 2023 paperconfig.yaml has MIT0801 entries pre-written, commented out.
  Re-enable on strain deployment via uncomment + adapter rebuild.
- **Estimated work:** ~1 day per `/deploy-strain` skill flow.
```

Future strains added below as they surface.

**Commented-out paperconfig pattern.** Paperconfigs reference all strains the paper covers, including undeployed ones; entries that depend on undeployed strains are YAML-commented with a `# BACKLOG:` marker:

```yaml
publication:
  experiments:
    kujawinski_metabolomics_9301:
      organism: "Prochlorococcus MIT9301"
      ...
    kujawinski_metabolomics_9313:
      organism: "Prochlorococcus MIT9313"
      ...
    # BACKLOG: requires MIT0801 deployment — see plans/strain_deployment_backlog.md
    # kujawinski_metabolomics_0801:
    #   organism: "Prochlorococcus MIT0801"
    #   ...

  supplementary_materials:
    metabolites_kegg_export:
      type: metabolite_assays_table
      ...
      assays:
        - id: kegg_export_9301_intracellular
          ...
        # BACKLOG: requires MIT0801 deployment
        # - id: kegg_export_0801_intracellular
        #   experiment: kujawinski_metabolomics_0801
        #   ...
```

Why YAML comments rather than a `disabled: true` field:
- Zero adapter complexity — comments are stripped by the YAML parser before validation runs
- `grep BACKLOG paperconfig.yaml` lists pending work directly
- Re-enabling is a simple uncomment, no field changes
- The validator + skill flows already see only the live entries

`tools/validate_paperconfig.py` gains a `--report-backlog` flag that scans for `# BACKLOG:` markers across all paperconfigs and prints them per file, keeping `strain_deployment_backlog.md` synced. The `/paperconfig` skill SKILL.md documents the convention.

## Out of scope (explicit)

- Paired-comparison FC/p-value edges — no current paper has them; if they appear, design a separate node/edge type
- `MassFeature` node type for unidentified peaks (Table S6 of biller 2022)
- Isotopologue tracing (Szul 2019)
- `Assay_classifies_metabolite` (no current paper needs categorical metabolite values)
- Auto-pivot of long-format CSVs (paperconfig is explicit about column → assay mapping)
- Per-replicate edges as a separate type (replicates are aggregated into one edge with `replicate_values` array)
- `/check-metabolite-ids` skill (the step-7 resolution report is sufficient for v1)
- MS-feature mass-window matching (m/z + RT placeholders reserved; deferred)
- `UnresolvedMetabolite` node type — unresolved rows produce no edge
- Lipid-class hierarchy (LipidMaps integration) — `compound_class: lipid` hint reserved

## Open questions to revisit after v1 lands

1. **Lipid-class resolution rate** — biller 2022 Table S1 lipids are LOBSTAHS notation. If most stay `unresolved`, evaluate adding a LipidMaps-aware resolver path or a `LipidClass` node type.
2. **Per-replicate query patterns** — if downstream queries frequently want per-replicate data and the `replicate_values` array isn't ergonomic, consider an opt-in per-replicate edge type.
3. **Three-tier metabolite_id_mapping** — activate `specific_lookup` and `multi_lookup` when papers introduce ID columns.
4. **Compartment vocab churn** — Phase 2 adds `intracellular`, `extracellular`, `exometabolome`, `column_extract`. Future papers may need more (e.g. `surface_DOM`, `dissolved_inorganic_pool`).
5. **Step 6 performance budget** — paperconfig walking + per-row resolver hits should stay well under existing step 6 budget. If degradation observed, batch resolver hits or precompute a once-per-build alias index.
