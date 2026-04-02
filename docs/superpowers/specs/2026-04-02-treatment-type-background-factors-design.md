# Treatment Type & Background Factors Redesign

## Problem

Experiment nodes currently have `treatment_type` as a single string. This loses two kinds of information:

1. **Multi-factor experiments** cannot be fully described. A coculture darkness experiment (Biller 2018) is tagged only `darkness` -- invisible when querying for coculture-related experiments.
2. **Background conditions** (light regime, culture status) are buried in free-text `experimental_context` and not queryable. Thompson 2016 tests viral infection in both light and dark, but you can't filter by light context.

## Design

### Two categorical list fields

| Field | Type | Required | Rule |
|---|---|---|---|
| `treatment_type` | `str[]` | yes (>=1 value) | What the DE comparison tests (treatment vs control) |
| `background_factors` | `str[]` | no (default `[]`) | Conditions present but not compared in the DE analysis. Axenic/coculture always goes here unless coculture IS the DE comparison |

Both fields use the same shared canonical vocabulary, validated in `validate_paperconfig.py`.

### Shared canonical vocabulary

Existing values (unchanged):
- `nitrogen_stress`, `phosphorus_stress`, `iron_stress`, `carbon_stress`, `salt_stress`, `oxygen_stress`, `temperature_stress`, `light_stress`, `darkness`, `plastic_stress`, `viral`, `coculture`, `growth_state`, `growth_medium`, `diel`

New values:
- `axenic` -- organism grown alone (not in coculture)
- `continuous_light` -- standard continuous illumination
- `diel_cycle` -- light:dark cycling regime

### Node changes

**Experiment node:**
- `treatment_type`: `str` -> `str[]`
- `background_factors`: new `str[]`

**ClusteringAnalysis node:**
- `treatment_type`: already `str[]` (unchanged)
- `background_factors`: new `str[]`

**GeneCluster node:**
- `treatment_type`: already `str[]` (unchanged)
- `background_factors`: new `str[]`

**Publication node (post-import computed):**
- `treatment_types`: unchanged (already `str[]`, aggregated from experiments)
- `background_factors`: new `str[]`, aggregated from experiments

**OrganismTaxon node (post-import computed):**
- `treatment_types`: unchanged (already `str[]`)
- `background_factors`: new `str[]`, aggregated from publications

### Classification rule

`treatment_type` = what is directly compared in the differential expression analysis.
`background_factors` = experimental conditions that are present and could affect interpretation, but are NOT the variable being compared.

**Special rule for axenic/coculture:** These always go in `background_factors` to indicate culture status, UNLESS `coculture` is itself the DE comparison (e.g., Aharonovich 2016 where DE tests coculture vs axenic).

### Examples

**Thompson 2016** (2x2 factorial: viral x light/dark):

| Experiment | treatment_type | background_factors |
|---|---|---|
| dark vs light, while infected | `[light_stress]` | `[viral]` |
| dark vs light, uninfected | `[light_stress]` | `[axenic]` |
| phage vs spent medium, in light | `[viral]` | `[continuous_light]` |
| phage vs spent medium, in dark | `[viral]` | `[darkness]` |

**Biller 2018** (darkness in axenic vs coculture):

| Experiment | treatment_type | background_factors |
|---|---|---|
| NATL2A darkness, axenic | `[darkness]` | `[axenic]` |
| NATL2A darkness, coculture | `[darkness]` | `[coculture]` |
| MIT1002 darkness, coculture | `[darkness]` | `[coculture]` |

**Weissberg 2025** (N-starvation in axenic vs coculture):

| Experiment | treatment_type | background_factors |
|---|---|---|
| MED4 N-starvation, axenic | `[nitrogen_stress]` | `[axenic]` |
| MED4 N-starvation, coculture | `[nitrogen_stress]` | `[coculture]` |
| HOT1A3 coculture effect | `[coculture]` | -- |

**Hennon 2017** (carbon stress under diel cycle):

| Experiment | treatment_type | background_factors |
|---|---|---|
| MIT9312 elevated CO2, coculture | `[carbon_stress]` | `[coculture, diel_cycle]` |
| EZ55 elevated CO2, coculture | `[carbon_stress]` | `[coculture, diel_cycle]` |
| EZ55 elevated CO2, axenic | `[carbon_stress]` | `[axenic, diel_cycle]` |

### Paperconfig format

```yaml
experiments:
  my_experiment:
    name: "..."
    organism: "..."
    treatment_type: [nitrogen_stress]        # list (required, >=1 value)
    background_factors: [axenic, continuous_light]  # list (optional)
    treatment_condition: "..."
    control_condition: "..."
    # ... rest unchanged
```

Gene clusters entry:

```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  treatment_type: [nitrogen_stress]          # already a list
  background_factors: [axenic, continuous_light]  # new
  # ... rest unchanged
```

### Validation

- Both fields validated against the same `CANONICAL_CONDITION_TYPES` set
- `treatment_type` required (>= 1 value); `background_factors` optional (defaults to `[]`)
- A value appearing in both fields for the same experiment triggers a warning (not error)
- Validator accepts both string and list for `treatment_type` (normalizes to list)

### Indexes

- Existing `experiment_treatment_type_idx` -- may need recreation for array property
- New `experiment_background_factors_idx`
- Existing `gene_cluster_treatment_type_idx` -- unchanged (already array)
- New index on ClusteringAnalysis/GeneCluster `background_factors` if needed

### Post-import aggregation

Publication level: flatten experiment array properties via `UNWIND` + `collect(DISTINCT)`, same pattern already used for `treatment_types` and `omics_types`. `background_factors` follows the identical aggregation.

OrganismTaxon level: same pattern, aggregated from publications via `apoc.coll.toSet(reduce(...))`.

## Implementation order

1. **Schema + validator + skill** -- define what's valid
2. **Adapters** -- omics_adapter.py, cluster_adapter.py (backward-compatible defaults)
3. **Paperconfig pass 1** -- mechanical: string->list, add obvious background_factors
4. **Paperconfig pass 2** -- manual review of flagged experiments (may need PDF reading)
5. **Post-import + indexes** -- aggregate background_factors, add indexes
6. **Tests** -- unit tests, KG validity tests
7. **Docs** -- CLAUDE.md updates

### Pass 1 vs Pass 2 paperconfig migration

**Pass 1 (mechanical):**
- Convert all `treatment_type` from string to list
- Add `background_factors` where derivable from existing fields:
  - Has `treatment_organism` -> `[coculture]` in background (unless treatment_type is coculture)
  - `light_condition` contains "continuous" -> `[continuous_light]`
  - No coculture indicators -> `[axenic]`
- Flag experiments needing manual review

**Pass 2 (manual, after review):**
- Fill in `background_factors` for flagged experiments
- May require reading paper PDFs to determine culture status or light regime
- Validate all paperconfigs pass

## Scope boundaries

**In scope:**
- Schema, adapter, validator, paperconfig skill, post-import, tests, docs
- All 23 paperconfig files updated

**Out of scope:**
- MCP tool changes (they query treatment_type which keeps working; background_factors filtering can be added later)
- Changing experiment node IDs
- Changing any DE data or expression edges
