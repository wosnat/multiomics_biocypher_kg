# KG Change: treatment_type → array, new background_factors field

**Date:** 2026-04-02

## What changed

### Experiment nodes

| Property | Before | After |
|---|---|---|
| `treatment_type` | `str` (single value, e.g. `"nitrogen_stress"`) | `str[]` (list, e.g. `["nitrogen_stress"]`) |
| `background_factors` | did not exist | `str[]` (list, e.g. `["axenic", "continuous_light"]`). May be null (empty list stored as null by Neo4j) |

### ClusteringAnalysis and GeneCluster nodes

| Property | Before | After |
|---|---|---|
| `treatment_type` | `str[]` (unchanged) | `str[]` (unchanged) |
| `background_factors` | did not exist | `str[]` |

### Publication nodes (post-import computed)

| Property | Before | After |
|---|---|---|
| `treatment_types` | `str[]` (unchanged) | `str[]` (unchanged, but aggregation logic changed for array source) |
| `background_factors` | did not exist | `str[]` (distinct values from experiments) |

### OrganismTaxon nodes (post-import computed)

| Property | Before | After |
|---|---|---|
| `treatment_types` | `str[]` (unchanged) | `str[]` (unchanged) |
| `background_factors` | did not exist | `str[]` (aggregated from publications) |

## Breaking Cypher changes

### Experiment.treatment_type is now an array

Any query comparing `treatment_type` as a scalar string will break or return wrong results.

| Before | After |
|---|---|
| `e.treatment_type = 'nitrogen_stress'` | `'nitrogen_stress' IN e.treatment_type` |
| `e.treatment_type IN ['coculture', 'viral']` | `ANY(t IN e.treatment_type WHERE t IN ['coculture', 'viral'])` |
| `NOT e.treatment_type IN ['coculture', 'viral']` | `NOT 'coculture' IN e.treatment_type AND NOT 'viral' IN e.treatment_type` |
| `collect(DISTINCT e.treatment_type)` | `apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.treatment_type, [])) \| s + t))` |

### New background_factors property

Available for filtering but may be null (when no background factors were specified, e.g. coculture-as-treatment experiments).

```cypher
// Find all experiments involving coculture in any capacity
MATCH (e:Experiment)
WHERE 'coculture' IN e.treatment_type OR 'coculture' IN coalesce(e.background_factors, [])
RETURN e

// Find nitrogen stress experiments, distinguish axenic vs coculture
MATCH (e:Experiment)
WHERE 'nitrogen_stress' IN e.treatment_type
RETURN e.id,
       CASE WHEN 'coculture' IN coalesce(e.background_factors, []) THEN 'coculture'
            WHEN 'axenic' IN coalesce(e.background_factors, []) THEN 'axenic'
            ELSE 'unknown' END AS culture_status

// Find all experiments under diel cycle
MATCH (e:Experiment)
WHERE 'diel_cycle' IN coalesce(e.background_factors, [])
RETURN e
```

## New indexes

- `experiment_background_factors_idx` on `Experiment.background_factors`

## Canonical vocabulary

Shared between `treatment_type` and `background_factors`:

`nitrogen_stress`, `phosphorus_stress`, `iron_stress`, `carbon_stress`, `salt_stress`, `oxygen_stress`, `temperature_stress`, `light_stress`, `darkness`, `plastic_stress`, `viral`, `coculture`, `growth_state`, `growth_medium`, `diel`, `axenic`, `continuous_light`, `diel_cycle`, `chemical_inhibitor`

## Counts

- 76 Experiment nodes, all have `treatment_type` as array (most single-element)
- 68 experiments have non-null `background_factors`
- 8 experiments have null `background_factors` (coculture-as-treatment with no additional context)
- 30 experiments discoverable as coculture-related (16 treatment + 14 background), up from 16
- 7 experiments with `diel_cycle` background (previously undiscoverable)
