# Phase 4: Edge Type Split + New Properties

**Status:** IN PROGRESS
**Date:** 2026-03-11
**Pre-condition snapshot:** `pre_phase4` — 188,501 total edges (170,904 EnvironmentalCondition + 17,597 OrganismTaxon)

---

## Overview

Split `Affects_expression_of` into 4 typed edges, add new properties, add Publication connectivity edges, and add `condition_category` to EnvironmentalCondition nodes.

---

## Task Checklist

### B — Schema + Adapters

#### B1: schema_config.yaml

**File:** `config/schema_config.yaml`

Replace the two `affects_expression_of` edge definitions (lines 31–63) and two `affects_expression_of_homolog` definitions (lines 71–111) with 4 new edge types:

```yaml
# NEW: organism (coculture) → gene
coculture to gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: coculture_changes_expression_of
  label_in_input: coculture_changes_expression_of
  source: organism_taxon
  target: gene
  properties:
    expression_direction: str
    control_condition: str
    experimental_context: str
    time_point: str
    log2_fold_change: float
    adjusted_p_value: float
    significant: boolean
    publications: str[]
    omics_type: str          # RNASEQ | PROTEOMICS | METABOLOMICS | MICROARRAY
    organism_strain: str     # which organism's genes are measured (e.g., "MED4")
    treatment_condition: str # what was applied (e.g., "Coculture with Alteromonas HOT1A3")
    statistical_test: str    # DESeq2, edgeR, etc.

# NEW: environmental condition → gene
condition to gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: condition_changes_expression_of
  label_in_input: condition_changes_expression_of
  source: environmental condition
  target: gene
  properties:
    expression_direction: str
    control_condition: str
    experimental_context: str
    time_point: str
    log2_fold_change: float
    adjusted_p_value: float
    significant: boolean
    publications: str[]
    omics_type: str
    organism_strain: str
    treatment_condition: str
    statistical_test: str

# NEW: organism (coculture) → gene (ortholog, inferred)
coculture to ortholog gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: coculture_changes_expression_of_ortholog
  label_in_input: coculture_changes_expression_of_ortholog
  source: organism_taxon
  target: gene
  properties:
    expression_direction: str
    control_condition: str
    experimental_context: str
    time_point: str
    log2_fold_change: float
    adjusted_p_value: float
    significant: boolean
    publications: str[]
    omics_type: str
    organism_strain: str
    treatment_condition: str
    statistical_test: str
    original_gene: str
    homology_source: str
    homology_cluster_id: str
    distance: str

# NEW: environmental condition → gene (ortholog, inferred)
condition to ortholog gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: condition_changes_expression_of_ortholog
  label_in_input: condition_changes_expression_of_ortholog
  source: environmental condition
  target: gene
  properties:
    expression_direction: str
    control_condition: str
    experimental_context: str
    time_point: str
    log2_fold_change: float
    adjusted_p_value: float
    significant: boolean
    publications: str[]
    omics_type: str
    organism_strain: str
    treatment_condition: str
    statistical_test: str
    original_gene: str
    homology_source: str
    homology_cluster_id: str
    distance: str
```

Also add to `environmental condition` node (after `condition_type`):
```yaml
    condition_category: str  # same value as condition_type (canonical)
```

Also add `Published_expression_data_about` edge (new, after the 4 expression edges):
```yaml
publication to expression source association:
  is_a: Association
  represented_as: edge
  label_as_edge: published_expression_data_about
  label_in_input: published_expression_data_about
  source: publication
  target: [environmental condition, organism taxon]
```

#### B2: omics_adapter.py

**File:** `multiomics_kg/adapters/omics_adapter.py`

1. **OMICSEdgeType enum** — replace `affects_expression_of` with:
   ```python
   condition_changes_expression_of = auto()
   coculture_changes_expression_of = auto()
   published_expression_data_about = auto()
   ```

2. **`set_edge_types`** — update default list to include new types.

3. **`get_nodes`** — in the environmental condition node creation loop, add `condition_category`:
   ```python
   # After the loop over env_data items:
   condition_type = env_data.get('condition_type')
   if condition_type:
       env_properties['condition_category'] = condition_type
   ```

4. **`get_edges`** — add collection of Publication→source edges (new). After processing all supp materials, emit one `published_expression_data_about` edge per distinct source node seen during analysis processing.

5. **`_load_and_create_edges`** —
   - Determine edge label based on source type:
     - `env_condition_id` → `condition_changes_expression_of`
     - organism → `coculture_changes_expression_of`
   - Add new properties to `edge_properties`:
     ```python
     omics_type = analysis.get('type', '')  # RNASEQ, PROTEOMICS, etc.
     organism_strain = analysis.get('organism', '').split()[-1] if analysis.get('organism') else ''
     treatment_condition = analysis.get('treatment_condition', '')
     statistical_test = analysis.get('test_type', '')
     ```
   - Return `edge_label` alongside edges so `get_edges` can collect unique source nodes for Publication edges.

#### B3: Organism strain extraction

For `organism_strain`, extract the last word of the organism name (e.g., `"Prochlorococcus MED4"` → `"MED4"`, `"Alteromonas macleodii HOT1A3"` → `"HOT1A3"`). This matches the strain field used in queries.

---

### C — Tests

**File:** `tests/test_omics_adapter_organism_gene.py`

Add tests for:
1. Config with `environmental_treatment_condition_id` → edge label `condition_changes_expression_of`
2. Config with `treatment_organism` → edge label `coculture_changes_expression_of`
3. New properties (`omics_type`, `organism_strain`, `treatment_condition`, `statistical_test`) present on edges
4. `condition_category` on EnvironmentalCondition nodes = `condition_type`
5. `published_expression_data_about` edges emitted from Publication → source nodes

---

### D — Skills

1. **`omics-edge-snapshot/SKILL.md`** — update to count `Condition_changes_expression_of` + `Coculture_changes_expression_of` instead of `Affects_expression_of`
2. **`.claude/skills/omics-edge-snapshot/omics_edge_snapshot.py`** — update Cypher queries
3. **`cypher-queries/SKILL.md`** — replace `Affects_expression_of` templates with new edge types
4. **`paperconfig/SKILL.md`** — document edge routing

---

### E — Docs

1. **`CLAUDE.md`** — Update "Actual Neo4j labels" section:
   - Replace `Affects_expression_of` with `Condition_changes_expression_of`, `Coculture_changes_expression_of`
   - Replace `Affects_expression_of_homolog` with `Condition_changes_expression_of_ortholog`, `Coculture_changes_expression_of_ortholog`
   - Add `Published_expression_data_about`

---

## Expected Counts After Rebuild

| Edge type | Expected count |
|---|---|
| `Condition_changes_expression_of` | ~170,904 (was EnvironmentalCondition total) |
| `Coculture_changes_expression_of` | ~17,597 (was OrganismTaxon/NamedThing total) |
| Total direct | ~188,501 (same as pre_phase4) |

---

## Gate Criteria

- [x] `pytest tests/test_omics_adapter_organism_gene.py -v` — all new edge label tests pass
- [x] `uv run python create_knowledge_graph.py --test` — no build errors
- [x] `pytest -m "not slow and not kg"` — all unit tests pass
- [x] No `affects_expression_of` string literals in adapter or schema
- [x] `condition_category` set on all EnvironmentalCondition nodes
- [x] `published_expression_data_about` edges emitted
- [x] All new properties populated on edges

---

## Gate Result

**PASSED 2026-03-11** — 188,501 → 188,501 edges (zero data loss). Split: 170,904 `Condition_changes_expression_of` + 17,597 `Coculture_changes_expression_of`. All 21 publications intact, zero per-gene regressions.
