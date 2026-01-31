# Graph Schema for Gene Expression Changes

The schema uses a **hybrid approach** optimized for LLM agent queries:

## Direct Edges (for simple LLM queries - 1 hop)

### Coculture Effects
**`organism to gene expression association`**: Represents how coculture with organism X changes expression of genes in organism Y.

Example: Alteromonas → affects_expression_of → Prochlorococcus gene PMM0001

- **Source**: `organism_taxon` (the coculture partner causing the effect)
- **Target**: `gene` (the affected gene, belongs to a different organism)
- **Key properties**: `log2_fold_change`, `adjusted_p_value`, `expression_direction`, `time_point`
- **Context property**: `environmental_treatment_condition_id` - references environmental condition node (what was held constant)

### Environmental Stress Effects
**`environmental condition to gene expression association`**: Represents how environmental perturbations affect gene expression.

Example: nitrogen_stress → affects_expression_of → gene PMM0001 (same label as organism edge)

- **Source**: `environmental condition`
- **Target**: `gene`
- **Key properties**: same as organism edge plus `control_condition`, `treatment_condition`
- **Context property**: `biological_context` - what biological conditions were held constant (e.g., "axenic", "coculture_with_alteromonas")

## Context Properties (Single-Factor Experiment Design)

Most experiments follow single-factor design: change one variable while holding others constant. The schema captures this with context properties:

| Edge Type | Varied Factor | Context Property | Description |
|-----------|---------------|------------------|-------------|
| `organism to gene expression` | Coculture partner | `environmental_treatment_condition_id` | References environmental condition node (what was held constant) |
| `environmental condition to gene expression` | Environmental stress | `biological_context` | String describing biological context (e.g., "axenic", "coculture_with_alteromonas") |

This design allows LLM agents to:
1. Query direct effects: "Genes affected by Alteromonas" → use organism edge
2. Filter by context: Join to environmental condition node via `environmental_treatment_condition_id`
3. Compare contexts: "Do the same genes respond to Alteromonas under different conditions?" → compare edges with different environmental_treatment_condition_id values

## Environmental Condition Node

The `environmental condition` node represents various environmental perturbations:

| Condition Type | Relevant Properties |
|----------------|---------------------|
| `gas_shock` | `oxygen_level`, `co2_level` (e.g., "0%", "depleted", "21%") |
| `nutrient_stress` | `nitrogen_source`, `nitrogen_level`, `phosphate_level` (e.g., "replete", "limited", "starved") |
| `light_stress` | `light_condition`, `light_intensity` (e.g., "darkness", "blue", "red", "low", "high") |

## Hub Model (for multi-factorial experiments)

For experiments combining multiple factors (e.g., coculture + low light), use the `statistical_test` node as a hub:

```
Alteromonas ──organism_treatment_in_test──▶ statistical_test ◀──environmental_condition_in_test── low_light
                                                  │
                                                  │ molecular_result_from_test
                                                  ▼
                                           Prochlorococcus gene
```

### Hub Model Edges

- **`organism treatment in test`**: Links organism_taxon → statistical_test (with `role` property: "coculture_partner", "predator", "symbiont", "host")
- **`environmental condition in test`**: Links environmental condition → statistical_test (with `role` property: "treatment", "control", "baseline")
- **`molecular result from test`**: Links gene/protein → statistical_test (with expression results)

## Query Patterns for LLM Agents

| Query | Recommended Pattern |
|-------|---------------------|
| "What affects gene PMM0001?" | `(factor)-[affects_expression_of]->(gene {id: 'PMM0001'})` - returns both organisms and environmental conditions |
| "Genes upregulated by Alteromonas" | `(Alteromonas)-[affects_expression_of {expression_direction: 'up'}]->(gene)` |
| "Genes affected by nitrogen stress" | `(nitrogen_stress:environmental_condition)-[affects_expression_of]->(gene)` |
| "Genes affected by coculture under low light" | Hub: Find statistical_test linked to both factors, then get gene results |
| "All experiments involving Alteromonas" | Hub: `(Alteromonas)-[organism_treatment_in_test]->(test)` |
