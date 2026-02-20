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

---

## Where Functional Information Lives

The graph has two node types that may seem redundant but have distinct roles:

- **`Gene`** = the genomic entity: has location, locus tag, expression edges, homology edges
- **`Protein`** = the functional/sequence entity: has UniProt annotations, GO terms, pathways, EC numbers

Every Gene (for all organisms) is linked to its Protein via `(Protein)-[:Gene_encodes_protein]->(Gene)`.

### Functional annotation coverage by organism

| Property | Prochlorococcus | Synechococcus | Alteromonas |
|---|---|---|---|
| `Gene.product` | ✓ (NCBI) | ✓ (NCBI) | ✓ (NCBI) |
| `Gene.function_description` | ✓ (copied from UniProt at import) | ✓ | ✓ |
| `Gene.go_biological_processes` | ✓ (Cyanorak + UniProt) | ✓ (UniProt) | ✓ (UniProt) |
| `Gene.cyanorak_Role`, `Gene.eggNOG` | ✓ (Cyanorak only) | ✗ | ✗ |
| `Protein.function_description` | ✓ (UniProt) | ✓ | ✓ |
| `BiologicalProcess` nodes (GO) | ✓ via Protein | ✓ via Protein | ✓ via Protein |

> **Note**: `Gene.function_description` and `Gene.go_biological_processes` are denormalized from the linked Protein node during post-import. They are the recommended starting point for function-based queries.

### Decision tree: "which genes are related to X?"

**If X is a biological function / process** (e.g., "oxidative stress", "photosynthesis"):

```cypher
// Option A — text search on Gene (fastest, works for all organisms):
MATCH (g:Gene)
WHERE g.function_description CONTAINS 'oxidative'
   OR any(t IN g.go_biological_processes WHERE t CONTAINS 'oxidative')
RETURN g

// Option B — structured GO term traversal (precise, requires knowing GO term name):
MATCH (g:Gene)<-[:Gene_encodes_protein]-(p:Protein)
      -[:protein_involved_in_biological_process]->(bp:BiologicalProcess)
WHERE bp.name CONTAINS 'oxidative stress'
RETURN g, bp

// Option C — KEGG / EC / pathway (for metabolic functions):
MATCH (g:Gene)<-[:Gene_encodes_protein]-(p:Protein)
      -[:protein_take_part_in_pathway]->(path:Pathway)
WHERE path.name CONTAINS 'oxidative'
RETURN g, path
```

**If X is an experimental treatment** (e.g., "nitrogen stress", "coculture with Alteromonas"):

```cypher
// Direct edge — 1 hop, most specific:
MATCH (c:EnvironmentalCondition)-[r:Affects_expression_of]->(g:Gene)
WHERE c.name CONTAINS 'nitrogen'
RETURN g, r.log2_fold_change, r.expression_direction

MATCH (o:OrganismTaxon)-[r:Affects_expression_of]->(g:Gene)
WHERE o.organism_name CONTAINS 'Alteromonas'
RETURN g, r.log2_fold_change, r.expression_direction
```

**If you want genes with BOTH a functional annotation AND expression evidence**:

```cypher
MATCH (g:Gene)
WHERE g.function_description CONTAINS 'oxidative'
AND (g)<-[:Affects_expression_of]-()
RETURN g
```

### What NOT to do

- Do **not** query `Gene.pathway_description` — pathways are on `Protein` nodes and linked via `protein_take_part_in_pathway` edges.
- Do **not** use `Gene.ec_numbers` — EC numbers are on `Protein` nodes and linked via `protein_catalyzes_ec_number` edges to the `EcNumber` hierarchy.
- For Prochlorococcus, `Gene.cyanorak_Role` and `Gene.eggNOG` provide additional functional text but are absent for Alteromonas/Synechococcus — prefer `Gene.function_description` for organism-agnostic queries.
