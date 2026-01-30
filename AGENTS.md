# Goal of this file

This file contain instructions for AI agents for creating and structuring the code of this repo

# Goal of the repo
This repo contains code for building a knowledge graph that describes the known genomic information on the cyanobacteria Prochlorococcus and marine bacteria Alteromonas. and the interaction between them as measured in coculture experiments. These include information on genes, proteins and metabolites from public databases (uniprot, kegg, ec, cyanorak, and more). As well as results of RNASEQ, affymetrix arrays, proteomics and metablomics from this study and other publications.

The resulting knowledge graph will be fed to LLM agents and used to analyze fresh omics data.
We are using biocypher infrastructure which can create neo4j graphs.

# directory structure of the repo:
* create_knowledge_graph.py - the main script that build the graph.
* config/schema_config.yaml - the graph schema
* multiomics_kg/adapters - adapters used to build the graph
* data/Prochlorococcus/papers_and_supp - existing publications. each subfolder contains the PDF of a publication as well as relevant supplemental informantion files and tables.
* data/Prochlorococcus/genomes/MED4 - genomic information on Prochlorococcus MED4 - the first organism to insert into the graph.

# implementation strategy
Start with a graph for a single Prochlorococcus strain - Prochlorococcus MED4. Reuse and tailor adapters from existing biocypher implementations where possible. 
Leverage unit tests whereever possible to check continuous correction as more data is added.


Use a combination of LLM based entity extraction and summarization, automatic downloads, and manual coding to convert the files to the required format. cache intermediate results in json file to avoid LLM costs and time consuming reruns.


# for the omics result adapter

## Overview
The omics adapter reads YAML configuration files from each publication folder and creates statistical test nodes and molecular result edges in the knowledge graph. This allows integration of differential expression analyses, proteomics results, metabolomics data, and other omics studies.

## Configuration File Format

Use a `paperconfig.yaml` file per publication, located in the paper's folder under `data/Prochlorococcus/papers_and_supp/<paper_name>/`

### YAML Structure

```yaml
publication:
  papername: <Paper Name and Year>
  papermainpdf: <path to PDF>

  # Define environmental conditions used in this paper (creates nodes)
  environmental_conditions:
    <condition_id>:  # unique ID for this condition
      condition_type: <type>  # growth_medium, nutrient_stress, light_stress, gas_shock
      name: <descriptive name>
      light_condition: <e.g., continuous_light, darkness>
      light_intensity: <e.g., 10 µmol photons m-2 s-1>
      temperature: <e.g., 24C>
      nitrogen_level: <e.g., replete, starved>
      # ... other condition-specific properties

  supplementary_materials:
    supp_table_<N>:
      type: csv  # or other file type
      filename: <path to data file>
      statistical_analyses:  # List of one or more analyses
        - type: RNASEQ  # or PROTEOMICS, METABOLOMICS, MICROARRAY, etc.
          name: <Descriptive test name>
          test_type: <Algorithm name, e.g., Rockhopper, DESeq2, limma>
          control_condition: <Baseline condition>
          treatment_condition: <Experimental condition>
          timepoint: <Time point, e.g., "20h", "day 7">
          reference_timepoint: <Baseline time if applicable>
          organism: <Organism name or strain>
          # Coculture experiment fields (for organism → gene edges)
          treatment_organism: <Name of coculture partner organism>
          treatment_taxid: <NCBI taxid of coculture partner>
          environmental_treatment_condition_id: <condition_id>  # references node defined above
          # Column mappings
          name_col: <Name of column with gene/protein IDs>
          logfc_col: <Name of column with log2 fold change>
          adjusted_p_value_col: <Name of column with adjusted p-value>
        - type: PROTEOMICS  # Additional analysis for same table
          # ... more fields ...
```

### Multiple Analyses Per Table

A single supplementary table can contain results from multiple statistical analyses. Specify them as a list under `statistical_analyses`:

```yaml
supp_table_1:
  type: csv
  filename: data/results.csv
  statistical_analyses:
    - type: RNASEQ
      name: "RNA-seq analysis"
      test_type: DESeq2
      # ... other fields ...
    - type: PROTEOMICS  
      name: "Proteomics analysis"
      test_type: limma
      # ... other fields ...
```

### Example Configuration

See [data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml](data/Prochlorococcus/papers_and_supp/Aharonovich%202016/paperconfig.yaml) for a complete working example.

## Creating a Configuration for a New Paper

1. **Identify supplementary tables**: Find which supplementary tables contain omics results (gene expression, protein levels, metabolite abundance, etc.)

2. **Create paperconfig.yaml** in the paper's folder with:
   - Publication metadata (name, PDF location)
   - One entry per supplementary table containing omics data
   - For each table:
     - File path to the CSV/Excel file with results
     - Column names in the data file (name_col, logfc_col, p_value_col, etc.)
     - Metadata about the statistical analysis:
       - Test type (e.g., RNA-seq method used)
       - Conditions being compared (control vs treatment)
       - Time point (if time series)
       - Organism studied
     - **For coculture experiments** (creates organism → gene edges):
       - `treatment_organism`: Name of coculture partner (e.g., "Alteromonas macleodii HOT1A3")
       - `treatment_taxid`: NCBI taxid of coculture partner (e.g., 28108)
       - `environmental_treatment_condition_id`: References an environmental condition node defined in the `environmental_conditions` section

3. **Data file requirements**:
   - CSV format recommended for compatibility
   - Must have columns for gene/protein identifiers and statistical results
   - Must have fold change and adjusted p-value columns

## OMICSAdapter Implementation

### Class: OMICSAdapter

Located in `multiomics_kg/adapters/omics_adapter.py`

#### Key Methods

- `load_config(config_file)`: Load and parse a paperconfig.yaml file
- `get_nodes()`: Generate statistical test nodes from the config
- `get_edges()`: Generate molecular result edges connecting genes/proteins to their test results

#### Generated Node Types

**statistical_test**: Represents a statistical analysis
- Properties:
  - `name`: Descriptive name of the test
  - `test_type`: Algorithm/method used (e.g., Rockhopper, DESeq2)
  - `control_condition`: Baseline condition
  - `treatment_condition`: Experimental condition
  - `timepoint`: When measurement was taken
  - `organism`: Organism studied

#### Generated Edge Types

**molecular_result_from_test**: Connects a gene/protein to a statistical test with results
- Source: gene or protein node
- Target: statistical_test node
- Properties:
  - `log2_fold_change`: Log2 of fold change between conditions
  - `adjusted_p_value`: Adjusted p-value from statistical test
  - `direction`: "up" for upregulation, "down" for downregulation

### Usage Example

```python
from multiomics_kg.adapters.omics_adapter import OMICSAdapter

# Initialize adapter with a config file
adapter = OMICSAdapter(
    config_file="data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml"
)

# Get nodes and edges
test_nodes = adapter.get_nodes()
result_edges = adapter.get_edges()

# Add to BioCypher graph
for node in test_nodes:
    bc.add_node(*node)
for edge in result_edges:
    bc.add_edge(*edge)
```

## Integration into create_knowledge_graph.py

To integrate omics data from a publication:

1. Create the paperconfig.yaml file in the publication folder
2. Add to create_knowledge_graph.py:
   ```python
   from multiomics_kg.adapters.omics_adapter import OMICSAdapter
   
   omics = OMICSAdapter(
       config_file="data/Prochlorococcus/papers_and_supp/<PaperName>/paperconfig.yaml"
   )
   
   # Add to graph
   for node in omics.get_nodes():
       bc.add_node(*node)
   for edge in omics.get_edges():
       bc.add_edge(*edge)
   ```

## Supported Omics Types

The adapter is flexible and can handle various omics data types. The `type` field in statistical_analysis specifies the data type:

- **RNASEQ**: Transcript abundance from RNA sequencing
- **MICROARRAY**: Gene expression from microarray platforms
- **PROTEOMICS**: Protein abundance from mass spectrometry or similar
- **METABOLOMICS**: Metabolite abundance measurements
- **PHOSPHOPROTEOMICS**: Phosphorylation site abundance
- Custom types can be added as needed

## Schema Integration

The adapter creates nodes and edges that conform to the BioCypher schema defined in `config/schema_config.yaml`:
- `statistical_test` node type with corresponding properties
- `molecular_result_from_test` edge type connecting to genes/proteins

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

# Agent Skills: Cypher Query Examples

This section provides ready-to-use Cypher query templates for LLM agents querying the Prochlorococcus-Alteromonas knowledge graph.

## Skill: Find What Affects a Gene

Find all factors (organisms or environmental conditions) that affect expression of a specific gene:

```cypher
// All factors affecting gene expression
MATCH (factor)-[r:affects_expression_of]->(g:gene {locus_tag: 'PMM0001'})
RETURN factor, r.expression_direction, r.log2_fold_change, r.adjusted_p_value
ORDER BY abs(r.log2_fold_change) DESC
```

## Skill: Find Genes Affected by Organism Coculture

```cypher
// Genes upregulated by Alteromonas coculture
MATCH (org:organism {organism_name: 'Alteromonas macleodii'})-[r:affects_expression_of {expression_direction: 'up'}]->(g:gene)
RETURN g.locus_tag, g.product, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.log2_fold_change DESC
LIMIT 20

// Genes downregulated by Alteromonas coculture
MATCH (org:organism)-[r:affects_expression_of {expression_direction: 'down'}]->(g:gene)
WHERE org.organism_name CONTAINS 'Alteromonas'
RETURN g.locus_tag, g.product, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.log2_fold_change ASC
LIMIT 20
```

## Skill: Find Genes Affected by Environmental Stress

```cypher
// Genes affected by nitrogen stress
MATCH (env:environmental_condition {condition_type: 'nutrient_stress'})-[r:affects_expression_of]->(g:gene)
WHERE env.nitrogen_level = 'starved'
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change
ORDER BY abs(r.log2_fold_change) DESC

// Genes affected by light/dark shift
MATCH (env:environmental_condition)-[r:affects_expression_of]->(g:gene)
WHERE env.light_condition = 'darkness'
RETURN g.locus_tag, r.expression_direction, r.log2_fold_change
```

## Skill: Compare Conditions

```cypher
// Genes affected by BOTH Alteromonas AND nitrogen stress
MATCH (org:organism)-[r1:affects_expression_of]->(g:gene)
WHERE org.organism_name CONTAINS 'Alteromonas'
MATCH (env:environmental_condition)-[r2:affects_expression_of]->(g)
WHERE env.nitrogen_level = 'starved'
RETURN g.locus_tag, g.product,
       r1.expression_direction AS coculture_direction,
       r2.expression_direction AS stress_direction,
       r1.log2_fold_change AS coculture_fc,
       r2.log2_fold_change AS stress_fc

// Genes with OPPOSITE response to coculture vs stress
MATCH (org:organism)-[r1:affects_expression_of]->(g:gene)
MATCH (env:environmental_condition)-[r2:affects_expression_of]->(g)
WHERE r1.expression_direction <> r2.expression_direction
RETURN g.locus_tag, g.product, r1.expression_direction, r2.expression_direction
```

## Skill: Filter by Context (Environmental Condition)

```cypher
// Genes affected by Alteromonas under specific environmental conditions
MATCH (org:organism)-[r:affects_expression_of]->(g:gene)
WHERE org.organism_name CONTAINS 'Alteromonas'
MATCH (env:environmental_condition {id: r.environmental_treatment_condition_id})
WHERE env.light_condition = 'continuous_light'
RETURN g.locus_tag, r.log2_fold_change, env.name

// Compare coculture response under different environmental conditions
MATCH (org:organism)-[r1:affects_expression_of]->(g:gene)
MATCH (env1:environmental_condition {id: r1.environmental_treatment_condition_id})
WHERE env1.light_condition = 'continuous_light'
MATCH (org)-[r2:affects_expression_of]->(g)
MATCH (env2:environmental_condition {id: r2.environmental_treatment_condition_id})
WHERE env2.light_condition = 'darkness'
RETURN g.locus_tag,
       r1.log2_fold_change AS light_fc,
       r2.log2_fold_change AS dark_fc,
       env1.name AS light_condition,
       env2.name AS dark_condition
ORDER BY abs(r1.log2_fold_change - r2.log2_fold_change) DESC
```

## Skill: Get Gene Details with Function

```cypher
// Gene with its protein and GO annotations
MATCH (g:gene {locus_tag: 'PMM0001'})
OPTIONAL MATCH (g)-[:Gene_encodes_protein]->(p:protein)
OPTIONAL MATCH (p)-[:protein_involved_in_biological_process]->(bp:biological_process)
OPTIONAL MATCH (p)-[:protein_enables_molecular_function]->(mf:molecular_function)
RETURN g.locus_tag, g.product, p.primary_protein_name,
       collect(DISTINCT bp.name) AS biological_processes,
       collect(DISTINCT mf.name) AS molecular_functions
```

## Skill: Find Functionally Related Affected Genes

```cypher
// Genes in same biological process that are affected by Alteromonas
MATCH (org:organism)-[r:affects_expression_of]->(g:gene)
WHERE org.organism_name CONTAINS 'Alteromonas'
MATCH (g)-[:Gene_encodes_protein]->(p:protein)-[:protein_involved_in_biological_process]->(bp:biological_process)
RETURN bp.name AS process,
       collect(g.locus_tag) AS affected_genes,
       avg(r.log2_fold_change) AS avg_fold_change,
       count(g) AS gene_count
ORDER BY gene_count DESC
LIMIT 10
```

## Skill: Multi-factorial Experiments (Hub Model)

```cypher
// Find experiments with both coculture and environmental stress
MATCH (org:organism)-[:organism_treatment_in_test]->(test:statistical_test)
MATCH (env:environmental_condition)-[:environmental_condition_in_test]->(test)
RETURN test.name, org.organism_name, env.name, test.timepoint

// Get gene results from multi-factorial experiment
MATCH (org:organism)-[:organism_treatment_in_test]->(test:statistical_test)
MATCH (env:environmental_condition)-[:environmental_condition_in_test]->(test)
MATCH (g:gene)-[r:molecular_result_from_test]->(test)
WHERE org.organism_name CONTAINS 'Alteromonas'
RETURN g.locus_tag, r.log2_fold_change, r.direction,
       org.organism_name, env.name
ORDER BY abs(r.log2_fold_change) DESC
```

## Skill: Significance Filtering

```cypher
// Significantly differentially expressed genes (p < 0.05, |FC| > 2)
MATCH (factor)-[r:affects_expression_of]->(g:gene)
WHERE r.adjusted_p_value < 0.05 AND abs(r.log2_fold_change) > 1
RETURN factor, g.locus_tag, r.expression_direction, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.adjusted_p_value ASC

// Top 10 most significantly affected genes
MATCH (org:organism)-[r:affects_expression_of]->(g:gene)
WHERE org.organism_name CONTAINS 'Alteromonas'
RETURN g.locus_tag, g.product, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.adjusted_p_value ASC
LIMIT 10
```

## Skill: Time Series Analysis

```cypher
// Gene expression across time points
MATCH (org:organism)-[r:affects_expression_of]->(g:gene {locus_tag: 'PMM0001'})
WHERE org.organism_name CONTAINS 'Alteromonas'
RETURN r.time_point, r.log2_fold_change, r.expression_direction
ORDER BY r.time_point

// Genes with consistent response across time
MATCH (org:organism)-[r:affects_expression_of]->(g:gene)
WITH g, collect(r.expression_direction) AS directions
WHERE size([d IN directions WHERE d = 'up']) = size(directions)
   OR size([d IN directions WHERE d = 'down']) = size(directions)
RETURN g.locus_tag, g.product, directions
```

## Skill: Publication Traceability

```cypher
// Find which publication reported the expression change
MATCH (factor)-[r:affects_expression_of]->(g:gene)
WHERE g.locus_tag = 'PMM0001'
RETURN factor, r.log2_fold_change, r.publications

// All genes from a specific study
MATCH (factor)-[r:affects_expression_of]->(g:gene)
WHERE 'doi:10.1234/example' IN r.publications
RETURN g.locus_tag, factor, r.expression_direction
```
