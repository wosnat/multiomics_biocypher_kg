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
