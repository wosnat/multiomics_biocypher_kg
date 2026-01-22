# OMICS Adapter Quick Start Guide

## Overview

The `OMICSAdapter` reads YAML configuration files (`paperconfig.yaml`) from publication folders and creates:
- **Statistical test nodes**: Metadata about differential expression/proteomics/metabolomics analyses
- **Molecular result edges**: Connections between genes/proteins and their test results with fold changes and p-values

## Quick Start

### 1. Create a Config File

Create `data/Prochlorococcus/papers_and_supp/<YourPaper>/paperconfig.yaml`:

```yaml
publication:
  papername: "Your Paper Name Year"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/YourPaper/paper.pdf"
  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "data/Prochlorococcus/papers_and_supp/YourPaper/results.csv"
      statistical_analyses:  # List allows multiple tests per table
        - type: RNASEQ
          name: "Your analysis description"
          test_type: "DESeq2"  # or Rockhopper, limma, etc.
          control_condition: "condition A"
          treatment_condition: "condition B"
          timepoint: "24h"
          organism: "Prochlorococcus MED4"
          name_col: "Gene"
          logfc_col: "log2FoldChange"
          adjusted_p_value_col: "padj"
        - type: PROTEOMICS
          name: "Proteomics results"
          test_type: "limma"
          control_condition: "condition A"
          treatment_condition: "condition B"
          timepoint: "24h"
          organism: "Prochlorococcus MED4"
          name_col: "ProteinID"
          logfc_col: "logFC"
          adjusted_p_value_col: "adj.P.Val"
```

### 2. Use in create_knowledge_graph.py

```python
from multiomics_kg.adapters.omics_adapter import OMICSAdapter

# Load adapter
omics = OMICSAdapter(
    config_file="data/Prochlorococcus/papers_and_supp/YourPaper/paperconfig.yaml"
)

# Add to graph
for node in omics.get_nodes():
    bc.add_node(*node)
for edge in omics.get_edges():
    bc.add_edge(*edge)
```

## Configuration Details

### Statistical Analysis Fields

Required:
- `name`: Human-readable test name
- `test_type`: Algorithm used (e.g., DESeq2, Rockhopper)
- `control_condition`: Baseline condition
- `treatment_condition`: Experimental condition
- `organism`: Organism name/strain

Column Mappings:
- `name_col`: Column containing gene/protein identifiers
- `logfc_col`: Column with log2 fold changes
- `adjusted_p_value_col`: Column with adjusted p-values

Optional:
- `type`: Omics data type (RNASEQ, PROTEOMICS, METABOLOMICS, etc.)
- `timepoint`: When measurement was taken
- `reference_timepoint`: Baseline timepoint if applicable

### Multiple Analyses Per Table

You can specify multiple statistical analyses for the same supplementary table using `statistical_analyses` (plural) as a list:

```yaml
supp_table_1:
  type: csv
  filename: "data/results.csv"
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

### Data File Requirements

- **Format**: CSV with headers
- **Gene ID column**: Must contain unique identifiers
- **Fold change column**: Numeric values (positive = upregulated)
- **P-value column**: Numeric values (0.0 to 1.0)

## Generated Nodes and Edges

### Nodes: `statistical_test`

Represents a statistical analysis. Properties include:
- `name`: Test name
- `test_type`: Algorithm/method
- `control_condition`: Baseline
- `treatment_condition`: Experimental
- `timepoint`: Time point
- `organism`: Organism studied

### Edges: `molecular_result_from_test`

Connects genes/proteins to test results:
- **Source**: Gene/protein node (e.g., `ncbigene:PMM0735`)
- **Target**: Statistical test node
- **Properties**:
  - `log2_fold_change`: Fold change value
  - `adjusted_p_value`: P-value
  - `direction`: "up" or "down"

## Example: Aharonovich 2016

See the complete working example:
- Config: `data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml`
- Data files: CSV files in same directory
- Results: 3 statistical tests, 6,251 molecular result edges

## Supported Omics Types

- `RNASEQ` - RNA sequencing
- `MICROARRAY` - Gene expression arrays
- `PROTEOMICS` - Protein abundance
- `METABOLOMICS` - Metabolite levels
- `PHOSPHOPROTEOMICS` - Phosphorylation sites
- Custom types as needed

## Troubleshooting

### YAML Parsing Errors
- Check indentation (use spaces, not tabs)
- Ensure colons are followed by spaces
- Validate YAML at: https://www.yamllint.com/

### Missing Data Files
- Verify file paths are correct
- Paths are relative to workspace root
- Check file extensions match type (csv, xlsx, etc.)

### No Edges Generated
- Verify column names in CSV match config
- Check that identifier column has non-empty values
- Ensure fold change and p-value columns have numeric data

## Integration with BioCypher Schema

The adapter creates nodes/edges conforming to:
- Node type: `statistical_test` (defined in `config/schema_config.yaml`)
- Edge type: `molecular_result_from_test`
- Automatically prefixes gene IDs with `ncbigene:` prefix

For more details, see [agents.md](agents.md#for-the-omics-result-adapter)
