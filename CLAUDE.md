# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioCypher-based knowledge graph for multi-omics data on marine cyanobacteria *Prochlorococcus* and *Alteromonas*. Integrates genomic data, proteins, and experimental omics results (RNA-seq, proteomics, metabolomics) from public databases and research publications. The graph feeds into LLM agents for omics data analysis.

## Build and Run Commands

```bash
# Install dependencies (uv recommended)
uv sync

# Build the knowledge graph
uv run python create_knowledge_graph.py

# Run tests
pytest

# Docker deployment (builds graph + Neo4j)
docker compose up -d
# Neo4j at localhost:7474 (HTTP), localhost:7687 (Bolt)
```

## Architecture

### Adapter Pattern
Each data source has a dedicated adapter in `multiomics_kg/adapters/` implementing:
- `get_nodes()` - Generator yielding `(node_id, label, properties)` tuples
- `get_edges()` - Generator yielding `(edge_id, source_id, target_id, edge_type, properties)` tuples
- `download_data(cache=True)` - Fetches data with caching support

### Key Adapters
- **uniprot_adapter.py** - Proteins, genes, organisms from UniProt API
- **cyanorak_ncbi_adapter.py** - Genomic data from GFF/GenBank files
- **omics_adapter.py** - Differential expression results from `paperconfig.yaml` files
- **go_adapter.py** - Gene Ontology annotations
- **ec_adapter.py** - Enzyme Commission number hierarchy

### Main Entry Point
`create_knowledge_graph.py` orchestrates the pipeline:
1. Initializes BioCypher instance
2. Configures adapters with organism ID (default: MED4 = 59919)
3. Calls `download_data()` then `write_nodes()`/`write_edges()` for each adapter
4. Generates Neo4j import scripts via `bc.write_import_call()`

### Configuration
- `config/schema_config.yaml` - Graph schema (node types, edge types, properties)
- `config/biocypher_config.yaml` - BioCypher framework settings

## Adding Omics Data from Publications

Create `paperconfig.yaml` in `data/Prochlorococcus/papers_and_supp/<PaperName>/`:

```yaml
publication:
  papername: "Author Year"
  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "path/to/results.csv"
      statistical_analyses:
        - type: RNASEQ  # or PROTEOMICS, METABOLOMICS, MICROARRAY
          test_type: "DESeq2"
          control_condition: "Axenic"
          treatment_condition: "Coculture"
          experimental_context: "in Pro99 medium under continuous light"
          organism: "Prochlorococcus MED4"
          name_col: "Gene"
          logfc_col: "log2FoldChange"
          adjusted_p_value_col: "padj"
```

Then add to `create_knowledge_graph.py`:
```python
omics = OMICSAdapter(config_file="path/to/paperconfig.yaml")
omics.download_data(cache=CACHE)
bc.write_nodes(omics.get_nodes())
bc.write_edges(omics.get_edges())
```

## Data Locations

- **Genomic data:** `data/Prochlorococcus/genomes/MED4/`
- **Publication data:** `data/Prochlorococcus/papers_and_supp/`

## Notes

- Uses custom fork of pypath-omnipath: `https://github.com/wosnat/pypath`
- Cache intermediate results in JSON to avoid redundant API calls
- Current focus: Prochlorococcus MED4 strain (NCBI taxid 59919)
- API keys for LLM extraction go in `.env` (OPENAI_API_KEY, ANTHROPIC_API_KEY)