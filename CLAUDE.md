# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioCypher-based knowledge graph for multi-omics data on marine cyanobacteria *Prochlorococcus* and *Alteromonas*. Integrates genomic data, proteins, and experimental omics results (RNA-seq, proteomics, metabolomics) from public databases and research publications. The graph feeds into LLM agents for omics data analysis.

## Build and Run Commands

```bash
# Install dependencies (uv required)
uv sync

# Build the knowledge graph
uv run python create_knowledge_graph.py

# Run adapter unit tests (excluding slow integration tests and KG tests)
pytest -m "not slow and not kg"

# Run a single test file
pytest tests/test_omics_adapter_organism_gene.py -v

# Run the full integration test (builds entire graph, ~1 hour)
pytest -m slow

# Run KG validity tests (requires running Docker graph at localhost:7687)
pytest tests/kg_validity/ -v
pytest -m kg                          # same thing via marker
pytest tests/kg_validity/ --neo4j-url bolt://localhost:7687  # explicit URL

# Docker deployment (builds graph + Neo4j)
docker compose up -d
# Neo4j at localhost:7474 (HTTP), localhost:7687 (Bolt, no auth)
# Biochatter UI at localhost:8501
```

## Development Flags in create_knowledge_graph.py

Two flags at the top of `create_knowledge_graph.py` control pipeline behavior:
- `TEST_MODE = True` — stops each adapter after 100 items; use for fast iteration
- `CACHE = True` — reuses cached JSON files from previous runs; set to `False` to re-fetch data

## Architecture

### Adapter Pattern
Each data source has a dedicated adapter in `multiomics_kg/adapters/` implementing:
- `get_nodes()` — Generator yielding `(node_id, label, properties)` tuples
- `get_edges()` — Generator yielding `(edge_id, source_id, target_id, edge_type, properties)` tuples
- `download_data(cache=True)` — Fetches/parses data with caching

Most adapters have a single-source class (e.g., `CyanorakNcbi`) and a **Multi*** wrapper (e.g., `MultiCyanorakNcbi`) that iterates over all configured strains/papers and aggregates nodes/edges. `create_knowledge_graph.py` only instantiates the Multi* wrappers.

### Key Adapters
- **uniprot_adapter.py** — Proteins, genes, organisms from UniProt API
- **cyanorak_ncbi_adapter.py** — Genomic data from GFF/GenBank files; also emits `organism taxon` nodes
- **omics_adapter.py** — Differential expression results from `paperconfig.yaml` files; creates `publication` nodes and `affects_expression_of` edges
- **go_adapter.py** — Gene Ontology terms and GO→protein edges
- **ec_adapter.py** — Enzyme Commission number hierarchy
- **pathway_adapter.py** — Metabolic/signaling pathway associations
- **pdf_publication_extraction.py** — LangChain-based LLM extraction of publication metadata from PDFs

### Main Entry Point
`create_knowledge_graph.py` orchestrates the pipeline:
1. Initializes BioCypher instance with `config/biocypher_config.yaml`
2. Reads three registry files:
   - `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` — genome → NCBI accession/taxid mapping
   - `data/Prochlorococcus/treatment_organisms.csv` — non-genomic organisms used as edge sources (e.g., phage, coculture partners)
   - `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` — list of paths to all `paperconfig.yaml` files
3. Calls `download_data()` then `write_nodes()`/`write_edges()` for each adapter
4. Generates Neo4j import scripts via `bc.write_import_call()`

### Configuration
- `config/schema_config.yaml` — Graph schema (node types, edge types, properties, ID conventions)
- `config/biocypher_config.yaml` — BioCypher framework settings (offline mode, delimiter, array delimiter `|`)
- `config/biocypher_docker_config.yaml` — Same but with Docker-specific paths; `scripts/build.sh` copies it over the main config

### Node ID Conventions
- Genes: NCBI GeneID (integer as string)
- Proteins: UniProt accession
- Organisms: INSDC GCF assembly accession (e.g., `GCF_000011465.1`)
- Publications: DOI
- GO terms: GO:XXXXXXX
- EC numbers: EC hierarchy string

### Docker Pipeline Stages
1. **build** — runs `create_knowledge_graph.py`, outputs CSVs
2. **import** — runs `neo4j-admin import` from generated script
3. **post-process** — executes `scripts/post-import.cypher`, which creates bidirectional `gene_is_homolog_of_gene` edges between genes sharing a Cyanorak cluster (needed because BioCypher can't infer these during import)
4. **deploy** — exposes Neo4j
5. **app** — Biochatter web UI

## Adding Omics Data from Publications

1. Create directory: `data/Prochlorococcus/papers_and_supp/<Author Year>/`
2. Create `paperconfig.yaml` there (see `/skill:paperconfig` for interactive creation):

```yaml
publication:
  papername: "Author Year"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Author Year/paper.pdf"
  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "data/.../results.csv"
      statistical_analyses:
        - id: "unique_id_within_publication"   # required; used as edge ID prefix
          type: RNASEQ  # or PROTEOMICS, METABOLOMICS, MICROARRAY
          test_type: "DESeq2"
          control_condition: "Axenic"
          treatment_condition: "Coculture"
          experimental_context: "in Pro99 medium under continuous light"
          organism: "Prochlorococcus MED4"
          name_col: "Gene"
          logfc_col: "log2FoldChange"
          adjusted_p_value_col: "padj"
```

3. Add the path to `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt`
4. Verify gene ID matching: use the `/check-gene-ids` skill; if mapping is needed, use `/fix-gene-ids`
5. The `MultiOMICSAdapter` picks it up automatically from the text file list

## Custom Claude Code Skills

The `.claude/skills/` directory provides project-specific skills:
- `/paperconfig` — interactive wizard to create paperconfig.yaml
- `/check-gene-ids` — validate gene ID match rates between CSV data and graph nodes
- `/fix-gene-ids` — map gene IDs to locus tags when mismatches are found
- `/cypher-queries` — run Cypher queries against Neo4j with ready-made templates

## Data Locations

- **Genomic data (raw):** `data/Prochlorococcus/genomes/<Strain>/` and equivalents for Synechococcus/Alteromonas
- **Genomic data (cached):** `cache/data/<Organism>/genomes/<Strain>/`
- **Publication data:** `data/Prochlorococcus/papers_and_supp/<Author Year>/`
- **Gene mapping CSVs:** `cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv`
- **UniProt cache:** `uniprot_raw_data.json`, `uniprot_preprocess_data.json` (root dir, not committed)
- **PDF extraction cache:** `pdf_extraction_cache.json`

## KG Validity Tests

`tests/kg_validity/` contains a fast (~14s) pytest suite that connects to the live deployed Neo4j instance and validates the output graph. Tests auto-skip if Neo4j is unreachable. All tests are marked `@pytest.mark.kg`.

### Test files

| File | What it validates |
|---|---|
| `test_structure.py` | Node type presence, minimum counts (>5K genes, ≥12 organisms), orphan detection (genes without organism, proteins without organism), Gene_encodes_protein edge validation, key property presence |
| `test_biology.py` | Ecotype/clade labels per strain, katG absence in *Prochlorococcus* (Black Queen Hypothesis), all expected strains present, locus-tag → UniProt spot checks |
| `test_expression.py` | `log2_fold_change` / `adjusted_p_value` are numeric, `adjusted_p_value` ∈ [0,1], `expression_direction` ∈ {up,down}, direction/sign consistency, required properties on all edges |
| `test_post_import.py` | Homolog edges exist and are bidirectional, `distance`/`cluster_id` properties present, `Affects_expression_of_homolog` propagated correctly |

### Actual Neo4j labels (BioCypher PascalCase output)

- Nodes: `Gene`, `Protein`, `OrganismTaxon`, `Publication`, `EnvironmentalCondition`, `Cyanorak_cluster`
- Relationships: `Gene_belongs_to_organism`, `Protein_belongs_to_organism`, `Gene_encodes_protein`, `Gene_in_cyanorak_cluster`, `Gene_is_homolog_of_gene`, `Affects_expression_of`, `Affects_expression_of_homolog`

### Key graph facts

- Gene↔Protein linkage: explicit `Gene_encodes_protein` edges (Protein→Gene) created by UniProt adapter via RefSeq protein_id join with gene_mapping.csv
- Expression sources: `EnvironmentalCondition` (~95K edges, stress experiments), `OrganismTaxon` (~12K edges, coculture experiments)
- `adjusted_p_value` may be null on expression edges (and propagated homolog edges) when the original study did not report it
- Strains in graph: MED4, AS9601, MIT9301, MIT9312, MIT9313, NATL1A, NATL2A, RSP50 (Prochlorococcus); CC9311 (Synechococcus); WH8102 (Parasynechococcus); MIT1002, EZ55, HOT1A3 (Alteromonas)

## Notes

- Uses custom fork of pypath-omnipath: `https://github.com/wosnat/pypath`
- Current organisms: 8 *Prochlorococcus* strains + 2 *Synechococcus/Parasynechococcus* + 3 *Alteromonas*; primary focus is MED4 (NCBI taxid 59919)
- API keys go in `.env`: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `ELSEVIER_API_KEY`, plus others for LangChain tracing and search
- `@pytest.mark.slow` tests run the full pipeline as a subprocess and take ~1 hour
- `@pytest.mark.kg` tests require a running Neo4j instance (Docker graph); skip automatically if not available
