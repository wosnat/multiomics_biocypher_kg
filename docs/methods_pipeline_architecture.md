# Methods: Pipeline Architecture and Deployment

## Overview

The knowledge graph is produced by a two-phase pipeline. A shell-orchestrated preprocessing phase (`scripts/prepare_data.sh`) downloads and integrates genomic, proteomic, and publication data from external databases into per-strain annotation tables. A Python build phase (`create_knowledge_graph.py`) consumes those tables through a set of BioCypher-compatible adapters that yield graph nodes and edges as tab-delimited CSV files. A five-stage Docker Compose pipeline then bulk-imports the CSVs into Neo4j, executes post-import graph enrichment queries, and exposes the database alongside a web query interface. The current graph spans 13 bacterial strains (8 *Prochlorococcus*, 2 *Synechococcus*/*Parasynechococcus*, 3 *Alteromonas*), 24 publications, and approximately ~188K differential expression edges (`Condition_changes_expression_of` + `Coculture_changes_expression_of`).

## Data Preprocessing Pipeline

### Pipeline orchestration

The preprocessing pipeline is implemented as a Bash script (`scripts/prepare_data.sh`) that executes five sequential steps (0–4), each invoking a dedicated Python module. Each step's output is consumed by subsequent steps, enforcing a strict dependency chain. Per-step logs are written to `logs/prepare_data_step{N}.log`. The script accepts CLI flags for selective execution: `--steps` restricts which steps run, `--strains` limits processing to specific strains, `--force` re-runs steps even when output files already exist, and `--skip-cyanorak` bypasses downloads from the intermittently unreachable Cyanorak server.

### Step descriptions

| Step | Script | Input | Output | Requires |
|------|--------|-------|--------|----------|
| 0 | `download_genome_data.py` | NCBI, Cyanorak, UniProt APIs | GFF, protein FASTA, GBK, UniProt JSON, `gene_mapping.csv` per strain | Network access |
| 1 | `build_protein_annotations.py` | UniProt JSON per taxid | `protein_annotations.json` per taxid | Step 0 |
| 2 | `build_gene_annotations.py` | `gene_mapping.csv` + eggNOG annotations + `protein_annotations.json` | `gene_annotations_merged.json` per strain | Step 1 |
| 3 | `build_gene_id_mapping.py` | `gene_annotations_merged.json` + paperconfig `id_translation`/`annotation_gff` entries | `gene_id_mapping.json` per strain | Step 2 |
| 4 | `resolve_paper_ids.py` | Paper CSVs + `gene_id_mapping.json` | `*_resolved.csv` per paper table | Step 3 |

Step 0 downloads raw genomic data from three external sources. NCBI Datasets provides GFF annotations, protein FASTA sequences, and GenBank flat files for each assembly accession listed in `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`. The Cyanorak v2.1 database provides supplementary GFF and GBK files for *Prochlorococcus* and *Synechococcus* strains. UniProt proteome records are fetched per NCBI taxonomy ID. The step concludes by building `gene_mapping.csv` for each strain, merging NCBI and Cyanorak gene coordinates and identifiers (see Methods: Position-Based Annotation Merge for the coordinate reconciliation algorithm used for MIT9313).

Step 1 restructures UniProt data into per-taxid protein annotation tables (`protein_annotations.json`) containing accessions, functional descriptions, EC numbers, GO terms, and RefSeq cross-references.

Step 2 merges three annotation sources — `gene_mapping.csv`, eggNOG-mapper output, and `protein_annotations.json` — into a unified `gene_annotations_merged.json` per strain. Merge rules are configured in `config/gene_annotations_config.yaml` and handle field-level priority (e.g., Cyanorak product descriptions override eggNOG for *Prochlorococcus* strains).

Step 3 builds the three-tier gene ID mapping (`gene_id_mapping.json`) that enables resolution of non-standard identifiers (JGI catalog IDs, microarray probesets, alternative locus tags) found in publication CSVs. The algorithm harvests additional ID sources from `id_translation`, `annotation_gff`, and `csv`-level `id_columns` entries declared in paperconfig files (see Methods: Three-Tier Gene Identifier Mapping for the full algorithm).

Step 4 pre-resolves gene identifiers in publication CSVs by applying the mappings from step 3. For each `csv` supplementary table whose `name_col` is not `locus_tag`, the script reads the source CSV, resolves each row's gene identifier to a locus tag, and writes a `*_resolved.csv` file alongside the original with two added columns: `locus_tag` and `resolution_method`. The omics adapter probes for these files during graph construction and uses pre-resolved locus tags instead of the original identifiers (see Methods: Differential Expression Integration for details).

### Shared resources and per-strain processing

UniProt data is downloaded once per NCBI taxonomy ID and shared across strains with the same taxid. The cache directory hierarchy reflects this: strain-specific files reside in `cache/data/<organism_group>/genomes/<strain>/`, while taxid-level resources are stored in `cache/data/<organism_group>/uniprot/<taxid>/`. This avoids redundant API calls for strains of the same species.

## Adapter Pattern and Graph Construction

### Adapter interface

The graph build phase uses the BioCypher adapter pattern. Each data source is encapsulated in a Python class that implements three methods: `download_data(cache=True)` fetches and parses source data with filesystem caching, `get_nodes()` yields `(node_id, label, properties)` tuples, and `get_edges()` yields `(edge_id, source_id, target_id, edge_type, properties)` tuples. Both generators yield one element at a time, avoiding the need to load entire datasets into memory.

### Multi-wrapper aggregation

Each adapter follows a two-level pattern. A single-source class (e.g., `CyanorakNcbi`) handles one strain or data source. A Multi* wrapper class (e.g., `MultiCyanorakNcbi`) reads a configuration file listing all strains or papers, instantiates per-source adapters, and aggregates their node and edge generators. The main entry point (`create_knowledge_graph.py`) instantiates only Multi* wrappers.

### Adapter execution order

Adapters are executed sequentially in a fixed order determined by data dependencies:

| Order | Adapter | Data source | Yields |
|-------|---------|-------------|--------|
| 1 | `MultiCyanorakNcbi` | `gene_annotations_merged.json`, `cyanobacteria_genomes.csv` | Gene, OrganismTaxon, Cyanorak_cluster nodes; Gene_belongs_to_organism, Gene_in_cyanorak_cluster edges |
| 2 | `MultiUniprot` | `protein_annotations.json`, `gene_mapping.csv` | Protein nodes; Gene_encodes_protein, Protein_belongs_to_organism edges |
| 3 | `MultiOMICSAdapter` | paperconfig YAMLs + publication CSVs | Publication, EnvironmentalCondition nodes; Condition_changes_expression_of and Coculture_changes_expression_of edges |
| 4 | `MultiGoAnnotationAdapter` | `gene_annotations_merged.json` + GO OBO | GO term nodes + Gene→GO edges + GO hierarchy |
| 5 | `MultiEcAnnotationAdapter` | `gene_annotations_merged.json` + Expasy | EC number nodes + Gene→EC edges + EC hierarchy |
| 6 | `MultiKeggAnnotationAdapter` | `gene_annotations_merged.json` + KEGG REST API | KEGG KO/Pathway/Subcategory/Category nodes + hierarchy edges |
| 7 | `MultiCogRoleAnnotationAdapter` | `gene_annotations_merged.json` + role definition files | COG category, Cyanorak Role, TIGR Role nodes + Gene→role edges |

The `MultiCyanorakNcbi` adapter must run first: it generates `gene_mapping.csv` files that the `MultiUniprot` adapter uses to join RefSeq protein IDs to gene nodes when creating `Gene_encodes_protein` edges. The `MultiUniprot` adapter must precede the omics adapter because expression edges target gene nodes that must already exist in the output. The four functional annotation adapters (orders 4–7) are independent of each other but depend on gene nodes from adapter 1. See Methods: Functional Annotation as Graph Edges for the schema modeling of these annotation types.

### String sanitization

BioCypher exports node and edge properties as tab-delimited CSV files destined for `neo4j-admin import`. Single-quote characters in string values break the field quoting, and pipe characters corrupt array fields (pipe is the configured array delimiter). Every adapter applies a sanitization function to all string properties before yielding: single quotes are replaced with carets (`^`) and pipes are removed. This is enforced locally in each adapter module rather than centralized, ensuring no unsanitized string reaches the CSV output.

## BioCypher Integration

### Configuration

The pipeline uses BioCypher in offline mode (`biocypher_config.yaml`: `offline: true`), bypassing online ontology validation against the Biolink Model. The graph schema is defined declaratively in `config/schema_config.yaml`, which maps semantic types (e.g., `gene`, `organism taxon`, `publication`) to Neo4j node labels and declares all properties with their data types. The schema currently defines 20 node types and 49 relationship types spanning genomic structure, functional annotation, expression effects, homology, and ontology hierarchies.

### CSV export and import script generation

`create_knowledge_graph.py` passes each adapter's generators to BioCypher's `write_nodes()` and `write_edges()` methods, which serialize tuples into CSV files using a tab field delimiter and a pipe array delimiter. After all adapters have been processed, `write_import_call()` generates a shell script (`neo4j-admin-import-call.sh`) that invokes `neo4j-admin import` with the correct file paths, node labels, and relationship types. The import is configured with `skip_duplicate_nodes` and `skip_bad_relationships` to tolerate edge cases (e.g., expression edges targeting genes absent from the current genome annotation) without aborting the entire import.

### Node ID conventions

Node identifiers follow CURIE conventions normalized via the Bioregistry:

| Node type | Prefix | Example |
|-----------|--------|---------|
| Gene | `ncbigene` | `ncbigene:1253526` |
| Protein | `uniprot` | `uniprot:Q7V6D4` |
| Organism | `insdc.gcf` | `insdc.gcf:GCF_000011465.1` |
| Publication | `doi` | `doi:10.1038/ismej.2014.63` |
| GO term | `go` | `go:GO:0006260` |
| EC number | `eccode` | `eccode:2.7.7.7` |
| KEGG KO | `kegg.orthology` | `kegg.orthology:K02338` |

These prefixed identifiers ensure uniqueness across node types and support interoperability with external databases.

## Docker Deployment Pipeline

### Five-stage service chain

The Docker Compose file (`docker-compose.yml`) defines five services that execute as a sequential pipeline:

| Stage | Image | Purpose | Output |
|-------|-------|---------|--------|
| build | `biocypher/base:1.2.0` | Runs `create_knowledge_graph.py` | CSV files + import script in shared volume |
| import | `neo4j:4.4-enterprise` | Executes generated `neo4j-admin import` script | Neo4j database store |
| post-process | `neo4j:4.4-enterprise` | Runs `scripts/post-import.cypher` | Homology + inferred expression edges |
| deploy | `neo4j:4.4-enterprise` | Serves Neo4j (ports 7474/7687) | Live database for queries and tests |
| app | `biocypher/biochatter-light:0.6.13` | Biochatter web UI (port 8501) | Natural-language query interface |

### Service dependencies and data flow

Each service declares a `depends_on` condition requiring the previous service to complete successfully before starting. Data passes between stages via a shared named Docker volume (`biocypher_neo4j_volume`). The build container copies project source code into its working directory and substitutes the Docker-specific BioCypher configuration (`config/biocypher_docker_config.yaml`), which sets the output path to the shared volume. The import container reads the generated CSV files and import script from this volume, produces the Neo4j database store, and copies the import report for diagnostic review. Authentication is disabled across all Neo4j services to simplify programmatic access by downstream tools and tests.

### Post-import graph enrichment

Certain edge types require graph traversal patterns that cannot be computed during CSV export. The post-process stage executes four Cypher queries from `scripts/post-import.cypher`:

1. **Cyanorak-based homology**: Creates bidirectional `Gene_is_homolog_of_gene` edges between genes sharing a Cyanorak cluster, with a `distance` property capturing phylogenetic proximity (same strain, same clade, same species, same genus, same order, or cross-order).
2. **Alteromonas within-family orthologs**: Creates homology edges between *Alteromonas* genes sharing an Alteromonadaceae-level eggNOG orthologous group.
3. **Cross-phylum orthologs**: Creates homology edges between *Alteromonas* and *Prochlorococcus*/*Synechococcus* genes sharing a Bacteria-level COG orthologous group.
4. **Expression propagation**: For every `Condition_changes_expression_of` or `Coculture_changes_expression_of` edge from a source to gene A, if A is a homolog of gene B, creates a `Condition_changes_expression_of_ortholog` or `Coculture_changes_expression_of_ortholog` edge (respectively) from the same source to gene B, carrying the original expression properties plus the homology provenance.

These queries use `MERGE` for homology edges (idempotent) and `CREATE` for expression propagation (each direct expression edge produces an independent homolog edge). See Methods: Homology Propagation (forthcoming) for the algorithmic rationale.

## Caching Strategy

The pipeline employs filesystem caching at three levels. Raw downloads from external databases (NCBI GFF/FASTA, Cyanorak GFF/GBK, UniProt JSON) are stored in `cache/data/<organism_group>/` and reused across pipeline runs. Preprocessed intermediate files (`gene_mapping.csv`, `gene_annotations_merged.json`, `protein_annotations.json`, `gene_id_mapping.json`) reside alongside raw data in the same cache hierarchy and are regenerated only when upstream inputs change or the `--force` flag is passed. BioCypher CSV output is written to `biocypher-log/` (local) or `/data/build2neo/` (Docker) and regenerated on every build.

Cache bypass mechanisms are provided at each level. The `--force` flag in `prepare_data.sh` re-runs preprocessing steps regardless of existing outputs. The `--no-cache` flag in `create_knowledge_graph.py` forces re-download within adapters. The `--skip-cyanorak` flag works around the intermittently unreachable Cyanorak server by skipping its download substep while reusing any previously cached files. In Docker, the `BUILD2NEO_CLEANUP` environment variable clears previous build output before a fresh run.

## Testing and Iteration

A `--test` flag in `create_knowledge_graph.py` limits each adapter to 100 items, enabling rapid iteration during development without running the full pipeline. Unit tests for individual adapters and utility modules run via `pytest -m "not slow and not kg"`. A full integration test (`pytest -m slow`) builds the entire graph as a subprocess and validates the output. After Docker deployment, a KG validity test suite (`tests/kg_validity/`) connects to the live Neo4j instance and verifies structural invariants, biological constraints, expression edge properties, and snapshot regression data. See Methods: Quality Assurance and Testing (forthcoming) for the full testing methodology.
