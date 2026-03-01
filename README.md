# Multiomics BioCypher Knowledge Graph

BioCypher-based knowledge graph for multi-omics data on marine cyanobacteria *Prochlorococcus* and *Alteromonas*.

## Build

```bash
# Install dependencies (uv required)
uv sync

# Full build
uv run python create_knowledge_graph.py

# Fast iteration (100 items per adapter)
uv run python create_knowledge_graph.py --test

# Include GO and/or EC hierarchy nodes/edges
uv run python create_knowledge_graph.py --go --ec

# Re-fetch all data (bypass cache)
uv run python create_knowledge_graph.py --no-cache

```

## Docker Deployment

Runs the full pipeline in containers: build → import → post-process → Neo4j + Biochatter UI.

```bash
docker compose up -d
```

- Neo4j browser: http://localhost:7474 (no auth)
- Neo4j Bolt: `bolt://localhost:7687`
- Biochatter UI: http://localhost:8501

The pipeline stages are:

| Stage | What it does |
|---|---|
| `build` | Runs `create_knowledge_graph.py`, outputs CSVs |
| `import` | Runs `neo4j-admin import` from generated import script |
| `post-process` | Runs `scripts/post-import.cypher` (creates bidirectional homolog edges) |
| `deploy` | Exposes Neo4j |
| `app` | Biochatter web UI |

### Options

Edit the `build` service environment in [docker-compose.yml](docker-compose.yml) before running:

| Variable | Default | Description |
|---|---|---|
| `BUILD2NEO_CLEANUP` | `yes` | Delete previous build artifacts before rebuilding |
| `BUILD_GO` | `no` | Include GO hierarchy nodes and edges |
| `BUILD_EC` | `no` | Include EC number hierarchy nodes and edges |

Example — enable GO and EC edges:

```yaml
environment:
  BUILD2NEO_CLEANUP: "yes"
  BUILD_GO: "yes"
  BUILD_EC: "yes"
```

## Prepare Genome Data

Run this before `create_knowledge_graph.py` when adding new genomes or refreshing annotation data.

```bash
# Download everything and build annotation tables
bash scripts/prepare_data.sh

# Force re-run all steps
bash scripts/prepare_data.sh --force

# Skip Cyanorak re-download (server can be slow/throttled; cached files are reused)
bash scripts/prepare_data.sh --force --skip-cyanorak

# Specific strains or steps only
bash scripts/prepare_data.sh --strains MED4 MIT9313 --force
bash scripts/prepare_data.sh --steps 1 2 --force
bash scripts/prepare_data.sh --steps 3 --strains MIT9301 --force  # rebuild gene ID mapping only
```

The script runs four steps:

| Step | Script | What it does |
|---|---|---|
| 0 | `download_genome_data.py` | NCBI genomes, Cyanorak GFF/GBK, UniProt, gene_mapping.csv |
| 1 | `build_protein_annotations.py` | Per-taxid protein annotation tables → `protein_annotations.json` |
| 2 | `build_gene_annotations.py` | Merges gene_mapping + eggNOG + UniProt → `gene_annotations_merged.json` |
| 3 | `build_gene_id_mapping.py` | Extended ID mappings from paperconfig id_translation/annotation_gff entries → `gene_id_mapping.json` |

Logs are written to `logs/prepare_data_step{0,1,2,3}.log`. Monitor with `tail -f logs/prepare_data_step0.log`.

## Organisms

# list of organisms mentioned in the papers:
Alteromonas macleodii EZ55
Alteromonas macleodii HOT1A3
Alteromonas macleodii MIT1002
Phage
Prochlorococcus AS9601
Prochlorococcus MED4
Prochlorococcus MIT9301
Prochlorococcus MIT9312
Prochlorococcus MIT9313
Prochlorococcus NATL1A
Prochlorococcus NATL2A
Prochlorococcus RSP50
Synechococcus CC9311
Synechococcus WH8102