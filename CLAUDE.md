# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BioCypher-based knowledge graph for multi-omics data on marine cyanobacteria *Prochlorococcus* and *Alteromonas*. Integrates genomic data, proteins, and experimental omics results (RNA-seq, proteomics, metabolomics) from public databases and research publications. The graph feeds into LLM agents for omics data analysis.

## Build and Run Commands

```bash
# Install dependencies (uv required)
uv sync

# Build the knowledge graph (full run)
uv run python create_knowledge_graph.py

# Fast iteration: test mode (100 items per adapter)
uv run python create_knowledge_graph.py --test

# Include GO and/or EC hierarchy nodes/edges
uv run python create_knowledge_graph.py --go --ec

# Re-fetch all data (bypass cache)
uv run python create_knowledge_graph.py --no-cache

# Custom output directory
uv run python create_knowledge_graph.py --output-dir ./my-output/

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

## CLI Flags for create_knowledge_graph.py

`create_knowledge_graph.py` accepts the following command-line options:

| Flag | Default | Description |
|---|---|---|
| `--test` | off | Test mode: stop each adapter after 100 items |
| `--go` | off | Download and write GO nodes/edges |
| `--ec` | off | Download and write EC nodes/edges |
| `--no-cache` | off | Re-fetch all data instead of reusing cached files |
| `--output-dir PATH` | `./biocypher-log/example_knowledge_graph/` | Output directory for CSV exports |

## Architecture

### Adapter Pattern
Each data source has a dedicated adapter in `multiomics_kg/adapters/` implementing:
- `get_nodes()` — Generator yielding `(node_id, label, properties)` tuples
- `get_edges()` — Generator yielding `(edge_id, source_id, target_id, edge_type, properties)` tuples
- `download_data(cache=True)` — Fetches/parses data with caching

#### String property sanitization
BioCypher wraps string fields in single quotes in the CSV output. Any string property value containing a single quote (`'`) or a pipe (`|`) will break `neo4j-admin import`. **All string properties must be sanitized before yielding from `get_nodes()` or `get_edges()`:**
```python
def _clean_str(value: str) -> str:
    return value.replace("'", "^").replace("|", "")
```
Define this helper locally in each adapter and apply it to every string field in node/edge property dicts. The array delimiter is also `|` (see `biocypher_config.yaml`), so pipe characters in values would corrupt array fields too.

Most adapters have a single-source class (e.g., `CyanorakNcbi`) and a **Multi*** wrapper (e.g., `MultiCyanorakNcbi`) that iterates over all configured strains/papers and aggregates nodes/edges. `create_knowledge_graph.py` only instantiates the Multi* wrappers.

### Key Adapters
- **uniprot_adapter.py** — Proteins, genes, organisms from UniProt API
- **cyanorak_ncbi_adapter.py** — Genomic data from GFF/GenBank files; also emits `organism taxon` nodes
- **omics_adapter.py** — Differential expression results from `paperconfig.yaml` files; creates `publication` nodes and expression edges (`Condition_changes_expression_of` / `Coculture_changes_expression_of`)
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
4. If the publication uses non-standard gene IDs (JGI catalog IDs, probesets, etc.), add `id_translation` and/or `annotation_gff` entries to the paperconfig (see "Gene ID Mapping" below), then run `bash scripts/prepare_data.sh --steps 3 4 --strains <Strain> --force`
5. Verify gene ID matching: use the `/check-gene-ids` skill; if mapping is needed, use `/fix-gene-ids`
6. The `MultiOMICSAdapter` picks it up automatically from the text file list

### paperconfig supplementary_materials entry types

| Type | Purpose |
|---|---|
| `csv` | Differential expression results table with `statistical_analyses` |
| `id_translation` | Pure ID mapping table (no DE data); bridges non-standard IDs to locus tags |
| `annotation_gff` | GFF3/GTF file; adds protein_id/Name bridges for a strain |

**`id_translation` example** — maps JGI catalog IDs via an annotated genome CSV:
```yaml
annotation_genome_9301:
  type: id_translation
  filename: "data/.../9301_annotated_genome.csv"
  organism: "Prochlorococcus MIT9301"
  id_columns:
    - column: "ID"
      id_type: jgi_id
    - column: "uniprot_gene_name"   # "dnaA P9301_05911" → whitespace-split anchor
      id_type: gene_name
    - column: "uniprot_entry_name"  # "DNAA_PROM0" → strip _ORGANISM suffix
      id_type: uniprot_entry_name
  product_columns:
    - column: "uniprot_protein_name"
```

**`id_translation` with `generate` block** — auto-generated via diamond protein matching for cross-assembly bridging:
```yaml
id_translation_ez55_author:
  type: id_translation
  filename: "data/.../aez55_to_ez55_id_translation.csv"
  organism: "Alteromonas EZ55"
  generate:
    method: diamond_protein_match
    source_fasta: "data/.../ez55_aa.fasta"
    source_id_col: aez55_id
    img_gff: "data/.../EZ55.gff"   # optional, for header remapping
  id_columns:
    - column: "locus_tag"
      id_type: locus_tag
    - column: "aez55_id"
      id_type: old_locus_tag
```

The `generate` block causes `build_gene_id_mapping.py` (step 3) to automatically run `scripts/map_img_to_ncbi_proteins.py` when the output file is missing or `--force` is given. Fields: `method` (`diamond_protein_match`), `source_fasta` (path to source protein FASTA), `source_id_col` (output column name), `img_gff` (optional GFF for header remapping). The `--ncbi-faa` and `--gene-mapping` args are derived from the organism's genome_dir.

**`annotation_gff` example** — extracts locus_tag/old_locus_tag/protein_id bridges:
```yaml
annotation_gff_mit9301:
  type: annotation_gff
  filename: "data/.../GCF_000015965.1_ASM1596v1_genomic.gff"
  organism: "Prochlorococcus MIT9301"
```

**`csv` with non-standard ID columns:**
```yaml
supp_table_1:
  type: csv
  organism: "Prochlorococcus MIT9301"
  id_columns:
    - column: "ID"
      id_type: jgi_id
  statistical_analyses:
    - id: "unique_id"
      ...
```

### Gene ID Mapping

`gene_id_mapping.json` (schema v2, per strain in `cache/data/<Organism>/genomes/<Strain>/`) stores a three-tier ID mapping built by `build_gene_id_mapping.py` (prepare_data step 3) from:
1. `gene_annotations_merged.json` (locus tags, old locus tags, UniProt accessions, protein IDs)
2. `id_translation` paperconfig entries (JGI IDs, probesets, etc.)
3. `annotation_gff` paperconfig entries (additional GFF bridges)
4. `csv` table `id_columns` from paperconfigs (extra identifier columns per paper)

**Three-tier ID classification:**

| Tier | ID types | Storage | Resolution |
|------|----------|---------|------------|
| 1 (gene-unique) | `locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `old_locus_tag`, `alternative_locus_tag`, `jgi_id`, `probeset_id`, `uniprot_entry_name` | `specific_lookup` (1:1) | Always used; conflicts = data errors |
| 2 (protein-level) | `protein_id_refseq`, `protein_id`, `uniprot_accession` | `multi_lookup` (1:many) | Used only when singleton for this organism (paralogs expected) |
| 3 (generic) | `gene_name`, `gene_synonym`, `gene_oln`, `em_preferred_name` | `multi_lookup` (1:many) | Used only when singleton |

**Iterative convergence**: All sources are processed in passes until no new ID-to-locus_tag mapping is added. Transitive closure is automatic (order-independent). Typically 2–3 passes.

**v2 schema** (top-level keys):
- `specific_lookup`: `{alt_id: locus_tag}` — fast 1:1 Tier 1 lookup
- `multi_lookup`: `{alt_id: [locus_tag, ...]}` — Tier 2+3, singletons usable for resolution
- `conflicts`: `{tier1_id: [lt1, lt2]}` — Tier 1 data errors for manual review
- `genes`: `{locus_tag: {tier1_ids, tier2_ids, tier3_ids}}` — per-gene ID catalogue

`gene_id_utils.py` exposes the v2 API:
- `load_mapping_v2(genome_dir)` → `MappingData`
- `resolve_row(row, name_col, id_columns, mapping_data)` → `(locus_tag | None, method_str)`
- `expand_list(raw_val)` → splits list-valued cells on `,` and `;`

**Resolution method strings**: `tier1:<col>`, `heuristic:<col>` (zero-pad / strip asterisk), `multi:<col>` (Tier 2+3 singleton), `tier1_conflict`, `ambiguous`, `unresolved`.

Diagnostic report after build: `gene_id_mapping_report.json` (per-ID-type stats, reclassification warnings).

See `docs/methods_gene_id_mapping.md` for the full scientific description.

## Custom Claude Code Skills

The `.claude/skills/` directory provides project-specific skills:
- `/paperconfig` — interactive wizard to create paperconfig.yaml
- `/check-gene-ids` — validate gene ID match rates between CSV data and graph nodes
- `/fix-gene-ids` — map gene IDs to locus tags when mismatches are found
- `/cypher-queries` — run Cypher queries against Neo4j with ready-made templates
- `/deploy-strain` — end-to-end checklist for deploying v2 gene ID mapping for a new strain (snapshot → rebuild mapping → resolve papers → verify → rebuild KG → compare)

## Genome Data Download Pipeline

`scripts/prepare_data.sh` orchestrates all genome annotation downloads and preprocessing. Run this before `create_knowledge_graph.py` when adding new genomes or refreshing data.

```bash
# Download everything + build annotation tables (skip eggNOG — run /eggnog-run skill separately)
bash scripts/prepare_data.sh

# Force re-run all steps
bash scripts/prepare_data.sh --force

# Cyanorak server throttles sometimes — skip it when files are already cached
bash scripts/prepare_data.sh --force --skip-cyanorak

# Specific strains or steps only
bash scripts/prepare_data.sh --strains MED4 MIT9313 --force
bash scripts/prepare_data.sh --steps 1 2 --force   # rebuild annotation tables only
```

Logs written to `logs/prepare_data_step0.log` … `logs/prepare_data_step4.log`. Monitor with `tail -f logs/prepare_data_step0.log`.

**Step 0** (`multiomics_kg/download/download_genome_data.py`) — sub-steps:
- 1: NCBI genome (GFF, protein FASTA, GBFF)
- 2: Cyanorak GFF/GBK (strains with `cyanorak_organism` only; server can be slow)
- 3: UniProt per unique taxid → `cache/data/<org_group>/uniprot/<taxid>/`
- 4: eggNOG-mapper (skipped by default; requires `EGGNOG_DATA_DIR` in `.env`)
- 5: Build `gene_mapping.csv`

**Step 1** (`multiomics_kg/download/build_protein_annotations.py`) — builds per-taxid protein annotation tables → `protein_annotations.json` per taxid. Requires step 0 (UniProt data must be cached first).

**Step 2** (`multiomics_kg/download/build_gene_annotations.py`) — merges gene_mapping.csv + eggNOG + UniProt (protein_annotations.json) → `gene_annotations_merged.json` per strain. Requires step 1.

**Step 3** (`multiomics_kg/download/build_gene_id_mapping.py`) — builds extended gene ID mappings (`gene_id_mapping.json` and backward-compat `gene_mapping_supp.csv`) per strain by harvesting `id_translation` and `annotation_gff` entries from all paperconfigs. Run as module: `uv run python -m multiomics_kg.download.build_gene_id_mapping`. Required for strains with non-standard IDs (JGI catalog IDs, probesets, etc.) in publication CSVs. Requires step 2.

**Step 4** (`multiomics_kg/download/resolve_paper_ids.py`) — pre-resolves gene IDs in paper CSVs: for each `csv` supplementary table whose `name_col` is not already `locus_tag`, loads the source CSV, resolves each `name_col` value via `build_id_lookup()` (uses `gene_id_mapping.json`), and writes `<stem>_resolved.csv` alongside the original with two extra columns: `locus_tag` (NaN when unresolved) and `resolution_method`. The `omics_adapter` probes for these files and uses the pre-resolved locus tags instead of the original `name_col`, avoiding dangling edges. Run as module: `uv run python -m multiomics_kg.download.resolve_paper_ids [--force] [--papers "Paper Name"]`. Requires step 3.

## Data Locations

- **Genomic data (raw):** `data/Prochlorococcus/genomes/<Strain>/` and equivalents for Synechococcus/Alteromonas
- **Genomic data (cached):** `cache/data/<Organism>/genomes/<Strain>/`
- **Publication data:** `data/Prochlorococcus/papers_and_supp/<Author Year>/`
- **Gene mapping CSVs:** `cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv`
- **Gene annotation tables:** `cache/data/<Organism>/genomes/<Strain>/gene_annotations_merged.json`
- **Extended gene ID mapping:** `cache/data/<Organism>/genomes/<Strain>/gene_id_mapping.json` (built by step 3; includes paper-derived IDs)
- **Backward-compat ID mapping:** `cache/data/<Organism>/genomes/<Strain>/gene_mapping_supp.csv` (generated alongside gene_id_mapping.json)
- **UniProt cache (taxid-keyed):** `cache/data/<org_group>/uniprot/<taxid>/uniprot_preprocess_data.json`
- **PDF extraction cache:** `pdf_extraction_cache.json`

## KG Validity Tests

`tests/kg_validity/` contains a fast (~14s) pytest suite that connects to the live deployed Neo4j instance and validates the output graph. Tests auto-skip if Neo4j is unreachable. All tests are marked `@pytest.mark.kg`.

### Test files

| File | What it validates |
|---|---|
| `test_structure.py` | Node type presence, minimum counts (>5K genes, ≥12 organisms), orphan detection (genes without organism, proteins without organism), Gene_encodes_protein edge validation, key property presence |
| `test_biology.py` | Ecotype/clade labels per strain, katG absence in *Prochlorococcus* (Black Queen Hypothesis), all expected strains present, locus-tag → UniProt spot checks |
| `test_expression.py` | `log2_fold_change` / `adjusted_p_value` are numeric, `adjusted_p_value` ∈ [0,1], `expression_direction` ∈ {up,down}, direction/sign consistency, required properties on all edges |
| `test_post_import.py` | Homolog edges exist and are bidirectional, `distance`/`cluster_id` properties present, `Condition_changes_expression_of_ortholog` and `Coculture_changes_expression_of_ortholog` propagated correctly |
| `test_snapshot.py` | Regression snapshot: verifies a sample of specific nodes and edges (with properties) still exist after rebuilds. Catches silent data loss. |

### Snapshot regression tests

`test_snapshot.py` loads `snapshot_data.json` (committed fixture) and checks that every sampled node/edge still exists with expected properties. To regenerate after intentional changes:

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

### Actual Neo4j labels (BioCypher PascalCase output)

- Nodes: `Gene`, `Protein`, `OrganismTaxon`, `Publication`, `EnvironmentalCondition`, `Cyanorak_cluster`, `BiologicalProcess`, `CellularComponent`, `MolecularFunction`
- Relationships: `Gene_belongs_to_organism`, `Protein_belongs_to_organism`, `Gene_encodes_protein`, `Gene_in_cyanorak_cluster`, `Gene_is_homolog_of_gene`, `Condition_changes_expression_of`, `Coculture_changes_expression_of`, `Condition_changes_expression_of_ortholog`, `Coculture_changes_expression_of_ortholog`, `Published_expression_data_about`, `Gene_involved_in_biological_process`, `Gene_located_in_cellular_component`, `Gene_enables_molecular_function`, `Biological_process_is_a_biological_process`, `Cellular_component_is_a_cellular_component`, `Molecular_function_is_a_molecular_function`

### Key graph facts

- Gene↔Protein linkage: explicit `Gene_encodes_protein` edges (Protein→Gene) created by UniProt adapter via RefSeq protein_id join with gene_mapping.csv
- Expression sources: `EnvironmentalCondition` → `Condition_changes_expression_of` (~170K edges, stress experiments); `OrganismTaxon` → `Coculture_changes_expression_of` (~17K edges, coculture experiments); total direct edges ~188K (TODO: update after build)
- Ortholog edges: `Condition_changes_expression_of_ortholog` (all homologs); `Coculture_changes_expression_of_ortholog` (same-phylum homologs only; cross-phylum filtered)
- New edge properties on all expression edge types: `omics_type` (RNASEQ | PROTEOMICS | METABOLOMICS | MICROARRAY), `organism_strain`, `treatment_condition`, `statistical_test`
- `condition_category` property on `EnvironmentalCondition` nodes (same value as `condition_type`)
- `preferred_name` property on `OrganismTaxon` nodes (e.g., `"Prochlorococcus MED4"`)
- `Published_expression_data_about` edges: `(Publication)→(EnvironmentalCondition)` and `(Publication)→(OrganismTaxon)` — one edge per distinct source node used in a publication's analyses
- `adjusted_p_value` may be null on expression edges (and propagated ortholog edges) when the original study did not report it
- Strains in graph: MED4, AS9601, MIT9301, MIT9312, MIT9313, NATL1A, NATL2A, RSP50 (Prochlorococcus); CC9311 (Synechococcus); WH8102 (Parasynechococcus); MIT1002, EZ55, HOT1A3 (Alteromonas)

## EggNOG Mapper Setup

EggNOG mapper is installed as a project dependency (`eggnog-mapper>=2.1.13` in `pyproject.toml`).
Set `EGGNOG_DATA_DIR` in `.env` to point to the database directory (e.g. `~/tools/eggnog-mapper`).

**Known bug (v2.1.13):** `download_eggnog_data.py` uses the wrong base URL (`eggnogdb.embl.de` → 404).
After `uv sync`, patch `.venv/bin/download_eggnog_data.py` lines 15 and 18:
```python
# Wrong (old domain, returns 404):
BASE_URL = f'http://eggnogdb.embl.de/download/emapperdb-{__DB_VERSION__}'
NOVEL_FAMS_BASE_URL = f'http://eggnogdb.embl.de/download/novel_fams-{__NOVEL_FAMS_DB_VERSION__}'
# Correct:
BASE_URL = f'http://eggnog5.embl.de/download/emapperdb-{__DB_VERSION__}'
NOVEL_FAMS_BASE_URL = f'http://eggnog5.embl.de/download/novel_fams-{__NOVEL_FAMS_DB_VERSION__}'
```
See: https://github.com/eggnogdb/eggnog-mapper/issues/575

## Known Issues

- **~46% of UniProt proteins are orphaned (no Protein_belongs_to_organism or Gene_encodes_protein edges)** — The UniProt adapter only creates these edges when a protein's RefSeq WP_ ID matches an entry in `gene_mapping.csv`. UniProt returns proteins for a taxid regardless of whether they have a WP_ cross-reference, so proteins absent from our NCBI gene_mapping are stranded. Investigation needed: check whether this fraction is a pre-existing data gap or a regression from the Feb 2026 UniProt adapter refactor (`fe5c2bb`). See `plans/orphan_proteins.md` for details. The KG validity tests `test_no_orphan_proteins` and `test_no_orphan_proteins_without_gene` are currently failing because of this.

- **Cyanorak web server is intermittently unavailable** — `bioinformatics.psb.ugent.be` returns connection errors or throttles requests without warning. This is a server-side issue outside our control. When it happens, use `--skip-cyanorak` (files already in cache are reused). The `TestCyanorakDownloadFile` tests in `tests/test_download_genome_data.py` may also be affected if the server is down during CI.

## Notes

- Uses custom fork of pypath-omnipath: `https://github.com/wosnat/pypath`
- Current organisms: 8 *Prochlorococcus* strains + 2 *Synechococcus/Parasynechococcus* + 3 *Alteromonas*; primary focus is MED4 (NCBI taxid 59919)
- API keys go in `.env`: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `ELSEVIER_API_KEY`, plus others for LangChain tracing and search
- `@pytest.mark.slow` tests run the full pipeline as a subprocess and take ~1 hour
- `@pytest.mark.kg` tests require a running Neo4j instance (Docker graph); skip automatically if not available
