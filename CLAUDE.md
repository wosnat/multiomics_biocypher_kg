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

# Docker deployment (builds graph + Neo4j, APOC plugin auto-installed)
docker compose up -d
# Neo4j at localhost:7474 (HTTP), localhost:7687 (Bolt, no auth)
# APOC core available in post-process and deploy containers
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

**This applies to computed/literal strings too, not just external data.** Never use `|` or `'` as a separator or literal character when constructing string property values in code. For example, use ` :: ` as a field separator in computed summary strings, not ` | `. `clean_text()` silently converts `|` → `,` and `'` → `^`, so using them as intentional separators would corrupt the value without any error.

Most adapters have a single-source class (e.g., `CyanorakNcbi`) and a **Multi*** wrapper (e.g., `MultiCyanorakNcbi`) that iterates over all configured strains/papers and aggregates nodes/edges. `create_knowledge_graph.py` only instantiates the Multi* wrappers.

### Key Adapters
- **uniprot_adapter.py** — Proteins, genes, organisms from UniProt API
- **cyanorak_ncbi_adapter.py** — Genomic data from GFF/GenBank files; also emits `organism taxon` nodes
- **omics_adapter.py** — Differential expression results from `paperconfig.yaml` files; creates `Publication` nodes, `Experiment` nodes, `Has_experiment` edges, `Tests_coculture_with` edges, and `Changes_expression_of` expression edges (Experiment → Gene)
- **go_adapter.py** — Gene Ontology terms and GO→protein edges
- **ec_adapter.py** — Enzyme Commission number hierarchy
- **pathway_adapter.py** — Metabolic/signaling pathway associations
- **functional_annotation_adapter.py** — Gene→GO edges (per-strain from `gene_annotations_merged.json`), Gene→Pfam edges + Pfam/PfamClan nodes (from `pfam_ids` in merged annotations + Pfam reference data)
- **ortholog_group_adapter.py** — OrthologGroup nodes and Gene_in_ortholog_group edges from pre-computed `ortholog_groups` field in `gene_annotations_merged.json`
- **metabolism_adapter.py** — Reaction + Metabolite nodes + Gene_catalyzes_reaction / Reaction_has_metabolite / Reaction_in_kegg_pathway / Metabolite_in_pathway edges (reads `cache/data/kegg/kegg_data.json` built by prepare-data step 6)
- **tcdb_adapter.py** — TCDB transport ontology + substrate→Metabolite bridge. Reads `cache/data/tcdb/{tcdb_hierarchy,tcdb_pruned}.json` built by step 6. Pruned to the gene-reachable subhierarchy (above + below)
- **cazy_adapter.py** — CAZy carbohydrate-active-enzyme ontology (observed-only; CAZY_CLASSES table from `multiomics_kg/utils/cazy_utils.py`)
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
- Experiments: `{doi}_{experiment_key}` (e.g., `10.1038/s41586-024-07246-9_n_limitation_med4_rnaseq`)
- GO terms: GO:XXXXXXX
- OrthologGroups: `ortholog_group:<source>:<name>` (e.g., `ortholog_group:cyanorak:CK_00000364`)
- Pfam domains: `pfam:PF00712`
- Pfam clans: `pfam.clan:CL0060`
- EC numbers: EC hierarchy string

### Docker Pipeline Stages
1. **build** — runs `create_knowledge_graph.py`, outputs CSVs
2. **import** — runs `neo4j-admin import` from generated script. Container runs as root with `/bin/bash` entrypoint (overriding the neo4j image's entrypoint that would drop to the `neo4j` user) so `scripts/import.sh` can chmod `/output` and write `import.status` + `import.report` artifacts the KG validity suite reads. Re-running `import` on an existing graph requires stopping `deploy` and `app` first (the deploy container holds a lock on the database files in `biocypher_neo4j_volume`)
3. **post-process** — executes `scripts/post-import.sh` (the authoritative post-import script run by Docker via `cypher-shell`). Creates indexes (scalar, full-text, OrthologGroup, Pfam, Experiment), derives `expression_status` on expression edges, and computes Gene routing signal properties (`annotation_types`, `expression_edge_count`, `significant_up_count`, `significant_down_count`, `reaction_count`, `metabolite_count`, `closest_ortholog_group_size`, `closest_ortholog_genera`) plus `rank_by_effect`, `rank_up`, and `rank_down` on expression edges. `scripts/post-import.cypher` is a reference copy kept in sync for non-Docker use -- both files must have identical Cypher logic. The script groups ~30 statements into three `cypher-shell` invocations (indexes / small aggregations / `CALL IN TRANSACTIONS` writes) to avoid per-call JVM startup overhead; each group is wrapped in `time` and prints a `[timing]` line to the container log. When modifying the cypher logic, run `scripts/post-import-validate.sh > baseline.txt` against the currently deployed graph first, then rebuild and `diff baseline.txt after.txt` — the dump is deterministic and must be byte-identical for a pure refactor. `pytest tests/kg_validity/ -v` (or `-m kg`) covers the same computed properties via explicit assertions.
4. **deploy** — exposes Neo4j
5. **app** — Biochatter web UI

## Adding Omics Data from Publications

1. Create directory: `data/Prochlorococcus/papers_and_supp/<Author Year>/`
2. Create `paperconfig.yaml` there (see `/skill:paperconfig` for interactive creation):

```yaml
publication:
  papername: "Author Year"
  doi: "10.NNNN/xxxxx"  # optional — overrides PDF-extracted DOI
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Author Year/paper.pdf"

  experiments:
    my_experiment:
      name: "MED4 nitrogen limitation RNA-seq"
      organism: "Prochlorococcus MED4"
      treatment_condition: "Nitrogen limitation"
      control_condition: "Replete medium"
      experimental_context: "in Pro99 medium under continuous light"
      omics_type: RNASEQ  # or PROTEOMICS, METABOLOMICS, MICROARRAY
      test_type: "DESeq2"
      treatment_type: ["nitrogen"]   # category for filtering -- list of categories
      background_factors: ["axenic", "continuous_light"]  # optional: experimental context factors
      medium: "Pro99"
      temperature: "24C"
      light_condition: "continuous light"
      light_intensity: ""
      # For coculture experiments, add:
      # treatment_organism: "Alteromonas macleodii HOT1A3"
      # treatment_taxid: 28108

  supplementary_materials:
    supp_table_1:
      type: csv
      filename: "data/.../results.csv"
      statistical_analyses:
        - id: "unique_id_within_publication"   # required; used as edge ID prefix
          experiment: my_experiment             # references experiments block
          name_col: "Gene"
          logfc_col: "log2FoldChange"
          adjusted_p_value_col: "padj"
          pvalue_threshold: 0.05
          # Optional time-series fields:
          # timepoint: "24h"
          # timepoint_hours: 24
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
| `gene_clusters` | Cluster assignment table; creates ClusteringAnalysis + GeneCluster nodes |
| `derived_metrics_table` | Column-level scalar summaries per gene (periodicity flags, classifiers, numeric scores); creates DerivedMetric nodes + one of 3 measurement edge types per metric (quantifies/flags/classifies) |

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

**`gene_clusters` example** — cluster assignment table with per-analysis metadata:
```yaml
med4_kmeans_nstarvation:
  type: gene_clusters
  name: "MED4 K-means N-starvation clusters"
  filename: "data/.../med4_kmeans_clusters.csv"
  organism: "Prochlorococcus MED4"
  gene_id_col: "gene_id"
  cluster_col: "cluster"
  cluster_type: "response_pattern"
  cluster_method: "K-means (K=9)"
  omics_type: MICROARRAY
  light_condition: "continuous light"
  treatment_type: ["nitrogen_stress"]
  treatment: "N-starvation time course (0, 3, 6, 12, 24, 48h)"
  experimental_context: "MED4 cells in Pro99; sampled at 0, 3, 6, 12, 24, 48h post N-deprivation"
  experiments:                  # optional: link to Experiment nodes
    - nitrogen_stress_experiment_key
```
The entry key (`med4_kmeans_nstarvation`) becomes the analysis key in node IDs. No `clusters:` block — cluster identities come from unique values in `cluster_col`. Per-cluster descriptions (`id`, `name`, `functional_description`, `temporal_pattern`, `expression_dynamics`) are supplied by a separately extracted `cluster_extraction_{entry_key}.json` file. Optional `score_col` and `p_value_col` fields carry membership scores onto `Gene_in_gene_cluster` edges.

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

# Refresh metabolism caches only (bypasses prepare_data.sh; avoids re-running NCBI/UniProt):
bash scripts/prepare_data.sh --steps 6 --force
```

Logs written to `logs/prepare_data_step0.log` … `logs/prepare_data_step6.log`. Monitor with `tail -f logs/prepare_data_step0.log`.

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

**Step 5** (`multiomics_kg/download/build_og_descriptions.py`) — extracts eggNOG ortholog group descriptions from the local `eggnog.db` SQLite database and writes a lightweight JSON cache at `cache/data/eggnog/og_descriptions.json` (~1 MB). This avoids needing the 39 GB database at KG build time (e.g., in Docker). The `ortholog_group_adapter` reads from this cache first, falling back to `eggnog.db` if the cache is missing. Run as module: `uv run python -m multiomics_kg.download.build_og_descriptions [--force]`. Requires step 2.

**Step 6** (`multiomics_kg/download/build_kegg_metabolism_xrefs.py`) — Phase 1.2.2 unified pruned KEGG cache + TCDB ontology assembly. Sub-steps:
1. Downloads TCDB reference TSVs → `cache/data/tcdb/raw/` (gitignored, mirrors `cache/data/kegg/raw/`).
2. Builds the TCDB hierarchy → `cache/data/tcdb/tcdb_hierarchy.json` (committed; parsing moved out of `build_mnx_resolver.py`).
3. Walks every strain's `gene_annotations_merged.json` to identify gene-reachable {KOs, reactions, compounds, pathways} AND gene-annotated TCDB IDs (bidirectionally — above + below the annotated level).
4. Prunes the TCDB hierarchy to the gene-reachable subhierarchy, then resolves TCDB substrate strings (`CHEBI:NNNN;name`) via the MNX resolver into Metabolite primary IDs.
5. Splits resolved substrate primaries into kegg.compound:* (extends the gene-reachable cpds set) and non-KEGG (chebi:*, mnx:*) buckets.
6. Computes transport-reachable pathways = ⋃ `compound→pathway` for every kegg.compound substrate; extends the gene-reachable pws set with these.
7. Prunes raw KEGG to the EXTENDED sets (catalysis ∪ transport-reachable), enriching reactions/compounds with MNX cross-refs. compound→pathway lists are filtered to extended pws (Option B).
8. Tags every compound entry with `evidence_sources` based on origin: `['metabolism']` if only catalysis-reachable, `['transport']` if only substrate-reachable, `['metabolism', 'transport']` if both. Non-KEGG primaries become `additional_compounds` entries with `['transport']`.
9. Writes `cache/data/tcdb/tcdb_pruned.json` (kept-node set + per-leaf metabolite primary IDs) and `cache/data/kegg/kegg_data.json` (~3-4 MB, indented JSON for git-friendly diffs).

Both `kegg_annotation_adapter` / `metabolism_adapter` (KEGG side) and `tcdb_adapter` (TCDB side) read these files — no per-adapter pruning at iteration time. The ~4 GB MNX resolver opens here only — it must be built first via `bash scripts/refresh_mnx.sh` (one-time setup; rerun only when MNX releases). The resolver location is configurable via the `MNX_DATA_DIR` env var (set in `.env`, e.g. `MNX_DATA_DIR=~/tools/mnx`) so a single ~4 GB cache can be shared between checkouts; defaults to `cache/data/mnx/` when unset. Run as module: `uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]`. Requires step 5 + the MNX resolver.

## Data Locations

- **Genomic data (raw):** `data/Prochlorococcus/genomes/<Strain>/` and equivalents for Synechococcus/Alteromonas
- **Genomic data (cached):** `cache/data/<Organism>/genomes/<Strain>/`
- **Publication data:** `data/Prochlorococcus/papers_and_supp/<Author Year>/`
- **Gene mapping CSVs:** `cache/data/<Organism>/genomes/<Strain>/gene_mapping.csv`
- **Gene annotation tables:** `cache/data/<Organism>/genomes/<Strain>/gene_annotations_merged.json`
- **Extended gene ID mapping:** `cache/data/<Organism>/genomes/<Strain>/gene_id_mapping.json` (built by step 3; includes paper-derived IDs)
- **Backward-compat ID mapping:** `cache/data/<Organism>/genomes/<Strain>/gene_mapping_supp.csv` (generated alongside gene_id_mapping.json)
- **UniProt cache (taxid-keyed):** `cache/data/<org_group>/uniprot/<taxid>/uniprot_preprocess_data.json`
- **Pfam reference cache:** `cache/data/pfam/pfam_reference.json`
- **eggNOG OG descriptions cache:** `cache/data/eggnog/og_descriptions.json` (built by step 5; ~1 MB)
- **KEGG unified pruned cache:** `cache/data/kegg/kegg_data.json` (built by step 6; ~3-4 MB indented; single source for `kegg_annotation_adapter` + `metabolism_adapter`)
- **TCDB committed outputs:** `cache/data/tcdb/{tcdb_hierarchy,tcdb_pruned}.json` (built by step 6; raw TSVs live under `cache/data/tcdb/raw/` and are gitignored)
- **MNX resolver:** location set by `MNX_DATA_DIR` env var (e.g. `~/tools/mnx`), default `cache/data/mnx/`. Holds 4 MNX TSVs + the ~4 GB SQLite resolver. Built by `bash scripts/refresh_mnx.sh`; shareable across checkouts via the env var.
- **PDF extraction cache:** `pdf_extraction_cache.json`

## KG Validity Tests

`tests/kg_validity/` contains a fast (~14s) pytest suite that connects to the live deployed Neo4j instance and validates the output graph. Tests auto-skip if Neo4j is unreachable. All tests are marked `@pytest.mark.kg`.

### Test files

| File | What it validates |
|---|---|
| `test_structure.py` | Node type presence, minimum counts (>50K genes, ≥25 organisms, ≥25 publications), orphan detection (genes without organism, proteins without organism), Gene_encodes_protein edge validation, key property presence |
| `test_biology.py` | Ecotype/clade labels per strain, katG absence in *Prochlorococcus* (Black Queen Hypothesis), all expected strains present, locus-tag → UniProt spot checks |
| `test_expression.py` | `log2_fold_change` / `adjusted_p_value` are numeric, `adjusted_p_value` ∈ [0,1], `expression_direction` ∈ {up,down}, direction/sign consistency, required properties on all edges |
| `test_organism.py` | No garbage taxonomy nodes, species on Alteromonas strains, all organisms bacterial (except Phage), organism count = 28 (23 strains + 5 treatment), precomputed stats consistency (gene_count, publication_count, experiment_count, treatment_types, omics_types), Publication.organisms ↔ OrganismTaxon.preferred_name alignment |
| `test_post_import.py` | OrthologGroup nodes exist with correct properties, Gene_in_ortholog_group membership edges exist, cross-strain bridging verified, no old homolog/ortholog edges remain |
| `test_snapshot.py` | Regression snapshot: verifies a sample of specific nodes and edges (with properties) still exist after rebuilds. Catches silent data loss. |
| `test_import_report.py` | Docker import artifact validation: reads `./output/import.status` + `./output/import.report`, fails if `neo4j-admin import` exited non-zero or skipped any relationships. Auto-skips when artifacts are missing. |

### Snapshot regression tests

`test_snapshot.py` loads `snapshot_data.json` (committed fixture) and checks that every sampled node/edge still exists with expected properties. To regenerate after intentional changes:

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

### Actual Neo4j labels (BioCypher PascalCase output)

- Nodes: `Gene`, `Protein`, `OrganismTaxon`, `Publication`, `Experiment`, `OrthologGroup`, `BiologicalProcess`, `CellularComponent`, `MolecularFunction`, `EcNumber`, `KeggTerm`, `CogFunctionalCategory`, `CyanorakRole`, `TigrRole`, `Pfam`, `PfamClan`, `TcdbFamily`, `CazyFamily`, `BriteCategory`, `ClusteringAnalysis`, `GeneCluster`, `Reaction`, `Metabolite`
- Relationships: `Gene_belongs_to_organism`, `Protein_belongs_to_organism`, `Gene_encodes_protein`, `Gene_in_ortholog_group`, `Og_has_cyanorak_role`, `Og_in_cog_category`, `Has_experiment`, `Tests_coculture_with`, `Changes_expression_of`, `Gene_involved_in_biological_process`, `Gene_located_in_cellular_component`, `Gene_enables_molecular_function`, `Biological_process_is_a_biological_process`, `Cellular_component_is_a_cellular_component`, `Molecular_function_is_a_molecular_function`, `Ec_number_is_a_ec_number`, `Gene_catalyzes_ec_number`, `Gene_has_kegg_ko`, `Kegg_term_is_a_kegg_term`, `Gene_in_cog_category`, `Gene_has_cyanorak_role`, `Gene_has_tigr_role`, `Cyanorak_role_is_a_cyanorak_role`, `Gene_has_pfam`, `Pfam_in_pfam_clan`, `Gene_has_tcdb_family`, `Tcdb_family_is_a_tcdb_family`, `Tcdb_family_transports_metabolite`, `Gene_has_cazy_family`, `Cazy_family_is_a_cazy_family`, `Brite_category_is_a_brite_category`, `Kegg_term_in_brite_category`, `Publication_has_clustering_analysis`, `Clustering_analysis_has_gene_cluster`, `Clusteringanalysis_belongs_to_organism`, `Experiment_has_clustering_analysis`, `Gene_in_gene_cluster`, `Gene_catalyzes_reaction`, `Reaction_has_metabolite`, `Reaction_in_kegg_pathway`, `Metabolite_in_pathway`, `Organism_has_metabolite`

### Key graph facts

- Gene↔Protein linkage: explicit `Gene_encodes_protein` edges (Gene→Protein) created by UniProt adapter via RefSeq protein_id join with gene_mapping.csv
- Experiment nodes: 167 nodes with properties: `name`, `organism_name`, `treatment_type` (str[]), `background_factors` (str[]), `treatment`, `control`, `experimental_context`, `coculture_partner`, `omics_type`, `statistical_test`, `is_time_course`, `medium`, `temperature`, `light_condition`, `light_intensity`, `table_scope`, `table_scope_detail`. Node IDs: `{doi}_{experiment_key}`. Each groups the analyses for one scientific comparison within a publication. `coculture_partner` is only present on coculture experiments (null/absent on non-coculture). `table_scope` describes what genes the source DE table contains: `all_detected_genes` | `significant_any_timepoint` | `significant_only` | `top_n` | `filtered_subset`. `table_scope_detail` provides free-text clarification for ambiguous cases. `treatment_type` canonical values: `nitrogen`, `phosphorus`, `iron`, `carbon`, `salt`, `light`, `temperature`, `plastic`, `darkness`, `diel`, `viral`, `coculture`, `growth_phase`, `compartment` (list — experiments may have multiple). `compartment` is for subcellular fraction comparisons (e.g. vesicle vs whole cell) where the actual fraction lives in the Experiment `compartment` field. `background_factors` canonical values: `axenic`, `light`, `diel`, `coculture`, `chemical`, `darkness`, `viral` (list of experimental context factors that aren't the primary treatment).
- Experiment computed properties (post-import): `gene_count` (int, **cumulative** count of `Changes_expression_of` edges — sums per-timepoint counts for time-course experiments, so for a paper measuring the same set of N genes at 6 TPs `gene_count = 6N`), `distinct_gene_count` (int, **distinct** gene count via `count(DISTINCT gene)` over the experiment's edges; this is the unique-gene scalar for detection-power / pathway-background reasoning; invariant `distinct_gene_count <= gene_count`, equality only for non-time-course or single-TP experiments), `significant_up_count` (int, edges where `expression_status = 'significant_up'`), `significant_down_count` (int, edges where `expression_status = 'significant_down'`), `time_point_count` (int), `time_point_labels` (str[], `""` = no label), `time_point_orders` (int[]), `time_point_hours` (float[], `-1.0` = unknown), `time_point_totals` (int[], per-tp gene counts), `time_point_significant_up` (int[], per-tp significant-up counts), `time_point_significant_down` (int[], per-tp significant-down counts). All time_point arrays are parallel (same length = `time_point_count`), ordered by `time_point_order`. Sentinel values used because Neo4j arrays cannot contain nulls. Non-time-course experiments get one entry with `time_point_labels = [""]`; experiments with 0 edges get `gene_count = 0`, `distinct_gene_count = 0`, and empty arrays. `not_significant_count` is not stored; derive as `gene_count - significant_up_count - significant_down_count`.
- Expression edges: ~227K `Changes_expression_of` edges (Experiment → Gene), unified single type replacing the old split (`Condition_changes_expression_of` + `Coculture_changes_expression_of`). Edge properties: `time_point` (str), `time_point_order` (int), `time_point_hours` (float), `log2_fold_change` (float), `adjusted_p_value` (float), `expression_direction` (str), `significant` (str), `expression_status` (str, post-import derived: `"significant_up"` | `"significant_down"` | `"not_significant"`), `rank_by_effect` (int, post-import: rank by |log2FC| among all genes per experiment x timepoint, 1 = strongest), `rank_up` (int|null, post-import: rank by |log2FC| among significant_up genes per experiment x timepoint; null if not significant_up), `rank_down` (int|null, post-import: rank by |log2FC| among significant_down genes per experiment x timepoint; null if not significant_down). Metadata (organism, omics type, conditions) lives on the Experiment node, not the edge.
- `Has_experiment` edges: 167 — `(Publication) → (Experiment)`, one per experiment node
- `Tests_coculture_with` edges: ~80 — `(Experiment) → (OrganismTaxon)`, present only for coculture/viral experiments, links to the treatment organism
- `rank_by_effect` property on expression edges: computed by post-import Cypher, ranks genes within each (experiment, time_point_order) group by descending |log2_fold_change| (1 = strongest effect)
- OrthologGroup nodes: ~45K nodes with `name` (raw OG ID, e.g. "CK_00000364" or "COG0592@2"), `source` ("cyanorak"|"eggnog"), `taxonomic_level` (e.g. "curated", "Prochloraceae", "Proteobacteria", "Bacteria"), `taxon_id` (int), `specificity_rank` (0=curated, 1=family, 2=intermediate, 3=domain), `consensus_product` (majority-vote product from members), `consensus_gene_name` (most frequent gene name), `description` (str, nullable — eggNOG functional narrative from local `eggnog.db`; null for Cyanorak groups), `functional_description` (str, nullable — semicolon-separated majority-vote CyanorakRole hierarchical names + CogFunctionalCategory names from member genes), `member_count` (int), `organism_count` (int), `genera` (str[], e.g. ["Prochlorococcus", "Alteromonas"]), `has_cross_genus_members` ("cross_genus"|"single_genus"). Sources: Cyanorak curated clusters (Pro/Syn), eggNOG family-level OGs (Prochloraceae, Synechococcus, Alteromonadaceae), eggNOG intermediate bridging OGs (Cyanobacteria for all cyanos, Proteobacteria for all proteobacteria), eggNOG Bacteria-level COGs. `ORGANISM_GROUP_LEVELS` in `ortholog_group_utils.py` maps 8 organism groups to their configured levels. Expression propagation is query-time via 2-hop (`gene->OG<-gene`) instead of materialized ortholog edges.
- `Gene_in_ortholog_group` edges: ~223K membership edges. A gene may have 2-4 memberships (Cyanorak + family + intermediate + bacteria-level). Intermediate bridging levels enable cross-genus comparison: Cyanobacteria (1117) bridges Prochlorococcus + Synechococcus + Thermosynechococcus; Proteobacteria (1224) bridges Alteromonas + Shewanella + Pseudomonas + Ruegeria. Meiothermus has only Bacteria-level OGs.
- `Og_has_cyanorak_role` edges: OrthologGroup → CyanorakRole, majority-vote (>50% of member genes). `Og_in_cog_category` edges: OrthologGroup → CogFunctionalCategory, same majority rule.
- Pfam domain nodes: ~5,500 `Pfam` nodes (only domains referenced by genes) with `name` (description), `short_name` (Pfam shortname); ~509 `PfamClan` nodes (superfamilies) with `name` (clan name). ~89K `Gene_has_pfam` edges, ~3.5K `Pfam_in_pfam_clan` edges. Node IDs use bioregistry CURIEs: `pfam:PF*` and `pfam.clan:CL*`. Pfam reference data downloaded from `Pfam-A.clans.tsv.gz`; eggNOG shortnames resolved to PF* accessions via reverse lookup
- TcdbFamily nodes: pruned to gene-reachable subhierarchy (above + below); 5 levels (tc_class | tc_subclass | tc_family | tc_subfamily | tc_specificity). Properties: name, tcdb_id, level (0-4), level_kind, superfamily (sparse, leaf-only), tc_class_id (post-import: pointer to root tc_class). Computed post-import: gene_count, organism_count, member_count, metabolite_count.
- CazyFamily nodes: observed-only; 3 levels (cazy_class | cazy_family | cazy_subfamily). Properties: name, cazy_id, level (0-2), level_kind. Computed post-import: gene_count, organism_count.
- Five new ontology edges: Gene_has_tcdb_family (Gene → TcdbFamily, attaches at exact eggNOG-annotated level), Tcdb_family_is_a_tcdb_family (child → parent), Tcdb_family_transports_metabolite (TcdbFamily-leaf → Metabolite, only on tc_specificity nodes), Gene_has_cazy_family (most specific observed level), Cazy_family_is_a_cazy_family.
- Metabolite.evidence_sources (str[]): values 'metabolism' (Reaction-reachable) and 'transport' (TCDB-substrate-reachable). Open-ended; 'metabolomics' reserved for a future metabolomics-DM spec.
- Metabolite.transporter_count (int, post-import): distinct TcdbFamily nodes pointing at this metabolite via Tcdb_family_transports_metabolite.
- Metabolite.gene_count / organism_count: now UNION across catalysis (Gene → Reaction → Metabolite) + transport (Gene → TcdbFamily → Metabolite) paths.
- Gene.tcdb_family_count / cazy_family_count: distinct ontology edges per gene (post-import). Gene.metabolite_count: UNION across catalysis + transport paths. Gene.annotation_types extended with 'tcdb' / 'cazy'.
- Organism_has_metabolite materialization: now includes the transport arm (Gene → TcdbFamily → tc_specificity → Metabolite). Single MERGE deduplicates by (organism, metabolite).
- Removed Gene properties (replaced by graph edges): transporter_classification (use Gene_has_tcdb_family), cazy_ids (use Gene_has_cazy_family).
- Reaction nodes (~3K) with primary ID `kegg.reaction:R*`. Properties: kegg_reaction_id, name, ec_numbers, kegg_pathway_ids, mnxr_id (nullable), rhea_ids, reaction_class ('transport'|'chemical'), mass_balance ('balanced'|'unbalanced'). Computed post-import: gene_count, organism_count, organisms[]. KEGG-native primary IDs (Spec 1.2 pivot from MNX-as-primary).
- Metabolite nodes (~1.5K) with primary ID `kegg.compound:C*` (with chebi: / mnx: fallback for compounds without KEGG entry). Properties: kegg_compound_id, name, formula, elements (str[], sparse — sorted unique element symbols Hill-parsed from formula via chemparse; absent when formula is null or contains no element symbols), mass, inchikey, smiles, mnxm_id (nullable), chebi_id (nullable), hmdb_id (nullable). Computed post-import: gene_count, organism_count.
- `Gene_catalyzes_reaction` edges from gene["kegg_reactions"]; `Reaction_has_metabolite` edges from KEGG /link/compound/reaction (direction-agnostic — Spec 1.2 explicitly does NOT model reaction direction); `Reaction_in_kegg_pathway` edges reuse existing KeggTerm pathway nodes (kegg.pathway:ko*); `Metabolite_in_pathway` edges (Metabolite → KeggTerm pathway) materialize compound→pathway annotations from KEGG, pruned to KO ∪ Reaction-reachable pathways (Option B; drops compound-only pathways like Leishmaniasis or MAPK signalling that are biologically inappropriate for marine bacteria); `Organism_has_metabolite` edges materialized post-import as a 2-hop save (Gene→Reaction→Metabolite, distinct).
- BriteCategory nodes: ~2,560 nodes from 12 KEGG BRITE functional hierarchy trees. Node IDs: `kegg.brite:{tree_code}.{A#}.{B#}...` (positional, 1-based, no tree-root node emitted). Properties: `name` (str), `tree` (snake_case tree name, e.g. `"transporters"`), `tree_code` (e.g. `"ko02000"`), `level` (int, 0=broadest/A-level=brite_class, 1=brite_subclass, 2=brite_family, 3=brite_subfamily), `level_kind` (str). Computed post-import: `member_ko_count` (int, recursive KO leaf count), `gene_count` (int, distinct genes reachable via KO leaves), `organism_count` (int). 12 trees: enzymes (ko01000), transporters (ko02000), peptidases (ko01002), transcription_factors (ko03000), secretion (ko02044), two_component (ko02022), defense (ko02048), chaperones (ko03110), ribosome (ko03011), translation_factors (ko03012), trna_biogenesis (ko03016), dna_replication (ko03032). Edges: ~2,500 `Brite_category_is_a_brite_category` (child→parent hierarchy within each tree) + ~4,100 `Kegg_term_in_brite_category` (KeggTerm KO → BriteCategory leaf). KO nodes connect to their immediate BRITE parent category. A-level entries are level=0 with no parent edge. BRITE output is pruned to the subhierarchy reachable from `KeggTerm` nodes — a `BriteCategory` node appears iff at least one KO in its subtree is referenced by a gene in the KG. Pruning is keyed off `MultiKeggAnnotationAdapter.all_ko_ids()` (see `docs/superpowers/specs/2026-04-14-brite-pruning-design.md`).
- ClusteringAnalysis nodes: intermediate layer between Publication and GeneCluster. Each groups clusters from one clustering analysis entry in a paperconfig. Properties: `name`, `organism_name`, `cluster_method`, `cluster_type`, `cluster_count`, `total_gene_count`, `omics_type`, `treatment_type` (str[]), `background_factors` (str[]), `treatment`, `light_condition`, `experimental_context`. Node IDs: `clustering_analysis:{doi_short}:{entry_key}`. Linked via `Publication_has_clustering_analysis` (Publication → ClusteringAnalysis), `Clusteringanalysis_belongs_to_organism` (ClusteringAnalysis → OrganismTaxon), and `Experiment_has_clustering_analysis` (Experiment → ClusteringAnalysis).
- GeneCluster nodes: one per unique cluster value in the CSV. Node IDs: `cluster:{doi_short}:{analysis_entry_key}:{csv_cluster_key}`. Properties: `id` (extracted short identifier), `name`, `organism_name`, `member_count`, `functional_description`, `temporal_pattern`, `expression_dynamics`. Per-cluster descriptions (`id`, `name`, `functional_description`, `temporal_pattern`, `expression_dynamics`) come from extraction JSON files (`cluster_extraction_{entry_key}.json`); only applied when stage3 validation verdict = "pass". Experimental metadata (cluster_method, cluster_type, treatment_type, omics_type, etc.) lives on the parent ClusteringAnalysis node, not on GeneCluster. Linked to analysis via `Clustering_analysis_has_gene_cluster`; linked to genes via `Gene_in_gene_cluster`.
- DerivedMetric nodes: 13 (Biller 2018: 6 boolean + 1 categorical; Waldbauer 2012: 6 numeric). Node IDs: `derived_metric:{doi_short}:{entry_key}:{metric_type}`. Properties: `name`, `metric_type`, `value_kind ∈ {numeric, boolean, categorical}`, `rankable` (str `"true"`/`"false"`), `has_p_value` (str `"true"`/`"false"`), `p_value_threshold` (float, nullable), `unit`, `allowed_categories` (str[]), `field_description`, `compartment`, plus denormalized parent-Experiment fields (`experiment_id`, `organism_name`, `publication_doi`, `omics_type`, `treatment_type` (str[]), `background_factors` (str[]), `treatment`, `light_condition`, `experimental_context`). Each DM emits exactly one of 3 measurement edge types: `Derived_metric_quantifies_gene` (numeric — 1,872 edges in Waldbauer 2012, 312 genes × 6 metrics), `Derived_metric_flags_gene` (boolean — 4,160 edges in Biller 2018), `Derived_metric_classifies_gene` (categorical — 258 edges in Biller 2018). Binding edges: `PublicationHasDerivedMetric`, `ExperimentHasDerivedMetric`, `DerivedMetricBelongsToOrganism`. See `docs/kg-changes/non-de-evidence-extension.md`.
- DerivedMetric computed properties (post-import): `total_gene_count` (int), `growth_phases` (str[] from parent Experiment).
- Numeric DM edge computed properties (post-import, only when parent `rankable="true"`): `rank_by_metric` (int, 1 = highest value), `metric_percentile` (float 0-100), `metric_bucket` (str: `top_decile`/`top_quartile`/`mid`/`low`, pinned thresholds ≥90/≥75/≥25). When parent `has_p_value="true"` AND threshold is set AND `adjusted_p_value` non-null: `significant` (str `"true"`/`"false"`).
- Experiment `compartment` property (str, adapter-emitted, default `"whole_cell"`). Vocab: `whole_cell`, `vesicle`, `exoproteome`, `secretome`. Computed Experiment DM rollups (post-import): `reports_fold_change` (str `"true"`/`"false"`), `reports_derived_metric_types` (str[]), `derived_metric_count` (int), `derived_metric_value_kinds` (str[]), `derived_metric_gene_count` (int), `derived_metric_search_text` (str, null on experiments with no DMs — aggregated DM `name` + tokenized `metric_type` + `field_description`; consumed by `experimentFullText`).
- Publication/OrganismTaxon DM rollups (post-import): `derived_metric_count` (int), `derived_metric_gene_count` (int), `compartments` (str[]), `derived_metric_types` (str[]), `derived_metric_value_kinds` (str[]). Publication also gets `derived_metric_search_text` (str, null on publications with no DMs — same aggregation as Experiment; consumed by `publicationFullText`).
- Gene DM routing signals (post-import): `numeric_metric_count`, `boolean_metric_count`, `categorical_metric_count` (ints, default 0); `numeric_metric_types_observed`, `boolean_metric_types_observed`, `categorical_metric_types_observed` (str[], default []); `compartments_observed` (str[], default []).
- EnvironmentalCondition nodes no longer exist — replaced by Experiment nodes (see above)
- `preferred_name` property on `OrganismTaxon` nodes (e.g., `"Prochlorococcus MED4"`)
- `organism_type` property on `OrganismTaxon` nodes: `"genome_strain"` (25 organisms — real genome assemblies from cultured isolates), `"treatment"` (5 organisms — non-genomic coculture partners), `"reference_proteome_match"` (2 organisms — identified by matching experimental data against a multi-organism reference DB). Reference proteome match organisms also carry `reference_database` (e.g., `"MarRef v6"`) and `reference_proteome` (assembly accession). See `docs/kg-changes/reference-proteome-match-organisms.md`.
- OrganismTaxon computed properties (post-import): `gene_count` (int), `publication_count` (int), `experiment_count` (int), `treatment_types` (str[]), `omics_types` (str[]), `background_factors` (str[]), `metabolite_count` (int), `reaction_count` (int). Computed after Publication.organisms alignment. `species` is derived from `preferred_name` when NCBI taxonomy is genus-level only (e.g., Alteromonas strains → "Alteromonas macleodii").
- Gene node computed fields for MCP gene lookup: `organism_name` (preferred organism name), `gene_summary` ("name :: product :: description"), `all_identifiers` (union of all alt IDs)
- Gene nodes carry ~27 properties. Kept properties: core (locus_tag, start, end, strand, product, protein_id), naming (gene_name, gene_name_synonyms, alternate_functional_descriptions, function_description), function (protein_family, catalytic_activities, transmembrane_regions, signal_peptide, bigg_reaction), computed (organism_name, gene_summary, all_identifiers), quality (annotation_quality, gene_category), routing signals (annotation_types, expression_edge_count, significant_up_count, significant_down_count, reaction_count, metabolite_count, tcdb_family_count, cazy_family_count, closest_ortholog_group_size, closest_ortholog_genera, cluster_membership_count, cluster_types — set by post-import Cypher, not adapter). reaction_count is the count of Gene_catalyzes_reaction edges; metabolite_count is the distinct count of Metabolite nodes reachable via any gene-reaching path (UNION across catalysis + transport, per TCDB-S3 / KG-A2). Redundant ID arrays removed: `gene_synonyms` (= alternative_locus_tags ∪ gene_name_synonyms), `alternative_locus_tags` (⊂ all_identifiers), `old_locus_tags` (⊂ alternative_locus_tags ⊂ all_identifiers). See `docs/gene_id_fields_spec.md` for full analysis.
- Publication computed properties (post-import): `experiment_count` (int), `treatment_types` (str[]), `omics_types` (str[]), `background_factors` (str[]), `organisms` (str[] — sorted distinct `organism_name` + `coculture_partner` from experiments). The adapter no longer sets organism; it is fully computed post-import as `organisms`.
- Post-import indexes: scalar (`gene_locus_tag_idx`, `gene_name_idx`, `gene_organism_name_idx`, `ortholog_group_id_idx`, `ortholog_group_name_idx`, `ortholog_group_level_idx`, `ortholog_group_rank_idx`, `pfam_name_idx`, `pfam_clan_name_idx`, `tcdb_family_level_idx`, `tcdb_family_level_kind_idx`, `tcdb_family_tcdb_id_idx`, `tcdb_family_tc_class_id_idx`, `cazy_family_level_idx`, `cazy_family_level_kind_idx`, `cazy_family_cazy_id_idx`, `brite_category_tree_idx`, `brite_category_level_idx`, `brite_category_name_idx`, `experiment_id_idx`, `experiment_organism_idx`, `experiment_treatment_type_idx`, `experiment_omics_type_idx`, `experiment_background_factors_idx`, `experiment_compartment_idx`, `clustering_analysis_organism_idx`, `clustering_analysis_method_idx`, `clustering_analysis_type_idx`, `derived_metric_metric_type_idx`, `derived_metric_value_kind_idx`, `derived_metric_compartment_idx`, `derived_metric_omics_type_idx`, `derived_metric_treatment_type_idx`, `derived_metric_organism_idx`, `derived_metric_experiment_idx`, `reaction_id_idx`, `reaction_kegg_id_idx`, `reaction_mnxr_idx`, `metabolite_id_idx`, `metabolite_kegg_id_idx`, `metabolite_mnxm_idx`, `metabolite_chebi_idx`) + full-text (`geneFullText` on gene_summary, all_identifiers, gene_name_synonyms, alternate_functional_descriptions; `orthologGroupFullText` on OrthologGroup consensus_product, consensus_gene_name, description, functional_description; `pfamFullText` on Pfam name, short_name; `pfamClanFullText` on PfamClan name; `tcdbFamilyFullText` on TcdbFamily name; `cazyFamilyFullText` on CazyFamily name; `briteCategoryFullText` on BriteCategory name; `experimentFullText` on Experiment name, treatment, control, experimental_context, light_condition, compartment, derived_metric_search_text; `publicationFullText` on Publication title, abstract, description, compartments, derived_metric_search_text; `clusteringAnalysisFullText` on ClusteringAnalysis name, treatment, experimental_context; `derivedMetricFullText` on DerivedMetric name, field_description; `reactionFullText` on Reaction name; `metaboliteFullText` on Metabolite name, formula; `biologicalProcessFullText`, `molecularFunctionFullText`, `cellularComponentFullText`, `ecNumberFullText`, `keggFullText`, `cogCategoryFullText`, `cyanorakRoleFullText`, `tigrRoleFullText` on ontology/role node `name` properties)
- Unified ontology `level` property: every term on `BiologicalProcess`, `MolecularFunction`, `CellularComponent`, `EcNumber`, `KeggTerm`, `CyanorakRole`, `TigrRole`, `CogFunctionalCategory`, `Pfam`, `PfamClan`, `BriteCategory` carries `level: int` with `0` = broadest/root. GO DAG terms whose min-path ≠ max-path to the namespace root additionally carry `level_is_best_effort = "true"` (sparse — property absent otherwise; check with `t.level_is_best_effort IS NOT NULL`). KEGG's original string `level` has been renamed to `level_kind` (`category|subcategory|pathway|ko`); the new integer `level` maps `category=0, subcategory=1, pathway=2, ko=3`. See `docs/kg-changes/ontology-level.md`.
- KeggTerm pathway-level nodes (level_kind = 'pathway', ~377 nodes) carry post-import `reaction_count: int` (incoming Reaction_in_kegg_pathway edges) and `metabolite_count: int` (incoming Metabolite_in_pathway edges). Sparse — only on pathway-level KeggTerm nodes; KOs and category nodes have neither.
- `adjusted_p_value` may be null on expression edges when the original study did not report it
- Strains in graph: MED4, AS9601, MIT9301, MIT9312, MIT9313, MIT9303, NATL1A, NATL2A, RSP50, **SS120** (Prochlorococcus); CC9311, WH7803, PCC7002, PCC7942, UTEX2973, **BL107** (Synechococcus); WH8102 (Parasynechococcus); BP1 (Thermosynechococcus); MIT1002, EZ55, HOT1A3 (Alteromonas); W3-18-1 (Shewanella); KT2440 (Pseudomonas); DSS-3 (Ruegeria); MruberA (Meiothermus). Reference proteome match organisms (not cultured strains): **Marinobacter (MarRef v6)** (HP15 reference proteome, GCF_000166295.1), **Alteromonas (MarRef v6)** (Alt_MarRef, GCA_003513035.1, DEH24_ locus tags). New as of 2026-04-15: SS120, BL107, HP15, AltMedDE (added for Domínguez 2017 + Fuszard 2012 + Moreno 2023 integration); HP15 renamed to Marinobacter (MarRef v6), AltMedDE renamed to Alt_MarRef as of 2026-04-16 — see `docs/community_proteomics_marref_saga.md` and `docs/kg-changes/reference-proteome-match-organisms.md`.

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
- `chemparse` — pure-Python Hill-notation formula parser (used in metabolism_adapter for Metabolite.elements)
- Current organisms: 10 *Prochlorococcus* (+SS120) + 7 *Synechococcus/Parasynechococcus/Thermosynechococcus* (+BL107) + 3 *Alteromonas* + 5 heterotrophs (Shewanella, Pseudomonas, Ruegeria, Meiothermus) + 2 reference proteome match organisms (**Marinobacter (MarRef v6)**, **Alteromonas (MarRef v6)**); 32 OrganismTaxon nodes total (25 genome strains + 5 treatment organisms + 2 reference proteome match); primary focus is MED4 (NCBI taxid 59919)
- API keys go in `.env`: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `ELSEVIER_API_KEY`, plus others for LangChain tracing and search
- `@pytest.mark.slow` tests run the full pipeline as a subprocess and take ~1 hour
- `@pytest.mark.kg` tests require a running Neo4j instance (Docker graph); skip automatically if not available

## Development Workflow

Every development task — feature, bug fix, refactor, or new data integration — must follow all five phases below. Do not skip phases or combine them implicitly.

### 1. Scope

Before writing any code, define the scope explicitly:
- **What** is changing (files, adapters, schema, config)
- **What is NOT changing** (boundaries prevent scope creep)
- **Acceptance criteria** — concrete, testable conditions for "done"
- For non-trivial tasks, create a plan file in `plans/` and get alignment before proceeding

### 2. Implement

- Make the code changes
- Follow existing patterns (adapter pattern, ID conventions, string sanitization)
- Keep changes minimal — only what the scope requires

### 3. Test

- Write or update unit tests for every code change
- Run the relevant test suite before considering the task complete:
  ```bash
  # Unit tests (fast, always run)
  pytest -m "not slow and not kg" -v

  # KG validity tests (when graph structure changes, requires running Neo4j)
  pytest -m kg -v
  ```
- For schema or adapter changes: use `/omics-edge-snapshot` before and after to verify no edge regressions
- For new paper integrations: use `/check-gene-ids` to validate match rates

### 4. Review

- Verify string literals are consistent across schema, adapters, and tests
- Check that no old vocabulary or dead code remains
- Confirm no unintended side effects on existing data (e.g., edge count changes)
- For schema changes: verify `schema_config.yaml`, adapter code, post-import scripts, and test fixtures all agree

### 5. Document

Update documentation to reflect the changes:
- **`CLAUDE.md`** — if the change affects build commands, architecture, key facts, or known issues
- **Plan files** (`plans/`) — mark completed phases, record decisions
- **Skill files** (`.claude/skills/`) — if the change affects skill behavior or templates
- **Memory** — if the change introduces new patterns or facts that future conversations need
