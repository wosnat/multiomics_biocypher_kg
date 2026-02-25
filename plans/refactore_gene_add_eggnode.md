# Plan: EggNOG Integration + Gene Annotation Pipeline

## Context

The KG draws gene/protein annotations from CyanorakNcbi (manually curated), UniProt (reviewed), and NCBI GFF (automated). All three are partial. EggNOG mapper has been run on all 13 strains and provides the broadest computational coverage. Goal: (0) refactor download logic into a standalone genome-data pipeline script; (1) merge all annotation sources into a per-genome wide table then collapse to canonical fields with priority Cyanorak > UniProt > NCBI > EggNOG; (2) clean up the Gene/Protein node schema using the full property mapping below; (3) optionally add an LLM-generated summary field per gene.

---

## Step 0 — Download Pipeline Script

**New file:** `multiomics_kg/download/download_genome_data.py`

Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` and for each genome:

| Step | Tool | Cache path | Notes |
|---|---|---|---|
| 1 | Download NCBI GFF | `cache/data/<org>/genomes/<strain>/<accession>.gff` | extracted from `cyanorak_ncbi_adapter._download_ncbi_genome()` |
| 2 | Download Cyanorak GFF/GBK | `cache/data/<org>/genomes/<strain>/cyanorak*.gff` | extracted from `cyanorak_ncbi_adapter._download_cyanorak_*()` |
| 3 | Download UniProt (per unique taxid) | `cache/data/<org_group>/uniprot/<taxid>/uniprot_raw_data.json` | extracted from `uniprot_adapter.download_uniprot_data()`; shared taxids downloaded once |
| 4 | Run eggnog-mapper | `cache/data/<org>/genomes/<strain>/eggnog/<strain>.emapper.annotations` | skip if exists |
| 5 | _(future)_ SignalP | `cache/data/<org>/genomes/<strain>/signalp/` | |
| 6 | _(future)_ DeepTMHMM | `cache/data/<org>/genomes/<strain>/tmhmm/` | |
| 7 | _(future)_ OrthoFinder | `cache/data/<org>/orthofinder/` | multi-genome, run once |

**UniProt taxid deduplication:** `cyanobacteria_genomes.csv` contains:
- 10 unique taxids across Prochlorococcus/Synechococcus strains
- Alteromonas MIT1002, EZ55, HOT1A3 **all share taxid 28108** → UniProt downloaded once to `cache/data/Alteromonas/uniprot/28108/`
- General pattern: `cache/data/<organism_group>/uniprot/<taxid>/` where `organism_group` is the top-level dir (Prochlorococcus, Synechococcus, Alteromonas)

**UniProt adapter changes:** `uniprot_adapter.download_uniprot_data()` is updated to:
- Accept an explicit `cache_dir` path instead of writing to the project root
- Read from `cache_dir` if it exists (skip download)
- The download pipeline passes the correct taxid-keyed path

**Design:** The adapter `download_data()` methods are refactored to be "load from cache only" — they assume the download pipeline has already populated the cache. The pipeline script becomes the single entry point for genome data acquisition. Running it with `--force <strain>` triggers re-download + re-annotation for a specific genome.

The `create_knowledge_graph.py` pipeline does NOT call the download pipeline automatically (it assumes cache is ready); the download pipeline is run manually when adding new genomes or re-annotating.

**CLI flags:**
- `--steps <N> [<N> ...]` — run only the specified step numbers (1–7 above); default: run all
- No explicit `--force`: by default all steps skip if the target cache file already exists; `--force` overrides this for specified strains
- Example: `uv run python multiomics_kg/download/download_genome_data.py --steps 1 2 3` (skip eggnog re-run)

---

## Step 1 — Merge Script: `multiomics_kg/download/build_gene_annotations.py`

For each strain in `cyanobacteria_genomes.csv`:

**Input sources and join strategy:**
1. `gene_mapping.csv` — base (38 cols: locus_tag + all Cyanorak/NCBI fields)
2. `eggnog/<strain>.emapper.annotations` — join on `protein_id` = `query` col
3. `uniprot_preprocess_data.json` — join via `xref_refseq` → `protein_id` (adds UniProt function_description, catalytic_activity, subcellular_location, is_reviewed, annotation_score, transmembrane_regions, signal_peptide)

**Step A — Build wide dict** (all source fields kept, source-prefixed), stored as JSON:
```
→ cache/data/<org>/genomes/<strain>/gene_annotations_wide.json
```
Format: `{ "PMM0001": { "ncbi_product": "...", "cyanorak_product": "...", "eggnog_Description": "...", "eggnog_GOs": [...], ... }, ... }`
Purpose: audit trail; inspect per-gene; allows changing merge rules without re-downloading.

**Step B — Apply merge rules → canonical fields**, also stored as JSON:
```
→ cache/data/<org>/genomes/<strain>/gene_annotations_merged.json
```
Format: `{ "PMM0001": { "gene_name": "katG", "product": "catalase-peroxidase", "go_terms": ["GO:0006979"], "product_source": "cyanorak", ... }, ... }`
Arrays stored as native JSON lists. Null fields omitted (sparse). Consistent with existing `uniprot_preprocess_data.json` convention.

**Step C — LLM summary** (optional flag `--llm-summary`):
For each gene with non-empty `function_description` or `product`, call Claude API to generate a 2–3 sentence plain-English gene summary incorporating: gene_name, product, function_description, cyanorak_role, cog_category, kegg_pathway.
- Store in `llm_summary` field in merged CSV
- Cache in `cache/data/<org>/genomes/<strain>/llm_summaries.json` (locus_tag → summary)
- Model: `claude-haiku-4-5-20251001` (cost-efficient for bulk annotation)
- Run only once per gene; skip if cached

**Step D — Coverage report:**
Print per-strain stats: gene count, eggnog match rate (%), fill rates for description, go_terms, cog_category, kegg_ko.

---

## Step 2 — Full Property Mapping Table

The following table maps every current property across all sources to the proposed new schema. Provenance is stored as an `{field}_source` field for key canonical fields (see bottom of table).

### A. Gene Identity & Locus

| Current (source) | New name | Action | Notes |
|---|---|---|---|
| `locus_tag` (gene_mapping) | `locus_tag` | KEEP | Primary node ID |
| `locus_tag_ncbi` (gene_mapping) | `locus_tag_ncbi` | KEEP | For ID mapping |
| `locus_tag_cyanoak` (gene_mapping) | `locus_tag_cyanorak` | KEEP + fix typo | |
| `protein_id` (gene_mapping) | `protein_id` | KEEP | RefSeq WP_ join key |
| `old_locus_tags` (gene_mapping) | `old_locus_tags` | KEEP | For MIT9313 remapping |
| `source`, `Note`, `exception`, `inference`, `gene` (gene_mapping) | — | DROP | Internal NCBI fields, usually empty |

### B. Genomic Location

| Current (source) | New name | Action |
|---|---|---|
| `start`, `end`, `strand` (NCBI) | `start`, `end`, `strand` | KEEP |
| `start_cyanorak`, `end_cyanorak`, `strand_cyanorak` | `start_cyanorak`, `end_cyanorak`, `strand_cyanorak` | KEEP |

### C. Gene Naming (hierarchical — increasing specificity)

| Level | New field | Priority merge rule |
|---|---|---|
| 1 | `gene_name` (str) | Cyanorak `gene_names_cyanorak`[0] > UniProt `gene_symbol` > NCBI `gene_names`[0] > EggNOG `Preferred_name` |
| 2 | `gene_synonyms` (str[]) | UNION of all alternative names from all sources not selected as `gene_name`; explicitly includes UniProt `protein_gene_names` (all UniProt gene name synonyms) |
| — | `gene_name_source` (str) | Which source provided `gene_name` ("cyanorak"/"uniprot"/"ncbi"/"eggnog") |

Replaces: `gene_names`, `gene_names_cyanorak`, `gene_synonym`, `gene_symbol` (on Protein). Cyanorak-specific names still in `gene_synonyms`.

### D. Functional Description (hierarchical — increasing detail)

| Level | New field | Priority merge rule |
|---|---|---|
| 1 | `product` (str) | Cyanorak `product_cyanorak` > NCBI `product` > EggNOG `Description` |
| 2 | `function_description` (str) | UniProt `function_description` > EggNOG `Description` (fallback when longer than product) |
| 3 | `llm_summary` (str) | LLM-generated (see Step 1C); null if not run |
| — | `product_source` (str) | Which source provided `product` |
| — | `function_description_source` (str) | Which source provided `function_description` |

Keep source-specific fields as-is on Gene node for provenance: `product_cyanorak` (Cyanorak curated), `product` (NCBI). Drop neither — they remain queryable.

### E. Categorical Functional Roles

| Current (source) | New name | Action | Notes |
|---|---|---|---|
| `cyanorak_Role`, `cyanorak_Role_description` | same | KEEP | Manually curated TIGR role categories |
| `tIGR_Role`, `tIGR_Role_description` | same | KEEP | TIGR functional role hierarchy |
| EggNOG `COG_category` | `cog_category` | NEW | COG letter codes, e.g. "C,G" |
| Protein `protein_family` | `protein_family` | DENORMALIZE to Gene | Useful for quick classification |

### F. Gene Ontology

| Current (source) | New name | Action | Notes |
|---|---|---|---|
| `go_component`, `go_function`, `go_process` (NCBI) | `go_terms` | MERGE | Consolidate all GO into one array |
| `Ontology_term_ncbi`, `Ontology_term_cyanorak` | `go_terms` + `ontology_terms` | MERGE GO ones into go_terms; keep non-GO in `ontology_terms` |
| `ontology_term_description` | `go_term_descriptions` | RENAME + KEEP |
| UniProt `go_cellular_components`, `go_biological_processes`, `go_molecular_functions` | `go_terms` | MERGE | Highest quality GO annotations |
| EggNOG `GOs` | `go_terms` | MERGE | Fill gaps for non-UniProt genes |
| Gene `go_biological_processes` (denormalized) | (removed, now in `go_terms`) | CONSOLIDATE |

Replaces 6 separate GO fields with one `go_terms` array (all namespaces) + `go_term_descriptions`. GO namespace categorization lives in the GO node hierarchy (existing `go_adapter`).

### G. EC Numbers

| Current (source) | New name | Action |
|---|---|---|
| `ec_numbers` (gene_mapping/NCBI) | `ec_numbers` | UNION all sources |
| Protein `ec_numbers` | — | merge |
| EggNOG `EC` | — | merge |

### H. KEGG

| Current (source) | New name | Action | Notes |
|---|---|---|---|
| `kegg`, `kegg_description` (gene_mapping) | `kegg_ko`, `kegg_ko_descriptions` | RENAME + MERGE | KO IDs (K00001 format) |
| Protein `kegg_ids`, `kegg_ko_ids` | — | merge into `kegg_ko` | Both contain KO IDs |
| EggNOG `KEGG_ko` | — | merge into `kegg_ko` | |
| EggNOG `KEGG_Pathway` | `kegg_pathway` | NEW | Pathway IDs (map.00010 format) |
| EggNOG `KEGG_Module` | `kegg_module` | NEW | Module IDs |
| EggNOG `KEGG_Reaction` | `kegg_reaction` | NEW | Reaction IDs (optional) |
| EggNOG `BRITE` | `kegg_brite` | NEW | BRITE hierarchy (optional, bulky) |
| EggNOG `KEGG_TC` | `transporter_classification` | NEW | TCDB TC numbers (e.g. 1.A.1.1.1) — keep for transporter analysis |
| EggNOG `BiGG_Reaction` | `bigg_reaction` | NEW | Metabolic model reaction IDs |

### I. EggNOG Orthology

| Current (source) | New name | Action |
|---|---|---|
| `eggNOG`, `eggNOG_description` (gene_mapping, from NCBI GFF) | `eggnog_ogs`, `eggnog_og_descriptions` | RENAME + MERGE with EggNOG mapper `eggNOG_OGs` |
| Protein `eggnog_ids` | — | merge into `eggnog_ogs` |
| EggNOG `seed_ortholog` | `seed_ortholog` | NEW |
| EggNOG `max_annot_lvl` | `max_annot_lvl` | NEW (deepest taxonomic level with annotation) |
| EggNOG `evalue` | `seed_ortholog_evalue` | NEW (confidence metric) |

### J. Protein Domains / PFAMs

| Current (source) | New name | Action |
|---|---|---|
| `protein_domains`, `protein_domains_description` (gene_mapping) | `pfam_ids`, `pfam_descriptions` | RENAME + MERGE |
| Protein `pfam_ids` | — | merge |
| EggNOG `PFAMs` | — | merge |
| Protein `domain_description` | `domain_description` | KEEP on Protein (longer text) |
| Protein `functional_motifs` | `functional_motifs` | KEEP on Protein |

### K. Carbohydrate-Active Enzymes (new)

| Current | New name | Action |
|---|---|---|
| EggNOG `CAZy` | `cazy_ids` | NEW |

### L. Protein-Level Properties (Protein node, some denormalized to Gene)

| Protein field | Protein node | Denormalize to Gene | Notes |
|---|---|---|---|
| `sequence_length`, `molecular_mass` | KEEP | no | Size metrics |
| `protein_name`, `protein_synonyms` | KEEP | no | UniProt name; already on Gene via gene_name |
| `amino_acid_sequence` | KEEP | no | |
| `prott5_embedding` | KEEP | no | |
| `catalytic_activity` | KEEP | **YES** | Key for function queries |
| `cofactors` | KEEP | no | |
| `subcellular_location` | KEEP | **YES** | Important for biology |
| `transmembrane_regions` | KEEP | **YES** | Structural classification |
| `signal_peptide` | KEEP | **YES** | Structural classification |
| `pathway_description` | KEEP | no | |
| `annotation_score` | KEEP | no | Quality indicator — stays on Protein only |
| `is_reviewed` | KEEP | no | SwissProt quality indicator — stays on Protein only |
| `caution_notes`, `interaction_notes` | KEEP | no | |
| `keywords`, `keyword_ids` | KEEP | no | UniProt keyword tags |
| `proteome_ids`, `refseq_ids`, `string_ids` | KEEP | no | Cross-references |
| `protein_family` | KEEP | **YES** | Useful for classification |

### M. Provenance / Source Tracking

Two complementary approaches:
- **Per-field source tags** for key canonical fields: `product_source`, `gene_name_source`, `function_description_source`, `go_terms_source` — simple strings ("cyanorak", "uniprot", "ncbi", "eggnog", "merged")
- **`annotation_quality`** score on Gene node: 0 = no annotation, 1 = eggnog-only, 2 = NCBI/Cyanorak, 3 = UniProt reviewed — quick filter for high-confidence genes

---

## Step 3 — Schema Update: `config/schema_config.yaml`

**Gene node — remove** (consolidated above):
`gene_names`, `gene_names_cyanorak`, `gene_synonym`, `kegg_ids`, `go_component`, `go_function`, `go_process`, `Ontology_term_ncbi`, `Ontology_term_cyanorak`, `eggNOG`, `eggNOG_description`, `kegg`, `kegg_description`, `protein_domains`, `protein_domains_description`, `go_biological_processes` (denormalized)

**Gene node — rename:**
`locus_tag_cyanoak` → `locus_tag_cyanorak`

**Gene node — add:**
`gene_name`, `gene_synonyms`, `llm_summary`, `cog_category`, `kegg_ko`, `kegg_ko_descriptions`, `kegg_pathway`, `kegg_module`, `kegg_reaction`, `kegg_brite`, `transporter_classification`, `bigg_reaction`, `cazy_ids`, `eggnog_ogs`, `eggnog_og_descriptions`, `seed_ortholog`, `max_annot_lvl`, `seed_ortholog_evalue`, `pfam_ids`, `pfam_descriptions`, `go_terms`, `go_term_descriptions`, `ontology_terms`, `annotation_quality`, `product_source`, `gene_name_source`, `catalytic_activity`, `subcellular_location`, `transmembrane_regions`, `signal_peptide`, `protein_family`

**Protein node** — largely unchanged except removing fields now also on Gene (avoids confusion for LLM queries).

---

## Step 4 — Adapter Update: `cyanorak_ncbi_adapter.py`

- `download_data()`: simplified to load-from-cache only; genome download extracted to Step 0 script
- `get_nodes()`: read `gene_annotations_merged.csv`; add new `GeneNodeField` enum members; populate all new canonical fields
- Keep backward compatibility: if merged CSV missing, fall back to gene_mapping.csv

---

## Step 5 — Tests

1. `tests/test_gene_annotations_merge.py` — unit tests for merge logic, priority rules, union dedup
2. `tests/kg_validity/test_structure.py` — add: MED4 COG category coverage > 80%; `gene_name` populated on >90% of genes
3. Regenerate `snapshot_data.json` after rebuild

---

## Critical Files

| File | Action |
|---|---|
| `multiomics_kg/download/download_genome_data.py` | **CREATE** (genome + UniProt download pipeline) |
| `multiomics_kg/download/build_gene_annotations.py` | **CREATE** (produces `gene_annotations_wide.json` + `gene_annotations_merged.json` per strain) |
| `config/schema_config.yaml` | **UPDATE** (schema cleanup per mapping table) |
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | **UPDATE** (load from cache only; read `gene_annotations_merged.json`) |
| `multiomics_kg/adapters/uniprot_adapter.py` | **UPDATE** (accept cache_dir; read/write taxid-keyed JSON) |
| `tests/test_gene_annotations_merge.py` | **CREATE** |
| `tests/kg_validity/test_structure.py` | **UPDATE** |

---

## Verification

```bash
# 1. Download pipeline (idempotent; skips already-cached)
uv run python multiomics_kg/download/download_genome_data.py

# 2. Build merged annotation tables for all strains
uv run python multiomics_kg/download/build_gene_annotations.py
# Optionally with LLM summaries (costs $):
uv run python multiomics_kg/download/build_gene_annotations.py --llm-summary

# 3. Inspect a few MED4 merged gene annotations
python -c "import json; d=json.load(open('cache/data/Prochlorococcus/genomes/MED4/gene_annotations_merged.json')); import pprint; pprint.pprint(list(d.items())[:3])"

# 4. Build KG (TEST_MODE=True for fast iteration)
uv run python create_knowledge_graph.py

# 5. Unit tests
pytest tests/test_gene_annotations_merge.py -v

# 6. Full build + KG validity
docker compose up -d
pytest tests/kg_validity/ -v
```
