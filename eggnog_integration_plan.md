# Plan: EggNOG Mapper Integration

## Context

Several functional annotation categories are absent or incomplete in the current KG, especially for **Alteromonas** strains (which have no Cyanorak annotations): COG categories, KEGG KO IDs, cross-strain OG identifiers, CAZy enzymes, and extended GO/EC coverage.

All 13 strain protein FASTA files are already cached at `cache/data/<Organism>/genomes/<Strain>/protein.faa`. The integration strategy is to **extend `cyanorak_ncbi_adapter.py`** to merge EggNOG TSV data into gene node properties during `download_data()` — same pattern as the existing EC number merge from GBFF files.

**Approach**: Install and run EggNOG mapper on a test strain first, inspect the real output, then develop the integration.

---

## Phase 1: Install EggNOG Mapper and Run Test

### Install
```bash
conda install -c bioconda eggnog-mapper
# OR
pip install eggnog-mapper
```

### Download databases (one-time, ~50 GB)
```bash
download_eggnog_data.py -y --data_dir /path/to/eggnog_db
```

### Test run on MED4
```bash
emapper.py \
  -i cache/data/Prochlorococcus/genomes/MED4/protein.faa \
  --output eggnog \
  --output_dir cache/data/Prochlorococcus/genomes/MED4/ \
  --cpu 8 \
  -m diamond
```

Output file: `cache/data/Prochlorococcus/genomes/MED4/eggnog.emapper.annotations`

**Inspect output before writing any code.** Verify:
1. What the `#query` column looks like (full FASTA header vs. just accession)
2. Which columns have good coverage for cyanobacteria
3. Coverage of CAZy, PFAMs, KEGG fields

---

## Phase 2: Preprocessing Script — `scripts/run_eggnog_mapper.py`

After inspecting the output, write a script that runs emapper on all 13 strains:
- Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` to get all organisms/strains
- Infers protein FASTA path: `cache/data/<Organism>/genomes/<Strain>/protein.faa`
- Output file: `cache/data/<Organism>/genomes/<Strain>/eggnog.emapper.annotations`
- **Skips** if output file already exists (caching)
- Logs per-strain progress and handles subprocess errors

---

## Phase 3: Adapter Changes — `cyanorak_ncbi_adapter.py`

**File:** [multiomics_kg/adapters/cyanorak_ncbi_adapter.py](multiomics_kg/adapters/cyanorak_ncbi_adapter.py)

### EggNOG mapper output columns (21 total)

| Column | Useful? | → Gene property |
|--------|---------|----------------|
| `#query` | join key | parse WP_xxx accession from header |
| `seed_ortholog` | maybe | skip (internal DB ref) |
| `evalue`, `score` | no | skip |
| `eggNOG_OGs` | **yes** | `eggnog_ogs` (str[], split on `,`) |
| `max_annot_lvl` | no | skip |
| `COG_category` | **yes** | `cog_category` (str[], split letter-by-letter) |
| `Description` | **yes** | update `eggNOG_description` where currently empty |
| `Preferred_name` | **yes** | `eggnog_preferred_name` |
| `GOs` | **yes** | `go_terms_emapper` (str[], raw — no deduplication) |
| `EC` | **yes** | merge into existing `ec_numbers` |
| `KEGG_ko` | **yes** | `kegg_ko_ids` (str[], strip `ko:` prefix) |
| `KEGG_Pathway` | **yes** | `kegg_pathway_ids` (str[]) |
| `KEGG_Module` | maybe | skip for now |
| `KEGG_Reaction` | maybe | skip for now |
| `KEGG_rclass` | no | skip |
| `BRITE` | maybe | skip for now |
| `KEGG_TC` | no | skip |
| `CAZy` | **yes** | `cazy_ids` (str[]) |
| `BiGG_Reaction` | no | skip |
| `PFAMs` | **yes** | `pfam_ids` (str[], complements existing `protein_domains`) |

### ID mapping

FASTA headers: `>WP_011306043.1 photosystem II protein D2 [Prochlorococcus...]`
EggNOG `#query` will likely be the full header or just the accession — verify in Phase 1.
Join: `query.split()[0].split('.')[0]` (or similar) → `protein_id` column in `gene_mapping.csv` → gene node.

### New `GeneNodeField` enum entries
```python
EGGNOG_OGS = "eggnog_ogs"
COG_CATEGORY = "cog_category"
KEGG_KO_IDS = "kegg_ko_ids"
KEGG_PATHWAY_IDS = "kegg_pathway_ids"
EGGNOG_PREFERRED_NAME = "eggnog_preferred_name"
GO_TERMS_EMAPPER = "go_terms_emapper"
CAZY_IDS = "cazy_ids"
PFAM_IDS = "pfam_ids"
```

### New `_split_field()` entries
Treat `eggnog_ogs`, `kegg_ko_ids`, `kegg_pathway_ids`, `cazy_ids`, `pfam_ids`, `go_terms_emapper` as arrays (split on `,`).
`cog_category`: split each character into an array element (e.g. `"KL"` → `["K","L"]`).

### New `_load_eggnog_annotations()` method
1. Check for `eggnog.emapper.annotations` in strain cache dir
2. Read TSV, skip `##` comment lines, parse relevant columns
3. Extract protein accession from `#query` (handle whether full header or accession-only)
4. Return DataFrame indexed by `protein_id`

### Extend `download_data()` (after EC number merge, ~line 428)
```python
eggnog_df = self._load_eggnog_annotations()
if eggnog_df is not None:
    self.data_df = self.data_df.merge(eggnog_df, on='protein_id', how='left')
    # Update eggNOG_description where currently empty
    mask = self.data_df['eggNOG_description'].isna()
    self.data_df.loc[mask, 'eggNOG_description'] = self.data_df.loc[mask, 'eggnog_description_emapper']
    # Merge EC: union of existing + emapper
    self.data_df['ec_numbers'] = (
        self.data_df['ec_numbers'].fillna('') + ',' +
        self.data_df['ec_emapper'].fillna('')
    ).str.strip(',').replace('', None)
```

---

## Phase 4: Schema — `config/schema_config.yaml`

**File:** [config/schema_config.yaml](config/schema_config.yaml)

Add to Gene node properties (after `eggNOG_description`):
```yaml
eggnog_ogs: str[]
cog_category: str[]
kegg_ko_ids: str[]
kegg_pathway_ids: str[]
eggnog_preferred_name: str
go_terms_emapper: str[]
cazy_ids: str[]
pfam_ids: str[]
```

---

## Phase 5: Pipeline Flag — `create_knowledge_graph.py`

**File:** [create_knowledge_graph.py](create_knowledge_graph.py)

```python
USE_EGGNOG = True   # Load EggNOG mapper annotations from cache if present
```

Pass to `CyanorakNcbi` / `MultiCyanorakNcbi` constructor; adapter loads TSV only when `True`.

---

## Critical Files

| File | Change type |
|------|-------------|
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | GeneNodeField enum + `_load_eggnog_annotations()` + merge in `download_data()` |
| `config/schema_config.yaml` | 8 new Gene node properties |
| `create_knowledge_graph.py` | `USE_EGGNOG` flag |
| `scripts/run_eggnog_mapper.py` | NEW preprocessing script |

---

## Verification

```bash
# After Phase 1+2: check output exists
ls cache/data/Prochlorococcus/genomes/MED4/eggnog.emapper.annotations

# Build graph in test mode
# TEST_MODE = True, CACHE = True, USE_EGGNOG = True
uv run python create_knowledge_graph.py

# Run unit tests
pytest -m "not slow and not kg" -v

# Check new columns appear in output CSVs
head -1 biocypher-out/*/Gene-*.csv | tr ',' '\n' | grep -E 'cog|kegg_ko|eggnog|cazy|pfam'

# After Docker deploy
# MATCH (g:Gene) WHERE g.cog_category IS NOT NULL RETURN g.locus_tag, g.cog_category LIMIT 10
```
