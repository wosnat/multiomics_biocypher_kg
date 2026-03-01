# Project Memory

## scripts/validate_annotations.py
- Validates functional annotations in paper CSVs against `gene_annotations_merged.json`
- Reuses `multiomics_kg/utils/gene_id_utils.py` for gene ID mapping
- Auto-detects annotation columns (product/annotation/description/definition/genbank/rast/role/function/note keywords)
- Skips gene-symbol columns (short values, avg len <12, no spaces)
- Token Jaccard overlap ≥0.5 = broad match; also reports exact match count
- LLM batch analysis (≤30 mismatch pairs per column) via LangChain init_chat_model
- Deduplicates CSVs within a paper (multiple analyses can share one CSV)
- skip_rows: uses same logic as omics_adapter (`pd.read_csv(..., skiprows=N)`)
- Usage: `uv run python scripts/validate_annotations.py [--no-llm] [--papers "Name"] [--llm-model MODEL]`
- Al-Hosani 2015 and Anjur 2025 have no annotation columns (by design — CSVs only have fold-change data)

## build_gene_id_mapping.py
- Location: `multiomics_kg/download/build_gene_id_mapping.py`
- Run as module: `uv run python -m multiomics_kg.download.build_gene_id_mapping [--strains X] [--force]`
- Phase 3 of gene_id_mapping_utility plan (implemented)
- Outputs `gene_id_mapping.json` and backward-compat `gene_mapping_supp.csv` per strain
- For MIT9301/Anjur 2025: resolves 1872/1939 genes from annotation CSV, 235/242 DE rows (97.1%)
- 7 unresolved DE rows: genuinely unannotated genes (all ID columns NaN in annotation CSV)
- Processing order: id_translation first → rebuild lookup → annotation_gff → csv tables

## Key patterns

### Organism → genome dir
- Defined in `gene_id_utils.ORGANISM_TO_GENOME_DIR`
- `get_genome_dir(organism_name, project_root)` does fuzzy matching

### gene_annotations_merged.json fields for annotation comparison
- `product` (NCBI) and `product_cyanorak` (Cyanorak) are the reference fields

### paperconfig_files.txt
- `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` lists all active papers

### skip_rows in paperconfigs
- Defined at the supp_table level (not per-analysis)
- Uses pandas `skiprows=N` (skip first N rows, next row is header)
