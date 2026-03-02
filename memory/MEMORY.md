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

## gene_id_mapping_utility.md plan status
- Phases 1-7: COMPLETE (build scripts, gene_id_utils extensions, omics_adapter, prepare_data.sh, paperconfig skill)
- Phase 8: COMPLETE (check/fix-gene-ids SKILL.md updated with mapping refresh note)
- MIT9313_resources paperconfig: CREATED — `data/Prochlorococcus/papers_and_supp/MIT9313_resources/paperconfig.yaml`
  - `id_translation` for MIT9313_genbank.tsv: locus_tag, PMTid (old_locus_tag), P9313name (old_locus_tag), uniprot_id
- barreto 2022: `annotation_gtf_ez55` entry added for `EZ55.exon.fixed2.gtf` (organism: Alteromonas macleodii EZ55)
- Phase 2 curation (20 remaining paperconfigs): PENDING — add id_columns/product_columns to each

## build_gene_id_mapping.py
- Location: `multiomics_kg/download/build_gene_id_mapping.py`
- Run as module: `uv run python -m multiomics_kg.download.build_gene_id_mapping [--strains X] [--force]`
- Phase 3 of gene_id_mapping_utility plan (implemented)
- Outputs `gene_id_mapping.json` and backward-compat `gene_mapping_supp.csv` per strain
- For MIT9301/Anjur 2025: resolves 1872/1939 genes from annotation CSV, 235/242 DE rows (97.1%)
- 7 unresolved DE rows: genuinely unannotated genes (all ID columns NaN in annotation CSV)
- Processing order: id_translation first → rebuild lookup → annotation_gff → csv tables

## resolve_paper_ids.py (Phase 6, prepare_data step 4)
- Location: `multiomics_kg/download/resolve_paper_ids.py`
- Run as module: `uv run python -m multiomics_kg.download.resolve_paper_ids [--force] [--papers "Name"]`
- `--papers` filter matches against papername OR directory name (handles papername≠dirname cases)
- For each csv table with non-locus_tag name_col: resolves via gene_id_mapping.json
- Writes `<stem>_resolved.csv` with `locus_tag` + `resolution_method` columns added
- Skips if _resolved.csv is newer than source (unless --force)
- Report: per publication / organism, shows unresolved row IDs
- Anjur 2025 result: 235/242 (97.1%), 7 unresolved (genuinely unannotated)

## omics_adapter.py (Phase 6 change)
- `_load_and_create_edges(filename, analysis, sep=',')` now has `sep` param
- Probes for `<stem>_resolved.csv` next to filename; uses it when present
- Uses `locus_tag` column from resolved CSV; NaN rows → skipped (no dangling edges)
- Falls back to original `name_col` when no resolved CSV (old behaviour)
- `get_edges()` passes sep + table-level organism into each stat_analysis

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
