---
name: check-gene-ids
description: Validate gene ID matching between paper CSV data and gene nodes in the knowledge graph. Reports per-paper match rates with full breakdowns (matched, RNA skipped, mismatched), cross-references Docker import report, and suggests fix strategies. Use to diagnose dangling edges.
argument-hint: [paper-name or "all"]
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *), Bash(docker *)
---

# Gene ID Validation Skill

Check whether gene IDs used in paper supplementary tables match existing gene nodes in the knowledge graph.

## Quick Start

```bash
# Check all papers against latest biocypher output (+ Docker import report if available)
uv run python .claude/skills/check-gene-ids/check_gene_ids.py

# Check a single paper
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"

# Specify biocypher output dir
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --biocypher-dir biocypher-out/20260205200505

# Use a specific import report file
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --import-report output/import.report

# Skip import report lookup
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --no-import-report
```

## What It Checks

For each publication and each supplementary table/analysis:

1. **Detects pre-resolved CSVs**: if `<stem>_resolved.csv` exists alongside the original file (created by `resolve_paper_ids.py` / prepare_data step 4), uses the `locus_tag` column from it — exactly what `omics_adapter` will use at build time
2. **Without resolved CSV**: reads the CSV file using `name_col` (and optional `secondary_name_col`) from the paperconfig, checks **ALL** gene IDs against Gene node `:ID` values, classifies each ID: **matched**, **matched via secondary column**, **RNA/tRNA/ncRNA** (expected skip), or **mismatched**
3. Cross-references mismatches against the **Docker import report** (if available) to confirm real dangling edges
4. Categorizes mismatched IDs by pattern (trailing asterisk, pseudo suffix, contig IDs, composite, etc.)
5. If not matched, determines **why** and suggests a **fix strategy**

## Key Features

### Full ID Checking
Checks ALL gene IDs in each CSV, not just a sample. This catches edge cases like Fang 2019 where 161/665 IDs had underscore mismatches that a 30-ID sample missed.

### Import Report Cross-Reference
Automatically loads the Neo4j import report to show which mismatches actually resulted in dangling edges:
- Tries Docker: `docker compose exec deploy cat /data/build2neo/import.report`
- Falls back to local: `output/import.report`
- Can be specified: `--import-report <path>`

### Secondary Gene ID Column
Papers may have a secondary gene identifier column (e.g., `locus_tag` alongside `gene_name`). Configure in paperconfig.yaml:
```yaml
statistical_analyses:
  - id: "my_analysis"
    name_col: "Gene"
    secondary_name_col: "locus_tag"  # fallback column
```
IDs are first checked via `name_col`. If not found, the corresponding `secondary_name_col` value is tried.

### RNA/tRNA/ncRNA Handling
IDs matching RNA patterns are counted separately and reported as "OK" — these features are intentionally not loaded as gene nodes. Patterns detected:
- `tRNA*`, `ncRNA*`, `rRNA*`, `Yfr*`, `tmRNA*`
- `RNA_<number>` (e.g., `RNA_48`)
- `PMT_ncRNA_*` (e.g., `PMT_ncRNA_Yfr7`)
- VIMSS tRNA/rRNA IDs (e.g., `A9601_tRNAAlaVIMSS1309073`)

### Pre-resolved CSV Detection
When `resolve_paper_ids.py` (prepare_data step 4) has run for a paper, it writes `<stem>_resolved.csv` alongside each original CSV. This file has a `locus_tag` column pre-computed from `gene_id_mapping.json`. The `omics_adapter` automatically uses it at build time, so check-gene-ids does too:

- **Rows with a locus_tag**: checked against gene nodes — should all match
- **Rows with NaN locus_tag** (`!Res` column): the ID could not be resolved; the omics_adapter will skip these rows (no dangling edges)
- Rows using the resolved CSV are marked with `*` in the summary table's Fix Strategy column

This means a paper can have non-standard IDs in its `name_col` (JGI IDs, probesets, etc.) **without** needing a paperconfig change — as long as `gene_id_mapping.json` covers those IDs and step 4 has been run.

### Mismatch Categorization
Mismatched IDs are grouped by pattern to help diagnose issues (only relevant when no `_resolved.csv` exists):
- `trailing_asterisk`: IDs ending with `*` (footnote artifacts)
- `trailing_whitespace`: IDs with trailing spaces
- `pseudo_suffix`: IDs with `_pseudo` suffix
- `unannotated`: literal `(unannotated)` entries
- `composite_ids`: IDs containing commas (multiple names in one cell)
- `contig_ids`: contig-based IDs (draft genome features)
- `other`: IDs that don't match known patterns

## Fix Strategies

| Strategy | Meaning | Action |
|----------|---------|--------|
| **MATCH** | IDs (or resolved locus tags) match gene nodes | No action needed |
| **RNA_ONLY** | All IDs are RNA features | No action needed |
| **PRE_UNRESOLVED** | `_resolved.csv` exists but all locus_tags are NaN | Re-run steps 3+4: `bash scripts/prepare_data.sh --steps 3 4 --force` |
| **REBUILD_RESOLVED** | `_resolved.csv` locus tags don't match gene nodes | Re-run steps 3+4 with `--force` to rebuild `gene_id_mapping.json` and `_resolved.csv` |
| **LOAD_ORGANISM** | Organism genome not loaded | Add genome to `cyanobacteria_genomes.csv` |
| **CHANGE_NAME_COL** | Another CSV column has matching IDs | Change `name_col` in paperconfig.yaml |
| **CREATE_MAPPING_CSV** | IDs match a column in `gene_mapping.csv` | Preferred: run steps 3+4 to build `_resolved.csv` automatically. Alternative: run `/fix-gene-ids` to create `_with_locus_tag.csv` and update paperconfig |
| **CREATE_MAPPING_SUPP** | IDs match `alt_id` in `gene_mapping_supp.csv` (cross-paper) | Preferred: run steps 3+4. Alternative: run `/fix-gene-ids` to resolve via supplementary mappings |
| **CREATE_MAPPING_GFF** | IDs found in GFF/GBK but NOT in `gene_mapping.csv` | Re-run `prepare_data.sh` to regenerate annotations, then steps 3+4 |
| **UNRELATED_IDS** | IDs don't match anything | Manual investigation needed |

### Preferred fix for non-locus_tag IDs (steps 3+4 pipeline)

When a paper uses non-standard gene IDs (JGI catalog IDs, probesets, numeric IDs, etc.), the preferred fix is:
1. Add `id_translation` and/or `annotation_gff` entries to the paperconfig if needed (for IDs not already in annotations)
2. Re-run step 3 to build/update `gene_id_mapping.json`:
   ```bash
   uv run python -m multiomics_kg.download.build_gene_id_mapping --strains <Strain> --force
   ```
3. Re-run step 4 to create/update `_resolved.csv`:
   ```bash
   uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Paper Name" --force
   ```
4. Re-run `/check-gene-ids` — the resolved CSV will be detected automatically and the paper should show MATCH

**No paperconfig change needed** — the original `name_col` and `filename` stay as-is; the omics_adapter uses `_resolved.csv` transparently.

## Annotation Validation

After gene IDs are confirmed to match, you can also validate that **functional annotations** in the paper CSV agree with the canonical annotations in `gene_annotations_merged.json`:

```bash
# All papers (token Jaccard ≥50% = match, no LLM)
uv run python scripts/validate_annotations.py --no-llm

# Single paper, with LLM mismatch analysis
uv run python scripts/validate_annotations.py --papers "Biller 2018"

# Options
uv run python scripts/validate_annotations.py --no-llm --papers "Name1" "Name2"
uv run python scripts/validate_annotations.py --llm-model gpt-4.1-nano
```

The script auto-detects annotation/product/description/definition columns in each CSV, maps gene IDs to locus tags (same logic as this skill), and reports per-column match rates, e.g. `(963/1794) 54% match`. It also optionally sends a batch of mismatches to an LLM for qualitative analysis (vocabulary differences, specificity, etc.).

Papers with no annotation columns (e.g. Al-Hosani 2015, Anjur 2025) are correctly reported as having no comparable data.

## Refreshing the Gene ID Mapping

This skill detects and uses `_resolved.csv` files (written by `resolve_paper_ids.py`, prepare_data step 4). When a `_resolved.csv` is present, it reflects the state of `gene_id_mapping.json` at the time it was generated. If match rates are unexpectedly low or a resolved CSV shows `PRE_UNRESOLVED`, rebuild:

```bash
# Rebuild gene_id_mapping.json for one strain (picks up new id_translation/annotation_gff entries)
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains MIT9301 --force

# Re-resolve paper CSVs for one paper (creates/updates _resolved.csv)
uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Anjur 2025" --force

# Or rebuild all strains + re-resolve all paper CSVs
bash scripts/prepare_data.sh --steps 3 4 --force
```

The fallback (when no `_resolved.csv` exists) uses `build_id_lookup()` from `gene_id_utils.py` directly, which also prefers `gene_id_mapping.json` over the older `gene_mapping_supp.csv`.

## Workflow

When invoked with a paper name (e.g., `/check-gene-ids "Biller 2018"`):

1. Find the paperconfig.yaml in `data/Prochlorococcus/papers_and_supp/$ARGUMENTS/`
2. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/$ARGUMENTS/paperconfig.yaml"`
3. Review the output — check whether a `_resolved.csv` was used (marked `*`), the Breakdown line, and mismatched ID categories
4. **If MATCH (with `*`)**: pre-resolution pipeline is working — no action needed
5. **If PRE_UNRESOLVED or REBUILD_RESOLVED**: re-run steps 3+4 with `--force` (see "Refreshing" above)
6. For CHANGE_NAME_COL: read the CSV to verify the suggested column, then update paperconfig.yaml
7. For CREATE_MAPPING_CSV / CREATE_MAPPING_SUPP: preferred fix is steps 3+4 (builds `_resolved.csv`); fallback is `/fix-gene-ids` which creates `_with_locus_tag.csv` and requires paperconfig update
8. For CREATE_MAPPING_GFF: re-run `prepare_data.sh` to regenerate annotations, then steps 3+4
9. For LOAD_ORGANISM: note the organism for future genome loading

When invoked with "all" or no arguments:

1. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py`
2. Review the summary table — columns: Match / 2nd (secondary col) / RNA / !Res (pre-unresolved) / Miss (mismatched) / IR (import report); `*` in Fix Strategy = using pre-resolved CSV
3. Prioritize fixes: PRE_UNRESOLVED/REBUILD_RESOLVED (run steps 3+4) > CHANGE_NAME_COL (easy) > CREATE_MAPPING_CSV (run steps 3+4 or `/fix-gene-ids`) > CREATE_MAPPING_GFF (needs annotation rebuild) > LOAD_ORGANISM (needs genome data)
