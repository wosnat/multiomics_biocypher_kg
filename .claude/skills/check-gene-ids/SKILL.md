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

1. Reads the CSV file using `name_col` (and optional `secondary_name_col`) from the paperconfig
2. Checks **ALL** gene IDs (not just a sample) against Gene node `:ID` values
3. Classifies each ID: **matched**, **matched via secondary column**, **RNA/tRNA/ncRNA** (expected skip), or **mismatched**
4. Cross-references mismatches against the **Docker import report** (if available) to confirm real dangling edges
5. Categorizes mismatched IDs by pattern (trailing asterisk, pseudo suffix, contig IDs, composite, etc.)
6. If not matched, determines **why** and suggests a **fix strategy**

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

### Mismatch Categorization
Mismatched IDs are grouped by pattern to help diagnose issues:
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
| **MATCH** | IDs match gene nodes | No action needed |
| **RNA_ONLY** | All IDs are RNA features | No action needed |
| **LOAD_ORGANISM** | Organism genome not loaded | Add genome to `cyanobacteria_genomes.csv` |
| **CHANGE_NAME_COL** | Another CSV column has matching IDs | Change `name_col` in paperconfig.yaml |
| **CREATE_MAPPING_CSV** | IDs match a column in `gene_mapping.csv` | Run `/fix-gene-ids` to create `_with_locus_tag.csv`, then update paperconfig |
| **CREATE_MAPPING_SUPP** | IDs match `alt_id` in `gene_mapping_supp.csv` (cross-paper) | Run `/fix-gene-ids` to resolve via supplementary mappings |
| **CREATE_MAPPING_GFF** | IDs found in GFF/GBK but NOT in `gene_mapping.csv` | Add the GFF attribute as a new column in gene_mapping.csv |
| **UNRELATED_IDS** | IDs don't match anything | Manual investigation needed |

## Workflow

When invoked with a paper name (e.g., `/check-gene-ids "Biller 2018"`):

1. Find the paperconfig.yaml in `data/Prochlorococcus/papers_and_supp/$ARGUMENTS/`
2. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/$ARGUMENTS/paperconfig.yaml"`
3. Review the output — check the Breakdown line and mismatched ID categories
4. For CHANGE_NAME_COL: read the CSV to verify the suggested column, then update paperconfig.yaml
5. For CREATE_MAPPING_CSV: run `/fix-gene-ids` on the paperconfig, then update filename and name_col
6. For CREATE_MAPPING_SUPP: run `/fix-gene-ids` (it uses gene_mapping_supp.csv automatically), then update filename and name_col
7. For CREATE_MAPPING_GFF: extract the GFF attribute and add it to gene_mapping.csv
7. For LOAD_ORGANISM: note the organism for future genome loading

When invoked with "all" or no arguments:

1. Run: `uv run python .claude/skills/check-gene-ids/check_gene_ids.py`
2. Review the summary table — columns show Match/2nd/RNA/Miss/IR counts per analysis
3. Prioritize fixes: CHANGE_NAME_COL (easy) > CREATE_MAPPING_CSV (moderate) > CREATE_MAPPING_GFF (needs GFF parsing) > LOAD_ORGANISM (needs genome data)
