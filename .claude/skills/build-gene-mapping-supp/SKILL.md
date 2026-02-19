---
name: build-gene-mapping-supp
description: For each organism, collect alternative gene ID columns from all paper supplementary CSVs and write gene_mapping_supp.csv to the organism's cache directory. Used to build extended ID mappings for check-gene-ids and fix-gene-ids.
argument-hint: [organism-name or "all"]
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *)
---

# Build Gene Mapping Supplement Skill

Scans all supplementary CSV files across all papers and collects alternative gene identifiers
for each organism, writing a `gene_mapping_supp.csv` file per organism to its cache directory.

## Quick Start

```bash
# Build for all organisms
uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py

# Build for a specific organism
uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py --organism "Prochlorococcus MED4"

# Dry run (show what would be written without writing files)
uv run python .claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py --dry-run
```

## What It Does

For each organism (strain):

1. Finds all papers in `paperconfig_files.txt` that have analyses for that organism
2. Loads the organism's `gene_mapping.csv` to resolve gene IDs → locus_tags
3. For each paper's supplementary CSV:
   - Reads the CSV using `name_col` and `skip_rows` from the paperconfig
   - For each row, resolves the `name_col` value to a locus_tag via `gene_mapping.csv`
   - Collects values from other "ID-like" columns (skipping fold change, p-value, and description columns)
4. Writes `gene_mapping_supp.csv` to the organism's cache directory

Skips tRNA/ncRNA/rRNA gene IDs automatically.

## Output Format

`cache/data/<Organism>/genomes/<Strain>/gene_mapping_supp.csv`:

| Column | Description |
|--------|-------------|
| `locus_tag` | Primary locus tag resolved from gene_mapping.csv |
| `alt_id` | Alternative identifier found in the supplementary CSV |
| `source_col` | Column name in the supplementary CSV where `alt_id` came from |
| `source_csv` | Basename of the supplementary CSV file |
| `paper` | Paper name (from paperconfig `papername`) |

Example:
```
locus_tag,alt_id,source_col,source_csv,paper
PMT9312_0001,WP_011375566.1,protein_id,DE genes MIT9312 ....csv,Tetu 2019
PMT9312_0001,PMT9312_RS00005,MIT9312 (2017 NCBI),DE genes MIT9312 ....csv,Tetu 2019
```

## Column ID Detection

The skill auto-detects which columns contain gene identifiers using these rules:
- **Include if** column name contains: locus, tag, gene, protein, orf, id, accession, homolog
- **Include if** values are short strings (avg < 50 chars) and not purely numeric
- **Exclude** columns specified in paperconfig: `name_col`, `logfc_col`, `adjusted_p_value_col`
- **Exclude** purely numeric columns (coordinates, fold changes, p-values)
- **Exclude** long-text columns (descriptions, GO terms)
- **Exclude** organism/taxonomy columns: `species`, `organism`, `strain`, `taxon`
- **Exclude** annotation/regulation columns: names containing `annotation`, `regulation`, `description`, `product`, etc. (e.g. "Gene annotation" is excluded despite containing "gene")
- **Exclude** low-cardinality columns: ≤10 unique values where each repeats ≥5× on average (e.g. "Up"/"Down" in a Regulation column)

## When to Run

- After adding new paperconfig files (to pick up their alternative ID columns)
- When `/check-gene-ids` reports `UNRELATED_IDS` (the IDs might be in another paper's columns)
- Periodically to refresh the mappings

## Workflow

When invoked (e.g., `/build-gene-mapping-supp MED4`):

1. Run the script for that organism
2. Check the output: `head cache/data/Prochlorococcus/genomes/MED4/gene_mapping_supp.csv`
3. Use the file as an additional lookup with `/check-gene-ids` and `/fix-gene-ids`
