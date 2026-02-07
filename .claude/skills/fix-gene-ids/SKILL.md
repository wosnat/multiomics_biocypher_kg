---
name: fix-gene-ids
description: Map gene IDs in paper CSVs to locus tags using gene_mapping.csv. Creates new CSV files with a locus_tag column. Use when /check-gene-ids reports CREATE_MAPPING_CSV fix strategy.
argument-hint: <paperconfig-path>
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *)
---

# Fix Gene IDs Skill

Map gene IDs in paper supplementary table CSVs to locus tags using the organism's `gene_mapping.csv`. Creates new `_with_locus_tag.csv` files with a `locus_tag` column that matches gene node primary IDs.

## Quick Start

```bash
# Fix a single paper
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/he 2022/paperconfig.yaml"

# Dry run (show mapping stats without writing files)
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --dry-run
```

## What It Does

For each analysis in a paperconfig:

1. Reads the paper CSV using `name_col` from the config
2. Loads the organism's `gene_mapping.csv` from the cache directory
3. Builds lookup dicts from columns: `locus_tag`, `gene_names`, `gene`, `gene_names_cyanorak`, `old_locus_tags`
4. For each gene ID in the paper CSV:
   - If it's already a valid locus_tag -> keeps it
   - If not, tries splitting by comma and looking up each part in gene_names/gene columns
   - Adds a `locus_tag` column with the mapped value (empty if unmapped)
5. Writes `<original_name>_with_locus_tag.csv`
6. Reports mapping stats

## After Running

Update the paperconfig.yaml:
- Change `filename` to point to the new `_with_locus_tag.csv`
- Change `name_col` to `"locus_tag"`

Then run `/check-gene-ids` to verify the fix.

## When to Use

When `/check-gene-ids` reports `CREATE_MAPPING_CSV` fix strategy, meaning the paper's gene IDs (gene names, old locus tags, etc.) can be mapped to locus tags via `gene_mapping.csv`.
