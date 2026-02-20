---
name: fix-gene-ids
description: Map gene IDs in paper CSVs to locus tags using gene_mapping.csv and gene_mapping_supp.csv. Creates new CSV files with a locus_tag column. Supports alternative column scanning and import report mismatch recovery.
argument-hint: <paperconfig-path>
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *), Bash(docker *)
---

# Fix Gene IDs Skill

Map gene IDs in paper supplementary table CSVs to locus tags using the organism's `gene_mapping.csv` and `gene_mapping_supp.csv` (cross-paper alternative IDs). Creates new `_with_locus_tag.csv` files with a `locus_tag` column that matches gene node primary IDs.

## Quick Start

```bash
# Fix a single paper (with alt-column scanning + import report)
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml"

# Dry run (show mapping stats without writing files)
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --dry-run

# With explicit import report for mismatch detection
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --import-report output/import.report

# Patch existing _with_locus_tag.csv files (only fix empty/mismatched rows)
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml" --patch --import-report output/import.report
```

## What It Does

For each analysis in a paperconfig:

### Phase 1: Primary mapping (via name_col)
1. Reads the paper CSV using `name_col` from the config
2. Loads the organism's `gene_mapping.csv` and `gene_mapping_supp.csv`
3. Maps each gene ID via: direct match, gene_names lookup, old_locus_tags, protein_id, zero-padding, composite splitting

### Phase 2: Alternative column recovery
4. For rows that are **unmapped** or where the mapped locus_tag is in the **import report mismatch set** (dangling edges):
   - Scans other CSV columns that look like gene identifiers (auto-detected)
   - Tries mapping each alternative value through the same lookup machinery
   - Uses the first successful mapping that produces a non-mismatched locus_tag
   - Tracks which alternative columns resolved each row

### Phase 3: Output
5. Writes `<original_name>_with_locus_tag.csv` (or overwrites in patch mode)
6. Writes `fix_gene_ids_report.md` in the paper's directory with full details

## Flags

| Flag | Description |
|------|-------------|
| `--paperconfig PATH` | Path to paperconfig.yaml (required) |
| `--dry-run` | Show stats without writing files |
| `--import-report PATH` | Path to import.report file (default: tries output/import.report, then Docker) |
| `--no-import-report` | Skip loading import report |
| `--patch` | Update existing `_with_locus_tag.csv` in place (only fix empty/mismatched rows) |

## After Running

Update the paperconfig.yaml:
- Change `filename` to point to the new `_with_locus_tag.csv`
- Change `name_col` to `"locus_tag"`

Then run `/check-gene-ids` to verify the fix.

## Report

Each run generates `fix_gene_ids_report.md` in the paper's directory with:
- Per-analysis primary mapping stats
- Alternative column recovery details (which column, which value, which locus_tag)
- Still-unmapped IDs with their values across all candidate columns

## When to Use

- When `/check-gene-ids` reports `CREATE_MAPPING_CSV` or `CREATE_MAPPING_SUPP` fix strategy
- When a paper's CSV has partial matches and you want to recover more IDs via alternative columns
- When patching existing `_with_locus_tag.csv` files to fill empty locus_tag values or fix import-report mismatches
