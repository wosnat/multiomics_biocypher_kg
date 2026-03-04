---
name: fix-gene-ids
description: Map gene IDs in paper CSVs to locus tags using gene_id_mapping.json (v2 three-tier mapping). Creates new CSV files with a locus_tag column. Supports alternative column scanning and import report mismatch recovery.
argument-hint: <paperconfig-path>
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *), Bash(docker *)
---

# Fix Gene IDs Skill

Map gene IDs in paper supplementary table CSVs to locus tags using the organism's `gene_id_mapping.json` (three-tier v2 mapping: Tier 1 specific lookup, Tier 2 protein-level multi-match, Tier 3 generic names). Creates new `_with_locus_tag.csv` files with a `locus_tag` column that matches gene node primary IDs.

> **Prefer the steps 3+4 pipeline for most cases.** Running `resolve_paper_ids.py` (prepare_data step 4) creates a `_resolved.csv` alongside the original file that the `omics_adapter` uses automatically — **no paperconfig change needed**. Use this skill (`fix-gene-ids`) only when the automated pipeline can't resolve IDs and you need manual multi-column scanning, or when you want a permanent `_with_locus_tag.csv` checked into the repo.
>
> Steps 3+4 pipeline (preferred):
> ```bash
> uv run python -m multiomics_kg.download.build_gene_id_mapping --strains <Strain> --force  # step 3
> uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Paper Name" --force   # step 4
> ```

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

### Phase 1: Primary mapping (via name_col + paperconfig id_columns)
1. Reads the paper CSV using `name_col` from the config
2. Loads the organism's `gene_id_mapping.json` (v2); falls back to `gene_annotations_merged.json` when not yet built
3. Uses `resolve_row()`: three-tier strategy — Tier 1 specific_lookup (direct + heuristics: zero-pad, strip asterisk), Tier 2+3 multi_lookup singletons. Also applies `id_columns` from paperconfig as fallback columns in the same pass.

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

> **Alternative (no paperconfig change):** if the IDs in the original CSV can be resolved via `gene_id_mapping.json`, run steps 3+4 instead — this creates `_resolved.csv` transparently and the omics_adapter uses it automatically without touching the paperconfig.

Optionally, validate that functional annotations in the CSV match the reference annotations:

```bash
uv run python scripts/validate_annotations.py --papers "Author Year" --no-llm
# or with LLM mismatch analysis:
uv run python scripts/validate_annotations.py --papers "Author Year"
```

This reports per-column match rates (e.g. `(963/1794) 54% match`) comparing annotation/product/description columns in the CSV against `gene_annotations_merged.json`. Useful to confirm that a newly added paper's data aligns with the reference genome annotation.

## Report

Each run generates `fix_gene_ids_report.md` in the paper's directory with:
- Per-analysis primary mapping stats
- Alternative column recovery details (which column, which value, which locus_tag)
- Still-unmapped IDs with their values across all candidate columns

## When to Use

Use this skill when the **steps 3+4 pipeline is insufficient**:

- `/check-gene-ids` reports `CREATE_MAPPING_CSV` or `CREATE_MAPPING_SUPP` and you need the mapping to be permanent (checked-in `_with_locus_tag.csv` rather than a generated `_resolved.csv`)
- A paper's CSV has partial matches and you need multi-column scanning to recover more IDs (tries all ID-like columns automatically)
- Patching existing `_with_locus_tag.csv` files to fill empty locus_tag values or fix import-report mismatches (`--patch` mode)
- The paper's IDs can't be resolved via `gene_id_mapping.json` and require custom lookup logic

**If the paper uses non-standard IDs** (JGI IDs, probesets, paper-specific annotation IDs), first add `id_translation` and/or `annotation_gff` entries to the paperconfig and rebuild, then decide: steps 3+4 (transparent `_resolved.csv`, no paperconfig change) or this skill (`_with_locus_tag.csv`, requires paperconfig update):

```bash
# Option A: transparent pipeline (preferred)
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains <Strain> --force
uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Paper Name" --force

# Option B: permanent mapping file (this skill)
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains <Strain> --force
uv run python .claude/skills/fix-gene-ids/fix_gene_ids.py --paperconfig "data/.../paperconfig.yaml"
# then update filename + name_col in paperconfig.yaml
```
