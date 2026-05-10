---
name: tcdb-diamond
description: Run diamond blastp vs. the curated TCDB FASTA per strain to generate per-protein TCDB classifications, with tiered confidence (5-part / 4-part / 3-part) plus consensus + eggNOG-agreement tags. Phase 1 — produces inspectable `<strain>.tcdb.calls.json` artifacts; KG integration deferred to Phase 2.
argument-hint: [--strain <name> | --force | --refresh-tcdb | --threads <n>]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(perl *), Bash(diamond *), Bash(git *)
---

# TCDB-Diamond Skill

Per-strain `diamond blastp` against the curated TCDB FASTA. Augments eggNOG's
family-level (3-part) `KEGG_TC` annotations with sequence-similarity-based calls
that can reach `tc_specificity` (5-part) when identity warrants it.

## Quick Start

```bash
# Run all genome strains (skip already-done)
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py

# Run a single strain
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MED4

# Force re-run even if calls.json exists
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --force

# Re-download TCDB FASTA + diamond DB
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --refresh-tcdb

# Use more threads
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --threads 8

# Smoke test: first 100 proteins of one strain (~10-30s)
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MIT1002 --limit 100
```

## One-Time Setup

```bash
# Saier Lab's TCDB tools (used to download the FASTA + build the diamond DB)
git clone https://github.com/SaierLaboratory/TCDBtools.git ~/tools/TCDB/TCDBtools
```

Optional `.env` entry to relocate the TCDB data dir (default: `~/tools/TCDB`):

```
TCDB_DATA_DIR=/path/to/TCDB
```

System tools required: `perl`, `wget`, `diamond` (already required by other parts of the project).

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`, filters to `organism_type='genome_strain'`.
2. Ensures `~/tools/TCDB/DB/tcdb.dmnd` exists, building it via `extractTCDB.pl -i tcdb -f diamond` if missing or `--refresh-tcdb`.
3. For each strain, runs `diamond blastp` of `<data_dir>/protein.faa` against the TCDB diamond DB, writing raw 8-column TSV to `<data_dir>/tcdb/<strain>.tcdb.tsv`.
4. Applies the per-hit tier policy + per-protein consensus collapse (Spec §6.3, §6.4-A).
5. Joins each protein's eggNOG `KEGG_TC` (from `<data_dir>/eggnog/<strain>.emapper.annotations`) and computes the agreement tag (Spec §6.4-B).
6. Tags class-9 calls (Spec §6.4-C).
7. Writes `<data_dir>/tcdb/<strain>.tcdb.calls.json` (per-protein records) and `<strain>.tcdb.skill_summary.json` (per-strain stats).
8. Prints a status table to stdout. Diamond's full per-strain stdout+stderr are captured to `logs/tcdb_diamond_<strain>.log` (auto-gitignored via `*.log`).

## Tier Policy

| Tier | Truncate to | Identity | Coverage rule | Notes |
|---|---|---|---|---|
| 1 | 5 parts (`tc_specificity`) | ≥ 70% | qcov ≥ 70% | High-confidence specificity call |
| 2 | 4 parts (`tc_subfamily`) | ≥ 40% | qcov ≥ 60% | Solid subfamily call |
| 3 | 3 parts (`tc_family`) | (no floor) | qcov ≥ 40% **OR** scov ≥ 40% | gblast3-style floor |

All tiers also require: e-value ≤ 0.001, HSP length ≥ 50 aa.

## Output Schema (`<strain>.tcdb.calls.json`)

Keyed by NCBI protein_id (WP_ accession):

```json
{
  "WP_011131900.1": {
    "tcid": "1.A.11.1.5",
    "level_kind": "tc_specificity",
    "tier": 1,
    "identity": 87.4,
    "qcov": 92.1,
    "scov": 89.7,
    "evalue": 1.2e-180,
    "length": 412,
    "consensus_n": 5,
    "consensus_agreement": "5_part",
    "egn_agreement": "refines",
    "egn_tcid": "1.A.11",
    "incompletely_characterized": false
  }
}
```

`egn_agreement` values: `confirms` (same TC) | `refines` (diamond deeper than eggNOG, same lineage) | `extends` (eggNOG had no TC) | `conflicts` (different family).

## Phase 2 (Future)

Phase 1 artifacts sit in the strain cache for inspection. Phase 2 (separate spec) will integrate them into `gene_annotations_merged.json` and the KG. See [docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md §7](../../../docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md).

## Workflow When Invoked

1. Verify one-time setup is done (`~/tools/TCDB/TCDBtools/bin/extractTCDB.pl` exists).
2. Run the skill: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py`.
3. Review the status table — note any FAILED strains.
4. Inspect `<data_dir>/tcdb/<strain>.tcdb.calls.json` for spot checks.
