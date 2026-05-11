---
name: tcdb-diamond
description: Run diamond blastp vs. the curated TCDB FASTA per strain to generate per-protein TCDB classifications, with tiered confidence (5-part / 4-part / 3-part) plus consensus + eggNOG-agreement tags. Phase 1 — produces inspectable `<strain>.tcdb.calls.json` artifacts; KG integration deferred to Phase 2.
argument-hint: "[--strain <name> | --force | --refresh-tcdb | --threads <n>]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(diamond *)
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

No manual setup required. The skill downloads two TCDB-curated resources on first run:

1. **TCDB FASTA** from `https://tcdb.org/public/tcdb` → rewritten + built into a diamond DB at `~/tools/TCDB/DB/tcdb.dmnd`.
2. **TCDB Pfam→TC map** from `https://www.tcdb.org/cgi-bin/projectv/public/pfam.py` → cleaned 3-column TSV at `~/tools/TCDB/DB/tcdb_pfam_map.tsv` (~600 KB, ~8.3K Pfam→TC pairs across ~1.3K Pfams). Used as a tiebreaker signal when diamond and eggNOG disagree.

Both are re-downloaded on `--refresh-tcdb`.

Optional `.env` entry to relocate the TCDB data dir (default: `~/tools/TCDB`):

```
TCDB_DATA_DIR=/path/to/TCDB
```

System tool required: `diamond` (already required by other parts of the project).

**Note:** Saier Lab's `extractTCDB.pl` performs the same FASTA download + header rewrite, but contains a bash-only redirect (`>&/dev/null`) that fails on Ubuntu's default dash-based `/bin/sh`. The skill replicates the operation in Python to sidestep the portability issue, so no `git clone` of TCDBtools is needed.

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`, filters to `organism_type='genome_strain'`.
2. Ensures `~/tools/TCDB/DB/tcdb.dmnd` and `~/tools/TCDB/DB/tcdb_pfam_map.tsv` exist, downloading them from TCDB if missing or `--refresh-tcdb`.
3. For each strain, runs `diamond blastp` of `<data_dir>/protein.faa` against the TCDB diamond DB, writing raw 8-column TSV to `<data_dir>/tcdb/<strain>.tcdb.tsv`.
4. Applies the per-hit tier policy + per-protein consensus collapse + best-tier promotion (Spec §6.3, §6.4-A).
5. Joins each protein's eggNOG `KEGG_TC` (from `<data_dir>/eggnog/<strain>.emapper.annotations`) and computes the agreement tag (Spec §6.4-B).
6. Joins each protein's Pfam annotations (from `<data_dir>/gene_annotations_merged.json`) against the TCDB Pfam→TC map and computes a per-protein `pfam_agreement` tag — an independent third signal alongside diamond and eggNOG (Spec §6.4-D).
7. Tags class-9 calls (Spec §6.4-C).
8. Writes `<data_dir>/tcdb/<strain>.tcdb.calls.json` (per-protein records, including `confidence_score` and the `pfam_*` fields) and `<strain>.tcdb.skill_summary.json` (per-strain stats — now includes `pfam_agreement_distribution`).
9. Prints a status table to stdout. Diamond's full per-strain stdout+stderr are captured to `logs/tcdb_diamond_<strain>.log` (auto-gitignored via `*.log`).

## Tier Policy

| Tier | Truncate to | Identity | Coverage rule | Notes |
|---|---|---|---|---|
| 1 | 5 parts (`tc_specificity`) | ≥ 70% | qcov ≥ 70% | High-confidence specificity call |
| 2 | 4 parts (`tc_subfamily`) | ≥ 40% | qcov ≥ 60% | Solid subfamily call |
| 3 | 3 parts (`tc_family`) | (no floor) | qcov ≥ 40% **OR** scov ≥ 40% | gblast3-style floor |

All tiers also require: e-value ≤ 0.001, HSP length ≥ 50 aa.

## Output Schema (`<strain>.tcdb.calls.json`)

Keyed by NCBI protein_id (WP_ accession). Each protein has shared fields plus a `calls` list — one candidate per 3-part TC family the protein hits, sorted by `confidence_score` descending. Multi-domain proteins (RND + MFP + OMF efflux components, etc.) produce multiple candidates instead of being rejected.

```json
{
  "WP_010951455.1": {
    "egn_tcids": ["8.A.1.2.1"],
    "pfam_ids":  ["PF_HlyD"],
    "pfam_tc_families": ["8.A.1"],
    "calls": [
      {
        "tcid": "2.A.6.1",
        "level_kind": "tc_subfamily",
        "tier": 1,
        "confidence_score": 0.85,
        "identity": 100.0,
        "qcov": 100.0,
        "scov": 100.0,
        "evalue": 1.2e-200,
        "length": 400,
        "consensus_n": 2,
        "consensus_agreement": "4_part",
        "egn_agreement": "conflicts",
        "pfam_agreement": "confirms_eggnog",
        "incompletely_characterized": false
      },
      {
        "tcid": "8.A.1.2",
        "level_kind": "tc_subfamily",
        "tier": 2,
        "confidence_score": 0.55,
        "...": "...",
        "egn_agreement": "confirms",
        "pfam_agreement": "confirms_both"
      }
    ]
  }
}
```

Field semantics:
- **Protein-level** (shared across all candidates): `egn_tcids`, `pfam_ids`, `pfam_tc_families` — facts about the gene, not about a specific call.
- **Per-candidate**:
  - `tier` (1/2/3) and `level_kind` come from `effective_tier = max(best_hit_tier, depth_tier)` within this family's hits — the more conservative of the strongest hit's identity-tier and the consensus depth's tier. `tcid` is truncated to the parts justified by `effective_tier`.
  - `confidence_score` ∈ [0,1] = `(identity / 100) × (qcov / 100) × agreement_weight`, where agreement_weight is 1.0 / 0.85 / 0.7 for 5/4/3-part consensus. The list is sorted by this descending — `calls[0]` is the strongest call.
  - `identity`, `qcov`, `scov`, `evalue`, `length` come from the **best (highest-identity)** hit in THIS family group.
  - `egn_agreement` / `pfam_agreement` are computed PER candidate — different candidates of the same protein can have different verdicts (the headline use case for the multi-call design: one candidate `confirms` eggNOG while another `conflicts`).
- `egn_tcids` is a list because eggNOG's `KEGG_TC` field is multi-valued (comma-separated in source TSV) — e.g. MreB-family proteins carry `["1.A.33.1", "9.B.157.1"]`. Diamond matching ANY value yields `confirms` / `refines`; only ALL-disagree counts as `conflicts`.
- `egn_agreement` values: `confirms` (any eggNOG TC matches diamond's family) | `refines` (any eggNOG TC is a strict ancestor of diamond's call) | `extends` (eggNOG had no TC) | `conflicts` (every eggNOG TC disagrees at family level).
- `pfam_ids`: the gene's Pfam annotations (from `gene_annotations_merged.json`). Independent of diamond/eggNOG.
- `pfam_tc_families`: 3-part TC families implied by the gene's Pfams via TCDB's curated Pfam→TC mapping (`~/tools/TCDB/DB/tcdb_pfam_map.tsv`). Empty when the gene has no Pfams or none of them are in TCDB's map.
- `pfam_agreement` values: `confirms_diamond` (Pfam-implied family matches diamond, not eggNOG) | `confirms_eggnog` (matches an eggNOG TC, not diamond) | `confirms_both` (matches both — multi-domain or shared Pfam) | `contradicts_both` (Pfam implies a TC family but neither side called it) | `neutral` (no Pfam signal — gene has no Pfams or its Pfams aren't in TCDB's map). The Pfam signal is independent from `egn_agreement`; combine both for stronger merge-rule decisions.

## Phase 2 (Future)

Phase 1 artifacts sit in the strain cache for inspection. Phase 2 (separate spec) will integrate them into `gene_annotations_merged.json` and the KG. See [docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md §7](../../../docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md).

## Workflow When Invoked

1. Verify `diamond` is on PATH (`which diamond`).
2. Run the skill: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py`. The TCDB diamond DB builds itself on first run (~30s download + 1min makedb).
3. Review the status table — note any FAILED strains.
4. Inspect `<data_dir>/tcdb/<strain>.tcdb.calls.json` for spot checks.
