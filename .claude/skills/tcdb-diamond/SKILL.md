---
name: tcdb-diamond
description: Run diamond blastp vs. the curated TCDB FASTA per strain to generate per-protein TCDB classifications, with tiered confidence (5-part / 4-part / 3-part) plus per-family consensus collapse. Pure sequence evidence — no eggNOG or gene-annotation inputs. Produces inspectable `<strain>.tcdb.calls.json` artifacts; **Phase 2 KG integration is DONE** — merged as a second evidence source on `Gene_has_tcdb_family` alongside eggNOG (`sources: ['eggnog','diamond']`). Re-running with --force requires `prepare_data.sh --steps 2 6 --force` + a Docker rebuild.
argument-hint: "[--strains <name> ... | --force | --refresh-tcdb | --threads <n>]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(diamond *)
---

# TCDB-Diamond Skill

Per-strain `diamond blastp` against the curated TCDB FASTA — a **second, independent
evidence source** for TCDB transporter classifications alongside eggNOG's `KEGG_TC`.

Measured across 42 strains, diamond's contribution is **breadth, not depth**: 15,074
proteins receive a TC call with no eggNOG TC at all, and 762 TC IDs appear that eggNOG
never mentions. (eggNOG is already the deeper source — 240 distinct 5-part IDs vs
diamond's 36.)

## Pure sequence evidence — by design

This runner reads **only** `protein.faa` and the TCDB diamond DB. It deliberately does
**not** read eggNOG annotations or `gene_annotations_merged.json`.

An earlier version did both, which created a cycle: `gene_annotations_merged.json` is
the *output* of prepare_data step 2, and step 2 is exactly where these calls get merged.
Re-running step 2 silently mutated `calls.json` even when the diamond blast was
byte-identical.

Cross-source reconciliation now happens downstream, where both sources are already present:

| Concern | Where it lives now |
|---|---|
| eggNOG agreement (confirms / refines / extends / conflicts) | Derivable from `Gene_has_tcdb_family.sources` + the `Tcdb_family_is_a_tcdb_family` hierarchy |
| Pfam corroboration | The `Pfam_in_tcdb_family` bridge edge, built from TCDB's curated Pfam→TC map in prepare_data step 6 |

There is also **no `filter_action`**. The old post-hoc filter chain made a candidate's
verdict depend on whether the protein happened to have *unrelated* sibling candidates.
The tier policy below is the quality gate, and it is sibling-independent by construction.

Full rationale: [`docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md`](../../../docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md) §3.

## Quick Start

```bash
# Run all genome strains (skip already-done)
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py

# Run a single strain
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strains MED4 MIT9313

# Force re-run even if calls.json exists
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --force

# Re-download TCDB FASTA + diamond DB
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --refresh-tcdb

# Use more threads
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --threads 8

# Smoke test: first 100 proteins of one strain (~10-30s)
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strains MIT1002 --limit 100
```

## One-Time Setup

No manual setup required. The skill downloads one resource on first run:

- **TCDB FASTA** from `https://tcdb.org/public/tcdb` → rewritten + built into a diamond DB at `~/tools/TCDB/DB/tcdb.dmnd`. Re-downloaded on `--refresh-tcdb`.

> TCDB's curated **Pfam→TC map** is no longer downloaded here — it moved to prepare_data
> step 6 (alongside the other TCDB reference TSVs) and backs the `Pfam_in_tcdb_family`
> bridge edge. See `multiomics_kg/utils/tcdb_utils.py`.

Optional `.env` entry to relocate the TCDB data dir (default: `~/tools/TCDB`):

```
TCDB_DATA_DIR=/path/to/TCDB
```

System tool required: `diamond` (already required by other parts of the project).

**Note:** Saier Lab's `extractTCDB.pl` performs the same FASTA download + header rewrite, but contains a bash-only redirect (`>&/dev/null`) that fails on Ubuntu's default dash-based `/bin/sh`. The skill replicates the operation in Python to sidestep the portability issue, so no `git clone` of TCDBtools is needed.

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`, filters to `organism_type='genome_strain'`.
2. Ensures `~/tools/TCDB/DB/tcdb.dmnd` exists, downloading + building it if missing or `--refresh-tcdb`.
3. For each strain, runs `diamond blastp` of `<data_dir>/protein.faa` against the TCDB diamond DB, writing raw 8-column TSV to `<data_dir>/tcdb/<strain>.tcdb.tsv`.
4. Groups each protein's hits by 3-part TC family, then applies the per-hit tier policy + per-family consensus collapse + best-tier promotion.
5. Tags class-9 (incompletely characterized) calls.
6. Writes `<data_dir>/tcdb/<strain>.tcdb.calls.json` and `<strain>.tcdb.skill_summary.json` (per-strain stats).
7. Prints a status table to stdout. Diamond's full per-strain stdout+stderr are captured to `logs/tcdb_diamond_<strain>.log` (auto-gitignored via `*.log`).

## Tier Policy

| Tier | Truncate to | Identity | Coverage rule | Notes |
|---|---|---|---|---|
| 1 | 5 parts (`tc_specificity`) | ≥ 70% | qcov ≥ 70% | High-confidence specificity call |
| 2 | 4 parts (`tc_subfamily`) | ≥ 40% | qcov ≥ 60% | Solid subfamily call |
| 3 | 3 parts (`tc_family`) | (no floor) | qcov ≥ 40% **OR** scov ≥ 40% | gblast3-style floor |

All tiers also require: e-value ≤ 0.001, HSP length ≥ 50 aa.

**Tier is both a confidence band and the claim depth** — weaker similarity yields a
deliberately *broader* claim. This makes tier 3 conservative, not overconfident: across
42 strains its median identity is 34% but its median e-value is 2.3e-30, and 3.9% of
tier-3 calls have ≥70% identity, demoted because the top hits scattered across
subfamilies (genuine subfamily ambiguity, correctly resolved by claiming only the family).

Effective tier is `max(best_hit_tier, depth_tier)` within a family group — the more
conservative of the strongest hit's identity-tier and the consensus depth's tier.

## Output Schema (`<strain>.tcdb.calls.json`)

Keyed by NCBI protein_id (WP_ accession). Each protein maps to a `calls` list — one
candidate per 3-part TC family the protein hits, sorted by `confidence_score`
descending. Multi-domain proteins (RND + MFP + OMF efflux components, etc.) produce
multiple candidates instead of being rejected.

Every field is derived from the diamond alignment alone. There are no protein-level
cross-source fields and no filter verdict — see "Pure sequence evidence" above.

```json
{
  "WP_010951455.1": {
    "calls": [
      {
        "tcid": "2.A.6.1",
        "level_kind": "tc_subfamily",
        "tier": 2,
        "confidence_score": 0.85,
        "identity": 100.0,
        "qcov": 100.0,
        "scov": 100.0,
        "evalue": 1.2e-200,
        "length": 400,
        "consensus_n": 2,
        "consensus_agreement": "4_part",
        "incompletely_characterized": false
      },
      {
        "tcid": "8.A.1.2",
        "level_kind": "tc_subfamily",
        "tier": 2,
        "confidence_score": 0.55,
        "...": "..."
      }
    ]
  }
}
```

Field semantics:
- `tier` (1/2/3) and `level_kind` come from `effective_tier = max(best_hit_tier, depth_tier)` within this family's hits. `tcid` is truncated to the parts justified by `effective_tier`.
- `confidence_score` ∈ [0,1] = `(identity / 100) × (qcov / 100) × agreement_weight`, where agreement_weight is 1.0 / 0.85 / 0.7 for 5/4/3-part consensus. The list is sorted by this descending — `calls[0]` is the strongest call.
- `identity`, `qcov`, `scov`, `evalue`, `length` come from the **best (highest-identity)** hit in THIS family group.
- `consensus_n` is how many hits backed this family; `consensus_agreement` is the depth at which they agreed.
- `incompletely_characterized` flags TCDB class 9. No demotion is applied — downstream consumers decide.

## Consuming the artifact

`multiomics_kg/utils/tcdb_diamond.py` exposes: `load_calls_json`, `iter_candidates`,
`best_call`, `call_tc_families`, `call_tcids`.

Nothing is pre-filtered. Consumers that want a stricter cut should threshold on `tier`,
`identity`, or `confidence_score` **explicitly**, so the choice is visible at the point of
use rather than frozen into the artifact.

## Refreshing TCDB names (T6)

`TcdbFamily` node names (all levels) come from a committed artifact,
`cache/data/tcdb/tcdb_names.json`, built by a standalone, slow-and-polite
scraper — TCDB publishes no bulk file for subclass/subfamily/specificity names:

```bash
# Full run (~352 requests at 2.5 s spacing ≈ 15-20 min; resumable)
uv run python -m multiomics_kg.download.scrape_tcdb_names

# Top-up after onboarding a strain that introduced unscraped families
# (step 6 logs a warning naming the exact families + this command):
uv run python -m multiomics_kg.download.scrape_tcdb_names --families 2.A.130 9.B.400
```

**When to run:** (1) step 6 warns `kept TCDB families not covered by
tcdb_names.json` — a new strain's annotations reached families the scrape
hasn't seen; (2) a TCDB upstream release worth picking up (`--force` re-fetches
cached pages). After a run, rebuild the hierarchy: `bash scripts/prepare_data.sh
--steps 6 --force`, then commit both `tcdb_names.json` and
`tcdb_hierarchy.json`.

**Hygiene (by design):** single-threaded, 2.5 s/request, exponential backoff,
aborts after 3 consecutive failures (tcdb.org has multi-day maintenance
windows — resume later; already-fetched pages under the gitignored
`cache/data/tcdb/raw/pages/` are skipped). Scope = the ~351 families carrying
kept 4/5-part IDs in `tcdb_pruned.json`, plus the full browse.php layer
(classes/subclasses/families). Subfamilies with no upstream name stay as bare
TC ids on purpose (`t.name = t.tcdb_id` is the fallback marker). Design:
`docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md`.

## Phase 2 — KG integration ✅ DONE

These artifacts are merged into `gene_annotations_merged.json` and the KG as a
**second evidence source on the existing `TcdbFamily` ontology** — one
`Gene_has_tcdb_family` edge per (gene, TC ID) carrying `sources: ['eggnog','diamond']`
provenance plus the diamond evidence fields (`tier`, `confidence_score`, `identity`,
`qcov`, `evalue`, `consensus_n`).

Coverage across 42 strains: genes with a TC call went 11,103 → **30,076**
(10,278 corroborated by both sources, 18,973 diamond-only).

**Re-running this skill changes the KG.** After a `--force` run, rebuild the caches
and the graph:

```bash
bash scripts/prepare_data.sh --steps 2 6 --force   # merge + re-seed TCDB pruning
docker compose down && docker compose up -d --build
```

Two consumer-facing notes:

- **No `filter_action`.** All candidates become edges; the tier policy is the quality
  gate. Post-import folds `'tcdb'` into `annotation_types` / `annotation_quality` only
  for eggNOG-sourced or **tier ≤ 2** edges, so tier-3 calls are findable without
  inflating annotation quality. Filter on `r.tier` explicitly if you want a stricter cut.
- **`transporter_classification` is a UNION** of both sources; `tcdb_eggnog_ids` /
  `tcdb_diamond_ids` preserve attribution.

Design: [`docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md`](../../../docs/superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md)
(supersedes the Phase-2 sketch in the 2026-05-10 spec §7, which assumed a specificity
win the real data does not show).
Release notes: [`docs/kg-changes/tcdb-cazy-ontologies.md`](../../../docs/kg-changes/tcdb-cazy-ontologies.md).

## Workflow When Invoked

1. Verify `diamond` is on PATH (`which diamond`).
2. Run the skill: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py`. The TCDB diamond DB builds itself on first run (~30s download + 1min makedb).
3. Review the status table — note any FAILED strains.
4. Inspect `<data_dir>/tcdb/<strain>.tcdb.calls.json` for spot checks.
