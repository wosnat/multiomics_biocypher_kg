---
name: tcdb-diamond
description: Run diamond blastp vs. the curated TCDB FASTA per strain to generate per-protein TCDB classifications, with tiered confidence (5-part / 4-part / 3-part) plus per-family consensus collapse. Pure sequence evidence — no eggNOG or gene-annotation inputs. Produces inspectable `<strain>.tcdb.calls.json` artifacts; **Phase 2 KG integration is DONE** — emitted as `Gene_has_tcdb_family` edges of their own (`tcdb_call` diamond_tier and diamond_best_hit, alongside eggnog), each stamped with `tcdb_gate`. Re-running with --force requires `prepare_data.sh --steps 2 6 --force` + a Docker rebuild.
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
| eggNOG agreement (confirms / refines / extends / conflicts) | A query: compare a gene's `Gene_has_tcdb_family` edges of different `tcdb_call` (eggnog / diamond_tier / diamond_best_hit) over the `Tcdb_family_is_a_tcdb_family` hierarchy. Nothing is stored (project 3b reconciles) |
| Pfam corroboration | The `Tcdb_family_has_pfam_domain` bridge edge, built from TCDB's curated Pfam→TC map in prepare_data step 6 |

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
| 1 | 5 parts (`tc_system`) | ≥ 70% | qcov ≥ 70% | High-confidence system call |
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

### `best_hit` block (best_hit_v1, 2026-10-19)

Alongside `calls`, a protein may carry a `best_hit` block — its **actual best-hit 5-part
TCDB system(s)**, which the tier/consensus collapse truncates away (e.g. MED4 cynD is called
`3.A.1` but its best hit is `3.A.1.16.2`, the cyanate/nitrite system). Since 2026-10-23 the
block carries **per-system statistics** in `hits` (best first, aligned with `tcids`), each from
that system's highest-bitscore subject row — no consensus step:

```json
"WP_011131703.1": {
  "calls": [ ... ],
  "best_hit": {
    "rule": "best_hit_v1",
    "tcids": ["3.A.1.203.8", "3.A.1.203.12"],
    "identity": 34.3,
    "attach_tcid": "3.A.1",
    "hits": [
      {"tcid": "3.A.1.203.8",  "identity": 34.3, "qcov": 80.9, "scov": 85.4, "evalue": 1.09e-100, "bitscore": 322.0},
      {"tcid": "3.A.1.203.12", "identity": 32.2, "qcov": 79.7, "scov": 91.8, "evalue": 2.05e-90,  "bitscore": 292.0}
    ]
  }
}
```

Rule (`select_best_hits`, independent of the tier policy): qcov ≥ 75 **and** scov ≥ 75
(scov from gapseq's reference-coverage rule), TC class 1–4, no identity or length floor;
**reciprocal best hit** (the top hit's subject row must pick this protein as its best gated
match in the strain); window = gated hits with bitscore ≥ 0.9 × top, up to 3 distinct 5-part
TCIDs, any family, ordered bitscore desc → TCID → subject. Subject TCIDs longer than 5 parts
are truncated (best-hit path only). `identity` = the top system's identity (kept for
compatibility; `hits[0].identity`). `attach_tcid` = the call in the top hit's 3-part family —
informational only since the KG gives every best-hit system its **own** edge. ~188 proteins
corpus-wide have a best hit but **no** call (top HSP < 50 aa): `{"calls": [], "best_hit":
{...}}` with no `attach_tcid`; they still get best-hit edges. `best_hit_tier(identity)` gives a
best hit's identity band (≥ 70 → 1, 40–70 → 2, < 40 → 3), the `tier` on its KG edge. Summary
keys: `proteins_with_best_hit`, `best_hit_multi_tcid`, `best_hit_without_call`,
`best_hit_fallback_attach` (expected 0). Design:
`docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md` (rule) and
`docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md` §5.3 (per-system stats).

**Regenerating without diamond:** `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --from-tsv [--strains ...]`
rebuilds calls.json + summary from the committed `.tcdb.tsv` (no FASTA, no diamond, always
rewrites). Use it when only the derivation changed — this is how the 2026-10-23 `hits` block and
the `tc_specificity` → `tc_system` rename were applied to the committed calls.json (0 tier-call
diffs). A real diamond rerun against a newer TCDB FASTA changes the existing tier calls.

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
# (step 6 logs a warning naming the exact families + this command). Merges into
# the existing artifact — other families' names are kept:
uv run python -m multiomics_kg.download.scrape_tcdb_names --families 2.A.130 9.B.400
```

**Merge-preserving (2026-10-19):** every run starts from the committed
`tcdb_names.json` and only *adds* names for ids not yet named — existing entries
stay verbatim, so a re-scrape cannot rename an existing `TcdbFamily`. (Before this,
`--families` rewrote the artifact from the listed families alone, dropping every
other family's subfamily/system names.) To pick up new upstream wording for
existing ids, delete those entries first.

**When to run:** (1) step 6 warns `kept TCDB families not covered by
tcdb_names.json` — a new strain's annotations reached families the scrape
hasn't seen, or `kept TCDB systems have no name` (kept systems — gene-anchored,
best-hit and range-only — render as the bare TCID until named); (2) a TCDB upstream release worth picking up (`--force` re-fetches
cached pages). After a run, rebuild the hierarchy: `bash scripts/prepare_data.sh
--steps 6 --force`, then commit both `tcdb_names.json` and
`tcdb_hierarchy.json`.

**Hygiene (by design):** single-threaded, 2.5 s/request, exponential backoff,
aborts after 3 consecutive failures (tcdb.org has multi-day maintenance
windows — resume later; already-fetched pages under the gitignored
`cache/data/tcdb/raw/pages/` are skipped — the cache is per clone; copy it from
another clone before a run to avoid re-fetching). Scope = the families carrying
any kept id in `tcdb_pruned.json` (gene-anchored, best-hit and range-only systems; a
family-only kept id is fetched via its smallest 4-part child) **plus** the diamond best-hit
systems (merged `tcdb_best_hit_ids`), plus the full browse.php layer
(classes/subclasses/families) — 592 family pages, 9,622 names as of 2026-10-23 (134 kept
systems unnamed upstream). Run order after a kept-set change: step 6 → scraper → step 6
`--force` (names + `direction`). Subfamilies with no upstream name stay as bare
TC ids on purpose (`t.name = t.tcdb_id` is the fallback marker). Design:
`docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md`.

## Phase 2 — KG integration ✅ DONE

These artifacts are merged into `gene_annotations_merged.json` (`tcdb_diamond_ids`, folded
into `transporter_classification`; best hits as `tcdb_best_hit_ids`) and the KG as evidence on
the `TcdbFamily` ontology. Since the TCDB cleanup (2026-10-23,
`docs/kg-changes/tcdb-cleanup.md`) each source is **its own edge**, one per (gene, TC node,
`tcdb_call`):

- `diamond_tier` — the tier-truncated call (`tier`, `identity`, `qcov`, `scov`, `evalue`;
  `confidence_score` / `consensus_n` stay in the artifact only). 41,242 edges.
- `diamond_best_hit` — each `best_hit` system with its own `hits` statistics (`rank`, `tier` =
  identity band, `identity`, `qcov`, `scov`, `evalue`, `bitscore`, `bitscore_ratio`). 17,512 edges.
- `eggnog` — eggNOG `KEGG_TC`, alongside (17,069).

Every edge carries `tcdb_gate` (`direct` | `annotated` | `below_gate`); best-hit systems also
seed TCDB pruning in step 6. `qcov` / `scov` / `identity` are percent.

**Re-running this skill changes the KG.** After a `--force` run, rebuild the caches
and the graph:

```bash
bash scripts/prepare_data.sh --steps 2 6 --force   # merge + re-seed TCDB pruning
# if the kept set gained systems: re-run the name scraper, then --steps 6 --force again
docker compose down && docker compose up -d --build
```

Two consumer-facing notes:

- **No `filter_action`.** All candidates become edges, tier 3 included. Downstream reads
  `tcdb_gate`: `'tcdb'` in `annotation_types` / `annotation_quality` and the TCDB gene counts use
  `tcdb_gate IN ['direct','annotated']` (eggNOG, tier ≤ 2 calls, rank-1 best hits with tier ≤ 2);
  the transport-arm metabolite counts use `direct` only. Tier-3 calls and near-tie best hits are
  `below_gate` — findable, but they feed no count.
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
