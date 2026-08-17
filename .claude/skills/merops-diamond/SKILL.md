---
name: merops-diamond
description: Run diamond blastp vs. the MEROPS peptidase scan library per strain to classify proteases/peptidases and their inhibitors into MEROPS families (clan / family / subfamily / identifier), with tiered confidence following the tcdb-diamond policy plus non-peptidase-homolog and inhibitor flags. Pure sequence evidence. Produces inspectable `<strain>.merops.calls.json` artifacts; Phase 1 only — no KG integration yet. Triggers: "run merops", "classify proteases/peptidases", "protease families for strain X", "exoprotease classification".
argument-hint: "[--strains <name> ... | --force | --limit N | --threads N | --refresh-merops]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(diamond *)
---

# MEROPS-Diamond Skill

Per-strain `diamond blastp` of `protein.faa` against the curated **MEROPS scan
library** (`merops_scan.lib`, 5,009 representative peptidase-unit sequences —
one per MEROPS identifier, release 12.5) to classify proteolytic enzymes and
their inhibitors into the MEROPS hierarchy: clan → family (S08) → subfamily
(S08A) → identifier (S08.001). Phase 1 — committed per-strain JSON artifacts,
no KG coupling. A tcdb-diamond clone: same tier policy, same consensus
collapse, same pure-sequence-evidence contract.

Design: [`docs/superpowers/specs/2026-08-17-merops-diamond-design.md`](../../../docs/superpowers/specs/2026-08-17-merops-diamond-design.md)

## Pure sequence evidence — by design

Reads **only** `protein.faa`, the MEROPS diamond DB, and the static
family→clan reference (`family.txt`). Deliberately does NOT read eggNOG
annotations or `gene_annotations_merged.json` (the tcdb-diamond lesson: that
creates a cycle with prepare_data step 2, where these calls would eventually
merge). Unlike TCDB there is **no second sequence source** for MEROPS — eggNOG
does not emit MEROPS ids and InterProScan has no MEROPS member DB — so any
Phase-2 evidence score starts single-source; corroboration options are the
MEROPS-published `interpro.txt` Pfam/InterPro map and `GO_annotation.txt`.

## Quick Start

```bash
# Run all genome strains (skip already-done)
uv run python .claude/skills/merops-diamond/run_merops_diamond.py

# Run one or more strains by name
uv run python .claude/skills/merops-diamond/run_merops_diamond.py --strains MED4 MIT9313

# Force re-run even if calls.json exists
uv run python .claude/skills/merops-diamond/run_merops_diamond.py --force

# Re-download the MEROPS scan library + reference tables, rebuild diamond DB
uv run python .claude/skills/merops-diamond/run_merops_diamond.py --refresh-merops

# Smoke test: first 100 proteins of one strain (~seconds; output gitignored)
uv run python .claude/skills/merops-diamond/run_merops_diamond.py --strains MED4 --limit 100
```

## One-Time Setup

No manual setup required. On first run (or `--refresh-merops`) the skill
downloads three files from `https://ftp.ebi.ac.uk/pub/databases/merops/current_release/`
into `~/tools/MEROPS/DB/`:

- `merops_scan.lib` (**1.9 MB**, 5,009 sequences) → headers rewritten to
  `>MERNUM|merops_id|subfam_token` → `diamond makedb` → `merops.dmnd` (~2 MB)
- `database_files/family.txt` + `database_files/clan.txt` (family → clan /
  name / peptidase-vs-inhibitor reference; a few hundred KB)

Total disk footprint: **< 10 MB** (the smallest tool DB in the project).

Optional `.env` entry to relocate (default `~/tools/MEROPS`):

```
MEROPS_DATA_DIR=/path/to/MEROPS
```

System tool required: `diamond` on PATH (already required by tcdb-diamond).

Gotchas:
- `merops_scan.lib` is **latin-1 encoded** (0x96 en-dashes in peptidase
  names), not UTF-8 — the runner opens it with `encoding="latin-1"`.
- MEROPS has been in maintenance mode since release 12.5 (Sept 2023); the
  `current_release` URL is stable. `--refresh-merops` is only needed if EBI
  ships a new release.

Verification: `ls ~/tools/MEROPS/DB/` should show `merops.dmnd`,
`merops.faa`, `merops_scan.raw.lib`, `family.txt`, `clan.txt`.

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` via
   `load_genome_rows`, optionally filtered by `--strains`.
2. Ensures `~/tools/MEROPS/DB/merops.dmnd` + `family.txt` exist (downloads +
   builds if missing or `--refresh-merops`).
3. For each strain, runs `diamond blastp -q <data_dir>/protein.faa -d merops.dmnd
   --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length evalue bitscore
   --evalue 0.001 --max-target-seqs 25 --more-sensitive` → raw TSV at
   `<data_dir>/merops/<strain>.merops.tsv`.
4. Groups each protein's accepted hits by MEROPS **family**, applies the
   per-hit tier policy + per-family consensus collapse + best-tier promotion
   (`multiomics_kg/utils/merops_diamond.py`, unit-tested).
5. Tags each candidate with clan / catalytic type / entry type (peptidase vs
   inhibitor) and non-peptidase-homolog evidence from the identifier tail.
6. Writes `<strain>.merops.calls.json` + `<strain>.merops.skill_summary.json`
   via `multiomics_kg/utils/tool_calls_io.py`.
7. Prints a status table; full diamond stdout+stderr go to
   `logs/merops/<strain>.log` (gitignored).

## Tier Policy

Thresholds are exact tcdb-diamond parity; only the claim vocabulary differs.

| Tier | Claim (level_kind) | Identity | Coverage rule |
|---|---|---|---|
| 1 | MEROPS identifier (`merops_id`, e.g. S08.036) | ≥ 70% | qcov ≥ 70% |
| 2 | subfamily (`merops_subfamily`, e.g. S08A) | ≥ 40% | qcov ≥ 60% |
| 3 | family (`merops_family`, e.g. S08) | (no floor) | qcov ≥ 40 **OR** scov ≥ 40 |

All tiers also require: e-value ≤ 0.001, HSP length ≥ 50 aa.

Consensus collapse per family group: agreement depth `id` > `subfamily` >
`family`, weights 1.0 / 0.85 / 0.7;
`effective_tier = max(best_hit_tier, depth_tier)`.

**Ragged-hierarchy note (deviation from TCDB):** ~44% of scan-library entries
sit in families with no subfamilies. There, tier-2 confidence claims land at
family level — `tier: 2` with `level_kind: "merops_family"`. `tier` is the
confidence band; `level_kind` is the actual claim depth.

## Output Schema (`<strain>.merops.calls.json`)

Keyed by NCBI protein_id (WP_ accession). Each protein maps to a `calls` list —
one candidate per MEROPS family the protein hits, sorted by
`confidence_score` descending. Multi-domain proteins emit multiple candidates.

```json
{
  "WP_049586340.1": {
    "calls": [
      {
        "code": "S08A",
        "family": "S08",
        "subfamily": "S08A",
        "clan": "SB",
        "catalytic_type": "S",
        "entry_type": "peptidase",
        "level_kind": "merops_subfamily",
        "tier": 2,
        "confidence_score": 0.4335,
        "identity": 51.0,
        "qcov": 85.0,
        "scov": 80.0,
        "evalue": 1.2e-100,
        "length": 400,
        "consensus_n": 5,
        "consensus_agreement": "subfamily",
        "best_hit_id": "S08.036",
        "best_hit_mernum": "MER0001234",
        "best_hit_kind": "holotype",
        "homolog_hit_fraction": 0.0
      }
    ]
  }
}
```

Field semantics:
- `code` is the **called** classification, truncated to the depth
  `effective_tier` justifies: identifier / subfamily / family. `family` is
  always present; `subfamily` is null for family-level claims and for
  families without subfamilies; `clan` is null when MEROPS has not assigned
  the family to a clan (clan codes ending `-` are the unassigned sentinel).
- `catalytic_type` is the family's first letter (S/C/M/A/T/G/N/U/P where P =
  mixed); **null for inhibitor families**. `entry_type` is
  `peptidase` | `inhibitor` (I-families are kept — protease inhibitors are
  biology — and flagged, never silently dropped).
- `confidence_score` ∈ [0,1] = `(identity/100) × (qcov/100) × agreement_weight`
  (1.0 / 0.85 / 0.7 for id/subfamily/family consensus). `calls[0]` is the
  strongest call.
- `identity`, `qcov`, `scov`, `evalue`, `length` come from the best
  (highest-identity) hit in this family group; `consensus_n` is how many hits
  backed it; `consensus_agreement` the depth they agreed at.
- `best_hit_id` / `best_hit_mernum` identify that best hit in MEROPS;
  `best_hit_kind` decodes its identifier tail: `holotype` (characterized
  peptidase, tails 001–899), `putative` (letter tails, e.g. S01.A33),
  `nonpeptidase_homolog` (9xx tails — **catalytically dead** relatives).
  `homolog_hit_fraction` is the fraction of the group's hits that are 9xx
  entries. A diamond hit alone cannot distinguish an active peptidase from a
  dead homolog — these two fields carry the evidence; nothing is filtered.

`<strain>.merops.skill_summary.json` — per-strain QC: `strain`,
`tool_version` (diamond version line), `input_proteins`, `wallclock_s`,
`raw_hit_lines`, `parse_failures`, `proteins_with_hits`,
`proteins_with_call`, `total_candidates`,
`candidates_per_protein_distribution`, `tier_distribution`,
`consensus_agreement_distribution`, `catalytic_type_distribution` (keys =
catalytic letters + `inhibitor`), `entry_type_distribution`,
`best_hit_kind_distribution`.

## QC

### Per-strain QC — what to inspect in `<strain>.merops.skill_summary.json`

- `parse_failures == 0` — non-zero means diamond's output or the rewritten
  subject-ID format drifted from what the parser expects.
- `proteins_with_call / input_proteins` ≈ **2–6%** — peptidases + inhibitors
  are a small fraction of any proteome. Near-zero means a broken DB; ≫10%
  means something is wrong with the tier floor.
- `tier_distribution` — tier 3 dominates (remote homology), tier 1 is a small
  minority. A strain with zero tier-1 calls is unusual but not alarming.
- `best_hit_kind_distribution` — `nonpeptidase_homolog` should be a minority
  of candidates; if it rivals `holotype`, inspect before trusting counts.
- `entry_type_distribution` — inhibitors ≪ peptidases everywhere.
- `wallclock_s` — seconds per strain (5K-sequence DB). 10× outliers signal a
  bad input.

(Calibrated against the 42-strain 2026-08-17 batch — see "Observed batch
results" below.)

### Cross-strain QC narrative

Observed across the 42-strain batch (2026-08-17):

- **Lifestyle gradient, exactly as predicted.** Streamlined *Prochlorococcus*
  HL strains: 57–63 peptidase-carrying proteins (~3.1% of the proteome);
  larger cyanobacteria (MIT9313, CC9311, PCC7002): 72–113; heterotrophs run
  2–3×: Alteromonas 136–182, KT2440 173, W3-18-1 153. A heterotroph landing
  in the Pro range (or vice versa) is a red flag.
- **Secretion cross-check (vs signalp calls.json):** 64 of MIT1002's 154
  merops-called proteins (42%) carry a signal peptide vs 6 of MED4's 58
  (10%) — heterotroph exoproteases vs cyanobacterial housekeeping proteases.
  All 6 of MIT1002's S08 (subtilisin) calls are SignalP `SP` — the canonical
  secreted-subtilisin signature.
- **Tier shape:** tier 3 dominates everywhere (~92% of candidates), tier 2
  ~8%, tier 1 nearly absent (6 of 4,254 — the scan library holds few marine
  bacterial holotypes, so ≥70% identity is rare). This is expected
  remote-homology behavior, not a defect.
- **Homolog load:** `nonpeptidase_homolog` best-hits are ~18% of candidates
  globally — concentrated in known metabolic-enzyme-adjacent families (C26
  gamma-glutamyl-hydrolase-like, C44 amidotransferases). Treat family counts
  in those families with care.
- Top families are unspectacular housekeeping sets (S33, S09, C26, M38, C44,
  S01, M24, M41, S14…), consistent with whole-proteome scanning rather than
  secretome enrichment.

### Spot checks

| Strain | Protein ID | Expected | Why this is the ground truth |
|---|---|---|---|
| MED4 | WP_011132947.1 | family S14, clan SK | PMM1313, ATP-dependent Clp protease proteolytic subunit — ClpP is the S14 type peptidase |
| MED4 | WP_011131865.1 | family M41, clan MA | PMM0226 FtsH1 — FtsH is the M41 type peptidase |
| WH7803 | WP_011932280.1 | family M41, clan MA | SynWH7803_0355 FtsH1 — Synechococcus representative of the M41 check |
| MIT1002 | WP_049586340.1 | family S08, clan SB | MIT1002_02011, annotated "S8 family serine peptidase" (subtilisin); also SignalP `SP` — secreted subtilisin |

```bash
jq '."WP_011132947.1".calls[0] | {code, family, clan, tier}' \
  cache/data/Prochlorococcus/genomes/MED4/merops/MED4.merops.calls.json
# Expected: family "S14", clan "SK"
```

## Observed batch results (42-strain run, 2026-08-17)

4,194 proteins classified (4,254 candidates — 60 proteins hit ≥2 families);
zero FAILED strains, zero parse failures. Cross-strain distributions
(candidates):

| Dimension | Distribution |
|---|---|
| Tier | T3 3,917 (92.1%) · T2 331 (7.8%) · T1 6 |
| Catalytic type | S 1,728 · M 1,527 · C 662 · T 124 · A 101 · inhibitor 52 · N 33 · U 23 · P 4 |
| Best-hit kind | holotype 2,523 (59%) · putative 954 (22%) · nonpeptidase_homolog 777 (18%) |
| Top families | S33 412 · S09 298 · C26 272 · M38 174 · C44 169 · S01 130 · M24 129 · M41 127 · S14 125 · M23 120 |

Per-strain proteins-with-call: **57** (Pro HL strains, ~1,850 proteins) →
**182** (Alt_MarRef, ~4,300 proteins). Per-strain wallclock: **0.5–0.7 s**
(the 5,009-sequence DB makes this the fastest tool in the project; the full
42-strain batch is under a minute).

## Phase 2 (Future)

NOT yet wired into `gene_annotations_merged.json` or any KG adapter — the
artifacts sit in the strain cache for inspection. The Phase-2 sketch
(spec §"Phase 2 sketch"): CAZy-shaped `MeropsFamily` ontology (clan → family →
subfamily, observed-only; catalytic type as node *property* — clans mix
catalytic types) + scored `Gene_has_merops_family` edges via a new
`merops_diamond` logical source, plus MEROPS-published Pfam/InterPro/GO
bridges. No substrate→Metabolite arm (peptidase substrates are proteins);
cleavage specificity becomes node properties. Integrate via
`/integrate-a-tool` with a fresh design spec.

## Workflow When Invoked

1. Verify `diamond` is on PATH (`diamond version`).
2. Run the batch: `uv run python .claude/skills/merops-diamond/run_merops_diamond.py`
   (DB self-installs on first run, < 1 min).
3. Review the status table — chase any FAILED rows via `logs/merops/<strain>.log`.
4. Verify the spot checks above with the `jq` one-liner.
5. Commit the per-strain `calls.json` + `skill_summary.json` artifacts
   (smoke-test `*.limited_*` files are gitignored).
