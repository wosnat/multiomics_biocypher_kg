# TCDB-Diamond Augmentation — Design

**Date:** 2026-05-10 (amended 2026-05-11: §3/§5/§6.1 — Python-native FASTA download + `diamond makedb`, no `extractTCDB.pl`; §6.2 — `--max-target-seqs 25` (was 5); §6.4-A — `effective_tier = max(best_hit_tier, depth_tier)`, plus **multi-candidate emission** — group hits by 3-part TC family, emit one candidate per family instead of forcing a single global consensus; §6.4-B — multi-valued eggNOG `KEGG_TC` aggregation; §6.4-D NEW — Pfam-corroboration tag using TCDB's curated Pfam→TC mapping; §6.5 — schema now wraps candidates in a `calls: list`, with `egn_tcids`, `pfam_ids`, `pfam_tc_families` at protein level and per-candidate `egn_agreement`/`pfam_agreement`)
**Status:** Draft, pending real-data refinement
**Related:** `cache/data/tcdb/{tcdb_hierarchy,tcdb_pruned}.json` (built by prepare_data step 6); `multiomics_kg/adapters/tcdb_adapter.py`; `multiomics_kg/download/build_gene_annotations.py`; `.claude/skills/eggnog-run/`; `.claude/skills/refresh-mnx/`.

## 1 — Problem

Today the only source of TCDB classifications in the KG is **eggNOG's `KEGG_TC` column**, which is family-level (3-part TCIDs like `1.A.11`, `3.A.1.27`) only. Consequences:

1. **Specificity loss.** The `tc_specificity` (5-part) leaves of the curated TCDB hierarchy are essentially unused by genes. `Tcdb_family_transports_metabolite` edges *are* materialized at every kept ancestor (via the step-6 `subtree_substrates` rollup), so substrate→metabolite paths still work — but the `gene_count` / `organism_count` of leaf nodes is 0, and per-gene substrate inference is one ancestor coarser than the underlying evidence supports.
2. **Coverage gaps.** Genes whose seed ortholog has no `KEGG_TC` get no TCDB call at all, even when a direct sequence comparison against the curated TCDB FASTA would find a strong hit.
3. **No verification surface.** A second independent source could cross-check eggNOG's family-level call (paralogs, ambiguous transporters, mis-annotated multi-domain proteins).

## 2 — Goal (this spec — Phase 1 only)

Build a self-contained `/tcdb-diamond` skill that runs **diamond blastp vs. the curated TCDB FASTA** per strain, applies a tiered confidence policy + lightweight post-steps, and writes a per-strain `<strain>.tcdb.calls.json` artifact suitable for inspection.

The artifact is **not yet wired into `gene_annotations_merged.json` or the KG**. That integration is **deferred to Phase 2** (separate spec), to be designed once we can inspect real diamond↔eggNOG agreement patterns on a few strains. Premature merge-rule design risks locking in choices before we know what the data looks like.

Phase 1 deliverable: a reproducible per-strain JSON that addresses all three goals (specificity + coverage + verification — goal D from the brainstorm) at the data-production layer, leaving the integration layer for a follow-up.

Tier thresholds, agreement-tag taxonomy, and the eventual merge rule are explicitly intended to be **revisited once we can inspect real data**. Build first, calibrate after.

## 3 — Architecture (Phase 1)

Single component:

| Component | Where | When it runs | Input | Output |
|---|---|---|---|---|
| **Heavy compute** — diamond run + tier policy + post-steps | `.claude/skills/tcdb-diamond/` (skill) | Opt-in (`/tcdb-diamond`) | `protein.faa` per strain | `<strain>.tcdb.tsv` (raw diamond) + `<strain>.tcdb.calls.json` (per-protein tiered call) + `<strain>.tcdb.skill_summary.json` (per-strain stats) |

This mirrors `/eggnog-run`: opt-in heavy compute that produces inspectable per-strain artifacts. No coupling to `build_gene_annotations.py` or the KG in this phase — the artifacts sit in the strain cache for inspection and downstream Phase 2 work.

**FASTA + DB build:** the skill downloads the curated TCDB FASTA from `https://tcdb.org/public/tcdb` and builds the diamond DB itself in pure Python (download → header rewrite → `diamond makedb`). Saier Lab's [TCDBtools](https://github.com/SaierLaboratory/TCDBtools) `extractTCDB.pl` performs the same steps, but contains a bash-only redirect (`>&/dev/null`) that fails on Ubuntu's default dash-based `/bin/sh`; replicating in Python sidesteps the portability issue and removes the perl dependency. The post-step logic (tier policy, consensus collapse, eggNOG agreement) is custom hierarchy-aware code distinct from `gblast3.py`'s TMS-based pipeline.

**Why Phase 2 is a separate spec:** the merge rule between diamond and eggNOG (union vs. resolve-conflicts vs. confidence-aware), the YAML config-driven extension shape, and the KG adapter implications all hinge on what real-data agreement patterns look like. Designing them up-front from intuition would lock in choices we'd likely revise after first inspection. Phase 1 produces the inspection artifact; Phase 2 designs the integration.

## 4 — File layout

**Shared resource directory** (per-machine, gitignored, sharable across checkouts):

```
~/tools/TCDB/                          # default; override via TCDB_DATA_DIR env var
└── DB/                                # downloaded data (skill-managed)
    ├── tcdb.raw.faa                   # raw download from https://tcdb.org/public/tcdb (TCDB-canonical headers)
    ├── tcdb.faa                       # rewritten headers (>lcl|<accession>-<tcid>); diamond makedb input
    ├── tcdb.dmnd                      # produced by `diamond makedb --in tcdb.faa -d tcdb`
    ├── tcdb_acc2tcid.tsv              # parsed from FASTA headers; sanity vs cache/data/tcdb/raw/acc2tcid.tsv
    └── refresh.log                    # mtime + record count of last refresh
```

`~/tools/TCDB/DB` can be safely deleted to force a clean re-download (or use `--refresh-tcdb`).

The `TCDB_DATA_DIR` env var follows the same pattern as `MNX_DATA_DIR` / `EGGNOG_DATA_DIR`:
- Read via `dotenv.load_dotenv(REPO_ROOT / ".env")` then `os.environ.get("TCDB_DATA_DIR")`.
- Default: `~/tools/TCDB/`.
- The skill prints the resolved location on every run.

**Per-strain cache** (in-repo, gitignored, mirrors `eggnog/`):

```
cache/<organism>/genomes/<strain>/tcdb/
├── <strain>.tcdb.tsv                  # raw diamond blastp output (8-column tab)
├── <strain>.tcdb.calls.json           # post-filtered, tiered TC calls keyed by protein_id
└── <strain>.tcdb.skill_summary.json   # per-strain skill output (tier dist, agreement dist)
```

Estimated size: ~0.5 MB per cyano strain, ~1 MB per heterotroph; ~10–15 MB total across all 25 strains. Raw TSV retained so the tier policy can be tuned without re-running diamond.

## 5 — One-time setup

No manual setup required. The skill creates `~/tools/TCDB/DB/` and downloads the TCDB FASTA on first run.

Optional `.env` entry:
```
TCDB_DATA_DIR=~/tools/TCDB
```

Required system tool: `diamond` (already required for `scripts/map_img_to_ncbi_proteins.py`).

## 6 — Skill: `/tcdb-diamond`

**Files:**
- `.claude/skills/tcdb-diamond/SKILL.md` — frontmatter + workflow doc (mirrors `/eggnog-run`)
- `.claude/skills/tcdb-diamond/run_tcdb_diamond.py` — orchestrator

**CLI (mirrors `/eggnog-run`):**
```bash
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py             # all strains, skip done
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MED4
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --force     # re-run even if calls.json present
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --refresh-tcdb  # re-download TCDB FASTA + diamond DB
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --threads 8   # default: os.cpu_count() or 4
```

### 6.1 — TCDB DB build (skill, run lazily)

Triggered when `~/tools/TCDB/DB/tcdb.dmnd` is missing OR `--refresh-tcdb` is given. The skill creates `~/tools/TCDB/DB/` if needed, then performs three steps in pure Python:

1. **Download** the curated TCDB FASTA via `urllib.request.urlretrieve("https://tcdb.org/public/tcdb", DB/tcdb.raw.faa)`.
2. **Rewrite headers** from TCDB-canonical `>gnl|TC-DB|<acc>|<tcid> <description>` to parser-friendly `>lcl|<acc>-<tcid>` (matches `parse_tcdb_subject_id` in `multiomics_kg/utils/tcdb_diamond.py`). Writes `DB/tcdb.faa`. Sequences with malformed headers are silently skipped.
3. **Build the diamond DB** via `diamond makedb --in DB/tcdb.faa -d DB/tcdb` → `DB/tcdb.dmnd`.

After build, the skill parses `DB/tcdb.faa` headers into `DB/tcdb_acc2tcid.tsv` (two columns: accession, 5-part TCID). Sanity-check: count of **distinct 5-part TCIDs** in the new file should be ≥ the count of distinct values in column 2 of `cache/data/tcdb/raw/acc2tcid.tsv` from prepare_data step 6 (TCDB grows; never shrinks meaningfully). Mismatch → warn, do not fail.

### 6.2 — Per-strain diamond invocation

```bash
diamond blastp \
  -q cache/<org>/genomes/<strain>/protein.faa \
  -d ~/tools/TCDB/DB/tcdb.dmnd \
  -o cache/<org>/genomes/<strain>/tcdb/<strain>.tcdb.tsv \
  --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length evalue bitscore \
  --evalue 0.001 \
  --max-target-seqs 5 \
  --more-sensitive \
  --threads <N>
```

**Key choices:**
- **No `--id` floor at the diamond step.** Saier's `gblast3.py` runs BLAST without an identity threshold; identity tiering happens in code. We do the same.
- **No `--query-cover` floor at the diamond step.** Diamond can only AND query/subject coverage; gblast3's rule is OR. We apply the OR rule in Python.
- **`--max-target-seqs 5`** preserves enough hits to compute the consensus collapse (post-step A) without bloating the TSV.
- **`--more-sensitive`**, not `--very-sensitive` — marginal gain at this identity range, ~5× the runtime.

### 6.3 — Per-hit tier policy

Each row in the diamond TSV is classified into a tier based on identity, coverage, length, and e-value. The TCID from the matched TCDB FASTA header is **truncated to `tier_parts` segments** before being recorded.

| Tier | Truncate to | Identity | Coverage rule | Length | E-value | Maps to ontology level |
|---|---|---|---|---|---|---|
| 1 | 5 parts | ≥ 70% | qcov ≥ 70% | ≥ 50 | ≤ 0.001 | `tc_specificity` |
| 2 | 4 parts | ≥ 40% | qcov ≥ 60% | ≥ 50 | ≤ 0.001 | `tc_subfamily` |
| 3 (gblast3 floor) | 3 parts | (none) | qcov ≥ 40% **OR** scov ≥ 40% | ≥ 50 | ≤ 0.001 | `tc_family` |

Hits failing tier 3 are dropped. Tier 3 is identical to gblast3's acceptance rule (`bin/gblast3.py:316,320`); tiers 1–2 ratchet confidence using identity. Thresholds are **explicitly subject to revision** once we have real-data agreement statistics.

### 6.4 — Post-diamond steps (per-protein, in `run_tcdb_diamond.py`)

After parsing the raw TSV but before writing `<strain>.tcdb.calls.json`:

**A — Multi-candidate emission + best-tier promotion.** Diamond's `--max-target-seqs 25` produces up to 25 hits per query. We do **not** force a single per-protein call. Instead:

1. **Group hits by 3-part TC family.** All hits with the same `1.A.11.*` or `2.A.6.*` prefix go to the same group.
2. **One candidate per family.** Run consensus_collapse + best-tier promotion **within each group**. A protein with hits in N distinct 3-part families produces N candidates. Single-family proteins are unaffected (one candidate, same as old behavior).
3. **No global rejection.** The old "disagree at 3 parts → reject the protein" branch is gone — the rejected proteins (~16% of hits) are exactly the multi-family case we now want to recover.

Within a single-family group, the consensus rules are unchanged:

| Hits in the group | Consensus depth | depth_tier |
|---|---|---|
| All agree at 5 parts | `5_part` | 1 |
| Agree at 4 parts but disagree at 5 | `4_part` | 2 |
| Agree at 3 parts but disagree at 4 | `3_part` | 3 |

Two signals drive the final tier and TCID truncation depth (per candidate):

- `depth_tier` (above): how confident the **shared prefix** is.
- `best_tier`: the per-hit identity-tier of the **strongest** hit in the group (per §6.3 thresholds).

`effective_tier = max(best_tier, depth_tier)` — the more conservative of the two. Using best (not worst) honors the strongest evidence: when 3 hits agree at 4-part but vary in identity (e.g. one tier-1 hit + two tier-3 hits), the strong hit's tier-1 promotes to tier-2 via the consensus floor rather than being dragged down to tier-3 by the weakest hit. The TCID is then truncated to the parts justified by `effective_tier` (5/4/3 for tier 1/2/3) — not the consensus depth alone. Metadata fields (`identity`, `qcov`, `scov`, `evalue`, `length`) are sourced from the best (highest-identity) hit.

Output records `consensus_n` (number of hits considered) and `consensus_agreement` (deepest level at which **all** top-N hits share the same prefix: `5_part` | `4_part` | `3_part`). Note: TCDB curates UniProt accessions only at the 5-part `tc_specificity` leaves, so every hit in the diamond TSV carries a 5-part TCID; consensus is therefore always meaningful (no need to handle short-TCID hits).

**Confidence score.** A continuous complement to the discrete `tier`:

```
confidence_score = (best_identity / 100) × (best_qcov / 100) × agreement_weight
```

with `agreement_weight = 1.0 / 0.85 / 0.7` for `5_part` / `4_part` / `3_part` consensus. Lets downstream consumers do their own thresholding without losing the underlying gradient. Range [0, 1], rounded to 4 decimals.

**B — eggNOG agreement tag.** Look up the gene's existing eggNOG `KEGG_TC` value(s) in `<strain>.emapper.annotations` and tag the diamond call. The `KEGG_TC` field is **multi-valued** (comma-separated in the source TSV) — e.g. MreB-family proteins carry `1.A.33.1,9.B.157.1` (the legacy Hsp70-cation-channel call plus the correct MreBCD-family call). All values are inspected; the strongest match wins.

| Tag | Condition |
|---|---|
| `confirms` | diamond call is identical to ANY eggNOG TC, OR a strict descendant of any eggNOG TC's family |
| `refines` | diamond call is a strict descendant of ANY eggNOG TC (e.g. eggNOG `1.A.11`, diamond `1.A.11.1.5`) — the headline specificity win |
| `extends` | eggNOG had no TC values; diamond produced one |
| `conflicts` | EVERY eggNOG TC disagrees with diamond at family level (different first 3 parts) |
| (`egn_only`) | not emitted in Phase 1 — only the diamond-side calls are recorded here; eggNOG-only genes are visible by their absence from the JSON |

Aggregation across multi-valued eggNOG TCs uses precedence `confirms > refines > conflicts` — any single confirming or refining match overrides per-pair conflicts.

**C — Class-9 tag.** TCDB class `9.*` = "Incompletely Characterized Transport Systems". Recorded as `incompletely_characterized: true` in the JSON. **No demotion** — let merge / downstream consumers decide.

**D — Pfam-corroboration tag.** Cross-check the gene's Pfam domain annotation against TCDB's curated Pfam→TC mapping (`https://www.tcdb.org/cgi-bin/projectv/public/pfam.py`, ~1.3K Pfams covering ~8.3K (Pfam, TC) pairs, cached as `~/tools/TCDB/DB/tcdb_pfam_map.tsv`). The mapping returns one or more 3-part TC families implied by each Pfam; we union across all of the gene's Pfams to get the set of TC families the Pfam evidence supports.

| Tag | Condition |
|---|---|
| `confirms_diamond` | Pfam-implied family includes diamond's call's family but no eggNOG family — Pfam supports diamond |
| `confirms_eggnog` | Pfam-implied family includes an eggNOG family but not diamond's — Pfam supports eggNOG |
| `confirms_both` | Pfam-implied family includes diamond's family AND at least one eggNOG family — multi-domain protein OR a Pfam that TCDB curates to multiple TC families |
| `contradicts_both` | Pfam-implied family is non-empty but matches neither diamond nor any eggNOG TC — both sequence-based calls may be wrong |
| `neutral` | gene has no Pfam annotations OR none of its Pfams appear in TCDB's curated map — no Pfam signal available |

This is independent from `egn_agreement`: a `conflicts` egn_agreement combined with a `confirms_diamond` pfam_agreement is a strong "diamond wins" signal; a `conflicts` + `confirms_eggnog` is a strong "eggNOG wins" signal. Phase 2's merge rule will combine both signals.

Recorded fields per protein:
- `pfam_ids`: list of the gene's Pfam IDs (passthrough from `gene_annotations_merged.json`)
- `pfam_tc_families`: sorted list of unique 3-part TC families implied by those Pfams via the TCDB map; empty when no signal
- `pfam_agreement`: one of the 5 tags above

Note: TCDB's Pfam→TC map is a *curated* resource and has known coverage gaps — e.g. PF04193 (the SWEET-family Pfam) is absent from the map even though TCDB curates the SWEET family at 2.A.123. Genes whose Pfams aren't in the map get `pfam_agreement = 'neutral'` regardless of how informative the Pfam might be in principle. We accept these gaps rather than augment the map.

**Explicitly NOT included:**
- HMMTOP / TMS topology check (gblast3 uses HMMTOP; we have UniProt `transmembrane_regions` already, but coverage is too partial — 5–16% of genes by strain — to be useful as a gate or even a confidence signal; deferred until UniProt coverage is improved or a future skill runs DeepTMHMM).
- CDD domain check (heavyweight; marginal value over hierarchy + consensus).
- Substrate compatibility check vs. KG metabolites (belongs in downstream analysis, not the per-strain calls).
- Vendoring `gblast3.py` (different problem — find novel transporters from scratch with TMS gates; we are augmenting curated proteomes).

### 6.5 — `<strain>.tcdb.calls.json` shape

Keyed by **NCBI protein_id (WP_ accession)** — same join key eggNOG uses, so a future Phase-2 merge step is a direct dict lookup.

```json
{
  "WP_010951455.1": {
    "egn_tcids": ["8.A.1.2.1"],
    "pfam_ids":  ["PF_HlyD"],
    "pfam_tc_families": ["8.A.1"],
    "calls": [
      {
        "tcid": "2.A.6.1", "level_kind": "tc_subfamily", "tier": 1,
        "confidence_score": 0.85,
        "identity": 100.0, "qcov": 100.0, "scov": 100.0,
        "evalue": 1.2e-200, "length": 400,
        "consensus_n": 2, "consensus_agreement": "4_part",
        "egn_agreement": "conflicts",
        "pfam_agreement": "confirms_eggnog",
        "incompletely_characterized": false
      },
      {
        "tcid": "8.A.1.2", "level_kind": "tc_subfamily", "tier": 2,
        "confidence_score": 0.55, "...": "...",
        "egn_agreement": "confirms",
        "pfam_agreement": "confirms_both"
      }
    ]
  }
}
```

**Per-protein** (shared across candidates):
- `egn_tcids` is the list of all eggNOG `KEGG_TC` values for this protein (`[]` when eggNOG had no TC). Multi-valued because eggNOG's column is comma-separated; see §6.4-B for aggregation semantics.
- `pfam_ids` / `pfam_tc_families`: see §6.4-D. Empty when the gene has no Pfam annotation OR none of its Pfams are in TCDB's curated Pfam→TC map.
- `calls`: list of candidates sorted by `confidence_score` descending. One candidate per distinct 3-part TC family the protein hits. Single-family proteins have `len(calls) == 1`; multi-domain proteins (RND + MFP partners, etc.) emit multiple candidates instead of being rejected by global consensus.

**Per-candidate**:
- `tier` and `level_kind` reflect `effective_tier = max(best_hit_tier, depth_tier)` within this family's hits only (§6.4-A); `tcid` is truncated to the parts justified by `effective_tier`.
- `confidence_score` ∈ [0, 1] = `(identity / 100) × (qcov / 100) × agreement_weight` (§6.4-A).
- `identity`, `qcov`, `scov`, `evalue`, `length` are from the **best (highest-identity)** hit in this family's hit group.
- `egn_agreement` / `pfam_agreement` are computed PER candidate against the protein-level `egn_tcids` / `pfam_tc_families`. Different candidates of the same protein can have different verdicts — the headline use case for the multi-call design (one candidate `confirms` eggNOG, another `conflicts`).
- `level_kind` mirrors the existing TCDB hierarchy vocabulary (`tc_class | tc_subclass | tc_family | tc_subfamily | tc_specificity`).

### 6.6 — Skill output summary (stdout)

After all strains complete, print a status table:

```
Strain     Status   Hits  T1   T2   T3   confirms  refines  extends  conflicts
MED4       OK       187    47   91   49        88       47       41          8
MIT9301    OK       194    52   92   50        91       48       43          7
MIT1002    OK       402   118  201   83       192      102       95         13
EZ55       skipped  -      -    -    -          -        -        -          -
...
```

`skipped` = `<strain>.tcdb.calls.json` already present and `--force` not given.

## 7 — Phase 2 (deferred): merge into `gene_annotations_merged.json`

**Out of scope for this spec.** Will be designed as a separate spec after Phase 1 ships and we can inspect real diamond↔eggNOG agreement patterns on at least 3-5 strains.

**Sketch of the likely Phase-2 shape** (recorded here for context, not committed):

- Implement as an additive extension of `config/gene_annotations_config.yaml` rather than bespoke code (mirroring how eggNOG / UniProt / gene_mapping are wired). Likely needs a new source `type: json_calls` (nested-record JSON), a new transform `expand_tcdb_ancestors` (walks the TCDB hierarchy), and an `optional: true` source flag (silently skip when calls.json is missing).
- The structured provenance sidecar (`transporter_classification_sources`) can't fit the YAML field-rule shape; implement it as a small post-merge function alongside the existing `enrich_pfam_fields()`.
- The merge rule itself (union vs. resolve-conflicts vs. confidence-aware) is the open question Phase 2 needs to answer. **The data Phase 1 produces is what unlocks that decision.**
- Phase 2 will also extend the per-strain `step2_metabolism_report.json` with a `tcdb_diamond_merge` block summarizing what flowed into the merge.

Until Phase 2 ships, the per-strain `<strain>.tcdb.calls.json` artifacts sit in the cache for manual inspection and ad-hoc queries. The KG itself is unchanged.

## 8 — Adapter / KG impact (Phase 1)

**None.** Phase 1 produces per-strain artifacts in `cache/<org>/genomes/<strain>/tcdb/` only. `multiomics_kg/adapters/tcdb_adapter.py`, `gene_annotations_merged.json`, and the KG build are untouched. KG-side impact is a Phase 2 concern.

## 9 — Future enhancements (out of scope for Phase 1)

1. **Phase 2: merge into `gene_annotations_merged.json`** — see §7 sketch. Separate spec, designed after Phase 1 produces real data.
2. **Calibration against gblast3.** Run gblast3 (full BLAST + HMMTOP + CDD) on one strain offline and compare its TC assignments against the diamond+tier output. Use agreement rate to validate (or retune) tier thresholds. One-time analysis, not a runtime dependency.
3. **TMD coverage closure.** A future skill running DeepTMHMM (or HMMTOP) on `protein.faa` to fill the UniProt TMD-curation gap. Once coverage exceeds ~80%, re-add a step D that flags membrane-class TC calls (TC class 1/2/4/5) where the gene has zero TMDs as `low_confidence`.
4. **Substrate compatibility cross-check.** Cross-reference the substrate of each tier-1 specificity call against KG metabolites; flag mismatches as a downstream analysis report.

## 10 — Acceptance criteria (Phase 1)

A working Phase 1 must satisfy:

1. `/tcdb-diamond` produces `<strain>.tcdb.tsv`, `<strain>.tcdb.calls.json`, and `<strain>.tcdb.skill_summary.json` for every strain with a `protein.faa` in the registry, without error.
2. The skill's stdout output includes a per-strain status table with hit count, tier distribution (T1/T2/T3), and agreement distribution (confirms / refines / extends / conflicts).
3. `<strain>.tcdb.calls.json` is keyed by protein_id and follows the §6.5 schema.
4. The skill is idempotent: re-running without `--force` skips strains whose calls.json already exists; with `--force` it regenerates.
5. `--refresh-tcdb` re-downloads the TCDB FASTA + diamond DB into `~/tools/TCDB/DB/` even if `tcdb.dmnd` already exists.
6. **No regression in existing files** — `gene_annotations_merged.json`, `step2_metabolism_report.json`, and the KG build are unchanged by Phase 1. (Phase 2 will integrate.)
7. Existing tests (`pytest -m "not slow and not kg"` and `pytest -m kg`) continue to pass — Phase 1 should not affect them since it touches only new files.

## 11 — Open questions Phase 1 answers (input to Phase 2)

To be answered after the first Phase 1 run on all strains; these drive the Phase 2 merge-rule design:

- What fraction of eggNOG-TC genes receive a `confirms` diamond tag? (Confidence sanity check.)
- What fraction receive `refines`? (Specificity headline win.)
- What fraction `extends`? (Coverage headline win.)
- What fraction `conflicts`? (Drives whether the v1 union rule needs replacement.)
- Tier 1/2/3 distribution — do the 70% / 40% identity bands produce roughly the expected ratio (~30% / ~40% / ~30%), or are they too strict / too lenient?
- Do `consensus_agreement: 3_part` rejections (top-N hits scattered across families) correlate with multi-domain proteins (suggests they're real ambiguity) or with low-quality genes (suggests we should drop them entirely)?
