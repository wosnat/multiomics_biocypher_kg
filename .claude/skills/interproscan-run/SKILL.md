---
name: interproscan-run
description: Run InterProScan 5 (via the interpro/interproscan Docker image) on each strain's protein.faa to predict protein domains/families across all member databases (Pfam, NCBIfam/TIGRFAM, Hamap, PROSITE, SFLD, PANTHER, Gene3D, SUPERFAMILY, CDD, PRINTS, SMART, …), integrated into InterPro entries with GO term attribution. Emits a faceted per-protein `calls.json` (per-library match rows, InterPro entry rollups, GO — no pathways). Phase 1 — produces inspectable `<strain>.interproscan.calls.json` artifacts; **Phase 2 KG integration is DONE** (`InterproEntry` ontology via `interpro_adapter` + split-out `NcbifamFamily` ontology via `ncbifam_adapter` — see `docs/kg-changes/interpro-multi-ontology.md`). Triggers on "run interproscan", "predict protein domains", "InterPro / domain annotation for all strains", "functional domains for the new strain".
argument-hint: "[--strains <name> ... | --force | --limit N | --threads N | --applications APP,APP | --no-xrefs | --normalize | --prepare-image | --refresh-data]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(docker *), Bash(jq *)
---

# InterProScan-Run Skill

Runs **InterProScan 5.78-109.0** over each configured strain's `protein.faa`
and writes per-protein domain calls. InterProScan is an orchestrator: it runs
every protein against ~15 member databases (Pfam, NCBIfam/TIGRFAM, Hamap,
PROSITE, SFLD, PANTHER, Gene3D, SUPERFAMILY, CDD, PRINTS, SMART, PIRSF,
MobiDBLite, Coils) and maps each hit into a unified **InterPro entry**
(`IPRxxxxxx`) with GO term attribution. This is **Phase 1** — it produces
inspectable `<strain>.interproscan.calls.json` artifacts (keyed by RefSeq WP_
accession) committed under git; KG integration is a separate Phase-2 spec.
Pathway xrefs are **not** retained in this artifact (see Output Schema) —
entry-level pathway lookup lives centrally in `interpro_reference.json`.

Default runs **all default free apps** (the licensed SignalP/TMHMM/Phobius are
deactivated in the free distribution). Restrict per-run with `--applications`.

## Quick Start

```bash
# One-time install (see One-Time Setup) — pull the image, then download data:
uv run python .claude/skills/interproscan-run/run_interproscan.py --prepare-image
uv run python .claude/skills/interproscan-run/run_interproscan.py --refresh-data

# Run all genome strains (skips strains that already have calls.json)
uv run python .claude/skills/interproscan-run/run_interproscan.py

# Run one or more strains by name
uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4
uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4 MIT9313

# Force re-run even if calls.json exists
uv run python .claude/skills/interproscan-run/run_interproscan.py --force

# Restrict to a fast, high-signal subset (overrides the all-apps default)
uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4 \
    --applications Pfam,NCBIFAM,Hamap,PROSITEPatterns,PROSITEProfiles,SFLD

# Smoke test: first 100 proteins of one strain (output is gitignored)
uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4 --limit 100

# GO + pathway xrefs are ON by default; opt out for a smaller/faster artifact
# (pathways are requested from IPS either way but discarded during parsing —
# see Output Schema)
uv run python .claude/skills/interproscan-run/run_interproscan.py --strains MED4 --no-xrefs

# Full batch (multi-day wallclock with all apps) — run detached:
nohup uv run python .claude/skills/interproscan-run/run_interproscan.py \
    > logs/interproscan/batch.log 2>&1 &

# Re-parse cached raw.json into the current faceted calls.json — no Docker,
# no re-scan (e.g. after a parser change). Requires the raw.json to still be
# on disk (gitignored, kept only on the machine that ran the scan).
uv run python .claude/skills/interproscan-run/run_interproscan.py --normalize
uv run python .claude/skills/interproscan-run/run_interproscan.py --normalize --strains MED4
```

`--normalize` needs `<strain>.interproscan.raw.json` locally — it is never
committed. A strain missing its raw.json reports `NO_RAW` and is left
untouched; on a machine that never ran the scan (a fresh checkout, CI, a
teammate's laptop), keep using the committed `<strain>.interproscan.calls.json`
as-is rather than trying to normalize.

## One-Time Setup

Two pieces: the Docker image (software) and the data directory (member-DB
models), downloaded separately and mounted read-only.

**1. Pull the image** (~2–3 G):

```bash
uv run python .claude/skills/interproscan-run/run_interproscan.py --prepare-image
# → docker pull interpro/interproscan:5.78-109.0 ; prints the RepoDigest and exits.
docker images | grep interproscan   # verify
```

**2. Download the data** (6.4 G compressed → extracted under `~/tools/InterProScan/`):

```bash
uv run python .claude/skills/interproscan-run/run_interproscan.py --refresh-data
# Downloads interproscan-data-5.78-109.0.tar.gz + .md5 from EBI, verifies the
# checksum, extracts to $INTERPROSCAN_DATA_DIR/interproscan-5.78-109.0/data,
# and exits. Skips the download if the data dir already exists (use --force).
ls $INTERPROSCAN_DATA_DIR/interproscan-5.78-109.0/data   # verify member-DB subdirs
```

- **`.env` entry (optional):** `INTERPROSCAN_DATA_DIR=~/tools/InterProScan`.
  Defaults to `~/tools/InterProScan` when unset — the canonical outside-the-repo
  location (never committed). Override per-run with `--data-dir`.
- **Disk footprint:** image ~2–3 G; data tarball 6.4 G download; extracted data
  **34 G** (measured for 5.78-109.0). Member-DB subdirs present: antifam, cdd,
  funfam, gene3d, hamap, ncbifam, panther, pfam, phobius, pirsf, pirsr, prints,
  prosite, sfld, smart, superfamily, tmhmm. (phobius/tmhmm data ship but the
  analyses are deactivated in the free distribution; SignalP is absent.)
- **Prereqs:** Docker on PATH (`docker --version`).
- **Portability note:** the runner passes `--user $(id -u):$(id -g)` and mounts
  a writable `--tempdir` under the strain's `interproscan/temp/` so IPS output
  is host-owned (not root). `--disable-precalc` is passed so IPS computes
  matches locally instead of querying the EBI precalc lookup service (novel
  bacterial proteins mostly miss the lookup anyway, and it removes a network
  dependency for reproducible batches).

## What It Does

`run_interproscan.py` loads genome rows from `cyanobacteria_genomes.csv`
(`load_genome_rows`). For each strain it:

1. Resolves `<data_dir>/protein.faa` (RefSeq WP_ accessions as FASTA ids);
   emits a `MISSING_INPUT` status row and continues if absent.
2. Skips the strain if `<data_dir>/interproscan/<strain>.interproscan.calls.json`
   already exists, unless `--force`.
3. Runs the Docker image, mounting the strain dir at `/input` and the data dir
   read-only at `/opt/interproscan/data`:

   ```
   docker run --rm --user UID:GID \
     -v <data_dir>:/input \
     -v <INTERPROSCAN_DATA_DIR>/interproscan-5.78-109.0/data:/opt/interproscan/data \
     interpro/interproscan:5.78-109.0 \
     --input /input/protein.faa \
     --output-file-base /input/interproscan/<strain>.interproscan.raw \
     --formats JSON --disable-precalc --tempdir /input/interproscan/temp \
     --cpu <threads> [--applications <apps>] --goterms --pathways
   ```

   `--goterms` / `--pathways` switch on the GO and pathway xrefs carried by
   integrated InterPro entries (each IMPLIES `-iprlookup`). Both are **local**
   lookups against the bundled entry data, so they stay compatible with
   `--disable-precalc` and add no network dependency and no measurable
   wallclock. Without them IPS emits `goXRefs: []` / `pathwayXRefs: []` on
   every entry — which is exactly what the original 2026-07-23 batch produced.
   `--no-xrefs` restores that older, smaller output. **`--pathways` is still
   requested from IPS here but its output is discarded at parse time** (see
   Output Schema) — only `--goterms` output survives into `calls.json`.

   Full stdout+stderr goes to `logs/interproscan/<strain>.log`.
4. Parses the raw JSON (`multiomics_kg.utils.interproscan.parse_interproscan_json`)
   into a WP_-keyed faceted calls dict and writes `<strain>.interproscan.calls.json`
   + `<strain>.interproscan.skill_summary.json` via `tool_calls_io`. There is no
   `entry_xrefs.json` sidecar any more (see Output Schema).
5. Prints a status table: `strain | n_proteins | n_calls | matches | wallclock_s | status`.

`--normalize` is a separate, faster mode that skips steps 1–3 entirely: it
re-parses the strain's already-cached `<strain>.interproscan.raw.json` through
the current `parse_interproscan_json` + `summarize` and rewrites `calls.json` +
`skill_summary.json` in place (deleting any stale `entry_xrefs.json` sidecar on
success). Use it after a parser change, instead of re-running Docker.

**Stale-output guard.** InterProScan refuses to overwrite an existing output
file — it silently writes `<base>_1.json`, `<base>_2.json`, … instead. The
runner therefore deletes the canonical raw JSON *and* any `_N` leftovers before
each run. Without that, a `--force` re-run re-parses the **previous** run's raw
JSON and reports success with unchanged data (this is how the first
`--goterms` re-run appeared to produce zero GO terms).

`--limit N` runs only the first N sequences of `protein.faa` (written to a
truncated FASTA); outputs get a `.limited_<N>.` infix and are auto-gitignored.

## Output Schema

**Faceted format (2026-08-17 multi-ontology redesign).** This supersedes the
old flat `matches[]` shape — see
`docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`
(which supersedes the output-schema section of the original 2026-07-22 phase-1
spec). Rewritten by `parse_interproscan_json`
(`multiomics_kg/utils/interproscan.py`); regenerate a strain's artifacts from
cached raw output with `--normalize` (see Quick Start) rather than hand-editing.

**`<strain>.interproscan.calls.json`** — top-level dict keyed by WP_ accession.
Proteins processed but matching nothing are kept (`match_count: 0`, all three
facets `{}`) — so a **missing key** means "not processed" and a present key
with `match_count == 0` means "processed, no domain found":

```jsonc
{
  "<WP_accession>": {
    "md5": "<str>",
    "match_count": 0,
    "libraries": {
      "<LIBRARY>": [
        {"accession": "<str>", "name": "<str|null>", "ipr": "<str|null>",
         "start": 0, "end": 0, "evalue": 0.0, "score": 0.0}
      ]
    },
    "interpro_entries": {
      "<IPR*>": {"type": "<str>", "libraries": ["<str>", "..."],
                 "match_count": 0, "start": 0, "end": 0,
                 "evalue": 0.0, "evalue_library": "<str|null>"}
    },
    "go_terms": {"<GO:NNNNNNN>": ["<IPR*>", "..."]}
  }
}
```

| Field | Notes |
|---|---|
| `md5` | InterProScan sequence MD5 (dedup key) |
| `match_count` | count of (signature × location) rows across all libraries; `0` = sentinel (no domain found) |
| `libraries` | **sparse** — a key exists only for member DBs that matched. One row per (signature × location); sorted by (start, accession). Each row's `accession` is version-stripped (`NF002735.2` → `NF002735`); `ipr` is that row's InterPro entry, or `null` if the signature isn't integrated |
| `interpro_entries` | rollup keyed by `IPR*`, one entry per distinct InterPro id the protein's matches carry; **sparse** (`{}` if no matches integrated) |
| `go_terms` | GO id → sorted list of the `IPR*` entries that donated it (attribution, not a flat list); **sparse** (`{}` if `--no-xrefs` or nothing carries GO) |

**Rules:**
- **Facets are sparse.** A `libraries` key exists only if that DB produced a
  hit; `interpro_entries` / `go_terms` are `{}`, not omitted, when empty.
- **`ipr: null`** on a `libraries` row means that signature is not integrated
  into any InterPro entry — the member-DB hit stands alone (common for
  PANTHER/PIRSF/MobiDBLite and some NCBIfam rows).
- **No `pathways` field anywhere** — not on the record, not on a `libraries`
  row, not on an `interpro_entries` entry. Pathway xrefs are entry-level and
  fully derivable from the central
  `cache/data/interpro/interpro_reference.json` (built once by `prepare_data`
  step 9); duplicating them per-strain per-protein was pure bloat (~162K
  Reactome refs on this proteome, mostly species-projection noise) with zero
  KEGG content — InterPro carries no KEGG pathway xrefs at all. For KEGG
  pathways, keep using the existing KO layer (`Gene_has_kegg_ko` →
  `KeggTerm`), which this redesign doesn't touch.
- **`interpro_entries` has no cross-library `score`** — "count, don't
  combine": each entry tracks the **min `evalue`** across its contributing
  rows plus `evalue_library` (which library reported it), never a synthesized
  blended score. Cross-library corroboration is read from
  `size(libraries)` (independent model families agreeing) vs. `match_count`
  (total rows, including within-library repeats), not from a combined number.
- **`<strain>.interproscan.entry_xrefs.json` no longer exists.** The old
  per-strain sidecar (`{IPR*: {go_terms, pathways}}`) is gone — GO provenance
  is now inline in each record's `go_terms` facet above, and entry-level
  lookups (names, types, hierarchy, GO, pathways, EC, CAZy) live centrally in
  `cache/data/interpro/interpro_reference.json`, not per strain. `--normalize`
  deletes any stale `entry_xrefs.json` it finds once it successfully rewrites
  a strain's `calls.json`.

**`<strain>.interproscan.skill_summary.json`** — per-strain QC:
`strain`, `tool_version`, `image_digest`, `applications`, `input_proteins`,
`calls_made`, `proteins_no_match`, `parse_failures`, `total_matches`,
`interpro_integrated_matches`, `distribution` (matches per member DB),
`sentinel_rate`, `wallclock_s`, `proteins_with_go_terms`, `distinct_go_terms`,
plus `xrefs_requested` (bool) when the caller passes it (the scan path does;
`--normalize` does not, since it isn't re-requesting anything from IPS). The
pathway QC fields — `proteins_with_pathways`, `distinct_pathways`,
`pathway_databases` — are **gone**; everything else is unchanged from the
pre-redesign shape.

## QC

### Per-strain QC — what to inspect in `<strain>.interproscan.skill_summary.json`

- `parse_failures == 0` — non-zero means the raw JSON drifted from what the
  parser expects (IPS version bump?).
- `calls_made == input_proteins` — every input protein should appear in the
  output; a shortfall means IPS dropped sequences (non-standard residues, empty
  seqs). Check `logs/interproscan/<strain>.log`.
- `sentinel_rate` (proteins with no match): calibrated to the 42-strain batch
  = **0.06–0.21** (mean 0.12). *Prochlorococcus* runs high (~0.15 — streamlined
  genomes carry more uncharacterized ORFs); *Alteromonas*/heterotrophs run low
  (~0.08). Flag anything **> 0.30** — that usually means the data dir is missing
  member-DB subdirs or the wrong `--applications` was passed.
- `interpro_integrated_matches / total_matches` — the fraction of hits that map
  to an `IPR` entry; expect a healthy majority. A collapse signals a data-dir
  version mismatch.
- `distribution` — member-DB hit counts. Pfam + NCBIfam should dominate; if a
  big DB (PANTHER/Gene3D) shows `0`, its data subdir wasn't downloaded.
- **GO coverage** — `xrefs_requested: true` with `distinct_go_terms: 0` means
  the flags were passed but the output is stale (see the stale-output guard) or
  the entry data is missing. On the MED4 smoke, 71/100 proteins carried a GO
  term; expect `proteins_with_go_terms / calls_made` in the **0.5–0.8** band.
  Both counters read `0` on a `--no-xrefs` run — that is correct, not a
  failure. (Pathway coverage is no longer tracked here at all — see Output
  Schema.)

### Cross-strain QC narrative

Observed across the 42-strain batch (2026-07-23). These are raw-scan facts
(matches, integration rate, sentinel rate) — the 2026-08-17 `--normalize`
re-parse of all 42 strains into the faceted schema re-read the same cached
`raw.json` files (no Docker re-run), so `calls_made`, `total_matches`, and
`interpro_integrated_matches` are unchanged per strain by construction; only
the on-disk shape changed. Pathway-specific numbers elsewhere in this doc
predate that redesign and no longer apply to the current artifact — see
Output Schema.

- **Domain density** (matches/protein) is tight and lifestyle-consistent:
  *Prochlorococcus* 7.45, *Synechococcus* 7.68, *Alteromonas* 8.86,
  other heterotrophs 8.71 — the larger heterotroph proteomes carry more
  multi-domain proteins, exactly as expected.
- **InterPro-integration fraction** is stable at **66–69%** across every genus
  (67.0% overall) — a member-DB-version-mismatch collapse would show up as a
  genus dropping well below this band.
- **Sentinel rate** cleanly tracks genome novelty: *Prochlorococcus* 15.3% (the
  streamlined-genome ORFs that resist all curated models) vs heterotrophs ~7.7%.
  No strain was an outlier; `parse_failures == 0` for all 42; `calls_made ==
  input_proteins` for all 42.
- **Cross-tool sanity:** Pfam is the single largest contributor (15.6% of all
  matches), consistent with the KG's existing eggNOG-derived Pfam layer — but
  IPS adds coordinates + the InterPro id + the genuinely-new NCBIfam (6.9%),
  Hamap (2.5%), PROSITE (6.8%), SFLD (0.4%) and structural Gene3D/SUPERFAMILY
  (26.5% combined) evidence that eggNOG does not provide.

### Spot checks

Verified against the committed, faceted `calls.json` (2026-08-17 re-normalize
of all 42 strains). Checks 1–3 confirm the per-row/per-entry/per-GO-term
shapes on one protein; checks 4–5 are a cross-strain biological sanity check
(TIGR00198/TIGR00357 presence — mirrors the katG Black Queen Hypothesis check
used elsewhere in this repo for *Prochlorococcus* genome streamlining).

| # | Strain / Protein | Check | Expected |
|---|---|---|---|
| 1 | MED4 `WP_002805124.1` | `libraries.PFAM[0]` accession + `ipr` | `{"accession": "PF02532", "ipr": "IPR003686"}` |
| 2 | MED4 `WP_002805124.1` | `libraries.NCBIFAM[0].ipr` | `null` (NCBIfam hit not integrated into an InterPro entry) |
| 3 | MED4 `WP_002805124.1` | `go_terms."GO:0015979"` | includes `"IPR003686"` |
| 4 | EZ55 `WP_156086936.1` / MIT1002 `WP_014977393.1` vs. MED4 (whole file) | `libraries.NCBIFAM[].accession` | `TIGR00198` present in EZ55 + MIT1002, **absent** across all of MED4 |
| 5 | MED4 `WP_011131653.1` | `libraries.NCBIFAM[].accession` | includes `TIGR00357` (MsrB) |

```bash
# 1. PFAM row + its InterPro attribution:
jq '."WP_002805124.1".libraries.PFAM[0] | {accession, ipr}' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# -> {"accession": "PF02532", "ipr": "IPR003686"}

# 2. NCBIfam row that isn't integrated:
jq '."WP_002805124.1".libraries.NCBIFAM[0].ipr' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# -> null

# 3. GO term carries entry attribution:
jq '."WP_002805124.1".go_terms."GO:0015979"' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# -> includes "IPR003686"

# 4. Cross-strain TIGR00198 (present in EZ55 + MIT1002, absent from MED4):
jq -r '."WP_156086936.1".libraries.NCBIFAM[].accession' \
  cache/data/Alteromonas/genomes/EZ55/interproscan/EZ55.interproscan.calls.json
jq -r '."WP_014977393.1".libraries.NCBIFAM[].accession' \
  cache/data/Alteromonas/genomes/MIT1002/interproscan/MIT1002.interproscan.calls.json
jq -r '[.[] | .libraries.NCBIFAM[]?.accession] | index("TIGR00198")' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# -> "TIGR00198" in both EZ55/MIT1002 output; null (absent) for MED4

# 5. MED4 carries TIGR00357 (MsrB):
jq -r '."WP_011131653.1".libraries.NCBIFAM[].accession' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# -> includes "TIGR00357"
```

## Observed batch results (42-strain run, 2026-07-23; re-run with xrefs 2026-08-09)

*These totals predate the 2026-08-17 faceted-schema redesign, which
`--normalize` applied to all 42 strains against the same cached `raw.json`
inputs — so `input_proteins` / `total_matches` / integration % are unchanged
(see Task 3, `.superpowers/sdd/2026-08-17-interpro-multi-ontology-redesign/task-3-report.md`).
The pathway numbers below are historical only — the current `calls.json` no
longer carries a pathway facet (see Output Schema); the GO numbers still
apply via the `go_terms` facet.*

**120,343 proteins** processed (all appeared in output; `parse_failures = 0`).
**986,526 matches**, of which **67.0% integrate to an InterPro entry**. All 42
strains `OK`, no failures. Total compute ~26.6 CPU-wallclock hours.

**The 2026-08-09 `--goterms --pathways` re-run reproduced all three of those
totals exactly** — 120,343 proteins / 986,526 matches / 67.0% integrated, with
`parse_failures = 0`, `calls_made == input_proteins` on all 42, and sentinel
rates spanning 0.064–0.211 (inside the documented band). The xref flags add
annotation without perturbing matching at all; that identity is the regression
check to repeat after any future IPS version bump. Compute 28.0 h (+5%).

New xref coverage from that re-run:

| Metric | Value |
|---|---|
| proteins carrying ≥1 GO term | 69,759 / 120,343 (**58.0%**) |
| proteins carrying ≥1 pathway | 76,377 / 120,343 (63.5%) |
| pathway xrefs by source DB | Reactome 10,985,362 · MetaCyc 275,755 |

Per-strain GO coverage is tight (0.52–0.61) and tracks lifestyle the same way
domain density does — heterotrophs (EZ55, HOT1A3, MIT1002) ~0.61, streamlined
*Prochlorococcus* ~0.53–0.58. The Reactome:MetaCyc ratio (~40:1) is the
species-projection effect described under Output Schema, already collapsed to
stable ids; it is not a QC problem.

Cross-strain member-DB distribution (matches, % of all matches):

| Member DB | Matches | % |
|---|---:|---:|
| PFAM | 153,669 | 15.6 |
| GENE3D | 145,734 | 14.8 |
| SUPERFAMILY | 115,911 | 11.7 |
| PANTHER | 82,894 | 8.4 |
| PRINTS | 77,138 | 7.8 |
| NCBIFAM | 67,918 | 6.9 |
| PIRSR | 61,822 | 6.3 |
| CDD | 56,071 | 5.7 |
| PROSITE_PROFILES | 41,919 | 4.2 |
| FUNFAM | 37,487 | 3.8 |
| SMART | 33,245 | 3.4 |
| MOBIDB_LITE | 28,440 | 2.9 |
| PROSITE_PATTERNS | 25,530 | 2.6 |
| HAMAP | 24,532 | 2.5 |
| COILS | 15,192 | 1.5 |
| PIRSF | 15,092 | 1.5 |
| SFLD | 3,931 | 0.4 |
| ANTIFAM | 1 | 0.0 |

Per-strain wallclock: **21 min** (SB, 1,839 proteins) → **67 min** (KT2440,
5,452 proteins), scaling ~linearly with proteome size at all apps.

## Phase 2 — DONE (KG integration)

Integrated via `/integrate-a-tool` (2026-07-26). The calls.json now flow into the
KG on **both** surfaces the earlier sketch weighed: an `interproscan` logical
source merged into `gene_annotations_merged.json` (light `interpro_entries` list),
**and** an `InterproEntry` ontology (`interpro_adapter`) — hierarchical nodes +
scored `Gene_has_interpro_entry` edges (coords/evalue/`evalue_library`/libraries — **no `score`**; 0 of
397,342 edges carry one as of the 2026-08-17 multi-ontology redesign, which dropped the synthesized
cross-library score in favor of `evalue_library` naming which member DB reported the min evalue) +
`Interpro_entry_is_a_interpro_entry` hierarchy + `Pfam_in_interpro_entry` bridge. That same redesign also
split NCBIfam/TIGRFAM out into its own `NcbifamFamily` ontology (`ncbifam_adapter`) — its
`Gene_has_ncbifam_family` edges DO keep both `evalue` and `score` (single-library HMMER scale, not a
cross-library rollup). Node names/types/hierarchy come from `prepare_data` **step 9**
(`build_interpro_reference.py` + `build_ncbifam_reference.py` → committed
`cache/data/interpro/interpro_reference.json` + `cache/data/ncbifam/ncbifam_reference.json`).

- **Integration design:** `docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md` (original);
  `docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md` (NCBIfam split + faceted format)
- **What-changed / MCP contract:** `docs/kg-changes/interpro-multi-ontology.md` (current — supersedes
  `docs/kg-changes/interproscan-extension.md` and `docs/kg-changes/interpro-two-layer.md`, both kept for history)

Re-running `/interproscan-run` on a new strain + `prepare_data.sh --steps 2` (and
a KG rebuild) is all that's needed for that strain's domains to appear in the graph.

**GO/pathway xrefs (2026-08-07):** no longer deferred — `--goterms --pathways`
are on by default and the 42-strain batch was re-run to populate them. The
adapter does not consume them yet; wiring InterPro→GO into the graph is a
follow-up (it would need a decision on how these GO terms relate to the existing
eggNOG/UniProt-derived `Gene_involved_in_biological_process` &co. edges — most
likely a new evidence source on those edges rather than a parallel edge type).

## Workflow When Invoked

1. Verify Docker: `docker --version`.
2. If `docker images | grep interproscan` is empty → run `--prepare-image`.
   If `$INTERPROSCAN_DATA_DIR/interproscan-5.78-109.0/data` is absent → run
   `--refresh-data`. Both exit after setup; they never start a batch.
3. Run the batch (single strain, or all). Inspect the status table.
4. Verify the spot checks with the `jq` one-liners above.
5. Chase any `FAILED` row via `logs/interproscan/<strain>.log`.
