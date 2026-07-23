---
name: interproscan-run
description: Run InterProScan 5 (via the interpro/interproscan Docker image) on each strain's protein.faa to predict protein domains/families across all member databases (Pfam, NCBIfam/TIGRFAM, Hamap, PROSITE, SFLD, PANTHER, Gene3D, SUPERFAMILY, CDD, PRINTS, SMART, …), integrated into InterPro entries with GO + pathway xrefs. Emits per-protein domain calls. Phase 1 — produces inspectable `<strain>.interproscan.calls.json` artifacts; KG integration deferred to Phase 2. Triggers on "run interproscan", "predict protein domains", "InterPro / domain annotation for all strains", "functional domains for the new strain".
argument-hint: "[--strains <name> ... | --force | --limit N | --threads N | --applications APP,APP | --prepare-image | --refresh-data]"
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(docker *), Bash(jq *)
---

# InterProScan-Run Skill

Runs **InterProScan 5.78-109.0** over each configured strain's `protein.faa`
and writes per-protein domain calls. InterProScan is an orchestrator: it runs
every protein against ~15 member databases (Pfam, NCBIfam/TIGRFAM, Hamap,
PROSITE, SFLD, PANTHER, Gene3D, SUPERFAMILY, CDD, PRINTS, SMART, PIRSF,
MobiDBLite, Coils) and maps each hit into a unified **InterPro entry**
(`IPRxxxxxx`) with GO + pathway xrefs. This is **Phase 1** — it produces
inspectable `<strain>.interproscan.calls.json` artifacts (keyed by RefSeq WP_
accession) committed under git; KG integration is a separate Phase-2 spec.

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

# Full batch (multi-day wallclock with all apps) — run detached:
nohup uv run python .claude/skills/interproscan-run/run_interproscan.py \
    > logs/interproscan/batch.log 2>&1 &
```

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
     --cpu <threads> [--applications <apps>]
   ```

   Full stdout+stderr goes to `logs/interproscan/<strain>.log`.
4. Parses the raw JSON (`multiomics_kg.utils.interproscan.parse_interproscan_json`)
   into a WP_-keyed calls dict and writes `<strain>.interproscan.calls.json` +
   `<strain>.interproscan.skill_summary.json` via `tool_calls_io`.
5. Prints a status table: `strain | n_proteins | n_calls | matches | wallclock_s | status`.

`--limit N` runs only the first N sequences of `protein.faa` (written to a
truncated FASTA); outputs get a `.limited_<N>.` infix and are auto-gitignored.

## Output Schema

**`<strain>.interproscan.calls.json`** — top-level dict keyed by WP_ accession.
Proteins processed but matching nothing are kept with an empty `matches` list,
so a missing key means "not processed" and `match_count == 0` means "no domain
found":

| Field | Type | Notes |
|---|---|---|
| `md5` | str | InterProScan sequence MD5 (dedup key) |
| `match_count` | int | number of (signature × location) records; `0` = sentinel |
| `interpro_entries` | str[] | distinct sorted `IPR*` ids |
| `go_terms` | str[] | distinct sorted `GO:*` from integrated entries |
| `pathways` | str[] | distinct sorted `DB:id` (KEGG/Reactome/MetaCyc) |
| `libraries` | str[] | distinct member DBs that produced a hit |
| `matches` | obj[] | one per (signature × location), sorted by (start, evalue, accession) |

Each `matches[]` record: `library`, `signature_accession`,
`signature_description`, `interpro_accession` (null if unintegrated),
`interpro_description` (null), `interpro_type` (null), `start`, `end`,
`evalue` (null for pattern/profile hits), `score` (null when N/A),
`go_terms` (str[]), `pathways` (str[]).

**`<strain>.interproscan.skill_summary.json`** — per-strain QC:
`strain`, `tool_version`, `image_digest`, `applications`, `input_proteins`,
`calls_made`, `proteins_no_match`, `parse_failures`, `total_matches`,
`interpro_integrated_matches`, `distribution` (matches per member DB),
`sentinel_rate`, `wallclock_s`.

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

### Cross-strain QC narrative

Observed across the 42-strain batch (2026-07-23):

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

Highly conserved MED4 proteins that must produce a canonical Pfam hit (all
within the first 3 sequences, so the `--limit 100` smoke covers them):

| Strain | Protein ID | Expected | Why this is the ground truth |
|---|---|---|---|
| MED4 | WP_002805854.1 | Pfam `PF00016` (RuBisCO_large) → `IPR000685` | rbcL, RuBisCO form I large subunit — universally conserved. ✅ verified (smoke, 2026-07-22) |
| MED4 | WP_002805169.1 | Pfam `PF00137` (ATP-synt_C) → `IPR000454` | atpH, ATP synthase F0 subunit c — universally conserved. ✅ verified (smoke, 2026-07-22) |
| WH8102 (Syn) | WP_011128579.1 | Pfam `PF00016`+`PF02788` → `IPR000685` | rbcL, RuBisCO large subunit. ✅ verified (batch, 2026-07-23) |
| KT2440 (Pseudomonas) | WP_010952474.1 | Pfam `PF00118` (TCP-1/cpn60) → `IPR002423` | groEL 60 kDa chaperonin — universally conserved. ✅ verified (batch, 2026-07-23) |

```bash
# rbcL should carry a PFAM RuBisCO_large signature:
jq -r '."WP_002805854.1".matches[] | select(.library=="PFAM") | .signature_accession' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# Expected to include: PF00016

# ...and integrate to an InterPro family entry:
jq -r '."WP_002805854.1".interpro_entries[]' \
  cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# Expected to include: IPR000685
```

## Observed batch results (42-strain run, 2026-07-23)

**120,343 proteins** processed (all appeared in output; `parse_failures = 0`).
**986,526 matches**, of which **67.0% integrate to an InterPro entry**. All 42
strains `OK`, no failures. Total compute ~26.6 CPU-wallclock hours.

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

## Phase 2 (Future)

Deferred to a separate spec — see
`docs/superpowers/specs/2026-07-22-interproscan-domains-design.md` (Phase 2
sketch). Phase-1 artifacts sit in each strain's `interproscan/` cache dir for
inspection and are **not** wired into `gene_annotations_merged.json` or any KG
adapter yet. Two candidate surfaces: an `interproscan` logical source merged
into `gene_annotations_merged.json`, or a new `InterProEntry` node type +
`Gene_has_interpro_entry` edges via `/integrate-a-tool`.

## Workflow When Invoked

1. Verify Docker: `docker --version`.
2. If `docker images | grep interproscan` is empty → run `--prepare-image`.
   If `$INTERPROSCAN_DATA_DIR/interproscan-5.78-109.0/data` is absent → run
   `--refresh-data`. Both exit after setup; they never start a batch.
3. Run the batch (single strain, or all). Inspect the status table.
4. Verify the spot checks with the `jq` one-liners above.
5. Chase any `FAILED` row via `logs/interproscan/<strain>.log`.
