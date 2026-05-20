---
name: add-a-strain
description: End-to-end onboarding for a brand-new bacterial strain in the KG AND re-deploy of an existing strain's v2 gene-ID mapping after a paperconfig fix. Identify the NCBI/Cyanorak assembly from the paper that motivated the addition, register the strain in all three hardcoded sites, run prepare_data, run every per-strain tool (eggnog, psortb, tcdb-diamond, signalp, …) in the background, verify gene-ID resolution against the new paper's CSVs, rebuild the Docker KG, refresh snapshot/unit tests, then loop back through prepare_data so tool outputs land in `gene_annotations_merged.json`. Use this whenever you're integrating a paper that needs a strain that's not in `cyanobacteria_genomes.csv` yet, when the user says "add SS120" / "onboard <strain>" / "we need <strain> for <paper>", OR when you need to re-deploy v2 gene-ID mapping for an existing strain ("deploy <strain>", "re-resolve <strain> papers" — skip the brand-new-strain prefix and jump to step 4).
argument-hint: <strain-name> [--paper "<Author Year>"] [--accession <GCF...>] [--redeploy]
user-invocable: true
allowed-tools: Read, Edit, Grep, Glob, Bash(uv *), Bash(python *), Bash(bash *), Bash(docker *), Bash(pytest *), Bash(git *), WebFetch
---

# Add a Strain

Brand-new strain → fully-integrated KG row, with every per-strain tool run and
every test/snapshot updated. Also covers re-deploying an existing strain's v2
gene-ID mapping after a paperconfig fix (skip steps 0–3, jump to step 4).

The skill is opinionated about ordering because several steps are slow (Docker
rebuild, eggNOG-mapper, PSORTb) and cheap to parallelise but expensive to redo.

## When to use

- A new paper requires expression data on a strain not yet in `cyanobacteria_genomes.csv`.
- The user names a strain ("SS120", "BL107", "WH7805", "AD45", "Pseudohoeflea") and asks to onboard / add / integrate it.
- An existing strain needs a re-deploy after a paperconfig change (new `id_translation`, fixed `id_type`, new annotation source) — skip steps 0–3.
- Reference-proteome-match strains for community/MarRef-style proteomics: set `organism_type=reference_proteome_match` in the CSV row and skip the Cyanorak fields.

Brand-new strain additions in 2026-04 followed this workflow — see
`docs/superpowers/specs/2026-04-14-csv-ready-papers-batch-design.md` for the
SS120 / BL107 / HP15 / Alt_MarRef worked example, and `plans/mit1002_deploy.md`
+ `plans/ez55_deploy.md` for the cross-assembly cases.

## Multi-strain batches

The workflow handles N strains in one pass — the placeholder `<NEW_STRAIN>`
below stands for either a single name or a space-separated list. Diffs from
the single-strain flow:

| Step | Multi-strain diff |
|---|---|
| 0 | Identify each accession individually — only step that's per-strain by nature. |
| 1 | Append N rows to `cyanobacteria_genomes.csv`; add N alias bundles to `ORGANISM_TO_GENOME_DIR`; add N entries to `CANONICAL_GENOMIC_ORGANISMS`. One pytest run validates all of them. |
| 2 | `bash scripts/prepare_data.sh --strains S1 S2 … SN --steps 0 1 2 3 4 5 6 7` — accepts a space-separated list. |
| 3 | **Drop the `--strain` flag entirely** in each runner. eggnog/psortb/signalp/tcdb-diamond all default to iterating `cyanobacteria_genomes.csv` and skipping strains that already have output — so the new strains run, the old ones no-op. The Wave 2 `wait $EGGNOG_PID` gate still works since each runner only exits after its full CSV walk finishes. See the multi-strain block inline in Step 3. |
| 4–9 | One snapshot, one Docker rebuild, one `/check-gene-ids` and one snapshot-compare cover the whole batch. |
| 10 | `bash scripts/prepare_data.sh --strains S1 S2 … SN --steps 1 2 3 4 5 6 7 --force`. |
| 11 | Final sweep is unchanged. |

Recent worked example: 4 strains (SS120 + BL107 + HP15 + Alt_MarRef) added
together in 2026-04 — see the csv-ready-papers spec linked above.

## Quick map

| Stage | File / command | Brand-new only? |
|---|---|---|
| Identify accession | NCBI Datasets + Cyanorak organism tables | yes |
| Register site 1 | `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` | yes |
| Register site 2 | `multiomics_kg/utils/gene_id_utils.py` (`ORGANISM_TO_GENOME_DIR`) | yes |
| Register site 3 | `scripts/validate_paperconfig.py` (`CANONICAL_GENOMIC_ORGANISMS`) | yes |
| Initial download + ID resolution | `bash scripts/prepare_data.sh --strains <S> --steps 0 1 2 3 4 5 6 7` | yes |
| Per-strain tools (background) | `eggnog-run`, `psortb-run`, `tcdb-diamond`, `signalp-run`, … | yes (rerun on paperconfig change only if id_translation changes) |
| Snapshot KG | `omics-edge-snapshot --save before_<S>` | both |
| Docker rebuild | `docker compose up -d --build` | both |
| Verify match rates | `/check-gene-ids` (uses Docker import report) | both |
| Compare snapshot | `omics-edge-snapshot --compare before_<S>` | both |
| KG validity + tests | `pytest -m kg`, regenerate `snapshot_data.json` | both |
| Loop back | re-run `prepare_data.sh --strains <S> --steps 1 2 3 4 5 6 7 --force` once tool outputs land | brand-new |

## Step 0 — Identify the assembly (brand-new only)

Read the paper's Methods + supplementary tables and locate:

1. The exact strain designation as written by the authors (e.g., `Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)`).
2. The NCBI assembly accession the paper used. Cross-check against:
   - **NCBI Datasets** — `https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=<taxid>`
   - **Cyanorak organism table** — `data/Cyanorak Organism Table prochlorococcus.csv` / `synechococcus.csv` for the cyanorak short name, clade, and the legacy nucleotide accession
   - **Paper's UniProt/RefSeq dump** (if cited) — e.g., Domínguez 2017 explicitly cites taxid 167539

If multiple assemblies exist for the same strain (MIT1002 has a draft
`GCF_001077695.1` and a complete `GCF_901457835.2`), pick the one the **paper
used** and document the choice — switching assemblies later requires
re-resolving every paper that references the strain.

Heterotrophs / reference-proteome-match strains skip Cyanorak — leave
`cyanorak_organism` and `clade` blank.

## Step 1 — Register the strain in the three hardcoded sites (brand-new only)

Until `genomes_registry.py` lands (deferred; see "Open follow-ups"), each
addition needs three edits.

### Site 1 — `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`

Schema: `ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade,preferred_name,organism_type,reference_database,reference_proteome`

Three flavors:
- **Genome strain with Cyanorak** (Pro / Syn) — both `cyanorak_organism` and `clade` populated; `organism_type=genome_strain`.
- **Genome strain without Cyanorak** (Alteromonas, Shewanella, Pseudomonas, Marinobacter as a strain) — `cyanorak_organism` and `clade` blank.
- **Reference-proteome-match** (community proteomics referencing MarRef etc.) — `organism_type=reference_proteome_match`, `reference_database` and `reference_proteome` populated. See `docs/kg-changes/reference-proteome-match-organisms.md`.

Example rows from the recent batch:

```csv
GCF_000007925.1,Pro_SS120,167539,SS120,cache/data/Prochlorococcus/genomes/SS120/,LLII,Prochlorococcus marinus subsp. marinus CCMP1375 (SS120),genome_strain,,
GCF_000153805.1,Syn_BL107,313625,BL107,cache/data/Synechococcus/genomes/BL107/,5.1 IV,Synechococcus sp. BL107,genome_strain,,
GCF_000166295.1,,225937,HP15,cache/data/Alteromonas/genomes/HP15/,,Marinobacter (MarRef v6),reference_proteome_match,MarRef v6,GCF_000166295.1
```

### Site 2 — `multiomics_kg/utils/gene_id_utils.py`

Add lowercase entries to `ORGANISM_TO_GENOME_DIR` for the canonical
`preferred_name` **AND** every alias variant a paper might use (with/without
subspecies, with/without strain prefix). Single typo = silent unresolved-IDs.

```python
"prochlorococcus marinus subsp. marinus ccmp1375 (ss120)": "cache/data/Prochlorococcus/genomes/SS120",
"prochlorococcus marinus ss120":                          "cache/data/Prochlorococcus/genomes/SS120",
"prochlorococcus ss120":                                  "cache/data/Prochlorococcus/genomes/SS120",
```

For heterotrophs include both the with-species and without-species form —
papers frequently shorten "Alteromonas macleodii MIT1002" to "Alteromonas
MIT1002".

### Site 3 — `scripts/validate_paperconfig.py`

Add the **exact** `preferred_name` string from site 1 to
`CANONICAL_GENOMIC_ORGANISMS`. Case-sensitive, byte-for-byte.

```python
"Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
```

These three sites transitively cover everything:
- `.claude/skills/check-gene-ids/check_gene_ids.py` imports `ORGANISM_TO_GENOME_DIR` → updated automatically
- `tests/test_gene_id_utils.py` iterates the dict → updated automatically
- `tests/test_paperconfig_validation.py` imports `CANONICAL_GENOMIC_ORGANISMS` → updated automatically

Run the registration tests:

```bash
uv run pytest tests/test_gene_id_utils.py tests/test_paperconfig_validation.py -v
```

## Step 2 — Download + per-strain annotation tables + gene-ID resolution (brand-new only)

### Pre-flight: configure paperconfig if the paper uses non-standard IDs

prepare_data step 3 harvests `id_translation` and `annotation_gff` entries from
every paperconfig — so set them up **before** running prepare_data. When the
paper's IDs aren't standard NCBI/Cyanorak locus_tags, add to the paper's
`paperconfig.yaml`:

**Annotation CSV (`id_translation`):**

```yaml
id_translation_pro_9312:
  type: id_translation
  filename: "data/Prochlorococcus/papers_and_supp/<Author Year>/pro_9312_anot.csv"
  organism: "Prochlorococcus MIT9312"
  id_columns:
    - column: "GID"
      id_type: old_locus_tag      # Tier 1 — e.g. PMT9312_1938
    - column: "SYMBOL"
      id_type: gene_name          # Tier 3
    - column: "uniprot_acc"
      id_type: uniprot_accession  # Tier 2
  product_columns:
    - column: "PRODUCT"
```

**GFF annotation file (`annotation_gff`):**

```yaml
annotation_gff_mit9312:
  type: annotation_gff
  filename: "data/Prochlorococcus/papers_and_supp/<Author Year>/annotation.gff"
  organism: "Prochlorococcus MIT9312"
```

The GFF `locus_tag` attribute → `locus_tag_ncbi` (Tier 1), `old_locus_tag` →
`old_locus_tag` (Tier 1), `protein_id` → `protein_id_refseq` (Tier 2).

For paper IDs that come from a different assembly entirely (draft genome,
IMG/RAST/custom annotation), use **cross-assembly protein sequence bridging** —
see [references/cross_assembly_bridging.md](references/cross_assembly_bridging.md).

#### Check for legacy `_with_locus_tag.csv` workarounds

Earlier strains were fixed using the `/fix-gene-ids` skill, which created
`_with_locus_tag.csv` copies with an added `locus_tag` column and changed the
paperconfig to `name_col: "locus_tag"`. With v2 mapping, this workaround is
unnecessary — the resolver handles mixed ID columns natively.

**Before deploying a strain**, check if its paperconfig was previously modified
by fix-gene-ids:

1. Check if `filename` points to a `_with_locus_tag.csv` file
2. Check if `name_col` was changed to `"locus_tag"` from the original column
3. Check if `id_columns` references a `locus_tag` column that doesn't exist in the original CSV

**To revert**: use `git show <commit>:<paperconfig path>` to find the original,
then change `filename`, `name_col`, and remove the synthetic `id_columns`
entries. Keep valid additions like `prefiltered`, `pvalue_threshold`,
`product_columns`.

### Run prepare_data

```bash
bash scripts/prepare_data.sh --strains <NEW_STRAIN> --steps 0 1 2 3 4 5 6 7
```

| Step | Output |
|---|---|
| 0.1 | NCBI GFF + protein.faa + GBFF → `cache/data/<Organism>/genomes/<Strain>/` |
| 0.2 | Cyanorak GFF + GBK (only if `cyanorak_organism` set; server may throttle — `--skip-cyanorak` resumes) |
| 0.3 | UniProt by taxid → `cache/data/<org_group>/uniprot/<taxid>/` |
| 0.5 | `gene_mapping.csv` (locus_tag ↔ protein_id ↔ gene names) |
| 1 | `protein_annotations.json` per taxid |
| 2 | `gene_annotations_merged.json` per strain |
| 3 | `gene_id_mapping.json` v2 (3-tier) |
| 4 | `<paper>/<csv>_resolved.csv` for any paper already using the strain (gene-ID resolution) |
| 5 | `cache/data/eggnog/og_descriptions.json` (OG description cache) |
| 6 | `cache/data/kegg/kegg_data.json` + `cache/data/tcdb/tcdb_pruned.json` + `cache/data/metabolomics/metabolite_id_mapping.json` (KEGG/TCDB pruning + metabolite alias index) |
| 7 | `<paper>/<csv>_resolved.csv` for `metabolite_assays_table` entries (metabolite-ID resolution) |

Logs land in `logs/prepare_data_step{N}.log`. Cyanorak is the flakiest server in
the chain — if step 0.2 hangs, `--skip-cyanorak` and retry later.

## Step 3 — Per-strain tools (background, parallel) (brand-new only)

Every integrated tool with a Phase-1 SKILL must be run for the new strain so
its output is available before Phase-2 integrations re-merge into
`gene_annotations_merged.json`. Kick these off as background jobs and proceed
to Step 4 immediately — they don't depend on each other, and the rest of the
workflow (snapshot, Docker rebuild, verification) runs in parallel with them.
Step 10's loop-back is the rendezvous point where you wait for them to finish.

```bash
mkdir -p logs/{eggnog,psortb,tcdb,signalp}

# Wave 1 — independent tools, run in parallel:

# eggNOG — functional annotation (read by prepare_data step 2 when present;
#          tcdb-diamond also reads its output → must finish before Wave 2)
nohup uv run python .claude/skills/eggnog-run/run_eggnog.py --strain <NEW_STRAIN> > logs/eggnog/<NEW_STRAIN>.log 2>&1 &
EGGNOG_PID=$!

# PSORTb — Gram-negative subcellular localization
nohup uv run python .claude/skills/psortb-run/run_psortb.py --strain <NEW_STRAIN> > logs/psortb/<NEW_STRAIN>.log 2>&1 &

# SignalP — signal peptide prediction
nohup uv run python .claude/skills/signalp-run/run_signalp.py --strain <NEW_STRAIN> > logs/signalp/<NEW_STRAIN>.log 2>&1 &

# Wave 2 — depends on eggNOG (reads <strain>.emapper.annotations for egn_agreement
#          + gene_annotations_merged.json for pfam_agreement). Subshell waits for
#          eggNOG to finish, then launches tcdb — the outer shell doesn't block,
#          so you can proceed to Step 4 immediately.
( wait $EGGNOG_PID && \
    nohup uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strains <NEW_STRAIN> > logs/tcdb/<NEW_STRAIN>.log 2>&1 ) &

# Add new tools here as they're created via /add-a-tool — note the wave they belong in
```

If you launch tcdb-diamond before eggnog has finished, the runner still
completes but every call's `egn_agreement = "extends"` (it can't see the
eggNOG TC assignments) — silently degraded output that's hard to spot
after-the-fact. The subshell + `wait $EGGNOG_PID` gate is what prevents this
without stalling the foreground workflow.

### Multi-strain variant

For a batch (e.g. all 8 missing Prochlorococcus strains), **omit the
`--strain` flag entirely**. Each runner then iterates `cyanobacteria_genomes.csv`
and skips strains that already have output — new strains run, old strains
no-op. The Wave 2 gate still works because eggnog only exits once its full
CSV walk finishes:

```bash
mkdir -p logs/{eggnog,psortb,tcdb,signalp}

# Wave 1 — same three tools, no --strain flag → iterate the whole CSV
nohup uv run python .claude/skills/eggnog-run/run_eggnog.py > logs/eggnog/batch.log 2>&1 &
EGGNOG_PID=$!
nohup uv run python .claude/skills/psortb-run/run_psortb.py > logs/psortb/batch.log 2>&1 &
nohup uv run python .claude/skills/signalp-run/run_signalp.py > logs/signalp/batch.log 2>&1 &

# Wave 2 — tcdb-diamond, after eggnog has annotated every strain in the CSV
( wait $EGGNOG_PID && \
    nohup uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py > logs/tcdb/batch.log 2>&1 ) &
```

Pass `--force` to a specific runner only if you want to re-process strains
that already have output for that tool.

> **Note on deferred conventions** — two divergences between this snippet
> and `/add-a-tool`'s canonical patterns, both pending opportunistic
> migration on the existing runners:
>
> - **`--strain` vs `--strains`** — psortb / eggnog / signalp accept
>   `--strain` (singular); tcdb-diamond accepts `--strains` (plural,
>   `nargs='+'`). `/add-a-tool`'s canonical is `--strains` plural.
> - **Log location** — the canonical convention is `logs/<tool>/<strain>.log`
>   (subfolder per tool); existing runners write the redirect target
>   shown above. The `mkdir -p` ahead of the runs creates the subfolders
>   so the redirect lands cleanly; once the runners themselves migrate
>   they'll create the subfolder internally.
>
> New `/add-a-tool` additions should follow the canonical patterns from
> day one.

Per-strain wallclock ranges from ~10 min (eggNOG with cached DB) to ~30 min
(PSORTb on a 5K-protein heterotroph). **The `nohup` block above (and its
multi-strain variant) is the canonical list** of per-strain tools — when
[`/add-a-tool`](../add-a-tool/SKILL.md) ships a new runner, its Step 6
appends a line to both blocks so future strains pick it up automatically.

eggNOG output is consumed by `prepare_data.sh --steps 1 2` (step 2 picks up
`<strain>.emapper.annotations` and merges it into
`gene_annotations_merged.json`). PSORTb / tcdb-diamond / signalp are Phase-1
only — their `<data_dir>/<tool>/<strain>.calls.json` artifacts sit in the
strain cache until each tool's Phase-2 integration spec lands.

## Step 4 — Snapshot the KG

```bash
STRAIN=<your-strain-name>   # REQUIRED — single-strain: the strain you're deploying (e.g., MIT9312, SS120, BL107).
                            # Multi-strain batch: use a batch tag (e.g., missing_pro_strains, batch_2026_05_20).
                            # The rest of the steps reference ${STRAIN}; set it once here.
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save before_${STRAIN}
```

Neo4j must be running. The snapshot is what step 9 compares against to detect
silent edge loss.

## Step 5 — Spot-check resolved CSVs (pre-Docker)

Open the `<csv>_resolved.csv` files prepare_data step 4 just wrote and skim the
`resolution_method` column. Authoritative match-rate verification happens
post-Docker in step 8 — that pass cross-references the import report (which
only exists after `neo4j-admin import` runs).

For known per-strain ID quirks see
[references/known_id_formats.md](references/known_id_formats.md).

## Step 6 — Validate paperconfigs + add unit tests for any new fix pattern

If step 5 surfaced a paperconfig fix, validate + add a regression test:

```bash
uv run pytest tests/test_paperconfig_validation.py -v   # all paperconfigs must pass
# Add tests to tests/test_gene_id_graph.py covering:
# - the failure mode (what broke and why)
# - the fix (what the correct id_type / config is)
uv run pytest tests/test_gene_id_graph.py -q
```

The unit test should describe **why** the fix was needed — future strains
hitting the same pattern will read these tests as documentation.

## Step 7 — Rebuild the Docker KG

```bash
docker compose down deploy app  # release the Neo4j lock
docker compose up -d --build
```

Rebuild takes ~30 min. Proceed to step 8 in parallel.

## Step 8 — Verify gene-ID resolution + update KG-validity / snapshot tests

After Docker is up, run the canonical gene-ID check (now has the import
report available):

```bash
uv run python .claude/skills/check-gene-ids/check_gene_ids.py
```

Expect ≥ 80% match rate per paper. Lower → use the diagnostic playbook at
[references/id_resolution_diagnostics.md](references/id_resolution_diagnostics.md)
(Phase 1/2/3 anchor logic, footnote artifacts, MIT9313 multi-annotation trap,
Alteromonas dual-assembly traps, etc.). If you fix a paperconfig here, re-run
`prepare_data.sh --strains <NEW_STRAIN> --steps 3 4 5 6 7 --force`, then
rebuild Docker.

Then run the KG validity suite:

```bash
pytest -m kg -v
```

If `test_organism.py::test_all_expected_strains_present` fails because the new
strain isn't listed, add it to the strains list there. **Don't expand the
assertion silently** — the test is intentionally an explicit "did you mean to
add this?" gate.

Regenerate the regression snapshot once all expected paper edges have landed:

```bash
uv run python tests/kg_validity/generate_snapshot.py
git add tests/kg_validity/snapshot_data.json
```

The committed snapshot is what catches silent data loss on future rebuilds.
A strain without snapshot coverage is invisible to regression.

## Step 9 — Compare snapshot

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before_${STRAIN}
```

Detects per-paper regressions (lost edges) vs improvements (gained edges from
better ID resolution). Lost edges that aren't explained by intentional
paperconfig changes = investigate before merging.

## Step 10 — Loop back through prepare_data once tool outputs are in (brand-new only)

Some tools (eggNOG today; SignalP / InterProScan as Phase-2 integrations land)
are consumed by `prepare_data.sh --steps 1 2` when their output exists. Wait
for every background tool from step 3 to finish, then re-run steps 1–7:

```bash
bash scripts/prepare_data.sh --strains <NEW_STRAIN> --steps 1 2 3 4 5 6 7 --force
```

Then redo step 7 (Docker rebuild), step 8 (post-Docker check + KG validity),
step 9 (snapshot compare) so the KG picks up the enriched annotations.

Tools that are still Phase-1-only (PSORTb, tcdb-diamond as of 2026-05) don't
strictly need this loop — their outputs sit in the strain cache until their
Phase-2 integration spec lands. The loop becomes mandatory once a tool's Phase
2 wires into `build_gene_annotations.py` or an adapter.

## Step 11 — Final sanity sweep

```bash
pytest -m "not slow and not kg" -v    # full unit suite
pytest -m kg -v                       # KG validity (Docker must be up)
uv run python .claude/skills/check-gene-ids/check_gene_ids.py
```

Then update the memory / docs trail:

- Add a line to `memory/MEMORY.md` under recent strain additions.
- If the strain has a quirk worth remembering (dual assembly, RAST/IMG annotation, MarRef proxy, cross-assembly bridging), one-line entry in `memory/known_bugs.md` or under "Common Gene ID Mismatches" in MEMORY.
- If a follow-up is open (e.g., DEH24→MADE_RS diamond mapping deferred), record it in `plans/<strain>_deploy.md`.
- Append the per-strain ID format to [references/known_id_formats.md](references/known_id_formats.md).

## Worked examples to copy from

| Strain | Reference | Lesson |
|---|---|---|
| SS120 | `docs/superpowers/specs/2026-04-14-csv-ready-papers-batch-design.md` | Standard Pro strain with Cyanorak; gene-ID resolution natively at 100% via `Pro_NNNN` locus tags |
| BL107 | same | Synechococcus draft assembly (WGS scaffolds); cyanorak `RefSeq` column blank — verify against NCBI Datasets, not Cyanorak |
| HP15 (Marinobacter) | same | Heterotroph reference-proteome-match (no Cyanorak); reference DB MarRef v6 |
| Alt_MarRef (*A. mediterranea* DE) | same | Non-NCBI `DEH24_*` locus prefix → required diamond-based id_translation bridge |
| MIT1002 | `plans/mit1002_deploy.md` | Dual-assembly trap; pick the assembly the paper actually used; document the choice |
| EZ55 | `plans/ez55_deploy.md` | Cross-assembly protein bridging via `scripts/map_img_to_ncbi_proteins.py` |

## Reference material (loaded on demand)

- [`references/id_resolution_diagnostics.md`](references/id_resolution_diagnostics.md) — playbook for when `/check-gene-ids` returns < 80%. Covers Phase 1/2/3 anchor logic, footnote artifacts, multi-annotation traps (MIT9313), Alteromonas dual-assembly diagnostics.
- [`references/cross_assembly_bridging.md`](references/cross_assembly_bridging.md) — diamond-based protein sequence matching for papers whose gene IDs come from a different assembly (draft / RAST / IMG). The EZ55 + MIT1002 workflow.
- [`references/known_id_formats.md`](references/known_id_formats.md) — per-strain table of which ID formats are seen in which paper, where they come from, and how they resolve. Append-only history.

## Open follow-ups

- **`genomes_registry.py` refactor** — collapses sites 2 and 3 above into "just edit the CSV". Designed but deferred (see `plans/gene_id_mapping_v2_status.md`). Until landed, the three-site dance is mandatory.
- **MarRef-style community proteomics** — `reference_proteome_match` organisms still need extra paperconfig conventions; see `docs/community_proteomics_marref_saga.md`.
- **Undeployed strains referenced by paperconfigs**: Pseudohoeflea + Thalassospira (ziegler 2025) and the strains tracked in `docs/superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md`.

## Workflow summary

Brand-new strain: 0 → 1 → 2 → 3 (kick off in background) → 4 → 5 → 6 → 7 → 8 → 9 → 10 (when tools done) → 11.

Re-deploy of existing strain: edit paperconfig → `bash scripts/prepare_data.sh --strains <S> --steps 3 4 5 6 7 --force` → 4 → 5 → 6 → 7 → 8 → 9 → 11.
