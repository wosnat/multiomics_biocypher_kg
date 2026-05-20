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

## Quick map

| Stage | File / command | Brand-new only? |
|---|---|---|
| Identify accession | NCBI Datasets + Cyanorak organism tables | yes |
| Register site 1 | `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` | yes |
| Register site 2 | `multiomics_kg/utils/gene_id_utils.py` (`ORGANISM_TO_GENOME_DIR`) | yes |
| Register site 3 | `scripts/validate_paperconfig.py` (`CANONICAL_GENOMIC_ORGANISMS`) | yes |
| Initial download | `bash scripts/prepare_data.sh --strains <S> --steps 0 1 2 3 4` | yes |
| Per-strain tools (parallel) | `eggnog-run`, `psortb-run`, `tcdb-diamond`, `signalp-run`, … | yes (rerun on paperconfig change only if id_translation changes) |
| Snapshot KG | `omics-edge-snapshot --save before_<S>` | both |
| Rebuild `gene_id_mapping.json` | `python -m multiomics_kg.download.build_gene_id_mapping --strains <S> --force` | both |
| Re-resolve paper CSVs | `python -m multiomics_kg.download.resolve_paper_ids --force` | both |
| Verify match rates | `/check-gene-ids` | both |
| Docker rebuild | `docker compose up -d --build` | both |
| Compare snapshot | `omics-edge-snapshot --compare before_<S>` | both |
| KG validity + tests | `pytest -m kg`, regenerate `snapshot_data.json` | both |
| Loop back | re-run `prepare_data.sh --steps 1 2 3 4` once tool outputs land | brand-new |

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

## Step 2 — Download + per-strain annotation tables (brand-new only)

```bash
bash scripts/prepare_data.sh --strains <NEW_STRAIN> --steps 0 1 2 3 4
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
| 4 | `<paper>/<csv>_resolved.csv` for any paper already using the strain |

Logs land in `logs/prepare_data_step{N}.log`. Cyanorak is the flakiest server in
the chain — if step 0.2 hangs, `--skip-cyanorak` and retry later.

## Step 3 — Per-strain tools (background, parallel) (brand-new only)

Every integrated tool with a Phase-1 SKILL must be run for the new strain so
its output is available before Phase-2 integrations re-merge into
`gene_annotations_merged.json`. Run these as background jobs — they don't
depend on each other.

```bash
# eggNOG — functional annotation (read by prepare_data step 2 when present)
nohup uv run python .claude/skills/eggnog-run/run_eggnog.py --strain <NEW_STRAIN> > logs/eggnog_<NEW_STRAIN>.log 2>&1 &

# PSORTb — Gram-negative subcellular localization
nohup uv run python .claude/skills/psortb-run/run_psortb.py --strain <NEW_STRAIN> > logs/psortb_<NEW_STRAIN>.log 2>&1 &

# TCDB-diamond — transporter classification
nohup uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strains <NEW_STRAIN> > logs/tcdb_<NEW_STRAIN>.log 2>&1 &

# SignalP — signal peptide prediction
nohup uv run python .claude/skills/signalp-run/run_signalp.py --strain <NEW_STRAIN> > logs/signalp_<NEW_STRAIN>.log 2>&1 &

# Add new tools here as they're created via /add-a-tool
```

> **Note on `--strain` vs `--strains`:** today's runners diverge — psortb /
> eggnog / signalp accept `--strain` (singular); tcdb-diamond accepts
> `--strains` (plural, `nargs='+'`). `/add-a-tool`'s canonical CLI standardizes
> on `--strains` plural going forward. The three singular runners will be
> migrated opportunistically to also accept `--strains`; until then, this
> block reflects each runner's current actual flag. New `/add-a-tool`
> additions should use `--strains` from day one.

Per-strain wallclock ranges from ~10 min (eggNOG with cached DB) to ~30 min
(PSORTb on a 5K-protein heterotroph). For the full list of tools that should
run on every new strain, see [`/add-a-tool`](../add-a-tool/SKILL.md) — adding a
tool to the list there is what makes future strains pick it up.

eggNOG output is consumed by `prepare_data.sh --steps 1 2` (step 2 picks up
`<strain>.emapper.annotations` and merges it into
`gene_annotations_merged.json`). PSORTb / tcdb-diamond / signalp are Phase-1
only — their `<data_dir>/<tool>/<strain>.calls.json` artifacts sit in the
strain cache until each tool's Phase-2 integration spec lands.

## Step 4 — Snapshot the KG

```bash
STRAIN=<your-strain-name>   # REQUIRED — replace with the strain you're deploying (e.g., MIT9312, SS120, BL107).
                            # The rest of the steps reference ${STRAIN}; set it once here.
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save before_${STRAIN}
```

Neo4j must be running. The snapshot is what step 11 compares against to detect
silent edge loss.

## Step 5 — Rebuild `gene_id_mapping.json`

```bash
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains $STRAIN --force
```

Review the diagnostic report:

```
cache/data/<Organism>/genomes/<Strain>/gene_id_mapping_report.json
```

What this build auto-includes beyond `gene_annotations_merged.json` (no
paperconfig needed):

| Source | File | What it adds |
|---|---|---|
| GCF `cds_from_genomic.fna` | `cache/.../genomes/<Strain>/genomic_gca.gff` | `cds_fna_id` Tier 1 — full `lcl|<chr>_cds_<protein>_<n>` IDs; also `locus_tag_ncbi`, `old_locus_tag`, `protein_id_refseq`, `gene_name` from FNA headers |
| GCA GFF | `cache/.../genomes/<Strain>/genomic_gca.gff` | Old locus_tag format (e.g. `PMT9312_1938`), old protein IDs (e.g. `ABS83190.1`), Cyanorak locus tag cross-refs, `Alternative locus ID` from Note field (ProPortal-era IDs like `P9313_NNNNN`, `PMED4_NNNNN`) |

Both files are downloaded by `scripts/prepare_data.sh` step 0 NCBI sub-steps —
no manual action needed.

## Step 6 — Re-resolve paper CSVs

```bash
uv run python -m multiomics_kg.download.resolve_paper_ids --force
```

This walks every `csv` supplementary table in every paperconfig, loads the
source CSV, resolves each `name_col` value via the v2 mapper, and writes
`<stem>_resolved.csv` alongside the original with two extra columns
(`locus_tag`, `resolution_method`). The omics adapter consumes the resolved
CSV, not the raw one.

### Adding non-standard ID sources (paperconfig)

When a paper uses IDs not in the NCBI/Cyanorak annotations, add entries to the
paper's `paperconfig.yaml` **before** running step 5.

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

For paper IDs that come from a different assembly entirely
(draft genome, IMG/RAST/custom annotation), use **cross-assembly protein
sequence bridging** — see [references/cross_assembly_bridging.md](references/cross_assembly_bridging.md).

## Step 7 — Verify gene-ID resolution

```bash
uv run python .claude/skills/check-gene-ids/check_gene_ids.py
```

Expect ≥ 80% match rate per paper. Lower → use the diagnostic playbook at
[references/id_resolution_diagnostics.md](references/id_resolution_diagnostics.md)
(Phase 1/2/3 anchor logic, footnote artifacts, MIT9313 multi-annotation trap,
Alteromonas dual-assembly traps, etc.).

For known per-strain ID quirks see
[references/known_id_formats.md](references/known_id_formats.md).

### Removing `_with_locus_tag.csv` workarounds

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

## Step 8 — Validate paperconfigs + add unit tests for any new fix pattern

If step 7 surfaced a paperconfig fix, validate + add a regression test:

```bash
uv run pytest tests/test_paperconfig_validation.py -v   # all paperconfigs must pass
# Add tests to tests/test_gene_id_graph.py covering:
# - the failure mode (what broke and why)
# - the fix (what the correct id_type / config is)
uv run pytest tests/test_gene_id_graph.py -q
```

The unit test should describe **why** the fix was needed — future strains
hitting the same pattern will read these tests as documentation.

## Step 9 — Rebuild the Docker KG

```bash
docker compose down deploy app  # release the Neo4j lock
docker compose up -d --build
```

Rebuild takes ~30 min. Proceed to step 10 in parallel.

## Step 10 — Update KG-validity + snapshot tests

After Docker is up:

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

## Step 11 — Compare snapshot

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before_${STRAIN}
```

Detects per-paper regressions (lost edges) vs improvements (gained edges from
better ID resolution). Lost edges that aren't explained by intentional
paperconfig changes = investigate before merging.

## Step 12 — Loop back through prepare_data once tool outputs are in (brand-new only)

Some tools (eggNOG today; SignalP / InterProScan as Phase-2 integrations land)
are consumed by `prepare_data.sh --steps 1 2` when their output exists. Once
every background tool from step 3 has finished, re-run steps 1–4:

```bash
bash scripts/prepare_data.sh --strains <NEW_STRAIN> --steps 1 2 3 4 --force
```

Then redo step 6 (re-resolve), step 7 (check), step 9 (Docker rebuild), step 11
(snapshot compare) so the KG picks up the enriched annotations.

Tools that are still Phase-1-only (PSORTb, tcdb-diamond as of 2026-05) don't
strictly need this loop — their outputs sit in the strain cache until their
Phase-2 integration spec lands. The loop becomes mandatory once a tool's Phase
2 wires into `build_gene_annotations.py` or an adapter.

## Step 13 — Final sanity sweep

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

Brand-new strain: 0 → 1 → 2 → 3 (kick off) → 4 → 5 → 6 → 7 → 8 → 9 → 10 → 11 → 12 (when tools done) → 13.

Re-deploy of existing strain: 4 → 5 → 6 → 7 → 8 → 9 → 10 → 11 → 13.
