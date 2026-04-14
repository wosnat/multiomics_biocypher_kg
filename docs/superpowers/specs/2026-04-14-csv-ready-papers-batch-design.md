---
date: 2026-04-14
status: design
topic: Add Domínguez 2017, Fuszard 2012, Moreno 2023 paperconfigs (CSV-ready batch)
---

# CSV-Ready Papers Batch — Domínguez 2017, Fuszard 2012, Moreno 2023

## Background

Eight new Prochlorococcus papers had supplementary tables committed to the repo this week with no paperconfig yet. After triage, the following batch of three is ready to integrate because their per-strain CSVs are already extracted and the data fits the existing `Changes_expression_of` DE schema:

- **Domínguez 2017** — quantitative proteomics, *Prochlorococcus marinus* SS120, control vs azaserine (N-stress proxy)
- **Fuszard 2012** — iTRAQ proteomics, MIT9312 + NATL2A + SS120, Pi replete vs Pi deplete
- **Moreno 2023** — proteomics + metabolomics, MED4 + SS120 + WH7803 + WH8102 + BL107, glucose × light/dark in non-axenic cultures

The remaining five are deferred:
- **Szul 2019, McDonagh 2012, Pandhal 2007** — only docx/xls supp; need table extraction
- **Biller 2022** — extracellular vesicle composition, not DE
- **Kujawinski 2023** — metabolite diversity catalog, not DE

## Goals

1. Add four new genome strains (SS120, BL107, *Marinobacter adhaerens* HP15, *Alteromonas mediterranea* DE) to the pipeline as full data-producing strains.
2. Create three new paperconfigs that validate against `validate_paperconfig.py` and resolve gene IDs at ≥ 80% match rate.
3. Generate ~64 new `Experiment` nodes (1 from Domínguez + 3 from Fuszard + ~60 from Moreno) and corresponding `Changes_expression_of` edges.

## Out of Scope

- Metabolomics extracted from Moreno S1 (`spectrum.03275-22-s0001.xlsx`) — proteomics only in this batch.
- The five deferred papers above.
- Schema changes (the `fold_change_type` field already exists from the Synechococcus integration; `adjusted_p_value` already accepts null).

## Hard Prerequisites — New Strains

All four strains must be added to `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` and pushed through `scripts/prepare_data.sh --strains <strain> --steps 0 1 2` before any of the three paperconfigs can resolve gene IDs.

| Strain | Cyanorak name | Clade | Assembly accession | NCBI taxid | Locus prefix | Notes |
|---|---|---|---|---|---|---|
| **SS120** | `Pro_SS120` | LLII | `GCF_000007925.1` (complete; chromosome `NC_005042` / `AE017126`) | 167539 | `Pro_NNNN` | Used by Domínguez, Fuszard, Moreno |
| **BL107** | `Syn_BL107` | sub 5.1 / clade IV | `GCF_000153805.1` (WGS draft scaffolds; cyanorak RefSeq column blank) | 313625 | `BL107_NNNNN` | Used by Moreno only |
| **Marinobacter adhaerens HP15** | n/a (no cyanorak entry) | — | `GCF_000166295.1` (complete; chromosome `NC_017506` + 2 plasmids) | 225937 | `HP15_RS#####` | Used by Moreno S4 |
| **Alteromonas mediterranea DE** | n/a (no cyanorak entry) | — | `GCF_000020585.3` (complete; chromosome `NC_011138.3`) | 1774373 (strain-level; species `314275`; genus `232` is what the paper CSV reports) | RefSeq uses `MADE_RS#####` — paper data uses non-NCBI `DEH24_NNNNN` (see ID-bridge note below) | Used by Moreno S3 |

For Marinobacter and Alteromonas, the cyanorak `Pigment type` / `SubCluster` columns do not apply — leave blank as for existing heterotrophs (W3-18-1, KT2440, DSS-3, MruberA).

### Verification trail

- **NCBI**: assemblies and taxids confirmed via NCBI Datasets (April 2026). Earlier draft of this spec listed `GCF_000153825.1` for BL107 (actually *Synechococcus* RS9916) and `GCF_001886215.1` for *A. mediterranea* DE (does not exist) — both corrected above.
- **Cyanorak**: SS120 + BL107 metadata (clade, locus prefix) cross-matched against `data/Cyanorak  Organism Table  prochlorococcus.csv` and `synechococcus.csv`; cyanorak gives nucleotide accessions (`AE017126` / `NC_005042` for SS120; `AATZ00000000` WGS base for BL107) which point to the NCBI assemblies above.
- **Paper methods**:
  - **Domínguez 2017** explicitly cites SS120 = CCMP1375 (Roscoff Culture Collection) and used UniProt SS120 entries (1,881 entries, 18 May 2015) — taxid 167539 confirmed.
  - **Fuszard 2012** Materials & Methods confirms strains MIT9312, NATL2A, SS120 grown in Pro99; n=3 iTRAQ replicates; significance cutoffs ±1.6 / 0.6 fold-change (no per-protein p-values reported).
  - **Moreno 2023** Materials & Methods confirms 5 strains used for proteomics (SS120, MED4, MIT9313, WH8102, WH7803, BL107) but **MIT9313 was excluded from results** because cultures stopped growing during the proteomics campaign — supplementary tables therefore cover 5 strains. Reference DB for protein matching is **MarRef v6** (per Fig 1 caption) — this explains the `DEH24_*` Alteromonas locus prefix (MarRef custom annotation, not NCBI), and motivates the id-translation bridge described under "Open Items".

## Per-Paper Plans

### Domínguez 2017

- **Files**: `paper.pdf`, `table s3 Upregulated prots rel quant`, `table s3 Downregulated prots rel quant`, `table s4 Carbohydrate metabolism prots`
- **Strain**: SS120
- **Experiment** (1):
  - `azaserine_ss120_proteomics`: control vs azaserine (10 µM); LC-MS/MS Progenesis label-free; `treatment_type: ["nitrogen"]` (azaserine inhibits glutamine synthase → effective N-stress); `omics_type: PROTEOMICS`
- **Statistical analyses** (1): merge `table s3 Upregulated prots rel quant sys003172107st8.csv` + `table s3 Downregulated prots rel quant sys003172107st8.csv` into a new file `table s3 Combined_modified.csv` (leave the two originals untouched — see "CSV transformation convention" below) with a signed `log2_fold_change` derived from the `Max fold change` column and the sign from `Highest mean condition` / `Lowest mean condition` (negative when `Highest mean condition == Control`); `adjusted_p_value` from `Anova (p)`; `pvalue_threshold: 0.05` (matches the legend's "P value < 0.05" filter).
- **Skip**: table s4 is a hand-curated subset of S3 — no new information.
- **ID resolution**: `name_col: Description`, parse `GN=Pro_NNNN` from the description string. Add an `id_translation` entry built from the `Description` column (`id_type: locus_tag` for the Pro_NNNN extraction). Most rows also have UniProt entry name (`Q7VDF9_PROMA`) and accession (in row index `1::Q7VDF9_PROMA`).

### Fuszard 2012

- **Files**: 3 strain-specific CSVs (`table s1 MIT9312`, `table s1 NATL2A`, `table s1 SS120`)
- **Strains**: MIT9312 (existing), NATL2A (existing), SS120 (new)
- **Experiments** (3, one per strain):
  - `pi_limitation_<strain>_itraq`: Pi replete (50 µM NaH₂PO₄) vs Pi deplete (10 µM NaH₂PO₄); iTRAQ; `treatment_type: ["phosphorus"]`; `omics_type: PROTEOMICS`
  - `growth_phase`: `["stationary"]` for MIT9312 + SS120; `["exponential"]` for NATL2A (paper notes NATL2A was harvested mid-exponential because stock cultures collapsed at OD > 0.4)
  - `medium`: `Pro99`; `temperature: 21°C` growth, `23°C` experimental; `light_condition: 13:11 L:D`; `light_intensity: 30/20/10 µE m⁻² s⁻¹` per strain
- **Statistical analyses** (3): one per CSV. Compute `log2_fold_change = log2(Ratio of means)`. `adjusted_p_value`: leave **null** (paper does not report per-protein p-values; only fold-change cutoffs of 1.6 / 0.6 from n=3 iTRAQ replicates).
  - `fold_change_type`: use the existing fold-change-threshold significance mode (added during Synechococcus integration); thresholds `up: 1.6`, `down: 0.6` (or `±log2(1.6) ≈ ±0.678`, depending on the field's existing semantics — verify in `omics_adapter.py`).
- **ID resolution**: `name_col: ORF` (column "ORF") — `Pro_NNNN` for SS120, `PMT9312_NNNN` for MIT9312 (matches existing `old_locus_tag` mapping), `PMN2A_NNNN` for NATL2A. No id_translation entry expected.

### Moreno 2023

The most complex of the three: 5 strains, 3 organism perspectives (cyano + Alt + Marino), 4 comparisons each ⇒ up to **~60 experiments**.

- **Files** in scope: 15 CSVs (5 strains × 3 perspectives); skip Moreno S1 metabolomics (`spectrum.03275-22-s0001.xlsx`).
- **Strains as primary organism**:
  - Cyano S2: MED4, SS120, WH7803, WH8102, BL107
  - Alt S3: *Alteromonas mediterranea* DE (5 cocultures, one per cyano partner)
  - Marino S4: *Marinobacter adhaerens* HP15 (5 cocultures, one per cyano partner)
- **Experiment naming convention** (avoids ambiguity between "high light" and "high glucose" by ordering as `<light_condition>_<glucose_dose>_glucose`):
  - Cyano: `<comparison>_<strain>_proteomics` — e.g. `light_low_glucose_med4_proteomics`, `light_high_glucose_med4_proteomics`, `dark_low_glucose_med4_proteomics`, `dark_high_glucose_med4_proteomics`
  - Heterotroph: `<comparison>_<heterotroph_short>_in_<cyano_partner>_proteomics` — e.g. `light_high_glucose_alt_in_med4_proteomics`, `dark_low_glucose_marino_in_ss120_proteomics`
- **Comparison set** (4 per CSV, matching the columns observed):
  - `light_low_glucose`: light vs light + 100 nM glucose
  - `light_high_glucose`: light vs light + 5 mM glucose
  - `dark_low_glucose`: dark vs dark + 100 nM glucose
  - `dark_high_glucose`: dark vs dark + 5 mM glucose
  - (Skip the L vs D overall comparison — only present in compiled S1, not per-strain CSVs.)
- **Common experiment metadata**:
  - `omics_type: PROTEOMICS`
  - `treatment_type: ["carbon"]` for the four glucose comparisons
  - `background_factors`:
    - cyano S2: `["coculture"]` + `["light"]` or `["darkness"]` per condition (heterotrophs are uncontrolled community contamination — no `treatment_organism`)
    - Alt S3 + Marino S4: `["coculture"]` + light/darkness; `treatment_organism` = the cyano partner (organism node already in KG); `treatment_taxid` = cyano taxid
  - `medium`: cyano media (Pro99 / SN), specifics per strain
  - `light_condition`: `"continuous light"` or `"continuous darkness"` per condition
- **Statistical analyses**: 4 per CSV (one per `TTEST.X.Y` / `FC.X.Y` column pair).
  - `logfc_col`: `FC.L.LGnM` etc. — but these are linear ratios, not log2; compute `log2(FC)` upstream during ID-resolution step or note in the adapter that this paper's `logfc_col` is linear (Synechococcus integration added `fold_change_type` for cases like this — apply the same).
  - `adjusted_p_value_col`: `TTEST.L.LGnM` etc.
  - `pvalue_threshold: 0.05`
- **ID resolution**:
  - Cyano S2: `name_col: Accession` (UniProt accession), with optional `id_columns: [{column: KEGG, id_type: locus_tag}]` since the `KEGG` column carries `Pro_NNNN` / `PMM_NNNN` / `BL107_NNNNN` etc. as a fallback.
  - Alt S3 + Marino S4: parse `GN=` from `Description` column to extract locus tag (`DEH24_NNNNN` / `HP15_NNNN`); fall back to the gene-name when GN is just a 3-letter symbol (`typA`, `gltA`, etc.). Will require an `id_translation` entry per heterotroph.
  - **Alteromonas DEH24 → MADE_RS bridge**: NCBI RefSeq for *A. mediterranea* DE uses locus prefix `MADE_RS#####`, but the Moreno data uses `DEH24_NNNNN` (a non-NCBI prefix from a re-annotation, likely IMG/PATRIC/MarRef). An `id_translation` entry is needed to map `DEH24_*` → `MADE_RS*`. Same pattern as the EZ55 (`AEZ55_*` → NCBI locus tags) and MIT1002 (`fig|*` → `MIT1002_*`) deployments — likely a `generate: diamond_protein_match` block against a DEH24-annotated protein FASTA. Source FASTA for DEH24_ may need to be sourced from IMG (Genome ID for *A. mediterranea* DE) or from the Moreno authors directly.

## CSV Transformation Convention

Whenever a paper's source CSV needs to be transformed before ingestion (merging up/down halves, computing `log2(FC)` from a linear ratio, splitting a multi-strain table, signing a fold-change column, etc.), **always write the result to a new file with a `_modified.csv` suffix and leave the originals untouched**. The paperconfig's `filename:` field then points at the `_modified.csv`. This preserves provenance — the original supp file matches what the paper deposited, and the transformation is reproducible from a script.

Concrete cases in this batch:
- Domínguez 2017: combine S3 up + S3 down → `table s3 Combined_modified.csv`
- Fuszard 2012: add `log2_fold_change` column from `Ratio of means` → `table s1 <strain> Quantitation data_modified.csv`
- Moreno 2023: each `table sN <strain> <whatever>.csv` may need a `_modified.csv` if `log2(FC)` columns must be precomputed before ingestion.

Each transformation should be done with a small reusable script committed under `scripts/` (or notebook cell), so the `_modified.csv` is regeneratable.

## Schema & Decision Rationale

### Why `fold_change_type` for Fuszard

Fuszard's Materials and Methods (Supplementary Material 1) state n=3 biological triplicates per condition (iTRAQ tags 113-115 replete, 116-118 deplete), but the publication-derived CSV has only the per-condition mean and SD (computed across replicates), not per-protein p-values. The authors' own significance criterion is fold-change-only (`> 1.6` or `< 0.6`). Computing a t-test from mean+SD without raw replicate values would be a hack. Leaving `adjusted_p_value` null and using the existing `fold_change_type` mechanism cleanly mirrors the authors' methodology.

### Why Moreno S2 has no `treatment_organism`

Moreno explicitly designed the experiment as **non-axenic** ("we utilized nonaxenic cultures as a proxy"), with naturally-occurring heterotrophs (Alteromonas, Marinobacter, Labrenzia, Salinimonas, …) varying *uncontrolled* across cultures. From the cyano perspective, there is no defined coculture partner. Setting `treatment_organism` would imply a controlled binary culture, which would be misleading. The `["coculture"]` background factor flags the non-axenic condition.

### Why Moreno S3/S4 *do* have `treatment_organism`

Once we shift perspective to Alt/Marino as the primary organism, the cyano *is* a defined partner from the heterotroph's standpoint (the heterotroph populations are characterized in cocultures with each named cyano).

### Why merge Domínguez S3 up/down

The two "Upregulated" / "Downregulated" CSVs are halves of a single Anova-tested comparison split by direction. Combining them into one CSV with a signed `log2_fold_change` aligns with the schema's "one statistical analysis = one comparison" convention and avoids the need for two redundant Experiment nodes.

## Acceptance Criteria

1. **Strains added**:
   - 4 new rows in `cyanobacteria_genomes.csv`
   - `cache/data/<org>/genomes/<strain>/` populated for all 4 (NCBI GFF, GBFF, protein FASTA; UniProt for each unique taxid; gene_mapping.csv; gene_annotations_merged.json)
   - eggNOG annotations available (or a documented decision to defer eggNOG for the new strains until the next batch)
2. **Paperconfigs validate**: all three pass `python .claude/skills/paperconfig/validate_paperconfig.py <path>`
3. **Gene-ID match rates** (via `/check-gene-ids`): each statistical analysis ≥ 80% of paper rows resolved to locus tags. Lower rates allowed only for heterotroph S3/S4 with documented reason.
4. **KG build smoke test**:
   - 3 new `Publication` nodes
   - ~64 new `Experiment` nodes (1 + 3 + ~60), within ±5 of expected
   - All new `Changes_expression_of` edges have valid (non-dangling) `Gene` targets per `tests/kg_validity/test_import_report.py`
5. **KG validity tests** (`pytest -m kg`) still pass; `/omics-edge-snapshot` shows only additions, no regressions on existing papers.
6. **Documentation**:
   - Update `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` with 3 new entries
   - Update `CLAUDE.md` strain list and node-count summary
   - Update `MEMORY.md` to reflect the new strains + papers

## Execution Order

Single PR, four sequential phases (revised — paperconfigs land before `prepare_data.sh` so any obvious config errors are caught early without re-running downloads):

1. **Add all 4 strains**: append 4 rows to `cyanobacteria_genomes.csv`. No data downloads yet.
2. **Add 3 paperconfigs + write any `*_modified.csv` files; validate**: write paperconfigs for Domínguez 2017, Fuszard 2012, Moreno 2023; run `validate_paperconfig.py` on each; sanity-check `_modified.csv` outputs by spot-reading a few rows. Validation should pass before download starts.
3. **Run `prepare_data.sh --strains SS120 BL107 HP15 AltMedDE --steps 0 1 2 3 4`**: NCBI + cyanorak (where applicable) downloads, gene_annotations builds, gene_id_mapping builds, paper-CSV resolution. Use `--skip-cyanorak` only as a fallback if the cyanorak server throttles.
4. **Build the KG**: `uv run python create_knowledge_graph.py`; verify CSV outputs land under `biocypher-log/example_knowledge_graph/`.

### Out of scope for this PR

- Manual gene-ID fixes beyond what `prepare_data.sh --steps 3 4` resolves automatically (no hand-curated `id_translation` entries; no DEH24 bridge work). Generate a **gene-ID status report** at the end summarizing per-paper match rates so the gaps are visible — fixes happen in a follow-up PR.
- The DEH24 → MADE_RS bridge for Moreno S3 (Alteromonas perspective). Skip those 5 CSVs in the Moreno paperconfig for this batch and note the deferral.
- Docker rebuild, KG validity tests against live Neo4j (`pytest -m kg`), and `/omics-edge-snapshot` comparison — done by the user post-merge once they verify the CSV outputs look reasonable.

### Nice-to-have (parallel, non-gating)

- eggNOG-mapper for HP15 + AltMedDE in the background after `prepare_data.sh --steps 0` finishes their NCBI protein FASTA. Functional annotation enrichment, not required for DE edges.

## Open Items (resolve during implementation)

- **`fold_change_type` semantics**: confirm whether the existing field expects raw fold-change thresholds (1.6 / 0.6) or log2-space thresholds (±0.678). Read `omics_adapter.py` and existing Synechococcus paperconfigs that use this field before writing the Fuszard config.
- **Moreno cyano `KEGG` column**: confirm it carries the locus tag in all 5 strain CSVs (verified for MED4 = `PMM1436` and SS120 = `Pro_1591`; need to spot-check WH7803, WH8102, BL107).

### Resolved during brainstorming

- ~~NCBI accessions/taxids for SS120, BL107, HP15, AltMedDE~~ — verified against NCBI Datasets, cyanorak tables, and paper methods sections (see "Verification trail").
- ~~Moreno S3 Alteromonas DEH24 bridge~~ — deferred (out of scope for this PR; see "Out of scope" in Execution Order).
- ~~Heterotroph eggNOG~~ — runs as a non-gating background "nice-to-have" parallel agent.
