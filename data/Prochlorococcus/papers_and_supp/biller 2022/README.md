# Biller 2022

**Citation:** Biller SJ, Lundeen RA, Hmelo LR, Becker KW, Arellano AA, Dooley K, Heal KR, Carlson LT, Van Mooy BAS, Ingalls AE, Chisholm SW. *Prochlorococcus* extracellular vesicles: molecular composition and adsorption to diverse microbes. Environmental Microbiology, 2022, 24(1): 420–435.
**DOI:** 10.1111/1462-2920.15834
**Organism(s):** *Prochlorococcus* MIT9312 (HL-adapted) and MIT9313 (LL-adapted)
**Topic:** Characterizes the lipid, protein, pigment, and metabolite content of extracellular vesicles (EVs) released by two *Prochlorococcus* ecotypes, and compares vesicle content to parent cells. Label-free shotgun proteomics identifies hundreds of proteins packaged in vesicles vs. cells (Tables S2-S4). Table S4 reports relative enrichment of proteins in vesicles vs. cells -- candidate per-gene quantitative compartment metric. Tables S5-S7 cover vesicle metabolites; Tables S8-S9 are method parameters.

## Classification

**Bucket A - integrated** (proteomics + Phase 2 metabolomics S7)

Proteomics S3 + S4 wired via 4 `derived_metrics_table` entries (per-strain biovolume-normalized vesicle/cell abundance + log2 vesicle/cell enrichment). Phase 2 metabolomics S7 wired as of 2026-05-06 via 4 `metabolite_assays_table` entries (cells + vesicles × MIT9312 + MIT9313, `cell_format: embedded_mean_sd_n` parsing the paper's `mean (sd), n=N` cells; 70/70 = 100% resolution after 10 aliases for paper-style spacing/punctuation variants). Both strains deployed.

**Still unintegrated** (lower-ROI follow-ups):
- **Table S1** (lipids/pigments/plastoquinones) — lipid-class names need LipidMaps integration; pigment names should resolve via current resolver. Defer until LipidMaps lands.
- **Table S5/S6** (untargeted metabolite features, ~4000 rows without compound IDs) — out of scope for v1; would need a `MassFeature` node type. See `docs/kg-changes/metabolomics-extension.md` § "Out of scope".
- **Table S2** raw spectral counts — sample-level, superseded by S3 biovolume-normalized abundances.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `Environmental Microbiology - 2021 - Biller - Prochlorococcus extracellular vesicles  molecular composition and adsorption.pdf` | PDF | Main article | reference | — |
| `emi15834-sup-0001-supinfo.pdf` | PDF | Appendix S1 — supplementary figures & methods | reference | — |
| `legend.txt` | TXT | Supplementary table legends | reference | — |
| `emi15834-sup-0002-tables1.xlsx` | XLSX | Table S1: complete lipid, pigment, and plastoquinone dataset for cells and vesicles of MIT9312 and MIT9313 | skip | Out of scope (lipidomics/metabolomics; not per-gene) |
| `emi15834-sup-0003-tables2.xlsx` | XLSX | Table S2: complete proteomics dataset — protein/peptide identification and raw spectral counts for MIT9312 and MIT9313 cells and vesicles | add | Raw spectral counts per sample → too sample-level; however per-protein detection flags (detected in vesicle? detected in cell?) could be distilled into a `derived_metrics_table` with boolean metrics |
| `emi15834-sup-0004-tables3.xlsx` | XLSX | Table S3: relative protein abundances (biovolume-normalized spectral counts) for MIT9312 and MIT9313 cells and vesicles | add | `derived_metrics_table` with `value_kind: numeric`, `rankable: true` — biovolume-normalized abundance per protein in each compartment (split into per-strain / per-compartment Experiments with `compartment: vesicle` vs `whole_cell`) |
| `emi15834-sup-0005-tables4.csv` / `emi15834-sup-0005-tables4.xlsx` | CSV / XLSX | Table S4: relative enrichment of proteins in vesicles vs. cells (fold-change + significance per gene, per strain) | add | **Primary candidate.** Add as `csv` with `statistical_analyses` — Experiment `compartment: vesicle`, `omics_type: PROTEOMICS`, treatment="vesicle fraction" vs control="whole cell"; two analyses (MIT9312, MIT9313). CSV is the preferred input. XLSX is a duplicate → skip the XLSX. |
| `emi15834-sup-0006-tables5.xlsx` | XLSX | Table S5: untargeted metabolite diversity in vesicle and cell samples | defer | MS feature list without compound IDs — out of scope until `MassFeature` node type lands |
| `emi15834-sup-0007-tables6.xlsx` | XLSX | Table S6: full biovolume-normalized peak areas for metabolites in vesicles and cells | defer | Same blocker as S5 — untargeted features without compound IDs |
| `emi15834-sup-0008-tables7.csv` / `emi15834-sup-0008-tables7.xlsx` | CSV / XLSX | Table S7: targeted metabolites detected in vesicles and cells (mean (sd), n=N out of 3 reps; analytical fraction per row) | already in | 4 `metabolite_assays_table` entries (cells + vesicles × MIT9312 + MIT9313); `cell_format: embedded_mean_sd_n`; CSV is the input, XLSX is duplicate. |
| `metabolite_aliases.yaml` | YAML | Per-paper free-text-name → primary-ID overrides (10 entries for paper-style spacing/punctuation) | reference | — |
| `emi15834-sup-0008-tables7_resolved.csv` / `emi15834-sup-0008-tables7_resolution_report.json` | CSV / JSON | Step 7 outputs (consumed by `metabolite_assay_adapter`); 70/70 = 100% resolved | — | auto-generated |
| `emi15834-sup-0009-tables8.xlsx` | XLSX | Table S8: experimental conditions for metabolite fractions | skip | Methods metadata |
| `emi15834-sup-0010-tables9.xlsx` | XLSX | Table S9: MS-DIAL parameters per analytical fraction | skip | Methods metadata |

## Current paperconfig summary

- Experiments defined: 6 — 2 PROTEOMICS (`vesicle_proteomics_mit9312`, `vesicle_proteomics_mit9313`) + 4 METABOLOMICS (cells + vesicles × MIT9312 + MIT9313). All `compartment: vesicle | whole_cell`, `treatment_type: [compartment]`.
- Statistical analyses (DE edges): 0 — non-DE quantitative proteomics + targeted metabolomics.
- Supplementary materials entry types: 4× `derived_metrics_table` (proteomics S3 + S4) + 4× `metabolite_assays_table` (metabolomics S7).
- Organisms covered: Prochlorococcus MIT9312, Prochlorococcus MIT9313.
- Non-DE evidence: 6 derived metrics (proteomics) + 4 metabolite assays (Phase 2 metabolomics).
- ID resolution (proteomics): MIT9312 uses `Gene Number` = `locus_tag`; MIT9313 uses `Gene Number` = `old_locus_tag`. Tables S4 use `Protein ID` (UniProt accession).
- Metabolite resolution (Phase 2): 70/70 = 100% via 60 `name_match` + 10 `alias_override` (paper-style spacing/parens variants for SAM, SAH, glycine betaine, sugar phosphates, glycerol-3P, Nε-acetyl-lysine, trans-hydroxyproline, retinal).

## Recommended actions

1. **No action** — proteomics (Tables S3, S4) and Phase 2 metabolomics (Table S7) are integrated.
2. **Defer** — Table S2 raw spectral counts are sample-level and superseded by S3 biovolume-normalized abundances; no separate boolean detection metric currently warranted.
3. **Defer (blocked on schema work)** — Tables S5, S6 are untargeted MS-feature lists without compound IDs; would need a `MassFeature` node type.
4. **Defer (blocked on LipidMaps integration)** — Table S1 lipids/pigments — pigments may resolve via current resolver but lipid-class names need LipidMaps.
5. **Skip** — Tables S8 and S9 are methods metadata.

## Notes

- Vesicle/cell compartment split is a natural use case for the `compartment` Experiment property (`vesicle` vs `whole_cell`) already supported by the non-DE evidence and Phase 2 metabolomics extensions. See `docs/kg-changes/non-de-evidence-extension.md` and `docs/kg-changes/metabolomics-extension.md`.
- Table S4 (relative enrichment = vesicle/cell fold-change with significance) maps cleanly to a differential-expression-style `csv` entry even though the "treatment" is a subcellular fraction rather than an environmental condition — DerivedMetric would also work but the log2FC + p-value structure better fits `statistical_analyses`.
- Table S7 was the first paper to exercise `cell_format: embedded_mean_sd_n` (parser was reserved in the spec but not wired to the row loop until 2026-05-06). Vesicle compartment is sparse: ~3-4 compounds detected (Carotene, trans Retinal, Vitamin K1, Tryptophan partial); cells compartment ~60 compounds per strain — that contrast is the paper's main signal.
- Both strains (MIT9312, MIT9313) are already in the KG — no organism work needed.
- MIT9312 uses PMT9312_#### locus tags; MIT9313 uses PMT#### — verify via `/check-gene-ids` before integration.
