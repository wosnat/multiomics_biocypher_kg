# Biller 2022

**Citation:** Biller SJ, Lundeen RA, Hmelo LR, Becker KW, Arellano AA, Dooley K, Heal KR, Carlson LT, Van Mooy BAS, Ingalls AE, Chisholm SW. *Prochlorococcus* extracellular vesicles: molecular composition and adsorption to diverse microbes. Environmental Microbiology, 2022, 24(1): 420–435.
**DOI:** 10.1111/1462-2920.15834
**Organism(s):** *Prochlorococcus* MIT9312 (HL-adapted) and MIT9313 (LL-adapted)
**Topic:** Characterizes the lipid, protein, pigment, and metabolite content of extracellular vesicles (EVs) released by two *Prochlorococcus* ecotypes, and compares vesicle content to parent cells. Label-free shotgun proteomics identifies hundreds of proteins packaged in vesicles vs. cells (Tables S2–S4). Table S4 reports relative enrichment of proteins in vesicles vs. cells — candidate per-gene quantitative compartment metric. Tables S5–S7 cover vesicle metabolites; Tables S8–S9 are method parameters.

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
| `emi15834-sup-0006-tables5.xlsx` | XLSX | Table S5: untargeted metabolite diversity in vesicle and cell samples | skip | Metabolomics — no gene targets |
| `emi15834-sup-0007-tables6.xlsx` | XLSX | Table S6: full biovolume-normalized peak areas for metabolites in vesicles and cells | skip | Metabolomics — no gene targets |
| `emi15834-sup-0008-tables7.xlsx` | XLSX | Table S7: targeted metabolites detected in vesicles and cells | skip | Metabolomics — no gene targets |
| `emi15834-sup-0009-tables8.xlsx` | XLSX | Table S8: experimental conditions for metabolite fractions | skip | Methods metadata |
| `emi15834-sup-0010-tables9.xlsx` | XLSX | Table S9: MS-DIAL parameters per analytical fraction | skip | Methods metadata |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Add** — Create `paperconfig.yaml` with two Experiments (one per strain), each with `compartment: vesicle`, `omics_type: PROTEOMICS`, `treatment_type: ["growth_phase"]` (or a new vocabulary value for vesicle-enrichment), and one `csv` supplementary_materials entry per strain pointing to `emi15834-sup-0005-tables4.csv`. Verify the CSV's gene ID column — likely PMT9312_#### / PMT9313_#### locus tags (resolvable directly against MIT9312/MIT9313 NCBI annotations).
2. **Add** — Add Table S3 as a `derived_metrics_table` with numeric `rankable: true` metrics: `biovolume_normalized_abundance_cell`, `biovolume_normalized_abundance_vesicle` per strain. Each yields `Derived_metric_quantifies_gene` edges.
3. **Add** — Add Table S2 as a `derived_metrics_table` with boolean flags (`detected_in_vesicles`, `detected_in_cells`) per strain using `Derived_metric_flags_gene` edges. Alternatively, fold these into the Table S3 numeric set.
4. **Skip** — Table S5–S9 and S1 (metabolomics / lipidomics / methods metadata). If a future metabolite-to-gene annotation is added, revisit Table S5/S7 for candidate `DerivedMetric`-style gene linkages via biosynthesis pathway annotations.
5. **Check gene IDs** — run `/check-gene-ids` on Table S4 CSV against MIT9312 and MIT9313; if the column is protein_id not locus_tag, add an `id_columns` entry with `id_type: protein_id_refseq`.

## Notes

- Vesicle/cell compartment split is a natural use case for the `compartment` Experiment property (`vesicle` vs `whole_cell`) already supported by the non-DE evidence extension. See `docs/kg-changes/non-de-evidence-extension.md`.
- Table S4 (relative enrichment = vesicle/cell fold-change with significance) maps cleanly to a differential-expression-style `csv` entry even though the "treatment" is a subcellular fraction rather than an environmental condition — DerivedMetric would also work but the log2FC + p-value structure better fits `statistical_analyses`.
- Both strains (MIT9312, MIT9313) are already in the KG — no organism work needed.
- MIT9312 uses PMT9312_#### locus tags; MIT9313 uses PMT#### — verify via `/check-gene-ids` before integration.
