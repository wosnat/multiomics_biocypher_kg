# Szul 2019

**Citation:** Szul MJ, Dearth SP, Campagna SR, Zinser ER. Carbon Fate and Flux in *Prochlorococcus* under Nitrogen Limitation. mSystems, 2019, 4(1): e00254-18.
**DOI:** 10.1128/mSystems.00254-18
**Organism(s):** *Prochlorococcus* sp. VOL29 (eMED4/HL-I ecotype; axenic isolate). **For KG integration: mapped to MIT9312** (HL-I ecotype, conspecific).
**Topic:** Compares carbon fixation rates, intracellular metabolite pools, and ¹³C-isotope labeling between N-limited (chemostat) and N-replete (batch) VOL29 cultures; includes cross-comparison with field metabolomics data. This is a metabolomics + physiology study — no per-gene transcriptomics or proteomics data.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `mSystems.00254-18.pdf` | PDF | Main article | reference | — |
| `mSystems.00254-18-st001.pdf` | PDF | Rendered PDF of Supplementary Table S1 | reference | — |
| `mSystems.00254-18-st001.docx` | DOCX | Supp Table S1 (machine-readable): likely growth rates / physiological measurements per culture | skip | Sample/culture-level, not gene-level |
| `mSystems.00254-18-st002.docx` | DOCX | Supp Table S2: metabolite pool concentrations per condition (intracellular) | skip | Metabolomics — no gene targets |
| `mSystems.00254-18-st003.docx` | DOCX | Supp Table S3: ¹³C-isotope labeling measurements per metabolite / timepoint | skip | Metabolomics / isotope tracing — no per-gene evidence |

## Classification

**Bucket B — metabolites actionable** (Metabolite + MetaboliteAssay support has landed; Phase 2 metabolomics). Map VOL29 → MIT9312.

Paper is purely metabolomics + ¹³C tracing on a Prochlorococcus VOL29 isolate. Tables S1-S3 supply intracellular metabolite pools, ¹³C isotopologue labeling, and field cell-counts/sampling metadata; no transcriptomics, no proteomics, no per-gene evidence of any kind. Both blockers have cleared: (a) Metabolite/MetaboliteAssay nodes exist, and (b) the VOL29 isolate is treated as MIT9312 for KG purposes (eMED4/HL-I ecotype, conspecific). Table S2 (intracellular pools) maps cleanly to `metabolite_assays_table` with `organism: "Prochlorococcus MIT9312"` and `compartment: whole_cell` → `Assay_quantifies_metabolite`. ¹³C isotopologue labeling (S3) does NOT fit the current MetaboliteAssay schema (no flux/labeling fields) and is deferred.

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Add Table S2** — Intracellular metabolite pool concentrations → `metabolite_assays_table` with `organism: "Prochlorococcus MIT9312"` and `compartment: whole_cell`. Single Experiment per (N-status × timepoint). VOL29 is mapped to MIT9312 (eMED4/HL-I ecotype, conspecific).
2. **Defer Table S3** — ¹³C isotopologue labeling has no native fit in the current MetaboliteAssay schema (no flux/labeling fields). Skip for now or open a separate spec.
3. **Skip Table S1** — culture-level physiological measurements; not gene-level and not metabolite-quant. Out of scope.

## Notes

- Three DOCX files: st001 is also supplied as a PDF (redundant). The st002 and st003 tables are metabolite-level and do not link to biosynthesis genes in the source files.
- Paper cites related work on N-regulation gene expression (Tolonen 2006, Grzymski/Dussaq 2011, Domínguez-Martín 2014) — those integrations belong to their own paper directories, not here.
- VOL29 → MIT9312 mapping: VOL29 is described as eMED4/HL-I; MIT9312 is the deployed HL-I representative. Pool concentrations from VOL29 will be attached to the MIT9312 OrganismTaxon node — note this in the paperconfig `organism` and in the Experiment's `experimental_context` for traceability.
