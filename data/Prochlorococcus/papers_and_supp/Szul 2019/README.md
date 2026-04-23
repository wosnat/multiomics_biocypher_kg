# Szul 2019

**Citation:** Szul MJ, Dearth SP, Campagna SR, Zinser ER. Carbon Fate and Flux in *Prochlorococcus* under Nitrogen Limitation. mSystems, 2019, 4(1): e00254-18.
**DOI:** 10.1128/mSystems.00254-18
**Organism(s):** *Prochlorococcus* sp. VOL29 (eMED4/HL-I ecotype; axenic isolate)
**Topic:** Compares carbon fixation rates, intracellular metabolite pools, and ¹³C-isotope labeling between N-limited (chemostat) and N-replete (batch) VOL29 cultures; includes cross-comparison with field metabolomics data. This is a metabolomics + physiology study — no per-gene transcriptomics or proteomics data.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `mSystems.00254-18.pdf` | PDF | Main article | reference | — |
| `mSystems.00254-18-st001.pdf` | PDF | Rendered PDF of Supplementary Table S1 | reference | — |
| `mSystems.00254-18-st001.docx` | DOCX | Supp Table S1 (machine-readable): likely growth rates / physiological measurements per culture | skip | Sample/culture-level, not gene-level |
| `mSystems.00254-18-st002.docx` | DOCX | Supp Table S2: metabolite pool concentrations per condition (intracellular) | skip | Metabolomics — no gene targets |
| `mSystems.00254-18-st003.docx` | DOCX | Supp Table S3: ¹³C-isotope labeling measurements per metabolite / timepoint | skip | Metabolomics / isotope tracing — no per-gene evidence |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Skip** — Paper has no per-gene differential expression, no proteomics, no gene-level evidence of any kind. The data is physiological rates + metabolite pool sizes + ¹³C labeling — none of which maps to Gene or Protein nodes in the current schema.
2. **No action** — Do not create a paperconfig. Keep the files in-place for provenance.
3. **Add organism check** — VOL29 is **not** in the KG (not listed among current strains in CLAUDE.md). The paper describes it as an eMED4/HL-I ecotype isolate; if future integration were desired, VOL29 would need a genome entry in `cyanobacteria_genomes.csv`. Since there is no gene-level data to integrate here, no urgency.

## Notes

- Three DOCX files: st001 is also supplied as a PDF (redundant). The st002 and st003 tables are metabolite-level and do not link to biosynthesis genes in the source files.
- Paper cites related work on N-regulation gene expression (Tolonen 2006, Grzymski/Dussaq 2011, Domínguez-Martín 2014) — those integrations belong to their own paper directories, not here.
- Similar metabolomics-only status as Kujawinski 2023 and McDonagh 2012 (for non-redox content). If a future schema extension adds Metabolite nodes with Gene-to-Metabolite biosynthesis links, revisit.
