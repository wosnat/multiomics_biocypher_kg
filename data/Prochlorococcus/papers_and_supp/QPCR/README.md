# QPCR (reference directory)

**Citation:** Domínguez-Martín MA, López-Lozano A, Rangel-Zúñiga OA, Díez J, García-Fernández JM. Distinct features of C/N balance regulation in *Prochlorococcus* sp. strain MIT9313. FEMS Microbiology Letters, 2018, 365(3): fnx278.
**DOI:** 10.1093/femsle/fnx278
**Organism(s):** *Prochlorococcus* sp. MIT9313
**Topic:** Semi-quantitative RT-qPCR measurement of four genes (`icd`, `ntcA`, `glnB`, `pipX`) under N/P/Fe starvation, GS/GOGAT/electron-transport inhibitors, and culture ageing in MIT9313. No per-gene transcriptomics or DE table — the study reports relative expression fold-changes (2^−ΔΔCt) for a hand-picked 4-gene panel across several conditions, presented only in the main-text figures of the PDF.

## Classification

**Bucket D — defer / nothing to do**

PDF-only paper. No machine-readable supplementary tables exist; the entire dataset is four genes (`icd`, `ntcA`, `glnB`, `pipX`) charted across ~10 stress conditions in main-text figures. Integration would require manual figure curation of ~40 bar-chart values into a CSV. Low ROI — the same regulators are already represented in MIT9313 by transcriptomics from Tolonen 2006 / Voigt 2014, and MIT9313 is in the KG so the genes themselves are present. If a downstream consumer ever needs RT-qPCR comparisons specifically, this could graduate to bucket B as a small hand-curated DE CSV (~40 edges, one Experiment per stress); for now, defer.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `fnx278.pdf` | PDF | Main article — RT-qPCR of 4 regulators (icd/ntcA/glnB/pipX) under multiple stress conditions in MIT9313 | reference | — |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Identify the paper** — Done: `fnx278.pdf` is Domínguez-Martín et al. 2018 (FEMS Microbiol Lett, doi:10.1093/femsle/fnx278). Rename the directory from the generic `QPCR/` to `Dominguez-Martin 2018/` to match other paper directories, or leave in place if keeping as a reference.
2. **No action (skip integration)** — Paper reports RT-qPCR of a 4-gene targeted panel with no accompanying data table — all numeric values are in the figures only. Extraction would require manual transcription of ~40 bar-chart values (4 genes × ~10 conditions) from figures into a CSV. Low ROI: the genes studied (icd, ntcA, glnB, pipX) are already covered by transcriptomics papers in the KG (Tolonen 2006 for MED4/MIT9313 N-limitation; Voigt 2014 for MIT9313 nitrogen). Do **not** integrate unless a specific downstream consumer needs this data.
3. **Alternative (if integration is desired later)** — Manually curate a CSV with columns `gene`, `condition`, `log2_fold_change`, `adjusted_p_value` (values from Student t-test vs. control, recoverable from paper figure captions). Add as a small `csv` with `statistical_analyses`, one Experiment per stress condition (N, P, Fe starvation; GS/GOGAT/PET inhibition; culture ageing). Only 4 genes × ~10 conditions = ~40 edges — low volume but targeted on master regulators.

## Notes

- This is unparsed reference material: a PDF with no supplementary tables. Integration is possible only via manual figure curation.
- MIT9313 is already deployed in the KG. The 4 genes (`icd`, `ntcA`, `glnB`, `pipX`) all have stable MIT9313 locus tags (PMT_####) and will resolve cleanly if a CSV is ever hand-built.
- Directory name `QPCR/` is non-standard for this repo — consider moving to `Dominguez-Martin 2018/` for consistency with the `<Author Year>/` convention used elsewhere in `data/Prochlorococcus/papers_and_supp/`.
