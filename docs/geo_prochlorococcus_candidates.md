# GEO *Prochlorococcus* datasets — KG coverage triage

Source: `gds_result.txt` (GEO search dump, 694 records — only records 1–28 are Series; 29+ are Platforms/Samples).
Triaged 2026-08-16 against `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt`.

**28 series → 10 already covered · 8 out of scope · 10 series / 8 studies are candidates.**

---

## Candidates (not in the KG)

| GSE | First author / year | Title | DOI |
|---|---|---|---|
| GSE335434 | Bohutskyi / Johnson (PNNL), public 2026 | P-SSP7 infection of *Prochlorococcus* MED4 under standard light, high light, and darkness | — unpublished |
| GSE314951 | Johnson CGM 2026 (bioRxiv) | Spatiotemporal 4D whole-cell modeling of a minimal autotroph… | 10.64898/2026.01.21.700212 |
| GSE150732 | Huang S 2021, *MicrobiologyOpen* | Temporal transcriptomes of a marine cyanopodovirus and its *Synechococcus* host during infection | 10.1002/mbo3.1150 |
| GSE74921 + GSE63690 (SuperSeries GSE74922) | Doron S 2016, *ISME J* | Transcriptome dynamics of a broad host-range cyanophage and its hosts | 10.1038/ismej.2015.210 |
| GSE171435 | Hackl T 2023, *Cell* | Novel integrative elements and genomic plasticity in ocean ecosystems | 10.1016/j.cell.2022.12.006 |
| GSE17075 | Steglich C 2010, *Genome Biology* | Short RNA half-lives in the slow-growing marine cyanobacterium *Prochlorococcus* | 10.1186/gb-2010-11-5-r54 |
| GSE154594 | Muñoz-Marín MDC 2022, *Microbiology Spectrum* | Differential timing for glucose assimilation in *Prochlorococcus* and coexistent microbial populations in the North Pacific Subtropical Gyre | 10.1128/spectrum.02466-22 |
| GSE53065 | Voigt K 2014, *ISME J* | Comparative transcriptomics of two environmentally relevant cyanobacteria reveals unexpected transcriptome diversity | 10.1038/ismej.2014.57 |

### DOIs only (cut & paste)

```
10.64898/2026.01.21.700212
10.1002/mbo3.1150
10.1038/ismej.2015.210
10.1016/j.cell.2022.12.006
10.1186/gb-2010-11-5-r54
10.1128/spectrum.02466-22
10.1038/ismej.2014.57
```

### Notes per candidate

- **GSE335434** — 111 samples (MOI 10; 0/0.5/1/2/4/8/24 h + 26/32 h dark recovery; triplicates; uninfected controls). MED4 deployed. Raw counts CSV only → we run our own DE. Unpublished.
- **GSE314951** — 44 samples, 2 h resolution over a full cycle, 12:12 L:D, 22 °C. MED4 deployed. Counts CSV. The DOI is the *linked* citation (a whole-cell-modeling preprint using the data), not a dataset paper.
- **GSE150732** — WH7803 deployed; supplement already ships a WH7803 DE file → lowest friction. `treatment_type: [viral]`.
- **GSE74921 / GSE63690** — syn9 + P-TIM40 across WH7803, WH8102, WH8109 (+ NATL2A in the array arm). WH8109 would need onboarding.
- **GSE171435** — MIT0604 mitomycin C / UV. Unfiltered sense + antisense DE tables (better `table_scope` than most). Blockers: MIT0604 not in KG (`/add-a-strain`) + custom `ThomasGenomeGFF` IDs (`annotation_gff` + `id_translation`). GEO lists no citation; Hackl 2023 is inferred from the shared submitters/tycheposon topic.
- **GSE17075** — rifampicin RNA decay, MED4 + MIT9313, 7 timepoints. Not DE — per-gene scalar → `derived_metrics_table` (like Waldbauer 2012). GEO has CEL files only; the published half-life table is the easier source.
- **GSE154594** — field samples at station ALOHA on a custom Agilent array. Needs probe→gene mapping + a reference-proteome-match-style organism. Borderline.
- **GSE53065** — TSS / operons / UTRs / ncRNAs, MED4 vs MIT9313. No schema slot today. Skip unless a TSS layer is wanted.

**Suggested order:** GSE150732 → GSE314951 → GSE335434 → GSE74921/GSE63690 → GSE171435.

---

## Already covered

| GSE | KG paper |
|---|---|
| GSE264347 | Coe 2024 |
| GSE195946 | He 2022 |
| GSE118155 | Tetu 2019 |
| GSE93197 | Biller 2018 |
| GSE73511 | Biller 2016 |
| GSE79359 | Thompson 2016 |
| GSE65684 | Bagby and Chisholm 2015 |
| GSE49517 | Wang 2014 |
| GSE26533 | Thompson 2011 |
| GSE8382 | Lindell 2007 |

## Out of scope

- GSE65290 / GSE65278 / GSE65277 — mouse autoimmune peptide arrays (*Prochlorococcus* MIT9202 only as a peptide-antigen source).
- GSE149511 — *Synechocystis* PCC 6803 *nblD* mutant (strain not in KG).
- GSE130464 / GSE44448 / GSE21502 / GSE9384 — environmental community microarrays (MicroTOOLs, genome-proxy); no strain-level gene IDs to resolve.

---

**Scope caveat:** the GEO query was *Prochlorococcus*-anchored — *Alteromonas*- or *Synechococcus*-only datasets are not represented.
