# Beliaev 2014

**Paper:** Beliaev AS et al. "Inference of interactions in cyanobacterial-heterotrophic co-cultures via transcriptome sequencing." *The ISME Journal* (2014) 8, 2243-2255. doi:10.1038/ismej.2014.69

**GEO accession:** GSE53360

## Organisms

- **Synechococcus sp. PCC 7002** (NCBI taxid 32049, GenBank NC_010475) -- euryhaline unicellular cyanobacterium
- **Shewanella putrefaciens W3-18-1** (NCBI taxid 351745, GenBank NC_008750) -- marine facultative aerobe heterotroph

Both organisms are NEW to the KG and require genome setup in Phase 3.

## Experimental Design

Chemostat co-cultures grown under continuous light at 30C, carbon-limited steady state.
Three conditions compared:
1. **Axenic controls** -- each organism grown alone (Syn 7002 on 8 mM NaHCO3; Shew W3-18-1 on 5 mM D,L-Na lactate)
2. **Co-culture with lactate** -- 5 mM D,L-Na lactate as carbon source, 1720 umol photons m-2 s-1
3. **Co-culture with HCO3-** -- 8 mM NaHCO3 as carbon source, 1720 umol photons m-2 s-1

Sequencing: SOLiD 5500XL, 50-base reads, mapped with SOLiDTM LifeScope v2.5.
Normalization: RPKM. No statistical test/p-values reported -- fold change threshold of >=2-fold used.

## Supplementary Tables

All tables have identical column structure with 3 fold change columns (lactate vs axenic, HCO3- vs axenic, HCO3- vs lactate). Tables are filtered by >=2-fold threshold on specific comparisons.

| Table | Organism | Filter criteria | Data rows | Used |
|-------|----------|----------------|-----------|------|
| S1 | Synechococcus 7002 | >=2-fold ALL conditions vs axenic | 475 | No (subset of S3 and S5) |
| S2 | Shewanella W3-18-1 | >=2-fold ALL conditions vs axenic | 236 | No (subset of S4 and S6) |
| S3 | Synechococcus 7002 | >=2-fold lactate vs axenic | 1186 | Yes |
| S4 | Shewanella W3-18-1 | >=2-fold lactate vs axenic | 458 | Yes |
| S5 | Synechococcus 7002 | >=2-fold HCO3- vs axenic | 692 | Yes |
| S6 | Shewanella W3-18-1 | >=2-fold HCO3- vs axenic | 1343 | Yes |
| S7 | Synechococcus 7002 | >=2-fold HCO3- vs lactate co-culture | 702 | Yes |
| S8 | Shewanella W3-18-1 | >=2-fold HCO3- vs lactate co-culture | 1278 | Yes |

## Paperconfig Design

6 experiments (3 comparisons x 2 organisms), 6 supplementary tables used (S3-S8).

Tables S1/S2 are skipped because they are strict subsets:
- S1 genes are almost entirely contained in both S3 and S5 (473/474 overlap)
- S2 genes are fully contained in both S4 and S6 (235/235)

## Key Characteristics

- **Fold changes are LINEAR** (not log2) -- values like 0.46, 3.32, 16.667. Uses `fold_change_type: linear`.
- **No p-values** -- tables only report fold changes with a >=2-fold cutoff.
- **Gene ID format:** `SYNPCC7002_A####` (Synechococcus), `SputW3181_####` (Shewanella)
- **Locus tag column** has a leading space: `" Locus tag"` (note space before L)
- **Both organisms have DE data** -- this is a bilateral coculture transcriptome study

## Notes

- The paper also includes NMR metabolomics data (Table 2) and pathway enrichment (Tables 3-6) which are not captured in the paperconfig.
- A light-limited condition (640 umol photons m-2 s-1, HCO3-) is mentioned in Table 1 but no separate DE table is provided for it.
