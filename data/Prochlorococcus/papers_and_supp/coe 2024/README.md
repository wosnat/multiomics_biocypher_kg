# Coe 2024

**Paper:** Coe et al. (2024) — Dark-tolerant Prochlorococcus
**Journal:** ISME Journal
**DOI:** See paperconfig.yaml / pdf_extraction_cache.json

## Summary

This study investigates how Prochlorococcus NATL2A evolves dark tolerance through long-term co-culture with Alteromonas macleodii MIT1002. The authors evolved NATL2A for extended dark survival and performed RNA-seq time-course experiments comparing the dark-tolerant evolved strain to the parental (wild-type) strain, both in co-culture with MIT1002 under a 13:11 light:dark diel cycle. Expression was sampled every ~4 hours over a full 24h diel cycle (0, 4, 8, 13, 16, 20, 24h).

## Organisms

- **Prochlorococcus NATL2A** — primary organism (both dark-tolerant evolved and parental strains)
- **Alteromonas macleodii MIT1002** — co-culture partner

## Data in Knowledge Graph

### Differential Expression (Supplemental Tables 2 and 4)

Two experiments comparing dark-tolerant vs parental strains across 7 diel timepoints:

1. **NATL2A RNA-seq** (Supplemental Table 2): Dark-tolerant vs parental NATL2A gene expression at 0, 4, 8, 13, 16, 20, 24h. Uses DESeq2 with p-value threshold 0.1 (asterisk-in-logFC encoding). Gene IDs: NCBI old_locus_tag (PMN2A_RS*), locus_tag_ncbi, locus_tag_cyanorak, gene_name.

2. **MIT1002 RNA-seq** (Supplemental Table 4): Dark-tolerant vs parental co-culture effect on Alteromonas MIT1002 gene expression at the same 7 timepoints. Gene IDs: NCBI old_locus_tag.

### Gene Expression Clusters (Supplemental Table 3)

Two clusterings of NATL2A diel expression patterns (2081 genes, 15 clusters each):

1. **Dark-tolerant strain clusters** (`dark-tolerant_cluster` column): Diel periodicity clusters from the dark-tolerant evolved NATL2A strain, identified using RAIN periodicity detection.

2. **Parental strain clusters** (`parental_cluster` column): Diel periodicity clusters from the parental (wild-type) NATL2A strain, same method.

The CSV also contains a `same_periodicity_pattern` boolean column indicating whether each gene falls into the same cluster pattern in both strains.

Gene IDs in the cluster table are mixed format: gene names (e.g., AtpT, Yfr1), cds-prefixed locus tags (cds-PMN2A_RS*), cds-prefixed protein IDs (cds-WP_*), fig| IDs, and rna-prefixed IDs.

## Files

| File | Description |
|------|-------------|
| `paperconfig.yaml` | Full configuration for KG integration |
| `coe_isme_ycae131.pdf` | Main paper PDF |
| `Supplemental Table 2.csv` | NATL2A DE results (dark-tolerant vs parental, 7 timepoints) |
| `Supplemental Table 3.csv` | NATL2A diel expression cluster assignments (2081 genes, 2 clusterings) |
| `Supplemental Table 4.csv` | MIT1002 DE results (dark-tolerant vs parental co-culture, 7 timepoints) |

## Notes

- The cluster gene IDs use mixed formats requiring ID translation (cds- prefix stripping, fig| IDs, gene names)
- The `same_periodicity_pattern` column is metadata about cluster correspondence, not used directly by the KG
- Both clusterings reference the same experiment node since they derive from the same RNA-seq dataset comparing dark-tolerant vs parental strains
