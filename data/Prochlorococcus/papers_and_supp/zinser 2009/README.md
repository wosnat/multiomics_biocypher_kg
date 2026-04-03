# Zinser 2009

**Citation:** Zinser ER, Lindell D, Johnson ZI, Futschik ME, Steglich C, et al. (2009)
Choreography of the Transcriptome, Photophysiology, and Cell Cycle of a Minimal
Photoautotroph, *Prochlorococcus*. PLoS ONE 4(4): e5135.
doi:10.1371/journal.pone.0005135

## Summary

Landmark diel transcriptome study of Prochlorococcus MED4 over a simulated
14:10 light:dark cycle. 80% of annotated genes showed cyclic expression.
Sixteen expression clusters were identified by soft clustering (Mfuzz),
plus two special clusters: aperiodic/expressed (17) and non-expressed (18).

## Organism

- Prochlorococcus MED4 (axenic)

## Experimental Conditions

- Medium: Sargasso Seawater-based Pro99, 10 mM HEPES pH 7.5
- Temperature: 24 +/- 0.2 C
- Light: 14h light / 10h dark with gradual dawn/dusk transitions
- Max light intensity: ~232 umol quanta m-2 s-1 (10:00-14:00 experimental time)
- Doubling time: ~1 day (1.0-1.1 days for replicate cultures)
- Sampling: every 2 hours over 48 hours (50 timepoints), duplicate cultures sampled alternately

## Platform

- Custom Affymetrix MD4-9313 GeneChip (includes MED4 and MIT9313 probes)
- Normalization: Robust Multi-Array Average (RMA)

## Data in Table_S1.csv

3610 rows (1697 ORFs + 1855 intergenic regions + 58 non-coding RNAs).
Columns:
- Gene or region: PMM#### (ORFs), PMMIG... (intergenic), PMM_... (ncRNAs)
- New gene name: PMED4_##### (updated nomenclature)
- Function: annotation
- Fourier: Fourier score for periodicity
- FDR: false discovery rate for Fourier score
- Peak: peak expression time (hours)
- Peak R value: Pearson correlation with shifted cosine
- Cluster: cluster assignment (1-16 periodic, 17 aperiodic, 18 non-expressed, NA)
- Cluster membership score: soft clustering membership (0-1)
- Cyanobase categories: functional category and sub-category
- 50 pairs of Normalized/StdDev columns (timepoints 0-48h at 2h intervals)

## Clusters

18 clusters total:
- Clusters 1-16: periodic expression patterns ordered by peak time
  - Cluster 1 (57 genes): peak ~08:20, enriched in Photosystem I and II
  - Cluster 13 (110 genes): peak ~03:24, enriched in ribosomal proteins
  - Cluster 16 (125 genes): peak ~05:30, enriched in ATP synthase and CO2 metabolism
- Cluster 17 (173 genes): aperiodic but expressed
- Cluster 18 (180 genes): non-expressed

## KG Integration

- gene_clusters entry only (no DE edges -- data are raw expression values, not fold-changes)
- Gene IDs are PMM#### format (standard MED4 locus tags) -- direct match expected
- Intergenic (PMMIG) and ncRNA (PMM_) probes will not match gene nodes (expected)
- Cluster "NA" rows have no cluster assignment and will be skipped by the adapter
