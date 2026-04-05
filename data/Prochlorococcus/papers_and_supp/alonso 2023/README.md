# Alonso-Saez 2023

**Title:** Transcriptional Mechanisms of Thermal Acclimation in *Prochlorococcus*

**Authors:** Laura Alonso-Saez, Antonio S. Palacio, Ana M. Cabello, Semidan Robaina-Estevez, Jose M. Gonzalez, Laurence Garczarek, Angel Lopez-Urrutia

**Journal:** mBio, Volume 14, Issue 3, May/June 2023

**DOI:** 10.1128/mbio.03425-22

## Summary

Quantitative transcriptome (RNA-Seq) analysis of *Prochlorococcus marinus* MIT9301 under long-term thermal acclimation across its thermal niche (17-30C). Cultures were synchronized to a 12h:12h light/dark cycle and acclimated to five temperatures (17, 20, 22, 25, 30C). RNA samples were collected 3h after subjective sunrise (day) and sunset (night).

Five gene clusters (A-E) were identified by fuzzy c-means clustering based on day/night expression patterns across the thermal gradient:
- **Cluster A** (273 genes): Core daytime genes -- carbon fixation, Calvin cycle, glycogen synthesis, ATP synthesis, some PSII components. Consistently expressed during daytime.
- **Cluster B** (198 genes): Photosystems I & II, oxygen-evolving complex, regulatory proteins. Expression decreases from Topt toward Tmin during daytime, paralleling growth rate decline.
- **Cluster C** (432 genes): Stress response, regulatory proteins. Strongly upregulated at Tmin during daytime.
- **Cluster D** (340 genes): Stress response, glycogen degradation, protein synthesis. Upregulated at Tmin during both daytime and nighttime.
- **Cluster E** (437 genes): Pentose P pathway, respiration, DNA replication, cell division. Upregulated at nighttime across all temperatures.

## Organism

- *Prochlorococcus marinus* MIT9301 (HLII clade representative)

## Data files

| File | Description |
|------|-------------|
| `TABLE S5 softcluster membership.csv` | Soft cluster assignments and per-cluster probability scores for 1831 MIT9301 protein-coding genes |
| `mbio.03425-22-s0005.xlsx` | Table S3: mRNA transcripts per cell |
| `mbio.03425-22-s0006.xlsx` | Table S4: mRNA transcripts per biovolume |
| `soft_clusters_mbio.03425-22-s0008.xlsx` | Original Excel for Table S5 |

## Gene IDs

- `RefSeq Locus Tag` column: P9301_NNNNN format (MIT9301 locus tags) -- 1760/1831 genes (96.1%)
- `BV-BRC_ID` column: fig|167546.4.peg.NNN format (PATRIC/BV-BRC IDs) -- all 1831 genes
- 71 genes have BV-BRC IDs only (no RefSeq locus tag); these will be unmatched without additional id_translation

## Notes

- No differential expression data in this paper (no pairwise comparisons with fold changes/p-values). Data is absolute transcript counts across conditions, not DE results.
- The soft clustering uses fuzzy c-means, so each gene has a probability score for each of the 5 clusters. The `score_col` is not set in the paperconfig because there is no single score column; instead there are 5 per-cluster probability columns. A preprocessing step could derive the max probability as a single score column.
- The CSV has 2 extra header rows (title + blank line) before the actual column headers; `skip_rows: 2` handles this.
- Tables S3 and S4 contain raw transcript abundance data (not DE), so they are not configured as `csv` statistical_analyses entries.
