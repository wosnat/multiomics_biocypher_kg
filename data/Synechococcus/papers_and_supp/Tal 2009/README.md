# Tai 2009 -- Synechococcus-Vibrio Interactions

## Citation

Tai V, Paulsen IT, Phillippy K, Johnson DA, Palenik B. (2009)
Whole-genome microarray analyses of Synechococcus-Vibrio interactions.
Environmental Microbiology. doi:10.1111/j.1462-2920.2009.01997.x

## Classification

**Bucket D -- defer / nothing to do (already integrated)**

WH8102 vs Vibrio parahaemolyticus LM5312 microarray coculture study; WH8102 is deployed in the KG (as Parasynechococcus WH8102) and Vibrio parahaemolyticus is set up as a treatment organism (taxid 670). The two split-direction tables (S1 up, S2 down) feed a single significant_only microarray experiment via paperconfig.yaml; both have `*_resolved.csv` siblings using the SYNW#### locus tag column. SAM scores are stored in the table but no adjusted p-values are reported, which is captured by the `significant_only` table_scope. Nothing more to add: the paper provides no metabolomics or non-DE per-gene tables.

## Summary

Whole-genome microarray analysis of Synechococcus sp. WH8102 gene expression
changes when co-cultured with the heterotroph Vibrio parahaemolyticus strain
LM5312. 285 genes significantly upregulated and 295 significantly downregulated
(SAM analysis, score thresholds >3.11 and <-3.09 respectively).

Key findings:
- Phosphate acquisition genes strongly upregulated despite phosphate not being depleted
- Cell wall biosynthesis and glycosyltransferase genes upregulated (potential operon)
- Zinc transporter genes upregulated
- Photosystem II, cytochrome b6f, ferredoxins, and antioxidant genes downregulated
  (possibly linked to Fur family iron regulation)
- NADH dehydrogenase genes downregulated
- Clp protease genes downregulated, secreted proteases upregulated
- Gene expression patterns share similarities with ammonia-grown WH8102

## Data Files

- `table s1.csv` — 285 significantly upregulated genes (SAM score > 3.11)
- `table s2.csv` — 295 significantly downregulated genes (SAM score < -3.09)

## Column Mapping

Note: The column names in the CSVs are confusingly swapped:
- `Gene Name` column contains the SYNW#### locus tags (used as name_col)
- `Gene ID` column contains gene product descriptions
- `Score` is the SAM statistic (not a p-value)
- `Log2 fold change` is log2 FC (coculture vs monoculture)

## Experimental Conditions

- Organism: Synechococcus sp. WH8102
- Coculture partner: Vibrio parahaemolyticus strain LM5312 (NCBI taxid 670)
- Medium: Artificial seawater with 9.0 mM NaNO3, vitamins
- Temperature: 25C
- Light: 30 uEinstein m-2 s-1 continuous white light
- RNA harvested at mid-log phase (~1x10^8 cells/mL)
- Two biological replicates, 7 and 3 technical replicates (dye-swapped)
- Statistical method: SAM (Significance Analysis of Microarrays)
- No adjusted p-values reported; significance based on SAM score thresholds

## Notes

- Table S2 has blank rows between functional categories — the omics adapter
  should handle these gracefully (rows without a valid Gene Name will be skipped)
- Vibrio parahaemolyticus (taxid 670) must be added to treatment_organisms.csv
  before this paper can be loaded into the KG
- Gene IDs are SYNW#### format — standard WH8102 locus tags, no mapping needed
