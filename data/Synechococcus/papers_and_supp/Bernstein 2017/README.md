# Bernstein 2017

**Paper:** Indirect Interspecies Regulation: Transcriptional and Physiological Responses of a Cyanobacterium to Heterotrophic Partnership

**Journal:** mSystems 2:e00181-16 (2017)

**DOI:** 10.1128/mSystems.00181-16

## Organisms

- **Primary:** Thermosynechococcus elongatus BP-1 (cyanobacterium)
- **Coculture partner:** Meiothermus ruber strain A (heterotroph)

**Note:** T. elongatus BP-1 is NOT currently in the KG genomes (`cyanobacteria_genomes.csv`). Gene IDs use NCBI locus tags (`SY28_RS#####` format). The organism would need to be added to the genome registry before this data can be integrated.

## Data available

This paper has clustering data only -- no differential expression (fold-change/p-value) tables.

| Data type | Status |
|-----------|--------|
| Differential expression | Not available (raw abundances only in data set s1) |
| Gene clustering (light) | Available -- 4 clusters (A-D) + unclustered, 4492 genes |
| Gene clustering (oxygen) | Available -- 4 clusters (E-H) + unclustered, 4492 genes |

## Files

| File | Description |
|------|-------------|
| `Synechococcus-Meiothermus co-culture mSystems 2017.pdf` | Full paper PDF |
| `data set s1 Exp_Data_BP1_AX_CC sys002172092sd1.csv` | Raw transcript abundances for axenic and coculture conditions (no fold-change; skipped) |
| `data set s2 clustered_Avg_Cond_Light sys002172092sd2.csv` | Clustered averaged transcript abundances by light condition (LLLO, MLLO, HLLO). Columns: ClustID_light (0-4), Clust name_light (A-D) |
| `data set s2 clusterd_Avg_Cond_Oxygen sys002172092sd2.csv` | Clustered averaged transcript abundances by oxygen condition (HLLO, HLMO, HLHO). Columns: ClustID_ox (0-4), Clust name_ox (E-H). Note filename typo: "clusterd" |
| `data set s2 Key for sample names sys002172092sd2.csv` | Key mapping condition IDs (HLLO, HLMO, etc.) to physical sample IDs and descriptions |
| `sys002172092sd1.xlsx` | Excel version of data set s1 |
| `sys002172092sd2.xlsx` | Excel version of data set s2 |
| `sys002172092s3.docx` | Supplementary text |
| `sys002172092sf4.pdf` | Supplementary figure 4 |
| `sys002172092st9.pdf` | Supplementary table 9 |
| `sys002172092st10.pdf` | Supplementary table 10 |
| `legends.txt` | Brief legends for data sets s1 and s2 |

## Paperconfig

The `paperconfig.yaml` defines two `gene_clusters` entries:
- **bp1_light_clusters** -- 4 clusters by light condition (low/medium/high light at low oxygen)
- **bp1_oxygen_clusters** -- 4 clusters by oxygen condition (low/medium/high oxygen at high light)

Cluster 0 in both cases represents unclustered genes (no cluster name assigned).

## Experimental conditions

Conditions are combinations of light intensity and dissolved oxygen:

| Abbreviation | Light | Oxygen | Used in |
|-------------|-------|--------|---------|
| LLLO | Low (30% power) | Low (pO2=0) | Light clustering |
| MLLO | Medium (68.64% power) | Low (pO2=0) | Light clustering |
| HLLO | High (100% power) | Low (pO2=0) | Light + oxygen clustering |
| HLMO | High (100% power) | Medium (pO2=0.30) | Oxygen clustering |
| HLHO | High (100% power) | High (pO2=0.59) | Oxygen clustering |

All conditions are binary coculture (T. elongatus BP-1 + M. ruber strain A).
