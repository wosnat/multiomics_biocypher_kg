# Wang 2014

**Paper:** Wang D, Ning K, Li J, Hu J, Han D, Wang H, Zeng X, Jing X, Zhou Q, Su X, Chang X, Wang A, Wang W, Jia J, Wei L, Xin Y, Qiao Y, Huang R, Chen J, Han B, Yoon K, Hill RT, Zohar Y, Chen F, Hu Q, Xu J. (2014) "Nannochloropsis Genomes Reveal Evolution of Microalgal Oleaginous Traits." *BMC Genomics* 14:11.

**Note:** Despite being filed under `papers_and_supp/`, the supplementary data extracted here pertains to *Prochlorococcus* MED4 gene expression levels and operon predictions from an earlier related analysis by the same group. The BMC article number (1471-2180-14-11) indicates *BMC Microbiology*.

## Organism

*Prochlorococcus marinus* MED4

## Available Data Files

| File | Description |
|------|-------------|
| `GeneClassification 12866_2013_2178_MOESM3_ESM.csv` | Classification of all MED4 CDS genes by expression level (VEG/MEG/HEG), with Core/Flexible designation, COG annotations |
| `operons 12866_2013_2178_MOESM1_ESM.csv` | Predicted operon structure for MED4 (operon ID, coordinates, member genes) |
| `12866_2013_2178_MOESM1_ESM.xlsx` | Supplementary Table 1 (xlsx version of operon predictions) |
| `12866_2013_2178_MOESM2_ESM.xlsx` | Supplementary Table 2 |
| `12866_2013_2178_MOESM3_ESM.xlsx` | Supplementary Table 3 (xlsx version of gene classification) |
| `12866_2013_2178_MOESM4_ESM.xlsx` | Supplementary Table 4 |
| `1471-2180-14-11.pdf` | Paper PDF |

### Gene Classification CSV columns

`sysName`, `locus tag`, `start`, `stop`, `strand`, `Core/Flexible`, `ExpressionLevel`, `name`, `desc`, `COG`, `COGFun`, `COGDesc`

- **ExpressionLevel** categories: VEG (very highly expressed genes), MEG (moderately expressed genes), HEG (highly expressed genes)
- Gene IDs use PMM-style sysNames (e.g., PMM0001) and PMED4_ locus tags

### Operons CSV columns

`Operon_ID`, `start`, `end`, `strand`, `length`, `Num. of genes`, followed by variable-width gene member columns (sysName format)

## KG Integration

**Gene expression level classification** integrated as `gene_clusters` entry via `GeneClassification_clean.csv` (preprocessed: title row removed, footnote suffixes stripped from column names).

- `gene_id_col: sysName` (PMM#### format)
- `cluster_col: ExpressionLevel` — VEG (1081), MEG (445), HEG (291), LEG (82), NEG (66), -- (89), empty (5)
- `cluster_type: expression_level`
- `cluster_method`: Quartile-based classification from RNA-seq RPKM across 10 growth conditions (Pro99 and AMP media, 21C, continuous light 28 umol quanta m-2 s-1)

**Operons** not integrated — operon structure predictions don't map to DE edges or gene clusters.

**Note:** `GeneClassification_clean.csv` is derived from the original CSV with cleaned headers (title row removed, footnote suffixes stripped). Original file preserved.
