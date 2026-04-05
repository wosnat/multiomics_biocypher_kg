# Ma 2022

**Title:** Cross-feeding between cyanobacterium *Synechococcus* and *Escherichia coli* in an artificial autotrophic-heterotrophic coculture system revealed by integrated omics analysis

**Journal:** Biotechnology for Biofuels and Bioproducts (2022) 15:69

**DOI:** 10.1186/s13068-022-02163-5

**Authors:** Jiajia Ma, Taohong Guo, Meijin Ren, Lei Chen, Xinyu Song, Weiwen Zhang

## Organism

- **Target organism:** *Synechococcus elongatus* UTEX 2973 (cscB+ sucrose-secreting variant)
- **Treatment organism:** *Escherichia coli* ABKm (3-HP producing strain)
- **Gene ID format:** `M744_#####` — old locus tags from the UTEX 2973 genome
- **NOTE:** *S. elongatus* UTEX 2973 is NOT currently in the KG genome list (`cyanobacteria_genomes.csv`). The genome must be added before this paper can be integrated.

## Experimental Design

Coculture of sucrose-secreting *S. elongatus* cscB+ and *E. coli* ABKm in CoBG-11 medium, separated by a dialysis bag (14 kDa cutoff). *E. coli* was inside the bag, *S. elongatus* outside. Compared to axenic *S. elongatus* cscB+ grown under the same conditions.

- **Medium:** CoBG-11 (BG-11 + 150 mM NaCl + 4 mM NH4Cl + 3 g/L TES, pH 8.3)
- **Temperature:** 30C
- **Light:** 100 umol photons m-2 s-1, continuous illumination
- **Duration:** 4 days coculture, sampled at 48h (transcriptomics/proteomics)
- **Replicates:** 3 biological replicates per condition

## Supplementary Tables

### Table S1 — Up-regulated transcripts (RNA-seq)
- 78 genes up-regulated in cocultured vs axenic *S. elongatus* cscB+
- Columns: Gene ID, Description/Gene name, Pathway, log2FoldChange
- Cutoff: fold change > 1.5 (log2FC > ~0.585) AND padj (Q-value) < 0.05
- No p-values in table (pre-filtered)
- Transcriptomic analysis by GENEWIZ (Suzhou, China)

### Table S2 — Down-regulated transcripts (RNA-seq)
- 55 genes down-regulated in cocultured vs axenic *S. elongatus* cscB+
- Same columns and cutoffs as S1
- No p-values in table (pre-filtered)

### Table S4 — Selected differentially expressed proteins (proteomics)
- ~40 selected proteins associated with cross-feeding and metabolite exchange
- Columns: Protein, Description, Mean Ratio
- Method: iTRAQ quantitative proteomics via LC-MS/MS (BGI, Shenzhen)
- Cutoff: fold change > 1.2 AND P-value < 0.05 (t-test of 9 pairwise comparison groups)
- No p-values in table (pre-filtered)
- **WARNING:** Mean Ratio uses non-standard sign convention for down-regulated proteins (negative values like -0.75 instead of 0.75 or log2-transformed values). Up-regulated values are standard linear ratios (>1). This ambiguity should be resolved before integration.

## Data Limitations

1. **No p-values in any table** — all tables are pre-filtered for significance
2. **Split direction tables (S1/S2)** — RNA-seq results split into up and down tables, both referencing the same experiment
3. **Proteomics fold change ambiguity** — Table S4 Mean Ratio has negative values for down-regulated proteins, which is not standard for linear fold changes
4. **Organism not in KG** — *S. elongatus* UTEX 2973 needs to be added to `cyanobacteria_genomes.csv` with its NCBI assembly accession before this data can be loaded
5. **Only S. elongatus DE data** — no differential expression data for the *E. coli* side
