# Wang 2014 — med4_expression_level

## Cluster -- | mixed | high

**Name:** Prochlorococcus MED4 cluster -- (mixed, expression level)
**Enrichment:** not described (p=None)
**Functional:** This cluster consists of 89 genes that are unclassified in terms of expression level based on RNA-seq RPKM quartiles across 10 growth conditions. These genes do not fit into the defined expression categories (VEG, HEG, MEG, LEG, NEG) and thus represent a heterogeneous group with mixed expression patterns. Their functional roles are not specifically described in the paper.

**Behavioral:** Genes in this cluster show variable or inconsistent expression levels across the tested conditions, lacking a clear pattern of high or low expression.

**Notes:** Cluster is defined by exclusion from other expression level categories; no specific functional or behavioral description provided.

**Quotes:**
- [Page 4] All MED4 CDS genes were classified into five subgroups (HEG, MEG, LEG, NEG, and VEG).
- [Analysis description] Cluster --: 89 genes (sample: RNA_43, RNA_37, asRNA_00641)

## Cluster HEG | up | high

**Name:** Prochlorococcus MED4 cluster HEG (up, expression level)
**Enrichment:** COG categories C, J, O (p=0.001)
**Functional:** This cluster includes 291 highly expressed genes (HEG) consistently showing high transcript abundance across multiple growth conditions. These genes are significantly enriched in the core genome and are involved in essential cellular functions such as energy production and conversion, translation and ribosomal structure, and protein modification, folding, and turnover. They include ribosomal proteins, photosynthetic apparatus components, and protein folding machinery, reflecting their central metabolic roles and evolutionary conservation.

**Behavioral:** Genes in this cluster are consistently and highly expressed, evolve slowly, and are often organized in operons. They tend to have rapid mRNA turnover, which may enhance protein fidelity and cellular economy.

**Notes:** Strong evidence from RNA-seq data and functional enrichment analyses supports the characterization of this cluster.

**Quotes:**
- [Page 5] Among these core HEG genes, several functional categories were more prominent than others. These included the “C” (energy production and conversion), “J” (translation and ribosomal structure), and “O” (protein modification, folding and turnover) categories.
- [Page 4] HEG had a significantly lower nonsynonymous substitution rate (Ka) than MEG or LEG (Kruskal-Wallis Test, two-tailed P < 0.001).
- [Page 6] Operons are strikingly enriched in HEG and MEG.
- [Page 7] Highly expressed genes were more likely degraded fast.

## Cluster LEG | down | high

**Name:** Prochlorococcus MED4 cluster LEG (down, expression level)
**Enrichment:** COG category R (general function) (p=0.023)
**Functional:** This cluster comprises 82 lowly expressed genes (LEG) that show consistently low transcript levels across the tested growth conditions. These genes are less represented in the core genome and more common in the flexible genome. Functional enrichment indicates a slight overrepresentation of genes with general function prediction (COG category R). They tend to evolve faster and have longer mRNA half-lives compared to highly expressed genes.

**Behavioral:** Genes in this cluster are expressed at low levels, evolve relatively faster, and have slower mRNA degradation rates compared to HEG and MEG clusters.

**Notes:** Expression and evolutionary data support the low expression and faster evolution; functional enrichment is modest.

**Quotes:**
- [Page 5] Additionally, category “R” (general function) was slightly enriched in both LEG and NEG (P = 0.023 and 0.055).
- [Page 4] LEG had a significantly higher nonsynonymous substitution rate than HEG (Kruskal-Wallis Test, two-tailed P < 0.001).
- [Page 7] Lowly expressed genes were more likely slowly degraded.

## Cluster MEG | mixed | medium

**Name:** Prochlorococcus MED4 cluster MEG (mixed, expression level)
**Enrichment:** Essential genes (DEG-hit) enriched (p=0.001)
**Functional:** This cluster contains 445 moderately expressed genes (MEG) with intermediate transcript abundance across growth conditions. These genes are enriched in the core genome relative to the flexible genome and include a significant proportion of essential genes. Functional categories are diverse, reflecting a range of cellular processes. MEG genes evolve at intermediate rates and have mRNA turnover rates between those of HEG and LEG clusters.

**Behavioral:** Genes in this cluster show moderate and relatively stable expression, slower evolution than lowly expressed genes but faster than HEG, and intermediate mRNA degradation rates.

**Notes:** Moderate confidence based on expression, essentiality, and evolutionary data.

**Quotes:**
- [Page 4] Although the MEG subclass had a significantly higher rate of DEG-hit genes (P < 0.001).
- [Page 4] MEG had a significantly lower nonsynonymous substitution rate (Ka) than LEG (Kruskal-Wallis Test, two-tailed P < 0.001).
- [Page 6] Operons are strikingly enriched in HEG and MEG.

## Cluster NEG | down | high

**Name:** Prochlorococcus MED4 cluster NEG (down, expression level)
**Enrichment:** Phage-related genes (p=None)
**Functional:** This cluster includes 66 non-expressed genes (NEG) that show no detectable transcript levels under the tested growth conditions. These genes are predominantly found in the flexible genome and are often associated with phage-related functions. They evolve faster and have longer mRNA half-lives compared to expressed genes. Functional enrichment is low due to lack of expression and annotation.

**Behavioral:** Genes in this cluster are transcriptionally silent under tested conditions, evolve rapidly, and are mostly flexible genome or phage-related genes.

**Notes:** Expression data clearly define this cluster; functional annotation limited by lack of expression.

**Quotes:**
- [Page 5] As expected, phage-related genes displayed the lowest expression levels in this study, as phage infection conditions were not tested.
- [Page 4] The core genome had fewer NEG than the flexible genome (1.5% < 6.6%).

## Cluster VEG | mixed | medium

**Name:** Prochlorococcus MED4 cluster VEG (mixed, expression level)
**Enrichment:** not described (p=None)
**Functional:** This cluster consists of 1081 variably expressed genes (VEG) that show fluctuating expression levels across different growth conditions. These genes are more prevalent in the flexible genome and include genes with diverse functions. They tend to evolve faster than constantly expressed genes and have variable mRNA turnover rates. This cluster reflects genes with condition-dependent or variable regulation.

**Behavioral:** Genes in this cluster exhibit variable expression patterns, faster evolutionary rates, and are less likely to be organized in operons compared to constantly expressed genes.

**Notes:** Defined by variable expression; functional and behavioral heterogeneity expected.

**Quotes:**
- [Page 4] CEG subclass had a lower Ka than VEG (Mann–Whitney U Test, two-tailed P < 0.001).
- [Page 4] The core genome had fewer VEG than the flexible genome (49.6% < 64.6%).
- [Page 6] Operons are strikingly enriched in HEG and MEG, less so in VEG.
