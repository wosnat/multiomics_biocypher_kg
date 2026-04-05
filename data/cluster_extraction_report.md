# Cluster Extraction Report

## Alonso-Saez 2023 / mit9301_softclusters_thermal_acclimation

### Cluster Cluster A | mixed | high
**Name:** Prochlorococcus cluster A (mixed, core metabolism)
**Enrichment:** Core metabolic pathways (p=None, sig=False)
**Functional:** Includes genes related to carbon fixation and assimilation such as RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), Calvin cycle enzymes (gap2, tktA, glpX, pgk, cbbA), glycogen synthesis (glgABC), ATP synthesis genes (atpADE), and some photosystem II components (psbA, psbC, psbD). These genes represent essential metabolic pathways maintained across the thermal niche.
**Behavioral:** Genes are consistently expressed during daytime across all temperatures with stable expression levels, indicating core daytime metabolic activity.
**Notes:** High confidence due to clear gene annotations and consistent expression patterns.
**Sources:** Figure 2, Table S4, Table S5

### Cluster Cluster B | down | high
**Name:** Prochlorococcus cluster B (down, photosynthesis)
**Enrichment:** Photosynthesis (p=None, sig=False)
**Functional:** Enriched for photosynthetic genes including photosystem I components (psaABDEFKL), some photosystem II components (psbBJH, psbO), and photosynthetic electron transport genes (petACGNM). Also includes phosphate binding protein genes (pstS).
**Behavioral:** Genes show a gradual decrease in expression from optimum to minimum temperature during daytime, paralleling growth rate decline. This cluster exhibits transcriptional suppression of photosynthetic apparatus under cold conditions and altered circadian expression patterns.
**Notes:** Strong evidence from expression patterns and correlation with growth rates.
**Sources:** Figure 2, Figure 4, Figure 5

### Cluster Cluster C | up | high
**Name:** Prochlorococcus cluster C (up, cold stress response)
**Enrichment:** Stress response and protein synthesis (p=None, sig=False)
**Functional:** Enriched for global stress response genes including cellular chaperones (groES/groEL, dnaK, clpBCP), fatty acid desaturases (desA, desC), oxidative damage protection genes (recA, ruvB, sod), carotenoid synthesis (pds, crtBH), and rubredoxin (rub). Also includes amino acid synthesis genes (glyA, serA, leuA), translation initiation factors (infABC), and nitrogen acquisition genes.
**Behavioral:** Genes are strongly upregulated at the minimum temperature (17°C) during both daytime and nighttime, with some chaperones showing prioritization during daytime. This reflects activation of cold stress and oxidative stress responses.
**Notes:** Clear functional enrichment and strong temperature-dependent expression.
**Sources:** Figure 2, Figure 3, Table S4, Table S5

### Cluster Cluster D | up | medium
**Name:** Prochlorococcus cluster D (up, cold stress response night)
**Enrichment:** Stress response and cell cycle (p=None, sig=False)
**Functional:** Includes genes involved in stress response, protein synthesis, glycogen degradation (glgP), and cell division (ftsZYQ). Also contains DNA replication genes (dnaA, nrdJ, gyrB) and pentose phosphate pathway genes (tal, gnd, zwf).
**Behavioral:** Genes are upregulated at the minimum temperature during both daytime and nighttime, with some showing increased expression at night, reflecting cold stress adaptation and increased demand for energy and protein synthesis.
**Notes:** Functional roles inferred from gene annotations and expression timing.
**Sources:** Figure 2, Figure 3, Table S4

### Cluster Cluster E | mixed | high
**Name:** Prochlorococcus cluster E (mixed, core nighttime metabolism)
**Enrichment:** Core metabolic pathways (p=None, sig=False)
**Functional:** Includes genes related to catabolic consumption (cyoB, ndhD), DNA replication (dnaA, nrdJ, gyrB), cell division (ftsZYQ), and pentose phosphate pathway (tal, gnd, zwf). These genes represent essential pathways typically upregulated at nighttime.
**Behavioral:** Genes are upregulated at nighttime across the thermal gradient, maintaining core nighttime metabolic functions.
**Notes:** Consistent with known nighttime expression patterns.
**Sources:** Figure 2, Table S5

## Bernstein 2017 / bp1_light_clusters

### Cluster 0 | mixed | low
**Name:** Thermosynechococcus elongatus cluster 0 (mixed, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 1 | up | low
**Name:** Thermosynechococcus elongatus cluster 1 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 2 | down | low
**Name:** Thermosynechococcus elongatus cluster 2 (down, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 3 | up | low
**Name:** Thermosynechococcus elongatus cluster 3 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 4 | down | low
**Name:** Thermosynechococcus elongatus cluster 4 (down, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

## Bernstein 2017 / bp1_oxygen_clusters

### Cluster 0 | up | low
**Name:** Thermosynechococcus elongatus cluster 0 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 1 | up | low
**Name:** Thermosynechococcus elongatus cluster 1 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 2 | up | low
**Name:** Thermosynechococcus elongatus cluster 2 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 3 | up | high
**Name:** Thermosynechococcus elongatus cluster 3 (up, photosynthesis and metabolism)
**Enrichment:** Photosynthesis and metabolism (p=0.05, sig=True)
**Functional:** Enriched for photosystem II components (psbV, psbX, psbV2), carboxysome genes (ccmK1, ccmL), and pyruvate metabolism genes (pdhB, pdhA). Also includes genes involved in amino acid metabolism (cysE, metE) and ROS detoxification (grxD, sodB).
**Behavioral:** Genes in this cluster show increased expression with increasing irradiance and decreased expression with increasing oxygen tension, corresponding to increased specific growth and photosynthesis rates. Expression is sustained across steady states.
**Notes:** Strong functional enrichment and consistent expression patterns with physiological data.
**Sources:** Figure 3A, Figure 6D

### Cluster 4 | down | high
**Name:** Thermosynechococcus elongatus cluster 4 (down, metabolic exchange and oxidative stress)
**Enrichment:** Metabolic exchange and oxidative stress (p=0.05, sig=True)
**Functional:** Includes genes involved in organic acid synthesis and export (acs, ackA, ldhA, gltA), nitrogen metabolism (nrtABD, glnA, nirA), vitamin B12 biosynthesis and salvage (cobWNT, btuCD), and ROS detoxification (flv4, prxQ-B2, sodB). Also includes M. ruber genes for carbon and nitrogen uptake and metabolism, indicating metabolic coupling.
**Behavioral:** Genes in this cluster generally increase expression with increasing irradiance and decrease with increasing oxygen tension, correlating with higher growth and photosynthesis rates. This pattern suggests coordinated metabolic exchange and acclimation to heterotrophic partnership.
**Notes:** Strong evidence from coexpression of genes from both species and functional enrichment.
**Sources:** Figure 6D, Figure 6G, Figure 3C

## Bernstein 2017 / mruber_light_clusters

### Cluster 0 | mixed | low
**Name:** Meiothermus ruber cluster 0 (mixed, response_pattern)
**Enrichment:**  (p=None, sig=False)
**Functional:** Cluster 0 contains 857 genes including SY28_RS00050, SY28_RS00500, and SY28_RS04985. Not specifically discussed in the paper for functional enrichment or highlighted genes.
**Behavioral:** Not discussed in paper.
**Notes:** No specific functional or behavioral description provided in the paper for this cluster.

### Cluster 1 | up | high
**Name:** Meiothermus ruber cluster 1 (up, response_pattern)
**Enrichment:** Carbon and nitrogen metabolism (p=None, sig=False)
**Functional:** Cluster 1 includes 324 genes such as SY28_RS00505, SY28_RS05035, and SY28_RS05040. Genes involved in uptake and metabolism of cyanobacterium-derived organic carbon and nitrogen, including acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), xylose isomerase (xylA), and branched-chain amino acid uptake (livGH).
**Behavioral:** Genes in this cluster show increased expression with increasing irradiance (I_i) and decreased expression with increasing oxygen tension (pO2), corresponding to increased specific growth and photosynthesis rates. This pattern suggests coordination with cyanobacterial physiology to accommodate resource availability.
**Notes:** Functional roles inferred from gene annotations and expression patterns correlated with irradiance and oxygen treatments.
**Sources:** Figure 6

### Cluster 2 | down | medium
**Name:** Meiothermus ruber cluster 2 (down, response_pattern)
**Enrichment:** Vitamin B12 and amino acid metabolism (p=None, sig=False)
**Functional:** Cluster 2 contains 188 genes including SY28_RS05090, SY28_RS05091, and SY28_RS05110. Includes genes involved in vitamin B12 uptake/scavenging (btuCD) and methionine biosynthesis (metHX), which decrease with increasing irradiance. Also includes genes related to cysteine biosynthesis (cysK, cysE) indicating M. ruber likely synthesizes cysteine via a vitamin B12-independent pathway.
**Behavioral:** Genes in this cluster show decreased expression with increasing irradiance and increased expression with increasing oxygen tension, opposite to some cyanobacterial genes, suggesting a coordinated exchange of vitamins and amino acids with T. elongatus.
**Notes:** Functional roles inferred from gene annotations and opposing expression patterns with cyanobacterial homologs.
**Sources:** Figure 6

### Cluster 3 | up | medium
**Name:** Meiothermus ruber cluster 3 (up, response_pattern)
**Enrichment:** Oxidative stress response and electron transfer (p=None, sig=False)
**Functional:** Cluster 3 includes 259 genes such as SY28_RS05100, SY28_RS05115, and SY28_RS05120. Contains genes involved in ROS detoxification and electron transfer processes, including peroxidases (bcp), superoxide dismutase (sod2), and electron transfer flavoprotein (fixAB).
**Behavioral:** Genes in this cluster generally increase with increasing oxygen tension (pO2) and decrease with increasing irradiance, corresponding to increased oxidative stress and electron transfer activity under higher oxygen conditions.
**Notes:** Functional roles inferred from gene annotations and expression patterns correlated with oxygen treatments.
**Sources:** Figure 6

### Cluster 4 | mixed | low
**Name:** Meiothermus ruber cluster 4 (mixed, response_pattern)
**Enrichment:**  (p=None, sig=False)
**Functional:** Cluster 4 contains 444 genes including SY28_RS05005, SY28_RS05010, and SY28_RS05045. Not specifically discussed in the paper for functional enrichment or highlighted genes.
**Behavioral:** Not discussed in paper.
**Notes:** No specific functional or behavioral description provided in the paper for this cluster.

## Bernstein 2017 / mruber_oxygen_clusters

### Cluster 0 | up | low
**Name:** Meiothermus ruber cluster 0 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 1 | up | medium
**Name:** Meiothermus ruber cluster 1 (up, carbon and nitrogen metabolism)
**Enrichment:** Carbon and nitrogen metabolism (p=None, sig=False)
**Functional:** Includes genes involved in carbon uptake and metabolism such as acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), xylose isomerase (xylA), and amino acid uptake (livGH). Also includes nitrogen metabolism genes such as glutamate dehydrogenase (gdhA) and amino acid uptake systems.
**Behavioral:** Genes in this cluster increase with increasing specific growth and photosynthesis rates, corresponding to increasing irradiance and decreasing oxygen tension.
**Notes:** Descriptions based on coexpression patterns and gene annotations.
**Sources:** Figure 6, Table S1, Table S2

### Cluster 2 | up | medium
**Name:** Meiothermus ruber cluster 2 (up, vitamin B12 and methionine metabolism)
**Enrichment:** Vitamin B12 and methionine metabolism (p=None, sig=False)
**Functional:** Contains genes related to vitamin B12 uptake/scavenging (btuCD) and methionine biosynthesis and degradation (metH, metE, metK, mtnA). Also includes genes involved in cysteine biosynthesis (cysK, cysE).
**Behavioral:** Genes in this cluster generally decrease with increasing irradiance and increase with increasing oxygen tension, showing opposite patterns to some T. elongatus genes, suggesting exchange of vitamins and amino acids between species.
**Notes:** Based on coexpression patterns and gene function annotations.
**Sources:** Figure 6, Table S1, Table S2

### Cluster 3 | down | medium
**Name:** Meiothermus ruber cluster 3 (down, oxidative stress and electron transfer)
**Enrichment:** Oxidative stress and electron transfer (p=None, sig=False)
**Functional:** Includes genes encoding peroxidases (bcp), superoxide dismutase (sod2), and electron transfer flavoproteins (fixAB), NADH-dehydrogenase, and NADH-quinone oxidoreductase components (nuoDFGHIJKN).
**Behavioral:** Genes in this cluster generally increase with oxygen tension and decrease with irradiance, showing coordinated expression with T. elongatus genes involved in ROS detoxification.
**Notes:** Based on gene function and coexpression patterns.
**Sources:** Figure 6, Table S1, Table S2

### Cluster 4 | down | low
**Name:** Meiothermus ruber cluster 4 (down, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

## Biller 2018 / mit1002_periodicity

### Cluster coculture_LD | mixed | medium
**Name:** Alteromonas macleodii MIT1002 cluster coculture_LD (mixed, periodicity)
**Enrichment:** Not discussed in paper. (p=None, sig=False)
**Functional:** This cluster includes 529 genes showing 24-h periodicity only in coculture under the diel light:dark (L:D) cycle. Functional enrichment is not explicitly detailed in the paper for this cluster, but it includes genes such as MIT1002_0002, MIT1002_0007, and MIT1002_0010 as examples.
**Behavioral:** Genes in this cluster exhibit 24-h periodic oscillations in coculture during the diel L:D cycle but do not show periodicity under extended darkness or in axenic culture conditions.
**Notes:** Functional enrichment not described; gene examples provided.
**Sources:** Figure 4A, Table S4A

### Cluster coculture_LD+coculture_darkness | mixed | low
**Name:** Alteromonas macleodii MIT1002 cluster coculture_LD+coculture_darkness (mixed, periodicity)
**Enrichment:** Not discussed in paper. (p=None, sig=False)
**Functional:** This cluster contains 1 gene (sample MIT1002_2235) that exhibits 24-h periodicity in coculture under both diel L:D and extended darkness conditions. Specific functional enrichment or gene roles are not discussed in the paper.
**Behavioral:** The gene in this cluster maintains 24-h periodic oscillations in coculture both during the diel L:D cycle and under extended darkness, indicating sustained rhythmic expression in the presence of Alteromonas even without light.
**Notes:** Single gene cluster with no detailed functional annotation provided.
**Sources:** Figure 4A, Table S4A

### Cluster coculture_darkness | mixed | low
**Name:** Alteromonas macleodii MIT1002 cluster coculture_darkness (mixed, periodicity)
**Enrichment:** Not discussed in paper. (p=None, sig=False)
**Functional:** This cluster includes 1 gene (sample MIT1002_2228) showing 24-h periodicity only in coculture under extended darkness. Functional roles or enrichment are not described in the paper.
**Behavioral:** The gene exhibits 24-h periodic oscillations exclusively in coculture during extended darkness, not under diel L:D or axenic conditions, suggesting heterotroph influence on maintaining periodicity in the dark.
**Notes:** Minimal functional information available; periodicity pattern noted.
**Sources:** Figure 4A, Table S4A

## Biller 2018 / natl2a_darkness_survival

### Cluster darkness_axenic+darkness_coculture | mixed | medium
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+darkness_coculture (mixed, extended darkness)
**Enrichment:**  (p=None, sig=False)
**Functional:** This cluster includes 96 genes with transcripts present in both axenic and coculture conditions during extended darkness. Functional enrichment is not explicitly detailed for this cluster in the paper, but it includes genes such as PMN2A_RS00075, PMN2A_RS00090, and PMN2A_RS00180 as examples. The cluster likely represents genes maintaining expression across both conditions during extended darkness.
**Behavioral:** Transcripts are present in both axenic and coculture conditions during extended darkness (72-144 h), indicating sustained expression under prolonged dark stress in both culture types.
**Notes:** Functional enrichment and specific gene functions for this cluster are not explicitly described in the paper; descriptions are based on transcript presence patterns.
**Sources:** Figure 1A, Figure 2, Figure 3

### Cluster darkness_axenic+unique_axenic | down | high
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+unique_axenic (down, extended darkness)
**Enrichment:** Biosynthesis, DNA repair, fermentation, stress response (p=0.02, sig=True)
**Functional:** This cluster of 76 genes is unique to axenic cultures during extended darkness and is enriched for transcripts involved in biosynthetic pathways, NAD metabolism, ATP synthase subunits, and DNA repair genes such as recA, ruvB, mutS, and radA. It also shows enrichment of fermentation-associated transcripts and stress response genes including groES, clpB, and dnaJ, as well as activation of the stringent response pathway.
**Behavioral:** Transcripts in this cluster decrease or are depleted during extended darkness in axenic cultures, reflecting cessation of growth and metabolic downregulation. DNA repair and stress response genes are enriched early in extended darkness, indicating stress. This cluster is not enriched or is depleted in cocultures, indicating a unique axenic response.
**Notes:** Strong evidence from transcriptome data and pathway enrichment analyses supports these functional assignments.
**Sources:** Figure 3, Table 1, Table 2

### Cluster darkness_coculture+unique_coculture | up | high
**Name:** Prochlorococcus NATL2A cluster darkness_coculture+unique_coculture (up, extended darkness)
**Enrichment:** Organic compound degradation and salvage, terpenoid and tetrapyrrole biosynthesis (p=0.006, sig=True)
**Functional:** This cluster of 87 genes is unique to cocultured Prochlorococcus during extended darkness and is enriched for transcripts involved in organic compound degradation and salvage pathways, including amino acid degradation, adenine and adenosine salvage (apt), lipoate salvage (lplA), purine nucleotide degradation (truB), terpenoid biosynthesis (dxs, ispG, ispE, lytB), and tetrapyrrole biosynthesis (hemB, hemD). It also includes genes related to the Entner-Doudoroff pathway (eda) suggesting mixotrophic metabolism.
**Behavioral:** Transcripts in this cluster are enriched or maintained during extended darkness in cocultured Prochlorococcus but are depleted or not enriched in axenic cultures, indicating sustained metabolic activity and utilization of heterotroph-derived organic compounds under dark stress.
**Notes:** Strong transcriptomic evidence and pathway enrichment support the functional description of mixotrophic metabolism in cocultures.
**Sources:** Figure 3, Table 2

## Biller 2018 / natl2a_periodicity

### Cluster axenic_LD | mixed | medium
**Name:** Prochlorococcus NATL2A cluster axenic_LD (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** This cluster includes 68 genes showing 24-h periodicity only in axenic L:D conditions. Not discussed in paper specifically for this cluster.
**Behavioral:** Genes exhibit 24-h periodic oscillations in axenic L:D cultures but not in coculture or extended darkness conditions.
**Notes:** Functional description not specifically discussed for this cluster; behavior inferred from periodicity data.
**Sources:** Figure 4A, Table S4A

### Cluster axenic_LD+axenic_darkness+coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster axenic_LD+axenic_darkness+coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster axenic_LD+coculture_LD | mixed | high
**Name:** Prochlorococcus NATL2A cluster axenic_LD+coculture_LD (mixed, periodicity)
**Enrichment:** Photosynthesis and central metabolism (p=0.001, sig=True)
**Functional:** This large cluster of 902 genes includes many transcripts showing 24-h periodicity in both axenic and cocultured L:D conditions but not under extended darkness. Includes genes involved in photosynthesis (psbA, psbB, psaA, psaB), Calvin cycle, glycolysis, fatty acid biosynthesis, and glycogen metabolism.
**Behavioral:** Genes exhibit 24-h periodic oscillations in both axenic and cocultured L:D conditions but lose periodicity under extended darkness.
**Notes:** Strong evidence from periodicity and transcript abundance data.
**Sources:** Figure 4A, Table S4A

### Cluster axenic_LD+coculture_LD+axenic_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster axenic_LD+coculture_LD+axenic_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster axenic_LD+coculture_LD+coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster axenic_LD+coculture_LD+coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster axenic_LD+coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster axenic_LD+coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster coculture_LD | mixed | medium
**Name:** Prochlorococcus NATL2A cluster coculture_LD (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** This cluster of 328 genes shows 24-h periodicity only in coculture L:D conditions. Not specifically discussed in paper for this cluster.
**Behavioral:** Genes exhibit 24-h periodic oscillations in coculture L:D but not in axenic or extended darkness conditions.
**Notes:** Functional description not specifically discussed for this cluster; behavior inferred from periodicity data.
**Sources:** Figure 4A, Table S4A

### Cluster coculture_LD+axenic_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster coculture_LD+axenic_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster coculture_LD+coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster coculture_LD+coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster coculture_darkness | mixed | low
**Name:** Prochlorococcus NATL2A cluster coculture_darkness (mixed, periodicity)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific information available for this cluster.

### Cluster not_periodic | mixed | low
**Name:** Prochlorococcus NATL2A cluster not_periodic (mixed, not periodic)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Genes do not exhibit 24-h periodicity under any tested conditions.
**Notes:** No specific information available for this cluster.

## Coe 2024 / supp_table_3_darktolerant_clusters

### Cluster 1 | up | low
**Name:** Prochlorococcus cluster 1 (up, 183 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 10 | up | low
**Name:** Prochlorococcus cluster 10 (up, 173 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 11 | up | low
**Name:** Prochlorococcus cluster 11 (up, 93 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 12 | up | low
**Name:** Prochlorococcus cluster 12 (up, 94 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 13 | up | low
**Name:** Prochlorococcus cluster 13 (up, 140 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 14 | up | low
**Name:** Prochlorococcus cluster 14 (up, 169 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 15 | up | low
**Name:** Prochlorococcus cluster 15 (up, 130 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 2 | up | low
**Name:** Prochlorococcus cluster 2 (up, 133 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 3 | up | low
**Name:** Prochlorococcus cluster 3 (up, 130 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 4 | up | low
**Name:** Prochlorococcus cluster 4 (up, 169 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 5 | up | low
**Name:** Prochlorococcus cluster 5 (up, 171 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 6 | up | low
**Name:** Prochlorococcus cluster 6 (up, 181 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 7 | up | low
**Name:** Prochlorococcus cluster 7 (up, 118 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 8 | up | low
**Name:** Prochlorococcus cluster 8 (up, 83 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 9 | up | low
**Name:** Prochlorococcus cluster 9 (up, 114 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

## Coe 2024 / supp_table_3_parental_clusters

### Cluster 1 | mixed | low
**Name:** Prochlorococcus cluster 1 (mixed, Yfr106, Yfr2_2, Yfr2_4)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 10 | down | low
**Name:** Prochlorococcus cluster 10 (down, AtpT, cds-PMN2A_RS01135, cds-WP_011293591.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 11 | mixed | low
**Name:** Prochlorococcus cluster 11 (mixed, Yfr1, Yfr2_3, cds-WP_009788866.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 12 | mixed | low
**Name:** Prochlorococcus cluster 12 (mixed, cds-PMN2A_RS05260, cds-PMN2A_RS09840, cds-WP_011293542.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 13 | mixed | low
**Name:** Prochlorococcus cluster 13 (mixed, cds-WP_011293552.1, cds-WP_011293553.1, cds-WP_011293557.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 14 | mixed | low
**Name:** Prochlorococcus cluster 14 (mixed, cds-WP_011293597.1, cds-WP_011293636.1, cds-WP_011293647.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 15 | mixed | low
**Name:** Prochlorococcus cluster 15 (mixed, Yfr10.like_1, Yfr22_1, cds-WP_011125327.1-2)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 2 | mixed | low
**Name:** Prochlorococcus cluster 2 (mixed, Yfr10.like_4, cds-WP_011293588.1, cds-WP_011293608.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 3 | mixed | low
**Name:** Prochlorococcus cluster 3 (mixed, cds-WP_011293551.1, cds-WP_011293554.1, cds-WP_011293576.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 4 | mixed | low
**Name:** Prochlorococcus cluster 4 (mixed, cds-WP_011293580.1, cds-WP_011293592.1, cds-WP_011293599.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 5 | mixed | low
**Name:** Prochlorococcus cluster 5 (mixed, cds-WP_011293563.1, cds-WP_011293645.1, cds-WP_011293646.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 6 | mixed | low
**Name:** Prochlorococcus cluster 6 (mixed, Yfr10.like_2, Yfr103_1, Yfr103_2)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 7 | mixed | low
**Name:** Prochlorococcus cluster 7 (mixed, Yfr108, cds-WP_011293543.1, cds-WP_011293544.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 8 | mixed | low
**Name:** Prochlorococcus cluster 8 (mixed, Yfr107, cds-WP_011293537.1, cds-WP_011293538.1)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 9 | mixed | low
**Name:** Prochlorococcus cluster 9 (mixed, Yfr13, Yfr22_2, Yfr22_3)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

## Lindell 2007 / med4_phage_transcription_groups

### Cluster 1 | up | high
**Name:** Prochlorococcus MED4 cluster 20 genes (up, phage infection response)
**Enrichment:** Host stress response and metabolism (p=None, sig=False)
**Functional:** Includes 20 genes transiently upregulated immediately after infection, such as high-light-inducible stress response genes (hli), carbon metabolism gene rbcLS, transcription genes rpoC2 and rpoD, and ribosomal protein genes rpl5, rpl6, rps8, rps11, and rps17. These genes are associated with host stress response and metabolism.
**Behavioral:** Genes show a rapid transient peak in expression immediately after phage infection, indicating an early stress response by the host.
**Notes:** Gene identities and timing are well supported by microarray and RT-PCR data.
**Sources:** Figure 3, Supplementary Table 3

### Cluster 2 | up | high
**Name:** Prochlorococcus MED4 cluster 25 genes (up, phage infection response)
**Enrichment:** RNA degradation, protein turnover, and stress response (p=None, sig=False)
**Functional:** Contains 25 genes upregulated starting 2 hours post infection, including RNA degradation and modification genes (rne, rnhB, dus, sun), protein turnover genes (clpS, AAA ATPase family), stress response genes (umuD, phoH), and genes of unknown function. Many are located in hypervariable genome islands and have homologues in phage genomes.
**Behavioral:** Genes show a delayed but sustained upregulation beginning 2 hours after infection, suggesting involvement in RNA processing, protein turnover, and stress response during phage infection.
**Notes:** Functional categories inferred from gene annotations; timing supported by expression data.
**Sources:** Figure 3, Supplementary Table 3

## Tolonen 2006 / med4_kmeans_nstarvation

### Cluster 1 | up | high
**Name:** Prochlorococcus cluster 1 (up, transport and binding)
**Enrichment:** transport and binding (p=0.01, sig=True)
**Functional:** Contains nitrogen transport genes such as urtA, cynA, and the highly upregulated unknown gene PMM0958 with a top-ranking NtcA binding site, indicating a role in nitrogen transport and assimilation during N starvation.
**Behavioral:** Most rapidly and highly upregulated cluster, with genes responding within 6 hours of N deprivation and peaking early in the time course.
**Notes:** Strong functional enrichment and early rapid upregulation.
**Sources:** Figure 3A, Table I, Table III

### Cluster 2 | up | medium
**Name:** Prochlorococcus cluster 2 (up, transport and binding)
**Enrichment:** transport and binding (p=0.12, sig=False)
**Functional:** Contains a subset of hli genes (hli10, hli21, hli22) and genes with putative NtcA binding sites, suggesting involvement in nitrogen stress response and photoprotection.
**Behavioral:** Rapid and high upregulation early in the N starvation time course, distinct from a later responding hli gene subset.
**Notes:** Functional enrichment not statistically significant but contains notable genes.
**Sources:** Figure 3A, Figure 4E

### Cluster 3 | up | medium
**Name:** Prochlorococcus cluster 3 (up, regulation)
**Enrichment:** regulation (p=0.07, sig=False)
**Functional:** Contains a lesser degree of upregulated hli genes and regulatory genes including sigma factors, indicating a role in transcriptional regulation during N stress.
**Behavioral:** Genes respond later and to a lesser degree compared to cluster 2, showing a more gradual upregulation pattern.
**Notes:** Functional enrichment marginal; includes regulatory genes and sigma factors.
**Sources:** Figure 3A, Figure 4B

### Cluster 4 | down | low
**Name:** Prochlorococcus cluster 4 (down, amino acid synthesis)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** Paper does not discuss this cluster.

### Cluster 5 | down | low
**Name:** Prochlorococcus cluster 5 (down, regulation)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** Paper does not discuss this cluster.

### Cluster 6 | down | high
**Name:** Prochlorococcus cluster 6 (down, translation)
**Enrichment:** translation (p=0.001, sig=True)
**Functional:** Strongly enriched for translation-related genes (p=0.001), including ribosomal proteins, indicating repression of protein synthesis during N starvation.
**Behavioral:** Highly downregulated cluster with rapid onset of repression during N starvation.
**Notes:** Strong functional enrichment and clear repression pattern.
**Sources:** Figure 3A

### Cluster 7 | down | high
**Name:** Prochlorococcus cluster 7 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=2.6e-09, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes (p=2.6e-9), including photosystem I genes (psaBDEIJKLM), indicating repression of photosynthetic machinery during N starvation.
**Behavioral:** Strongly downregulated cluster with rapid repression during N starvation.
**Notes:** Strong functional enrichment and clear repression pattern.
**Sources:** Figure 3A

### Cluster 8 | down | high
**Name:** Prochlorococcus cluster 8 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=1.4e-05, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes (p=1.4e-5), including ATP synthase subunits and carbon fixation genes (rbcLS), indicating repression of energy metabolism during N starvation.
**Behavioral:** Strongly downregulated cluster with rapid repression during N starvation.
**Notes:** Strong functional enrichment and clear repression pattern.
**Sources:** Figure 3A

## Tolonen 2006 / mit9313_kmeans_nstarvation

### Cluster 1 | up | high
**Name:** MIT9313 cluster 1 (up, transport and binding)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Enriched for transport and binding (p=0.04). Contains nitrogen transport genes urtA and the nitrite permease, and hli genes hliS and hli7. Also includes genes with high-ranking NtcA binding sites such as ntcA and sigma factors.
**Behavioral:** Most rapidly and strongly upregulated cluster during nitrogen starvation, with genes responding within the first hours of the time course and sustained upregulation.
**Notes:** Strong functional enrichment and early rapid upregulation.
**Sources:** Figure 3, Table I, Table II

### Cluster 2 | up | high
**Name:** MIT9313 cluster 2 (up, regulation)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** Contains regulatory genes including two sigma factors (PMT2246 and PMT0346) that are upregulated during nitrogen starvation. One sigma factor (PMT2246) has a strong NtcA binding site, suggesting indirect regulation by NtcA.
**Behavioral:** Upregulated during nitrogen starvation, likely contributing to the regulation of gene expression in response to nitrogen stress.
**Notes:** Sigma factors upregulated with evidence of NtcA regulation for one factor.
**Sources:** Figure 4B, Table I

### Cluster 3 | up | medium
**Name:** MIT9313 cluster 3 (up, regulation)
**Enrichment:** regulation (p=0.09, sig=False)
**Functional:** Contains regulatory genes with weaker enrichment (p=0.09). Specific genes not detailed in paper.
**Behavioral:** Upregulated during nitrogen starvation with a slower and sustained pattern compared to cluster 1.
**Notes:** Functional enrichment is weak and specific genes not discussed.
**Sources:** Figure 3

### Cluster 4 | down | high
**Name:** MIT9313 cluster 4 (down, transport and binding)
**Enrichment:** transport and binding (p=0.39, sig=False)
**Functional:** Contains genes linking nitrogen and carbon metabolism including glnB, icd, acnB, and rbcLS. These genes are repressed only at the final time point, likely reflecting general transcriptional shutdown rather than specific nitrogen stress response.
**Behavioral:** Genes unchanged until repressed at 48h, coinciding with severe starvation and general transcriptional shutdown.
**Notes:** Repression likely due to general shutdown, not specific N stress response.
**Sources:** Figure 3, Figure 6

### Cluster 5 | down | medium
**Name:** MIT9313 cluster 5 (down, fatty acid and phospholipid metabolism)
**Enrichment:** fatty acid and phospholipid metabolism (p=0.004, sig=True)
**Functional:** Enriched for fatty acid and phospholipid metabolism (p=0.004). Specific genes not detailed in paper.
**Behavioral:** Downregulated during nitrogen starvation with increasing repression over time.
**Notes:** Strong functional enrichment but limited gene details.
**Sources:** Figure 3

### Cluster 6 | down | high
**Name:** MIT9313 cluster 6 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=0.00049, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes (p=4.9e-4). Contains genes for photosystem I and II and phycoerythrin gene cpeB.
**Behavioral:** Genes gradually repressed during nitrogen starvation, with repression increasing over time and peaking at the final time point.
**Notes:** Strong enrichment and consistent repression pattern.
**Sources:** Figure 3

### Cluster 7 | down | high
**Name:** MIT9313 cluster 7 (down, translation)
**Enrichment:** translation (p=7.9e-10, sig=True)
**Functional:** Enriched for translation genes (p=7.9e-10). Contains ribosomal protein genes and other translation-related genes.
**Behavioral:** Strongly downregulated during nitrogen starvation with sustained repression.
**Notes:** Strong functional enrichment and clear repression pattern.
**Sources:** Figure 3

## Wang 2014 / med4_expression_level

### Cluster -- | mixed | low
**Name:** Prochlorococcus cluster -- (mixed, expression level)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific functional or behavioral description provided for this cluster in the paper.

### Cluster HEG | up | high
**Name:** Prochlorococcus cluster HEG (up, expression level)
**Enrichment:** Energy production and conversion, Translation and ribosomal structure, Protein modification and folding (p=0.04, sig=True)
**Functional:** Highly expressed genes (HEG) are overrepresented in the core genome and enriched for central metabolic functions including energy production and conversion (COG C), translation and ribosomal structure (COG J), and protein modification, folding and turnover (COG O). Genes such as ribosomal proteins, photosynthetic apparatus components, and protein folding genes are included.
**Behavioral:** Genes in this cluster are consistently highly expressed across all growth conditions, indicating stable and abundant expression.
**Notes:** Strong functional enrichment and consistent expression across conditions.
**Sources:** Figure 4c, Figure 5b

### Cluster LEG | down | low
**Name:** Prochlorococcus cluster LEG (down, expression level)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific functional or behavioral description provided for this cluster in the paper.

### Cluster MEG | mixed | medium
**Name:** Prochlorococcus cluster MEG (mixed, expression level)
**Enrichment:** Essential genes (p=0.001, sig=True)
**Functional:** Moderately expressed genes (MEG) are enriched in the core genome and have a significantly higher rate of essential genes (DEG-hit). They include genes involved in central metabolism but with less enrichment than HEG.
**Behavioral:** Genes in this cluster are moderately and consistently expressed across all growth conditions.
**Notes:** Moderate functional enrichment and consistent expression.
**Sources:** Figure 3c, Figure 4b

### Cluster NEG | down | low
**Name:** Prochlorococcus cluster NEG (down, expression level)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.
**Notes:** No specific functional or behavioral description provided for this cluster in the paper.

### Cluster VEG | up | medium
**Name:** Prochlorococcus cluster VEG (up, expression level)
**Enrichment:** Variable expression (p=None, sig=False)
**Functional:** Variably expressed genes (VEG) are less conserved, more frequent in the flexible genome, and have higher nonsynonymous substitution rates. They include genes with variable expression across conditions and are less likely to be essential.
**Behavioral:** Genes in this cluster show variable expression levels across different growth conditions rather than consistent expression.
**Notes:** Described as variably expressed with less conservation and higher evolutionary rates.
**Sources:** Figure 3c

## Zinser 2009 / med4_diel_clusters

### Cluster 1 | up | high
**Name:** Prochlorococcus cluster 1 (up, photosystem I and II)
**Enrichment:** Photosystem I and II (p=1.5e-09, sig=True)
**Functional:** Enriched for photosystem I (p=1.5e-9) and photosystem II (p=2.5e-5) genes. Includes reaction center genes psbA, psbD, and psaA involved in photosynthetic light reactions.
**Behavioral:** Genes peak in expression near mid-morning to mid-day (8.3 h), coinciding with high light intensity and photosynthetic activity, with 24h periodicity.
**Sources:** Figure 4A, Table 1, Table S4

### Cluster 10 | up | low
**Name:** Prochlorococcus cluster 10 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 11 | up | low
**Name:** Prochlorococcus cluster 11 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 12 | up | low
**Name:** Prochlorococcus cluster 12 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 13 | up | high
**Name:** Prochlorococcus cluster 13 (up, ribosomal proteins)
**Enrichment:** Ribosomal proteins (p=3.7e-41, sig=True)
**Functional:** Strongly enriched for ribosomal proteins (p=3.7e-41), including 45 of 53 ribosomal protein genes.
**Behavioral:** Genes peak near dawn (3.4 h), consistent with preparation for protein synthesis early in the day, with 24h periodicity.
**Sources:** Table 1

### Cluster 14 | up | low
**Name:** Prochlorococcus cluster 14 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 15 | up | high
**Name:** Prochlorococcus cluster 15 (up, Calvin cycle and ATP synthase)
**Enrichment:** ATP synthase and CO2 metabolism (p=8.4e-08, sig=True)
**Functional:** Enriched for ATP synthase (p=8.4e-8) and CO2 metabolism (p=1.8e-5) genes, including Calvin cycle genes rbcL, csoS3, fbaA, and ATPase subunits.
**Behavioral:** Genes peak near dawn (4.7 h), consistent with initiation of carbon fixation and energy production at the start of the light period.
**Sources:** Table 1, Figure 8C

### Cluster 16 | up | high
**Name:** Prochlorococcus cluster 16 (up, ATP synthase and CO2 metabolism)
**Enrichment:** ATP synthase and CO2 metabolism (p=8.4e-08, sig=True)
**Functional:** Enriched for ATP synthase (p=8.4e-8) and CO2 metabolism (p=1.8e-5) genes, including ATPase subunits and Calvin cycle genes.
**Behavioral:** Genes peak near dawn (5.5 h), supporting energy production and carbon fixation early in the light cycle.
**Sources:** Table 1

### Cluster 17 | down | low
**Name:** Prochlorococcus cluster 17 (down, aperiodic expressed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 18 | down | medium
**Name:** Prochlorococcus cluster 18 (down, menaquinone and ubiquinone)
**Enrichment:** Menaquinone and ubiquinone biosynthesis (p=0.00045, sig=True)
**Functional:** Enriched for menaquinone and ubiquinone biosynthesis genes (p=0.00045).
**Behavioral:** Not discussed in paper.
**Notes:** Behavioral description not provided in paper.
**Sources:** Table 1

### Cluster 2 | up | low
**Name:** Prochlorococcus cluster 2 (up, not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Not discussed in paper.
**Behavioral:** Not discussed in paper.

### Cluster 3 | up | high
**Name:** Prochlorococcus cluster 3 (up, cytochrome b6/f and photosystem II)
**Enrichment:** Cytochrome b6/f and Photosystem II (p=3.7e-08, sig=True)
**Functional:** Enriched for cytochrome b6/f (p=0.0059) and photosystem II (p=3.7e-8) genes, including petA, petB, petD, and psb genes.
**Behavioral:** Genes peak near mid-day (12.5 h), coinciding with peak photosynthetic electron transport activity.
**Sources:** Table 1, Figure 4C

### Cluster 4 | up | medium
**Name:** Prochlorococcus cluster 4 (up, cytochrome b6/f)
**Enrichment:** Cytochrome b6/f (p=0.073, sig=False)
**Functional:** Enriched for cytochrome b6/f genes (p=0.073).
**Behavioral:** Genes peak in mid-afternoon (15.8 h), slightly later than cluster 3, possibly reflecting continued photosynthetic electron transport activity.
**Notes:** Marginal enrichment significance.
**Sources:** Table 1

### Cluster 5 | up | high
**Name:** Prochlorococcus cluster 5 (up, respiratory terminal oxidases and protein degradation)
**Enrichment:** Respiratory terminal oxidases and protein degradation (p=0.019, sig=True)
**Functional:** Enriched for respiratory terminal oxidases (p=0.019) and degradation of proteins, peptides, and glycopeptides (p=0.071).
**Behavioral:** Genes peak near early night (17.5 h), consistent with increased respiration and protein turnover at night.
**Notes:** Protein degradation enrichment is marginal.
**Sources:** Table 1

### Cluster 6 | up | high
**Name:** Prochlorococcus cluster 6 (up, purine ribonucleotide biosynthesis)
**Enrichment:** Purine ribonucleotide biosynthesis (p=0.0049, sig=True)
**Functional:** Enriched for purine ribonucleotide biosynthesis genes (p=0.0049).
**Behavioral:** Genes peak near early night (18.6 h), consistent with nucleotide synthesis during DNA replication phase.
**Sources:** Table 1

### Cluster 7 | up | medium
**Name:** Prochlorococcus cluster 7 (up, nitrogen metabolism)
**Enrichment:** Nitrogen metabolism (p=0.087, sig=False)
**Functional:** Includes nitrogen metabolism genes such as amt1, glnA, glsF, and icd, with marginal enrichment (p=0.087).
**Behavioral:** Genes peak near sunset (20.1 h), consistent with nitrogen uptake and assimilation during the night.
**Notes:** Marginal enrichment; limited description in paper.
**Sources:** Figure 7B

### Cluster 8 | up | high
**Name:** Prochlorococcus cluster 8 (up, chaperones)
**Enrichment:** Chaperones (p=0.00023, sig=True)
**Functional:** Enriched for chaperone genes (p=0.00023).
**Behavioral:** Genes peak near late night (21.0 h), possibly reflecting protein folding and stress response during night.
**Sources:** Table 1

### Cluster 9 | up | high
**Name:** Prochlorococcus cluster 9 (up, RNA synthesis, modification, and DNA transcription)
**Enrichment:** RNA synthesis, modification, and DNA transcription (p=0.035, sig=True)
**Functional:** Enriched for RNA synthesis, modification, and DNA transcription genes (p=0.035).
**Behavioral:** Genes peak near late night (22.4 h), consistent with preparation for transcription during night.
**Sources:** Table 1

## Warnings

- [Biller 2018 / natl2a_periodicity / cluster not_periodic] low confidence but behavioral_description is not 'Not discussed in paper.': Genes do not exhibit 24-h periodicity under any tested condi...
- [Biller 2018 / mit1002_periodicity / cluster coculture_LD] locus tag in functional_description: MIT1002_0002
- [Biller 2018 / mit1002_periodicity / cluster coculture_LD+coculture_darkness] locus tag in functional_description: MIT1002_2235
- [Biller 2018 / mit1002_periodicity / cluster coculture_darkness] locus tag in functional_description: MIT1002_2228
- [Biller 2018 / mit1002_periodicity / cluster coculture_LD+coculture_darkness] low confidence but functional_description is not 'Not discussed in paper.': This cluster contains 1 gene (sample MIT1002_2235) that exhi...
- [Biller 2018 / mit1002_periodicity / cluster coculture_LD+coculture_darkness] low confidence but behavioral_description is not 'Not discussed in paper.': The gene in this cluster maintains 24-h periodic oscillation...
- [Biller 2018 / mit1002_periodicity / cluster coculture_darkness] low confidence but functional_description is not 'Not discussed in paper.': This cluster includes 1 gene (sample MIT1002_2228) showing 2...
- [Biller 2018 / mit1002_periodicity / cluster coculture_darkness] low confidence but behavioral_description is not 'Not discussed in paper.': The gene exhibits 24-h periodic oscillations exclusively in ...
- [Biller 2018 / natl2a_darkness_survival / cluster darkness_axenic+darkness_coculture] locus tag in functional_description: PMN2A_RS00075
- [Biller 2018 / natl2a_darkness_survival / cluster darkness_axenic+darkness_coculture] filler phrase 'likely' in functional_description
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 1] locus tag in functional_description: PMM0958
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 8] near-identical functional_description as cluster 7
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 2] locus tag in functional_description: PMT2246
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 2] filler phrase 'likely' in behavioral_description
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 4] filler phrase 'likely' in functional_description
- [Zinser 2009 / med4_diel_clusters / cluster 4] filler phrase 'possibly' in behavioral_description
- [Zinser 2009 / med4_diel_clusters / cluster 8] filler phrase 'possibly' in behavioral_description
- [Zinser 2009 / med4_diel_clusters / cluster 16] near-identical functional_description as cluster 15
- [Wang 2014 / med4_expression_level / cluster VEG] filler phrase 'likely' in functional_description
- [Bernstein 2017 / mruber_light_clusters / cluster 0] low confidence but functional_description is not 'Not discussed in paper.': Cluster 0 contains 857 genes including SY28_RS00050, SY28_RS...
- [Bernstein 2017 / mruber_light_clusters / cluster 4] low confidence but functional_description is not 'Not discussed in paper.': Cluster 4 contains 444 genes including SY28_RS05005, SY28_RS...
- [Bernstein 2017 / mruber_light_clusters / cluster 2] filler phrase 'likely' in functional_description
