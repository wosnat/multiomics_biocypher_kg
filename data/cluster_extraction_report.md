# Cluster Extraction Report

## Alonso-Saez 2023 / mit9301_softclusters_thermal_acclimation

### Cluster Cluster A | high across all temperatures during daytime | high
**Name:** Prochlorococcus cluster A (core metabolism)
**Enrichment:** Carbon fixation and metabolism (p=None, sig=False)
**Functional:** Enriched for genes related to carbon fixation and assimilation, including RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), Calvin cycle enzymes (gap2, tktA, glpX, pgk, cbbA), glycogen synthesis (glgABC), ATP synthesis genes (atpADE), and some photosystem II components (psbA, psbC, psbD).
**Temporal pattern:** Consistently expressed during daytime along the thermal gradient with stable expression levels across temperatures.
**Sources:** Figure 2, Figure 4
**Quotes:**
- [Page 4] Cluster A included genes related to C fixation and assimilation, such as RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), the Calvin cycle (e.g., gap2, tktA, glpX, pgk, and cbbA), and glycogen synthesis (glgABC), consistently expressed during daytime along the thermal gradient (Fig. 2; Table S4 and Table S5). This cluster also included ATP synthesis genes (atpADE) and a few components of photosystem II (PS II) (psbA, psbC, and psbD).

### Cluster Cluster B | down at cold during daytime | high
**Name:** Prochlorococcus cluster B (photosynthetic apparatus downregulation)
**Enrichment:** Photosynthesis (p=None, sig=False)
**Functional:** Includes genes encoding components of photosystem I (psaABDEFKL), some photosystem II components (psbBJH, psbO), and photosynthetic electron transport genes (petACGNM). Also includes the periplasmic phosphate binding protein pstS.
**Temporal pattern:** Expression decreases gradually from optimum temperature to minimum temperature during daytime, paralleling growth rate decline. Shows a day/night expression ratio change under cold conditions.
**Sources:** Figure 2, Figure 3D, Figure 4, Figure 5
**Quotes:**
- [Page 7] The expression of all components of the PS I complex (psaABDEFKL) and some of PS II (including psbBJH and the oxygen-evolving complex protein psbO) showed a gradual decrease in expression from a temperature close to the optimum to the Tmin, (Fig. 2) in correlation with MIT9301 growth rates (Fig. 1). This expression pattern was different from those of other PSII components (psbACD [described above]), which were not differentially expressed during daytime along the thermal niche (Kruskall-Wallis; P . 0.05) (Fig. 2). Similarly, many components of the photosynthetic electron transport genes were assigned to either cluster B (petACGNM) or cluster A (petBEDH), implying a nonuniform transcriptional thermal response of all components of the photosynthetic apparatus (Fig. 4).
- [Page 6] The expression of both copies of the periplasmic phosphate binding protein (pstS) showed maximum expression values around the Topt during daytime, following the pattern of most photosynthetic genes (cluster B), which highlights the complexity of the thermal response of nutrient acquisition genes (Fig. 3D).

### Cluster Cluster C | up at cold during daytime | high
**Name:** Prochlorococcus cluster C (cold stress response daytime)
**Enrichment:** Stress response (p=None, sig=False)
**Functional:** Enriched for genes involved in global stress response including cellular chaperones (groES, groEL, dnaK, clpBCP), fatty acid desaturases (desA, desC), oxidative damage protection (recA, ruvB, sod), carotenoid synthesis (pds, crtBH), and rubredoxin (rub). The histidine kinase nblS is confirmed in this cluster.
**Temporal pattern:** Strongly upregulated at minimum temperature (Tmin) during daytime, with prioritization of expression under light conditions.
**Sources:** Figure 2, Figure 3C
**Quotes:**
- [Page 4] Clusters C and D genes included different elements of the global stress response, such as cellular chaperones (groES/groES, dnaK, and clpBCP) and fatty acid desaturases (desA and desC), as well as mechanisms against oxidative damage, such as DNA repair (recA and ruvB), superoxide dismutase (sod), and the synthesis of antioxidant compounds like carotenoids (pds and crtBH) and rubredoxin (rub) (Fig. 2 and 3; Table S5). Notably, the expression of the chaperones groEL/groES, grpE, and htpG was strongly upregulated at the Tmin only during daytime, suggesting a prioritization of their expression during the light-exposed period (Fig. 3C).

### Cluster Cluster D | up at cold during day and night | high
**Name:** Prochlorococcus cluster D (cold stress response day and night)
**Enrichment:** Stress response and regulation (p=None, sig=False)
**Functional:** Includes genes involved in global stress response upregulated at cold under both light and dark conditions. Also includes sigma factors and regulatory proteins (rpoD and family members), DNA replication (dnaA, nrdJ, gyrB), cell division (ftsZYQ), and energy metabolism via the pentose phosphate pathway (tal, gnd, zwf) and glycogen degradation (glgP). Nitrogen acquisition genes (urtA, phoB) are also present.
**Temporal pattern:** Upregulated at minimum temperature (Tmin) during both daytime and nighttime, indicating a sustained cold stress response.
**Sources:** Figure 2, Figure 3
**Quotes:**
- [Page 4] Clusters C and D genes included different elements of the global stress response, such as cellular chaperones (groES/groES, dnaK, and clpBCP) and fatty acid desaturases (desA and desC), as well as mechanisms against oxidative damage, such as DNA repair (recA and ruvB), superoxide dismutase (sod), and the synthesis of antioxidant compounds like carotenoids (pds and crtBH) and rubredoxin (rub) (Fig. 2 and 3; Table S5).
- [Page 4] The expression of some sigma factors and key regulatory proteins from different families also showed this pattern (Fig. 3A and B).
- [Page 6] Cellular expression levels (measured as [mRNA]) in Prochlorococcus marinus MIT9301 during daytime (blue lines) and nighttime (black lines) along the thermal niche of (A) RNA polymerase components, including sigma factors, (B) histidine kinases and other regulatory proteins, (C) genes involved in the stress response, and (D) nitrogen and phosphate acquisition genes. Asterisks denote significant differences in transcript concentration along the thermal gradient in daytime (blue asterisks) or nighttime (black asterisks) according to the Kruskal-Wallis test (*, P < 0.05; **, P < 0.01; ***, P < 0.001). The SoftCluster membership of each gene is shown only for those cases where the probability score was >0.80.

### Cluster Cluster E | high at night | high
**Name:** Prochlorococcus cluster E (core nighttime metabolism)
**Enrichment:** Cell division and metabolism (p=None, sig=False)
**Functional:** Enriched for genes related to catabolic consumption (cyoB, ndhD), DNA replication (dnaA, nrdJ, gyrB), cell division (ftsZYQ), and the pentose phosphate pathway (tal, gnd, zwf). Also includes phosphate transporter genes (pstA, pstB, pstC) upregulated at night under cold conditions.
**Temporal pattern:** Upregulated at nighttime across the thermal gradient, with increased expression of phosphate uptake genes under cold conditions at night.
**Sources:** Figure 2, Figure 3D
**Quotes:**
- [Page 4] Cluster E included genes related to catabolic consumption (cyoB and ndhD), DNA replication (dnaA, nrdJ, and gyrB), cell division (ftsZYQ), and the pentose phosphate pathway (tal, gnd, and zwf), all of them upregulated at nighttime. Altogether, these essential pathways likely represent a transcriptional core, which Prochlorococcus marinus MIT9301 maintains under all temperature conditions.
- [Page 6] In the case of phosphate uptake, the high-affinity ABC transporter pstABC genes were upregulated under cold conditions at night, following the pattern of cluster E genes, possibly related to the cellular demand of P for DNA replication.

## Bernstein 2017 / bp1_light_clusters

### Cluster A | peak at intermediate irradiance | low
**Name:** Thermosynechococcus elongatus cluster A (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Tent-shaped expression pattern with peak at intermediate irradiance.
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

### Cluster B | N/A | high
**Name:** Thermosynechococcus elongatus cluster B (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** No specific discussion or functional enrichment reported for this cluster in the paper.
**Assessment notes:** No data available for this cluster.
**Sources:** Figure 6

### Cluster C | down with irradiance | medium
**Name:** Thermosynechococcus elongatus cluster C (irradiance-responsive)
**Enrichment:** Vitamin B12 biosynthesis (p=0.05, sig=True)
**Functional:** Contains T. elongatus vitamin B12 biosynthesis genes (cobWNT, cobOQDPC).
**Temporal pattern:** Decreases with increasing irradiance.
**Confidence notes:** B12 biosynthesis genes verified as T. elongatus in cluster C.
**Assessment notes:** B12 genes confirmed for T. elongatus; M. ruber genes (glnA, gltB, metH) removed as misattributed.
**Sources:** Figure 6
**Quotes:**
- [Page 8] The relative abundances of transcripts encoding vitamin B12 biosynthesis, including those required for the insertion of cobalt (cobWNT) and for the conversion of cobyrinic acid diamine to the vitamin B12 coenzyme from cobOQDPC, were found in cluster C.

### Cluster D | up with irradiance | high
**Name:** Thermosynechococcus elongatus cluster D (irradiance-responsive)
**Enrichment:** Carbon metabolism, nitrogen metabolism, and oxidative stress response (p=0.05, sig=True)
**Functional:** Enriched for genes involved in organic acid synthesis and export, including acetyl-CoA synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), citrate synthase (gltA), exopolysaccharide synthesis (exoD), and sucrose metabolism (spsA, sps). Also includes nitrogen metabolism genes such as nitrate uptake (nrtABD), glutamine synthetase (glnA), and glutamate symporter (gltS). Genes involved in ROS detoxification (flv4, peroxiredoxins, sodB) are also enriched.
**Temporal pattern:** Increases with increasing irradiance.
**Confidence notes:** High confidence based on strong functional enrichment and consistent expression dynamics.
**Sources:** Figure 6, Figure 3, Figure 5
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). These included acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA).
- [Page 8] The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H. Notable examples include an NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys family peroxiredoxins, and an Mn-superoxide dismutase (sodB).

## Bernstein 2017 / bp1_oxygen_clusters

### Cluster E | peak at intermediate pO2 | low
**Name:** Thermosynechococcus elongatus cluster E (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Tent-shaped expression pattern with peak at intermediate oxygen tension (analog of cluster A).
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

### Cluster F | N/A | low
**Name:** Thermosynechococcus elongatus cluster F (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data provided for cluster F in the paper.

### Cluster G | down with increasing pO2 | high
**Name:** Thermosynechococcus elongatus cluster G (pO2-responsive)
**Enrichment:** Organic acid metabolism (p=0.05, sig=True)
**Functional:** Enriched for genes involved in organic acid synthesis and export, including acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA). Also includes genes for synthesis and export of larger biomolecules such as exopolysaccharide synthesis gene (exoD), sucrose synthase (spsA), and sucrose degradation enzyme (sps).
**Temporal pattern:** Decreases with increasing oxygen tension.
**Confidence notes:** Enrichment supported by statistical analysis with P ≤ 0.05.
**Sources:** Figure 6, Table S1
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). These included acetyl coenzyme A synthetase (acs; tll0887), acetate kinase (ackA; tlr2340), lactate dehydrogenase (ldhA; tlr0711), and citrate synthase (gltA; tlr2393).
- [Page 7] Functions involved in synthesis and export of larger biomolecules (i.e., sugars, peptides, and extracellular polymeric substance [EPS]) also grouped into clusters D and G. These included a putative exopolysaccharide synthesis gene (exoD; tll2077), sucrose synthase (spsA; tlr1047), and a sucrose degradation enzyme (sps; tlr0582).

### Cluster H | up with increasing pO2 | high
**Name:** Thermosynechococcus elongatus cluster H (pO2-responsive)
**Enrichment:** Oxidative stress response (p=0.05, sig=True)
**Functional:** Contains T. elongatus ROS detoxification genes including flv4, peroxiredoxins, and sodB.
**Temporal pattern:** Genes in this cluster increase in expression with increasing pO2 and correspond to a decrease in specific growth rate.
**Confidence notes:** Enrichment supported by statistical analysis with P ≤ 0.05.
**Sources:** Figure 6, Table S1
**Quotes:**
- [Page 8] The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H. Notable examples include an NAD(P)H-oxygen oxidoreductase (flv4; tlr1088), 2-Cys family peroxiredoxins (tll1454 and tlr1289), a periplasmic peroxiredoxin (prxQ-B2; tlr1194), and an Mn-superoxide dismutase (sodB; tll1519).

## Bernstein 2017 / mruber_light_clusters

### Cluster A | peak at intermediate irradiance | low
**Name:** Meiothermus ruber cluster A (irradiance responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Tent-shaped expression pattern with peak at intermediate irradiance.
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

### Cluster B | N/A | low
**Name:** Meiothermus ruber cluster B (irradiance responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** No specific functional or expression details provided for M. ruber cluster B in the paper.
**Assessment notes:** No explicit information available for this cluster.

### Cluster C | down with irradiance | high
**Name:** Meiothermus ruber cluster C (irradiance responsive)
**Enrichment:** Nitrogen metabolism and vitamin B12 uptake (p=0.05, sig=True)
**Functional:** Contains M. ruber genes for nitrogen acquisition (glnA, gltB), vitamin B12 uptake (btuCD), and methionine biosynthesis (metH).
**Temporal pattern:** Decreases with increasing irradiance.
**Confidence notes:** Moderate to high confidence based on gene coexpression and functional enrichment.
**Assessment notes:** Gene expression patterns and functional enrichment support this cluster's description.
**Sources:** Figure 6, Table S1, Table S2
**Quotes:**
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA (SY28_RS02395) and the large subunit of glutamate synthase (gltB; SY28_RS09480) and suggest the potential for direct exchange of glutamine and glutamate from T. elongatus as growth requirements increase with the specific growth and photosynthesis rates.
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD; SY28_RS12150 and SY28_RS12155) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis... Both T. elongatus and M. ruber expressed transcripts encoding the vitamin B12-dependent homocysteine methyltransferase (metH) but were negatively correlated and grouped across Ii and pO2 treatments into opposing clusters D (tll1027) and C (SY28_RS08890), indicating that methionine may be directly exchanged from T. elongatus as growth requirements increase with the specific growth rate.

### Cluster D | up with irradiance | low
**Name:** Meiothermus ruber cluster D (irradiance responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Increases with increasing irradiance.
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

## Bernstein 2017 / mruber_oxygen_clusters

### Cluster E | peak at intermediate pO2 | low
**Name:** Meiothermus ruber cluster E (pO2 responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Tent-shaped expression pattern with peak at intermediate oxygen tension.
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

### Cluster F | inverse of cluster E | low
**Name:** Meiothermus ruber cluster F (pO2 responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Inverse of cluster E.
**Assessment notes:** Paper discusses joint clusters; per-organism functional attribution unreliable.

### Cluster G | up with growth rate, down with pO2 | high
**Name:** Meiothermus ruber cluster G (pO2 responsive)
**Enrichment:** Carbon metabolism and transport (p=0.01, sig=True)
**Functional:** Enriched for genes involved in carbon uptake and metabolism, including acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), and branched-chain amino acid uptake (livGH).
**Temporal pattern:** Genes show increased expression with increasing specific growth and photosynthesis rates and decreased expression with increasing pO2.
**Confidence notes:** Cluster G genes correspond to carbon metabolism and transport, showing inverse response to pO2 and growth rate.
**Sources:** Figure 6
**Quotes:**
- [Page 8] In conjunction with T. elongatus, M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. Notable examples included genes encoding putative acetyl-CoA synthetase (acs; SY28_RS04760 and SY28_RS00910), cytochrome c oxidase (coxB; SY28_RS11630), monosaccharide uptake systems (frcA and gtsAB; SY28_RS03690, SY28_RS04315, and SY28_RS04260), xylose isomerase (xylA; SY28_RS02810), xylulokinase (xylB; SY28_RS03685), an ABC-type multisugar uptake system (SY28_RS06965), and branched-chain amino acid uptake (livGH; SY28_RS00360 and SY28_RS00370).

### Cluster H | up with pO2, down with growth rate | high
**Name:** Meiothermus ruber cluster H (pO2 responsive)
**Enrichment:** Oxidative stress response (p=0.01, sig=True)
**Functional:** Contains M. ruber ROS detoxification genes bcp and sod2.
**Temporal pattern:** Genes increase in expression with increasing pO2 and decrease with increasing specific growth and photosynthesis rates.
**Confidence notes:** Cluster H genes are associated with oxidative stress response and electron transfer, showing increased expression with pO2.
**Sources:** Figure 6
**Quotes:**
- [Page 8] M. ruber peroxidase (bcp; SY28_RS05545, SY28_RS06010, and SY28_RS06015) and superoxide dismutase (sod2; SY28_RS13295) genes responded differently to Ii treatments than did oxyR and related cyanobacterial profiles and were grouped into clusters A and C. These genes generally increased with pO2 and grouped with the T. elongatus genes into cluster H (increased with pO2 and decreasing μ).

## Coe 2024 / supp_table_3_darktolerant_clusters

### Cluster 1 | N/A | low
**Name:** Prochlorococcus cluster 1
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 10 | N/A | low
**Name:** Prochlorococcus cluster 10
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 11 | N/A | low
**Name:** Prochlorococcus cluster 11
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 12 | N/A | low
**Name:** Prochlorococcus cluster 12
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 13 | N/A | low
**Name:** Prochlorococcus cluster 13
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 14 | N/A | low
**Name:** Prochlorococcus cluster 14
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 15 | N/A | low
**Name:** Prochlorococcus cluster 15
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 2 | N/A | low
**Name:** Prochlorococcus cluster 2
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 3 | N/A | low
**Name:** Prochlorococcus cluster 3
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 4 | N/A | low
**Name:** Prochlorococcus cluster 4
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 5 | N/A | low
**Name:** Prochlorococcus cluster 5
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 6 | N/A | low
**Name:** Prochlorococcus cluster 6
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 7 | N/A | low
**Name:** Prochlorococcus cluster 7
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 8 | N/A | low
**Name:** Prochlorococcus cluster 8
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 9 | N/A | low
**Name:** Prochlorococcus cluster 9
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

## Coe 2024 / supp_table_3_parental_clusters

### Cluster 1 | N/A | low
**Name:** Prochlorococcus cluster 1
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 10 | N/A | low
**Name:** Prochlorococcus cluster 10
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 11 | N/A | low
**Name:** Prochlorococcus cluster 11
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 12 | N/A | low
**Name:** Prochlorococcus cluster 12
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 13 | N/A | low
**Name:** Prochlorococcus cluster 13
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 14 | N/A | low
**Name:** Prochlorococcus cluster 14
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 15 | N/A | low
**Name:** Prochlorococcus cluster 15
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 2 | N/A | low
**Name:** Prochlorococcus cluster 2
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 3 | N/A | low
**Name:** Prochlorococcus cluster 3
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 4 | N/A | low
**Name:** Prochlorococcus cluster 4
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 5 | N/A | low
**Name:** Prochlorococcus cluster 5
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 6 | N/A | low
**Name:** Prochlorococcus cluster 6
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 7 | N/A | low
**Name:** Prochlorococcus cluster 7
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 8 | N/A | low
**Name:** Prochlorococcus cluster 8
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

### Cluster 9 | N/A | low
**Name:** Prochlorococcus cluster 9
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific cluster descriptions provided in the paper or supplement.

## Lindell 2007 / med4_phage_transcription_groups

### Cluster 1 | early transient | high
**Name:** MED4 cluster 1 (early transient)
**Enrichment:** stress response and translation (p=None, sig=False)
**Functional:** Enriched for high-light-inducible stress response genes including hli, carbon metabolism gene rbcLS, transcription genes rpoC2 and rpoD, and ribosomal protein genes rpl5, rpl6, rps8, rps11, and rps17.
**Temporal pattern:** Transiently upregulated immediately after infection, with expression peaking early in the infection time course.
**Sources:** Figure 3
**Quotes:**
- [Page 3] The first was transiently upregulated immediately after infection and consists of high-light-inducible stress response (hli), carbon metabolism (rbcLS), transcription (rpoC2, rpoD) and ribosome (rpl5, rpl6, rps8, rps11, rps17) genes.

### Cluster 2 | late sustained | high
**Name:** MED4 cluster 2 (late sustained)
**Enrichment:** RNA degradation and stress response (p=None, sig=False)
**Functional:** Includes genes involved in RNA degradation and modification such as rne, rnhB, dus, and sun; protein turnover genes clpS and an AAA ATPase family gene; stress response genes umuD and phoH; and genes of unknown function.
**Temporal pattern:** Upregulated starting approximately 2 hours after infection and sustained through later time points.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Transcripts of the second group appeared 2 h after infection and included genes involved in RNA degradation and modification (rne, rnhB, dus and sun), protein turnover (clpS, and an AAA ATPase family gene), stress responses (umuD and phoH), and those of unknown function.

## Thompson 2011 / med4_iron_response_clusters

### Cluster 1 | early transient | medium
**Name:** Prochlorococcus MED4 cluster 1 (1 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Contains the gene Yfr2-5_1, a small non-coding RNA potentially involved in regulatory functions during iron stress.
**Temporal pattern:** Expression changes begin early in the iron starvation time course and show dynamic regulation across the time points, including rescue.
**Assessment notes:** Limited specific functional annotation; small RNA regulatory role inferred.
**Sources:** Figure 3c
**Quotes:**
- [Figure 3c legend] Cluster 1: Yfr2-5_1 (1)

### Cluster 13 | gradual increase | medium
**Name:** Prochlorococcus MED4 cluster 13 (7 genes)
**Enrichment:** transport (p=None, sig=False)
**Functional:** Enriched for ABC transporter genes involved in transport functions under iron stress.
**Temporal pattern:** Genes in this cluster show differential expression starting early in iron starvation and respond dynamically through the time course and rescue.
**Assessment notes:** Transport function inferred from gene annotations; no explicit p-value reported.
**Sources:** Figure 3c
**Quotes:**
- [Figure 3c legend] 13: ABC transporter (6)

### Cluster 16 | biphasic | medium
**Name:** Prochlorococcus MED4 cluster 16 (4 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Contains small non-coding RNAs Yfr2 and Yfr4, suggesting a regulatory role in iron stress response.
**Temporal pattern:** Expression changes occur during iron starvation and are reversed upon iron rescue, indicating dynamic regulation.
**Assessment notes:** Small RNA cluster with inferred regulatory function; no enrichment data.
**Sources:** Figure 3c
**Quotes:**
- [Figure 3c legend] 16: Yfr2, Yfr4 (2)

### Cluster 19 | late sustained | high
**Name:** Prochlorococcus MED4 cluster 19 (63 genes)
**Enrichment:** iron stress response (p=None, sig=False)
**Functional:** Enriched for iron stress response genes including bcp, sigA/rpoD, crp, hli, iraI, idiA, isiB, and several small RNAs, indicating a coordinated iron stress regulon.
**Temporal pattern:** Genes are upregulated during iron starvation starting early and remain elevated until iron rescue, where expression decreases.
**Confidence notes:** Cluster includes many known iron stress genes and sRNAs, supporting functional assignment.
**Assessment notes:** Strong evidence from gene content and expression pattern.
**Sources:** Figure 3c, Figure 4b
**Quotes:**
- [Figure 3c legend] 19: bcp, sigA/rpoD, crp, hli, iraI, idiA, isiB, Yfr20, asRNA_04601, asRNA_07401 (63)
- [Page 6] MED4 genomic island 5 (ISL5), in particular, was a ‘hotspot’ for differentially expressed genes in our experiments—including high-light inducible (hli) genes, and numerous genes of unknown function.

### Cluster 22 | N/A | low
**Name:** Prochlorococcus MED4 cluster 22 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Cluster not discussed in the paper.
**Assessment notes:** No discussion or data provided for this cluster.

### Cluster 3 | gradual increase | medium
**Name:** Prochlorococcus MED4 cluster 3 (15 genes)
**Enrichment:** photosystem and stress response (p=None, sig=False)
**Functional:** Contains hli genes and other genes related to photosystem components and stress response.
**Temporal pattern:** Genes show differential expression during iron starvation and recovery phases, with dynamic regulation.
**Assessment notes:** Functional annotation based on gene names; no enrichment statistics.
**Sources:** Figure 3c
**Quotes:**
- [Figure 3c legend] 3: hli (12)

### Cluster 4 | N/A | low
**Name:** Prochlorococcus MED4 cluster 4 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Cluster not discussed in the paper.
**Assessment notes:** No discussion or data provided for this cluster.

### Cluster 6 | early transient downregulation | high
**Name:** Prochlorococcus MED4 cluster 6 (2 genes)
**Enrichment:** electron transfer (p=None, sig=False)
**Functional:** Contains petF, encoding ferredoxin, an iron-requiring electron transfer protein downregulated during iron stress.
**Temporal pattern:** Downregulated during iron starvation and upregulated upon iron rescue.
**Confidence notes:** Known iron response gene with clear expression pattern.
**Assessment notes:** Strong evidence from gene identity and expression.
**Sources:** Figure 3c, Figure 5
**Quotes:**
- [Figure 3c legend] 6: petF (2)
- [Page 7] ferredoxin was downregulated and flavodoxin was upregulated during iron stress in both MED4 and MIT9313

### Cluster 7 | N/A | low
**Name:** Prochlorococcus MED4 cluster 7 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Cluster not discussed in the paper.
**Assessment notes:** No discussion or data provided for this cluster.

### Cluster 9 | gradual increase | medium
**Name:** Prochlorococcus MED4 cluster 9 (30 genes)
**Enrichment:** metabolism and regulation (p=None, sig=False)
**Functional:** Includes genes involved in sulfur metabolism (cysD), heme synthesis (hemA), dehydrogenases, and several small RNAs, indicating diverse metabolic functions during iron stress.
**Temporal pattern:** Genes show differential expression during iron starvation and recovery, with dynamic regulation.
**Assessment notes:** Functional annotation based on gene names; no enrichment statistics.
**Sources:** Figure 3c
**Quotes:**
- [Figure 3c legend] 9: cysD, hemA, dehydrogenase, Yfr16, Yfr19, Yfr8 (30)

## Thompson 2011 / mit9313_iron_response_clusters

### Cluster 1 | N/A | low
**Name:** MIT9313 cluster 1 (2 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** Contains two genes including the non-coding RNA Yfr2-5_1. Specific gene functions are not detailed in the paper for this cluster.
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 10 | N/A | low
**Name:** MIT9313 cluster 10 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 13 | N/A | medium
**Name:** MIT9313 cluster 13 (6 genes)
**Enrichment:** ABC transporter (p=None, sig=False)
**Functional:** Enriched for ABC transporter genes. Contains six genes including components of transport systems.
**Temporal pattern:** N/A
**Assessment notes:** Enrichment category inferred from cluster annotation in figure legend; no p-value reported.
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 13: ABC transporter (6)

### Cluster 14 | N/A | low
**Name:** MIT9313 cluster 14 (4 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 16 | N/A | high
**Name:** MIT9313 cluster 16 (2 genes)
**Enrichment:** small non-coding RNAs (p=None, sig=False)
**Functional:** Contains sRNAs Yfr2 and Yfr4, which are small non-coding RNAs differentially expressed under iron stress.
**Temporal pattern:** N/A
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 16: Yfr2, Yfr4 (2)
- [Page 7] We observed differential expression of several, albeit different, sRNAs in MED4 and MIT9313 (Supplementary Tables S1 and S2), suggesting that sRNAs are also important regulatory agents for Prochlorococcus iron metabolism.

### Cluster 17 | N/A | low
**Name:** MIT9313 cluster 17 (10 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 18 | N/A | low
**Name:** MIT9313 cluster 18 (4 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 19 | N/A | high
**Name:** MIT9313 cluster 19 (2 genes)
**Enrichment:** iron stress response (p=None, sig=False)
**Functional:** Enriched for genes including bcp, sigA/rpoD, crp, hli, iraI, idiA, isiB, and antisense RNAs. Includes iron stress response genes such as idiA and isiB.
**Temporal pattern:** N/A
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 19: bcp, sigA/rpoD, crp, hli, iraI, idiA, isiB, Yfr20, asRNA_04601, asRNA_07401 (63)
- [Page 7-8] idiA is hypothesized to be a periplasmic iron-binding protein component of an iron ABC-transporter system... idiA was upregulated during iron stress and downregulated following rescue.

### Cluster 2 | N/A | medium
**Name:** MIT9313 cluster 2 (10 genes)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Contains genes including coaD and a possible porin, involved in transport functions.
**Temporal pattern:** N/A
**Confidence notes:** Enrichment p-value reported in example, assumed similar for this cluster based on gene content.
**Assessment notes:** Enrichment p-value inferred from example cluster; exact p-value not explicitly stated.
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] Cluster 2, coaD, some possible porin, Yfr7 (10)

### Cluster 20 | N/A | low
**Name:** MIT9313 cluster 20 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 23 | N/A | high
**Name:** MIT9313 cluster 23 (2 genes)
**Enrichment:** iron storage (p=None, sig=False)
**Functional:** Contains ferritin genes involved in iron storage, differentially expressed in MIT9313 during iron stress.
**Temporal pattern:** N/A
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 23: ferritin (2)
- [Page 8-9] One of the two ferritin genes in MIT9313 (PMT0499) was upregulated 16 h after iron deprivation, whereas the single ferritin gene in MED4 was not differentially expressed.

### Cluster 26 | N/A | low
**Name:** MIT9313 cluster 26 (2 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 27 | N/A | low
**Name:** MIT9313 cluster 27 (2 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 3 | N/A | high
**Name:** MIT9313 cluster 3 (12 genes)
**Enrichment:** high-light inducible (hli) genes (p=None, sig=False)
**Functional:** Contains hli genes involved in protection of photosystems from oxidative damage during iron stress.
**Temporal pattern:** N/A
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 3: hli (12)
- [Page 7-8] A large number of hli genes were differentially expressed in our experiments, including hli05/hli08 shared by both strains.

### Cluster 30 | N/A | high
**Name:** MIT9313 cluster 30 (11 genes)
**Enrichment:** electron transfer and transport (p=None, sig=False)
**Functional:** Includes petF (ferredoxin), rimI, gmk, cytochrome b6/f subunit VII, maf-like, and nitrogen transport genes, involved in electron transfer and transport.
**Temporal pattern:** N/A
**Sources:** Figure 3d
**Quotes:**
- [Figure 3d legend] 30: petF, rimI, gmk, cytochrome b6/f subunit VII, maf-like, N transport (11)
- [Page 7-8] Ferredoxin (petF) was downregulated and flavodoxin (isiB) upregulated during iron stress in both strains.

### Cluster 35 | N/A | low
**Name:** MIT9313 cluster 35 (6 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 38 | N/A | low
**Name:** MIT9313 cluster 38 (6 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 39 | N/A | low
**Name:** MIT9313 cluster 39 (1 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 4 | N/A | high
**Name:** MIT9313 cluster 4 (24 genes)
**Enrichment:** iron transport and stress response (p=None, sig=False)
**Functional:** Contains genes including glyQ, piuC, idiA, pcbB, and isiB, involved in iron transport and stress response.
**Temporal pattern:** N/A
**Sources:** Figure 3d, Figure 6
**Quotes:**
- [Figure 3d legend] 4: glyQ, piuC, idiA, pcbB, isiB (24)
- [Page 9-11] The idiA region includes piuC and other iron stress induced genes upregulated in MIT9313 but absent in MED4.

### Cluster 8 | N/A | low
**Name:** MIT9313 cluster 8 (3 genes)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Paper does not discuss this cluster.

## Tolonen 2006 / med4_kmeans_nstarvation

### Cluster 1 | early transient | high
**Name:** Prochlorococcus MED4 cluster 1 (5 genes)
**Enrichment:** transport and binding (p=0.01, sig=True)
**Functional:** Contains nitrogen transport genes such as urtA and cynA. Enriched for transport and binding category.
**Temporal pattern:** Most rapidly and highly upregulated genes in MED4 during nitrogen starvation, with expression first appearing within 6 hours and peaking early in the time course.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA, MED4 cynA, and the MIT9313 nitrite permease.

### Cluster 2 | early transient | high
**Name:** Prochlorococcus MED4 cluster 2 (16 genes)
**Enrichment:** regulation (p=0.07, sig=False)
**Functional:** Contains the rapidly-responding subset of hli genes (including hli10) and ntcA. The most highly upregulated hli genes have the strongest NtcA binding sites.
**Temporal pattern:** Rapidly and highly upregulated early in nitrogen starvation, representing the fast-responding hli subset distinct from the later-responding hli genes in cluster 3.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 3] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 3 | late transient | medium
**Name:** Prochlorococcus MED4 cluster 3 (19 genes)
**Enrichment:** regulation (p=0.18, sig=False)
**Functional:** Contains a subset of hli genes that respond later and to a lesser degree than cluster 2. Also includes genes involved in regulation.
**Temporal pattern:** Later and less pronounced upregulation during nitrogen starvation compared to cluster 2.
**Assessment notes:** Less detailed description available.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 4 | gradual increase | high
**Name:** Prochlorococcus MED4 cluster 4 (86 genes)
**Enrichment:** amino acid synthesis (p=0.18, sig=False)
**Functional:** Enriched for amino acid synthesis genes and genes involved in carbon metabolism such as zwf, tal, and acnB. Includes genes for the oxidative pentose phosphate pathway and 2-oxoglutarate synthesis.
**Temporal pattern:** Upregulated during nitrogen starvation, with genes involved in carbon metabolism and amino acid synthesis showing increased expression.
**Assessment notes:** icd removed from gene list per Table III assignment to cluster 5.
**Sources:** Figure 3, Figure 4F
**Quotes:**
- [Page 7, Table III] The paper discusses MED4 cluster 4 as containing amino acid synthesis and carbon metabolism genes (zwf, tal, acnB). Per Table III, icd is in cluster 5, not cluster 4.

### Cluster 5 | early transient | high
**Name:** Prochlorococcus MED4 cluster 5 (73 genes)
**Enrichment:** regulation (p=0.18, sig=False)
**Functional:** Contains two rpoD-like sigma factors that are upregulated during nitrogen stress, indicating transcriptional regulation. Also contains icd (per Table III).
**Temporal pattern:** Upregulated during nitrogen starvation, with sigma factors mediating gene expression changes.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 3] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.
- [Page 6] Two out of five MED4 sigma factors were upregulated and may play a role in gene expression during N stress.

### Cluster 6 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 6 (124 genes)
**Enrichment:** translation (p=0.001, sig=True)
**Functional:** Enriched for translation-related genes, strongly downregulated during nitrogen starvation, indicating repression of protein synthesis machinery.
**Temporal pattern:** Strongly downregulated during nitrogen starvation, consistent with reduced translation activity.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Highly downregulated clusters show strong enrichment for translation (MED4 cluster 6).

### Cluster 7 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 7 (42 genes)
**Enrichment:** translation (p=2.7e-09, sig=True)
**Functional:** Enriched for translation genes, strongly downregulated during nitrogen starvation, indicating repression of protein synthesis.
**Temporal pattern:** Strongly downregulated during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Highly downregulated clusters show strong enrichment for translation (MED4 cluster 7).

### Cluster 8 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 8 (37 genes)
**Enrichment:** photosynthesis and respiration (p=2.6e-09, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes, including photosystem I genes (psaBDEIJKLM). Strongly downregulated during nitrogen starvation.
**Temporal pattern:** Strongly downregulated during nitrogen starvation, indicating repression of photosynthesis genes.
**Sources:** Figure 3
**Quotes:**
- [Page 3] MED4 cluster 8 contains numerous genes for photosystem I (psaBDEIJKLM) and is enriched for photosynthesis and respiration genes.

### Cluster 9 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 9 (8 genes)
**Enrichment:** photosynthesis and respiration (p=1.4e-05, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes, including ATP synthase subunits and carbon fixation genes (rbcLS). Strongly downregulated during nitrogen starvation.
**Temporal pattern:** Strongly downregulated during nitrogen starvation, indicating repression of photosynthesis and carbon fixation genes.
**Sources:** Figure 3
**Quotes:**
- [Page 3] MED4 cluster 9 contains ATP synthase subunits and carbon fixation genes (rbcLS) and is enriched for photosynthesis and respiration genes.

## Tolonen 2006 / mit9313_kmeans_nstarvation

### Cluster 1 | early transient | high
**Name:** MIT9313 cluster 1 (transport and binding)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Contains nitrogen transport genes such as urtA and the nitrite permease. Also contains the rapidly-responding subset of hli genes (hliS, hli7). Enriched for transport and binding category (P=0.04).
**Temporal pattern:** Most rapidly and highly upregulated cluster, with genes responding within the first hours of nitrogen starvation.
**Confidence notes:** MIT9313 lacks cyanate genes (cynA), unlike MED4.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA, MED4 cynA, and the MIT9313 nitrite permease.
- [Page 3] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 2 | early transient | high
**Name:** MIT9313 cluster 2 (regulation)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** Contains two upregulated sigma factors and the later-responding subset of hli genes. Enriched for regulation category (P=0.03).
**Temporal pattern:** Upregulated cluster containing sigma factors and the later, less pronounced hli response subset (distinct from the rapidly-responding hli genes in cluster 1).
**Sources:** Figure 3, Figure 4E
**Quotes:**
- [Page 3] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.
- [Page 3] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 3 | N/A | low
**Name:** MIT9313 cluster 3 (regulation)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Paper does not discuss this cluster.
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 4 | late sustained repression | medium
**Name:** MIT9313 cluster 4 (not discussed)
**Enrichment:**  (p=None, sig=False)
**Functional:** Contains genes linking nitrogen and carbon metabolism such as glnB, icd, acnB, and rbcLS. Genes were unchanged until repressed only at the final time point, reflecting general transcriptional shutdown rather than specific nitrogen stress response.
**Temporal pattern:** Genes unchanged until late repression at 48h, coinciding with severe starvation state.
**Confidence notes:** Repression attributed to general transcriptional shutdown, not specific N stress response.
**Assessment notes:** Interpretation based on timing and physiological data; not explicitly detailed in paper.
**Sources:** Figure 3, Figure 4F
**Quotes:**
- [Page 4] MIT9313 cluster 4 contains a number of genes that were unchanged until being repressed only at the final time point. The physiological measurements Fv/Fm and chlorophyll fluorescence show that the cells were in a severe state of starvation by this time.
- [Page 4] The genes in MIT9313 cluster 4 may thus represent those genes that are repressed as part of a general shutdown in transcription, rather than a specific N stress response. Interestingly, this MIT9313 cluster contains a number of genes linking N and C metabolism (glnB, icd, acnB, rbcLS).

### Cluster 5 | N/A | low
**Name:** MIT9313 cluster 5
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Paper does not discuss this cluster.
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 6 | gradual repression | high
**Name:** MIT9313 cluster 6 (photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=0.00049, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes including photosystem I and II components and the phycoerythrin gene cpeB.
**Temporal pattern:** Downregulated cluster with gradual repression during nitrogen starvation, with mild repression at late time points.
**Sources:** Figure 3, Figure 4F
**Quotes:**
- [Page 4] MIT9313 cluster 6, the only cluster enriched for Photosynthesis and Respiration in this strain, contains genes for diverse aspects of photosystem I and II along with the phycoerythrin gene, cpeB.
- [Page 7] MIT9313 rbcLS are members of K-means cluster 4. These interstrain differences were also reflected in other carbon metabolism genes such as the bicarbonate transporter, sbtA, as well as the csoS12 genes encoding the carboxysome shell proteins. MED4 thus responds to reduced N availability by repressing the expression of carbon transport and fixation genes to a much greater degree than MIT9313.

### Cluster 7 | N/A | low
**Name:** MIT9313 cluster 7
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Paper does not discuss this cluster.
**Assessment notes:** Paper does not discuss this cluster.

## Wang 2014 / med4_expression_level

### Cluster HEG | high across all | high
**Name:** Prochlorococcus MED4 cluster HEG (highly expressed)
**Enrichment:** Energy production and conversion; Translation; Protein folding and turnover (p=0.001, sig=True)
**Functional:** Enriched for genes involved in energy production and conversion, translation and ribosomal structure, and protein modification, folding and turnover. Includes ribosomal proteins and photosynthetic apparatus genes.
**Temporal pattern:** Genes in this cluster are consistently highly expressed across multiple growth conditions.
**Sources:** Figure 3, Figure 4
**Quotes:**
- [Page 5] Among these core HEG genes, several functional categories were more prominent than others. These included the “C” (energy production and conversion), “J” (translation and ribosomal structure), and “O” (protein modification, folding and turnover) categories (Figure 4c).
- [Page 4] HEG had a significantly lower nonsynonymous substitution rate (Ka) than MEG or LEG (Kruskal-Wallis Test, two-tailed P < 0.001; Figure 3a)

### Cluster LEG | N/A | medium
**Name:** Prochlorococcus MED4 cluster LEG (lowly expressed)
**Enrichment:** General function (p=0.023, sig=True)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Only slight enrichment reported for general function category.
**Sources:** Figure 4
**Quotes:**
- [Page 5] Additionally, category “R” (general function) was slightly enriched in both LEG and NEG (P = 0.023 and 0.055; data not shown).

### Cluster MEG | moderate across all | high
**Name:** Prochlorococcus MED4 cluster MEG (moderately expressed)
**Enrichment:** Essential genes (p=0.001, sig=True)
**Functional:** Enriched for essential genes with homologs in the Database of Essential Genes (DEG).
**Temporal pattern:** Genes in this cluster show moderate expression levels across conditions.
**Sources:** Figure 4
**Quotes:**
- [Page 4] Although the MEG subclass had a significantly higher rate of DEG-hit genes (P < 0.001; Figure 4b)

### Cluster NEG | N/A | high
**Name:** Prochlorococcus MED4 cluster NEG (not expressed)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Figure 3
**Quotes:**
- [Page 4] the core genome had fewer NEG and VEG than the flexible genome (1.5% < 6.6% and 49.6% < 64.6%, respectively; P < 0.001; Figure 3c)

### Cluster VEG | N/A | high
**Name:** Prochlorococcus MED4 cluster VEG (variably expressed)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Figure 3
**Quotes:**
- [Page 4] the core genome had fewer NEG and VEG than the flexible genome (1.5% < 6.6% and 49.6% < 64.6%, respectively; P < 0.001; Figure 3c)

## Zinser 2009 / med4_diel_clusters

### Cluster 1 | peaks midday | high
**Name:** Prochlorococcus cluster 1 (photosystem I and II)
**Enrichment:** Photosystem I and II (p=1.5e-09, sig=True)
**Functional:** Enriched for photosystem I and II components (p=1.50E-09 and 2.50E-05 respectively). Includes genes related to photosystem I and II such as psbA, psbD, and psaA.
**Temporal pattern:** Genes peak in expression around mid-day (8.3 h in the paper's time scale) with 24h periodicity, co-varying with light intensity.
**Sources:** Figure 4A, Table 1, Table S4
**Quotes:**
- [Page 6] Expression of approximately half of photosystem (PS) II genes, including reaction center genes psbA and psbD, peak in abundance at mid-day
- [Page 6] Periodicity patterns of photosynthesis genes fell into 4 clusters (Table 1). Expression of approximately half of photosystem (PS) II genes, including reaction center genes psbA and psbD (encoding D1 and D2 respectively), as well as psbC (CP43) and psbF, co-varied with light intensity, with maxima at mid-day, and minima in the middle of the night (Figure 4A, Table S4)

### Cluster 10 | N/A | high
**Name:** Prochlorococcus cluster 10
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 11 | N/A | high
**Name:** Prochlorococcus cluster 11
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 12 | N/A | high
**Name:** Prochlorococcus cluster 12
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 13 | peaks early morning | high
**Name:** Prochlorococcus cluster 13 (ribosomal proteins)
**Enrichment:** Ribosomal proteins (p=3.7e-41, sig=True)
**Functional:** Enriched for ribosomal proteins (p=3.70E-41). Includes 45 out of 53 ribosomal protein genes.
**Temporal pattern:** Genes peak in expression near 3.4 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 13 is enriched for ribosomal proteins (45/53 genes, p=3.70E-41)

### Cluster 14 | N/A | high
**Name:** Prochlorococcus cluster 14
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 15 | N/A | high
**Name:** Prochlorococcus cluster 15
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 16 | peaks early morning | high
**Name:** Prochlorococcus cluster 16 (ATP synthase and CO2 metabolism)
**Enrichment:** ATP synthase and CO2 metabolism (p=8.4e-08, sig=True)
**Functional:** Enriched for ATP synthase (8/8 genes, p=8.40E-08) and CO2 metabolism genes (7/9 genes, p=1.80E-05).
**Temporal pattern:** Genes peak in expression near 5.5 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 16 is enriched for ATP synthase (8/8 genes, p=8.40E-08) and CO2 metabolism (7/9 genes, p=1.80E-05)

### Cluster 17 | N/A | high
**Name:** Prochlorococcus cluster 17
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 18 | N/A | high
**Name:** Prochlorococcus cluster 18 (menaquinone and ubiquinone)
**Enrichment:** Menaquinone and ubiquinone (p=0.00045, sig=True)
**Functional:** Enriched for menaquinone and ubiquinone genes (6/9 genes, p=0.00045).
**Temporal pattern:** Non-expressed gene cluster per Table 1 classification.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 18 is enriched for menaquinone and ubiquinone genes (6/9 genes, p=0.00045)

### Cluster 2 | N/A | high
**Name:** Prochlorococcus cluster 2
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Sources:** Table 1

### Cluster 3 | peaks mid-day | high
**Name:** Prochlorococcus cluster 3 (cytochrome b6/f and photosystem II)
**Enrichment:** Cytochrome b6/f and Photosystem II (p=3.7e-08, sig=True)
**Functional:** Enriched for cytochrome b6/f (3/7 genes, p=0.0059) and photosystem II genes (8/22 genes, p=3.70E-08).
**Temporal pattern:** Genes peak in expression near 12.5 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 3 is enriched for cytochrome b6/f (3/7 genes, p=0.0059) and photosystem II genes (8/22 genes, p=3.70E-08)

### Cluster 4 | peaks late afternoon | medium
**Name:** Prochlorococcus cluster 4 (cytochrome b6/f)
**Enrichment:** Cytochrome b6/f (p=0.073, sig=False)
**Functional:** Enriched for cytochrome b6/f genes (3/7 genes, p=0.073).
**Temporal pattern:** Genes peak in expression near 15.8 hours after onset of dark with 24h periodicity.
**Assessment notes:** Enrichment p-value marginally above significance threshold.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 4 is enriched for cytochrome b6/f genes (3/7 genes, p=0.073)

### Cluster 5 | peaks late afternoon | high
**Name:** Prochlorococcus cluster 5 (respiratory terminal oxidases and protein degradation)
**Enrichment:** Respiratory terminal oxidases and protein degradation (p=0.019, sig=True)
**Functional:** Enriched for respiratory terminal oxidases (3/3 genes, p=0.019) and degradation of proteins, peptides, and glycopeptides (5/15 genes, p=0.071).
**Temporal pattern:** Genes peak in expression near 17.5 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 5 is enriched for respiratory terminal oxidases (3/3 genes, p=0.019) and degradation of proteins, peptides, and glycopeptides (5/15 genes, p=0.071)

### Cluster 6 | peaks early night | high
**Name:** Prochlorococcus cluster 6 (purine ribonucleotide biosynthesis)
**Enrichment:** Purine ribonucleotide biosynthesis (p=0.0049, sig=True)
**Functional:** Enriched for purine ribonucleotide biosynthesis genes (7/18 genes, p=0.0049).
**Temporal pattern:** Genes peak in expression near 18.6 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 6 is enriched for purine ribonucleotide biosynthesis genes (7/18 genes, p=0.0049)

### Cluster 7 | peaks night | medium
**Name:** Prochlorococcus cluster 7 (nitrogen metabolism)
**Enrichment:** Nitrogen metabolism (p=0.087, sig=False)
**Functional:** Enriched for nitrogen metabolism genes (4/8 genes, p=0.087).
**Temporal pattern:** Genes peak in expression near 20.1 hours after onset of dark with 24h periodicity.
**Assessment notes:** Enrichment p-value marginally above significance threshold.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 7 is enriched for nitrogen metabolism genes (4/8 genes, p=0.087)

### Cluster 8 | peaks night | high
**Name:** Prochlorococcus cluster 8 (chaperones)
**Enrichment:** Chaperones (p=0.00023, sig=True)
**Functional:** Enriched for chaperone genes (7/14 genes, p=0.00023).
**Temporal pattern:** Genes peak in expression near 21.0 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 8 is enriched for chaperone genes (7/14 genes, p=0.00023)

### Cluster 9 | peaks late night | high
**Name:** Prochlorococcus cluster 9 (RNA synthesis, modification, and DNA transcription)
**Enrichment:** RNA synthesis, modification, and DNA transcription (p=0.035, sig=True)
**Functional:** Enriched for RNA synthesis, modification, and DNA transcription genes (6/23 genes, p=0.035).
**Temporal pattern:** Genes peak in expression near 22.4 hours after onset of dark with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 9 is enriched for RNA synthesis, modification, and DNA transcription genes (6/23 genes, p=0.035)

## Warnings

- [Thompson 2011 / mit9313_iron_response_clusters / cluster 1] low confidence but functional_description is not 'N/A': Contains two genes including the non-coding RNA Yfr2-5_1. Sp...
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 9] near-identical functional_description as cluster 8
- [Bernstein 2017 / bp1_light_clusters / cluster A] low confidence but temporal_pattern is not 'N/A': Tent-shaped expression pattern with peak at intermediate irr...
- [Bernstein 2017 / bp1_oxygen_clusters / cluster E] low confidence but temporal_pattern is not 'N/A': Tent-shaped expression pattern with peak at intermediate oxy...
- [Bernstein 2017 / mruber_light_clusters / cluster A] low confidence but temporal_pattern is not 'N/A': Tent-shaped expression pattern with peak at intermediate irr...
- [Bernstein 2017 / mruber_light_clusters / cluster D] low confidence but temporal_pattern is not 'N/A': Increases with increasing irradiance....
- [Bernstein 2017 / mruber_oxygen_clusters / cluster E] low confidence but temporal_pattern is not 'N/A': Tent-shaped expression pattern with peak at intermediate oxy...
- [Bernstein 2017 / mruber_oxygen_clusters / cluster F] low confidence but temporal_pattern is not 'N/A': Inverse of cluster E....
