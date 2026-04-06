# Cluster Extraction Report

## Alonso-Saez 2023 / mit9301_softclusters_thermal_acclimation

### Cluster Cluster A | high across all temperatures during daytime | high
**Name:** Prochlorococcus cluster A (core metabolism)
**Enrichment:** Carbon fixation and metabolism (p=None, sig=False)
**Functional:** Includes genes related to carbon fixation and assimilation such as RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), Calvin cycle enzymes (gap2, tktA, glpX, pgk, cbbA), glycogen synthesis (glgABC), ATP synthesis genes (atpADE), and some photosystem II components (psbA, psbC, psbD).
**Temporal pattern:** Consistently expressed during daytime across the thermal gradient with stable expression levels.
**Confidence notes:** High confidence based on gene function and expression pattern consistency.
**Sources:** Figure 2
**Quotes:**
- [Page 4] Cluster A included genes related to C fixation and assimilation, such as RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), the Calvin cycle (e.g., gap2, tktA, glpX, pgk, and cbbA), and glycogen synthesis (glgABC), consistently expressed during daytime along the thermal gradient.

### Cluster Cluster B | down at cold during daytime | high
**Name:** Prochlorococcus cluster B (photosynthesis down at cold)
**Enrichment:** Photosynthesis (p=None, sig=False)
**Functional:** Contains photosynthetic genes including components of photosystem I (psaABDEFKL), some photosystem II components (psbBJH, psbO), and photosynthetic electron transport genes (petACGNM). These genes show a decreasing expression trend from optimum to minimum temperature during daytime.
**Temporal pattern:** Gradual decrease in expression from optimum to minimum temperature during daytime, correlating with growth rate decline.
**Confidence notes:** Strong evidence from expression patterns and correlation with growth rates.
**Sources:** Figure 2, Figure 4
**Quotes:**
- [Page 7] The expression of all components of the PS I complex (psaABDEFKL) and some of PS II (including psbBJH and the oxygen-evolving complex protein psbO) showed a gradual decrease in expression from a temperature close to the optimum to the Tmin, in correlation with MIT9301 growth rates.

### Cluster Cluster C | up at cold during daytime | high
**Name:** Prochlorococcus cluster C (cold stress daytime)
**Enrichment:** Stress response (p=None, sig=False)
**Functional:** Includes genes involved in global stress response such as cellular chaperones (groEL/groES, dnaK, clpBCP), fatty acid desaturases (desA, desC), oxidative damage protection (recA, ruvB, sod), carotenoid synthesis (pds, crtBH), rubredoxin (rub), amino acid synthesis (glyA, serA, leuA), translation initiation factors (infABC), and nitrogen acquisition genes.
**Temporal pattern:** Strong upregulation at minimum temperature during daytime, indicating prioritization of stress response under light conditions.
**Confidence notes:** High confidence due to multiple stress-related genes and clear expression pattern.
**Sources:** Figure 2, Figure 3
**Quotes:**
- [Page 4] Clusters C and D genes included different elements of the global stress response, such as cellular chaperones (groES/groES, dnaK, and clpBCP) and fatty acid desaturases (desA and desC), as well as mechanisms against oxidative damage, such as DNA repair (recA and ruvB), superoxide dismutase (sod), and the synthesis of antioxidant compounds like carotenoids (pds and crtBH) and rubredoxin (rub). Notably, the expression of the chaperones groEL/groES, grpE, and htpG was strongly upregulated at the Tmin only during daytime, suggesting a prioritization of their expression during the light-exposed period.

### Cluster Cluster D | up at cold during day and night | high
**Name:** Prochlorococcus cluster D (cold stress all time)
**Enrichment:** Stress response (p=None, sig=False)
**Functional:** Contains genes involved in global stress response similar to cluster C, including chaperones and oxidative stress protection genes, with upregulation at minimum temperature during both daytime and nighttime.
**Temporal pattern:** Upregulated at minimum temperature during both daytime and nighttime, indicating sustained cold stress response.
**Confidence notes:** High confidence based on expression pattern and gene functions.
**Sources:** Figure 2, Figure 3
**Quotes:**
- [Page 4] Clusters C and D were associated with mechanisms of cold stress response, as they were characterized by a strong upregulation at the Tmin either during daytime or during both daytime and nighttime, respectively.

### Cluster Cluster E | high at night across all temperatures | high
**Name:** Prochlorococcus cluster E (core night metabolism)
**Enrichment:** DNA replication and cell division (p=None, sig=False)
**Functional:** Includes genes related to catabolic consumption (cyoB, ndhD), DNA replication (dnaA, nrdJ, gyrB), cell division (ftsZYQ), and pentose phosphate pathway (tal, gnd, zwf), typically upregulated at nighttime.
**Temporal pattern:** Upregulated at nighttime across the thermal gradient, representing essential metabolic pathways maintained under all temperature conditions.
**Confidence notes:** High confidence based on gene functions and consistent nighttime expression.
**Sources:** Figure 2
**Quotes:**
- [Page 4] Cluster E included genes related to catabolic consumption (cyoB and ndhD), DNA replication (dnaA, nrdJ, and gyrB), cell division (ftsZYQ), and the pentose phosphate pathway (tal, gnd, and zwf), all of them upregulated at nighttime.

## Bernstein 2017 / bp1_light_clusters

### Cluster A | peak at mid irradiance | high
**Name:** Thermosynechococcus elongatus cluster A (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint irradiance (1,190 µmol photons m-2 s-1). Genes in cluster A responded inversely to those in cluster B. Specific genes are not individually named for this cluster in the paper.
**Temporal pattern:** Maximum transcript abundance at intermediate irradiance (1,190 µmol photons m-2 s-1), with a tent-shaped pattern showing increase then decrease as irradiance changes.
**Sources:** Figure 6
**Quotes:**
- [Page 7] Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint Ii (1,190 mol photons m2 s1).

### Cluster B | minimum at mid irradiance | high
**Name:** Thermosynechococcus elongatus cluster B (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Cluster B genes responded inversely to those in cluster A. Specific genes are not individually named for this cluster in the paper.
**Temporal pattern:** Transcript abundance decreases at intermediate irradiance and shows an inverse pattern to cluster A.
**Sources:** Figure 6
**Quotes:**
- [Page 7] Genes in cluster B responded inversely to those in cluster A.

### Cluster C | down with light, up with oxygen | high
**Name:** Thermosynechococcus elongatus cluster C (irradiance-responsive)
**Enrichment:** Vitamin B12 biosynthesis and methionine metabolism (p=0.05, sig=True)
**Functional:** Cluster C contains genes that show a relative decrease with increasing irradiance. Genes involved in vitamin B12 biosynthesis (cobWNT, cobOQDPC), vitamin B12 uptake/scavenging (btuCD), and methionine biosynthesis (metHX) are included. M. ruber genes encoding vitamin B12-dependent methionine synthesis also group here. This suggests potential exchange of vitamin B12 and methionine between species.
**Temporal pattern:** Transcript abundance decreases with increasing irradiance and increases with increasing oxygen tension, showing opposite patterns to clusters D and G.
**Confidence notes:** Enrichment p-value reported as significant for vitamin B12 biosynthesis and related functions.
**Sources:** Figure 6
**Quotes:**
- [Pages 7-8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis... Both T. elongatus and M. ruber expressed transcripts encoding the vitamin B12-dependent homocysteine methyltransferase (metH) but were negatively correlated and grouped across Ii and pO2 treatments into opposing clusters D (tll1027) and C (SY28_RS08890), indicating that methionine may be directly exchanged from T. elongatus as growth requirements increase with the specific growth rate.

### Cluster D | up with light, down with oxygen | high
**Name:** Thermosynechococcus elongatus cluster D (irradiance-responsive)
**Enrichment:** Carbon metabolism, nitrogen metabolism, ROS detoxification, methionine metabolism (p=0.05, sig=True)
**Functional:** Cluster D contains genes that show a relative increase with increasing irradiance. Genes involved in organic acid synthesis (acs, ackA, ldhA, gltA), exopolysaccharide synthesis (exoD), sucrose metabolism (spsA, sps), nitrogen metabolism (nrtABD, glnA, glnB, nirA, gltS), ROS detoxification (flv4, 2-Cys peroxiredoxins, prxQ-B2, sodB), and methionine degradation and salvage (ahcY, metK, mtaP, mtnA) are included. M. ruber genes involved in uptake and metabolism of related compounds also group here. This cluster corresponds to increased specific growth and photosynthesis rates.
**Temporal pattern:** Transcript abundance increases with increasing irradiance and decreases with increasing oxygen tension.
**Confidence notes:** Enrichment p-values significant for multiple metabolic functions including carbon fixation, nitrogen metabolism, and ROS detoxification.
**Sources:** Figure 6
**Quotes:**
- [Pages 7-9] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D)... Functions involved in synthesis and export of larger biomolecules (i.e., sugars, peptides, and extracellular polymeric substance [EPS]) also grouped into clusters D and G... T. elongatus nitrogen metabolism genes found in these clusters include glutamine synthetase gene glnA, nitrogen regulatory protein gene glnB, assimilatory ferredoxin-nitrate reductase gene nirA, and glutamate symporter gene gltS... The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H... The relative abundances of T. elongatus transcripts encoding the degradation and salvage of methionine (ahcY, metK, mtaP, and mtnA) increased with Ii but decreased with pO2 treatments and grouped within clusters D and G, respectively.

## Bernstein 2017 / bp1_oxygen_clusters

### Cluster E | down with pO2 | high
**Name:** Thermosynechococcus elongatus cluster E (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Contains 211 genes including those involved in nitrogen metabolism such as glnA, gltB, and vitamin B12-independent homocysteine methyltransferase (metE). Genes in this cluster generally decrease with increasing pO2 and increase with specific growth and photosynthesis rates. This cluster also includes T. elongatus genes involved in vitamin B12 biosynthesis and methionine metabolism that are negatively correlated with M. ruber genes in this cluster.
**Temporal pattern:** Transcript abundances decrease with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles analogous to cluster A (tent-shaped eigen-gene with maximum at midpoint irradiance).
**Sources:** Figure 6, Figure 3
**Quotes:**
- [Page 7] Clusters E to H contain transcripts that trend with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles that are analogous to clusters A to D.
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis.
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA (SY28_RS02395) and the large subunit of glutamate synthase (gltB; SY28_RS09480).

### Cluster F | down with pO2 | medium
**Name:** Thermosynechococcus elongatus cluster F (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Contains 187 genes including M. ruber genes involved in nitrogen uptake and metabolism, such as glutamate dehydrogenase (gdhA) and amino acid uptake systems. This cluster shows expression patterns that decrease with increasing pO2 and increase with specific growth and photosynthesis rates, indicating coordination of nitrogen metabolism between species.
**Temporal pattern:** Transcript abundances decrease with increasing pO2 treatments and increase with specific growth and photosynthesis rates, showing inverse patterns to clusters E and H.
**Assessment notes:** The paper discusses nitrogen metabolism genes in clusters D and G more explicitly; cluster F is less detailed but inferred from related gene groups.
**Sources:** Figure 6
**Quotes:**
- [Page 8] Several M. ruber nitrogen-associated genes also grouped into clusters D and/or G, including a glutamate dehydrogenase gene (gdhA) and amino acid uptake system genes.
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates.

### Cluster G | down with pO2 | high
**Name:** Thermosynechococcus elongatus cluster G (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Contains 913 genes including T. elongatus genes involved in organic acid synthesis and export (acs, ackA, ldhA, gltA), sucrose metabolism (spsA, sps), and exopolysaccharide synthesis (exoD). Also includes M. ruber genes involved in uptake and metabolism of related compounds such as acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), and amino acid uptake (livGH).
**Temporal pattern:** Transcript abundances decrease with increasing pO2 treatments and increase with specific growth and photosynthesis rates, showing profiles analogous to cluster D with respect to irradiance.
**Sources:** Figure 6
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). These included acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA).
- [Page 7-8] In conjunction with T. elongatus, M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. Notable examples included genes encoding putative acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA and gtsAB), xylose isomerase (xylA), xylulokinase (xylB), an ABC-type multisugar uptake system, and branched-chain amino acid uptake (livGH).

### Cluster H | up with pO2 | high
**Name:** Thermosynechococcus elongatus cluster H (pO2-responsive)
**Enrichment:** Oxidative stress response (p=0.01, sig=True)
**Functional:** Contains 888 genes including T. elongatus genes involved in oxidative stress responses such as flv4, 2-Cys peroxiredoxins (tll1454, tlr1289), periplasmic peroxiredoxin (prxQ-B2), and Mn-superoxide dismutase (sodB). M. ruber genes include the H2O2-responsive transcriptional regulator oxyR and manganese catalase. Genes associated with electron transfer processes potentiating ROS group into cluster G, which decreases with increasing pO2.
**Temporal pattern:** Transcript abundances increase with increasing pO2 treatments and decrease with specific growth and photosynthesis rates, showing profiles analogous to cluster C with respect to irradiance.
**Confidence notes:** Enrichment in oxidative stress response genes is statistically significant with p-value 0.01.
**Sources:** Figure 6, Figure 3
**Quotes:**
- [Page 8-9] The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H. Notable examples include an NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys family peroxiredoxins (tll1454 and tlr1289), a periplasmic peroxiredoxin (prxQ-B2), and an Mn-superoxide dismutase (sodB).
- [Page 9] M. ruber peroxidase (bcp) and superoxide dismutase (sod2) genes generally increased with pO2 and grouped with the T. elongatus genes into cluster H (increased with pO2 and decreasing growth rate).

## Bernstein 2017 / mruber_light_clusters

### Cluster A | peak at mid irradiance | high
**Name:** Meiothermus ruber cluster A (irradiance-responsive)
**Enrichment:** ROS/RNS detoxification (p=0.05, sig=True)
**Functional:** Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint irradiance (1,190 µmol photons m-2 s-1). Genes in this cluster show a relative decrease or increase with irradiance, corresponding to increasing specific growth and photosynthesis rates. M. ruber genes in this cluster include those involved in ROS/RNS detoxification such as peroxidase (bcp) and superoxide dismutase (sod2), which generally increase with oxygen tension and group with T. elongatus genes into cluster H (increased with pO2 and decreasing growth).
**Temporal pattern:** Maximum expression at intermediate irradiance (1,190 µmol photons m-2 s-1), with a tent-shaped pattern showing increase then decrease as irradiance changes.
**Confidence notes:** Cluster A includes genes with clear expression patterns correlated with irradiance and ROS detoxification, supported by transcriptomic data.
**Sources:** Figure 6, Figure 5
**Quotes:**
- [Page 6] Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint Ii (1,190 mol photons m2 s1).
- [Page 8] M. ruber peroxidase (bcp) and superoxide dismutase (sod2) genes responded differently to Ii treatments than did oxyR and related cyanobacterial profiles and were grouped into clusters A and C.

### Cluster B | N/A | low
**Name:** Meiothermus ruber cluster B (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data for cluster B genes in M. ruber.

### Cluster C | up with growth rate, down with irradiance | high
**Name:** Meiothermus ruber cluster C (irradiance-responsive)
**Enrichment:** Nitrogen metabolism, methionine and vitamin B12 metabolism (p=0.05, sig=True)
**Functional:** Genes in cluster C show expression patterns opposite to cluster D, increasing with specific growth and photosynthesis rates. Key genes include glnA (glutamine synthetase) and gltB (glutamate synthase large subunit), suggesting potential direct exchange of glutamine and glutamate from T. elongatus. M. ruber transcripts encoding methionine biosynthesis proteins (metHX) also decrease with irradiance and group into cluster C, indicating salvage of cyanobacterium-derived methionine. Vitamin B12 uptake/scavenging genes (btuCD) decrease with irradiance, correlating with decreased T. elongatus vitamin B12 biosynthesis genes in cluster C.
**Temporal pattern:** Increase with specific growth and photosynthesis rates, opposite to cluster D; decrease with irradiance.
**Confidence notes:** Strong evidence from coexpression and functional annotation supports metabolic coupling in nitrogen and vitamin metabolism.
**Sources:** Figure 6
**Quotes:**
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA and the large subunit of glutamate synthase (gltB).
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis.

### Cluster D | up with irradiance, down with oxygen | high
**Name:** Meiothermus ruber cluster D (irradiance-responsive)
**Enrichment:** Carbon metabolism, oxidative stress response (p=0.05, sig=True)
**Functional:** Cluster D contains genes that increase with irradiance and/or decrease with oxygen tension. M. ruber genes in this cluster include acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), xylose isomerase (xylA), xylulokinase (xylB), ABC-type multisugar uptake system, and branched-chain amino acid uptake (livGH). These genes are involved in uptake and metabolism of organic carbon compounds derived from T. elongatus, supporting heterotrophic growth. M. ruber genes encoding glutamate dehydrogenase (gdhA) and amino acid uptake systems also group here. Additionally, M. ruber oxyR (H2O2-responsive transcriptional regulator) groups into clusters D and G, adjacent to manganese catalase gene, indicating oxidative stress response.
**Temporal pattern:** Increase with irradiance and decrease with oxygen tension, corresponding to increased specific growth and photosynthesis rates.
**Confidence notes:** Cluster D genes show coordinated expression with T. elongatus genes involved in carbon metabolism and oxidative stress, supporting metabolic coupling.
**Sources:** Figure 6
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G).
- [Page 7] In conjunction with T. elongatus, M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G.
- [Page 8] M. ruber contains an H2O2-responsive transcriptional regulator (oxyR), which grouped into clusters D and G and is located adjacent to a putative manganese catalase.

## Bernstein 2017 / mruber_oxygen_clusters

### Cluster E | peak at intermediate pO2 | high
**Name:** Meiothermus ruber cluster E (pO2-responsive)
**Enrichment:** Nitrogen metabolism and vitamin B12-related functions (p=0.05, sig=True)
**Functional:** Cluster E contains 328 M. ruber genes that show transcriptional patterns trending with increasing pO2 treatments, analogous to cluster A for irradiance. Genes include those involved in nitrogen acquisition such as glnA and gltB, and vitamin B12-dependent methylmalonyl-CoA mutase subunits (mcmA1 and mcmA2).
**Temporal pattern:** Genes in cluster E show a tent-shaped pattern with maximum relative mRNA abundances at the midpoint pO2 treatment, decreasing at higher and lower pO2 levels.
**Confidence notes:** Based on coexpression patterns and functional enrichment analysis with significant p-values.
**Sources:** Figure 6, Table S1, Table S2
**Quotes:**
- [Page 6-7] Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint Ii (1,190 mol photons m2 s1). ... Clusters E to H contain transcripts that trend with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles that are analogous to clusters A to D.
- [Page 8] M. ruber contains the vitamin B12-dependent methylmalonyl-CoA mutase subunits (mcmA1 and -A2; SY28_RS14080 and SY28_RS02560) which grouped into clusters A and H.
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA (SY28_RS02395) and the large subunit of glutamate synthase (gltB; SY28_RS09480).

### Cluster F | down with increasing pO2 | high
**Name:** Meiothermus ruber cluster F (pO2-responsive)
**Enrichment:** Amino acid metabolism and transport (p=0.05, sig=True)
**Functional:** Cluster F contains 309 M. ruber genes with transcriptional profiles that decrease with increasing pO2, analogous to cluster B for irradiance. Genes include those involved in amino acid uptake and metabolism, such as amino acid uptake system genes and glutamate dehydrogenase (gdhA).
**Temporal pattern:** Genes in cluster F show a decrease in relative mRNA abundance with increasing pO2.
**Confidence notes:** Functional enrichment significant; gene examples support amino acid metabolism role.
**Sources:** Figure 6, Table S1, Table S2
**Quotes:**
- [Page 8] Several M. ruber nitrogen-associated genes also grouped into clusters D and/or G, including a glutamate dehydrogenase gene (gdhA; SY28_RS07480 and SY28_RS07475) and amino acid uptake system genes (SY28_RS09865 and SY28_RS06725).
- [Page 8] Genes involved in cysteine biosynthesis shared common transcriptional patterning between species, including cysteine synthases (cysK; tlr0504 and SY28_RS03805) and serine O-acetyltransferases (cysE; tlr0851 and SY28_RS05065), which cogrouped into clusters D and A, respectively.

### Cluster G | down with increasing pO2, up with growth rate | high
**Name:** Meiothermus ruber cluster G (pO2-responsive)
**Enrichment:** Carbon metabolism and transport (p=0.05, sig=True)
**Functional:** Cluster G contains 468 M. ruber genes that increase with specific growth and photosynthesis rates and decrease with increasing pO2. Genes include those involved in carbon uptake and metabolism such as acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA and gtsAB), and branched-chain amino acid uptake (livGH).
**Temporal pattern:** Genes in cluster G decrease in relative mRNA abundance with increasing pO2 and increase with specific growth and photosynthesis rates.
**Confidence notes:** Supported by coexpression with T. elongatus genes and functional enrichment analysis.
**Sources:** Figure 6, Table S1, Table S2
**Quotes:**
- [Page 7-8] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). ... In conjunction with T. elongatus, M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. Notable examples included genes encoding putative acetyl-CoA synthetase (acs; SY28_RS04760 and SY28_RS00910), cytochrome c oxidase (coxB; SY28_RS11630), monosaccharide uptake systems (frcA and gtsAB; SY28_RS03690, SY28_RS04315, and SY28_RS04260), xylose isomerase (xylA; SY28_RS02810), xylulokinase (xylB; SY28_RS03685), an ABC-type multisugar uptake system (SY28_RS06965), and branched-chain amino acid uptake (livGH; SY28_RS00360 and SY28_RS00370).
- [Page 9] M. ruber genes associated with electron transfer processes that are potentiators of ROS grouped into cluster G, which decreased with increasing pO2 treatments and increased with specific growth and photosynthesis rates. These include subunits for an electron transfer flavoprotein (fixAB; SY28_RS10190 and SY28_RS10195), NADH-dehydrogenase (SY28_RS04385), and principal components of the NADH-quinone oxidoreductase (nuoDFGHIJKN).

### Cluster H | up with increasing pO2 | high
**Name:** Meiothermus ruber cluster H (pO2-responsive)
**Enrichment:** Oxidative stress response and vitamin metabolism (p=0.05, sig=True)
**Functional:** Cluster H contains 387 M. ruber genes that increase with increasing pO2. Genes include those involved in ROS detoxification such as peroxidase (bcp) and superoxide dismutase (sod2), as well as vitamin B12-dependent methionine synthesis (metH).
**Temporal pattern:** Genes in cluster H increase in relative mRNA abundance with increasing pO2.
**Confidence notes:** Supported by coexpression with T. elongatus genes and functional enrichment analysis.
**Sources:** Figure 6, Table S1, Table S2
**Quotes:**
- [Page 9] M. ruber peroxidase (bcp; SY28_RS05545, SY28_RS06010, and SY28_RS06015) and superoxide dismutase (sod2; SY28_RS13295) genes responded differently to Ii treatments than did oxyR and related cyanobacterial profiles and were grouped into clusters A and C. These genes generally increased with pO2 and grouped with the T. elongatus genes into cluster H (increased with pO2 and decreasing μ).
- [Page 9] M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. ... M. ruber genes associated with electron transfer processes that are potentiators of ROS grouped into cluster G, which decreased with increasing pO2 treatments and increased with specific growth and photosynthesis rates.

## Biller 2018 / mit1002_periodicity

### Cluster coculture_LD | periodic in L:D only | high
**Name:** Alteromonas macleodii MIT1002 cluster coculture_LD (L:D periodic)
**Enrichment:** Photosynthesis and metabolic pathways (p=0.01, sig=True)
**Functional:** Enriched for transcripts exhibiting 24-h periodicity in coculture under diel light:dark conditions, including genes associated with photosynthesis, Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosystem components.
**Temporal pattern:** Genes show 24-h periodic oscillations in coculture during the normal diel light:dark cycle but lose periodicity under extended darkness when grown alone. Periodicity is maintained longer in coculture than in axenic culture.
**Confidence notes:** High confidence based on transcriptome periodicity analysis and pathway enrichment tests.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 1] More Prochlorococcus transcripts exhibited 24-h periodic oscillations in coculture than in pure culture, both over the normal diel cycle and after the shift to extended darkness.
- [Page 11] The largest group of Prochlorococcus transcripts (42% of all protein-encoding genes) showed 24-h periodicity in both axenic and cocultures under diel L:D conditions but did not continue to oscillate under extended darkness.
- [Page 11] Within these general trends in periodic behavior, we found that subsets of the Prochlorococcus transcriptome exhibited periodic oscillations under different combinations of culture conditions.

### Cluster coculture_LD+coculture_darkness | periodic across all conditions | high
**Name:** Alteromonas macleodii MIT1002 cluster coculture_LD+coculture_darkness (periodic across all conditions)
**Enrichment:** Metabolic and photosynthetic pathways (p=0.01, sig=True)
**Functional:** Includes transcripts that maintain 24-h periodicity in coculture both under diel light:dark cycles and extended darkness, associated with metabolic pathways such as Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis.
**Temporal pattern:** Genes exhibit 24-h periodic oscillations in coculture under both diel light:dark and extended darkness conditions, indicating sustained metabolic activity in the dark when in coculture.
**Confidence notes:** High confidence from transcriptomic periodicity data and pathway enrichment analysis.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 11] Looking more closely at the set of transcripts that continued oscillating under extended darkness in cocultured but not in axenic cells, we find transcripts associated with a variety of metabolic pathways, including the Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis.
- [Page 11] This further indicates that the presence of the heterotroph may help Prochlorococcus survive the stress of extended darkness by contributing to the maintenance of regular metabolic functionality.

### Cluster coculture_darkness | not periodic | high
**Name:** Alteromonas macleodii MIT1002 cluster coculture_darkness (not periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** N/A
**Sources:** Table 3
**Quotes:**
- [Page 11] Only 6% of transcripts exhibited 24-h periodicity through the first 13 h of extended darkness in axenic cultures, while 25% did so in cocultures; periodicity essentially vanished in Alteromonas under extended darkness.

## Biller 2018 / natl2a_darkness_survival

### Cluster darkness_axenic+darkness_coculture | constitutive during extended darkness in both axenic and coculture conditions, with partial retention of periodicity in coculture | high
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+darkness_coculture (darkness)
**Enrichment:** Photosynthesis, Fermentation, Organic compound degradation (p=0.001, sig=True)
**Functional:** Genes in this cluster are expressed during extended darkness in both axenic and coculture conditions. They include photosystem I and II components (psbA, psbB, psaA, psaB) and genes involved in light reactions of photosynthesis, adenylate kinase, and respiratory functions. Transcripts for DNA repair genes (recA, ruvB, mutS, radA) were enriched only in axenic cells, not in cocultures. This cluster shows depletion of transcripts for biosynthetic pathways, NAD metabolism, ATP synthase subunits, glycogen degradation (glgP), oxidative pentose phosphate pathway, and late glycolysis steps in axenic cultures. Enrichment of fermentation-associated transcripts and pyruvate kinase (pykF) was observed in axenic cultures. Cocultured cells showed enrichment in transcripts encoding organic compound degradation/salvage pathways and the Entner-Doudoroff pathway enzyme Eda, suggesting mixotrophic metabolism. Transcripts for transporters were depleted in axenic but not cocultured cells, indicating higher uptake potential in cocultures. The stringent response was activated in axenic but not cocultured cells. Overall, this cluster reflects metabolic shifts and stress responses during extended darkness, with heterotroph presence mitigating some stress effects.
**Temporal pattern:** Genes in this cluster have detectable transcripts during late extended darkness (72-144h) in both axenic and coculture conditions. Transcript abundance changes occur primarily within the first 5 hours past the missed sunrise and continue through extended darkness. Periodic oscillations in transcript abundance are lost in axenic cultures but partially retained in cocultures during extended darkness.
**Confidence notes:** High confidence based on multiple transcriptomic analyses and pathway enrichment tests.
**Sources:** Figure 1A-D, Figure 3, Table 1, Table 2, Table S3
**Quotes:**
- [Page 6] The bulk of the transcriptional response was consistent with the cessation of growth in the dark... Cells were depleted in transcripts encoding a number of biosynthetic pathways, NAD metabolism genes, and ATP synthase subunits relative to cells that experienced sunrise on schedule... The activity that continued appeared to be focused largely on recovering energy and/or nutrients... Some activity also appeared directed toward collecting and utilizing photons, given that the transcriptomes were enriched for genes involved in the light reactions of photosynthesis...
- [Pages 8-9] Transcripts for DNA repair genes such as recA, ruvB, mutS, and radA were enriched only in the axenic cells, not in cocultures, under extended darkness... cocultured, but not axenic, Prochlorococcus cells were enriched in transcripts encoding organic compound degradation/salvage pathways... transcripts for the key enzyme in the ED pathway, Eda, were enriched in coculture and consistently depleted in axenic cultures...
- [Page 9] The stringent response was activated in axenic cultures but not similarly in cocultures... cocultures did not exhibit any enrichment for hpf, a ribosome hibernation-promoting factor...

### Cluster darkness_axenic+unique_axenic | N/A | low
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+unique_axenic (darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data provided for this cluster in the paper.

### Cluster darkness_coculture+unique_coculture | constitutive during extended darkness in coculture, with partial retention of periodicity | high
**Name:** Prochlorococcus NATL2A cluster darkness_coculture+unique_coculture (darkness)
**Enrichment:** Organic compound degradation, Amino acid degradation, Terpenoid and tetrapyrrole biosynthesis (p=0.006, sig=True)
**Functional:** Genes in this cluster are expressed during extended darkness uniquely or predominantly in coculture with Alteromonas. They include transcripts encoding organic compound degradation/salvage pathways, such as amino acid degradation, adenine and adenosine salvage (apt), lipoate salvage (lplA), purine nucleotide degradation (truB), terpenoid biosynthesis (dxs, ispG, ispE, lytB), and tetrapyrrole biosynthesis (hemB, hemD). This cluster shows enrichment of transcripts for the Entner-Doudoroff pathway enzyme Eda, suggesting enhanced glycolytic metabolism of organic molecules. Transcripts for transporters were relatively maintained compared to axenic cultures, indicating higher uptake potential. The cluster is associated with maintenance of biosynthetic potential and metabolic activity during extended darkness in coculture, supporting mixotrophic metabolism.
**Temporal pattern:** Genes in this cluster have detectable transcripts during late extended darkness (72-144h) uniquely or predominantly in coculture conditions. Transcript abundance changes occur primarily within the first 5 hours past the missed sunrise and continue through extended darkness. Periodic oscillations in transcript abundance are partially retained in cocultures during extended darkness.
**Confidence notes:** High confidence based on transcriptomic enrichment and pathway analysis.
**Sources:** Figure 3, Table 2, Table S3
**Quotes:**
- [Pages 8-9] Cocultured, but not axenic, Prochlorococcus cells were enriched in transcripts encoding organic compound degradation/salvage pathways, such as those for amino acids... transcripts for the key enzyme in the ED pathway, Eda, were enriched in coculture and consistently depleted in axenic cultures...
- [Page 9] Transcripts involved in some amino acid and nucleotide biosynthetic pathways were depleted under extended darkness in the axenic cells, but not in the cocultures... The enrichment of salvage pathways in the cocultures indicates that the unknown compound(s) we propose were supplied by Alteromonas could also have been required biosynthetic intermediates...
- [Page 9] Transcripts of at least 10 transporters were relatively depleted in axenic but not cocultured Prochlorococcus transcriptomes, raising the possibility that cocultures could have maintained a relatively higher uptake potential...

## Biller 2018 / natl2a_periodicity

### Cluster axenic_LD | periodic in axenic L:D only | high
**Name:** NATL2A cluster axenic_LD (axenic L:D periodic)
**Enrichment:** Photosynthesis (p=0.021, sig=True)
**Functional:** Genes in this cluster exhibit 24-h periodicity only in axenic cultures under the diel light:dark cycle. This cluster includes genes involved in photosystem components and light-harvesting complexes that show oscillations tied to the light cycle in axenic conditions.
**Temporal pattern:** Genes show 24-h periodicity exclusively in axenic L:D conditions and lose periodicity in coculture and extended darkness conditions.
**Confidence notes:** High confidence based on transcriptome periodicity analysis and enrichment tests.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 11] The largest group of Prochlorococcus transcripts (42% of all protein-encoding genes) showed 24-h periodicity in both axenic and cocultures under diel L:D conditions but did not continue to oscillate under extended darkness.

### Cluster axenic_LD+axenic_darkness+coculture_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+axenic_darkness+coculture_darkness (periodic in axenic L:D and darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data for this cluster in the paper.

### Cluster axenic_LD+coculture_LD | periodic in L:D only | high
**Name:** NATL2A cluster axenic_LD+coculture_LD (periodic in both L:D conditions)
**Enrichment:** Metabolism and photosynthesis (p=0.001, sig=True)
**Functional:** This cluster contains the largest number of genes showing 24-h periodicity in both axenic and cocultured Prochlorococcus under diel L:D conditions. Genes include those involved in core metabolic pathways such as Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis.
**Temporal pattern:** Genes oscillate with 24-h periodicity in both axenic and cocultured L:D conditions but lose periodicity under extended darkness.
**Confidence notes:** High confidence due to large gene set and consistent periodicity patterns.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 11] The largest group of Prochlorococcus transcripts (42% of all protein-encoding genes) showed 24-h periodicity in both axenic and cocultures under diel L:D conditions but did not continue to oscillate under extended darkness.
- [Page 11] Within these general trends in periodic behavior, we found that subsets of the Prochlorococcus transcriptome exhibited periodic oscillations under different combinations of culture conditions.

### Cluster axenic_LD+coculture_LD+axenic_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+coculture_LD+axenic_darkness (periodic in axenic darkness and L:D)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No explicit mention or data for this cluster.

### Cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness | periodic across all conditions | medium
**Name:** NATL2A cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness (periodic across all conditions)
**Enrichment:** Core metabolism and photosynthesis (p=0.005, sig=True)
**Functional:** Genes in this cluster maintain 24-h periodicity across all four conditions, including both L:D cycles and extended darkness in axenic and coculture. These include genes involved in core metabolic and photosynthetic functions, indicating sustained rhythmicity even in darkness when cocultured.
**Temporal pattern:** Genes oscillate with 24-h periodicity in all conditions, including extended darkness in both axenic and coculture.
**Confidence notes:** Moderate to high confidence based on periodicity analysis and gene function.
**Assessment notes:** Some uncertainty about exact gene membership and functional enrichment.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 11] Other transcripts maintained periodic oscillations in different sets of culture conditions, including those that oscillated under all conditions.
- [Page 11] Looking more closely at the set of transcripts that continued oscillating under extended darkness in cocultured but not in axenic cells, we find transcripts associated with a variety of metabolic pathways, including the Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis.

### Cluster axenic_LD+coculture_LD+coculture_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+coculture_LD+coculture_darkness (periodic in coculture darkness and L:D)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data for this cluster.

### Cluster axenic_LD+coculture_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+coculture_darkness (periodic in axenic L:D and coculture darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No explicit mention or data for this cluster.

### Cluster coculture_LD | periodic in coculture L:D only | medium
**Name:** NATL2A cluster coculture_LD (coculture L:D periodic)
**Enrichment:** Regulatory and metabolic genes (p=0.02, sig=True)
**Functional:** Genes showing 24-h periodicity only in cocultured Prochlorococcus under diel L:D conditions. This cluster includes genes that maintain rhythmic expression in coculture but not in axenic cultures, suggesting heterotroph influence on timing mechanisms.
**Temporal pattern:** Genes oscillate with 24-h periodicity exclusively in coculture L:D conditions and lose periodicity in axenic cultures and extended darkness.
**Confidence notes:** Moderate confidence based on periodicity data and coculture effects.
**Sources:** Figure 4A, Figure 4B
**Quotes:**
- [Page 11] More Prochlorococcus transcripts retained their periodicity in cocultures versus axenic cultures, both under extended darkness and even during the normal diel L:D cycle.

### Cluster coculture_LD+axenic_darkness | N/A | low
**Name:** NATL2A cluster coculture_LD+axenic_darkness (periodic in coculture L:D and axenic darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific data or discussion for this cluster.

### Cluster coculture_LD+coculture_darkness | N/A | low
**Name:** NATL2A cluster coculture_LD+coculture_darkness (periodic in coculture L:D and darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No explicit mention or data for this cluster.

### Cluster coculture_darkness | N/A | low
**Name:** NATL2A cluster coculture_darkness (coculture darkness periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion or data for this cluster.

### Cluster not_periodic | not periodic | high
**Name:** NATL2A cluster not_periodic (not periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Genes that do not exhibit significant 24-h periodicity under any of the four conditions. These genes likely represent constitutively expressed or non-rhythmic functions in Prochlorococcus.
**Temporal pattern:** Genes lack 24-h periodic oscillations across all tested conditions, including axenic and coculture L:D and extended darkness.
**Confidence notes:** High confidence as these genes were identified as non-periodic by the analysis.
**Sources:** Figure 4A
**Quotes:**
- [Page 11] Other transcripts did not exhibit periodicity under any conditions.

## Coe 2024 / supp_table_3_darktolerant_clusters

### Cluster 1 | N/A | low
**Name:** Prochlorococcus cluster 1 (183 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text; only clustering method and general periodicity reported.
**Sources:** Supplemental Table 3
**Quotes:**
- [Supplemental Text 1] Genes in parental and dark-tolerant cells were clustered into 15 clusters. If a gene in both parental and dark-tolerant cells was found in the same cluster, the pattern was marked “true”, and if not, as “false”.

### Cluster 10 | N/A | low
**Name:** Prochlorococcus cluster 10 (173 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 11 | N/A | low
**Name:** Prochlorococcus cluster 11 (93 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 12 | N/A | low
**Name:** Prochlorococcus cluster 12 (94 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 13 | N/A | low
**Name:** Prochlorococcus cluster 13 (140 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 14 | N/A | low
**Name:** Prochlorococcus cluster 14 (169 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 15 | N/A | low
**Name:** Prochlorococcus cluster 15 (130 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 2 | N/A | low
**Name:** Prochlorococcus cluster 2 (133 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 3 | N/A | low
**Name:** Prochlorococcus cluster 3 (130 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 4 | N/A | low
**Name:** Prochlorococcus cluster 4 (169 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 5 | N/A | low
**Name:** Prochlorococcus cluster 5 (171 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 6 | N/A | low
**Name:** Prochlorococcus cluster 6 (181 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 7 | N/A | low
**Name:** Prochlorococcus cluster 7 (118 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 8 | N/A | low
**Name:** Prochlorococcus cluster 8 (83 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

### Cluster 9 | N/A | low
**Name:** Prochlorococcus cluster 9 (114 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of individual clusters in the text.
**Sources:** Supplemental Table 3

## Coe 2024 / supp_table_3_parental_clusters

### Cluster 1 | N/A | low
**Name:** Prochlorococcus cluster 1 (97 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 10 | N/A | low
**Name:** Prochlorococcus cluster 10 (70 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 11 | N/A | low
**Name:** Prochlorococcus cluster 11 (205 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 12 | N/A | low
**Name:** Prochlorococcus cluster 12 (144 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 13 | N/A | low
**Name:** Prochlorococcus cluster 13 (165 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 14 | N/A | low
**Name:** Prochlorococcus cluster 14 (89 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 15 | N/A | low
**Name:** Prochlorococcus cluster 15 (132 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 2 | N/A | low
**Name:** Prochlorococcus cluster 2 (110 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 3 | N/A | low
**Name:** Prochlorococcus cluster 3 (245 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 4 | N/A | low
**Name:** Prochlorococcus cluster 4 (100 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 5 | N/A | low
**Name:** Prochlorococcus cluster 5 (65 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 6 | N/A | low
**Name:** Prochlorococcus cluster 6 (62 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 7 | N/A | low
**Name:** Prochlorococcus cluster 7 (202 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 8 | N/A | low
**Name:** Prochlorococcus cluster 8 (174 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

### Cluster 9 | N/A | low
**Name:** Prochlorococcus cluster 9 (221 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific description of this cluster in the paper or supplement.
**Sources:** Supplemental Table 3

## Lindell 2007 / med4_phage_transcription_groups

### Cluster 1 | early transient | high
**Name:** MED4 cluster 1 (early transient peak)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** This cluster is transiently upregulated immediately after infection and consists of high-light-inducible stress response genes (hli), carbon metabolism genes (rbcLS), transcription genes (rpoC2, rpoD), and ribosomal protein genes (rpl5, rpl6, rps8, rps11, rps17).
**Temporal pattern:** Genes in this cluster are transiently upregulated immediately after infection, showing an early transient peak in expression around 1 hour post infection.
**Sources:** Figure 3
**Quotes:**
- [Page 3] The first was transiently upregulated immediately after infection and consists of high-light-inducible stress response (hli), carbon metabolism (rbcLS), transcription (rpoC2, rpoD) and ribosome (rpl5, rpl6, rps8, rps11, rps17) genes.

### Cluster 2 | late sustained | high
**Name:** MED4 cluster 2 (late sustained induction)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** This cluster is upregulated starting 2 hours after infection and includes genes involved in RNA degradation and modification (rne, rnhB, dus, sun), protein turnover (clpS, AAA ATPase family gene), stress responses (umuD, phoH), and genes of unknown function. Many of these genes are located in hypervariable genome islands and have homologues in phage genomes, suggesting co-evolutionary significance.
**Temporal pattern:** Genes in this cluster begin to be upregulated about 2 hours after infection and show sustained induction through the infection time course.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Transcripts of the second group appeared 2 h after infection and included genes involved in RNA degradation and modification (rne, rnhB, dus and sun), protein turnover (clpS, and an AAA ATPase family gene), stress responses (umuD and phoH), and those of unknown function.
- [Page 3] Perhaps the most compelling evidence that upregulated host genes are part of the co-evolutionary process in this system is that 34% of them (more than would occur by chance P , 0.001) are found in hypervariable host genome islands (Supplementary Table 3), which are thought to be mobilized by phages.

## Tolonen 2006 / med4_kmeans_nstarvation

### Cluster 1 | early transient | high
**Name:** MED4 cluster 1 (5 genes)
**Enrichment:** transport and binding (p=0.01, sig=True)
**Functional:** Most rapidly and highly upregulated genes containing nitrogen transport genes such as urtA, cynA, and the nitrite permease in MIT9313. Includes highly upregulated unknown genes with strong NtcA binding sites such as PMM0958 and PMM1462.
**Temporal pattern:** Genes respond rapidly and highly within the first 6 hours of nitrogen starvation, with peak upregulation at 6 hours and some decline thereafter.
**Confidence notes:** High confidence due to direct quotes and strong NtcA binding site evidence.
**Sources:** Figure 3, Table I, Figure 4E
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 urtA, MED4 cynA, and the MIT9313 nitrite permease.
- [Page 7] MED4 PMM0958 was the most upregulated gene at all time points and has the top-ranking NtcA binding site in the genome.

### Cluster 2 | early transient | medium
**Name:** MED4 cluster 2 (16 genes)
**Enrichment:** transport and binding (p=0.12, sig=False)
**Functional:** Contains a subset of hli genes (hli10, hli21, hli22) that respond rapidly and highly to nitrogen starvation, and two sigma factors that are upregulated.
**Temporal pattern:** Rapid and high response early in nitrogen starvation time course.
**Confidence notes:** Moderate confidence; functional enrichment not statistically significant but gene content is notable.
**Assessment notes:** Functional enrichment not significant but gene content supports description.
**Sources:** Figure 3, Figure 4B, Figure 4E
**Quotes:**
- [Page 4] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1).
- [Page 4] MED4 cluster 5 contains two sigma factors that are upregulated during N stress.

### Cluster 3 | late sustained | medium
**Name:** MED4 cluster 3 (19 genes)
**Enrichment:** regulation (p=0.07, sig=False)
**Functional:** Contains a second subset of hli genes that respond later and to a lesser degree to nitrogen starvation, and an upregulated cluster containing two sigma factors.
**Temporal pattern:** Later and less pronounced response compared to cluster 2.
**Confidence notes:** Moderate confidence based on gene content and expression patterns.
**Assessment notes:** Functional enrichment borderline significant.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 4] Those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).
- [Page 4] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.

### Cluster 4 | early sustained repression | high
**Name:** MED4 cluster 4 (86 genes)
**Enrichment:** photosynthesis and respiration (p=2.6e-09, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes including photosystem I genes (psaBDEIJKLM). Contains genes for carbon fixation (rbcLS) and ATP synthase in cluster 9.
**Temporal pattern:** Rapid repression beginning early in nitrogen starvation and sustained through the time course.
**Confidence notes:** High confidence due to strong enrichment and gene content.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MED4 cluster 8 contains numerous genes for photosystem I (psaBDEIJKLM).
- [Page 4] MED4 cluster 9 contains ATP synthase subunits and the carbon fixation genes (rbcLS).

### Cluster 5 | early sustained repression | high
**Name:** MED4 cluster 5 (73 genes)
**Enrichment:** translation (p=0.001, sig=True)
**Functional:** Enriched for translation genes, including ribosomal proteins and genes involved in protein synthesis.
**Temporal pattern:** Strong repression early in nitrogen starvation, sustained over time.
**Confidence notes:** High confidence due to strong enrichment and gene content.
**Assessment notes:** Cluster 7 in figure corresponds to translation; cluster 5 in prompt likely a numbering difference, assumed cluster 7 here.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MED4 cluster 7 is enriched for translation genes.

### Cluster 6 | early transient | high
**Name:** MIT9313 cluster 1 (42 genes)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Most rapidly and highly upregulated genes including nitrogen transport genes such as urtA, nitrite permease, and hli genes hliS and hli7.
**Temporal pattern:** Rapid and high upregulation early in nitrogen starvation, peaking around 6-12 hours.
**Confidence notes:** High confidence due to direct quotes and strong NtcA binding site evidence.
**Sources:** Figure 3, Table I, Figure 4E
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MIT9313 urtA and the nitrite permease.
- [Page 7] In MIT9313, hliS and hli7 genes were by far the most upregulated of all genes in the genome (approximately 70-fold).

### Cluster 7 | late sustained | medium
**Name:** MIT9313 cluster 2 (37 genes)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** Contains regulatory genes including two sigma factors that are upregulated and a subset of hli genes that respond later and to a lesser degree.
**Temporal pattern:** Later and sustained upregulation during nitrogen starvation.
**Confidence notes:** Moderate to high confidence based on gene content and enrichment.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 4] MIT9313 cluster 2 contains two sigma factors and a subset of hli genes that respond later and to a lesser degree.

### Cluster 8 | late repression | high
**Name:** MIT9313 cluster 3 (37 genes)
**Enrichment:** photosynthesis and respiration (p=0.00049, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes including photosystem I and II genes and phycoerythrin gene cpeB.
**Temporal pattern:** Repressed only at the final time point (48h) indicating a shutdown response rather than specific nitrogen stress response.
**Confidence notes:** High confidence due to strong enrichment and gene content.
**Assessment notes:** Cluster 6 in figure corresponds to photosynthesis; cluster 3 in prompt likely a numbering difference, assumed cluster 6 here.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 6, the only cluster enriched for Photosynthesis and Respiration in this strain, contains genes for diverse aspects of photosystem I and II along with the phycoerythrin gene, cpeB.
- [Page 4] MIT9313 cluster 4 contains a number of genes that were unchanged until being repressed only at the final time point.

### Cluster 9 | early sustained repression | high
**Name:** MIT9313 cluster 4 (8 genes)
**Enrichment:** translation (p=7.9e-10, sig=True)
**Functional:** Enriched for translation genes involved in protein synthesis and ribosomal proteins.
**Temporal pattern:** Strong repression early in nitrogen starvation, sustained over time.
**Confidence notes:** High confidence due to strong enrichment and gene content.
**Assessment notes:** Cluster 7 in figure corresponds to translation; cluster 4 in prompt likely a numbering difference, assumed cluster 7 here.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 7 is enriched for translation genes.

## Tolonen 2006 / mit9313_kmeans_nstarvation

### Cluster 1 | early transient | high
**Name:** MIT9313 cluster 1 (transport and binding)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Most rapidly and highly upregulated genes including nitrogen transport genes such as urtA, nitrite permease, and hli genes hliS and hli7. Enriched for transport and binding (p=0.04).
**Temporal pattern:** Genes respond rapidly and highly within the first hours of nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA
- [Page 4] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA, MED4 cynA, and the MIT9313 nitrite permease.

### Cluster 2 | early transient | high
**Name:** MIT9313 cluster 2 (regulation)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** Contains two upregulated sigma factors and a subset of hli genes (hliS, hli7). Enriched for regulation (p=0.03).
**Temporal pattern:** Genes respond rapidly and highly early in nitrogen starvation.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 4] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.
- [Page 6] Two out of seven MIT9313 sigma factors (PMT2246 and PMT0346) were upregulated and may play a role in the upregulation of gene expression during N stress in Prochlorococcus.

### Cluster 3 | N/A | low
**Name:** MIT9313 cluster 3 (regulation)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Paper does not discuss this cluster.
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 4 | late sustained repression | high
**Name:** MIT9313 cluster 4 (transport and binding)
**Enrichment:** transport and binding (p=0.39, sig=False)
**Functional:** Contains genes repressed only at the final time point, possibly representing general transcriptional shutdown rather than specific N stress response. Includes genes linking nitrogen and carbon metabolism such as glnB, icd, acnB, and rbcLS.
**Temporal pattern:** Genes unchanged until repression at the final time point (48h) during severe starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 4 contains a number of genes that were unchanged until being repressed only at the final time point.
- [Page 4] The genes in MIT9313 cluster 4 may thus represent those genes that are repressed as part of a general shutdown in transcription, rather than a specific N stress response.

### Cluster 5 | N/A | low
**Name:** MIT9313 cluster 5 (fatty acid and phospholipid)
**Enrichment:**  (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** Paper does not discuss this cluster.
**Assessment notes:** Paper does not discuss this cluster.

### Cluster 6 | gradual repression | high
**Name:** MIT9313 cluster 6 (photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=0.00049, sig=True)
**Functional:** Strongly enriched for photosynthesis and respiration (p=0.00049). Contains genes for photosystem I and II and phycoerythrin gene cpeB. Repressed during nitrogen starvation.
**Temporal pattern:** Genes gradually repressed during nitrogen starvation, especially at later time points.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 6, the only cluster enriched for Photosynthesis and Respiration in this strain, contains genes for diverse aspects of photosystem I and II along with the phycoerythrin gene, cpeB.

### Cluster 7 | sustained repression | high
**Name:** MIT9313 cluster 7 (translation)
**Enrichment:** translation (p=7.9e-10, sig=True)
**Functional:** Strongly enriched for translation (p=7.9e-10). Contains ribosomal and translation-related genes. Repressed during nitrogen starvation.
**Temporal pattern:** Genes strongly repressed during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Highly downregulated clusters show strong enrichment for specific functional categories, namely Translation (MIT9313 cluster 7).

## Wang 2014 / med4_expression_level

### Cluster -- | N/A | high
**Name:** Prochlorococcus MED4 cluster -- (--)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** No specific discussion or data provided for this cluster in the paper.
**Assessment notes:** No information available for this cluster.

### Cluster HEG | high across all | high
**Name:** Prochlorococcus MED4 cluster HEG (heg)
**Enrichment:** Energy production and conversion, Translation, Protein folding (p=0.001, sig=True)
**Functional:** Highly expressed genes (HEG) are enriched in the core genome and include genes involved in energy production and conversion, translation and ribosomal structure, and protein modification, folding and turnover. Examples include ribosomal proteins and photosynthetic apparatus genes.
**Temporal pattern:** HEG genes are consistently highly expressed across all tested growth conditions.
**Confidence notes:** Strong enrichment in core genome and functional categories with significant p-values.
**Sources:** Figure 3, Figure 4c, Figure 5b, Figure 6b
**Quotes:**
- [Page 5] Among these core HEG genes, several functional categories were more prominent than others. These included the “C” (energy production and conversion), “J” (translation and ribosomal structure), and “O” (protein modification, folding and turnover) categories (Figure 4c).
- [Page 4] HEG had a significantly lower nonsynonymous substitution rate (Ka) than MEG or LEG (Kruskal-Wallis Test, two-tailed P < 0.001; Figure 3a), indicating a strong negative correlation between gene expression level and evolutionary rate.

### Cluster LEG | low across all | medium
**Name:** Prochlorococcus MED4 cluster LEG (leg)
**Enrichment:** General function (p=0.023, sig=True)
**Functional:** Lowly expressed genes (LEG) are less enriched in the core genome and show slight enrichment in the general function category (COG category R).
**Temporal pattern:** LEG genes are consistently lowly expressed across the tested conditions.
**Confidence notes:** LEG genes show slight enrichment in general function category with borderline significance.
**Assessment notes:** Limited detailed functional description; enrichment only slight.
**Sources:** Figure 3, Figure 4c
**Quotes:**
- [Page 5] Additionally, category “R” (general function) was slightly enriched in both LEG and NEG (P = 0.023 and 0.055; data not shown).
- [Page 4] LEG had a significantly higher nonsynonymous substitution rate than HEG (Kruskal-Wallis Test, two-tailed P < 0.001; Figure 3a).

### Cluster MEG | moderate across all | high
**Name:** Prochlorococcus MED4 cluster MEG (meg)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Moderately expressed genes (MEG) are enriched in the core genome and have a significantly higher rate of essential genes (DEG-hit).
**Temporal pattern:** MEG genes are moderately expressed consistently across conditions.
**Confidence notes:** MEG subclass has higher essential gene rate but no significant difference in expression level between DEG-hit and DEG-miss genes.
**Sources:** Figure 3, Figure 4b
**Quotes:**
- [Page 5] Although the MEG subclass had a significantly higher rate of DEG-hit genes (P < 0.001; Figure 4b), the mean expression level of the DEG-hit genes (mean RPKM = 602.62) was not significantly different from that of the DEG-miss genes (mean RPKM = 874.81; Student's t-test, two-tailed P = 0.084).
- [Page 4] Next, we compared the five gene expression subclasses of the core genome to that of the flexible genome. Our analysis clearly indicates that the genes in the HEG and MEG subclasses were more enriched in the core genome than in the flexible genome (17.7% > 11.5% and 26.8% > 15.3%, respectively; P < 0.001; Figure 3c).

### Cluster NEG | not expressed | high
**Name:** Prochlorococcus MED4 cluster NEG (neg)
**Enrichment:** General function (p=0.055, sig=False)
**Functional:** Non expressed genes (NEG) are less represented in the core genome compared to the flexible genome and show the lowest expression levels, including phage-related genes which were not tested under phage infection conditions.
**Temporal pattern:** Genes in the NEG cluster are consistently not expressed across the 10 growth conditions tested.
**Confidence notes:** Phage-related genes in NEG cluster have lowest expression; no phage infection conditions tested, limiting interpretation.
**Sources:** Figure 3, Figure 5b
**Quotes:**
- [Page 5] As expected, phage-related genes displayed the lowest expression levels in this study, as phage infection conditions were not tested.
- [Page 4] Conversely, the core genome had fewer NEG and VEG than the flexible genome (1.5% < 6.6% and 49.6% < 64.6%, respectively; P < 0.001; Figure 3c).

### Cluster VEG | variable across conditions | high
**Name:** Prochlorococcus MED4 cluster VEG (veg)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** Variably expressed genes (VEG) are more prevalent in the flexible genome and have higher nonsynonymous substitution rates than constantly expressed genes (CEG).
**Temporal pattern:** VEG genes show variable expression across the 10 growth conditions.
**Confidence notes:** VEG genes have higher evolutionary rates and are less conserved.
**Sources:** Figure 3
**Quotes:**
- [Page 4] CEG subclass had a lower Ka than VEG (Mann–Whitney U Test, two-tailed P < 0.001; Figure 3b), even when HEG were excluded from the CEG because of their bias with the lowest evolutionary rate among all expression subclasses (data not shown).
- [Page 4] Conversely, the core genome had fewer NEG and VEG than the flexible genome (1.5% < 6.6% and 49.6% < 64.6%, respectively; P < 0.001; Figure 3c).

## Zinser 2009 / med4_diel_clusters

### Cluster 1 | peaks near dawn | high
**Name:** Prochlorococcus cluster 1 (57 genes)
**Enrichment:** Photosystem I and II (p=1.5e-09, sig=True)
**Functional:** Enriched for photosystem I (p=1.5e-09) and photosystem II (p=2.5e-05) genes, including reaction center genes psbA and psbD.
**Temporal pattern:** Genes peak in expression near 8.3 h (about 2.3 hours after dawn) with 24h periodicity.
**Sources:** Figure 4A, Table 1, Table S4
**Quotes:**
- [Page 6] Expression of approximately half of photosystem (PS) II genes, including reaction center genes psbA and psbD, peak in abundance at mid-day
- [Page 6] Periodicity patterns of photosynthesis genes fell into 4 clusters (Table 1). Expression of approximately half of photosystem (PS) II genes, including reaction center genes psbA and psbD (encoding D1 and D2 respectively), as well as psbC (CP43) and psbF, co-varied with light intensity, with maxima at mid-day, and minima in the middle of the night (Figure 4A, Table S4)

### Cluster 10 | peaks just after dark | medium
**Name:** Prochlorococcus cluster 10 (77 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 0.3 h (just after dark) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 10.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle, with the exception of clusters 12 and 13, and 14 and 15, which had peak expression times less than one half hour apart.

### Cluster 11 | peaks early morning | medium
**Name:** Prochlorococcus cluster 11 (91 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 1.6 h (early morning) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 11.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle, with the exception of clusters 12 and 13, and 14 and 15, which had peak expression times less than one half hour apart.

### Cluster 12 | peaks early morning | medium
**Name:** Prochlorococcus cluster 12 (111 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 3.0 h (early morning) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 12.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle, with the exception of clusters 12 and 13, and 14 and 15, which had peak expression times less than one half hour apart.

### Cluster 13 | peaks early morning | high
**Name:** Prochlorococcus cluster 13 (110 genes)
**Enrichment:** Ribosomal proteins (p=3.7e-41, sig=True)
**Functional:** Strongly enriched for ribosomal proteins (45/53 genes, p=3.7e-41).
**Temporal pattern:** Genes peak at 3.4 h (early morning) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Page 3] Cluster 13 is significantly enriched for ribosomal proteins (45/53 genes, p=3.7e-41)

### Cluster 14 | peaks early morning | medium
**Name:** Prochlorococcus cluster 14 (22 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 4.6 h (early morning) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 14.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle, with the exception of clusters 12 and 13, and 14 and 15, which had peak expression times less than one half hour apart.

### Cluster 15 | peaks early morning | medium
**Name:** Prochlorococcus cluster 15 (107 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 4.7 h (early morning) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 15.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle, with the exception of clusters 12 and 13, and 14 and 15, which had peak expression times less than one half hour apart.

### Cluster 16 | peaks early morning | high
**Name:** Prochlorococcus cluster 16 (125 genes)
**Enrichment:** ATP synthase and CO2 metabolism (p=8.4e-08, sig=True)
**Functional:** Enriched for ATP synthase (8/8 genes, p=8.4e-08) and CO2 metabolism (7/9 genes, p=1.8e-05).
**Temporal pattern:** Genes peak at 5.5 h (early morning) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Page 3] Cluster 16 is significantly enriched for ATP synthase (8/8 genes, p=8.4e-08) and CO2 metabolism (7/9 genes, p=1.8e-05)

### Cluster 17 | N/A | high
**Name:** Prochlorococcus cluster 17 (173 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** Cluster 17 is aperiodic expressed genes, no temporal pattern or functional enrichment.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 17 is aperiodic expressed genes

### Cluster 18 | N/A | high
**Name:** Prochlorococcus cluster 18 (180 genes)
**Enrichment:** Menaquinone and ubiquinone (p=0.00045, sig=True)
**Functional:** Enriched for menaquinone and ubiquinone genes (6/9 genes, p=0.00045).
**Temporal pattern:** N/A
**Assessment notes:** Cluster 18 is mostly non-expressed genes but shows enrichment for menaquinone and ubiquinone.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 18 is non-expressed genes but enriched for menaquinone and ubiquinone genes (6/9 genes, p=0.00045)

### Cluster 2 | peaks mid-morning | medium
**Name:** Prochlorococcus cluster 2 (52 genes)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes peak at 9.6 h (mid-morning) with 24h periodicity.
**Assessment notes:** No specific functional description or enrichment reported for cluster 2.
**Sources:** Table 1
**Quotes:**
- [Page 3] The size of the clusters ranged from 22 to 138 (average 88) genes and peak transcription levels of the clusters were spread fairly evenly over the photocycle.

### Cluster 3 | peaks mid-day | high
**Name:** Prochlorococcus cluster 3 (23 genes)
**Enrichment:** Cytochrome b6/f and Photosystem II (p=0.0059, sig=True)
**Functional:** Enriched for cytochrome b6/f (3/7 genes, p=0.0059) and photosystem II (8/22 genes, p=3.7e-08).
**Temporal pattern:** Genes peak at 12.5 h (mid-day) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 3 is enriched for cytochrome b6/f (3/7 genes, p=0.0059) and photosystem II (8/22 genes, p=3.7e-08)

### Cluster 4 | peaks afternoon | medium
**Name:** Prochlorococcus cluster 4 (62 genes)
**Enrichment:** Cytochrome b6/f (p=0.073, sig=False)
**Functional:** Marginal enrichment for cytochrome b6/f (3/7 genes, p=0.073).
**Temporal pattern:** Genes peak at 15.8 h (afternoon) with 24h periodicity.
**Assessment notes:** Enrichment is marginal and not statistically significant.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 4 shows marginal enrichment for cytochrome b6/f (3/7 genes, p=0.073)

### Cluster 5 | peaks evening | high
**Name:** Prochlorococcus cluster 5 (120 genes)
**Enrichment:** Respiratory terminal oxidases and protein degradation (p=0.019, sig=True)
**Functional:** Enriched for respiratory terminal oxidases (3/3 genes, p=0.019) and degradation of proteins, peptides, and glycopeptides (5/15 genes, p=0.071).
**Temporal pattern:** Genes peak at 17.5 h (evening) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 5 is enriched for respiratory terminal oxidases (3/3 genes, p=0.019) and degradation of proteins, peptides, and glycopeptides (5/15 genes, p=0.071)

### Cluster 6 | peaks evening | high
**Name:** Prochlorococcus cluster 6 (137 genes)
**Enrichment:** Purine ribonucleotide biosynthesis (p=0.0049, sig=True)
**Functional:** Enriched for purine ribonucleotide biosynthesis (7/18 genes, p=0.0049).
**Temporal pattern:** Genes peak at 18.6 h (evening) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 6 is enriched for purine ribonucleotide biosynthesis (7/18 genes, p=0.0049)

### Cluster 7 | peaks evening | medium
**Name:** Prochlorococcus cluster 7 (121 genes)
**Enrichment:** Nitrogen metabolism (p=0.087, sig=False)
**Functional:** Marginal enrichment for nitrogen metabolism (4/8 genes, p=0.087).
**Temporal pattern:** Genes peak at 20.1 h (evening) with 24h periodicity.
**Assessment notes:** Enrichment is marginal and not statistically significant.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 7 shows marginal enrichment for nitrogen metabolism (4/8 genes, p=0.087)

### Cluster 8 | peaks evening | high
**Name:** Prochlorococcus cluster 8 (89 genes)
**Enrichment:** Chaperones (p=0.00023, sig=True)
**Functional:** Enriched for chaperones (7/14 genes, p=0.00023).
**Temporal pattern:** Genes peak at 21.0 h (evening) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 8 is enriched for chaperones (7/14 genes, p=0.00023)

### Cluster 9 | peaks late evening | high
**Name:** Prochlorococcus cluster 9 (99 genes)
**Enrichment:** RNA synthesis, modification, and DNA transcription (p=0.035, sig=True)
**Functional:** Enriched for RNA synthesis, modification, and DNA transcription (6/23 genes, p=0.035).
**Temporal pattern:** Genes peak at 22.4 h (late evening) with 24h periodicity.
**Sources:** Table 1
**Quotes:**
- [Table 1] Cluster 9 is enriched for RNA synthesis, modification, and DNA transcription (6/23 genes, p=0.035)

## Warnings

- [Biller 2018 / natl2a_periodicity / cluster not_periodic] filler phrase 'likely' in functional_description
- [Biller 2018 / natl2a_darkness_survival / cluster darkness_coculture+unique_coculture] near-identical functional_description as cluster darkness_axenic+darkness_coculture
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 1] locus tag in functional_description: PMM0958
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 8] near-identical functional_description as cluster 4
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 4] filler phrase 'possibly' in functional_description
- [Bernstein 2017 / bp1_oxygen_clusters / cluster H] locus tag in functional_description: tll1454
