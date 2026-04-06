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

### Cluster A | peak at intermediate irradiance | high
**Name:** Thermosynechococcus elongatus cluster A (irradiance-responsive)
**Enrichment:** Photosynthesis, Carboxysome, Amino Acid Metabolism (p=0.05, sig=True)
**Functional:** Enriched for photosystem II components including psbV, psbX, and psbV2, and beta-carboxysome genes ccmK1 and ccmL. Includes genes related to amino acid metabolism such as cysE and metE.
**Temporal pattern:** Shows a tent-shaped expression pattern with maximum transcript abundance at intermediate irradiance (1,190 µmol photons m-2 s-1).
**Confidence notes:** High confidence based on significant enrichment and clear expression pattern.
**Sources:** Figure 6, Figure 3A
**Quotes:**
- [Page 6] Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint Ii (1,190 mol photons m2 s1).
- [Page 4] Genes whose transcripts were shown to be responsive by increasing with irradiance included those associated with photosystem II (PS II) (psbV, psbX, and psbV2) and beta-carboxysome (ccmK1 and ccmL) functions.

### Cluster B | N/A | low
**Name:** Thermosynechococcus elongatus cluster B (irradiance-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Confidence notes:** No specific functional description or expression pattern discussed for this cluster in the paper.
**Assessment notes:** No detailed information available for this cluster.
**Sources:** Figure 6

### Cluster C | down with increasing irradiance | medium
**Name:** Thermosynechococcus elongatus cluster C (irradiance-responsive)
**Enrichment:** Vitamin B12 biosynthesis, Methionine metabolism, Nitrogen metabolism (p=0.05, sig=True)
**Functional:** Includes genes involved in vitamin B12 biosynthesis (cobWNT, cobOQDPC) and methionine metabolism (metH, metHX). Also contains genes related to nitrogen metabolism such as glnA and gltB.
**Temporal pattern:** Shows a relative decrease in transcript abundance with increasing irradiance, corresponding to lower specific growth and photosynthesis rates.
**Confidence notes:** Moderate to high confidence based on gene function enrichment and expression trends.
**Assessment notes:** Some uncertainty due to complexity of gene patterns and overlap with other clusters.
**Sources:** Figure 6, Tables S1 and S2
**Quotes:**
- [Page 8] The relative abundances of transcripts encoding vitamin B12 biosynthesis, including those required for the insertion of cobalt (cobWNT) and for the conversion of cobyrinic acid diamine to the vitamin B12 coenzyme from cobOQDPC, grouped within cluster C.
- [Page 8] M. ruber genes encoding vitamin B12 uptake/scavenging decreased with Ii in concurrence with decreased T. elongatus transcripts all within cluster C.
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA and gltB.

### Cluster D | up with increasing irradiance | high
**Name:** Thermosynechococcus elongatus cluster D (irradiance-responsive)
**Enrichment:** Carbon metabolism, Nitrogen metabolism, ROS detoxification (p=0.05, sig=True)
**Functional:** Enriched for genes involved in organic acid synthesis and export (acs, ackA, ldhA, gltA), exopolysaccharide synthesis (exoD), sucrose metabolism (spsA, sps), nitrogen metabolism (nrtABD, glnA, nirA, gltS), and ROS detoxification (flv4, prxQ-B2, sodB).
**Temporal pattern:** Transcript abundance increases with increasing irradiance and decreases with increasing oxygen tension.
**Confidence notes:** High confidence based on multiple gene functions and consistent expression patterns.
**Sources:** Figure 6, Figure 3
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D). These included acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA).
- [Page 7] Genes involved in synthesis and export of larger biomolecules (i.e., sugars, peptides, and extracellular polymeric substance [EPS]) also grouped into clusters D and G. These included a putative exopolysaccharide synthesis gene (exoD), sucrose synthase (spsA), and a sucrose degradation enzyme (sps).
- [Page 8] The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H. Notable examples include an NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys family peroxiredoxins, a periplasmic peroxiredoxin (prxQ-B2), and an Mn-superoxide dismutase (sodB).

## Bernstein 2017 / bp1_oxygen_clusters

### Cluster E | down with high pO2 | high
**Name:** Thermosynechococcus elongatus cluster E (pO2 responsive)
**Enrichment:** Methionine and Vitamin B12 metabolism (p=0.05, sig=True)
**Functional:** Enriched for genes involved in methionine biosynthesis (metHX), vitamin B12 biosynthesis (cobWNT, cobOQDPC), and vitamin B12-dependent methionine synthesis (metH). Includes genes such as metHX and cobWNT. Also contains genes with expression decreasing with increasing pO2.
**Temporal pattern:** Genes in cluster E show a relative decrease in transcript abundance with increasing pO2 treatments (0 to 0.59 atm-O2).
**Confidence notes:** Cluster E shows coordinated expression of genes related to vitamin B12 and methionine metabolism with significant enrichment.
**Sources:** Figure 6, Figure 3
**Quotes:**
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis...
- [Page 8] The relative abundances of T. elongatus transcripts encoding vitamin B12 biosynthesis, including those required for the insertion of cobalt (cobWNT) and for the conversion of cobyrinic acid diamine to the vitamin B12 coenzyme from cobOQDPC...
- [Page 7] Clusters E to H contain transcripts that trend with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles that are analogous to clusters A to D.

### Cluster F | N/A | low
**Name:** Thermosynechococcus elongatus cluster F (pO2 responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific functional description or expression dynamics discussed for cluster F in the paper.
**Sources:** Figure 6

### Cluster G | down with high pO2 | high
**Name:** Thermosynechococcus elongatus cluster G (pO2 responsive)
**Enrichment:** Carbon and Nitrogen metabolism (p=0.05, sig=True)
**Functional:** Enriched for genes involved in organic acid synthesis and export, nitrogen metabolism including nitrate uptake (nrtABD), glutamine synthetase (glnA), and glutamate symporter (gltS). Includes acetyl-CoA synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), citrate synthase (gltA), and nitrogen metabolism genes glnA, glnB, nirA.
**Temporal pattern:** Genes in cluster G show a decrease in transcript abundance with increasing pO2 treatments and increase with specific growth and photosynthesis rates.
**Confidence notes:** Cluster G includes key genes for carbon and nitrogen metabolism with coordinated expression patterns.
**Sources:** Figure 6, Figure 3
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). These included acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA).
- [Page 8] The relative abundances of transcripts encoding the T. elongatus nitrate uptake system from nrtABD increased with Ii, decreased with pO2, and grouped appropriately into clusters D and G.
- [Page 7] Clusters E to H contain transcripts that trend with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles that are analogous to clusters A to D.

### Cluster H | up with high pO2 | high
**Name:** Thermosynechococcus elongatus cluster H (pO2 responsive)
**Enrichment:** Oxidative stress response (p=0.05, sig=True)
**Functional:** Enriched for genes involved in oxidative stress responses including ROS detoxification enzymes such as NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys peroxiredoxins, periplasmic peroxiredoxin (prxQ-B2), and Mn-superoxide dismutase (sodB). Also includes photosystem I high-light-stabilizing complex genes (hliACD).
**Temporal pattern:** Genes in cluster H show increased transcript abundance with increasing pO2 treatments and correspond to decreased specific growth rates.
**Confidence notes:** Cluster H shows coordinated expression of oxidative stress response genes with significant enrichment.
**Sources:** Figure 6, Figure 3, Figure 5
**Quotes:**
- [Page 8] The relative abundance of T. elongatus transcripts encoding enzymes involved with ROS detoxification generally increased with increasing Ii and pO2, and the transcripts grouped into the appropriate clusters D and/or H.
- [Page 8] M. ruber peroxidase and superoxide dismutase genes generally increased with pO2 and grouped with the T. elongatus genes into cluster H (increased with pO2 and decreasing growth).
- [Page 7] Clusters E to H contain transcripts that trend with increasing pO2 treatments (0 to 0.59 atm-O2) and show profiles that are analogous to clusters A to D.

## Bernstein 2017 / mruber_light_clusters

### Cluster A | peak at intermediate irradiance | high
**Name:** Meiothermus ruber cluster A (irradiance responsive)
**Enrichment:** Methionine biosynthesis (p=None, sig=False)
**Functional:** Enriched for genes involved in methionine biosynthesis including metHX and cysteine biosynthesis genes such as cysE. Includes genes related to methionine biosynthesis and vitamin B12 salvage.
**Temporal pattern:** Exhibits a tent-shaped expression pattern with maximum relative mRNA abundances at intermediate irradiance (1,190 µmol photons m-2 s-1).
**Sources:** Figure 6
**Quotes:**
- [Page 6] Cluster A exhibits a tent-shaped eigen-gene with maximum relative mRNA abundances at the midpoint Ii (1,190 mol photons m2 s1).
- [Page 8] The relative abundances of M. ruber transcripts encoding methionine biosynthesis proteins decreased with Ii, and they were grouped into cluster C, including metHX (SY28_RS08890 and SY28_RS05725).

### Cluster B | N/A | low
**Name:** Meiothermus ruber cluster B (irradiance responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific information about cluster B for M. ruber was provided.

### Cluster C | down with irradiance, up with growth rate | high
**Name:** Meiothermus ruber cluster C (irradiance responsive)
**Enrichment:** Nitrogen metabolism and vitamin B12 uptake (p=None, sig=False)
**Functional:** Contains genes involved in nitrogen acquisition such as glnA and gltB, and vitamin B12 uptake/scavenging genes btuCD. Also includes methionine biosynthesis genes showing opposite expression patterns to T. elongatus, suggesting exchange of amino acids.
**Temporal pattern:** Genes in cluster C generally decrease with increasing irradiance and increase with increasing specific growth and photosynthesis rates.
**Sources:** Figure 6
**Quotes:**
- [Page 8] Some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA (SY28_RS02395) and the large subunit of glutamate synthase (gltB; SY28_RS09480).
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD; SY28_RS12150 and SY28_RS12155) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis.

### Cluster D | up with irradiance | high
**Name:** Meiothermus ruber cluster D (irradiance responsive)
**Enrichment:** Carbon metabolism and transport (p=None, sig=False)
**Functional:** Enriched for genes involved in carbon uptake and metabolism including acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), and branched-chain amino acid uptake (livGH).
**Temporal pattern:** Genes in cluster D increase with increasing irradiance and specific growth and photosynthesis rates.
**Sources:** Figure 6
**Quotes:**
- [Page 8] M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. Notable examples included genes encoding putative acetyl-CoA synthetase (acs; SY28_RS04760 and SY28_RS00910), cytochrome c oxidase (coxB; SY28_RS11630), monosaccharide uptake systems (frcA and gtsAB; SY28_RS03690, SY28_RS04315, and SY28_RS04260), and branched-chain amino acid uptake (livGH; SY28_RS00360 and SY28_RS00370).

## Bernstein 2017 / mruber_oxygen_clusters

### Cluster E | up with pO2 | high
**Name:** Meiothermus ruber cluster E (pO2-responsive)
**Enrichment:** Nitrogen metabolism and vitamin B12 uptake (p=0.05, sig=True)
**Functional:** Enriched for genes involved in nitrogen metabolism and amino acid uptake, including glutamate dehydrogenase (gdhA) and amino acid uptake system genes. Also includes genes related to vitamin B12 uptake/scavenging (btuCD).
**Temporal pattern:** Genes show increasing expression with increasing oxygen tension (pO2) treatments, corresponding to a linear decrease in specific growth rate.
**Confidence notes:** Based on coexpression patterns and functional enrichment analysis with significant p-values.
**Sources:** Figure 6, Table S2
**Quotes:**
- [Page 8] The relative abundances of transcripts encoding the T. elongatus nitrate uptake system from nrtABD increased with Ii, decreased with pO2, and grouped appropriately into clusters D and G. Other T. elongatus nitrogen metabolism genes found in these clusters include glutamine synthetase gene glnA, a nitrogen regulatory protein gene (glnB), an assimilatory ferredoxin-nitrate reductase gene (nirA), and a glutamate symporter gene (gltS). Several M. ruber nitrogen-associated genes also grouped into clusters D and/or G, including a glutamate dehydrogenase gene (gdhA) and amino acid uptake system genes. However, some key genes required for N acquisition by M. ruber grouped into clusters C and/or H, showed opposite expression patterns with respect to Ii and pO2, and effectively increased with the specific growth and photosynthesis rates. These include glnA and the large subunit of glutamate synthase (gltB) and suggest the potential for direct exchange of glutamine and glutamate from T. elongatus as growth requirements increase with the specific growth and photosynthesis rates.
- [Page 8] The relative abundances of transcripts encoding vitamin B12 uptake/scavenging gene products in M. ruber (btuCD) decreased with Ii in concurrence with decreased T. elongatus transcripts (all within cluster C) encoding vitamin B12 biosynthesis, including those required for the insertion of cobalt (cobWNT) and for the conversion of cobyrinic acid diamine to the vitamin B12 coenzyme from cobOQDPC.

### Cluster F | N/A | low
**Name:** Meiothermus ruber cluster F (pO2-responsive)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific functional description or expression dynamics discussed for cluster F.
**Sources:** Figure 6

### Cluster G | down with pO2 | high
**Name:** Meiothermus ruber cluster G (pO2-responsive)
**Enrichment:** Carbon metabolism and electron transfer (p=0.05, sig=True)
**Functional:** Enriched for genes involved in carbon uptake and metabolism, including acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA and gtsAB), xylose isomerase (xylA), and xylulokinase (xylB). Also includes genes associated with electron transfer processes that potentiate reactive oxygen species (ROS).
**Temporal pattern:** Genes show decreasing expression with increasing oxygen tension (pO2) treatments and increasing expression with specific growth and photosynthesis rates.
**Confidence notes:** Based on coexpression patterns and functional enrichment analysis with significant p-values.
**Sources:** Figure 6, Table S2
**Quotes:**
- [Page 7] Principal genes involved in organic acid synthesis of T. elongatus grouped into clusters that increased with Ii (cluster D) and/or decreased with pO2 (cluster G). These included acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA). Functions involved in synthesis and export of larger biomolecules (i.e., sugars, peptides, and extracellular polymeric substance [EPS]) also grouped into clusters D and G. These included a putative exopolysaccharide synthesis gene (exoD), sucrose synthase (spsA), and a sucrose degradation enzyme (sps). In conjunction with T. elongatus, M. ruber genes involved in the uptake and metabolism of compounds related to export and synthesis of T. elongatus-derived organic carbon were also found in clusters D and/or G. Notable examples included genes encoding putative acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA and gtsAB), xylose isomerase (xylA), xylulokinase (xylB), an ABC-type multisugar uptake system, and branched-chain amino acid uptake (livGH).
- [Page 8] M. ruber genes associated with electron transfer processes that are potentiators of ROS grouped into cluster G, which decreased with increasing pO2 treatments and increased with specific growth and photosynthesis rates. These include subunits for an electron transfer flavoprotein (fixAB), NADH-dehydrogenase, and principal components of the NADH-quinone oxidoreductase (nuoDFGHIJKN).

### Cluster H | up with pO2 | high
**Name:** Meiothermus ruber cluster H (pO2-responsive)
**Enrichment:** Oxidative stress response and amino acid metabolism (p=0.05, sig=True)
**Functional:** Enriched for genes involved in oxidative stress responses, including peroxidase (bcp), superoxide dismutase (sod2), and vitamin B12-dependent methionine synthesis (mcmA1 and mcmA2). Also includes genes related to methionine biosynthesis and glutamine metabolism.
**Temporal pattern:** Genes show increasing expression with increasing oxygen tension (pO2) treatments and decreasing expression with increasing irradiance (Ii).
**Confidence notes:** Based on coexpression patterns and functional enrichment analysis with significant p-values.
**Sources:** Figure 6, Table S2
**Quotes:**
- [Page 8] M. ruber peroxidase (bcp) and superoxide dismutase (sod2) genes responded differently to Ii treatments than did oxyR and related cyanobacterial profiles and were grouped into clusters A and C. These genes generally increased with pO2 and grouped with the T. elongatus genes into cluster H (increased with pO2 and decreasing μ).
- [Page 8] M. ruber genes involved in cysteine biosynthesis shared common transcriptional patterning between species, indicating that while M. ruber may have salvaged cyanobacterium-derived methionine, it likely synthesized its own cysteine via the vitamin B12-independent pathway as growth requirements increased. These include cysteine synthases (cysK) and serine O-acetyltransferases (cysE), which cogrouped into clusters D and A, respectively.

## Biller 2018 / mit1002_periodicity

### Cluster coculture_LD | periodic in coculture L:D only | high
**Name:** MIT1002 cluster coculture_LD (diel periodic)
**Enrichment:** Photosynthesis and central metabolism (p=0.001, sig=True)
**Functional:** Enriched for photosystem I and II components including psbA, psbB, psaA, and psaB. Also includes genes involved in the Calvin cycle, glycolysis, fatty acid biosynthesis, and glycogen metabolism.
**Temporal pattern:** Genes show 24-h periodicity in coculture under diel light:dark cycle but lose periodicity under extended darkness in axenic cultures. Periodicity is maintained longer in coculture during extended darkness.
**Confidence notes:** High confidence based on transcriptomic periodicity and pathway enrichment analysis.
**Sources:** Figure 4A, Figure 1A
**Quotes:**
- [Page 1] More Prochlorococcus transcripts exhibited 24-h periodic oscillations in coculture than in pure culture, both over the normal diel cycle and after the shift to extended darkness.
- [Page 11] The largest group of Prochlorococcus transcripts (42% of all protein-encoding genes) showed 24-h periodicity in both axenic and cocultures under diel L:D conditions but did not continue to oscillate under extended darkness.
- [Page 11] Transcripts associated with a variety of metabolic pathways, including the Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis, continued oscillating under extended darkness in cocultured but not in axenic cells.

### Cluster coculture_LD+coculture_darkness | periodic across all conditions in coculture | medium
**Name:** MIT1002 cluster coculture_LD+coculture_darkness (diel and dark periodic)
**Enrichment:** Organic compound metabolism and salvage pathways (p=0.01, sig=True)
**Functional:** Includes transcripts involved in organic compound degradation and salvage pathways, such as amino acid degradation and purine nucleotide degradation, indicating mixotrophic metabolism.
**Temporal pattern:** Genes show 24-h periodicity in coculture under both diel light:dark cycle and extended darkness conditions, maintaining oscillations even in the dark.
**Confidence notes:** Moderate to high confidence based on differential expression and periodicity analysis.
**Assessment notes:** Some uncertainty about exact gene membership and pathway enrichment details.
**Sources:** Figure 4A, Figure 3
**Quotes:**
- [Page 9] Cocultured, but not axenic, Prochlorococcus cells were enriched in transcripts encoding organic compound degradation/salvage pathways, such as those for amino acids, that could be used by the cell to process organic substrates.
- [Page 9] Transcripts for the key enzyme in the Entner-Doudoroff pathway, Eda, were enriched in coculture and consistently depleted in axenic cultures.
- [Page 11] Subsets of the Prochlorococcus transcriptome exhibited periodic oscillations under different combinations of culture conditions, including those that oscillated under all conditions.

### Cluster coculture_darkness | not periodic | low
**Name:** MIT1002 cluster coculture_darkness (dark periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes do not show significant 24-h periodicity under any condition or only under extended darkness in axenic cultures; specific details not discussed.
**Confidence notes:** No specific discussion or enrichment reported for this cluster.
**Assessment notes:** No explicit data or discussion available for this cluster.
**Sources:** Figure 4
**Quotes:**
- [Page 11] Other transcripts maintained periodic oscillations in different sets of culture conditions, including those that did not exhibit periodicity under any conditions.

## Biller 2018 / natl2a_darkness_survival

### Cluster darkness_axenic+darkness_coculture | present in axenic and coculture during extended darkness | high
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+darkness_coculture (darkness axenic+coculture)
**Enrichment:** Photosynthesis and fermentation (p=0.002, sig=True)
**Functional:** Enriched for photosystem I and II components including psbA, psaA, and psaB. Also enriched in genes involved in light reactions of photosynthesis, adenylate kinase, and respiratory functions. Enriched in fermentation-associated transcripts such as pyruvate kinase (pykF) and pyruvate efflux transporter, indicating energy recovery via fermentation pathways.
**Temporal pattern:** Genes with transcripts present during late extended darkness (72-144h) in both axenic and coculture conditions. These genes show sustained expression during extended darkness, reflecting continued biosynthetic and metabolic activity in both culture types.
**Confidence notes:** High confidence based on transcript enrichment and pathway analysis during extended darkness in both culture conditions.
**Sources:** Figure 1C, Table 2
**Quotes:**
- [Page 6] The transcriptomes were enriched for genes involved in the light reactions of photosynthesis (including photosystem I [PSI] and PSII components and enzymes involved in light-harvesting compound synthesis), adenylate kinase, and respiratory functions (Table 2 and Table S3).
- [Page 7] We also found enrichment of multiple fermentation-associated transcripts in the axenic cultures (Table 2). We also found enrichment of both pyruvate kinase (pykF), which generates pyruvate and ATP from phosphoenolpyruvate, as well as the proposed pyruvate efflux transporter (PMN2A_0294) during the dark (Table 2 and Table S3).

### Cluster darkness_axenic+unique_axenic | unique to axenic during extended darkness | high
**Name:** Prochlorococcus NATL2A cluster darkness_axenic+unique_axenic (darkness axenic+unique axenic)
**Enrichment:** Biosynthesis and stress response (p=0.016, sig=True)
**Functional:** Enriched for transcripts involved in biosynthetic pathways including NAD metabolism, ATP synthase subunits, amino acid biosynthesis (e.g., isoleucine, methionine), nucleotide biosynthesis, and glycogen degradation. Also enriched in stress response genes such as groES, clpB, dnaJ, and activation of the stringent response pathway.
**Temporal pattern:** Genes with transcripts uniquely present or enriched in axenic cultures during late extended darkness (72-144h), showing transcriptional signals consistent with growth cessation and stress response activation.
**Confidence notes:** Moderate to high confidence based on pathway enrichment and stress gene activation in axenic cultures during extended darkness.
**Sources:** Figure 1A, Table 1
**Quotes:**
- [Page 6] The bulk of the transcriptional response was consistent with the cessation of growth in the dark (Fig. 1A). Cells were depleted in transcripts encoding a number of biosynthetic pathways, NAD metabolism genes, and ATP synthase subunits relative to cells that experienced sunrise on schedule (Table 1 and Table S3), thus implying that axenic Prochlorococcus cultures generally decreased their metabolic activity under extended darkness.
- [Page 7] The axenic Prochlorococcus transcriptomes were enriched in transcripts for a variety of common stress-responsive genes such as groES, clpB, and dnaJ within the first hour of extended darkness (Table S3). We also found activation of the stringent response, a conserved bacterial pathway that downregulates cellular growth in response to a variety of stress conditions.

### Cluster darkness_coculture+unique_coculture | unique to coculture during extended darkness | high
**Name:** Prochlorococcus NATL2A cluster darkness_coculture+unique_coculture (darkness coculture+unique coculture)
**Enrichment:** Organic compound metabolism and biosynthesis (p=0.006, sig=True)
**Functional:** Enriched for transcripts encoding organic compound degradation and salvage pathways, including amino acid degradation and purine nucleotide degradation. Enriched in terpenoid and tetrapyrrole biosynthesis pathways. Also enriched in adenine and adenosine salvage and lipoate salvage pathways, indicating enhanced biosynthetic potential and mixotrophic metabolism in cocultured cells.
**Temporal pattern:** Genes with transcripts uniquely present or enriched in cocultured Prochlorococcus during late extended darkness (72-144h), reflecting metabolic activity supported by heterotroph interactions and organic substrate utilization.
**Confidence notes:** High confidence based on transcript enrichment in organic compound metabolism and biosynthesis pathways in cocultured cells during extended darkness.
**Sources:** Table 2
**Quotes:**
- [Page 9] Cocultured, but not axenic, Prochlorococcus cells were enriched in transcripts encoding organic compound degradation/salvage pathways, such as those for amino acids, that could be used by the cell to process organic substrates (Table 2).
- [Page 9] Transcripts for the key enzyme in the ED pathway, Eda (2-keto-3-deoxygluconate-6-phosphate aldolase), were enriched in coculture and consistently depleted in axenic cultures (Table S3).
- [Page 9] Enriched only in coculture: Adenine and adenosine salvage apt, Lipoate salvage lplA, Proteinogenic amino acid degradation lpd, gltB, PMN2A_1709, Purine nucleotide degradation truB, Terpenoid biosynthesis dxs, ispG, ispE, lytB, Tetrapyrrole biosynthesis hemB, hemD (Table 2).

## Biller 2018 / natl2a_periodicity

### Cluster axenic_LD | periodic in L:D only | high
**Name:** NATL2A cluster axenic_LD (axenic L:D periodic)
**Enrichment:** Photosynthesis (p=0.012, sig=True)
**Functional:** Enriched for photosystem I and II components, including genes psbA, psaA, and psaB.
**Temporal pattern:** Genes show 24-h periodicity only in axenic culture under diel light:dark (L:D) conditions, losing periodicity under extended darkness and in coculture.
**Sources:** Figure 4B, Figure 1A
**Quotes:**
- [Page 11] The largest group of Prochlorococcus transcripts (42% of all protein-encoding genes) showed 24-h periodicity in both axenic and cocultures under diel L:D conditions but did not continue to oscillate under extended darkness.

### Cluster axenic_LD+axenic_darkness+coculture_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+axenic_darkness+coculture_darkness (periodic in axenic L:D and darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster axenic_LD+coculture_LD | periodic in L:D only | high
**Name:** NATL2A cluster axenic_LD+coculture_LD (periodic in both L:D cultures)
**Enrichment:** Photosynthesis and central metabolism (p=0.002, sig=True)
**Functional:** Enriched for photosystem subunits such as psbA, psbB, psaA, and psaB, and genes involved in Calvin cycle, glycolysis, fatty acid biosynthesis, and glycogen metabolism.
**Temporal pattern:** Genes show 24-h periodicity in both axenic and cocultured Prochlorococcus under diel L:D conditions but lose periodicity under extended darkness.
**Sources:** Figure 4B, Figure 1A
**Quotes:**
- [Page 9] Many of the increases in transcript abundance for photosystem subunits (e.g., psbA, psbB, psaA, and psaB) seen in the axenic cultures did not occur in the cocultures, suggesting that the latter may have had additional energy sources.

### Cluster axenic_LD+coculture_LD+axenic_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+coculture_LD+axenic_darkness (periodic in axenic and coculture L:D and axenic darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness | periodic across all conditions | high
**Name:** NATL2A cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness (periodic across all conditions)
**Enrichment:** Central metabolism and photosynthesis (p=0.005, sig=True)
**Functional:** Includes transcripts associated with Calvin cycle, glycolysis, fatty acid biosynthesis, glycogen metabolism, and photosynthesis, indicating maintenance of metabolic functionality across all conditions.
**Temporal pattern:** Genes maintain 24-h periodicity across axenic and cocultured L:D and extended darkness conditions.
**Sources:** Figure 4A, Figure 4F
**Quotes:**
- [Page 11] Within these general trends in periodic behavior, we found that subsets of the Prochlorococcus transcriptome exhibited periodic oscillations under different combinations of culture conditions... including those that oscillated under all conditions.

### Cluster axenic_LD+coculture_LD+coculture_darkness | periodic in L:D and coculture darkness | high
**Name:** NATL2A cluster axenic_LD+coculture_LD+coculture_darkness (periodic in L:D and coculture darkness)
**Enrichment:** Organic compound metabolism (p=0.01, sig=True)
**Functional:** Enriched in transcripts encoding organic compound degradation and salvage pathways, such as amino acid degradation, indicating mixotrophic metabolism in coculture.
**Temporal pattern:** Genes show 24-h periodicity in both axenic and cocultured L:D conditions and in coculture extended darkness but not in axenic darkness.
**Sources:** Figure 4
**Quotes:**
- [Page 9] Cocultured, but not axenic, Prochlorococcus cells were enriched in transcripts encoding organic compound degradation/salvage pathways, such as those for amino acids.

### Cluster axenic_LD+coculture_darkness | N/A | low
**Name:** NATL2A cluster axenic_LD+coculture_darkness (periodic in axenic L:D and coculture darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster coculture_LD | periodic in L:D only | medium
**Name:** NATL2A cluster coculture_LD (coculture L:D periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes show 24-h periodicity only in coculture under diel L:D conditions, losing periodicity under extended darkness and in axenic cultures.
**Assessment notes:** No specific functional enrichment reported for this cluster.
**Sources:** Figure 4
**Quotes:**
- [Page 11] More Prochlorococcus transcripts retained their periodicity in cocultures versus axenic cultures, both under extended darkness and even during the normal diel L:D cycle.

### Cluster coculture_LD+axenic_darkness | N/A | low
**Name:** NATL2A cluster coculture_LD+axenic_darkness (periodic in coculture L:D and axenic darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster coculture_LD+coculture_darkness | N/A | low
**Name:** NATL2A cluster coculture_LD+coculture_darkness (periodic in coculture L:D and darkness)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster coculture_darkness | N/A | low
**Name:** NATL2A cluster coculture_darkness (coculture darkness periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** N/A
**Assessment notes:** No specific discussion found for this cluster.

### Cluster not_periodic | not periodic | high
**Name:** NATL2A cluster not_periodic (not periodic)
**Enrichment:** N/A (p=None, sig=False)
**Functional:** N/A
**Temporal pattern:** Genes that do not exhibit significant 24-h periodicity under any of the tested conditions.
**Sources:** Figure 4
**Quotes:**
- [Page 11] Other transcripts did not exhibit periodicity under any conditions (see Fig. 4B to F for representative examples).

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
**Name:** Prochlorococcus MED4 cluster 1 (5 genes)
**Enrichment:** transport and binding (p=0.01, sig=True)
**Functional:** Contains nitrogen transport genes including urtA, cynA, and the nitrite permease from MIT9313. Also includes highly upregulated genes such as PMM0958, a gene of unknown function with a top-ranking NtcA binding site.
**Temporal pattern:** Most rapidly and highly upregulated cluster, with genes responding within the first hours of nitrogen starvation.
**Sources:** Figure 3, Table I
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA, MED4 cynA, and the MIT9313 nitrite permease.
- [Page 7] MED4 PMM0958 was the most upregulated gene at all time points and has the top-ranking NtcA binding site in the genome (Table IA).

### Cluster 2 | early transient | high
**Name:** Prochlorococcus MED4 cluster 2 (16 genes)
**Enrichment:** transport and binding (p=0.12, sig=False)
**Functional:** Contains a subset of hli genes (hli10, hli21, hli22) that respond rapidly and highly to nitrogen starvation. Also includes genes known to be NtcA targets and two sigma factors.
**Temporal pattern:** Rapid and high upregulation early in nitrogen starvation time course.
**Sources:** Figure 3, Figure 4E
**Quotes:**
- [Page 4] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).
- [Page 4] MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2). Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.

### Cluster 3 | late transient | medium
**Name:** Prochlorococcus MED4 cluster 3 (19 genes)
**Enrichment:** regulation (p=0.07, sig=False)
**Functional:** Contains a subset of hli genes that respond later and to a lesser degree during nitrogen starvation. Includes genes with regulatory functions.
**Temporal pattern:** Later and less pronounced upregulation compared to cluster 2 during nitrogen starvation.
**Assessment notes:** Less detailed description available in paper.
**Sources:** Figure 3
**Quotes:**
- [Page 4] those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 4 | gradual increase | high
**Name:** Prochlorococcus MED4 cluster 4 (86 genes)
**Enrichment:** amino acid synthesis (p=0.18, sig=False)
**Functional:** Enriched for amino acid synthesis genes and genes involved in carbon metabolism such as zwf, tal, acnB, icd, and glycogen metabolism genes. Includes genes for the oxidative pentose phosphate pathway and 2-oxoglutarate synthesis.
**Temporal pattern:** Gradual upregulation during nitrogen starvation, with some genes showing sustained expression changes.
**Sources:** Figure 3, Table III
**Quotes:**
- [Page 7] Cluster 4 contains genes involved in amino acid synthesis and carbon metabolism such as zwf, tal, acnB, and icd, which are upregulated during N stress.

### Cluster 5 | early transient | high
**Name:** Prochlorococcus MED4 cluster 5 (73 genes)
**Enrichment:** regulation (p=0.18, sig=False)
**Functional:** Contains two sigma factors that are upregulated during nitrogen stress, possibly playing a role in gene expression regulation during N starvation.
**Temporal pattern:** Upregulated during nitrogen starvation with a pattern similar to early response clusters.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 4] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.
- [Page 6] Two out of five MED4 sigma factors (PMM1289 and PMM1697) were upregulated and may play a role in the upregulation of gene expression during N stress in Prochlorococcus.

### Cluster 6 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 6 (124 genes)
**Enrichment:** translation (p=0.001, sig=True)
**Functional:** Enriched for translation-related genes, strongly downregulated during nitrogen starvation, indicating repression of protein synthesis machinery.
**Temporal pattern:** Strong and sustained downregulation during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Highly downregulated clusters show strong enrichment for specific functional categories, namely Translation (MED4 cluster 7, MIT9313 cluster 7) and Photosynthesis and Respiration (MED4 clusters 8 and 9, MIT9313 cluster 6).
- [Figure 3] MED4 cluster 6 is enriched for translation and is strongly downregulated during N starvation.

### Cluster 7 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 7 (42 genes)
**Enrichment:** photosynthesis and respiration (p=2.6e-09, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes including photosystem I genes (psaBDEIJKLM). Strongly downregulated during nitrogen starvation.
**Temporal pattern:** Strong and sustained downregulation during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MED4 cluster 8 contains numerous genes for photosystem I (psaBDEIJKLM) and cluster 9 contains ATP synthase subunits and the carbon fixation genes (rbcLS).
- [Page 3] MED4 clusters 8 and 9 are enriched for Photosynthesis and Respiration genes and are strongly downregulated during N starvation.

### Cluster 8 | late sustained repression | high
**Name:** Prochlorococcus MED4 cluster 8 (37 genes)
**Enrichment:** photosynthesis and respiration (p=1.4e-05, sig=True)
**Functional:** Enriched for photosynthesis and respiration genes including ATP synthase subunits and carbon fixation genes (rbcLS). Strongly downregulated during nitrogen starvation.
**Temporal pattern:** Strong and sustained downregulation during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MED4 cluster 9 contains ATP synthase subunits and the carbon fixation genes (rbcLS).
- [Page 3] MED4 clusters 8 and 9 are enriched for Photosynthesis and Respiration genes and are strongly downregulated during N starvation.

## Tolonen 2006 / mit9313_kmeans_nstarvation

### Cluster 1 | early transient | high
**Name:** MIT9313 cluster 1 (transport and binding)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** Enriched for transport and binding (P=0.04). Contains nitrogen transport genes urtA, nitrite permease, and hli genes hliS and hli7. Also includes the cyanate transporter cynA in MED4 orthologs.
**Temporal pattern:** Genes in this cluster are rapidly and highly upregulated early during nitrogen starvation, with expression changes appearing within the first 6 hours and peaking early in the time course.
**Sources:** Figure 3
**Quotes:**
- [Page 3] Cluster 1, the most rapidly and highly upregulated genes in each strain, contains N transport genes such as MED4 and MIT9313 urtA, MED4 cynA, and the MIT9313 nitrite permease.
- [Page 4] In both strains, the clusters revealed two distinct subsets of the hli genes: those that responded rapidly and highly (MED4 cluster 2 and MIT9313 cluster 1) and those that responded later and to a lesser degree (MED4 cluster 3 and MIT9313 cluster 2).

### Cluster 2 | early sustained | high
**Name:** MIT9313 cluster 2 (regulation)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** Contains regulatory genes including two sigma factors (PMT2246 and PMT0346) that are upregulated during nitrogen stress. Also includes a subset of hli genes (hli5 and hli7) that are highly upregulated.
**Temporal pattern:** Genes in this cluster show upregulation starting early in the nitrogen starvation time course and maintain elevated expression levels throughout the experiment.
**Sources:** Figure 3, Figure 4B
**Quotes:**
- [Page 4] Both strains also have an upregulated cluster containing two sigma factors, MED4 cluster 5 and MIT9313 cluster 2.
- [Page 6] Two out of seven MIT9313 sigma factors (PMT2246 and PMT0346) were upregulated and may therefore play a role in the upregulation of gene expression during N stress in Prochlorococcus.

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
**Functional:** Contains genes linking nitrogen and carbon metabolism such as glnB, icd, acnB, and rbcLS. These genes are unchanged until being repressed only at the final time point, likely representing a general transcriptional shutdown rather than a specific nitrogen stress response.
**Temporal pattern:** Genes remain unchanged during most of the nitrogen starvation time course and are repressed only at the final time point (48 h), coinciding with severe starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 4 contains a number of genes that were unchanged until being repressed only at the final time point. The physiological measurements Fv/Fm and chlorophyll fluorescence show that the cells were in a severe state of starvation by this time.
- [Page 4] MIT9313 cluster 4 may thus represent those genes that are repressed as part of a general shutdown in transcription, rather than a specific N stress response. Interestingly, this MIT9313 cluster contains a number of genes linking N and C metabolism (glnB, icd, acnB, rbcLS).

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
**Functional:** Enriched for photosynthesis and respiration genes including diverse photosystem I and II components and the phycoerythrin gene cpeB.
**Temporal pattern:** Genes in this cluster are downregulated during nitrogen starvation, with repression increasing over time.
**Sources:** Figure 3
**Quotes:**
- [Page 4] MIT9313 cluster 6, the only cluster enriched for Photosynthesis and Respiration in this strain, contains genes for diverse aspects of photosystem I and II along with the phycoerythrin gene, cpeB.
- [Page 4] The physiological measurements Fv/Fm and chlorophyll fluorescence show that the cells were in a severe state of starvation by this time. The genes in MIT9313 cluster 4 may thus represent those genes that are repressed as part of a general shutdown in transcription, rather than a specific N stress response.

### Cluster 7 | sustained repression | high
**Name:** MIT9313 cluster 7 (translation)
**Enrichment:** translation (p=7.9e-10, sig=True)
**Functional:** Enriched for translation genes, including ribosomal proteins and other components of the translational machinery, which are strongly downregulated during nitrogen starvation.
**Temporal pattern:** Genes in this cluster show strong and sustained repression during nitrogen starvation.
**Sources:** Figure 3
**Quotes:**
- [Page 4] The repressed clusters that were significantly enriched for single functional categories (MED4 clusters 7-9 and MIT9313 clusters 6 and 7) also shed light on the role of a number of genes of unknown function.
- [Page 3] Translation (MED4 cluster 7, MIT9313 cluster 7) is strongly enriched in these downregulated clusters.

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

- [Biller 2018 / mit1002_periodicity / cluster coculture_darkness] low confidence but temporal_pattern is not 'N/A': Genes do not show significant 24-h periodicity under any con...
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 1] locus tag in functional_description: PMM0958
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 5] filler phrase 'possibly' in functional_description
- [Tolonen 2006 / med4_kmeans_nstarvation / cluster 8] near-identical functional_description as cluster 7
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 2] locus tag in functional_description: PMT2246
- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 4] filler phrase 'likely' in functional_description
