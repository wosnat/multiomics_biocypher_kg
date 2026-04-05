# Cluster Extraction Report

## Alonso-Saez 2023 / mit9301_softclusters_thermal_acclimation

### Cluster Cluster A | up | high
**Name:** Prochlorococcus cluster A (up, core metabolism)
**Enrichment:** core metabolism and photosynthesis (p=None, sig=False)
**Functional:** Cluster A includes genes related to carbon fixation and assimilation such as RuBisCO (rbcLS), CO2 transporters (csoS2), carboxysome shell proteins (ccmK), Calvin cycle enzymes (gap2, tktA, glpX, pgk, cbbA), and glycogen synthesis (glgABC). It also contains ATP synthesis genes (atpADE) and some photosystem II components (psbA, psbC, psbD). These genes are consistently expressed during daytime across the thermal gradient, representing essential metabolic pathways maintained under all temperature conditions.
**Behavioral:** Genes in this cluster show stable expression during daytime across temperatures, indicating maintenance of core metabolic functions.
**Notes:** High confidence based on gene function and consistent expression pattern.

### Cluster Cluster B | down | high
**Name:** Prochlorococcus cluster B (down, photosynthesis and growth)
**Enrichment:** photosynthesis and growth-related metabolism (p=None, sig=False)
**Functional:** Cluster B contains genes involved in photosystem I and II components (psaABDEFKL, psbBJH, psbO), photosynthetic electron transport (petACGNM), and nutrient acquisition genes such as pstS. These genes show a decreasing expression trend from optimum to minimum temperature during daytime, paralleling the decline in growth rate, suggesting association with metabolic processes limiting growth under cold conditions.
**Behavioral:** Genes in this cluster are downregulated during daytime as temperature decreases, correlating with reduced growth rates.
**Notes:** Moderate to high confidence based on correlation with growth and gene functions.

### Cluster Cluster C | up | high
**Name:** Prochlorococcus cluster C (up, cold stress response day)
**Enrichment:** stress response and regulation (p=None, sig=False)
**Functional:** Cluster C includes genes involved in global stress response such as cellular chaperones (groES, groEL, dnaK, clpBCP), fatty acid desaturases (desA, desC), oxidative damage protection (recA, ruvB, sod, pds, crtBH, rub), and regulatory proteins including sigma factors and histidine kinases. These genes are strongly upregulated at the minimum temperature during daytime, indicating prioritization of stress response under light-exposed cold conditions.
**Behavioral:** Genes in this cluster are upregulated during daytime at cold temperatures, reflecting activation of stress response mechanisms.
**Notes:** High confidence due to clear upregulation and known stress-related gene functions.

### Cluster Cluster D | up | medium
**Name:** Prochlorococcus cluster D (up, cold stress response night)
**Enrichment:** stress response and metabolism (p=None, sig=False)
**Functional:** Cluster D contains genes related to global stress response, including chaperones, DNA repair, protein synthesis, glycogen degradation (glgP), amino acid synthesis (glyA, serA, leuA), translation initiation factors (infABC), and nitrogen acquisition genes. These genes are upregulated at the minimum temperature during both daytime and nighttime, reflecting a sustained stress response and increased demand for protein synthesis and energy mobilization under cold stress.
**Behavioral:** Genes in this cluster show upregulation at cold temperatures during both day and night, indicating a prolonged stress response.
**Notes:** Moderate to high confidence based on gene functions and expression patterns.

### Cluster Cluster E | up | high
**Name:** Prochlorococcus cluster E (up, core metabolism night)
**Enrichment:** core metabolism and cell cycle (p=None, sig=False)
**Functional:** Cluster E includes genes related to catabolic consumption (cyoB, ndhD), DNA replication (dnaA, nrdJ, gyrB), cell division (ftsZYQ), and the pentose phosphate pathway (tal, gnd, zwf). These genes are upregulated at nighttime across the thermal gradient, representing essential cellular and metabolic processes maintained under all temperature conditions.
**Behavioral:** Genes in this cluster show consistent upregulation at night, reflecting core metabolic and cell cycle functions.
**Notes:** High confidence as these are essential nighttime processes.

## Bernstein 2017 / bp1_light_clusters

### Cluster 0 | up | high
**Name:** Thermosynechococcus vestitus BP-1 cluster 0 (up, primary productivity)
**Enrichment:** primary productivity and carbon metabolism (p=0.05, sig=True)
**Functional:** This cluster contains 79 genes from T. elongatus BP-1 that show increased transcript abundance with increasing irradiance and decreased abundance with increasing oxygen tension. Genes include those involved in organic acid synthesis such as acetyl-CoA synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA). Also included are genes for synthesis and export of larger biomolecules like exopolysaccharides (exoD) and sucrose metabolism (spsA, sps). These genes are implicated in the production and export of reduced carbon compounds that support heterotrophic growth.
**Behavioral:** Genes in this cluster show a tent-shaped expression pattern with peak transcript abundance at intermediate irradiance and decreased expression with increasing oxygen. This pattern corresponds to increased primary productivity and carbon export supporting heterotrophic partner growth.
**Notes:** High confidence based on strong correlation with irradiance and oxygen treatments and functional enrichment.

### Cluster 1 | up | high
**Name:** Thermosynechococcus vestitus BP-1 cluster 1 (up, ROS detoxification)
**Enrichment:** oxidative stress response (p=0.05, sig=True)
**Functional:** This cluster contains 914 genes showing increased transcript abundance with increasing irradiance and oxygen tension. It includes genes involved in reactive oxygen species (ROS) detoxification such as alkylhydroperoxidase (ahpD), peroxiredoxins (BCP-A family, 2-Cys family), and Mn-superoxide dismutase (sodB). These genes are involved in mitigating oxidative stress induced by heterotrophic partnership.
**Behavioral:** Genes show increased expression with increasing irradiance and oxygen, corresponding to increased ROS detoxification activity in response to heterotrophic partnership.
**Notes:** High confidence based on consistent expression patterns and functional annotation.

### Cluster 2 | mixed | medium
**Name:** Thermosynechococcus vestitus BP-1 cluster 2 (mixed, nitrogen and vitamin exchange)
**Enrichment:** vitamin B12 and amino acid metabolism (p=0.05, sig=True)
**Functional:** This cluster contains 255 genes showing mixed expression patterns with respect to irradiance and oxygen. It includes genes involved in vitamin B12 biosynthesis and salvage (cobWNT, cobOQDPC), methionine metabolism (metH, metE), and amino acid biosynthesis and degradation. The cluster suggests coordinated exchange of vitamins and amino acids between T. elongatus and M. ruber.
**Behavioral:** Genes show complex expression patterns, with some increasing and others decreasing with irradiance and oxygen, reflecting metabolic coupling for vitamin and amino acid exchange.
**Notes:** Moderate confidence due to complexity of expression patterns and functional diversity.

### Cluster 3 | down | medium
**Name:** Thermosynechococcus vestitus BP-1 cluster 3 (down, growth related)
**Enrichment:** photosynthesis and growth regulation (p=0.05, sig=True)
**Functional:** This cluster contains 201 genes showing decreased transcript abundance with increasing irradiance and oxygen tension. Genes include those involved in photosystem II main subunits, type II restriction-modification system, and ornithine biosynthesis (argJ). These genes are associated with growth and photosynthetic functions that decrease under stress conditions.
**Behavioral:** Genes show decreased expression with increasing irradiance and oxygen, corresponding to reduced growth and photosynthetic activity under high oxygen stress.
**Notes:** Moderate confidence due to less detailed description in text.

### Cluster 4 | down | high
**Name:** Thermosynechococcus vestitus BP-1 cluster 4 (down, primary productivity)
**Enrichment:** nitrogen metabolism and oxidative stress response (p=0.05, sig=True)
**Functional:** This cluster contains 971 genes showing decreased transcript abundance with increasing irradiance and increased abundance with increasing oxygen tension. Genes include those involved in nitrogen uptake and metabolism such as nitrate uptake system (nrtABD), glutamine synthetase (glnA), nitrogen regulatory protein (glnB), assimilatory ferredoxin-nitrate reductase (nirA), and glutamate symporter (gltS). Also included are genes involved in ROS detoxification and photosystem I stabilization.
**Behavioral:** Genes in this cluster show decreased expression with increasing irradiance and increased expression with increasing oxygen tension, corresponding to decreased growth and photosynthesis rates under high oxygen.
**Notes:** High confidence due to significant enrichment and consistent physiological correlation.

## Bernstein 2017 / bp1_oxygen_clusters

### Cluster 0 | up | high
**Name:** Thermosynechococcus vestitus BP-1 cluster 0 (up, primary productivity and carbon exchange)
**Enrichment:** primary productivity and carbon exchange (p=0.05, sig=True)
**Functional:** This cluster contains genes from both T. elongatus and M. ruber that show increased transcript abundance with increasing irradiance and decreased abundance with increasing oxygen tension. Key T. elongatus genes include acetyl coenzyme A synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA), involved in organic acid synthesis. Genes involved in synthesis and export of larger biomolecules such as exopolysaccharide synthesis (exoD) and sucrose metabolism (spsA, sps) are also included. Corresponding M. ruber genes encode enzymes for uptake and metabolism of these compounds, including acetyl-CoA synthetase, cytochrome c oxidase, monosaccharide uptake systems, and branched-chain amino acid uptake.
**Behavioral:** Genes in this cluster show increased expression with irradiance and decreased expression with oxygen, reflecting coordinated metabolic coupling for carbon production and uptake between the cyanobacterium and heterotroph.
**Notes:** High confidence based on strong coexpression patterns and functional enrichment.

### Cluster 1 | down | medium
**Name:** Thermosynechococcus vestitus BP-1 cluster 1 (down, nitrogen and amino acid exchange)
**Enrichment:** nitrogen and amino acid exchange (p=0.05, sig=True)
**Functional:** This cluster includes genes from both species that show decreased expression with increasing irradiance and increased expression with increasing oxygen tension. T. elongatus genes include nitrate uptake system (nrtABD), glutamine synthetase (glnA), nitrogen regulatory protein (glnB), assimilatory ferredoxin-nitrate reductase (nirA), and glutamate symporter (gltS). M. ruber genes include glutamate dehydrogenase (gdhA) and amino acid uptake systems. Some M. ruber nitrogen acquisition genes such as glnA and glutamate synthase (gltB) show opposite expression patterns, suggesting exchange of glutamine and glutamate between species.
**Behavioral:** Genes in this cluster reflect nitrogen metabolism and amino acid exchange, with coordinated expression patterns indicating metabolic coupling for nitrogen assimilation and transfer.
**Notes:** Moderate confidence; some genes show opposite patterns indicating complex regulation.

### Cluster 2 | mixed | medium
**Name:** Thermosynechococcus vestitus BP-1 cluster 2 (mixed, methionine and vitamin B12 exchange)
**Enrichment:** methionine and vitamin B12 exchange (p=0.05, sig=True)
**Functional:** This cluster contains genes with mixed expression patterns related to methionine and vitamin B12 metabolism. M. ruber is B12 auxotrophic and expresses vitamin B12 uptake genes (btuCD) that decrease with irradiance, while T. elongatus expresses B12 biosynthesis genes (cobWNT, cobOQDPC) that also decrease with irradiance. Both species express vitamin B12-dependent and independent homocysteine methyltransferases (metH and metE) with opposing expression patterns, suggesting exchange of methionine and vitamin B12. Genes involved in methionine degradation and salvage show varied expression across clusters.
**Behavioral:** Genes in this cluster indicate coordinated exchange and metabolism of methionine and vitamin B12 between the cyanobacterium and heterotroph, reflecting metabolic interdependence.
**Notes:** Moderate confidence; complex and opposing expression patterns require cautious interpretation.

### Cluster 3 | up | high
**Name:** Thermosynechococcus vestitus BP-1 cluster 3 (up, oxidative stress responses)
**Enrichment:** oxidative stress responses (p=0.05, sig=True)
**Functional:** This cluster includes genes from both species involved in reactive oxygen species (ROS) detoxification and oxidative stress responses. T. elongatus genes include NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys family peroxiredoxins, periplasmic peroxiredoxin (prxQ-B2), and Mn-superoxide dismutase (sodB). M. ruber genes include an H2O2-responsive transcriptional regulator (oxyR) and manganese catalase, as well as peroxidases and superoxide dismutase genes. Expression of these genes generally increases with irradiance and oxygen tension, reflecting increased oxidative stress and protective responses during heterotrophic partnership.
**Behavioral:** Genes in this cluster show increased expression with increasing irradiance and oxygen, indicating enhanced oxidative stress protection and ROS detoxification in response to heterotrophic partnership.
**Notes:** High confidence based on consistent expression patterns and functional annotation.

### Cluster 4 | down | medium
**Name:** Thermosynechococcus vestitus BP-1 cluster 4 (down, electron transfer and ROS potentiators)
**Enrichment:** electron transfer and ROS potentiators (p=0.05, sig=True)
**Functional:** This cluster contains genes primarily from M. ruber associated with electron transfer processes that potentiate reactive oxygen species (ROS). These include electron transfer flavoprotein subunits (fixAB), NADH-dehydrogenase, and components of NADH-quinone oxidoreductase (nuoDFGHIJKN). Expression of these genes decreases with increasing oxygen tension and increases with specific growth and photosynthesis rates, suggesting modulation of electron transfer activity in response to environmental conditions and partnership.
**Behavioral:** Genes in this cluster show decreased expression with increasing oxygen, reflecting downregulation of electron transfer processes that may generate ROS, potentially as an acclimation to oxidative stress during partnership.
**Notes:** Moderate confidence; functional roles inferred from gene annotations and expression patterns.

## Bernstein 2017 / mruber_light_clusters

### Cluster 0 | up | high
**Name:** M. ruber cluster 0 (up, primary productivity)
**Enrichment:** carbon metabolism and transport (p=0.05, sig=True)
**Functional:** This cluster contains 857 genes from M. ruber that show increased transcript abundance with increasing irradiance and decreased abundance with increasing oxygen tension. Genes include those involved in carbon uptake and metabolism such as acetyl-CoA synthetase (acs), cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), xylose isomerase (xylA), and branched-chain amino acid uptake (livGH). These genes are coexpressed with T. elongatus genes involved in organic acid synthesis and export of biomolecules like exopolysaccharides and sucrose, indicating coordinated metabolic coupling for carbon exchange.
**Behavioral:** Genes in this cluster increase expression with irradiance and decrease with oxygen, reflecting enhanced primary productivity and carbon exchange supporting heterotrophic growth.
**Notes:** High confidence due to strong coexpression patterns and functional enrichment.

### Cluster 1 | down | medium
**Name:** M. ruber cluster 1 (down, nitrogen metabolism)
**Enrichment:** nitrogen metabolism (p=0.05, sig=True)
**Functional:** This cluster of 324 genes includes M. ruber genes involved in nitrogen acquisition and amino acid metabolism, such as glutamine synthetase (glnA) and glutamate synthase (gltB), which show decreased expression with increasing irradiance and increased expression with increasing oxygen tension. These genes exhibit opposite expression patterns to T. elongatus nitrogen metabolism genes, suggesting a coordinated exchange of reduced nitrogen compounds like glutamine and glutamate between the species.
**Behavioral:** Genes decrease expression with irradiance and increase with oxygen, indicating nitrogen metabolism adjustments in response to environmental and partnership cues.
**Notes:** Moderate confidence; gene functions and expression patterns support nitrogen exchange inference.

### Cluster 2 | mixed | high
**Name:** M. ruber cluster 2 (mixed, vitamin B12 and amino acid exchange)
**Enrichment:** vitamin and amino acid metabolism (p=0.05, sig=True)
**Functional:** This cluster of 188 genes includes M. ruber genes involved in vitamin B12 uptake and methionine biosynthesis and salvage pathways. M. ruber is B12 auxotrophic and depends on T. elongatus for vitamin B12. Genes encoding vitamin B12-dependent methionine synthesis enzymes show expression patterns coordinated with T. elongatus B12 biosynthesis genes, indicating metabolic coupling. Methionine and cysteine biosynthesis genes also show coordinated expression patterns, suggesting exchange and salvage of these amino acids between species.
**Behavioral:** Genes show mixed expression patterns reflecting vitamin B12 and amino acid exchange and salvage between M. ruber and T. elongatus.
**Notes:** High confidence due to clear metabolic dependency and coexpression evidence.

### Cluster 3 | up | medium
**Name:** M. ruber cluster 3 (up, oxidative stress response)
**Enrichment:** oxidative stress response (p=0.05, sig=True)
**Functional:** This cluster of 259 genes includes M. ruber genes involved in oxidative stress responses such as the H2O2-responsive transcriptional regulator oxyR and manganese catalase. Peroxidase and superoxide dismutase genes show expression patterns increasing with oxygen tension, indicating a response to oxidative stress. Genes associated with electron transfer processes that potentiate reactive oxygen species decrease with increasing oxygen tension, suggesting modulation of redox balance in partnership.
**Behavioral:** Genes increase expression with oxygen tension, reflecting oxidative stress response and redox regulation in M. ruber during partnership.
**Notes:** Moderate confidence; gene functions and expression patterns support oxidative stress response.

### Cluster 4 | down | low
**Name:** M. ruber cluster 4 (down, growth related)
**Enrichment:** growth and metabolism (p=0.05, sig=True)
**Functional:** This cluster of 444 genes includes M. ruber genes that show decreased expression with increasing irradiance and increased expression with increasing oxygen tension, opposite to clusters associated with growth and metabolism. These genes include those involved in amino acid biosynthesis and other growth-related functions that may be downregulated as growth rates increase, reflecting metabolic adjustments to environmental conditions and partnership.
**Behavioral:** Genes decrease expression with irradiance and increase with oxygen tension, indicating modulation of growth-related functions in response to environmental and partnership cues.
**Notes:** Low to medium confidence; functional annotations are broad and expression patterns less distinct.

## Bernstein 2017 / mruber_oxygen_clusters

### Cluster 0 | up | high
**Name:** Meiothermus ruber cluster 0 (up, primary productivity and carbon exchange)
**Enrichment:** primary productivity and carbon exchange (p=0.05, sig=True)
**Functional:** This cluster contains 580 genes from M. ruber and T. elongatus that show increased transcript abundance with increasing irradiance and decreased abundance with increasing oxygen tension. It includes genes involved in organic acid synthesis and export, such as acetyl-CoA synthetase (acs), acetate kinase (ackA), lactate dehydrogenase (ldhA), and citrate synthase (gltA). Genes for synthesis and export of larger biomolecules like exopolysaccharides (exoD), sucrose synthase (spsA), and sucrose degradation (sps) are also present. M. ruber genes in this cluster encode enzymes for uptake and metabolism of cyanobacterium-derived organic carbon, including acetyl-CoA synthetase, cytochrome c oxidase (coxB), monosaccharide uptake systems (frcA, gtsAB), xylose isomerase (xylA), xylulokinase (xylB), and ABC-type multisugar uptake systems.
**Behavioral:** Genes in this cluster increase expression with irradiance and decrease with oxygen, reflecting coordinated metabolic coupling for carbon acquisition and utilization between M. ruber and T. elongatus.
**Notes:** High confidence due to clear coexpression patterns and functional enrichment consistent with carbon metabolism.

### Cluster 1 | down | medium
**Name:** Meiothermus ruber cluster 1 (down, nitrogen and amino acid exchange)
**Enrichment:** nitrogen and amino acid exchange (p=0.05, sig=True)
**Functional:** This cluster of 328 genes includes those involved in nitrogen metabolism and amino acid exchange. T. elongatus genes encoding nitrate uptake (nrtABD), glutamine synthetase (glnA), nitrogen regulatory protein (glnB), assimilatory ferredoxin-nitrate reductase (nirA), and glutamate symporter (gltS) are present. M. ruber genes include glutamate dehydrogenase (gdhA) and amino acid uptake systems. Some key M. ruber nitrogen acquisition genes (glnA, gltB) show opposite expression patterns, increasing with growth and photosynthesis rates, suggesting direct exchange of glutamine and glutamate from T. elongatus.
**Behavioral:** Genes in this cluster decrease expression with irradiance and increase with oxygen, indicating nitrogen and amino acid exchange dynamics between the two species.
**Notes:** Moderate confidence; gene expression patterns support nitrogen exchange but some genes show complex patterns.

### Cluster 2 | mixed | medium
**Name:** Meiothermus ruber cluster 2 (mixed, methionine and vitamin B12 exchange)
**Enrichment:** methionine and vitamin B12 exchange (p=0.05, sig=True)
**Functional:** This cluster of 309 genes includes those involved in methionine and vitamin B12 metabolism and exchange. M. ruber is B12 auxotrophic and expresses vitamin B12 uptake/scavenging genes (btuCD) that decrease with irradiance, concurrent with decreased T. elongatus vitamin B12 biosynthesis genes (cobWNT, cobOQDPC). Both species express vitamin B12-dependent homocysteine methyltransferase (metH) with opposing expression patterns, suggesting methionine exchange. Vitamin B12-independent homocysteine methyltransferase (metE) transcripts are also expressed. Genes for methionine degradation and salvage show varied expression patterns between species.
**Behavioral:** Genes in this cluster show mixed expression patterns, reflecting complex exchange and salvage of methionine and vitamin B12 between M. ruber and T. elongatus.
**Notes:** Moderate confidence; complex and opposing expression patterns suggest nuanced exchange mechanisms.

### Cluster 3 | up | high
**Name:** Meiothermus ruber cluster 3 (up, oxidative stress responses)
**Enrichment:** oxidative stress responses (p=0.05, sig=True)
**Functional:** This cluster of 468 genes includes those involved in oxidative stress responses and ROS detoxification. T. elongatus transcripts encoding NAD(P)H-oxygen oxidoreductase (flv4), 2-Cys family peroxiredoxins, periplasmic peroxiredoxin (prxQ-B2), and Mn-superoxide dismutase (sodB) increase with irradiance and oxygen. M. ruber contains an H2O2-responsive regulator (oxyR) adjacent to manganese catalase, grouped in clusters increasing with irradiance and oxygen. M. ruber peroxidase (bcp) and superoxide dismutase (sod2) genes respond differently to irradiance but increase with oxygen, grouping with T. elongatus genes. M. ruber genes for electron transfer processes that potentiate ROS decrease with oxygen and increase with growth and photosynthesis rates.
**Behavioral:** Genes in this cluster increase expression with irradiance and oxygen, indicating coordinated oxidative stress response and ROS detoxification between species.
**Notes:** High confidence; consistent expression patterns and functional roles support oxidative stress response.

### Cluster 4 | down | low
**Name:** Meiothermus ruber cluster 4 (down, growth repression)
**Enrichment:** growth repression and stress response (p=None, sig=False)
**Functional:** This cluster of 387 genes includes those that decrease expression with increasing irradiance and increase with oxygen tension, showing inverse patterns to growth and photosynthesis rates. Specific gene functions are not detailed in the paper, but the cluster likely includes genes involved in growth repression or stress responses that are downregulated as growth conditions improve.
**Behavioral:** Genes in this cluster decrease expression with irradiance and increase with oxygen, possibly reflecting repression under favorable growth conditions.
**Notes:** Low confidence due to lack of detailed description and functional annotation in the paper.

## Biller 2018 / mit1002_periodicity

### Cluster coculture_LD | up | high
**Name:** MIT1002 cluster 1 (up, periodic coculture LD)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 529 genes in Alteromonas macleodii MIT1002 that exhibit 24-hour periodic oscillations in transcript abundance during the diel light:dark (L:D) cycle in coculture with Prochlorococcus. These genes show rhythmic expression patterns synchronized with the daily light cycle, reflecting regulation linked to photosynthetic activity and metabolic processes. The cluster is characterized by increased transcript levels during the light phase, consistent with active metabolism and growth.
**Behavioral:** Genes in this cluster show upregulated expression during the light period of the diel cycle in coculture conditions, indicating synchronization with the daily light-driven metabolic rhythms.
**Notes:** High confidence based on clear periodicity and replicate concordance during diel L:D cycle in coculture.

### Cluster coculture_LD+coculture_darkness | mixed | low
**Name:** MIT1002 cluster 2 (mixed, periodic coculture LD+darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster contains 1 gene from Alteromonas macleodii MIT1002 that exhibits 24-hour periodic oscillations in transcript abundance both during the diel L:D cycle and under extended darkness in coculture with Prochlorococcus. The gene maintains rhythmic expression despite the absence of light, suggesting regulation influenced by biotic interactions or internal metabolic cues rather than light alone.
**Behavioral:** The gene in this cluster shows periodic oscillations in transcript levels under both light:dark and extended darkness conditions in coculture, indicating sustained rhythmic regulation independent of light.
**Notes:** Low confidence due to single gene representation; periodicity detected but biological significance uncertain.

### Cluster coculture_darkness | down | low
**Name:** MIT1002 cluster 3 (down, periodic coculture darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 1 gene from Alteromonas macleodii MIT1002 that exhibits 24-hour periodic oscillations in transcript abundance only under extended darkness in coculture with Prochlorococcus. The gene's expression is rhythmic in the absence of light, possibly reflecting metabolic or stress-related regulation specific to dark conditions.
**Behavioral:** The gene shows periodic transcript oscillations exclusively during extended darkness in coculture, indicating a dark-specific rhythmic expression pattern.
**Notes:** Low confidence due to single gene and limited data; biological role unclear.

## Biller 2018 / natl2a_darkness_survival

### Cluster darkness_axenic+darkness_coculture | mixed | high
**Name:** NATL2A cluster 1 (mixed, darkness axenic+coculture)
**Enrichment:** Photosynthesis and energy metabolism (p=0.021, sig=True)
**Functional:** This cluster includes 96 genes with transcripts present in both axenic and coculture conditions during extended darkness (72-144h). Genes in this cluster show mixed directionality in expression changes, reflecting a combination of responses to extended darkness that are common to both culture conditions. Functional roles include photosynthesis, energy metabolism, and stress response, with genes such as PMN2A_RS00075, PMN2A_RS00090, and PMN2A_RS00180 as examples. The cluster reflects a core set of genes maintaining activity or regulation under darkness regardless of heterotrophic presence.
**Behavioral:** Genes in this cluster maintain transcript presence during late extended darkness in both axenic and coculture conditions, showing mixed expression changes without a clear uniform up- or down-regulation pattern.
**Notes:** High confidence due to consistent transcript presence in both conditions and supporting gene examples.

### Cluster darkness_axenic+unique_axenic | down | medium
**Name:** NATL2A cluster 2 (down, darkness axenic+unique axenic)
**Enrichment:** Biosynthesis and nucleotide metabolism (p=0.016, sig=True)
**Functional:** This cluster comprises 76 genes with transcripts present in axenic cultures during extended darkness, including those unique to axenic conditions. These genes generally show downregulation under extended darkness, reflecting metabolic downshifts such as decreased biosynthesis, nucleotide metabolism, and energy generation. Examples include PMN2A_RS00025, PMN2A_RS00125, and PMN2A_RS00885. The cluster highlights the transcriptional response of axenic Prochlorococcus to energy deprivation and stress without heterotrophic support.
**Behavioral:** Genes in this cluster are downregulated during extended darkness in axenic cultures, indicating reduced metabolic activity and stress response in the absence of heterotrophs.
**Notes:** Moderate to high confidence based on differential expression and pathway enrichment analyses.

### Cluster darkness_coculture+unique_coculture | up | high
**Name:** NATL2A cluster 3 (up, darkness coculture+unique coculture)
**Enrichment:** Organic compound degradation and salvage pathways (p=0.006, sig=True)
**Functional:** This cluster contains 87 genes with transcripts present in coculture conditions during extended darkness, including those unique to coculture. These genes tend to be upregulated or maintained in expression, reflecting sustained metabolic activity supported by heterotrophic interactions. Functional categories include organic compound degradation, amino acid metabolism, and salvage pathways. Representative genes include PMN2A_RS00520, PMN2A_RS00570, and PMN2A_RS00890. The cluster suggests mixotrophic metabolism enabled by heterotroph-derived substrates during darkness.
**Behavioral:** Genes in this cluster are upregulated or stable during extended darkness in coculture, indicating continued biosynthetic and metabolic activity facilitated by heterotrophic support.
**Notes:** High confidence supported by transcriptomic data and functional enrichment consistent with mixotrophy.

## Biller 2018 / natl2a_periodicity

### Cluster axenic_LD | up | high
**Name:** NATL2A cluster 1 (up, axenic_LD)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 68 genes that exhibit 24-h periodicity exclusively in axenic light:dark (L:D) conditions. Genes such as NATL2_00231, NATL2_01421, and NATL2_01431 are representative. These genes are involved in typical diel metabolic and biosynthetic processes that are synchronized with the light cycle in the absence of heterotrophs.
**Behavioral:** Genes in this cluster show increased transcript abundance peaking during the day under axenic L:D conditions, reflecting light-driven regulation.
**Notes:** Cluster defined by periodicity only in axenic L:D; direction and timing inferred from diel expression patterns.

### Cluster axenic_LD+axenic_darkness+coculture_darkness | mixed | low
**Name:** NATL2A cluster 2 (mixed, axenic_LD+axenic_darkness+coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster contains a single gene (NATL2_17641) that exhibits 24-h periodicity in axenic L:D, axenic extended darkness, and coculture extended darkness conditions. The gene's function is not specifically described in the paper.
**Behavioral:** The gene shows mixed expression patterns across conditions, with periodicity maintained in both light and some dark conditions.
**Notes:** Single gene cluster with limited description; periodicity across multiple conditions.

### Cluster axenic_LD+coculture_LD | up | high
**Name:** NATL2A cluster 3 (up, axenic_LD+coculture_LD)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This large cluster of 902 genes exhibits 24-h periodicity in both axenic and cocultured L:D conditions. Representative genes include NATL2_00001, NATL2_00011, and NATL2_00031. These genes are involved in core diel metabolic processes, including photosynthesis, carbon fixation, and biosynthesis, synchronized with the light cycle.
**Behavioral:** Genes peak during the day in both axenic and coculture L:D conditions, reflecting light-driven diel regulation.
**Notes:** Cluster represents core diel genes oscillating in both culture types under L:D.

### Cluster axenic_LD+coculture_LD+axenic_darkness | up | medium
**Name:** NATL2A cluster 4 (up, axenic_LD+coculture_LD+axenic_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 38 genes showing 24-h periodicity in axenic L:D, coculture L:D, and axenic extended darkness conditions. Representative genes include NATL2_00951, NATL2_01161, and NATL2_02051. These genes may be involved in metabolic pathways that maintain some rhythmicity even in axenic extended darkness.
**Behavioral:** Genes maintain diel oscillations in light conditions and some periodicity in axenic darkness, indicating partial persistence of rhythmic expression.
**Notes:** Periodic expression extends into axenic darkness for this subset.

### Cluster axenic_LD+coculture_LD+axenic_darkness+coculture_darkness | up | medium
**Name:** NATL2A cluster 5 (up, axenic_LD+coculture_LD+axenic_darkness+coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 82 genes exhibits 24-h periodicity across all four conditions: axenic L:D, coculture L:D, axenic extended darkness, and coculture extended darkness. Representative genes include NATL2_00601, NATL2_00891, and NATL2_01071. These genes likely represent core metabolic and regulatory functions that maintain rhythmicity even under extended darkness and in coculture.
**Behavioral:** Genes show sustained oscillations across all conditions, indicating robust periodic regulation.
**Notes:** Strong periodicity across all conditions suggests key regulatory roles.

### Cluster axenic_LD+coculture_LD+coculture_darkness | up | medium
**Name:** NATL2A cluster 6 (up, axenic_LD+coculture_LD+coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster contains 377 genes with 24-h periodicity in axenic L:D, coculture L:D, and coculture extended darkness. Representative genes include NATL2_00021, NATL2_00061, and NATL2_00091. These genes may be involved in metabolic pathways supported by heterotroph interactions during darkness.
**Behavioral:** Genes maintain diel oscillations in light and coculture darkness, suggesting heterotroph-mediated maintenance of rhythmicity.
**Notes:** Periodic expression in coculture darkness indicates heterotroph influence.

### Cluster axenic_LD+coculture_darkness | mixed | low
**Name:** NATL2A cluster 7 (mixed, axenic_LD+coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This small cluster of 5 genes shows 24-h periodicity in axenic L:D and coculture extended darkness conditions. Representative genes include NATL2_07321, NATL2_07931, and NATL2_08291. These genes may have mixed expression patterns reflecting adaptation to heterotroph presence during darkness.
**Behavioral:** Genes show mixed oscillation patterns, with periodicity maintained in axenic L:D and coculture darkness but not consistently across all conditions.
**Notes:** Limited data; mixed directionality and timing.

### Cluster coculture_LD | up | medium
**Name:** NATL2A cluster 8 (up, coculture_LD)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 328 genes exhibiting 24-h periodicity exclusively in coculture L:D conditions. Representative genes include NATL2_00051, NATL2_00071, and NATL2_00251. These genes may be involved in metabolic pathways enhanced or stabilized by heterotroph interactions under light conditions.
**Behavioral:** Genes peak during the day in coculture L:D but do not show periodicity in axenic cultures, indicating heterotroph influence on diel regulation.
**Notes:** Periodic only in coculture L:D, suggesting heterotroph-driven regulation.

### Cluster coculture_LD+axenic_darkness | mixed | low
**Name:** NATL2A cluster 9 (mixed, coculture_LD+axenic_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster contains 2 genes showing 24-h periodicity in coculture L:D and axenic extended darkness. Representative genes include NATL2_13951 and NATL2_18631. These genes may have mixed expression patterns reflecting partial rhythmicity in these conditions.
**Behavioral:** Genes show mixed oscillation patterns, with periodicity in coculture L:D and axenic darkness but not consistently across all conditions.
**Notes:** Very small cluster with limited data.

### Cluster coculture_LD+coculture_darkness | up | medium
**Name:** NATL2A cluster 10 (up, coculture_LD+coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 71 genes exhibiting 24-h periodicity in coculture L:D and coculture extended darkness. Representative genes include NATL2_00131, NATL2_00461, and NATL2_00621. These genes likely participate in metabolic pathways supported by heterotroph interactions during darkness.
**Behavioral:** Genes maintain diel oscillations in coculture light and darkness, indicating heterotroph-mediated maintenance of rhythmicity.
**Notes:** Periodic expression in coculture darkness supports heterotroph influence.

### Cluster coculture_darkness | mixed | low
**Name:** NATL2A cluster 11 (mixed, coculture_darkness)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This small cluster of 5 genes shows 24-h periodicity exclusively in coculture extended darkness. Representative genes include NATL2_00291, NATL2_01401, and NATL2_08455. These genes may be involved in metabolic functions maintained by heterotroph interactions during darkness.
**Behavioral:** Genes show oscillations only in coculture darkness, indicating heterotroph-driven rhythmicity under dark stress.
**Notes:** Limited data; periodicity only in coculture darkness.

### Cluster not_periodic | not_described | high
**Name:** NATL2A cluster 12 (not_periodic)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 202 genes that do not exhibit significant 24-h periodicity in any of the four conditions. Representative genes include NATL2_00181, NATL2_00201, and NATL2_00321. These genes may be involved in constitutive or stress-related functions not tied to diel cycles.
**Behavioral:** Genes show no significant oscillations across any condition, indicating non-periodic expression.
**Notes:** Non-periodic genes with stable or variable expression.

## Coe 2024 / supp_table_3_darktolerant_clusters

### Cluster 1 | mixed | medium
**Name:** Prochlorococcus NATL2A cluster 1 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster includes 183 genes with expression patterns reflecting the dark-tolerant phenotype of Prochlorococcus NATL2A. Key genes include Yfr1 and Yfr106, which are small RNAs involved in regulatory processes. The cluster shows a mixed direction of expression changes, with some genes upregulated and others downregulated in dark-tolerant cells relative to parental cells. Functional roles relate to metabolic adjustments supporting survival in extended darkness, including shifts in carbon metabolism and stress responses.
**Behavioral:** Genes in this cluster exhibit altered diel expression patterns in dark-tolerant cells, with changes in timing and amplitude reflecting metabolic reprogramming for dark survival.
**Notes:** Cluster composition inferred from gene examples and overall study context; specific enrichment not detailed.

### Cluster 10 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 10 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster contains 173 genes including Yfr13 and other genes with locus tags. The expression pattern is mixed, with genes showing both up- and down-regulation in dark-tolerant cells. The cluster likely represents genes involved in metabolic and regulatory processes adapting to extended darkness, but specific functions are not detailed in the paper.
**Behavioral:** Genes show altered diel expression with mixed directionality, reflecting complex metabolic adjustments in dark-tolerant cells.
**Notes:** No specific functional enrichment or detailed gene roles provided.

### Cluster 11 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 11 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 11 includes 93 genes with locus tag samples. The expression pattern is mixed, with genes variably regulated in dark-tolerant cells. Functional roles are not explicitly described, but likely relate to metabolic and regulatory changes supporting dark tolerance.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with metabolic adaptation to darkness.
**Notes:** No detailed functional or enrichment data available.

### Cluster 12 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 12 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 94 genes includes locus tag samples and shows mixed expression changes in dark-tolerant cells. Specific functional roles are not described but likely involve metabolic and regulatory processes related to dark adaptation.
**Behavioral:** Genes show variable diel expression changes consistent with metabolic shifts in dark-tolerant cells.
**Notes:** No explicit functional or enrichment details provided.

### Cluster 13 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 13 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 13 contains 140 genes with locus tag samples. Expression changes are mixed in dark-tolerant cells, likely reflecting metabolic and regulatory adaptations to extended darkness, though specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

### Cluster 14 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 14 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 169 genes includes Yfr10-like and Yfr2 small RNAs and locus tag samples. Expression changes are mixed, likely involving regulatory and metabolic genes contributing to dark tolerance, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression changes consistent with metabolic and regulatory adaptation to darkness.
**Notes:** No explicit functional or enrichment data.

### Cluster 15 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 15 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 15 includes 130 genes with Yfr10-like and Yfr22 small RNAs and locus tag samples. Expression is mixed, likely reflecting regulatory and metabolic changes supporting dark tolerance, but specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

### Cluster 2 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 2 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 133 genes includes Yfr2 small RNA and locus tag samples. Expression changes are mixed in dark-tolerant cells, likely involving metabolic and regulatory genes related to dark adaptation, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression changes consistent with metabolic shifts in dark-tolerant cells.
**Notes:** No explicit functional or enrichment data.

### Cluster 3 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 3 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 3 contains 130 genes with locus tag samples. Expression changes are mixed in dark-tolerant cells, likely reflecting metabolic and regulatory adaptations to extended darkness, though specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

### Cluster 4 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 4 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 169 genes includes AtpT and Yfr107 small RNAs and locus tag samples. Expression changes are mixed, likely involving regulatory and metabolic genes contributing to dark tolerance, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression changes consistent with metabolic and regulatory adaptation to darkness.
**Notes:** No explicit functional or enrichment data.

### Cluster 5 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 5 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 5 includes 171 genes with locus tag samples. Expression is mixed, likely reflecting regulatory and metabolic changes supporting dark tolerance, but specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

### Cluster 6 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 6 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 181 genes includes Yfr103 small RNAs and locus tag samples. Expression changes are mixed, likely involving metabolic and regulatory genes related to dark adaptation, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression changes consistent with metabolic shifts in dark-tolerant cells.
**Notes:** No explicit functional or enrichment data.

### Cluster 7 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 7 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 7 contains 118 genes with Yfr10-like and Yfr108 small RNAs and locus tag samples. Expression changes are mixed in dark-tolerant cells, likely reflecting metabolic and regulatory adaptations to extended darkness, though specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

### Cluster 8 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 8 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** This cluster of 83 genes includes locus tag samples. Expression changes are mixed, likely involving metabolic and regulatory genes related to dark adaptation, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression changes consistent with metabolic shifts in dark-tolerant cells.
**Notes:** No explicit functional or enrichment data.

### Cluster 9 | mixed | low
**Name:** Prochlorococcus NATL2A cluster 9 (mixed, dark tolerance metabolism)
**Enrichment:** not described in paper (p=None, sig=False)
**Functional:** Cluster 9 includes 114 genes with Yfr22 and Yfr23 small RNAs. Expression is mixed, likely reflecting regulatory and metabolic changes supporting dark tolerance, but specific functions are not detailed.
**Behavioral:** Genes exhibit mixed diel expression changes consistent with dark tolerance adaptations.
**Notes:** No detailed functional or enrichment information.

## Coe 2024 / supp_table_3_parental_clusters

### Cluster 1 | mixed | medium
**Name:** Prochlorococcus cluster 1 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 97 genes with diel expression patterns featuring Yfr106, Yfr2_2, and Yfr2_4 as representative genes. The genes show mixed directionality in expression over the diel cycle, with no clear peak time or period described. Functional roles are not explicitly detailed in the paper, but the presence of Yfr family RNAs suggests involvement in regulatory processes related to diel rhythms.
**Behavioral:** Genes in this cluster exhibit mixed expression patterns across the diel cycle without a clear peak or trough, indicating complex regulation possibly linked to circadian or metabolic processes.
**Notes:** Cluster composition inferred from gene samples; functional roles not explicitly described in the paper.

### Cluster 10 | up | medium
**Name:** Prochlorococcus cluster 10 (up, ATP synthesis)
**Enrichment:** ATP synthesis (p=None, sig=False)
**Functional:** This cluster contains 70 genes including AtpT, which is involved in ATP synthesis. Genes in this cluster show upregulated expression, particularly related to energy metabolism and ATP synthesis pathways. The cluster likely reflects diel regulation of energy production processes.
**Behavioral:** Genes peak in expression during periods requiring increased ATP synthesis, consistent with diel energy demands.
**Notes:** AtpT presence suggests ATP synthesis function; exact timing not detailed.

### Cluster 11 | mixed | medium
**Name:** Prochlorococcus cluster 11 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 205 genes with representative genes Yfr1, Yfr2_3, and cds-WP_009788866.1. The genes show mixed diel expression patterns, possibly involved in regulatory and metabolic functions related to diel cycles.
**Behavioral:** Expression patterns are mixed without a clear peak, indicating complex regulation possibly linked to circadian rhythms or metabolic adjustments.
**Notes:** No explicit functional description; inference based on gene samples.

### Cluster 12 | mixed | low
**Name:** Prochlorococcus cluster 12 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster has 144 genes including cds-PMN2A_RS05260, cds-PMN2A_RS09840, and cds-WP_011293542.1. Functional roles are not described in the paper. Expression patterns are mixed over the diel cycle.
**Behavioral:** Genes show mixed diel expression without clear peak or trough, suggesting diverse or complex regulation.
**Notes:** No functional or behavioral description provided.

### Cluster 13 | mixed | low
**Name:** Prochlorococcus cluster 13 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster contains 165 genes including cds-WP_011293552.1, cds-WP_011293553.1, and cds-WP_011293557.1. Functional roles and diel expression patterns are not described in the paper.
**Behavioral:** Genes exhibit mixed expression patterns over the diel cycle with no clear peak or directionality.
**Notes:** No functional or behavioral details provided.

### Cluster 14 | mixed | low
**Name:** Prochlorococcus cluster 14 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 89 genes such as cds-WP_011293597.1, cds-WP_011293636.1, and cds-WP_011293647.1. Functional roles and diel expression patterns are not described in the paper.
**Behavioral:** Genes show mixed diel expression without clear timing or directionality.
**Notes:** No functional or behavioral description provided.

### Cluster 15 | mixed | medium
**Name:** Prochlorococcus cluster 15 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster has 132 genes including Yfr10.like_1, Yfr22_1, and cds-WP_011125327.1-2. The genes show mixed diel expression patterns, likely involved in regulatory functions related to diel rhythms.
**Behavioral:** Genes exhibit mixed expression patterns over the diel cycle without a clear peak or trough.
**Notes:** No explicit functional description; inference based on gene samples.

### Cluster 2 | mixed | medium
**Name:** Prochlorococcus cluster 2 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 110 genes with Yfr10.like_4, cds-WP_011293588.1, and cds-WP_011293608.1 as representative genes. Functional roles are not described, but Yfr RNAs suggest regulatory functions related to diel cycles.
**Behavioral:** Genes show mixed diel expression patterns without clear peak or trough.
**Notes:** No detailed functional description available.

### Cluster 3 | mixed | low
**Name:** Prochlorococcus cluster 3 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster contains 245 genes including cds-WP_011293551.1, cds-WP_011293554.1, and cds-WP_011293576.1. Functional roles and diel expression patterns are not described.
**Behavioral:** Genes exhibit mixed diel expression without clear timing or directionality.
**Notes:** No functional or behavioral details provided.

### Cluster 4 | mixed | low
**Name:** Prochlorococcus cluster 4 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 100 genes such as cds-WP_011293580.1, cds-WP_011293592.1, and cds-WP_011293599.1. Functional roles and diel expression patterns are not described.
**Behavioral:** Genes show mixed diel expression without clear peak or trough.
**Notes:** No functional or behavioral description provided.

### Cluster 5 | mixed | low
**Name:** Prochlorococcus cluster 5 (mixed, unknown)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster has 65 genes including cds-WP_011293563.1, cds-WP_011293645.1, and cds-WP_011293646.1. Functional roles and diel expression patterns are not described.
**Behavioral:** Genes exhibit mixed diel expression without clear timing or directionality.
**Notes:** No functional or behavioral details provided.

### Cluster 6 | mixed | medium
**Name:** Prochlorococcus cluster 6 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster contains 62 genes including Yfr10.like_2, Yfr103_1, and Yfr103_2. The genes likely have regulatory roles related to diel rhythms, but specific functions are not described.
**Behavioral:** Genes show mixed diel expression patterns without clear peak or trough.
**Notes:** No explicit functional description; inferred from gene samples.

### Cluster 7 | mixed | medium
**Name:** Prochlorococcus cluster 7 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 202 genes with Yfr108, cds-WP_011293543.1, and cds-WP_011293544.1 as representative genes. Functional roles are not described, but Yfr RNAs suggest regulatory functions related to diel cycles.
**Behavioral:** Genes exhibit mixed diel expression patterns without clear peak or trough.
**Notes:** No detailed functional description available.

### Cluster 8 | mixed | medium
**Name:** Prochlorococcus cluster 8 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster contains 174 genes including Yfr107, cds-WP_011293537.1, and cds-WP_011293538.1. Functional roles are not described, but Yfr RNAs suggest regulatory functions related to diel rhythms.
**Behavioral:** Genes show mixed diel expression patterns without clear peak or trough.
**Notes:** No explicit functional description; inferred from gene samples.

### Cluster 9 | mixed | medium
**Name:** Prochlorococcus cluster 9 (mixed, Yfr)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster includes 221 genes with Yfr13, Yfr22_2, and Yfr22_3 as representative genes. Functional roles are not described, but Yfr RNAs suggest regulatory functions related to diel cycles.
**Behavioral:** Genes exhibit mixed diel expression patterns without clear peak or trough.
**Notes:** No detailed functional description available.

## Lindell 2007 / med4_phage_transcription_groups

### Cluster 1 | up | high
**Name:** MED4 cluster 1 (up, early transient)
**Enrichment:** host stress response and metabolism (p=0.001, sig=True)
**Functional:** This cluster consists of host genes that are transiently upregulated immediately after phage infection. It includes genes involved in high-light-inducible stress response (hli), carbon metabolism (rbcLS), transcription (rpoC2, rpoD), and ribosomal proteins (rpl5, rpl6, rps8, rps11, rps17). These genes likely represent an early host stress response or phage-induced activation.
**Behavioral:** Genes in this cluster show a sharp transient peak in expression around 1 hour post infection, then decline.
**Notes:** High confidence based on microarray data and RT-PCR validation; gene functions well characterized.

### Cluster 2 | up | medium
**Name:** MED4 cluster 2 (up, late sustained)
**Enrichment:** RNA degradation, protein turnover, stress response (p=0.001, sig=True)
**Functional:** This cluster includes host genes upregulated starting about 2 hours post infection and sustained thereafter. It contains genes involved in RNA degradation and modification (rne, rnhB, dus, sun), protein turnover (clpS, AAA ATPase family gene), stress responses (umuD, phoH), and genes of unknown function. These genes may be involved in host RNA degradation and protein turnover, potentially exploited by phage for nucleotide synthesis and stress adaptation.
**Behavioral:** Genes in this cluster show increased expression starting 2 hours post infection and remain elevated through the infection period.
**Notes:** Moderate to high confidence; gene functions inferred from annotations and expression timing; some genes of unknown function included.

## Tolonen 2006 / med4_kmeans_nstarvation

### Cluster 1 | up | high
**Name:** MED4 cluster 1 (up, transport and binding)
**Enrichment:** transport and binding (p=0.01, sig=True)
**Functional:** This cluster contains the most rapidly and highly upregulated genes in MED4 during nitrogen starvation, including nitrogen transport genes such as urtA and cynA, and the nitrite permease in MIT9313 orthologs. It also includes the highly upregulated hli10 gene and other genes with putative NtcA binding sites.
**Behavioral:** Genes in this cluster are rapidly and strongly upregulated early in the nitrogen starvation time course.
**Notes:** Strong functional enrichment and early rapid upregulation support confidence in this cluster.

### Cluster 2 | up | medium
**Name:** MED4 cluster 2 (up, transport and binding)
**Enrichment:** transport and binding (p=0.12, sig=False)
**Functional:** This cluster contains genes involved in transport and binding, including a subset of hli genes that respond rapidly and highly to nitrogen starvation. It also includes regulatory genes such as sigma factors.
**Behavioral:** Genes in this cluster are rapidly upregulated during nitrogen starvation but to a lesser degree than cluster 1.
**Notes:** Functional enrichment is weaker but includes important transport and hli genes.

### Cluster 3 | up | medium
**Name:** MED4 cluster 3 (up, regulation)
**Enrichment:** regulation (p=0.07, sig=False)
**Functional:** This cluster contains regulatory genes including sigma factors and other genes involved in regulation during nitrogen starvation.
**Behavioral:** Genes in this cluster are upregulated during nitrogen starvation with moderate timing and magnitude.
**Notes:** Regulatory theme but weaker enrichment.

### Cluster 4 | up | low
**Name:** MED4 cluster 4 (up, amino acid synthesis)
**Enrichment:** amino acid synthesis (p=0.18, sig=False)
**Functional:** This cluster is enriched for genes involved in amino acid synthesis, which are upregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show upregulation during nitrogen starvation, reflecting increased amino acid synthesis activity.
**Notes:** Weak enrichment for amino acid synthesis genes.

### Cluster 5 | up | medium
**Name:** MED4 cluster 5 (up, regulation)
**Enrichment:** regulation (p=0.18, sig=False)
**Functional:** This cluster contains regulatory genes including sigma factors that are upregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster are upregulated during nitrogen starvation and may play roles in transcriptional regulation.
**Notes:** Regulatory genes present but weak enrichment.

### Cluster 6 | down | high
**Name:** MED4 cluster 6 (down, translation)
**Enrichment:** translation (p=0.001, sig=True)
**Functional:** This cluster is enriched for genes involved in translation, including ribosomal proteins, which are strongly downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show strong downregulation during nitrogen starvation, indicating reduced protein synthesis under stress.
**Notes:** Strong enrichment for translation-related genes.

### Cluster 7 | down | high
**Name:** MED4 cluster 7 (down, translation)
**Enrichment:** translation (p=2.7e-09, sig=True)
**Functional:** This cluster is also enriched for translation genes and shows strong downregulation during nitrogen starvation.
**Behavioral:** Genes in this cluster are strongly downregulated during nitrogen starvation, consistent with reduced translation activity.
**Notes:** Highly significant enrichment for translation genes.

### Cluster 8 | down | high
**Name:** MED4 cluster 8 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=2.6e-09, sig=True)
**Functional:** This cluster is enriched for photosynthesis and respiration genes, including photosystem I genes, which are downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show strong downregulation during nitrogen starvation, reflecting reduced photosynthetic activity.
**Notes:** Strong enrichment for photosynthesis genes.

### Cluster 9 | down | high
**Name:** MED4 cluster 9 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=1.4e-05, sig=True)
**Functional:** This cluster is enriched for photosynthesis and respiration genes, including ATP synthase subunits and carbon fixation genes (rbcLS), which are strongly downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show the strongest downregulation during nitrogen starvation, indicating a major reduction in photosynthetic and respiratory activity.
**Notes:** Strongest downregulated cluster with clear photosynthesis gene enrichment.

## Tolonen 2006 / mit9313_kmeans_nstarvation

### Cluster 1 | up | high
**Name:** MIT9313 cluster 1 (up, transport and binding)
**Enrichment:** transport and binding (p=0.04, sig=True)
**Functional:** This cluster contains the most rapidly and highly upregulated genes in MIT9313 during nitrogen starvation, including nitrogen transport genes such as urtA and the nitrite permease, as well as hli genes (hliS, hli7). It also includes the transcription factor ntcA and other genes with putative NtcA binding sites.
**Behavioral:** Genes in this cluster are rapidly and strongly upregulated early in the nitrogen starvation time course.
**Notes:** Strong functional enrichment and early rapid upregulation support confidence in this cluster.

### Cluster 2 | up | high
**Name:** MIT9313 cluster 2 (up, regulation)
**Enrichment:** regulation (p=0.03, sig=True)
**Functional:** This cluster contains genes involved in regulation, including two sigma factors (PMT2246 and PMT0346) and other regulatory genes. One of the sigma factors (PMT2246) has a strong NtcA binding site, suggesting indirect regulation by NtcA.
**Behavioral:** Genes in this cluster are upregulated during nitrogen starvation, but with a slower and more sustained pattern than cluster 1.
**Notes:** Functional enrichment for regulatory genes and presence of sigma factors support confidence.

### Cluster 3 | up | medium
**Name:** MIT9313 cluster 3 (up, regulation)
**Enrichment:** regulation (p=0.09, sig=False)
**Functional:** This cluster contains genes involved in regulation, but with less clear functional enrichment. It may include additional regulatory genes responding to nitrogen starvation.
**Behavioral:** Genes in this cluster show upregulation during nitrogen starvation with moderate timing and magnitude.
**Notes:** Weaker functional enrichment than cluster 2, but still regulatory theme.

### Cluster 4 | down | medium
**Name:** MIT9313 cluster 4 (down, transport and binding)
**Enrichment:** transport and binding (p=0.39, sig=False)
**Functional:** This cluster contains genes that were unchanged until being repressed only at the final time point, possibly representing genes repressed as part of a general transcriptional shutdown. It includes genes linking nitrogen and carbon metabolism such as glnB, icd, acnB, and rbcLS.
**Behavioral:** Genes in this cluster are repressed late in the nitrogen starvation time course, likely reflecting severe starvation and general shutdown.
**Notes:** Late repression and mixed functional roles suggest general shutdown rather than specific nitrogen response.

### Cluster 5 | down | high
**Name:** MIT9313 cluster 5 (down, fatty acid and phospholipid metabolism)
**Enrichment:** fatty acid and phospholipid metabolism (p=0.004, sig=True)
**Functional:** This cluster is enriched for genes involved in fatty acid and phospholipid metabolism, which are downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show downregulation during nitrogen starvation, possibly reflecting reduced lipid metabolism under stress.
**Notes:** Significant enrichment for fatty acid and phospholipid metabolism genes.

### Cluster 6 | down | high
**Name:** MIT9313 cluster 6 (down, photosynthesis and respiration)
**Enrichment:** photosynthesis and respiration (p=0.00049, sig=True)
**Functional:** This cluster contains genes involved in photosynthesis and respiration, including photosystem I and II genes and the phycoerythrin gene cpeB, which are downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster are strongly downregulated during nitrogen starvation, reflecting reduced photosynthetic activity under stress.
**Notes:** Strong enrichment and consistent downregulation of photosynthesis genes.

### Cluster 7 | down | high
**Name:** MIT9313 cluster 7 (down, translation)
**Enrichment:** translation (p=7.9e-10, sig=True)
**Functional:** This cluster is enriched for genes involved in translation, including ribosomal proteins, which are downregulated during nitrogen starvation.
**Behavioral:** Genes in this cluster show strong downregulation during nitrogen starvation, indicating reduced protein synthesis under stress.
**Notes:** Highly significant enrichment for translation-related genes.

## Wang 2014 / med4_expression_level

### Cluster -- | mixed | high
**Name:** Prochlorococcus MED4 cluster -- (mixed, expression level)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster consists of 89 genes that are unclassified in terms of expression level based on RNA-seq RPKM quartiles across 10 growth conditions. These genes do not fit into the defined expression categories (VEG, HEG, MEG, LEG, NEG) and thus represent a heterogeneous group with mixed expression patterns. Their functional roles are not specifically described in the paper.
**Behavioral:** Genes in this cluster show variable or inconsistent expression levels across the tested conditions, lacking a clear pattern of high or low expression.
**Notes:** Cluster is defined by exclusion from other expression level categories; no specific functional or behavioral description provided.

### Cluster HEG | up | high
**Name:** Prochlorococcus MED4 cluster HEG (up, expression level)
**Enrichment:** COG categories C, J, O (p=0.001, sig=True)
**Functional:** This cluster includes 291 highly expressed genes (HEG) consistently showing high transcript abundance across multiple growth conditions. These genes are significantly enriched in the core genome and are involved in essential cellular functions such as energy production and conversion, translation and ribosomal structure, and protein modification, folding, and turnover. They include ribosomal proteins, photosynthetic apparatus components, and protein folding machinery, reflecting their central metabolic roles and evolutionary conservation.
**Behavioral:** Genes in this cluster are consistently and highly expressed, evolve slowly, and are often organized in operons. They tend to have rapid mRNA turnover, which may enhance protein fidelity and cellular economy.
**Notes:** Strong evidence from RNA-seq data and functional enrichment analyses supports the characterization of this cluster.

### Cluster LEG | down | high
**Name:** Prochlorococcus MED4 cluster LEG (down, expression level)
**Enrichment:** COG category R (general function) (p=0.023, sig=True)
**Functional:** This cluster comprises 82 lowly expressed genes (LEG) that show consistently low transcript levels across the tested growth conditions. These genes are less represented in the core genome and more common in the flexible genome. Functional enrichment indicates a slight overrepresentation of genes with general function prediction (COG category R). They tend to evolve faster and have longer mRNA half-lives compared to highly expressed genes.
**Behavioral:** Genes in this cluster are expressed at low levels, evolve relatively faster, and have slower mRNA degradation rates compared to HEG and MEG clusters.
**Notes:** Expression and evolutionary data support the low expression and faster evolution; functional enrichment is modest.

### Cluster MEG | mixed | medium
**Name:** Prochlorococcus MED4 cluster MEG (mixed, expression level)
**Enrichment:** Essential genes (DEG-hit) enriched (p=0.001, sig=True)
**Functional:** This cluster contains 445 moderately expressed genes (MEG) with intermediate transcript abundance across growth conditions. These genes are enriched in the core genome relative to the flexible genome and include a significant proportion of essential genes. Functional categories are diverse, reflecting a range of cellular processes. MEG genes evolve at intermediate rates and have mRNA turnover rates between those of HEG and LEG clusters.
**Behavioral:** Genes in this cluster show moderate and relatively stable expression, slower evolution than lowly expressed genes but faster than HEG, and intermediate mRNA degradation rates.
**Notes:** Moderate confidence based on expression, essentiality, and evolutionary data.

### Cluster NEG | down | high
**Name:** Prochlorococcus MED4 cluster NEG (down, expression level)
**Enrichment:** Phage-related genes (p=None, sig=False)
**Functional:** This cluster includes 66 non-expressed genes (NEG) that show no detectable transcript levels under the tested growth conditions. These genes are predominantly found in the flexible genome and are often associated with phage-related functions. They evolve faster and have longer mRNA half-lives compared to expressed genes. Functional enrichment is low due to lack of expression and annotation.
**Behavioral:** Genes in this cluster are transcriptionally silent under tested conditions, evolve rapidly, and are mostly flexible genome or phage-related genes.
**Notes:** Expression data clearly define this cluster; functional annotation limited by lack of expression.

### Cluster VEG | mixed | medium
**Name:** Prochlorococcus MED4 cluster VEG (mixed, expression level)
**Enrichment:** not described (p=None, sig=False)
**Functional:** This cluster consists of 1081 variably expressed genes (VEG) that show fluctuating expression levels across different growth conditions. These genes are more prevalent in the flexible genome and include genes with diverse functions. They tend to evolve faster than constantly expressed genes and have variable mRNA turnover rates. This cluster reflects genes with condition-dependent or variable regulation.
**Behavioral:** Genes in this cluster exhibit variable expression patterns, faster evolutionary rates, and are less likely to be organized in operons compared to constantly expressed genes.
**Notes:** Defined by variable expression; functional and behavioral heterogeneity expected.

## Zinser 2009 / med4_diel_clusters

### Cluster 1 | up | high
**Name:** Prochlorococcus cluster 1 (up, photosynthesis)
**Enrichment:** Photosystem I and II (p=1.5e-09, sig=True)
**Functional:** Cluster 1 contains 57 genes enriched for photosystem I and photosystem II components. Genes such as psbA, psbD, and psaA show peak expression near sunrise (8.3 h), coinciding with the onset of light and photosynthetic activity. This cluster reflects the upregulation of photosynthetic machinery preparing the cell for daytime energy capture.
**Behavioral:** Genes in this cluster peak in expression shortly after dawn, increasing transcript abundance to support photosynthesis during the day.
**Notes:** Strong enrichment and clear peak timing near dawn support high confidence in this cluster's functional assignment.

### Cluster 10 | down | low
**Name:** Prochlorococcus cluster 10 (down, not described)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 10 contains 77 genes with peak expression at 0.3 h (just after dark). No significant functional enrichment or detailed description is provided in the paper.
**Behavioral:** Genes peak just after dark but lack clear functional annotation or description.
**Notes:** No functional enrichment or description available.

### Cluster 11 | down | low
**Name:** Prochlorococcus cluster 11 (down, not described)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 11 contains 91 genes peaking at 1.6 h after dark. No significant functional enrichment or detailed description is provided in the paper.
**Behavioral:** Genes peak early in the night but lack clear functional annotation or description.
**Notes:** No functional enrichment or description available.

### Cluster 12 | down | low
**Name:** Prochlorococcus cluster 12 (down, not described)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 12 contains 111 genes peaking at 3.0 h after dark. No significant functional enrichment or detailed description is provided in the paper.
**Behavioral:** Genes peak in the early night but lack clear functional annotation or description.
**Notes:** No functional enrichment or description available.

### Cluster 13 | down | high
**Name:** Prochlorococcus cluster 13 (down, ribosomal proteins)
**Enrichment:** Ribosomal proteins (p=3.7e-41, sig=True)
**Functional:** Cluster 13 includes 110 genes highly enriched for ribosomal proteins. Genes peak at 3.4 h after dark, consistent with ribosome biogenesis and protein synthesis during the night.
**Behavioral:** Genes peak early in the night, supporting ribosomal protein synthesis and preparation for translation.
**Notes:** Very strong enrichment and clear timing support high confidence.

### Cluster 14 | down | low
**Name:** Prochlorococcus cluster 14 (down, not described)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 14 contains 22 genes peaking at 4.6 h after dark. No significant functional enrichment or detailed description is provided in the paper.
**Behavioral:** Genes peak in the early night but lack clear functional annotation or description.
**Notes:** No functional enrichment or description available.

### Cluster 15 | down | medium
**Name:** Prochlorococcus cluster 15 (down, not described)
**Enrichment:** Calvin cycle and ATP synthase (p=0.071, sig=False)
**Functional:** Cluster 15 contains 107 genes peaking at 4.7 h after dark. This cluster is significantly enriched for Calvin cycle and ATP synthase genes, indicating involvement in carbon fixation and energy metabolism.
**Behavioral:** Genes peak early in the night, supporting Calvin cycle activity and ATP synthesis preparation for the day.
**Notes:** Enrichment is marginal; functional role inferred from gene annotations and timing.

### Cluster 16 | down | high
**Name:** Prochlorococcus cluster 16 (down, ATP synthase and CO2 metabolism)
**Enrichment:** ATP synthase and CO2 metabolism (p=8.4e-08, sig=True)
**Functional:** Cluster 16 contains 125 genes highly enriched for ATP synthase and CO2 metabolism genes. Genes peak at 5.5 h after dark, consistent with energy metabolism and carbon fixation processes preparing for the day.
**Behavioral:** Genes peak early in the morning, supporting ATP synthesis and CO2 metabolism at dawn.
**Notes:** Strong enrichment and timing consistent with energy metabolism support confidence.

### Cluster 17 | not_described | medium
**Name:** Prochlorococcus cluster 17 (aperiodic, expressed)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 17 contains 173 genes expressed but lacking significant periodicity. These genes do not show clear diel cycling and may represent housekeeping or constitutively expressed functions.
**Behavioral:** Genes are expressed without clear diel periodicity, maintaining steady expression levels throughout the cycle.
**Notes:** Lack of periodicity and functional enrichment limits interpretation.

### Cluster 18 | not_described | high
**Name:** Prochlorococcus cluster 18 (not expressed)
**Enrichment:** Menaquinone and ubiquinone metabolism (p=0.00045, sig=True)
**Functional:** Cluster 18 contains 180 genes not expressed under the experimental conditions. These genes include those involved in menaquinone and ubiquinone metabolism, which may be inactive under optimal growth conditions.
**Behavioral:** Genes are not expressed or expressed below detection threshold during the diel cycle.
**Notes:** Significant enrichment despite lack of expression suggests these genes are inactive under tested conditions.

### Cluster 2 | up | medium
**Name:** Prochlorococcus cluster 2 (up, not described)
**Enrichment:** not described (p=None, sig=False)
**Functional:** Cluster 2 contains 52 genes with peak expression around 9.6 h. The cluster is not significantly enriched for any functional category described in the paper.
**Behavioral:** Genes peak in the mid-morning but lack clear functional enrichment or description in the paper.
**Notes:** No significant enrichment or functional description available.

### Cluster 3 | up | high
**Name:** Prochlorococcus cluster 3 (up, photosynthesis)
**Enrichment:** Cytochrome b6/f and Photosystem II (p=0.0059, sig=True)
**Functional:** Cluster 3 includes 23 genes enriched for cytochrome b6/f complex and photosystem II components. Genes such as petA, petB, and psbD peak around 12.5 h, coinciding with mid-day photosynthetic activity.
**Behavioral:** Genes peak in mid-day, supporting photosynthetic electron transport and photosystem II function during peak light intensity.
**Notes:** Strong enrichment for photosynthesis-related genes supports functional assignment.

### Cluster 4 | up | medium
**Name:** Prochlorococcus cluster 4 (up, photosynthesis)
**Enrichment:** Cytochrome b6/f (p=0.073, sig=False)
**Functional:** Cluster 4 contains 62 genes enriched for cytochrome b6/f complex. Genes peak at 15.8 h, in the afternoon, supporting photosynthetic electron transport during the latter part of the day.
**Behavioral:** Genes peak in the afternoon, likely supporting photosynthetic electron transport chain activity.
**Notes:** Enrichment is marginally significant; functional role inferred from gene annotations.

### Cluster 5 | up | high
**Name:** Prochlorococcus cluster 5 (up, respiration and protein degradation)
**Enrichment:** Respiratory terminal oxidases and protein degradation (p=0.019, sig=True)
**Functional:** Cluster 5 includes 120 genes enriched for respiratory terminal oxidases and degradation of proteins, peptides, and glycopeptides. Genes peak at 17.5 h, coinciding with the transition to night and increased respiration and protein turnover.
**Behavioral:** Genes peak in the early night, supporting respiratory metabolism and protein degradation during dark periods.
**Notes:** Significant enrichment and timing consistent with night-time respiration support confidence.

### Cluster 6 | up | high
**Name:** Prochlorococcus cluster 6 (up, purine biosynthesis)
**Enrichment:** Purine ribonucleotide biosynthesis (p=0.0049, sig=True)
**Functional:** Cluster 6 contains 138 genes enriched for purine ribonucleotide biosynthesis. Genes peak at 18.6 h, during the night, consistent with nucleotide biosynthesis supporting DNA replication and cell division.
**Behavioral:** Genes peak in the night, supporting nucleotide biosynthesis during DNA synthesis phase.
**Notes:** Strong enrichment and timing consistent with cell cycle phase support confidence.

### Cluster 7 | up | medium
**Name:** Prochlorococcus cluster 7 (up, nitrogen metabolism)
**Enrichment:** Nitrogen metabolism (p=0.087, sig=False)
**Functional:** Cluster 7 includes 121 genes enriched for nitrogen metabolism. Genes peak at 20.1 h, near sunset, consistent with nitrogen uptake and assimilation during the night.
**Behavioral:** Genes peak near sunset, supporting nitrogen metabolism and assimilation during the night phase.
**Notes:** Enrichment is marginal; functional role inferred from gene annotations and timing.

### Cluster 8 | up | high
**Name:** Prochlorococcus cluster 8 (up, chaperones)
**Enrichment:** Chaperones (p=0.00023, sig=True)
**Functional:** Cluster 8 contains 90 genes enriched for chaperones. Genes peak at 21.0 h, in the late evening, likely supporting protein folding and stress responses during the night.
**Behavioral:** Genes peak in the late evening, supporting chaperone activity and protein maintenance at night.
**Notes:** Strong enrichment and clear timing support functional assignment.

### Cluster 9 | up | high
**Name:** Prochlorococcus cluster 9 (up, RNA synthesis and transcription)
**Enrichment:** RNA synthesis, modification, and DNA transcription (p=0.035, sig=True)
**Functional:** Cluster 9 includes 99 genes enriched for RNA synthesis, modification, and DNA transcription. Genes peak at 22.4 h, late at night, consistent with preparation for transcriptional activity before dawn.
**Behavioral:** Genes peak late at night, supporting RNA synthesis and transcriptional regulation in preparation for the next day.
**Notes:** Significant enrichment and timing consistent with transcriptional preparation support confidence.

## Warnings

- [Tolonen 2006 / mit9313_kmeans_nstarvation / cluster 2] locus tag in functional_description: PMT2246
- [Zinser 2009 / med4_diel_clusters / cluster 11] duplicate id: med4_down_not_described (also used by cluster 10)
- [Zinser 2009 / med4_diel_clusters / cluster 12] duplicate id: med4_down_not_described (also used by cluster 11)
- [Zinser 2009 / med4_diel_clusters / cluster 14] duplicate id: med4_down_not_described (also used by cluster 12)
- [Zinser 2009 / med4_diel_clusters / cluster 15] duplicate id: med4_down_not_described (also used by cluster 14)
