# Cluster Extraction Review: Tolonen 2006
**Organism:** Prochlorococcus MED4  
**Extracted:** 2026-04-01T10:32:54

## Cluster 1: med4_up_transport_amino_acid_biosynthesis — PASS

- **functional:** Enrichment: Transport and binding proteins (table very_high, p=0.0133533; visual high, p=0.01; both significant). Amino acid biosynthesis is indicated in semantic data (p=0.127730, not significant). Key genes include PMM0970, PMM1462, PMM0958, PMM0920, PMM0370; visual data specifically highlights N-transport genes such as urtA and cynA. Overall, cluster 1 comprises transport/binding gene members with a subset linked to amino acid biosynthesis signals, skewed toward rapid nitrogen-stress responsiveness.
- **behavioral:** Upregulated rapidly and highly in response to nitrogen starvation; early timepoints (rapid response 0–3 hours) with strong induction of the listed genes (e.g., PMM0970, PMM1462, PMM0370).

## Cluster 2: med4_up_hli_transport — WARN

- **functional:** Enrichment: Transport and binding proteins (table very_high, p=0.124545; visual high, p=0.12; not significant). Notably contains high-light-inducible (hli) genes; visual data notes rapid, high induction of hli genes. Additional context from semantic data includes various transport-related and auxiliary functions (e.g., PMT0601 glutamate family-related gene, PMT1314). Representative genes in the cluster include PMM0690, PMM0337, PMM1390, PMM0374, PMM1041, PMM0336, PMM0030, PMM0365, PMM1391, PMM0689, PMM1463, PMM0371, PMM1623, PMM0687, PMM0246, PMM0359, plus hli-related entries PMT1154, PMT2246, PMT0993.
- **behavioral:** Upregulated in response to nitrogen starvation; rapid response with high-level induction of hli genes noted in visual data.
- **issue:** While the genes are correctly listed and the behavioral description is accurate, the enrichment for transport and binding proteins is not clearly supported as significant in the paper.

## Cluster 3: med4_up_hli_regulation_late — WARN

- **functional:** Enrichment: Regulatory functions (table very_high, p=0.0741993; visual high, p=0.07; not significant). Contains high-light-inducible (hli) genes. Genes include PMM0936, PMM1118, PMM1404, PMM0817, PMM1385, PMM0997, PMM1562, PMM0684, PMM1135, PMM0548, PMM1672, PMM0818, PMM1262, PMM0819, PMM0338, PMM1396, PMM1397, PMM0988, PMM1384. The semantic note explicitly points to a later, lesser hli response compared with cluster 2.
- **behavioral:** Upregulated, but later and to a lesser degree than cluster 2; represents a slower regulatory response to nitrogen starvation.
- **issue:** The regulatory function enrichment is not clearly significant in the paper, but the genes and behavioral description are correct.

## Cluster 4: med4_up_nstarvation_mixed_cluster4 — WARN

- **functional:** Enrichment: Amino acid biosynthesis (table very_high, p=0.180753; visual high, p=0.18; both not significant). The cluster contains a broad mix of functions: hypothetical proteins, transport, photosynthesis/respiration (carboxysome), translation, and metabolism. Semantic data lists additional diverse categories (PMT0411, PMT0493, PMT2123, PMT1355, PMT1454, PMT1566, PMT1203, PMT1275, PMT1781, PMT2130, PMT1333, PMT0453). Representative genes highlighted in supporting quotes include PMT0411, PMT0493, PMT2123, PMT1355, PMT1454, PMT1566, PMT1203, PMT1275, PMT1781, PMT2130, PMT1333, PMT0453.
- **behavioral:** Upregulated in response to nitrogen starvation over the 0–48 h time course, but the cluster comprises diverse functions with no single overrepresented pathway.
- **issue:** The enrichment for amino acid biosynthesis is not supported as significant. The genes listed are diverse, matching the paper's description of mixed functions.

## Cluster 5: med4_up_nstarvation_mixed_cluster5 — WARN

- **functional:** Enrichment: Regulatory functions (table very_high, p=0.178981; visual high, p=0.18; semantic medium, p=0.985297). Not significant overall. Notable entries include sigma factors; various genes span photosynthesis/respiration (e.g., PsbD D2), DNA replication/repair (DNA polymerase III alpha subunit), amino acid biosynthesis, and regulatory/hypothetical proteins. Supporting notes list PMM0911, PMM0012, PMM1157 (PsbD), PMM0380, PMM0719, PMM0945, PMM0169, PMM0733, PMM1272, PMM0908, among others.
- **behavioral:** Upregulated across the nitrogen starvation time course (0, 3, 6, 12, 24, 48 h) with a mixed functional signature including regulatory and photosynthesis-related genes; some sigma factors noted.
- **issue:** The regulatory function enrichment is not significant. The genes listed include a mix of functions, consistent with the paper.

## Cluster 6: med4_down_translation_cluster6 — WARN

- **functional:** Enrichment: Translation (table very_high, p=0.00141159; visual high, p=0.001). Strongly significant downregulation of translation-related genes. Representative components include numerous ribosomal proteins and translation-associated factors distributed across the long gene list (e.g., PMM1142, PMM1182, PMM1664, PMM0584, PMM1264, PMM1482, PMM0747, PMM0215, PMM1335, PMM1345, PMM0536, PMM0251, PMM0273, PMM0283, PMM1070, PMM0642, PMM0060, PMM1524, PMM1336, PMM1150, PMM0346, PMM0553, PMM0265, PMM0027, PMM0428, PMM0056, PMM1483, PMM0753, PMM0134, PMM1661, PMM1509, PMM0500, PMM0288, PMM1511, PMM0203, PMM0298, PMM1489, PMM0829, PMM0046, PMM1183, PMM0844, PMM1229, PMM1251, PMM0510, PMM0494, PMM0295, PMM0619, PMM0605, PMM0206, PMM1532, PMM1109, PMM0483, PMM0609, PMM1546, PMM0213, PMM0907, PMM1644, PMM0580, PMM1440, PMM1055, PMM1191, PMM0429, PMM1510, PMM0784, PMM0200, PMM0769, PMM0026, PMM1643, PMM0253, PMM1437, PMM0902, PMM1662, PMM1543, PMM0204, PMM1607, PMM0299, PMM0869, PMM0930, PMM0296, PMM1547, PMM1522, PMM0622, PMM1506, PMM1600, PMM1542, PMM0137, PMM1079, PMM0751, PMM1485, PMM0511, PMM1663, PMM1617, PMM0987, PMM0929, PMM1011, PMM0017, PMM0599, PMM1608, PMM1247, PMM0754, PMM1534, PMM0741, PMM0870, PMM1530, PMM1618, PMM0474, PMM0560, PMM0700, PMM0312, PMM0142, PMM0347, PMM0348, PMM0016, PMM1375, PMM0135, PMM0864, PMM1060, PMM1544, PMM1585, PMM1548, PMM0924, PMM0356, PMM0876, PMM0297
- **behavioral:** Direction: down; Temporal pattern: not described in paper; Treatment: nitrogen starvation.
- **issue:** The enrichment for translation is strongly supported, and the genes are correct. However, the temporal pattern is not described in the paper.

## Cluster 7: med4_down_translation_cluster7 — PASS

- **functional:** Enrichment: Translation (table very_high, p=2.70545e-09; visual high, p=2.7e-9; semantic medium, p=7.8986e-10). Strong, coordinated downregulation of translation-related genes, including ribosomal proteins and translation-associated products. Representative genes cited in supporting notes include PMM1553, PMM1610, PMM0202, PMM1537, PMM1655, PMM1557, PMM1706, PMM0583, PMM1550, PMM1545, PMM0781, PMM1298, PMM1536, PMM1556, PMM0552, PMM1549, PMM0528, PMM1028, PMM1629, PMM1558, PMM1540, PMM0673, PMM1609, PMM1456, PMM1555, PMM1297, PMM1539, PMM1457, PMM1551, PMM1552, PMM1554, PMM1402, PMM0201, PMM1538, PMM0041, PMM0541, PMM0023, PMM1541, PMM0410, PMM1507, PMM1285, PMM1439.
- **behavioral:** Downregulated under nitrogen starvation; timing across 0–48 h; pattern consistent with a broad shutdown of translation machinery.

## Cluster 8: med4_down_photosynthesis_respiration_cluster8 — PASS

- **functional:** Enrichment: Photosynthesis and respiration (table very_high, p=2.61704e-09; visual high, p=2.6e-9; semantic medium, p=0.985297). Contains many photosystem I components; MED4 cluster 8 enriched for Photosynthesis and Respiration, with psa genes as prominent members (psaBDEIJKLM). Representative genes include PSA subunits psaB, psaD, psaE, psaI, psaJ, psaK, psaL, psaM, as well as other photosynthesis/respiration-related genes.
- **behavioral:** Downregulated in response to nitrogen starvation, across the 0–48 h time course, with repression of photosynthetic/respiration genes.

## Cluster 9: med4_conflicting_photosynthesis_cluster9 — FAIL

- **functional:** Enrichment: Photosynthesis and respiration (table very_high, p=1.46201e-06; visual high, p=1.4e-8; semantic medium, p=1.46e-06). Contains ATP synthase subunits and carbon fixation genes (rbcLS). Key genes include PMM1451, PMM1455, PMM1453, PMM1452, PMM0550, PMM1450, PMM1454, PMM0551, plus rbcL and rbcS.
- **behavioral:** Conflicting directional signals for nitrogen starvation response: visual data indicate downregulation, while semantic data indicate upregulation for some Photosynthesis/Respiration genes (e.g., rbcLS and ATP synthase components) across the 0–48 h time course.
- **issue:** The enrichment is correct, but there are conflicting directional signals in the paper. The behavioral description does not match the paper's findings.
