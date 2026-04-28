# Christie-Oleza et al. 2017

**Title:** Nutrient recycling facilitates long-term stability of marine microbial phototroph-heterotroph interactions

**Journal:** Nature Microbiology 2, 17100 (2017)

**DOI:** 10.1038/nmicrobiol.2017.100

**Authors:** Joseph A. Christie-Oleza, Despoina Sousoni, Matthew Lloyd, Jean Armengaud, David J. Scanlan

## Classification

**Bucket D -- defer / nothing to do (already integrated)**

Synechococcus WH7803 + Ruegeria pomeroyi DSS-3 long-term coculture proteomics; both organisms are deployed in the KG. The three FC-bearing tables -- S2A (Syn early/ASW), S5A (Syn late/seawater), S7 (Ruegeria long-term/seawater) -- are wired to three experiments in paperconfig.yaml with `*_resolved.csv` siblings, using the documented signed-linear ratio convention plus p-value/q-value statistics. Tables S2B and S5B (R. pomeroyi abundance only) are intentionally skipped because they have no FC. Nothing further to add unless DerivedMetric ingestion is wired up, in which case the abundance-only side tables could move to bucket C.

## Summary

Long-term coculture study of the marine cyanobacterium *Synechococcus* sp. WH7803 with the heterotrophic bacterium *Ruegeria pomeroyi* DSS-3. Demonstrates that nutrient recycling between phototroph and heterotroph enables stable mutualistic interactions in both nutrient-rich ASW medium and nutrient-poor natural seawater. Comparative proteomics reveals functional specialization: *Synechococcus* dedicates resources to photosynthesis and CO2 fixation while relying on the heterotroph for remineralization of leaked organic matter.

## Organisms

- **Synechococcus sp. WH7803** (NCBI taxid: 32051) -- phototroph, NOT in KG yet (Phase 3 needed)
- **Ruegeria pomeroyi DSS-3** (NCBI taxid: 89184) -- heterotroph, NOT in KG yet (Phase 3 needed)

## Experiments

| Experiment | Organism | Medium | Time | Table |
|---|---|---|---|---|
| Early coculture (ASW) | Syn WH7803 | ASW (nutrient-rich) | 35 days | S2A |
| Late coculture (seawater) | Syn WH7803 | Natural seawater | 10 days | S5A |
| Long-term coculture (seawater) | R. pomeroyi DSS-3 | Natural seawater | Long-term | S7 |

## Methods

- **Proteomics:** Shotgun proteomics via nanoLC-MS/MS
  - ASW experiments: LTQ-Orbitrap XL, SDS-PAGE (6 bands), NSAF quantification
  - Seawater experiments: Orbitrap Fusion, SDS-PAGE (1 band), NAF quantification
- **Statistics:** t-test (S2A), ANOVA (S5A), Student's t-test via Perseus with q-value (S7)
- **Replicates:** 3 biological replicates per condition
- **Growth conditions:** 22C, 10 umol photons m-2 s-1 continuous light, 140 rpm shaking

## Supplementary Tables

| Table | Contents | Rows | FC? | Used? |
|---|---|---|---|---|
| S2A | Syn WH7803 proteins, ASW axenic vs coculture | 564 | Yes (signed linear) | Yes |
| S2B | R. pomeroyi proteins, ASW coculture | 37 | No (abundance only) | No |
| S5A | Syn WH7803 proteins, seawater axenic vs coculture | 1050 | Yes (signed linear) | Yes |
| S5B | R. pomeroyi proteins, seawater coculture | ~640 | No (abundance only) | No |
| S7 | R. pomeroyi proteins, seawater mono vs coculture | 1724 | Yes (signed linear + q-value) | Yes |

## Fold Change Convention

All tables use a **signed linear ratio** convention (not standard linear FC):
- **Positive values** = higher abundance in coculture (upregulated in treatment)
- **Negative values** = higher abundance in axenic/mono-culture (downregulated in treatment)
- No values between 0 and 1 (unlike standard linear FC where downregulation = 0 < FC < 1)

The sign convention is correct for the KG (positive = up in treatment), so values are passed through as-is without `fold_change_type` conversion. This matches the Ma 2022 precedent for non-standard signed ratios.

## Gene ID Format

- **Syn WH7803:** `SynWH7803_NNNN` (e.g., SynWH7803_0366) -- locus tag format, in "Other reference" column
- **R. pomeroyi DSS-3:** `SPONNNN` (e.g., SPO1796) -- annotation reference, in "Annotation reference" column

## Phase 3 Requirements

Both organisms need genome deployment before this paper can produce valid KG edges:
- Synechococcus WH7803: NCBI accession GCF_000063505.1, taxid 32051
- Ruegeria pomeroyi DSS-3: NCBI accession GCF_000011965.2, taxid 89184
