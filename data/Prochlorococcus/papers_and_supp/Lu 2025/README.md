# Lu 2025

**Citation:** Lu Z, Plummer S, Kizziah J, Biller SJ, Morris JJ. "Enzymatically active exudates from *Alteromonas* facilitate *Prochlorococcus* survival in stationary phase." *bioRxiv* preprint, posted 2025-05-28.
**DOI:** 10.1101/2025.05.28.656624
**Organism(s):** *Alteromonas macleodii* EZ55 (and three LTPE-experiment evolved derivatives: LTPE26, LTPE397, LTPE403); coculture context is *Prochlorococcus* MIT9312.
**Topic:** Proteomic characterization of *Alteromonas* exudates and extracellular membrane vesicles (EVs) that promote *Prochlorococcus* MIT9312 survival in stationary phase. Identifies enzymatic activities (alkaline phosphatase, peroxide degradation, etc.) consistent with leaky Black Queen functions packaged into EVs. Supplementary table 1 reports per-protein detection in the >50 kDa supernatant size fraction and in crude cell lysates across three EZ55-evolved derivatives.

## Classification

**Bucket B — new metrics — integrated 2026-04-29 via EZ55 roll-up.**

Per-protein detection in the >50 kDa cell-free supernatant (exoproteome) and crude cell lysate (whole-cell) fractions of *Alteromonas macleodii* EZ55. Three EZ55-evolved sublines from the LTPE experiment (LTPE26 = ancestor, LTPE397 = 400 ppm pCO2 evolved 500 generations, LTPE403 = 800 ppm pCO2 evolved 500 generations) are treated as biological replicates of the EZ55 parent — the per-subline distinction is collapsed at integration. The first paper in this project to use the `EXOPROTEOMICS` / `exoproteome` compartment vocabulary.

## Integration

`paperconfig.yaml` defines two non-DE Experiments (compartment-split per the schema rule), each with one numeric `derived_metrics_table` entry:

| Experiment | omics_type | compartment | DerivedMetric (numeric, 0–3) | Source column |
|---|---|---|---|---|
| `ez55_exudate_proteome` | `EXOPROTEOMICS` | `exoproteome` | `exoproteome_detection_replicates` | "number of replicates where found in Supernatant ..." |
| `ez55_whole_cell_proteome` | `PROTEOMICS` | `whole_cell` | `whole_cell_detection_replicates` | "number of replicates where found in Lysate ..." |

Each metric is a 0–3 integer = "number of LTPE-evolved sublines (out of 3) in which this protein was detected by LC-MS/MS in this fraction". `rankable: "false"` (mass ties on a 4-value ordinal make percentile/bucket meaningless). `has_p_value: "false"`. Both metric_types are novel — see `multiomics_kg/vocab/non_de_evidence.py KNOWN_METRIC_TYPES` (validator warning, expected for first user).

`treatment_type: [compartment]` on both Experiments — the paper's biological framing IS exudate-vs-lysate proteome comparison, even though no formal DE table is provided. `background_factors: [axenic]` — EZ55 was grown alone in Pro99 + 0.1% glucose, separate from coculture/CM growth experiments described elsewhere in the paper.

## Available data inventory

| File | Type | Content | KG status |
|------|------|---------|-----------|
| `2025.05.28.656624v1.full.pdf` | PDF | Main bioRxiv preprint | reference |
| `table s1.csv` | CSV | 602 proteins detected by LC-MS/MS. Identifier column `Accession` = RefSeq WP_ (Tier-2 protein_id resolution against EZ55). Used columns: `Accession`, `Name` (product), the two replicate-count columns (`number of replicates where found in Supernatant ...` and `number of replicates where found in Lysate ...`). | integrated |

## What was dropped at integration

- **Six per-LTPE boolean columns** (`LTPE26_Supernatant`, `LTPE397_Supernatant`, `LTPE403_Supernatant`, `LTPE26_Lysate`, `LTPE397_Lysate`, `LTPE403_Lysate`): under the EZ55 roll-up these are subsumed by the replicate-count columns (the count is exactly the sum of the LTPE booleans for that fraction). Keeping them would emit redundant DM edges.
- **Ten KEGG-category boolean columns** (Amino Acid Metabolism, Carbohydrate Metabolism, ...): functional-category presence flags derived from KO assignments. Duplicative with the KG's existing `KeggTerm` / `BriteCategory` ontology coverage; per-protein KEGG facts arrive via the standard annotation pipeline, not via paper-supplied flags.

## Notes

- MIT9312 is the *Prochlorococcus* host in the paper's growth-phenotype experiments. The proteomics samples are NOT from coculture — EZ55 was grown axenically in Pro99 + 0.1% glucose at 30°C / 120 rpm shaking, scaled to 2 L, then size-fractionated by tangential flow filtration (0.22 µm → 50 kDa).
- Preprint, not peer-reviewed; treat citation as bioRxiv-only.
- The per-subline distinction (LTPE26 vs LTPE397 vs LTPE403, capturing pCO2-evolution effects) is intentionally lost under the roll-up. If a future analysis needs it, deploying the three sublines as separate OrganismTaxon nodes is reversible — would require re-resolving the paperconfig with three separate per-subline boolean DerivedMetrics instead of the two count-based numeric ones.

## Known limitation: ~39% gene ID resolution against EZ55

The paper's "NCBI *Alteromonas* EZ55 protein database" (per Methods) does NOT correspond to a single deployable assembly. **Of 602 protein accessions in Table S1, only 234 (38.9%) resolve against our deployed EZ55 (`GCF_901457815.2`).** The 368 unresolved break down as:

| Source | Count | Disposition |
|---|---:|---|
| `WP_*` accessions in `GCF_000808635.1` (ASM80863v1, an alternate EZ55 RefSeq assembly) but NOT in `GCF_901457815.2` | ~88 | Could be recovered via `id_translation` + `generate: method: diamond_protein_match` against ASM80863v1 protein FASTA |
| `WP_*` accessions in NEITHER deployed nor ASM80863v1 (PGAP-retired or from yet another assembly) | ~180 | Likely unrecoverable without identifying the exact protein-DB version the authors used |
| `OES*` accessions (41) — INSDC GenBank prefix from BioProject `PRJNA338983` (Cusick et al. 2016, US Naval Research Lab — *"Draft Genome Sequence of four Alteromonas macleodii strains isolated from copper coupons"*) | 41 | A separate GenBank-only submission with no RefSeq counterpart. None match `GCF_901457815.2` or `GCF_000808635.1`. Recovery would require sequence-level bridging against the `MIPX*` WGS contigs |

**Experiment tried 2026-04-29**: added `annotation_gff` entry pointing at `GCF_000808635.1`'s GFF, hoping `build_gene_id_mapping.py` would absorb its WP_ accessions into EZ55's mapping. Result: zero improvement. The two assemblies use disjoint locus-tag schemes (`RJ44_*` in ASM80863v1 vs `EZ55_*` / `ALTBGP6_*` in `GCF_901457815.2`), so the build's locus-tag-based anchoring couldn't bridge them — all ~5,000 GFF rows unanchorable. Anchoring would require **protein-sequence-level bridging** (`generate: method: diamond_protein_match`), the same machinery used for AltMedDE in the community-proteomics saga.

The 38.9% rate is the cost of integrating now vs. waiting for an EZ55 reference-proteome-match deployment (analogous to Marinobacter / Alt_MarRef per `docs/community_proteomics_marref_saga.md` and `docs/kg-changes/reference-proteome-match-organisms.md`). Decision was to ship at 38.9% since the integration is otherwise low-risk (proteomics-only, no DE, no time-course) and the paper is a bioRxiv preprint.
