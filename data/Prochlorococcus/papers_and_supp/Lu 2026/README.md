# Lu 2026

**Citation:** Lu Z, Plummer S, Kizziah J, Biller SJ, Morris JJ. "Vesicle-associated exudates from *Alteromonas* enhance growth and survival of *Prochlorococcus* in batch culture." *Applied and Environmental Microbiology*, 2026.
**DOI:** 10.1128/aem.00798-26
**Organism(s):** *Alteromonas macleodii* EZ55, as three LTPE-experiment strains (LTPE26, LTPE397, LTPE403); coculture context is *Prochlorococcus* MIT9312.
**Topic:** Proteomic characterization of *Alteromonas* exudates and extracellular membrane vesicles (EVs) that promote *Prochlorococcus* MIT9312 growth and survival across all phases of batch culture. Identifies enzymatic activities (alkaline phosphatase, peroxide degradation, siderophore capability) consistent with leaky Black Queen functions packaged into EVs. Table S1 reports per-protein detection in the >50 kDa supernatant size fraction and in crude cell lysates, separately for each of the three LTPE strains.

## Publication history — supersedes the bioRxiv preprint

This paper was previously integrated from its preprint at `Lu 2025/` (bioRxiv `10.1101/2025.05.28.656624`, titled "Enzymatically active exudates from *Alteromonas* facilitate *Prochlorococcus* survival in stationary phase"). The published AEM version replaces it entirely; the preprint directory was deleted on integration of this one.

**Node IDs changed** with the DOI: the `Publication`, both `Experiment` (`{doi}_{experiment_key}`), and all `DerivedMetric` (`derived_metric:{doi_short}:...`) node IDs differ from the preprint's. Nothing outside this directory referenced the old DOI except `cache/pdf_extraction_cache.json`.

**Data table changed shape.** The published `aem.00798-26-s0002.csv` is the preprint's `table s1.csv` minus two author-derived columns (`number of replicates where found in Supernatant/Lysate`). The 602 accessions merge 1:1 and no other column was added or removed.

## Classification

**Bucket B — new metrics — re-integrated 2026-08-12 as six per-strain boolean DerivedMetrics.**

Per-protein detection in the >50 kDa cell-free supernatant (exoproteome) and crude cell lysate (whole-cell) fractions, reported per LTPE strain.

## Integration

`paperconfig.yaml` defines two non-DE Experiments (compartment-split per the schema rule), each with one `derived_metrics_table` entry carrying three boolean metrics:

| Experiment | omics_type | compartment | DerivedMetrics (boolean) | Source columns |
|---|---|---|---|---|
| `ez55_exudate_proteome` | `EXOPROTEOMICS` | `exoproteome` | `exoproteome_detected_ancestor`, `exoproteome_detected_evolved_400ppm`, `exoproteome_detected_evolved_800ppm` | `LTPE26_Supernatant`, `LTPE397_Supernatant`, `LTPE403_Supernatant` |
| `ez55_whole_cell_proteome` | `PROTEOMICS` | `whole_cell` | `whole_cell_detected_ancestor`, `whole_cell_detected_evolved_400ppm`, `whole_cell_detected_evolved_800ppm` | `LTPE26_Lysate`, `LTPE397_Lysate`, `LTPE403_Lysate` |

Each emits `Derived_metric_flags_gene` edges. `true_tokens: ["yes"]` / `false_tokens: ["no"]` (the CSV's literal cell values); `blank_policy: skip`. All six `metric_type`s are registered in `multiomics_kg/vocab/non_de_evidence.py KNOWN_METRIC_TYPES`.

`treatment_type: [compartment]` on both Experiments — the paper's framing IS exudate-vs-lysate proteome comparison, even though no formal DE table is provided. `background_factors: [axenic]` — EZ55 was grown alone in Pro99 + 0.1% glucose, separate from the coculture growth experiments described elsewhere in the paper.

### Why per-strain booleans and not a 0–3 detection count

**The three LTPE strains are not replicates.** LTPE26 is the ancestral lineage; LTPE397 and LTPE403 evolved 500 generations at 400 ppm and 800 ppm pCO₂ respectively. The preprint integration (2026-04-29) modelled them as replicates and collapsed detection into two 0–3 numeric counts. That was wrong on the paper's own terms:

- The abstract's claim is that *"both the composition and activity of the secreted material changed after 500 generations of experimental evolution, suggesting genetic control."*
- The paper quantifies exactly this: the 800 ppm strain's exudate proteome shared only **1.7% and 4.2%** of proteins with the ancestor and the 400 ppm strain respectively, versus a much larger overlap between the ancestor and the 400 ppm strain.
- Three supplemental figures are devoted to the axis (S2 extracted proteins before/after evolution, S12 proteome before/after evolution, S13 intracellular/extracellular overlap).

Under a count, `2` is ambiguous between {ancestor, 400 ppm}, {400 ppm, 800 ppm} and {ancestor, 800 ppm} — i.e. "lost during evolution", "gained during evolution" and "lost only at 400 ppm" become indistinguishable. The count is recoverable from the booleans (`sum(detected)`); the booleans are not recoverable from the count. Counts are therefore not emitted.

**Cross-paper join.** `Lu 2025/` (Lu et al., *ISME J* 19(1) wrae259, `10.1093/ismejo/wrae259`) is the same group's report on the same LTPE experiment, listing convergently mutated EZ55 genes binned by pCO₂ treatment. Keeping per-strain resolution here is what makes a genotype (mutations) → phenotype (exudate proteome) join possible on the same lineages.

## Available data inventory

| File | Type | Content | KG status |
|------|------|---------|-----------|
| `lu-et-al-2026-...-in-batch.pdf` | PDF | Main AEM article | reference |
| `aem.00798-26-s0001.pdf` | PDF | Supplemental figures S1–S20 only; no additional data tables | reference |
| `aem.00798-26-s0002.csv` | CSV | Table S1 — 602 proteins detected by LC-MS/MS. Identifier column `Accession` = RefSeq WP_ (Tier-2 protein_id resolution against EZ55). Used columns: `Accession`, `Name` (product), the six `LTPE*_{Supernatant,Lysate}` detection columns. | integrated |

## What was dropped at integration

- **Ten KEGG-category boolean columns** (Amino Acid Metabolism, Carbohydrate Metabolism, …): functional-category presence flags derived from KO assignments. Duplicative with the KG's existing `KeggTerm` / `BriteCategory` ontology coverage; per-protein KEGG facts arrive via the standard annotation pipeline, not via paper-supplied flags.
- **`KO`, `MW`, `Localization`, `Score` columns**: KO and localization are covered by the annotation pipeline (eggNOG, PSORTb); MW and the ProteinPilot confidence score are per-run artifacts.

## Notes

- MIT9312 is the *Prochlorococcus* host in the paper's growth-phenotype experiments. The proteomics samples are NOT from coculture — EZ55 was grown axenically in Pro99 + 0.1% glucose at 30 °C / 120 rpm shaking, scaled to 2 L, then size-fractionated by tangential flow filtration (0.22 µm → 50 kDa).
- `light_condition` / `light_intensity` are absent (validator warns): the paper does not report a light regime for these shaken heterotroph cultures.
- Both Experiments trigger the validator's "no analyses referencing it" warning — expected for a DerivedMetric-only paper with no `statistical_analyses`.

## Known limitation: ~39% gene ID resolution against EZ55

**Of 602 protein accessions in Table S1, only 234 (38.9%) resolve against our deployed EZ55 (`GCF_901457815.2`).** Unchanged from the preprint — the accession set is identical. The paper's "NCBI *Alteromonas* EZ55 protein database" (per Methods) is in fact a multi-strain *A. macleodii* protein set, not an EZ55-specific reference.

Investigation 2026-04-30 (cross-checked the 368 unresolved IDs against `cache/data/Alteromonas/uniprot/28108/uniprot_preprocess_data.json`, then UniProt's proteome metadata):

| Bucket | Count | What it is |
|---|---:|---|
| Sister-strain proteins (taxid 28108, not EZ55) | 245 | 214 in `UP000095392` (strain **KCP01**, GCA_001750365.1) + 26 in `UP000063991` (strain **D7**, GCA_001562235.1) + 5 in `UP000509458`. EZ55 has no UniProt proteome of its own; UniProt's species-level taxid lumps 20 different *A. macleodii* strain proteomes together |
| `WP_*` not in any cached UniProt proteome | 113 | Likely from an assembly we don't cache (PGAP-retired or older) |
| `OES*` GenBank-only accessions (BioProject `PRJNA338983`) | 39 | No RefSeq counterpart; no match in our cache |

These IDs are not recoverable as EZ55 locus tags via `build_gene_id_mapping.py` — they aren't EZ55 proteins. ~67% (245/368) belong to identifiable sister strains; the rest aren't in our cache at all.

**Earlier experiment 2026-04-29** (now superseded): tried adding an `annotation_gff` for `GCF_000808635.1` (AD037), which yielded zero improvement. The two assemblies use disjoint locus-tag schemes, so locus-tag-based anchoring can't bridge them. Protein-sequence bridging (`generate: method: diamond_protein_match`) would only help the small subset that IS actually EZ55 sequence — most unresolved IDs aren't.
