# Metabolite scaffold — design overview

**Date:** 2026-04-28
**Status:** Phase 1 split into 3 sub-specs (each independently shippable); Phases 2 & 3 sketched
**Driver goal:** identify nutrient exchanges between *Prochlorococcus* and *Alteromonas*

## Goal

Add metabolites and the surrounding metabolism scaffold so we can answer, in increasing depth:

1. Which compounds are predicted as Pro→Alt or Alt→Pro nutrient exchanges from genome content alone? (Phase 1)
2. Of those, which have actually been measured in published experiments — especially in cocultures and exometabolome samples? (Phase 2)
3. What does a query layer over both look like? (Phase 3)

## Architecture overview

Three independent phases, each independently shippable + testable. Phase 1 is decomposed into 3 sub-specs.

| Phase / spec | Scope | Output |
|---|---|---|
| **1.1 Foundation** ([spec](2026-04-28-metabolite-foundation-design.md)) | MNX/TCDB/CAZy reference data + SQLite resolver + step 4 eggNOG enrichment | enriched `gene_annotations_merged.json`; **no KG-visible changes yet** |
| **1.2 Reactions scaffold** ([spec](2026-04-28-metabolite-reactions-scaffold-design.md)) | Metabolite + Reaction nodes; intra-scaffold edges; `Gene_catalyzes_reaction`; `Organism_has_metabolite` (substrate/product columns) | Pro↔Alt driver query works |
| **1.3 Transport + CAZy** ([spec](2026-04-28-metabolite-transport-cazy-design.md)) | TransporterClass + CazyFamily nodes + hierarchies; gene edges; backfills `Organism_has_metabolite` TC/CAZy fields | enriched substrate-class lookups |
| **2 (sketch below)** | Paper measurement evidence (MetaboliteAssay) | spec deferred until Phase 1 lands |
| **3 (sketch below)** | MCP tooling + augmentation (KofamScan/gapseq/...) | spec deferred until Phase 2 lands |

Identifier backbone (used by all three Phase-1 sub-specs): **MetaNetX MNXref** (release 4.5+) — Rosetta-stone for chemistry IDs, cross-referencing KEGG, ChEBI, BiGG, MetaCyc, ModelSEED, HMDB, Reactome, Rhea, LipidMaps, SwissLipids, VMH, enviPath. IDs: `mnx:MNXM<NNNNNN>` (compounds), `mnx:MNXR<NNNNNN>` (reactions).

## Phase 2 — Paper measurement evidence (sketch)

Detailed design deferred until Phase 1 lands. Key shape decisions made:

- **`MetaboliteAssay`** as its own node type (separate from DerivedMetric — different terminology, different domain), but reusing every property name that maps cleanly: `name`, `value_kind`, `unit`, `field_description`, `compartment`, `metric_type`, `rankable`, `has_p_value`, `p_value_threshold`, denormalized block (`experiment_id`, `organism_name`, `publication_doi`, `omics_type`, `treatment_type`, `background_factors`, `treatment`, `light_condition`, `experimental_context`).
- **Three measurement edge types** mirroring the DerivedMetric pattern: `Assay_quantifies_metabolite` (numeric), `Assay_flags_metabolite` (boolean / presence), `Assay_classifies_metabolite` (categorical).
- **Binding edges**: `Experiment_has_metabolite_assay`, `Publication_has_metabolite_assay`, `MetaboliteAssayBelongsToOrganism`.
- **Edge properties on `Assay_quantifies_metabolite`** parallel to `Changes_expression_of`: `value`, `value_se`, `condition_label`, `time_point`, `time_point_order`, `time_point_hours`, `log2_fold_change`, `adjusted_p_value`. Post-import `rank_by_metric`, `metric_percentile`, `metric_bucket`, `significant` (parallel to DM numeric edges).
- **Augments existing `Organism_has_metabolite` edges** with measurement properties (`measured_assay_count`, `measured_compartments`, `measured_value_kinds`, `measured_paper_count`).
- **Resolution** uses the Phase-1 `metabolite_resolver.db`. Unresolved compounds logged in a per-paper `metabolite_resolution_report.json`. Policy for unresolved compounds (drop / keep as `unresolved_metabolite_feature` node / manual override) is decided per paper.
- **Compartment vocab** decided here once we see actual paper data; likely additions: `exometabolome` and possibly `particulate` / `lipid_fraction`.

Marquee paper candidates for the first integration:

- **Capovilla 2023** — chitin utilization, MIT9303/MIT9313 + GlcNAc, paired DE genes + intracellular metabolites. CSV-ready.
- **Kujawinski 2023** — Pro ecotype metabolite diversity. Profiling-style data, exercises `Assay_flags_metabolite`. PDF + docx supplementary; needs extraction.
- **Szul 2019** — Pro N-limitation metabolomics. Multi-condition `condition_label` + isotopologue time-course (deferred or treated specially).
- **Fang 2019** — viral lysis DOM products. Exometabolome compartment. Needs PDF extraction.

## Phase 3 — Tooling and augmentation (sketch)

- **MCP tools**: `metabolite_overview`, `predicted_exchange_candidates`, `measured_exchange_candidates`, `metabolite_response_profile`, `genes_predicted_to_produce_metabolite`, `genes_predicted_to_consume_metabolite`, `metabolites_in_compartment`.
- **OrthologGroup inheritance for reactions**: post-import `Reaction_catalyzed_by_ortholog_group` edges with `support` property (fraction of OG members with the catalyzing KO). Lifts effective gene→reaction reachability from ~52% to ~85–95% via OG bridging.
- **KofamScan**: parallel KO source via KEGG-official HMM scoring. Adds ~5–15% KO coverage. New step in `prepare_data.sh`.
- **gapseq**: full genome-scale model reconstruction per strain. Output is per-strain SBML with explicit gene-reaction associations and predicted exchange reactions. Heaviest, but directly answers "predicted Pro–Alt nutrient exchanges" with a curated network. Output parsed into additional `Gene_catalyzes_reaction` edges with `evidence_source ∈ {gapseq_curated, gapseq_gapfilled}`.
- **CarveMe** (alternative to gapseq): faster ML-based GEM reconstruction, BiGG namespace, smaller compound coverage.
- **Pfam → reaction bridges**: derived from UniProt-curated Pfam ↔ Rhea links. Adds another precision tier (~78% gene coverage). Less precise than KO; included only if Phase 2 reveals coverage gaps.
- **antiSMASH**: biosynthetic gene clusters + secondary-metabolite class predictions. Useful if marine NRPS/PKS compounds become a topic.
- **BiGG model GPR ingestion**: parse BiGG model JSON files to extract gene-reaction rules; resolves the eggNOG `BiGG_Reaction` column. Likely subsumed by gapseq.
- **TC → Metabolite direct edges**: for TC entries with explicit substrate compounds in TCDB or KEGG TRANSPORTER. Supplements the `substrate_classes` text property with compound-level edges.

## Open questions (revisit based on data)

These are explicitly left open and revisited after each phase lands:

1. **Direct vs KO-mediated `Gene → Reaction` edges.** Phase 1 keeps both: `Gene_catalyzes_reaction` (direct, from eggNOG `KEGG_Reaction` column, ~25% gene coverage, per-gene precision) AND `Reaction_catalyzed_by → KeggTerm` (MNX-derived, ~52% gene coverage via existing `Gene_has_kegg_ko`, class-level). Trade-off: redundancy for KEGG-sourced reactions vs. precision/coverage tier flexibility. Revisit after observing query-pattern usage.
2. **OrthologGroup inheritance for reactions.** Deferred from Phase 1; Phase 3 candidate.
3. **Compartment vocabulary for Phase 2.** Existing Experiment values: `whole_cell` (default), `vesicle`, `exoproteome`, `secretome`. Phase 2 likely adds `exometabolome`. Final vocab decided once paper data is examined.
4. **Phase 2 marquee paper choice.** TBD when Phase 2 starts; Capovilla 2023 leading candidate.
5. **Unresolved-name policy for Phase 2.** Per-paper decision: drop / keep as `unresolved_metabolite_feature` node / manual override. Default proposal: drop with logging.
6. **TransporterClass → Metabolite compound-level edges.** Phase 1 keeps TC substrate info as text property (`substrate_classes`). Phase 2/3 may add direct edges where TCDB/KEGG TRANSPORTER has explicit compound data.

## Out of scope for Phase 1

- **`BiGG_Reaction` ingestion.** eggNOG column has gene refs (`<model>.<gene>`) not reaction IDs; needs BiGG model JSON GPR parsing. Coverage marginal beyond `KEGG_Reaction` (~14 BiGG-only genes for MED4, ~88 for EZ55). Phase 3 candidate, likely subsumed by gapseq.
- **`KEGG_Pathway`, `KEGG_Module` nodes.** Functional groupings; partial overlap with BRITE.
- **`KEGG_rclass`.** Abstract reaction patterns.
- **Lipid-specific extensions.** Lipid metabolites reachable via MNX (LipidMaps + SwissLipids), but lipid-class hierarchies and LM-specific properties not modelled.
- **Thermodynamics, charge balancing, sub-integer stoichiometry.** Phase 1 stoichiometry is the integer coefficient from MNX.
- **Paper measurement integration of any kind.** Strictly Phase 2.
