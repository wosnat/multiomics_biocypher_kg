# Lu 2025

**Citation:** Lu Z, Plummer S, Kizziah J, Biller SJ, Morris JJ. "Enzymatically active exudates from *Alteromonas* facilitate *Prochlorococcus* survival in stationary phase." *bioRxiv* preprint, posted 2025-05-28.
**DOI:** 10.1101/2025.05.28.656624
**Organism(s):** *Alteromonas macleodii* EZ55 (and three LTPE-experiment evolved derivatives: LTPE26, LTPE397, LTPE403); coculture context is *Prochlorococcus* MIT9312.
**Topic:** Proteomic characterization of *Alteromonas* exudates and extracellular membrane vesicles (EVs) that promote *Prochlorococcus* MIT9312 survival in stationary phase. Identifies enzymatic activities (alkaline phosphatase, peroxide degradation, etc.) consistent with leaky Black Queen functions packaged into EVs. Supplementary table 1 reports per-protein detection in the >50 kDa supernatant size fraction and in crude cell lysates across three EZ55-evolved derivatives.

## Classification

**Bucket B — new metrics / DE / resolution (want to add) — strain question to resolve first**

Per-protein detection (yes/no in supernatant + lysate fractions) across three EZ55-evolved derivatives (LTPE26, LTPE397, LTPE403). Fits the boolean `derived_metrics_table` pattern shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`). **Open decision before integrating:** treat LTPE26/397/403 as new strains (would require deploying three near-identical genomes) or roll the data up to the parent EZ55 strain (already deployed) and lose subline distinction. The simpler path is roll-up; the higher-fidelity path is full strain deployment. No paperconfig until that call is made.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `2025.05.28.656624v1.full.pdf` | PDF | Main bioRxiv preprint | reference | — |
| `table s1.csv` | CSV | Per-protein boolean detection (yes/no) per LTPE derivative x {Supernatant, Lysate}, replicate-counts (0-3) per fraction, KEGG metabolic-category boolean flags. ~hundreds of proteins. Identifier column is `Accession` = RefSeq WP_ protein IDs (e.g. `WP_014950802`). | not in | Pending strain decision (LTPE roll-up vs deploy). Once decided, wire as `derived_metrics_table` with `value_kind: boolean` (detection flags) + `value_kind: numeric` (replicate counts) per `docs/kg-changes/non-de-evidence-extension.md`; Tier-2 `protein_id` resolution against EZ55 (parent). |

## Notes

- Identifier column is `Accession` = RefSeq WP_ protein IDs. Tier-2 protein_id lookup against EZ55 will resolve cleanly once integration is unblocked, but only at the EZ55-parent level — the per-derivative distinction (LTPE26 vs LTPE397 vs LTPE403) requires the three sublines to exist as separate OrganismTaxon nodes.
- MIT9312 is the *Prochlorococcus* host in this study; it is already deployed in the KG.
- The KEGG category columns (Amino Acid Metabolism, Carbohydrate Metabolism, etc.) are functional-category presence flags derived from KO assignments — duplicative with the KG's existing `KeggTerm` ontology coverage; should NOT be re-imported as gene properties when this paper is finally integrated.
- Preprint, not peer-reviewed; treat citation as bioRxiv-only.
