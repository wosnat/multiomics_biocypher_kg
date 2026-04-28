# McDonagh 2012

**Citation:** McDonagh B, Domínguez-Martín MA, Gómez-Baena G, López-Lozano A, Diez J, Bárcena JA, García-Fernández JM. Nitrogen starvation induces extensive changes in the redox proteome of *Prochlorococcus* sp. strain SS120. Environmental Microbiology Reports, 2012, 4(2): 257–267.
**DOI:** 10.1111/j.1758-2229.2012.00329.x
**Organism(s):** *Prochlorococcus* SS120
**Topic:** Shotgun redox proteomics (biotin-HPDP labeling of reversibly oxidized cysteines) in SS120 under 24-h nitrogen starvation vs. replete control. Identifies ~80 proteins with Cys modifications; Table 1 = 7 proteins uniquely detected in one condition only; Table 2 = 20 proteins with significant ANOVA fold-change between conditions. Supp Table S1 holds the full identification list across all three replicates of both conditions.

## Classification

**Bucket D — defer / nothing actionable**

S1 supplement contains per-replicate Sequest scores (low-value DM). Tables 1 and 2 in the paper carry real fold-change values but are PDF-only — extraction effort exceeds expected KG value relative to other backlog items. The DerivedMetric loader (shipped 2026-04-21, see `docs/kg-changes/non-de-evidence-extension.md`) could in principle ingest extracted FC values, but until somebody manually extracts those PDF tables there is no machine-readable input.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `Environ Microbiol Rep - 2012 - McDonagh - Nitrogen starvation induces extensive changes in the redox proteome of.pdf` | PDF | Main article | reference | — |
| `legend.txt` | TXT | One-line legend describing Table S1 | reference | — |
| `emi4_329_sm_ts1.xls` | XLS | Supp Table S1: all proteins identified by MS/MS in at least 2 of 3 independent cultures per strain per condition; columns per replicate (Sequest score, % coverage, unique peptides); grouped by control (3 reps) and N-starved (3 reps) | add | Export Table S1 to CSV. Status depends on content: if it includes cross-condition fold-change or p-values per protein → add as `csv` with `statistical_analyses`; if it is purely per-replicate identification counts → add as `derived_metrics_table` with boolean `detected_in_control` / `detected_in_Nstarved` + numeric `seqeust_score`, `percent_coverage`, `unique_peptides` per condition |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Add** — Convert `emi4_329_sm_ts1.xls` to CSV and inspect column structure. Two plausible integrations:
   - **If** the XLS contains a computed log2 fold-change + p-value column per protein → add as a `csv` entry with `statistical_analyses` and `omics_type: PROTEOMICS`, treatment="nitrogen starvation", control="replete", `test_type: "ANOVA (Progenesis LC-MS)"`. Approx. 10 proteins with ANOVA-significant Cys-peptide fold-change (Table 2 of the paper) would be integrated as DE edges.
   - **Otherwise** (per-replicate raw scores only): add as a `derived_metrics_table` with numeric metrics (`sequest_score_control`, `sequest_score_nstarved`, `coverage_control`, `coverage_nstarved`, `unique_peptides_control`, `unique_peptides_nstarved`) + boolean flags (`detected_in_control`, `detected_in_nstarved`, `cys_modified`).
2. **Add (manual curation)** — The 20 proteins in paper's Table 2 (reversibly oxidized cysteines with significant fold-change between conditions) and the 7 proteins in Table 1 (detected only under one condition) are small, high-value sets. Consider hand-curating these ~27 proteins with fold-change + p-value (Table 2) or boolean detection flags (Table 1) as a small DE-style CSV even if the large S1 table turns out to be too raw.
3. **Skip** — Do not attempt integration until the XLS has been examined; the legend hints strongly that S1 is per-replicate raw scores, not aggregated DE.
4. **Add organism check** — SS120 is already deployed in the KG (added 2026-04-15, see `docs/community_proteomics_marref_saga.md`). No organism action needed.

## Notes

- Gene IDs in the paper's Table 1 use legacy NCBI GI numbers (e.g., `157413499`). These map to RefSeq WP_ / protein_id via `gene_id_mapping.json` Tier-2 `protein_id` ID type — match rate should be high via `multi_lookup`, but verify with `/check-gene-ids`.
- SS120 uses the Pro_#### locus-tag prefix in its NCBI annotation. The XLS may list protein names or GI numbers; an `id_columns` entry with `id_type: protein_id` (+ fall-back `gene_name`) will likely be needed.
- Protein Cys oxidation (reversible glutathionylation) is a post-translational signal and doesn't fit neatly into either DE or DerivedMetric schemas — but the fold-change summary in paper Table 2 is still a legitimate protein-level abundance comparison.
