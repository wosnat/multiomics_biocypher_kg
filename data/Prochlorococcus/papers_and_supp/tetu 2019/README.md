# Tetu 2019

**Citation:** Tetu SG, Sarker I, Schrameyer V, Pickford R, Elbourne LDH, Moore LR, Paulsen IT. "Plastic leachates impair growth and oxygen production in *Prochlorococcus*, the ocean's most abundant photosynthetic bacteria." *Communications Biology* 2, 184 (2019).
**DOI:** 10.1038/s42003-019-0410-x
**Organism(s):** *Prochlorococcus* MIT9312 (HLII) and NATL2A (LLI)
**Topic:** Short-term exposure (90–120 min, AMP1 medium, 22 °C, continuous light) of axenic MIT9312 and NATL2A to leachate from HDPE (50%) or PVC (2%) consumer plastics. RNA-seq DESeq2 analysis; supplementary datasets S3/S4 list pre-filtered DE genes (significant for HDPE, PVC, or both).

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `s42003-019-0410-x.pdf` | PDF | Main paper | reference | — |
| `42003_2019_410_MOESM1_ESM.pdf` | PDF | Supplementary figures/tables document | reference | — |
| `42003_2019_410_MOESM2_ESM.pdf` | PDF | Supplementary data 2 (likely summary tables) | reference | — |
| `DE genes MIT9312 42003_2019_410_MOESM6_ESM.csv` | CSV | MIT9312 DE: original locus, 2017 NCBI locus, protein_id, product, HDPE log2FC, PVC log2FC | already in | — |
| `DE genes MIT9312 42003_2019_410_MOESM6_ESM.xlsx` | XLSX | Excel original of MOESM6 | skip | Skip (duplicate of csv already in). |
| `DE genes Natl2a 42003_2019_410_MOESM7_ESM.csv` | CSV | NATL2A DE: original locus, 2017 NCBI locus, protein_id, product, HDPE log2FC, PVC log2FC | already in | — |
| `DE genes Natl2a 42003_2019_410_MOESM7_ESM.xlsx` | XLSX | Excel original of MOESM7 | skip | Skip (duplicate of csv already in). |
| `DE genes MIT9312 42003_2019_410_MOESM6_ESM_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `DE genes Natl2a 42003_2019_410_MOESM7_ESM_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `*_resolved_report.txt` (×2) | text | Resolution reports | — | — (auto-generated; skip) |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Superseded draft | — | — (skip / historical) |

## Classification

**Bucket D — defer / nothing to do (fully integrated)**

Four DE experiments (HDPE/PVC × MIT9312/NATL2A) are all wired up; both strains deployed; tri-column ID resolution is robust. Supplementary PDFs are pre-filtered summary tables, no untapped per-gene measurements remain. Action #3 in this README ("verify MOESM1/MOESM2 don't carry hidden per-gene physiology") is a low-priority optional sanity check, not a blocker — the paper's growth/PSII assays are sample-level by design. No remaining integration work.

## Current paperconfig summary

- Experiments defined: 4 — HDPE-vs-control and PVC-vs-control for each of MIT9312 and NATL2A (all RNASEQ / DESeq2 / treatment_type=[plastic]).
- Statistical analyses (DE edges): 4 (one per experiment). Two analyses share each CSV (HDPE log2FC and PVC log2FC columns).
- Supplementary materials entry types: 2× `csv` (MOESM6 and MOESM7).
- Organisms covered: *Prochlorococcus* MIT9312, NATL2A.
- Table scope: `filtered_subset` — "genes significant for HDPE or PVC or both" (pre-filtered by the authors; `prefiltered: true` on each analysis).
- Non-DE evidence: none.
- ID resolution: multi-tier — `MIT9312 original locus tag` / `NATL2A original locus tag` → `old_locus_tag` (primary name_col); fallback columns `... locus tag (2017 NCBI)` → `locus_tag_ncbi` and `protein_id` → `protein_id_refseq`. Triple bridging covers both legacy locus tags and reassembly-era NCBI IDs.

## Recommended actions

1. **No action** on DE — four-experiment integration is complete.
2. **Skip** — XLSX originals are format duplicates.
3. **Verify** — MOESM1/MOESM2 supplementary PDFs should be scanned once more for any per-gene physiology measurements (photosynthetic quantum yield per gene, etc.) that could be encoded as `derived_metrics_table`; based on the abstract they appear to be summary figures/tables only (skip unless a per-gene row table is found).
4. **Note on p-values** — no p-value column in either CSV; `prefiltered: true` is the correct signal and DE edges get only log2FC (treated as significant-by-default per prefilter).

## Notes

- Gene IDs: tetu 2019 is the canonical "legacy locus + 2017 NCBI locus + protein_id" tri-column pattern. Resolution falls back through all three via `gene_id_mapping.json` tier-1/tier-2 logic.
- HDPE and PVC exposures have different durations (120 min vs 90 min) and different leachate dilutions — correctly encoded as separate experiments with distinct timepoints.
- The paper's growth and PSII-yield assays are sample-level and correctly absent from the KG.
