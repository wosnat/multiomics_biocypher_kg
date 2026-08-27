# He 2022

**Citation:** He X, Liu H, Long L, Dong J, Huang S (2022). Acclimation and stress response of *Prochlorococcus* to low salinity. *Frontiers in Microbiology* 13:1038136.
**DOI:** 10.3389/fmicb.2022.1038136
**GEO:** GSE195946
**Organism(s):** *Prochlorococcus* NATL1A, MED4
**Topic:** RNA-seq of cultures acclimated through five rounds of transfer in Pro99 medium with reduced salinity (28 psu) vs normal seawater salinity (34 psu); 21C, continuous light (10 uE/m2/s).

## Classification

**Integrated — full GEO tables (upgraded 2026-08 from the significant-only supplements)**

The DE source was **replaced** in the 2026-08 GEO processed-supplements pass
(plan item 1.1 of `plans/geo_paperconfig_updates.md`): the paper's S1/S2
"highly DE" subsets (MED4 31 rows / NATL1A 82 rows, `prefiltered: true`,
`significant_only`) are unwired; the full GEO `GSE195946_Group_*_DE_anno.txt`
tables (MED4 2,042 / NATL1A 2,239 rows, `all_detected_genes`,
`prefiltered: false`, real `FDR` as the adjusted-p column) are wired instead.
The old S1/S2 CSVs stay on disk for provenance but nothing references them.
Applying the paper's own cut (FDR <= 0.05 AND |logFC| >= 1) to the full tables
reproduces 30 / 81 significant genes vs S1/S2's 31 / 82.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `fmicb-13-1038136.pdf` | PDF | Main paper | reference | — |
| `paperconfig.yaml` | YAML | KG integration config | **wired** | — |
| `paperconfig_orig.yaml` | YAML | Pre-replacement config kept for provenance | reference | — |
| `geo_source.txt` | TXT | GEO provenance + column notes | reference | — |
| `GSE195946_Group_MED4-34-VS-MED4-28_DE_anno.txt` | TSV | Full MED4 DE table (2,042 rows; logFC/PValue/FDR) | **wired** | — |
| `GSE195946_Group_NATL1A-34-VS-NATL1A-28_DE_anno.txt` | TSV | Full NATL1A DE table (2,239 rows) | **wired** | — |
| `GSE195946_Group_*_DE_significant_anno.txt.gz` | TSV | GEO's own significant subsets (54 / 139 rows) | superseded | — |
| `GSE195946_*-all.fpkm.txt.gz` | TSV | Per-sample FPKM matrices | reference (Tier 4, D1: no DE computed by us) | — |
| `supp table legends.csv` | CSV | Legends text | reference | — |
| `Data Sheet 1.docx` | DOCX | Supplementary narrative | reference | — |
| `table s1 ... NATL1A .csv` / `table s2 ... MED4 .csv` | CSV | Old S1/S2 significant subsets | superseded (unwired) | — |

Auto-generated `*_resolved.txt` / `*_resolved_report.txt` are produced by prepare_data step 4.

## Current paperconfig summary

- Experiments defined: **2** (`salt_low_salinity_acclimation_28_natl1a_rnaseq`, `salt_low_salinity_acclimation_28_med4_rnaseq`)
- Statistical analyses (DE edges): **2** (one per strain; single `acclimated` timepoint, `growth_phase: acclimated_steady_state`)
- Supplementary materials entry types: `csv` × 2 (the GEO full tables)
- Table scope: `all_detected_genes`, `prefiltered: false`, `pvalue_threshold: 0.05`, `logfc_threshold: 1.0`
- ID resolution: `GeneID` (`gene-PMM0087` / `gene-NATL1_08781` form) resolves tier-1; MED4 1,948/2,042 (95.4%), NATL1A 2,190/2,239 (97.8%). Unresolved are dominated by GEO-custom RNA features with no Gene nodes; ~6 MED4 protein rows are recoverable via a strip-`gene-` heuristic (filed in `plans/backlog.md`).
- `adjusted_p_value_col: FDR` — a real BH-adjusted value, non-null on every row.

## Notes

- Sign convention verified against raw counts: positive logFC = up in the
  low-salinity (28 psu) treatment.
- One known duplicate resolution: `gene-RNA_41` → PMM0521 (frr), a
  pre-existing MED4 mapping quirk (5S-rRNA rrf/frr name collision), benign
  here (not_significant edge); filed in `plans/backlog.md`.
- The paper's Methods say DESeq while the GEO tables carry edgeR's `logCPM`
  signature; the config follows the paper (`test_type: DESeq`).
