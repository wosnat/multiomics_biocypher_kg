# Ziegler 2025 (bioRxiv preprint)

**Citation:** Ziegler CA, Mullet JI, Coe A, Vo NN, Salcedo E, Arrigan DM, Parker SM, Chisholm SW (2025). The shared and distinct roles of *Prochlorococcus* and co-occurring heterotrophic bacteria in regulating community dynamics. bioRxiv preprint.
**DOI:** 10.1101/2025.09.25.678681
**Organism(s):** *Prochlorococcus* MED4; heterotrophs *Marinobacter*, *Alteromonas*, *Thalassospira*, *Pseudohoeflea* (synthetic community members)
**Topic:** RNA-seq of axenic MED4, MED4 with each individual heterotroph co-culture, and MED4 in a four-member synthetic community (Pro99 medium, continuous light 18 µE, 24 °C). Sampled at Day 2 (mid-exp), Day 4 (late-exp), and Day 5 (stationary transition). Two DE contrast families: (a) MED4 in each coculture/community vs axenic; (b) each heterotroph in community vs in single coculture with MED4. Absolute RNA/DNA/cell-count quantification underlies the analysis.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `2025.09.25.678681v1.full.pdf` | PDF | Main preprint (bioRxiv) | reference | — |
| `media-1.pdf` | PDF | Supplementary figures | reference | — |
| `MED4 DE genes table 6 media-7.csv` | CSV | Table 6: MED4 DE (coculture/community vs axenic) — `ID` (genome-tagged locus, e.g. `1_1_MED4_genome`), eggNOG/KEGG annotations, per-partner `{Partner}_LFC` + `{Partner}_padj` columns for Marinobacter, Alteromonas, Thalassospira (×1+×2 replicates), Pseudohoeflea, Community, per Day 2/4/5 | already in | — (drives 15 DE entries in paperconfig) |
| `heterptroph DE genes supp table 8 media-9.csv` | CSV | Table 8: heterotroph DE (community vs single coculture with MED4) — same format; rows cover all four heterotrophs | already in | — (drives 12 DE entries in paperconfig) |
| `media-6.csv` | CSV | Alternative export of the MED4 DE data — contains `organism`, per-day LFC/p columns; same underlying measurements as `media-7.csv` but in different layout | skip | — (duplicate of `media-7.csv` content) |
| `gene counts media-10.xlsx` | XLSX | Raw gene counts matrix (samples × genes) | skip | — (sample-level raw counts, out of KG scope) |
| `transcript counts media-10.xlsx` | XLSX | Raw transcript counts matrix (identical bytes to `media-10.xlsx` — duplicate of the gene counts export under a different name) | skip | — |
| `media-10.xlsx` | XLSX | Same-bytes duplicate of the two above (raw counts) | skip | — |
| `media-11.xlsx` | XLSX | Large supplementary table — raw/normalized expression or sample metadata (very large, >7 MB) | skip | — (likely sample-level) |
| `media-2.xlsx` | XLSX | Small supp table — heterotroph isolate strain list / genome accession table (Supp Table 1 in paper text) | skip | — (metadata, not per-gene) |
| `media-3.xlsx` | XLSX | Small supp table — cell-count / flow cytometry table (Supp Table 2) | skip | — (sample-level measurements) |
| `media-4.xlsx` | XLSX | Small supp table — metagenomic read-abundance table (Supp Table 3) | skip | — (sample-level) |
| `media-8.xlsx` | XLSX | Medium supp table — Thalassospira SNV / 2 Mb inversion comparison (Supp Table 4–5) | skip | — (strain-level genomics, not per-gene DE) |
| `paperconfig.yaml` | YAML | Active paperconfig | already in | — |
| `paperconfig_orig.yaml` | YAML | Pre-resolution copy of paperconfig | reference | — |

Note: no `_resolved.csv` / `_resolved_report.txt` files appear — the paperconfig uses `id_type: other`, which the resolver handles via multi-column heuristics at build time rather than persisted resolved CSVs.

## Classification

**Bucket B — new metrics / DE / resolution (want to add)**

27 DE analyses across MED4 + 4 heterotrophs are wired up, but heterotroph-side edges will dangle until the partner genomes are deployed. Marinobacter (MarRef v6) and Alteromonas (MarRef v6) are now reference-proteome-match organisms in the KG (added 2026-04-15/16), so two of the four heterotrophs are partly addressable; **Pseudohoeflea** and **Thalassospira** remain undeployed and Table 8 expression edges for them will not resolve. The community-vs-single contrasts also use a placeholder `treatment_taxid: 2742` for the multi-organism community (TODO in paperconfig). Action: deploy Pseudohoeflea and Thalassospira (use `/deploy-strain`), verify the heterotroph-side `id_type: other` resolution, and audit the swap of subject/partner on Table 8 analyses.

## Current paperconfig summary

- Experiments defined: 0 declared at the `experiments:` block level — every `statistical_analyses` entry is a free-standing RNA-seq DE with inline `treatment_condition` / `control_condition` / `organism` / `treatment_organism` / `treatment_taxid`. The `omics_adapter` synthesizes Experiment nodes from these (~15 MED4 ones + ~12 heterotroph ones = ~27 Experiment nodes).
- Statistical analyses (DE edges): 27 — MED4 vs each partner/community × 3 days (15 entries from Table 6) + each heterotroph community-vs-single × 3 days (12 entries from Table 8)
- Supplementary materials entry types: 27× `csv` (all pointing to either `MED4 DE genes table 6 media-7.csv` or `heterptroph DE genes supp table 8 media-9.csv`, each selecting different LFC/padj column pairs)
- Organisms covered: *Prochlorococcus* MED4 (cultured), and heterotrophs "Marinobacter", "Alteromonas", "Thalassospira", "Pseudohoeflea" (referenced by taxid). See notes on strain deployment below.
- Table scope(s): `all_detected_genes` inferred (full DE output per contrast; not prefiltered)
- Non-DE evidence: none
- ID resolution: `id_type: other` — IDs have genome-prefix form `{N}_{N}_{Organism}_genome` (e.g. `1_1_MED4_genome`). These require custom parsing; resolution relies on the `organism` column in the CSV and downstream heuristics in `build_id_lookup()`.

## Recommended actions

1. **No action** on the two primary DE CSVs — both are fully wired up.
2. **Skip** — the four `media-10.xlsx` / `gene counts media-10.xlsx` / `transcript counts media-10.xlsx` are raw sample-level count matrices (and at least one pair is literally duplicate bytes). Raw counts are out of scope per template rule.
3. **Skip** — `media-6.csv` appears to be a redundant export of the MED4 DE table already consumed via `media-7.csv`; integrating would double-count edges.
4. **Skip** — `media-2.xlsx`, `media-3.xlsx`, `media-4.xlsx`, `media-8.xlsx`, `media-11.xlsx` are strain metadata, flow-cytometry counts, metagenomic abundance, SNV analyses, or sample-level measurements — none supply per-gene evidence the KG supports.
5. **Skip** — `media-1.pdf` is reference (supplementary figures).
6. **Add organism** — the heterotroph strains (*Marinobacter*, *Pseudohoeflea*, *Thalassospira*) used here are **not** currently deployed in the KG as cultured genome strains. Table 8 expression edges will dangle unless these genomes are added to `cyanobacteria_genomes.csv` and the `_genome`-style locus tags are bridged. Current KG only has *Alteromonas* genomes (MIT1002/EZ55/HOT1A3). Track via `/deploy-strain` workflow.
7. **Change upload (candidate)** — the `treatment_taxid: 59919` (MED4) appears on heterotroph-side Table 8 analyses, with `treatment_organism: "Prochlorococcus MED4"`. Verify that swap of subject/partner is intentional on these analyses; otherwise the `Tests_coculture_with` edges will point at the wrong entity.

## Notes

- The Table-6 and Table-8 CSVs each contain many contrast columns. One physical CSV drives many paperconfig `csv` entries, each selecting its own `logfc_col` / `adjusted_p_value_col` — this is normal for wide DE tables.
- Gene IDs use a custom `{row}_{rep}_{Organism}_genome` form (e.g. `1_1_MED4_genome`, `23_1_Marinobacter_genome`). `id_type: other` signals that the resolver must combine ID + `organism` column rather than treating `ID` as a locus tag. Heterotroph rows may fail to resolve until the corresponding genomes are deployed (see action 6).
- The synthetic-community "Community" contrast uses `treatment_taxid: 2742` with a `# TODO` comment in the paperconfig noting that a proper identifier for a multi-organism treatment does not exist. The community is currently attached to MED4 Experiments via this placeholder.
- Two Thalassospira strains were sequenced; only one is used in paperconfig DE columns (`Thalassospira_1_*`). The `Thalassospira_2_*` column family, if present in the CSV, is unintegrated.
- No `_resolved.csv` / `_resolved_report.txt` files are present because `id_type: other` bypasses the standard tier-1 pre-resolution step.
