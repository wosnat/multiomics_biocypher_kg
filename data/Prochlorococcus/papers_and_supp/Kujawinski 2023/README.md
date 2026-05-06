# Kujawinski 2023

**Citation:** Kujawinski EB, Braakman R, Longnecker K, Becker JW, Chisholm SW, Dooley K, Kido Soule MC, Swarr GJ, Halloran K. Metabolite diversity among representatives of divergent *Prochlorococcus* ecotypes. mSystems, 2023, 8(5): e01261-22.
**DOI:** 10.1128/msystems.01261-22
**Organism(s):** *Prochlorococcus* MIT9301 (HLII), MIT0801 (LLI), MIT9313 (LLIV)
**Topic:** Mass-spectrometry-based metabolomics profiling of intracellular and extracellular metabolite pools across three ecologically divergent *Prochlorococcus* ecotypes, including a phosphorus-limitation condition for MIT9301. Reports ~35 intracellular + 18 extracellular metabolites with per-cell concentration (attomole/cell). This is a purely metabolomics paper - no per-gene differential expression, no proteomics.

## Classification

**Bucket A — integrated** (Phase 2 metabolomics; Metabolite + MetaboliteAssay nodes via `metabolite_assays_table`)

Wired up in `paperconfig.yaml` with 7 `metabolite_assays_table` entries: 6 numeric quant entries from `ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv` (3 strains × {intracellular, extracellular}) + 1 presence-flag entry from `table s2.csv`. All three strains — MIT9301, MIT0801, MIT9313 — are deployed (MIT0801 = `Prochlorococcus MIT0801`, LLI ecotype, GCF_000757865.1). `metabolite_aliases.yaml` and `*_resolved.csv` artifacts are present, indicating step 7 has run successfully. Glycine-betaine table S3 is the only machine-readable table not yet wired (small, optional follow-up).

Secondary: a small bucket-B opportunity exists if someone manually curates the paper's pathway-absence discussion (e.g., glycine betaine biosynthesis genes missing in MIT9313) into a per-strain Y/N gene-presence DerivedMetric — but ROI is low and the data is locked in DOCX/PDF prose.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `kujawinski-et-al-2023-metabolite-diversity-among-representatives-of-divergent-prochlorococcus-ecotypes.pdf` | PDF | Main article | reference | — |
| `msystems.01261-22-s0001.docx` | DOCX | Single combined supplementary file (Tables S1-S2 + supplementary methods) - metabolite concentrations per strain and growth condition | reference | Source of truth for the CSV derivations |
| `table s2.csv` | CSV | 35 intracellular metabolites + extraction efficiency + intracellular/extracellular presence flag | already in | `metabolite_assays_table` entry `presence_flags_table_s2` (`Assay_flags_metabolite`) |
| `table s3.csv` | CSV | Per-strain glycine betaine concentration (fg/cell) by clade + light intensity (MIT9301 HLII, MIT0801 LLI, MIT9313 LLIV) | not yet | Optional follow-up: small numeric quant table; one `metabolite_assays_table` entry per strain × light |
| `ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv` | CSV | KEGG-IDed metabolites x 13 sample columns (strain x P-status x light); concentrations in mol/cell | already in | 6 `metabolite_assays_table` entries (3 strains × {intracellular, extracellular}); `*_resolved.csv` + `metabolite_aliases.yaml` present |
| `metabolite_aliases.yaml` | YAML | Per-paper alias overrides feeding step 6/7 metabolite resolver | reference | — |
| `*_resolved.csv` / `*_resolution_report.json` | CSV/JSON | Step 7 outputs (consumed by `metabolite_assay_adapter`) | — | auto-generated |
| `legends.txt` | TXT | Empty placeholder | reference | — |

## Current paperconfig summary

- Experiments defined: per-strain × per-compartment Phase 2 metabolomics Experiments (`omics_type: METABOLOMICS`).
- Statistical analyses (DE edges): 0 — paper has no per-gene DE.
- Supplementary materials entry types: 7× `metabolite_assays_table` (6 numeric quant on `ChisholmPro_cellSpecific_KEGGexport`, 1 boolean presence on `table s2`).
- Organisms covered: *Prochlorococcus* MIT9301, MIT0801, MIT9313 (all deployed via `cyanobacteria_genomes.csv`).
- Companion files: `metabolite_aliases.yaml` (alias overrides) + `*_resolved.csv` (step 7 outputs).
- Not yet wired: `table s3.csv` glycine-betaine quant (small follow-up).

## Recommended actions

1. **No action** — Phase 2 metabolomics is integrated for the `ChisholmPro_cellSpecific_KEGGexport` quant tables (6 entries) and the `table s2` presence flags (1 entry). All three strains resolve.
2. **Optional follow-up** — wire `table s3.csv` (glycine betaine per strain × light) as one or more `metabolite_assays_table` numeric quant entries. Small.
3. **Optional bucket-B follow-up** — pathway absences in the prose (e.g., glycine betaine biosynthesis genes missing in MIT9313) could be encoded as per-gene boolean DerivedMetrics, but require manual extraction from DOCX — low ROI.
4. **Reference** — Keep the PDF and supplementary DOCX in-place for provenance.

## Notes

- Strains covered: MIT9301, MIT0801, MIT9313 — all in KG (`cyanobacteria_genomes.csv`). MIT0801 is `Prochlorococcus MIT0801` (LLI, GCF_000757865.1).
- Related context: metabolite release patterns from this paper overlap conceptually with Biller 2022 (vesicle metabolites) and Becker 2014 exometabolites.
