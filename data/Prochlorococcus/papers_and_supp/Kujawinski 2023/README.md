# Kujawinski 2023

**Citation:** Kujawinski EB, Braakman R, Longnecker K, Becker JW, Chisholm SW, Dooley K, Kido Soule MC, Swarr GJ, Halloran K. Metabolite diversity among representatives of divergent *Prochlorococcus* ecotypes. mSystems, 2023, 8(5): e01261-22.
**DOI:** 10.1128/msystems.01261-22
**Organism(s):** *Prochlorococcus* MIT9301 (HLII), MIT0801 (LLI), MIT9313 (LLIV)
**Topic:** Mass-spectrometry-based metabolomics profiling of intracellular and extracellular metabolite pools across three ecologically divergent *Prochlorococcus* ecotypes, including a phosphorus-limitation condition for MIT9301. Reports ~35 intracellular + 18 extracellular metabolites with per-cell concentration (attomole/cell). This is a purely metabolomics paper — no per-gene differential expression, no proteomics.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `kujawinski-et-al-2023-metabolite-diversity-among-representatives-of-divergent-prochlorococcus-ecotypes.pdf` | PDF | Main article | reference | — |
| `msystems.01261-22-s0001.docx` | DOCX | Single combined supplementary file (likely Tables S1–S2 + supplementary methods) — metabolite concentrations per strain and growth condition | skip | Metabolite-level, not gene-level evidence; no direct gene targets in the data |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Skip** — Paper is metabolomics-only with no per-gene quantitation and no direct gene-to-metabolite mapping exported by the authors. The KG schema supports DE edges and per-gene DerivedMetrics; metabolite pool concentrations don't fit either pattern.
2. **No action** — Unless a future release maps the 35 intracellular + 18 extracellular metabolites back to biosynthesis genes (e.g., via KEGG pathway hits), there is nothing to add. The text identifies some pathway absences (e.g., glycine betaine biosynthesis genes missing in MIT9313) that could in principle be encoded as per-gene boolean flags, but that would require manual curation from the PDF — low ROI.
3. **Reference** — Keep the PDF and supplementary DOCX in-place for provenance. If future integration is desired (e.g., as metabolite nodes with strain presence/absence flags), it would require new node types the KG does not currently define.

## Notes

- Strains covered: MIT9301 (in KG), MIT9313 (in KG), MIT0801 (**not in KG** — LLI ecotype, NCBI taxid 1110371 or similar). If the team later decides to add MIT0801 and extract biosynthesis gene presence from the paper's pathway discussion, the DOCX would need manual parsing.
- No CSV/TSV tables — data is embedded in the single DOCX file, which would need extraction before any integration.
- Related context: metabolite release patterns from this paper overlap conceptually with Biller 2022 (vesicle metabolites) and Becker 2014 exometabolites.
