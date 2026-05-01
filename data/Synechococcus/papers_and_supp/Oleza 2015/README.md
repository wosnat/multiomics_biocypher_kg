# Oleza 2015 -- Synechococcus WH7803 + R. pomeroyi DSS-3 exoproteome

**Citation:** Christie-Oleza JA, Armengaud J, Bertrand P, Laez B, Boez B, Mariadassou S, Pelletier E (2015) Functional assessments of microbial interactions in the ocean using exoproteomics. *Proteomics* 15:3315-3325.

## Classification

**Bucket B — new metrics / DE / resolution (want to add)**

Per-protein NSAF-style abundance for WH7803 and Ruegeria pomeroyi DSS-3 across axenic vs coculture conditions. No fold-change table, but the NSAF values fit the numeric `derived_metrics_table` pattern shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`). Both organisms are already deployed (WH7803, DSS-3). Action: paperconfig with `derived_metrics_table` entries per strain x compartment (Synechococcus exoproteins, R. pomeroyi proteins). Companion paper Oleza 2017 already provides FC-bearing tables and is integrated; this 2015 dataset adds single-condition abundance context.

## Organism

- *Synechococcus* sp. WH7803 (cyanobacterium)
- *Ruegeria pomeroyi* DSS-3 (heterotrophic coculture partner)

## Data files

| File | Description |
|------|-------------|
| `pmic8102-sup-0001-tables1.xlsx` | Table S1 |
| `pmic8102-sup-0002-tables2.xlsx` | Table S2: Synechococcus WH7803 exoproteins (NSAF abundance) |
| `pmic8102-sup-0003-tables3.xlsx` | Table S3 |
| `pmic8102-sup-0004-tables4.xlsx` | Table S4 |
| `pmic8102-sup-0005-tables5.xlsx` | Table S5 |
| `Synechococcus coculture exoproteome Proteomics 2015.pdf` | Paper PDF |

## KG Integration Status

**Ready to wire** — no paperconfig yet, but unblocked.

**Reason:** Although the data reports NSAF / spectral-count abundance rather than fold-change (so no `Changes_expression_of` edges), the `derived_metrics_table` loader shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`) supports `value_kind: numeric` entries that fit this dataset directly. Action: add a paperconfig with one `derived_metrics_table` per (strain, compartment).

## Notes

- All supplementary data is in xlsx format only (no CSVs extracted)
- Table S2 has Synechococcus exoproteins with fold-change axenic vs co-culture and p-value for some entries, but the primary quantification is NSAF abundance
- Table S2b (R. pomeroyi proteins) has abundance only, no fold-change
- The companion paper Oleza 2017 extends this work with proper fold-change comparisons and is integrated separately
