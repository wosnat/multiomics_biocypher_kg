# Waldbauer 2012

**Paper:** Waldbauer JR, Rodrigue S, Coleman ML, Chisholm SW (2012) Transcriptome and Proteome Dynamics of a Light-Dark Synchronized Bacterial Cell Cycle. PLoS ONE 7(8): e43432. doi:10.1371/journal.pone.0043432

**Organism:** Prochlorococcus MED4

**Topic:** Combined transcriptomics (RNA-seq) and quantitative proteomics (LC-MS ¹⁴N/¹⁵N SILAC) time course over a light-dark synchronized cell cycle. Measured paired mRNA-protein abundances for 312 genes every 2 hours, showing that transcript oscillations are broadly damped at the protein level.

## Integration

Integrated as a single `PAIRED_RNASEQ_PROTEOME` experiment with 6 per-gene derived metrics from Table S2 (phase, amplitude, lag, damping for 312 cycling genes). No DE comparison — diel-cycle summary statistics only.

## Files

| File | Description |
|------|-------------|
| `Waldbauer  2012 file.pdf` | Main paper PDF |
| `paperconfig.yaml` | KG integration config |
| `Table_S2.pdf` | Original per-gene cycling metrics table (9 pages) |
| `Table_S2_modified.csv` | Extracted CSV (312 rows × 9 cols) — built by `scripts/build_modified_csv/build_waldbauer2012_modified_csv.py` |
| `Table_S1.pdf` | Global cycling-detection stats (not per-gene; not ingested) |
| `Table_S3.pdf` | CBB/PPP pathway subset of Table S2 (redundant; not ingested) |
| `Table_S4.pdf` | Per-timepoint transcript coverage stats (not per-gene; not ingested) |
| `Table_S5.pdf` | Proteomics filtering summary (not per-gene; not ingested) |
| `Figure_S1.pdf` – `Figure_S10.pdf` | Supplementary figures |
| `Text_S1.pdf` | Supplementary text |

## Notes

- Waldbauer uses `PMED4_xxxxx` locus tags (JGI IMG draft annotation). The gene ID resolution pipeline maps these to canonical `PMM####` NCBI locus tags via `cache/data/Prochlorococcus/genomes/MED4/gene_id_mapping.json`.
- The paper also reports a raw per-timepoint time course in separate deposits; only the Table S2 derived metrics are ingested (per the non-DE-evidence extension spec, invariant #4: paired-modality derived metrics attach to a single `PAIRED_RNASEQ_PROTEOME` experiment).
