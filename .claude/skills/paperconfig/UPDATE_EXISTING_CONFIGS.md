# Auditing the Fold-Change and Significance Fields of an Existing paperconfig

Use this when re-checking an existing paperconfig, for example after a corpus audit or when a
paper's significant count looks wrong. The 2026-10-31 audit of all DE papers
(`plans/signed_linear_fold_change.md`) found every class of error listed below at least once.

## The fields

| Field | Type | Meaning |
|-------|------|---------|
| `fold_change_type` | `log2` (default) · `linear` · `signed_linear` | How the `logfc_col` values are written. `linear` = unsigned ratio (0.5 = 2-fold down); `signed_linear` = +3 is 3-fold up, -3 is 3-fold down. Under `signed_linear` a value in (-1, 1) is an error |
| `prefiltered` | boolean | The table lists only the genes the authors called DE; every row is significant |
| `pvalue_asterisk_in_logfc` | boolean | Significance is marked by `*` in the fold-change cell |
| `asterisk_min_level` | positive int (default 1) | Minimum stars that count as significant; requires `pvalue_asterisk_in_logfc` |
| `pvalue_threshold` | float in (0, 1] (default 0.05) | The paper's p / q cutoff |
| `logfc_threshold` | float >= 0, **no default** | The paper's \|log2FC\| cutoff; `0` = no fold cutoff. Required unless `prefiltered` or `pvalue_asterisk_in_logfc` is set |

The adapter calls significance in this order: `prefiltered` → asterisks → an empty p-value cell
(not significant) → `pvalue_threshold` and `logfc_threshold`. A table with no p-value column at
all is judged on the fold change alone.

## Checklist, per analysis

1. **Scale.** Read the table legend and look at the values.
   - Values never between -1 and 1, legend "positive and negative values indicate increase and
     decline" → `fold_change_type: signed_linear`. The adapter warns when a `log2` table has no value
     in (-1, 1) and max |x| > 20, but a signed table with small effects (Lindell 2007, max 14) slips
     under that, so read the legend.
   - Only positive values with 1 in the middle → `linear`.
   - Any other convention needs a `scripts/build_modified_csv/` builder that writes a `_modified.csv`
     with a true log2 column (see that directory's README):
     - minus-sign-marked ratios (Ma 2022: -0.55 = ratio 0.55)
     - ratios inflated by a tiny pseudocount (Aharonovich HOT1A3: 1 read vs 0 = 101-fold)
     - a zero-read floor (Beliaev 2014)
     - a corrupt ratio column (Fuszard 2012 NATL2A)

     Then point `filename` at the new file and rerun
     `uv run python -m multiomics_kg.download.resolve_paper_ids --papers "<Paper>" --force`.
2. **Is the table prefiltered?** Count its rows against the genome size and the paper's reported
   DE count. A whole-genome table, or one filtered only at P < 0.05 when the paper's criterion is
   stricter, is *not* `prefiltered`.
3. **Criterion.** Copy the methods' criterion exactly:
   - q/p only → `logfc_threshold: 0`
   - "fold change > 2" → `1.0`
   - "|log2FC| > 0.8" → `0.8`
   - a p column rounded to 2 decimals with "p < 0.05" → `pvalue_threshold: 0.049`
4. **Asterisks.** If the legend grades the stars and the paper's DE call is the stricter band, set
   `asterisk_min_level: 2`.
5. **Zero-count rows.** If a table publishes fold changes against zero counts that cannot be
   recomputed (VST-transformed values, no raw counts), leave a comment in the paperconfig saying so,
   as Capovilla 2023, Anjur-Dietrich 2026 and Huang 2021 do.
6. **Verify.** Run the validator. Then build the paper's edges offline and compare the significant
   count with the number the paper reports. A gap of more than a few percent means one of the
   choices above is wrong.

```bash
uv run python scripts/validate_paperconfig.py "data/<Genus>/papers_and_supp/<Author Year>/paperconfig.yaml"
```

The validator errors on:
- a missing `logfc_threshold` where one is required
- an unknown `fold_change_type`
- a `signed_linear` column holding a value in (-1, 1)
- an `asterisk_min_level` that is not a positive int, or is set without `pvalue_asterisk_in_logfc`
- a `pvalue_threshold` outside (0, 1], or a negative `logfc_threshold`

It warns when `prefiltered` is combined with thresholds or asterisks.

Log every change to an existing paper in `CHANGELOG.md` → `### Data`: what the KG asserted before
and what it asserts now.
