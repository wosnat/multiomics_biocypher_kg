# DE table quality (2026-11-19)

Spec: `docs/superpowers/specs/2026-11-18-de-table-quality-design.md`. Follow-ups: `plans/backlog.md`, section "DE table quality - follow-ups".

## What changed

Differential-expression experiments now say how far their table can be trusted: curated named limitations, what the KG did to the published table, the fold-change scale the table used, and per-time-point table statistics. Stored values (`log2_fold_change`, `adjusted_p_value`, `significant`) are unchanged.

## Properties

| Where | Property | Meaning |
|---|---|---|
| Experiment | `data_caveats` (str[]) | Problems visible in the published data. Curated closed vocabulary (5 values, `config/controlled_vocabularies.yaml`). Dense on DE experiments; `[]` is valid. |
| Experiment | `fold_change_provenance` (str[]) | What the KG did to the published table: `computed_by_kg`, `sign_reversed_by_kg`, `significance_rule_differs_from_paper`. Not a problem with the data. Dense on DE experiments; `[]` = values, signs and significance as the paper reports them. |
| Experiment | `data_notes` | General free-text notes (any experiment, in `experimentFullText`): evidence for a caveat or a provenance value, and which time points or tables it covers; or a concern that rests only on the Methods of the paper, which is written here and not tagged. |
| Experiment | `data_caveats_reviewed_on` | Review date. `[]` with a date = reviewed, none apply; `[]` without = not reviewed. |
| Experiment | `fold_change_type` | `log2` / `linear` / `signed_linear`: scale of the table the KG read. Sparse (absent without a DE table). |
| Experiment | `low_signal_basis` | `condition_at_column_minimum` / `mean_in_bottom_decile`. Sparse; only where level columns are declared. |
| Changes_expression_of | `low_signal_reason` | Same values; sparse. |
| Experiment | six `time_point_*` arrays | Post-import, parallel to `time_point_orders`: `negative_fractions`, `median_log2fc`, `max_abs_log2fc`, `large_fc_not_significant` (|log2FC| > 5, not significant), `extreme_fc_count` (|log2FC| > 10), `low_signal_count`. |

## Reading the statistics by `table_scope`

The statistics are computed over the experiment's published rows on every scope; they are gene-agnostic.

- `all_detected_genes`: they describe the measured genome. A negative fraction far from 0.5 or a large median suggests unnormalized conditions.
- `significant_only`, `filtered_subset`, `top_n`: they describe the published subset. `time_point_large_fc_not_significant` is 0 by construction in a significant-only table. `time_point_extreme_fc_count` stays meaningful.
- `time_point_low_signal_count` is 0 wherever no level columns are declared (all papers until the level-column backlog item is done; `plans/backlog.md`, section "DE table quality - follow-ups").

## Curation and the allowlist

`scripts/validate_paperconfig.py` requires `data_caveats` + `data_caveats_reviewed_on` (and a note when a caveat is set) on every DE experiment, except those listed in `config/data_caveats_unreviewed.yaml`. That list only shrank: a review set the fields and removed the key. The backfill finished on 2026-11-19 and the file was deleted, so every DE experiment must now be reviewed (193 experiments in 44 papers: 17 with a data caveat, 29 with a provenance value, 16 with a Methods note and no tag).

The rule (2026-11-19): the KG presents the data and the statistics as the paper reports them. A tag needs a problem measured or visible in the data (zero-count fold changes, a constant offset between conditions). A Methods choice that may not be sound (imputation, normalization the authors describe) is written in `data_notes` only. What the KG itself did (computed or negated a fold change, applied its own significance rule) is `fold_change_provenance`, not a caveat. A builder in `scripts/build_modified_csv/` that computes or negates fold changes adds the matching provenance value to every experiment it feeds. The validator still reads the file if it exists (a missing file means an empty list), so a later vocabulary change can re-queue experiments for review. The validator also prints a QA report (statistics, candidate level columns, caveats to consider) from the shared `de_calls` code.

## Deferred

Backfill of the caveat review, level-column declarations for the ~11 papers that publish them, and the section 3.4 suspects: see `plans/backlog.md`, section "DE table quality - follow-ups".
