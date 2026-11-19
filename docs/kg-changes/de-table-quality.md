# DE table quality (2026-11-19)

Spec: `docs/superpowers/specs/2026-11-18-de-table-quality-design.md`. Follow-ups: `plans/backlog.md`, section "DE table quality - follow-ups".

## What changed

Differential-expression experiments now say how far their table can be trusted: curated named limitations, the fold-change scale the table used, and per-time-point table statistics. Stored values (`log2_fold_change`, `adjusted_p_value`, `significant`) are unchanged.

## Properties

| Where | Property | Meaning |
|---|---|---|
| Experiment | `data_caveats` (str[]) | Curated closed vocabulary (8 values, `config/controlled_vocabularies.yaml`). Dense on DE experiments; `[]` is valid. |
| Experiment | `data_caveat_note` | Evidence, and which time points or tables a caveat covers. |
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

`scripts/validate_paperconfig.py` requires `data_caveats` + `data_caveats_reviewed_on` (and a note when a caveat is set) on every DE experiment, except those listed in `config/data_caveats_unreviewed.yaml`. That list only shrinks: a review sets the fields and removes the key; delete the file and its check when empty. 29 experiments in 12 papers are curated; 164 experiments in 37 papers remain on the allowlist. The validator also prints a QA report (statistics, candidate level columns, caveats to consider) from the shared `de_calls` code.

## Deferred

Backfill of the caveat review, level-column declarations for the ~11 papers that publish them, and the section 3.4 suspects: see `plans/backlog.md`, section "DE table quality - follow-ups".
