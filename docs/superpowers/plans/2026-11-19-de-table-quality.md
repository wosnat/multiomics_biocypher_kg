# DE Table Quality Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish, on every DE `Experiment`, curated data caveats, per-time-point table statistics, a low-signal count and `fold_change_type`. Also publish a per-edge `low_signal_reason`, so the MCP can warn before an analysis reads values.

**Architecture:** New pure module `multiomics_kg/utils/de_calls.py` holds the per-row logic: reading a row, converting the scale, calling significance, the low-signal rule, and the per-time-point statistics.
- `omics_adapter.py` uses it to build edges and emit the new Experiment and edge properties.
- `scripts/validate_paperconfig.py` uses the same module for its QA report, plus new curation checks with an allowlist.
- Post-import Cypher computes the stored statistics. A KG test pins them to the Python statistics.

**Tech Stack:** Python 3 (pandas, numpy), BioCypher adapters, Neo4j Cypher (`scripts/post-import.sh` + `scripts/post-import.cypher`), pytest.

**Spec:** `docs/superpowers/specs/2026-11-18-de-table-quality-design.md`. Where sections disagree, the dated amendments win: §7.4–7.7, §9, §10.

## Global Constraints

- Line B: stored `log2_fold_change`, `adjusted_p_value` and `significant` values do **not** change. The expression-edge snapshot before and after the rebuild must match edge counts per paper (`/omics-edge-snapshot`).
- House rules:
  - R1: snake_case values.
  - R3: no stored threshold over a value the node or edge already carries. There is no edge flag for |log2FC| > 10.
  - R5: no bool, and no third `'true'`-only marker.
  - Strings go through `clean_text` / `_clean_str`. No `|` or `'` in any value.
- Every post-import property is declared in `config/schema_config.yaml` with a `# post-import:` comment. `scripts/post-import.sh` and `scripts/post-import.cypher` carry identical Cypher.
- Every new value set is an entry in `config/controlled_vocabularies.yaml`, with `value_descriptions` (researcher-facing, no build provenance in the text; provenance goes in a `# provenance:` comment).
- Thresholds, copied verbatim:
  - large FC = |log2FC| > 5;
  - extreme FC = |log2FC| > 10;
  - bottom decile = mean ≤ 10th percentile of the table's mean column;
  - rule 2 also needs |log2FC| > 5.
- Empty level cell = no signal (rule 1: at the column minimum; rule 2: in the bottom decile).
- The statistics are computed on **every** DE experiment, whatever its `table_scope`.
- Out of scope (spec §7.6 / §10):
  - `array_design`, `shares_probe_with`, `paralogs_share_one_array_probe`;
  - level-column declarations for the ~11 papers;
  - the §3.4 suspects;
  - the allowlist backfill;
  - F2–F6.
- Concurrent sessions edit this clone: `git add` explicit paths only, never `-A`.

## Review Focus

1. A level column holding non-numeric text (`"n.d."`, `"-"`): treat it as empty, so no signal. It must not crash or count as a value. Test in Task 2.
2. A level column with every cell empty or identical: every row is "at minimum". The rule must still return a mask, with no division by zero and no NaN percentile crash. Test in Task 2.
3. An analysis declaring `treatment_level_cols` without `control_level_cols`, or both per-arm and mean columns: the validator errors and the adapter raises `ValueError`, never a silent half-rule. Test in Tasks 2 and 6.
4. A time point whose edges are all `log2_fold_change = 0` (no direction): `negative_fraction = 0.0`, `median = 0.0`, and every such edge counts as not significant. Python and Cypher must agree. Test in Task 2, pinned live in Task 9.
5. An experiment present in the allowlist that also carries `data_caveats_reviewed_on`: the validator errors as a stale entry, so the list can only shrink. Test in Task 6.

---

### Task 1: Shared per-row DE calls, with the adapter moved onto them

**Files:**
- Create: `multiomics_kg/utils/de_calls.py`
- Modify: `multiomics_kg/adapters/omics_adapter.py`. Replace the bodies of `_convert_fold_change` (line 77) and `_check_significance` (line 673). Replace the table read and the per-row fold-change / p-value / significance block in `_load_and_create_edges` (lines 757–770 and 880–945).
- Test: `tests/test_de_calls.py`

**Interfaces:**
- Produces:
  - `read_de_table(filename: str, analysis: dict) -> tuple[pd.DataFrame, bool]`: the dataframe plus `used_resolved`. It probes `get_resolved_path(filename)` first, with the same logic as today.
  - `convert_fold_change(fc: float, fold_change_type: str | None) -> float | None`
  - `check_significance(log2fc: float | None, pvalue: float | None, asterisk_significant: bool | None, analysis: dict, *, default_pvalue_threshold: float | None = None, default_logfc_threshold: float | None = None, p_missing: bool = False) -> str`: returns `'significant' | 'not significant' | 'unknown'`.
  - `@dataclass(frozen=True) RowCall(raw_fc: float, log2_fold_change: float | None, adjusted_p_value: float | None, significant: str)`
  - `call_row(row: dict, analysis: dict, columns: set[str], *, default_pvalue_threshold=None, default_logfc_threshold=None) -> RowCall | None`. It returns `None` when the fold-change cell is empty, `NA` or non-numeric. It returns a `RowCall` with `log2_fold_change=None` when the conversion is invalid; the adapter appends `raw_fc` and then skips that row, as today.
  - `expression_status(significant: str, log2fc: float | None) -> str`: `'significant_up' | 'significant_down' | 'not_significant'`. Mirrors post-import: significant with log2FC > 0 → up, < 0 → down, else not significant.
- `omics_adapter._convert_fold_change` stays importable: it delegates to `convert_fold_change`, because existing tests import it. `OMICSAdapter._check_significance` stays as a one-line wrapper passing `self.default_*`.

- [ ] **Step 1: Write the failing tests** in `tests/test_de_calls.py`:
  - `test_call_row_asterisk_levels`: with `{"logfc_col": "fc", "pvalue_asterisk_in_logfc": True, "asterisk_min_level": 2, "pvalue_threshold": 0.05}`:
    - `"**1.5"` → `significant == "significant"`, `adjusted_p_value == 0.05`;
    - `"1.5 *"` → `"not significant"`, `adjusted_p_value == 1.0`.
  - `test_call_row_empty_pvalue_is_not_significant`: `adjusted_p_value_col: "padj"` present in `columns`, the row's `padj` is NaN, `logfc_threshold: 1`, fc 3.0 → `"not significant"`.
  - `test_call_row_unparseable_fc_returns_none`: fc `"NA"`, `""`, `"abc"` and `inf` → `None`.
  - `test_call_row_linear_conversion`: `fold_change_type: "linear"`, fc `0.25` → `log2_fold_change == -2.0`; fc `0` → `RowCall` with `raw_fc == 0.0` and `log2_fold_change is None`.
  - `test_check_significance_explicit_zero_threshold`: `logfc_threshold: 0`, `pvalue_threshold: 0.05`, fc 0.01, p 0.01 → `"significant"`.
  - `test_expression_status_zero_fc`: `expression_status("significant", 0.0) == "not_significant"`.
- [ ] **Step 2: Run them and confirm they fail.** Run `uv run pytest tests/test_de_calls.py -v`. Expected: FAIL (`ModuleNotFoundError: multiomics_kg.utils.de_calls`).
- [ ] **Step 3: Implement the five functions and `RowCall`** in `de_calls.py`. Move the existing code; don't rewrite it. Keep its behaviour exactly, including:
  - the asterisk strip;
  - the `is not None` threshold defaults;
  - the asterisk placeholder p-value `pvalue_threshold or default or 0.05`;
  - `p_missing = asterisk is None and pval is None and p_col declared and present`.
- [ ] **Step 4: Point `omics_adapter.py` at the module.** The row loop keeps the gene-id handling, `significance_mode`, direction, time point and edge-id logic.
- [ ] **Step 5: Run the new and existing adapter tests.** Run `uv run pytest tests/test_de_calls.py tests/test_omics_adapter.py tests/test_omics_adapter_fold_change_type.py tests/test_omics_duplicate_rows.py tests/test_omics_adapter_experiment.py -v`. Expected: all PASS, with no edits to the existing test files.
- [ ] **Step 6: Commit.**
```bash
git add multiomics_kg/utils/de_calls.py multiomics_kg/adapters/omics_adapter.py tests/test_de_calls.py
git commit -m "refactor(omics): move per-row DE calls into utils/de_calls (behaviour-preserving)"
```

### Task 2: Low-signal rule and per-time-point statistics in `de_calls`

**Files:**
- Modify: `multiomics_kg/utils/de_calls.py`
- Test: `tests/test_de_calls.py`

**Interfaces:**
- Produces:
  - `LOW_SIGNAL_CONDITION = "condition_at_column_minimum"`, `LOW_SIGNAL_MEAN = "mean_in_bottom_decile"`
  - `LARGE_FC = 5.0`, `EXTREME_FC = 10.0`
  - `level_basis(analysis: dict) -> str | None`. Returns `LOW_SIGNAL_CONDITION` when `treatment_level_cols` and `control_level_cols` are both present, `LOW_SIGNAL_MEAN` when only `mean_level_cols` is, and `None` when none are. It raises `ValueError` when only one arm is given, or when per-arm and mean columns are both present. A single string is accepted as a one-item list.
  - `low_signal_mask(df: pd.DataFrame, analysis: dict) -> list[bool] | None`: one entry per df row, in df order. `None` when `level_basis` is `None`.
  - `low_signal_reason(basis: str | None, in_mask: bool, log2fc: float | None) -> str | None`:
    - condition basis → `basis` if `in_mask`;
    - mean basis → `basis` if `in_mask` and `abs(log2fc) > LARGE_FC`;
    - else `None`.
  - `@dataclass(frozen=True) TimePointStats(total: int, negative_fraction: float, median_log2fc: float, max_abs_log2fc: float, large_fc_not_significant: int, extreme_fc_count: int, low_signal_count: int)`
  - `time_point_stats(rows: Iterable[tuple[float, str, str | None]]) -> TimePointStats`. Each row is `(log2_fold_change, expression_status, low_signal_reason)`. Definitions:
    - negative fraction = count(fc < 0) / total;
    - median = the interpolated median (equal to Cypher `percentileCont(fc, 0.5)`);
    - large-not-significant = count(|fc| > 5 and status == `not_significant`);
    - extreme = count(|fc| > 10);
    - low signal = count(reason is not None).

  An empty input returns `total=0` with 0.0 / 0 fields.

Algorithm for `low_signal_mask`. A cell counts as **no signal** when `pd.to_numeric(..., errors="coerce")` is NaN, or when the cell equals its column's minimum (the minimum taken over the column's numeric non-NaN values in the whole table).
- Rule 1: the row is True when, for either arm, every column of that arm is no-signal.
- Rule 2: average the mean columns per row (skipna). The row is True when that value is NaN, or is ≤ `np.nanpercentile(values, 10)` (linear).
- A column with no numeric values: every cell is no signal (rule 1), or every row is NaN and so in the bottom decile (rule 2).

- [ ] **Step 1: Write the failing tests:**
  - `test_level_basis_kinds_and_errors`:
    - both arms → condition;
    - mean only → mean;
    - none → `None`;
    - treatment only → `ValueError`;
    - arms plus mean → `ValueError`;
    - `"baseMean"` as a string → mean.
  - `test_mask_rule1_minimum_and_empty`: `t=[0.01, 5, 7, None]`, `c=[3, 0.01, 4, 2]` → `[True, True, False, True]`.
  - `test_mask_rule1_replicates_need_all_at_minimum`: treatment columns `r1=[0, 0, 2]`, `r2=[0, 3, 0]`, control `c=[5, 6, 7]` → `[True, False, False]`. Row 0: both replicates are at their minimum (and so is the control). Row 1: `r2` is not at its minimum. Row 2: `r1` is not at its minimum.
  - `test_mask_constant_column_is_all_minimum`: `t=[1, 2, 3]`, `c=[5, 5, 5]` → `[True, True, True]`, with no exception.
  - `test_mask_non_numeric_text_is_empty`: `t=["n.d.", 4, 5]`, `c=[3, 4, 5]` → row 0 True.
  - `test_mask_all_empty_column`: `t=[None, None]`, `c=[1, 2]` → `[True, True]`, with no exception.
  - `test_mask_rule2_bottom_decile`: mean `[1..20]` → True for values ≤ 2.9 (`nanpercentile` 10 = 2.9), i.e. rows with 1 and 2.
  - `test_low_signal_reason_rule2_needs_large_fc`: `(MEAN, True, 5.0) is None`; `(MEAN, True, -5.01) == MEAN`; `(CONDITION, True, 0.1) == CONDITION`.
  - `test_time_point_stats_values`: rows `(-1, "not_significant", None)`, `(2, "significant_up", None)`, `(-6, "not_significant", "condition_at_column_minimum")`, `(11, "significant_up", None)`, `(0, "not_significant", None)` → total 5, negative_fraction 0.4, median 0.0, max_abs 11.0, large_fc_not_significant 1, extreme 1, low_signal 1.
  - `test_time_point_stats_all_zero`: three `(0.0, "not_significant", None)` rows → negative_fraction 0.0, median 0.0.
- [ ] **Step 2: Run them and confirm they fail.** `uv run pytest tests/test_de_calls.py -v`. Expected: the new tests FAIL with ImportError.
- [ ] **Step 3: Implement** the constants and the four functions above.
- [ ] **Step 4: Run them and confirm they pass.** `uv run pytest tests/test_de_calls.py -v`. Expected: PASS.
- [ ] **Step 5: Commit.**
```bash
git add multiomics_kg/utils/de_calls.py tests/test_de_calls.py
git commit -m "feat(de_calls): low-signal rule and per-time-point table statistics"
```

### Task 3: Vocabularies and schema declarations

**Files:**
- Modify: `config/controlled_vocabularies.yaml` (new entries next to `Experiment.control_timing`, line 786, and `Changes_expression_of.expression_status`, line 987)
- Modify: `config/schema_config.yaml` (Experiment block, near line 71; `changes_expression_of` block, near line 448)
- Test: `tests/test_controlled_vocab.py` (existing generic tests) + `tests/test_de_quality_vocab.py`

**Interfaces:**
- Produces vocabulary ids, read later through `load_vocabularies()[id].values`:
  - `Experiment.data_caveats`: `string_array`, closed. Values are the 8 in spec §4 item 1, without `paralogs_share_one_array_probe`:
    - `fold_change_inflated_at_low_counts`
    - `no_biological_replicates`
    - `conditions_not_normalized_to_each_other`
    - `comparison_direction_unclear`
    - `table_may_be_a_different_comparison`
    - `significance_calls_differ_from_paper`
    - `fold_change_recomputed_by_kg`
    - `fold_change_sign_reversed_by_kg`

    `value_descriptions` are copied verbatim from the spec table. The description says: dense on DE experiments; `[]` with `data_caveats_reviewed_on` = reviewed and clean; `[]` without it = not yet reviewed.
  - `Experiment.fold_change_type`: string, closed, sparse; `[log2, linear, signed_linear]`. The description says it is the scale of the table the KG read; for papers whose table a `_modified.csv` builder rewrites, that is the derived table. The stored `log2_fold_change` is always log2.
  - `Experiment.low_signal_basis`: string, closed, sparse; `[condition_at_column_minimum, mean_in_bottom_decile]`. The `mean_in_bottom_decile` description states the |log2FC| > 5 cutoff and why it is needed.
  - `Changes_expression_of.low_signal_reason`: string, closed, sparse; same two values.
- Schema (`schema_config.yaml`):
  - Experiment: `data_caveats: str[]`, `data_caveat_note: str`, `data_caveats_reviewed_on: str`, `fold_change_type: str`, `low_signal_basis: str`.
  - Post-import (each with a `# post-import:` comment): `time_point_negative_fractions: float[]`, `time_point_median_log2fc: float[]`, `time_point_max_abs_log2fc: float[]`, `time_point_large_fc_not_significant: int[]`, `time_point_extreme_fc_count: int[]`, `time_point_low_signal_count: int[]`. Each comment says that in a significant-only, filtered or top-n table the numbers describe the published subset.
  - `changes_expression_of`: `low_signal_reason: str`.

- [ ] **Step 1: Write the failing test** `tests/test_de_quality_vocab.py::test_de_quality_vocabularies`. It asserts:
  - the four ids exist;
  - `set(v["Experiment.data_caveats"].values)` equals the 8 values above;
  - `"paralogs_share_one_array_probe" not in` them;
  - every value has a description.
- [ ] **Step 2: Run it and confirm it fails.** `uv run pytest tests/test_de_quality_vocab.py -v`. Expected: FAIL (KeyError).
- [ ] **Step 3: Add the four vocabulary entries and the schema declarations above.**
- [ ] **Step 4: Run it and the generic vocabulary tests.** `uv run pytest tests/test_de_quality_vocab.py tests/test_controlled_vocab.py tests/test_controlled_vocabulary_adapter.py -v`. Expected: PASS. That includes `test_descriptions_carry_no_build_provenance`.
- [ ] **Step 5: Commit.**
```bash
git add config/controlled_vocabularies.yaml config/schema_config.yaml tests/test_de_quality_vocab.py
git commit -m "feat(vocab): data_caveats, fold_change_type, low_signal_basis/_reason + schema"
```

### Task 4: Adapter emits the new Experiment and edge properties

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py`: `get_nodes` (Experiment props, lines 404–466) and `_load_and_create_edges`
- Test: `tests/test_omics_adapter_de_quality.py`

**Interfaces:**
- Consumes: `level_basis`, `low_signal_mask`, `low_signal_reason` (Task 2); the vocabularies (Task 3).
- Produces Experiment props, only when `exp_key in de_experiments`:
  - `data_caveats`: always a list (`[]` default), each value cleaned;
  - `fold_change_type`: from the experiment's analyses, default `"log2"`;
  - `low_signal_basis`: only when some analysis declares level columns.

  Props on any experiment, when present: `data_caveat_note`, `data_caveats_reviewed_on` (sparse, cleaned).

  Edge prop: `low_signal_reason`, set only when not `None`.
- `get_nodes` raises `ValueError(f"{experiment_id}: analyses declare different fold_change_type values {sorted(set)}")` when an experiment's analyses disagree. It raises the same way, listing the bases, when they disagree on a non-`None` `level_basis`. An analysis without level columns does not conflict with one that has them.

- [ ] **Step 1: Write the failing tests.** Use the `_run_adapter` / `tmp_path` paperconfig pattern from `tests/test_omics_adapter_fold_change_type.py`.
  - `test_experiment_carries_caveats_and_fold_change_type`: the experiment has `data_caveats: [fold_change_inflated_at_low_counts]`, `data_caveat_note: "x"`, `data_caveats_reviewed_on: "2026-11-19"`; the analysis has `fold_change_type: linear`. The props equal those values, and `fold_change_type == "linear"`.
  - `test_de_experiment_defaults`: no caveat fields → `data_caveats == []`, `fold_change_type == "log2"`, and no `data_caveats_reviewed_on` / `low_signal_basis` keys.
  - `test_non_de_experiment_has_no_de_quality_props`: an experiment with no analyses has none of `data_caveats`, `fold_change_type`.
  - `test_mixed_fold_change_type_raises`: two analyses on one experiment, `log2` and `linear` → `pytest.raises(ValueError, match="fold_change_type")`.
  - `test_low_signal_reason_on_edges`: a CSV with `t=[0.01, 5, 6]`, `c=[4, 5, 0.01]`, `fc=[−8.6, 0.0, 9.1]`, with `treatment_level_cols: [t]`, `control_level_cols: [c]` declared. Edges for rows 1 and 3 carry `low_signal_reason == "condition_at_column_minimum"`; row 2 has none. The experiment's `low_signal_basis == "condition_at_column_minimum"`.
  - `test_low_signal_survives_unresolved_rows`: as above, but row 1 has an empty gene id. Row 3's edge is still flagged, which shows the mask is indexed by df row, not by edge.
- [ ] **Step 2: Run them and confirm they fail.** `uv run pytest tests/test_omics_adapter_de_quality.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement.** In `_load_and_create_edges`:
  - compute `basis = level_basis(analysis)` and `mask = low_signal_mask(df, analysis)` once, before the loop;
  - inside the loop, use the loop index `idx` to read `mask[idx]`.
- [ ] **Step 4: Run the new and adapter tests.** `uv run pytest tests/test_omics_adapter_de_quality.py tests/test_omics_adapter*.py tests/test_de_calls.py -v`. Expected: PASS.
- [ ] **Step 5: Run the adapter over every real paperconfig, so no ValueError fires on today's data.** Run `uv run python -c "from multiomics_kg.adapters.omics_adapter import MultiOMICSAdapter as M; a=M(config_list_file=['data/Prochlorococcus/papers_and_supp/paperconfig_files.txt','data/Synechococcus/papers_and_supp/paperconfig_files.txt']); a.download_data(); print(sum(1 for n in a.get_nodes() if n[1]=='experiment'))"`. This mirrors `create_knowledge_graph.py:159`. Expected: an experiment count is printed, with no traceback.
- [ ] **Step 6: Commit.**
```bash
git add multiomics_kg/adapters/omics_adapter.py tests/test_omics_adapter_de_quality.py
git commit -m "feat(omics): Experiment data_caveats/fold_change_type/low_signal_basis + edge low_signal_reason"
```

### Task 5: Post-import statistics, low-signal rollup and the dense `data_caveats`

**Files:**
- Modify: `scripts/post-import.sh` (defaults block ~line 410, aggregation ~line 425, dense-list block ~line 356). Make the identical change in `scripts/post-import.cypher` (aggregation near line 370).
- Test: `tests/test_post_import_de_quality.py` (static check)

**Interfaces:**
- Produces Experiment arrays, parallel to `time_point_orders` and in the same order: the six `time_point_*` arrays declared in Task 3. Each defaults to `[]` in the existing defaults block.
- Produces `data_caveats = coalesce(e.data_caveats, [])` for every Experiment with `e.fold_change_type IS NOT NULL`. Only DE experiments carry `fold_change_type`, so non-DE experiments stay absent.

Cypher, to add to the existing per-time-point `WITH` (aggregation) and its matching `collect(...)`:
```cypher
     count(CASE WHEN r.log2_fold_change < 0 THEN 1 END) AS neg,
     percentileCont(r.log2_fold_change, 0.5) AS med,
     max(abs(r.log2_fold_change)) AS maxabs,
     count(CASE WHEN abs(r.log2_fold_change) > 5 AND r.expression_status = 'not_significant' THEN 1 END) AS large_ns,
     count(CASE WHEN abs(r.log2_fold_change) > 10 THEN 1 END) AS extreme,
     count(r.low_signal_reason) AS low_sig
```
Then collect `toFloat(neg) / total` for the negative fraction, plus `med`, `maxabs`, `large_ns`, `extreme`, `low_sig`, and SET the six arrays.

- [ ] **Step 1: Write the failing static test.** `test_post_import_sets_de_quality_arrays`: for both script files, assert that each of the six property names, plus `data_caveats = coalesce(e.data_caveats, [])`, occurs in the text. Then assert that the set of lines containing `time_point_` is identical in the two files.
- [ ] **Step 2: Run it and confirm it fails.** `uv run pytest tests/test_post_import_de_quality.py -v`. Expected: FAIL.
- [ ] **Step 3: Edit both scripts** as above.
- [ ] **Step 4: Run it and confirm it passes.** `uv run pytest tests/test_post_import_de_quality.py -v`. Expected: PASS.
- [ ] **Step 5: Commit.**
```bash
git add scripts/post-import.sh scripts/post-import.cypher tests/test_post_import_de_quality.py
git commit -m "feat(post-import): per-time-point table statistics + low-signal rollup + dense data_caveats"
```

### Task 6: Validator rules for caveats, the allowlist, level columns and scale

**Files:**
- Create: `config/data_caveats_unreviewed.yaml`
- Modify: `scripts/validate_paperconfig.py`. Add `validate_de_quality(pub, supp, config_path, allowlist) -> tuple[list, list]` next to `validate_control_timing` (line 293), and call it beside that call (line 1797).
- Test: `tests/test_validate_paperconfig_de_quality.py`

**Interfaces:**
- Consumes: `level_basis` (Task 2), the `Experiment.data_caveats` vocabulary (Task 3).
- The allowlist format, a mapping from repo-relative paperconfig path to unreviewed experiment keys:
```yaml
# DE experiments not yet reviewed for data_caveats (spec §7.4). Only shrinks:
# a review removes its keys. Delete this file and its check when empty.
"data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml":
  - coculture_alteromonas_hot1a3_med4_rnaseq
```
- `load_caveat_allowlist(path: Path | None = None) -> dict[str, set[str]]`. Keys are normalised with `Path(...).as_posix()`, relative to the project root.
- Rules, where DE experiment = has ≥ 1 `statistical_analyses` entry (the `validate_control_timing` definition):
  - errors:
    1. A `data_caveats` value not in the vocabulary.
    2. `paralogs_share_one_array_probe` or any value outside the vocabulary.
    3. A non-empty `data_caveats` without `data_caveat_note`.
    4. `data_caveat_note` containing `|` or `'`.
    5. `data_caveats_reviewed_on` not matching `^\d{4}-\d{2}-\d{2}$`, or not a valid date.
    6. A DE experiment **not** in the allowlist that lacks the `data_caveats` key or `data_caveats_reviewed_on`.
    7. A DE experiment **in** the allowlist that has `data_caveats_reviewed_on` (stale entry: "remove it from config/data_caveats_unreviewed.yaml").
    8. `level_basis` raising `ValueError`, with its message.
    9. A declared level column absent from the table's columns. Read the header with the existing `_read_csv_safe`, or the resolved CSV when present.
    10. Analyses of one experiment with different `fold_change_type`, or different non-`None` level bases.
  - Warnings:
    - `data_caveats` on a non-DE experiment.
    - An allowlist key naming an experiment that no longer exists.

- [ ] **Step 0: Record today's corpus result.** `uv run python scripts/validate_paperconfig.py --all | tail -1 > "$SCRATCH/validate_before.txt"`, where `$SCRATCH` is the session scratchpad.
- [ ] **Step 1: Write the failing tests.** Build configs in memory, the way `tests/test_validate_paperconfig_fold_change_type.py` does, and call `validate_de_quality(pub, supp, "x.yaml", allowlist)` directly. One test per rule:
  - `test_unknown_caveat_value`
  - `test_probe_caveat_rejected`
  - `test_caveats_need_note`
  - `test_note_forbids_pipe_and_quote`
  - `test_reviewed_on_must_be_iso_date` (`"2026-13-01"` errors)
  - `test_new_de_experiment_requires_review_fields`
  - `test_allowlisted_experiment_may_omit_review_fields`
  - `test_stale_allowlist_entry_errors`
  - `test_one_arm_level_cols_errors`
  - `test_missing_level_column_errors`
  - `test_mixed_fold_change_type_errors`
  - `test_caveats_on_non_de_experiment_warns`

  Each asserts that one error or warning substring is present and no unrelated errors are.
- [ ] **Step 2: Run them and confirm they fail.** `uv run pytest tests/test_validate_paperconfig_de_quality.py -v`. Expected: FAIL (ImportError).
- [ ] **Step 3: Implement `validate_de_quality`, `load_caveat_allowlist` and the call site.**
- [ ] **Step 4: Generate the initial allowlist.** List every DE experiment in every registered paperconfig that has no `data_caveats_reviewed_on`. Use a one-off script in the session scratchpad that calls `load_all_paperconfigs()` and writes the YAML above, sorted by path then key. Don't commit the script.
- [ ] **Step 5: Run the tests and the whole corpus.** `uv run pytest tests/test_validate_paperconfig_de_quality.py tests/test_paperconfig_validation.py -v && uv run python scripts/validate_paperconfig.py --all | tail -3`. Expected: the tests PASS, and the last line equals `validate_before.txt` from Step 0: the allowlist makes the new rules pass on every existing experiment.
- [ ] **Step 6: Commit.**
```bash
git add scripts/validate_paperconfig.py config/data_caveats_unreviewed.yaml tests/test_validate_paperconfig_de_quality.py
git commit -m "feat(validator): data_caveats rules, unreviewed allowlist, level columns, scale consistency"
```

### Task 7: Validator QA report and the `/paperconfig` QA step

**Files:**
- Modify: `scripts/validate_paperconfig.py`. Add `de_quality_report(config: dict, config_path: str) -> str`. Print it from `validate()` only in single-file mode, after `_print_results`; `--all` stays quiet. Pass a `report: bool` argument through `validate`.
- Modify: `.claude/skills/paperconfig/SKILL.md`. Add a "QA step" subsection after "Validation Script" (line 1285), and a checklist line in "Validation Checklist" (line 1178).
- Modify: `.claude/skills/paperconfig/UPDATE_EXISTING_CONFIGS.md`. This is the audit guide for existing paperconfigs, and it becomes the backfill checklist.
- Modify: `scripts/build_modified_csv/README.md` (builder convention).
- Test: `tests/test_validate_paperconfig_de_quality.py`

**Interfaces:**
- Consumes: `read_de_table`, `call_row`, `expression_status`, `level_basis`, `low_signal_mask`, `low_signal_reason`, `time_point_stats` (Tasks 1–2).
- Report format: one block per DE analysis, headed `QA <analysis id> (<experiment key>, tp <timepoint or ->)`. It lists:
  - `rows`, `negative_fraction`, `median_log2fc`, `max_abs_log2fc`, `large_fc_not_significant`, `extreme_fc_count`;
  - `low_signal_count` (or `level columns: not declared`);
  - when not declared, `candidate level columns: [...]`: headers matching, case-insensitive, `basemean|rpkm|nsaf|cpm|count|mean|expression|abundance`, excluding the declared fold-change and p-value columns. Peptide and SD columns are excluded on purpose (spec §4.4).
  - `caveats to consider: [...]`, with these prompt rules (never stored):
    - `extreme_fc_count > 0` or `low_signal_count > 0` → `fold_change_inflated_at_low_counts`;
    - a full table (`table_scope == all_detected_genes`) with `negative_fraction > 0.75 or < 0.25`, or `abs(median) > 0.5` → `conditions_not_normalized_to_each_other`;
    - the table's filename ends `_modified.csv` → both `fold_change_recomputed_by_kg` and `fold_change_sign_reversed_by_kg`, "(check the builder)";
    - `table_scope_detail` holds an integer N and the table's significant count differs from N by more than 20% → `significance_calls_differ_from_paper`.

  It ends with `This report is information for the curator; it never fails validation.`

- [ ] **Step 1: Write the failing tests:**
  - `test_report_lists_statistics_and_prompts`: a `tmp_path` CSV with a −11 log2FC row and no level columns. The report contains `extreme_fc_count: 1`, `fold_change_inflated_at_low_counts` and `candidate level columns: ['baseMean']` (the CSV has a `baseMean` column).
  - `test_report_never_adds_errors`: `validate_paperconfig_content` returns the same errors with and without the report.
  - `test_report_uses_source_csv_without_resolved`: no `_resolved.csv` present → the report is still produced.
- [ ] **Step 2: Run them and confirm they fail.** `uv run pytest tests/test_validate_paperconfig_de_quality.py -k report -v`. Expected: FAIL.
- [ ] **Step 3: Implement `de_quality_report`** and the single-file print.
- [ ] **Step 4: Write the SKILL.md QA step.** Four numbered steps, from spec §7.7 minus `array_design`:
  1. Propose level columns from the report's candidates and confirm them against the paper's Methods.
  2. Run the validator and read the QA report.
  3. Go through the caveat vocabulary with the paper open, starting from the prompts. Set `data_caveats` (`[]` allowed), a `data_caveat_note` when non-empty, and `data_caveats_reviewed_on`.
  4. Don't add new experiments to `config/data_caveats_unreviewed.yaml`.

  Add the new paperconfig fields to the experiment and statistical-analysis field tables in the same file, with one YAML example.

  In `UPDATE_EXISTING_CONFIGS.md`:
  - add a section "Reviewing data caveats (allowlist backfill)". It gives the same four steps, plus: remove the experiment's key from `config/data_caveats_unreviewed.yaml`, and delete the file and its validator check when it is empty;
  - add one line to its significance-order section saying that the calls are implemented in `multiomics_kg/utils/de_calls.py`, shared by the adapter and the validator.

  In `scripts/build_modified_csv/README.md`'s Convention: a builder that recomputes fold changes adds `fold_change_recomputed_by_kg`, and one that negates them adds `fold_change_sign_reversed_by_kg`, to every experiment it feeds, with the reason in `data_caveat_note`.
- [ ] **Step 5: Run the tests and a real paper.** `uv run pytest tests/test_validate_paperconfig_de_quality.py -v && uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml" | grep -A9 "^QA "`. Expected: PASS. The Aharonovich MED4 block shows `extreme_fc_count` > 0 and `fold_change_inflated_at_low_counts` among the prompts.
- [ ] **Step 6: Commit.**
```bash
git add scripts/validate_paperconfig.py .claude/skills/paperconfig/SKILL.md .claude/skills/paperconfig/UPDATE_EXISTING_CONFIGS.md scripts/build_modified_csv/README.md tests/test_validate_paperconfig_de_quality.py
git commit -m "feat(validator): DE QA report from shared de_calls; /paperconfig QA step"
```

### Task 8: Curate the named cases

**Files:**
- Modify: the paperconfigs below and `config/data_caveats_unreviewed.yaml` (remove each curated key).

Each experiment gets `data_caveats`, a `data_caveat_note` naming the evidence (and the affected time points or tables when partial) and `data_caveats_reviewed_on: "<today>"`. Before tagging, **check each assignment against the source**: the PDF text (`pdftotext`), the table, and for the `_kg` caveats the builder in `scripts/build_modified_csv/`. If the check doesn't confirm it, leave the caveat off and note why in the commit message. The starting assignments come from spec §4 item 1:

| Paperconfig | Experiment keys | Caveats |
|---|---|---|
| `Aharonovich 2016` | `coculture_alteromonas_hot1a3_med4_rnaseq` | `fold_change_inflated_at_low_counts`, `conditions_not_normalized_to_each_other`, `comparison_direction_unclear` |
| `Aharonovich 2016` | `coculture_alteromonas_hot1a3_mit9313_rnaseq_inoc_5e6`, `…_inoc_05e6` | `fold_change_inflated_at_low_counts` |
| `Aharonovich 2016` | `coculture_prochlorococcus_mit9313_hot1a3_rnaseq`, `coculture_prochlorococcus_mit9313_high_alteromonas_hot1a3_rnaseq` | `fold_change_recomputed_by_kg` |
| `Capovilla 2023` | `chitosan_addition_mit9303_rnaseq` | `fold_change_inflated_at_low_counts` |
| `Anjur 2026` | `growth_state_freeliving_planktonic_enrichment_mit9301_rnaseq` | `fold_change_inflated_at_low_counts` |
| `Huang 2020` dir (spec "Huang 2021") | `viral_s_sbp1_infection_wh7803_rnaseq` | `fold_change_inflated_at_low_counts`, `no_biological_replicates`, `significance_calls_differ_from_paper` |
| `barreto 2022` | the three experiments whose tables are S12–S14 (the `*_at_800_ez55_rnaseq` set; confirm from `supplementary_materials`) | `table_may_be_a_different_comparison` |
| `Pandhal 2007` | `light_HL_vs_LL_med4_proteomics` (check ML vs LL too) | `significance_calls_differ_from_paper` |
| `Oleza 2017` | the experiment fed by Table S7 | `significance_calls_differ_from_paper` |
| `Beliaev 2014` | the experiments whose tables the Beliaev builder recomputes | `fold_change_recomputed_by_kg` |
| `Fuszard 2012` | `pi_limitation_natl2a_itraq` | `fold_change_recomputed_by_kg` |
| `Ferrer 2020` | `coculture_tpseudonana_dss3_rnaseq` | `fold_change_sign_reversed_by_kg` |
| `Barber 2026` | the three experiments whose builder negates (`build_barber2026_rpom_proteomics.py`) | `fold_change_sign_reversed_by_kg` |
| `Koch 2019` | the DE experiments whose builder negates (`build_koch2019_de_nsaf.py`) | `fold_change_sign_reversed_by_kg` |

- [ ] **Step 1: Check and edit, one paper at a time.** After each paper, run `uv run python scripts/validate_paperconfig.py "<path>"`. Expected: `VALIDATION PASSED`.
- [ ] **Step 2: Run the whole corpus.** `uv run python scripts/validate_paperconfig.py --all | tail -1`. Expected: the same line as Task 6 Step 0.
- [ ] **Step 3: Commit,** with every edited paperconfig path listed explicitly, plus `config/data_caveats_unreviewed.yaml`. Message: `data(caveats): curate named DE limitations (spec §4.1)`. The body lists any assignment the source check did not confirm.

### Task 9: KG tests, docs, rebuild and verification

**Files:**
- Create: `tests/kg_validity/test_de_quality.py`
- Create: `docs/kg-changes/de-table-quality.md`
- Modify: `CLAUDE.md` (Experiment bullet in "Key graph facts"; paperconfig fields in "Adding Omics Data"), `CHANGELOG.md` (`## [Unreleased]` → `### Added` + `### Data`), `plans/backlog.md` (nothing to remove; check the follow-ups section still matches)
- Modify: `.claude/skills/cypher-queries/SKILL.md`. Add the new properties to the `Experiment` row (line 45) and the `Changes_expression_of` row (line 98) of its property tables, plus one ready-made query: "experiments with caveats, with their notes" (`MATCH (e:Experiment) WHERE size(e.data_caveats) > 0 RETURN e.id, e.data_caveats, e.data_caveat_note`).
- Modify: `docs/methods_omics_adapter.md`:
  - the significance cascade (line 108) now names `de_calls` as its single implementation;
  - add a "Table quality" subsection covering `fold_change_type`, the two low-signal rules with their thresholds, and the six statistics.
- Modify: `docs/kg_mcp_guide.md` (time-course arrays bullet, line 182). One bullet for alpha testers: read `data_caveats` (and `data_caveats_reviewed_on`) before trusting an experiment; the `time_point_*` statistics describe the published subset when `table_scope` is not `all_detected_genes`.
- Modify: `docs/superpowers/specs/2026-11-18-de-table-quality-design.md` line 114. Drop "`shares_probe_with` and" from the tied-groups note: the feature was removed in §10.

**Interfaces:**
- Consumes: `time_point_stats`, `expression_status` (Task 2); the `run_query` fixture from `tests/kg_validity/conftest.py`.

- [ ] **Step 1: Write the KG tests** (marker `kg`):
  - `test_de_experiments_have_dense_caveats_and_scale`: every Experiment with ≥ 1 `Changes_expression_of` edge has `data_caveats IS NOT NULL` and `fold_change_type IN ['log2', 'linear', 'signed_linear']`.
  - `test_statistics_arrays_parallel`: for every such experiment, `size()` of each of the six arrays equals `size(time_point_orders)`.
  - `test_statistics_match_python`: for these sample experiments, ordered by `time_point_order`, recompute `time_point_stats` from the live edges (`log2_fold_change`, `expression_status`, `low_signal_reason`) and compare with the stored arrays, using `pytest.approx` with `abs=1e-9`:
    - `10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq`
    - one multi-time-point experiment, the first returned by `MATCH (e:Experiment) WHERE e.time_point_count > 3 RETURN e.id ORDER BY e.id LIMIT 1`
    - one significant-only experiment (`table_scope = 'significant_only'`, first by id)
  - `test_curated_caveats_live`: `10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq` has `'fold_change_inflated_at_low_counts' IN data_caveats` and a `data_caveats_reviewed_on`.
  - `test_low_signal_count_is_edge_rollup`: for every experiment, `sum(time_point_low_signal_count)` equals the count of its edges with `low_signal_reason IS NOT NULL`. That count is 0 everywhere until level columns are declared.

  Coverage of the new properties by `test_schema_coverage.py` and `test_controlled_vocabularies.py` is automatic.
- [ ] **Step 2: Write `docs/kg-changes/de-table-quality.md`.** It covers: what changed, the property list, how to read each statistic by `table_scope`, the allowlist process, and what is deferred (link the backlog section). Add the CHANGELOG bullets:
  - `### Added`: the new properties and the validator QA report.
  - `### Data`: the curated experiments.

  Add the CLAUDE.md lines: one Experiment-bullet sentence, plus the new paperconfig fields in the paperconfig YAML example. Make the `cypher-queries`, `methods_omics_adapter.md` and `kg_mcp_guide.md` edits listed under Files.

  Skills and docs checked and left unchanged: `omics-edge-snapshot` (counts edges only), `check-gene-ids` / `fix-gene-ids` (gene IDs only), `release-kg` (the preflight reads CHANGELOG, which is updated here). The historical `docs/kg-changes/*` files are records and are not rewritten.
- [ ] **Step 3: Run the unit suite.** `uv run pytest -m "not slow and not kg" -q`. Expected: no failures.
- [ ] **Step 4: Commit,** then push for the Docker rebuild.
```bash
git add tests/kg_validity/test_de_quality.py docs/kg-changes/de-table-quality.md CLAUDE.md CHANGELOG.md docs/superpowers/specs/2026-11-18-de-table-quality-design.md .claude/skills/cypher-queries/SKILL.md docs/methods_omics_adapter.md docs/kg_mcp_guide.md
git commit -m "test+docs: KG checks and docs for DE table quality"
```
  The Docker build runs in the other clone (`multiomics_biocypher_kg`). Ask the user to pull and rebuild there. Before the rebuild, take the `/omics-edge-snapshot` baseline against the current graph.
- [ ] **Step 5: Verify after the rebuild** (once the user has rebuilt):
  - `/omics-edge-snapshot` compare. Expected: every paper's edge count is unchanged.
  - `uv run pytest tests/kg_validity/ -v`. Expected: PASS, including `test_de_quality.py`, `test_schema_coverage.py` and `test_controlled_vocabularies.py`.
  - Spot check: `MATCH (e:Experiment {id:'10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq'}) RETURN e.data_caveats, e.time_point_extreme_fc_count, e.time_point_large_fc_not_significant`. Expected: the three caveats, a count > 0, and 291.
