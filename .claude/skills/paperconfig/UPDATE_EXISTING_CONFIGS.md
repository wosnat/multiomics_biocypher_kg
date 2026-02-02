# Updating Existing paperconfig.yaml Files with Significance Fields

## New Optional Fields

Three new fields can be added to each `statistical_analyses` entry to describe how significance is determined in the publication data:

| Field | Type | Description |
|-------|------|-------------|
| `prefiltered` | boolean | `true` if the table only contains rows the authors deemed significant |
| `pvalue_threshold` | float (0, 1] | Adjusted p-value cutoff stated in the paper (e.g., `0.05`) |
| `logfc_threshold` | float >= 0 | Absolute log2 fold-change cutoff stated in the paper (e.g., `1.0`) |

These fields describe what the **publication** says about significance. They do not change existing behavior unless the adapter is constructed with `significance_mode="significant_only"`.

## How to Decide Which Fields to Add

For each statistical analysis entry, determine which category it falls into:

### Category 1: Pre-filtered tables
The supplementary table only lists genes that passed significance criteria. The authors already filtered out non-significant rows before publishing.

**How to identify:** The paper or table legend says something like "Table contains only significantly differentially expressed genes" or "Genes with adjusted p-value < 0.05 and |log2FC| > 1".

**Action:** Add `prefiltered: true`. Optionally also record the thresholds used for documentation:
```yaml
prefiltered: true
pvalue_threshold: 0.05    # optional, for documentation
logfc_threshold: 1.0      # optional, for documentation
```

### Category 2: Asterisk-marked significance (`pvalue_asterisk_in_logfc: true`)
The table contains all genes, with asterisks in the fold-change column marking significant ones.

**How to identify:** Already has `pvalue_asterisk_in_logfc: true` in the config.

**Action:** No changes required. The adapter already uses asterisks to determine significance. Optionally add thresholds if the paper states them:
```yaml
pvalue_asterisk_in_logfc: true
pvalue_threshold: 0.1     # if the paper states the asterisk threshold
```

### Category 3: Tables with p-value columns
The table includes a column with adjusted p-values. The adapter stores the p-value but currently does not filter.

**How to identify:** Already has `adjusted_p_value_col` in the config.

**Action:** Add `pvalue_threshold` and/or `logfc_threshold` if the paper states significance criteria:
```yaml
adjusted_p_value_col: padj
pvalue_threshold: 0.05
logfc_threshold: 1.0
```

### Category 4: No significance information
The table has fold-change values but no p-values, no asterisks, and is not stated as pre-filtered.

**How to identify:** No `adjusted_p_value_col`, no `pvalue_asterisk_in_logfc`, and the paper doesn't say the table is pre-filtered.

**Action:** Leave as-is. The adapter will not assign a `significant` property to edges from this analysis. Alternatively, if you believe the table is implicitly pre-filtered (common for microarray papers), add `prefiltered: true`.

## Existing Config Assessment

Below is each existing paperconfig with recommended updates. **Read the actual paper** to confirm the thresholds before adding them.

### Al-Hosani 2015 (1 analysis)
- Has `adjusted_p_value_col: padj`
- **Check paper for:** stated p-value and FC thresholds
- **Likely update:** Add `pvalue_threshold` and optionally `logfc_threshold`, or `prefiltered: true` if the table only lists DE genes

### Aharonovich 2016 (3 analyses)
- Has `adjusted_p_value_col: q Value`
- **Check paper for:** Rockhopper default significance criteria
- **Likely update:** Add `pvalue_threshold` (Rockhopper default is typically q < 0.01)

### Anjur-Dietrich 2025 (1 analysis)
- Has `adjusted_p_value_col: FDR`
- Already mentions thresholds in `experimental_context`: "log2FC threshold 0.8 and q-value threshold 0.05"
- **Update:** Add `pvalue_threshold: 0.05` and `logfc_threshold: 0.8`, and likely `prefiltered: true`

### Bagby and Chisholm 2015 (4 analyses)
- No p-value column, no asterisks
- Microarray data
- **Check paper for:** whether tables are pre-filtered
- **Likely update:** `prefiltered: true` if the paper states only significant genes are listed

### Barreto 2022 (9 analyses)
- Has `adjusted_p_value_col: PValue`
- **Check paper for:** edgeR significance thresholds
- **Likely update:** Add `pvalue_threshold` and `logfc_threshold`

### Biller 2016 (9 analyses)
- Has `pvalue_asterisk_in_logfc: true`
- **Check paper for:** what the asterisk threshold is (likely adjusted p < 0.1)
- **Optional update:** Add `pvalue_threshold: 0.1` for documentation

### Biller 2018 (8 analyses)
- Has `pvalue_asterisk_in_logfc: true`
- Same as Biller 2016
- **Optional update:** Add `pvalue_threshold` for documentation

### Coe 2024 (16 analyses)
- Has `pvalue_asterisk_in_logfc: true`
- **Check paper for:** asterisk threshold
- **Optional update:** Add `pvalue_threshold` for documentation

### Fang 2019 (9 analyses)
- Has `pvalue_asterisk_in_logfc: true`
- **Check paper for:** asterisk threshold
- **Optional update:** Add `pvalue_threshold` for documentation

### He 2022 (2 analyses)
- Has `adjusted_p_value_col: p-value`
- **Check paper for:** stated significance criteria
- **Likely update:** Add `pvalue_threshold` and optionally `logfc_threshold`

## Validation

After updating any paperconfig, run the validation script:

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/<Author Year>/paperconfig.yaml"
```

The validator checks that:
- `prefiltered` is a boolean
- `pvalue_threshold` is a number in (0, 1]
- `logfc_threshold` is a non-negative number
- `prefiltered: true` is not redundantly combined with thresholds (warning)
- `prefiltered: true` is not combined with `pvalue_asterisk_in_logfc: true` (warning)
