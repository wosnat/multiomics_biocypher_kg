"""Per-row differential-expression calls, shared by the omics adapter and the validator.

Pure functions: read a DE table, turn one row into (log2 fold change, p-value,
significance). Moved out of ``omics_adapter`` unchanged so every consumer judges
a row the same way.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import pandas as pd

from multiomics_kg.download.resolve_paper_ids import get_resolved_path


def read_de_table(filename: str, analysis: dict) -> tuple[pd.DataFrame, bool]:
    """Read a DE table, preferring the pre-resolved CSV. Returns ``(df, used_resolved)``.

    A ``_resolved.csv`` is always a clean comma-separated file with no extra header
    rows, so it is read with plain ``pd.read_csv`` whatever ``sep`` / ``skip_rows`` say.
    """
    resolved_path = get_resolved_path(filename)
    if resolved_path.exists():
        return pd.read_csv(str(resolved_path)), True
    skip_rows = analysis.get('skip_rows', 0)
    sep = analysis.get('sep', ',')
    if skip_rows:
        return pd.read_csv(filename, sep=sep, skiprows=skip_rows), False
    return pd.read_csv(filename, sep=sep), False


def convert_fold_change(fc: float, fold_change_type: str | None) -> float | None:
    """Convert fold-change value to log2 scale.

    Args:
        fc: The raw fold-change value from the CSV.
        fold_change_type: 'log2' (default), 'linear', 'signed_linear', or None
            (treated as log2). 'signed_linear' is the "+3 = 3-fold up, -3 = 3-fold
            down" convention: log2 = sign(x) * log2|x|.

    Returns:
        log2-scale fold-change, or None if conversion is invalid (linear FC <= 0;
        signed linear FC in (-1, 1)).
    """
    if fold_change_type is None or fold_change_type == "log2":
        return fc
    if fold_change_type == "linear":
        if fc <= 0:
            return None
        return math.log2(fc)
    if fold_change_type == "signed_linear":
        # Values in (-1, 1) cannot occur in this convention; the whole analysis is
        # refused up front (signed_linear_violations), so this is only a guard.
        if abs(fc) < 1:
            return None
        return math.copysign(math.log2(abs(fc)), fc)
    return fc


def check_significance(
    log2fc: float | None,
    pvalue: float | None,
    asterisk_significant: bool | None,
    analysis: dict,
    *,
    default_pvalue_threshold: float | None = None,
    default_logfc_threshold: float | None = None,
    p_missing: bool = False,
) -> str:
    """Determine whether a row is statistically significant.

    Priority order:
    1. prefiltered: true -> always "significant"
    2. pvalue_asterisk_in_logfc -> use asterisk_significant
    3. p_missing (the table has a p-value column but this row's cell is empty,
       e.g. a DESeq2 NA padj) -> "not significant": untested is not significant
    4. Thresholds (config overrides the defaults)
    5. "unknown" if no info available

    Returns:
        "significant", "not significant", or "unknown".
    """
    if analysis.get('prefiltered'):
        return 'significant'

    if asterisk_significant is not None:
        return 'significant' if asterisk_significant else 'not significant'

    if p_missing:
        return 'not significant'

    # `is not None`, not `or`: an explicit threshold of 0 must not fall back to the default.
    pval_thresh = analysis.get('pvalue_threshold')
    if pval_thresh is None:
        pval_thresh = default_pvalue_threshold
    logfc_thresh = analysis.get('logfc_threshold')
    if logfc_thresh is None:
        logfc_thresh = default_logfc_threshold

    if pval_thresh is None and logfc_thresh is None:
        return 'unknown'

    sig = True
    if logfc_thresh is not None and log2fc is not None:
        sig = sig and (abs(log2fc) >= logfc_thresh)
    if pval_thresh is not None and pvalue is not None:
        sig = sig and (pvalue <= pval_thresh)

    return 'significant' if sig else 'not significant'


@dataclass(frozen=True)
class RowCall:
    raw_fc: float
    log2_fold_change: float | None  # None when the conversion is invalid
    adjusted_p_value: float | None
    significant: str


def call_row(
    row: dict,
    analysis: dict,
    columns: set[str],
    *,
    default_pvalue_threshold: float | None = None,
    default_logfc_threshold: float | None = None,
) -> RowCall | None:
    """Fold change, p-value and significance for one table row.

    Returns None when the fold-change cell is empty, ``NA`` or non-numeric (the row
    carries no call). A RowCall with ``log2_fold_change=None`` means the value could
    not be converted to log2 (e.g. linear FC <= 0).
    """
    logfc_col = analysis.get('logfc_col')
    p_value_col = analysis.get('adjusted_p_value_col')
    pvalue_asterisk = analysis.get('pvalue_asterisk_in_logfc', False)
    asterisk_min_level = int(analysis.get('asterisk_min_level', 1))

    fc_val = row.get(logfc_col)
    if pd.isna(fc_val):
        return None
    fc_str = str(fc_val).strip()
    if fc_str == '' or fc_str == 'NA':
        return None
    asterisk_significant = None
    # Handle asterisk-based significance markers in logFC column
    # Supports trailing "1.1 *", leading "* 1.1", or bare "1.1*"
    if pvalue_asterisk:
        # Count the marker: "**1.1", "1.1 **". A paper whose "*" is a
        # weaker band than its DE criterion sets asterisk_min_level: 2.
        n_stars = max(len(fc_str) - len(fc_str.lstrip('*')),
                      len(fc_str) - len(fc_str.rstrip('*')))
        asterisk_significant = n_stars >= asterisk_min_level
    # Always strip stray whitespace and stars from logFC values
    fc_str = fc_str.strip('*').strip()
    try:
        raw_fc = float(fc_str)
    except (ValueError, TypeError):
        return None
    if not math.isfinite(raw_fc):
        return None

    log2fc = convert_fold_change(raw_fc, analysis.get('fold_change_type', None))

    pval = None
    if asterisk_significant is not None:
        # Asterisk encodes significance at pvalue_threshold; use threshold as placeholder
        pval_thresh = analysis.get('pvalue_threshold') or default_pvalue_threshold or 0.05
        pval = pval_thresh if asterisk_significant else 1.0
    elif p_value_col and p_value_col in columns and not pd.isna(row.get(p_value_col)):
        try:
            pval = float(row[p_value_col])
            if not math.isfinite(pval):
                pval = None
        except (ValueError, TypeError):
            pass

    # A declared p-value column with an empty / non-numeric cell (DESeq2 NA padj,
    # single-peptide proteins) means "not tested", never "significant on fold alone".
    p_missing = (asterisk_significant is None and pval is None
                 and bool(p_value_col) and p_value_col in columns)
    significant = check_significance(
        log2fc, pval, asterisk_significant, analysis,
        default_pvalue_threshold=default_pvalue_threshold,
        default_logfc_threshold=default_logfc_threshold,
        p_missing=p_missing,
    )
    return RowCall(raw_fc, log2fc, pval, significant)


def expression_status(significant: str, log2fc: float | None) -> str:
    """Mirror of the post-import ``expression_status``."""
    if significant == 'significant' and log2fc is not None:
        if log2fc > 0:
            return 'significant_up'
        if log2fc < 0:
            return 'significant_down'
    return 'not_significant'
