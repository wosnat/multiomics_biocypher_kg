# multiomics_kg/extraction/cluster/table.py
"""Path: table — CSV parsing (generic) + per-paper supplementary parsers.

Reads column names from the paperconfig. No hardcoded column names.
Two data sources:
1. Cluster CSV (always available) — membership + structured columns
2. Supplementary XLS/TXT (per-paper) — enrichment p-values
"""
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def extract_from_csv(table_config: dict,
                     project_root: Path = Path("."),
                     ) -> dict[str, dict]:
    """Extract per-cluster data from the cluster CSV.

    Reads gene_id_col, cluster_col from table_config.
    Auto-detects additional numeric columns and computes per-cluster aggregates.

    Returns {cluster_key: {genes: [...], gene_count: int, col_mean: float, ...}}
    """
    csv_path = Path(table_config["filename"])
    if not csv_path.is_absolute():
        csv_path = project_root / csv_path

    gene_id_col = table_config["gene_id_col"]
    cluster_col = table_config["cluster_col"]

    df = pd.read_csv(csv_path)

    # Columns to skip when auto-detecting numeric aggregates
    skip_cols = {gene_id_col, cluster_col}
    # Also skip any extra id/name columns
    for col in df.columns:
        if col.lower() in ("locus_tag", "gene_symbol", "gene_product",
                           "name", "desc", "bvbrc_id", "start", "stop",
                           "strand", "core_flexible"):
            skip_cols.add(col)

    numeric_cols = [
        c for c in df.columns
        if c not in skip_cols and pd.api.types.is_numeric_dtype(df[c])
    ]

    results: dict[str, dict] = {}
    for key, subset in df.groupby(cluster_col):
        # Convert integer-valued floats (e.g. 5.0 → "5") for consistent keys
        if isinstance(key, float) and key == int(key):
            str_key = str(int(key))
        else:
            str_key = str(key)
        entry: dict = {
            "genes": subset[gene_id_col].dropna().tolist(),
            "gene_count": len(subset),
            "source": "table",
            "confidence": "very_high",
        }

        for col in numeric_cols:
            vals = subset[col].dropna()
            if len(vals) == 0:
                continue
            entry[f"{col}_mean"] = round(float(vals.mean()), 4)
            entry[f"{col}_median"] = round(float(vals.median()), 4)

        results[str_key] = entry

    return results


def parse_tolonen_enrichment_xls(paper_dir: Path,
                                 cluster_keys: list[str],
                                 organism_hint: str = "",
                                 ) -> dict[str, dict]:
    """Tolonen 2006 specific: parse medPvals.xls / mitPvals.xls.

    Returns {cluster_key: {enrichment_category, enrichment_pvalue,
                           enrichment_significant, all_enrichments}}.
    """
    prefix = ""
    if "med4" in organism_hint.lower():
        prefix = "med"
    elif "mit9313" in organism_hint.lower():
        prefix = "mit"

    results: dict[str, dict] = {}
    for xls_path in sorted(paper_dir.rglob("*[Pp]val*.[Xx][Ll][Ss]*")):
        if not xls_path.is_file():
            continue
        if prefix and prefix not in xls_path.name.lower():
            continue
        try:
            df = pd.read_excel(xls_path, header=None)
            _parse_tolonen_pval_df(df, cluster_keys, results)
        except Exception as e:
            logger.debug("Could not parse %s: %s", xls_path, e)
    return results


def _parse_tolonen_pval_df(df: pd.DataFrame,
                           cluster_keys: list[str],
                           results: dict[str, dict]) -> None:
    """Parse Tolonen-format p-value table (rows=categories, cols=cluster numbers)."""
    header_idx = None
    for i in range(min(10, len(df))):
        cells = [str(v).strip() for v in df.iloc[i] if pd.notna(v)]
        has_fc = any("FUNCTIONAL CATEGORY" == c.upper() for c in cells)
        has_nums = any(c.replace(".0", "").isdigit() for c in cells)
        if has_fc and has_nums:
            header_idx = i
            break
    if header_idx is None:
        return

    header = df.iloc[header_idx]
    col_to_cluster: dict[int, str] = {}
    cat_col = 0
    for col_idx in range(len(header)):
        val = header.iloc[col_idx]
        if pd.isna(val):
            continue
        val_str = str(val).strip()
        if "FUNCTIONAL CATEGORY" in val_str.upper():
            cat_col = col_idx
            continue
        try:
            num = str(int(float(val_str)))
            if num in cluster_keys:
                col_to_cluster[col_idx] = num
        except (ValueError, TypeError):
            pass
    if not col_to_cluster:
        return

    per_cluster: dict[str, list[tuple[str, float, bool]]] = {
        k: [] for k in cluster_keys
    }
    for i in range(header_idx + 1, len(df)):
        category = df.iloc[i, cat_col]
        if pd.isna(category) or not str(category).strip():
            continue
        category = str(category).strip()
        for col_idx, ckey in col_to_cluster.items():
            val = df.iloc[i, col_idx]
            if pd.isna(val):
                continue
            val_str = str(val).strip()
            has_stars = "*" in val_str
            try:
                pval = float(val_str.replace("*", "").strip())
            except (ValueError, TypeError):
                continue
            per_cluster[ckey].append((category, pval, pval < 0.05 or has_stars))

    for ckey, entries in per_cluster.items():
        if not entries:
            continue
        entries.sort(key=lambda x: x[1])
        best_cat, best_pval, best_sig = entries[0]
        results[ckey] = {
            "enrichment_category": best_cat,
            "enrichment_pvalue": best_pval,
            "enrichment_significant": best_sig,
            "all_enrichments": [
                {"category": c, "pvalue": p, "significant": s}
                for c, p, s in entries if s
            ],
            "source": "table",
            "confidence": "very_high",
        }


def run_table(table_config: dict,
              paper_dir: Path,
              cluster_keys: list[str],
              organism_hint: str = "",
              project_root: Path = Path("."),
              ) -> dict[str, dict]:
    """Run Path: table — CSV extraction + supplementary parsers.

    Returns per-cluster dict with all fields extractable from data files.
    """
    results = extract_from_csv(table_config, project_root)

    # Paper-specific enrichment parsers
    enrichments = parse_tolonen_enrichment_xls(paper_dir, cluster_keys,
                                               organism_hint)
    for key, enr in enrichments.items():
        if key in results:
            results[key].update(enr)
        else:
            results[key] = enr

    return results
