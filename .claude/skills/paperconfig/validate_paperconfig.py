#!/usr/bin/env python
"""Validate a paperconfig.yaml file for the multiomics knowledge graph.

Usage:
    python .claude/skills/paperconfig/validate_paperconfig.py <path_to_paperconfig.yaml>

Checks:
    - YAML is parseable
    - Required top-level fields exist
    - PDF file exists
    - CSV files exist
    - Column names (name_col, logfc_col, adjusted_p_value_col) match CSV headers
    - Environmental condition references resolve
    - Each analysis has either organism or environmental condition edge source
    - All required analysis fields are present
    - Analysis IDs are unique
    - skip_rows is applied correctly when reading CSV headers
    - Sample data rows are readable and logfc_col contains numeric-like values
"""

import sys
import os
import yaml
import pandas as pd

REQUIRED_ANALYSIS_FIELDS = [
    "type", "name", "id", "test_type",
    "control_condition", "treatment_condition",
    "organism", "name_col", "logfc_col",
]

VALID_TYPES = {"RNASEQ", "MICROARRAY", "PROTEOMICS", "METABOLOMICS"}


def validate(config_path: str) -> bool:
    errors = []
    warnings = []

    # --- Parse YAML ---
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"FAIL: Cannot parse YAML: {e}")
        return False

    if not config or "publication" not in config:
        print("FAIL: Missing top-level 'publication' key")
        return False

    pub = config["publication"]

    # --- papername ---
    if "papername" not in pub:
        errors.append("Missing 'papername'")
    else:
        print(f"  papername: {pub['papername']}")

    # --- PDF ---
    pdf = pub.get("papermainpdf", "")
    if not pdf:
        warnings.append("Missing 'papermainpdf' (optional but recommended)")
    elif not os.path.exists(pdf):
        errors.append(f"PDF not found: {pdf}")
    else:
        print(f"  PDF exists: {pdf}")

    # --- Environmental conditions ---
    env_conds = pub.get("environmental_conditions", {})
    if env_conds:
        print(f"  Environmental conditions: {list(env_conds.keys())}")

    # --- Supplementary materials ---
    supp = pub.get("supplementary_materials")
    if not supp:
        errors.append("Missing 'supplementary_materials'")
        _print_results(errors, warnings)
        return len(errors) == 0

    all_ids = []

    for table_key, table in supp.items():
        print(f"\n  [{table_key}]")
        fn = table.get("filename", "")
        if not fn:
            errors.append(f"{table_key}: Missing 'filename'")
            continue

        if not os.path.exists(fn):
            errors.append(f"{table_key}: File not found: {fn}")
            continue

        print(f"    file: {fn}")

        # Read CSV headers
        skip = table.get("skip_rows", 0)
        try:
            if skip:
                df = pd.read_csv(fn, skiprows=range(1, skip + 1), nrows=5)
            else:
                df = pd.read_csv(fn, nrows=5)
            cols = list(df.columns)
        except Exception as e:
            errors.append(f"{table_key}: Cannot read CSV: {e}")
            continue

        print(f"    columns: {cols}")

        analyses = table.get("statistical_analyses", [])
        if not analyses:
            errors.append(f"{table_key}: No 'statistical_analyses' defined")
            continue

        for i, analysis in enumerate(analyses):
            aid = analysis.get("id", f"<missing-id-{i}>")
            all_ids.append(aid)
            print(f"    Analysis: {aid}")

            # Check required fields
            for field in REQUIRED_ANALYSIS_FIELDS:
                if field not in analysis:
                    errors.append(f"{aid}: Missing required field '{field}'")

            # Check type value
            atype = analysis.get("type", "")
            if atype and atype not in VALID_TYPES:
                warnings.append(f"{aid}: type '{atype}' not in {VALID_TYPES}")

            # Check column references
            name_col = analysis.get("name_col", "")
            logfc_col = analysis.get("logfc_col", "")
            pval_col = analysis.get("adjusted_p_value_col")

            if name_col and name_col not in cols:
                errors.append(f"{aid}: name_col '{name_col}' not found in CSV columns")
            elif name_col:
                print(f"      name_col '{name_col}': OK")

            if logfc_col and logfc_col not in cols:
                errors.append(f"{aid}: logfc_col '{logfc_col}' not found in CSV columns")
            elif logfc_col:
                print(f"      logfc_col '{logfc_col}': OK")
                # Check if values look numeric
                if len(df) > 0:
                    sample_vals = df[logfc_col].dropna().head(3).tolist()
                    non_numeric = []
                    for v in sample_vals:
                        sv = str(v).rstrip("*")
                        try:
                            float(sv)
                        except ValueError:
                            non_numeric.append(v)
                    if non_numeric:
                        warnings.append(
                            f"{aid}: logfc_col '{logfc_col}' has non-numeric values: {non_numeric}"
                        )

            if pval_col:
                if pval_col not in cols:
                    errors.append(
                        f"{aid}: adjusted_p_value_col '{pval_col}' not found in CSV columns"
                    )
                else:
                    print(f"      adjusted_p_value_col '{pval_col}': OK")

            # Check edge source: organism or environmental condition
            has_treatment_org = "treatment_organism" in analysis and "treatment_taxid" in analysis
            has_env_treat = "environmental_treatment_condition_id" in analysis

            if not has_treatment_org and not has_env_treat:
                errors.append(
                    f"{aid}: Must have either (treatment_organism + treatment_taxid) "
                    "or environmental_treatment_condition_id"
                )

            # Validate environmental condition references
            env_treat_id = analysis.get("environmental_treatment_condition_id")
            env_ctrl_id = analysis.get("environmental_control_condition_id")

            if env_treat_id:
                if env_treat_id in env_conds:
                    print(f"      env_treatment_id '{env_treat_id}': OK")
                else:
                    errors.append(
                        f"{aid}: environmental_treatment_condition_id '{env_treat_id}' "
                        "not found in environmental_conditions"
                    )

            if env_ctrl_id:
                if env_ctrl_id in env_conds:
                    print(f"      env_control_id '{env_ctrl_id}': OK")
                else:
                    errors.append(
                        f"{aid}: environmental_control_condition_id '{env_ctrl_id}' "
                        "not found in environmental_conditions"
                    )

    # --- Check ID uniqueness ---
    print(f"\n  All analysis IDs: {all_ids}")
    seen = set()
    dupes = []
    for aid in all_ids:
        if aid in seen:
            dupes.append(aid)
        seen.add(aid)
    if dupes:
        errors.append(f"Duplicate analysis IDs: {dupes}")
    else:
        print("  All IDs unique: OK")

    _print_results(errors, warnings)
    return len(errors) == 0


def _print_results(errors, warnings):
    print()
    if warnings:
        print(f"WARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  ⚠ {w}")
    if errors:
        print(f"ERRORS ({len(errors)}):")
        for e in errors:
            print(f"  ✗ {e}")
        print("\nVALIDATION FAILED")
    else:
        print("VALIDATION PASSED")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_paperconfig.yaml>")
        sys.exit(2)

    config_path = sys.argv[1]
    if not os.path.exists(config_path):
        print(f"File not found: {config_path}")
        sys.exit(1)

    print(f"Validating: {config_path}\n")
    ok = validate(config_path)
    sys.exit(0 if ok else 1)
