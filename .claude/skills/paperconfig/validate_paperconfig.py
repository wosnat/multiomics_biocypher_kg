#!/usr/bin/env python
"""Validate a paperconfig.yaml file for the multiomics knowledge graph.

Usage:
    python .claude/skills/paperconfig/validate_paperconfig.py <path_to_paperconfig.yaml>

Checks:
    - YAML is parseable
    - Required top-level fields exist (or strain-level resource config with no publication block)
    - PDF file exists
    - CSV files exist and are readable
    - Column names (name_col, logfc_col, adjusted_p_value_col) match CSV headers
    - id_columns column names match CSV headers for csv and id_translation entries
    - product_columns column names match CSV headers
    - annotation_gff entries have organism and a valid GFF/GTF/GTF extension
    - id_translation entries have organism and id_columns
    - Environmental condition references resolve
    - Each analysis has either organism or environmental condition edge source
    - All required analysis fields are present
    - Analysis IDs are unique
    - skip_rows is applied correctly when reading CSV headers
    - Sample data rows are readable and logfc_col contains numeric-like values
    - Canonical vocabulary: organism, condition_type, test_type must use allowed values
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

# ── New: Required statistical_analyses fields (enforced as errors) ──────────
# These three are required regardless of other fields:
REQUIRED_STATS_FIELDS = {
    "id": "unique identifier within publication",
    "type": "one of RNASEQ, PROTEOMICS, METABOLOMICS, MICROARRAY",
    "treatment_condition": "non-empty string describing the treatment condition",
}

VALID_TYPES = {"RNASEQ", "MICROARRAY", "PROTEOMICS", "METABOLOMICS"}
VALID_ID_TYPES = {
    "locus_tag", "locus_tag_ncbi", "locus_tag_cyanorak",
    "old_locus_tag", "alternative_locus_tag",
    "gene_name", "gene_synonym",
    "protein_id_refseq", "uniprot_accession", "uniprot_entry_name",
    "jgi_id", "probeset", "rast_id", "annotation_specific", "other",
}

# ── Canonical vocabulary ─────────────────────────────────────────────────────

# Canonical organism names for the 'organism' field in statistical_analyses
# and for the 'organism' field in supplementary table entries.
CANONICAL_GENOMIC_ORGANISMS = {
    "Prochlorococcus MED4",
    "Prochlorococcus AS9601",
    "Prochlorococcus MIT9301",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL1A",
    "Prochlorococcus NATL2A",
    "Prochlorococcus RSP50",
    "Synechococcus WH8102",
    "Synechococcus CC9311",
    "Alteromonas macleodii HOT1A3",
    "Alteromonas macleodii EZ55",
    "Alteromonas MIT1002",
}

# Canonical treatment organism names loaded from treatment_organisms.csv
# (genus-level organisms used as coculture partners with no loaded genome).
# These are loaded at runtime from the CSV; this set is the fallback default.
_TREATMENT_ORGANISMS_CSV_DEFAULT = {
    "Phage",
    "Marinobacter",
    "Thalassospira",
    "Pseudohoeflea",
    "Alteromonas",
}

# Canonical condition_type values for environmental_conditions entries.
CANONICAL_CONDITION_TYPES = {
    "growth_medium",
    "nitrogen_stress",
    "phosphorus_stress",
    "iron_stress",
    "salt_stress",
    "carbon_stress",
    "light_stress",
    "darkness",
    "plastic_stress",
    "viral",
    "coculture",
    "growth_state",
    "temperature_stress",
}

# Canonical test_type values for statistical_analyses entries.
CANONICAL_TEST_TYPES = {
    "DESeq2",
    "DESeq",
    "edgeR",
    "Rockhopper",
    "microarray",
    "microarray_Cyber-T",
    "microarray_LPE",
    "microarray_Goldenspike",
}


def _load_treatment_organisms(project_root: str) -> set:
    """Load treatment organism names from treatment_organisms.csv.

    Returns the default set if the file cannot be read.
    """
    csv_path = os.path.join(
        project_root,
        "data", "Prochlorococcus", "treatment_organisms.csv",
    )
    try:
        df = pd.read_csv(csv_path, comment="#")
        names = set(df["organism_name"].dropna().str.strip().tolist())
        return names
    except Exception:
        return set(_TREATMENT_ORGANISMS_CSV_DEFAULT)


def _canonical_organism_error(config_path: str, context: str, value: str, allowed: set) -> str:
    """Format a clear vocabulary error message for an organism field."""
    sorted_allowed = sorted(allowed)
    return (
        f"{config_path} | {context} | organism '{value}' is not in the canonical list | "
        f"allowed: {sorted_allowed}"
    )


def _canonical_field_error(config_path: str, context: str, field: str, value: str, allowed: set) -> str:
    """Format a clear vocabulary error message for a canonical enum field."""
    sorted_allowed = sorted(allowed)
    return (
        f"{config_path} | {context} | {field} '{value}' is not a canonical value | "
        f"allowed: {sorted_allowed}"
    )


def _read_csv_safe(fn, sep, skip, errors, table_key):
    """Read CSV headers + 5 rows. Returns (df, cols) or (None, None) on failure."""
    try:
        if skip:
            df = pd.read_csv(fn, sep=sep, skiprows=skip, nrows=5)
        else:
            df = pd.read_csv(fn, sep=sep, nrows=5)
        return df, list(df.columns)
    except Exception as e:
        errors.append(f"{table_key}: Cannot read CSV: {e}")
        return None, None


def _validate_id_columns(id_columns, cols, table_key, errors, warnings):
    """Check each declared id_column against actual CSV columns."""
    for id_col_def in id_columns:
        col_name = id_col_def.get("column", "")
        id_type = id_col_def.get("id_type", "")
        if not col_name:
            errors.append(f"{table_key}: id_columns entry missing 'column'")
        elif col_name not in cols:
            errors.append(f"{table_key}: id_columns column '{col_name}' not found in CSV")
        else:
            print(f"      id_col '{col_name}' ({id_type}): OK")
        if id_type and id_type not in VALID_ID_TYPES:
            warnings.append(
                f"{table_key}: id_type '{id_type}' not in known types; "
                f"valid: {sorted(VALID_ID_TYPES)}"
            )


def _validate_product_columns(product_columns, cols, table_key, warnings):
    """Check each declared product_column against actual CSV columns."""
    for prod_col_def in product_columns:
        col_name = prod_col_def.get("column", "")
        if not col_name:
            continue
        if col_name not in cols:
            warnings.append(f"{table_key}: product_columns column '{col_name}' not found in CSV")
        else:
            print(f"      product_col '{col_name}': OK")


def validate(config_path: str) -> bool:
    errors = []
    warnings = []

    # Derive project root from the config file path so we can load treatment_organisms.csv
    # regardless of where the script is invoked from.
    abs_config = os.path.abspath(config_path)
    # Walk up from the config file until we find a directory containing
    # data/Prochlorococcus/treatment_organisms.csv (i.e., the project root).
    _project_root = os.path.dirname(abs_config)
    for _ in range(10):
        candidate = os.path.join(
            _project_root, "data", "Prochlorococcus", "treatment_organisms.csv"
        )
        if os.path.exists(candidate):
            break
        _project_root = os.path.dirname(_project_root)
    else:
        _project_root = None

    treatment_organisms = _load_treatment_organisms(_project_root or "")
    # The full allowed set for the 'organism' field is the union of genomic + treatment organisms
    all_canonical_organisms = CANONICAL_GENOMIC_ORGANISMS | treatment_organisms

    # --- Parse YAML ---
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"FAIL: Cannot parse YAML: {e}")
        return False

    # Detect strain-level resource config (no 'publication' block)
    is_resource_config = bool(config) and "publication" not in config

    if is_resource_config:
        print("  [strain-level resource config — no publication block]")
        pub = {}
        supp = config.get("supplementary_materials")
        env_conds = {}
    else:
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
            # Validate condition_type in each environmental condition
            for cond_key, cond_val in env_conds.items():
                if not isinstance(cond_val, dict):
                    continue
                cond_type = cond_val.get("condition_type")
                if cond_type is not None and cond_type not in CANONICAL_CONDITION_TYPES:
                    errors.append(
                        _canonical_field_error(
                            config_path,
                            f"environmental_conditions.{cond_key}",
                            "condition_type",
                            cond_type,
                            CANONICAL_CONDITION_TYPES,
                        )
                    )
                elif cond_type is None:
                    warnings.append(
                        f"{config_path} | environmental_conditions.{cond_key} | "
                        "missing 'condition_type'"
                    )

        supp = pub.get("supplementary_materials")

    if not supp:
        errors.append("Missing 'supplementary_materials'")
        _print_results(errors, warnings)
        return len(errors) == 0

    all_ids = []

    for table_key, table in supp.items():
        table_type = table.get("type", "csv")
        print(f"\n  [{table_key}] (type: {table_type})")

        fn = table.get("filename", "")
        if not fn:
            errors.append(f"{table_key}: Missing 'filename'")
            continue

        if not os.path.exists(fn):
            errors.append(f"{table_key}: File not found: {fn}")
            continue

        print(f"    file: {fn}")

        # ── annotation_gff ──────────────────────────────────────────────────
        if table_type == "annotation_gff":
            organism = table.get("organism", "")
            if not organism:
                errors.append(f"{table_key}: annotation_gff requires 'organism'")
            else:
                print(f"    organism: {organism}")
                if organism not in all_canonical_organisms:
                    errors.append(
                        _canonical_organism_error(
                            config_path, table_key, organism, all_canonical_organisms
                        )
                    )
            ext = os.path.splitext(fn)[1].lower()
            if ext not in (".gff", ".gff3", ".gtf"):
                warnings.append(
                    f"{table_key}: expected .gff/.gff3/.gtf extension, got '{ext}'"
                )
            # No CSV parsing, no analyses needed for this type
            continue

        # ── id_translation ──────────────────────────────────────────────────
        if table_type == "id_translation":
            organism = table.get("organism", "")
            if not organism:
                errors.append(f"{table_key}: id_translation requires 'organism'")
            else:
                print(f"    organism: {organism}")
                if organism not in all_canonical_organisms:
                    errors.append(
                        _canonical_organism_error(
                            config_path, table_key, organism, all_canonical_organisms
                        )
                    )

            id_columns = table.get("id_columns", [])
            if not id_columns:
                errors.append(f"{table_key}: id_translation requires 'id_columns'")

            sep = table.get("sep", ",")
            skip = table.get("skip_rows", 0)
            df, cols = _read_csv_safe(fn, sep, skip, errors, table_key)
            if cols is None:
                continue
            print(f"    columns: {cols}")

            if id_columns:
                _validate_id_columns(id_columns, cols, table_key, errors, warnings)

            prod_cols = table.get("product_columns", [])
            _validate_product_columns(prod_cols, cols, table_key, warnings)

            # No statistical_analyses for id_translation
            continue

        # ── csv (default) ───────────────────────────────────────────────────
        sep = table.get("sep", ",")
        skip = table.get("skip_rows", 0)
        df, cols = _read_csv_safe(fn, sep, skip, errors, table_key)
        if cols is None:
            continue

        print(f"    columns: {cols}")

        # Validate id_columns if declared
        id_columns = table.get("id_columns", [])
        if id_columns:
            _validate_id_columns(id_columns, cols, table_key, errors, warnings)

        # Validate product_columns if declared
        prod_cols = table.get("product_columns", [])
        _validate_product_columns(prod_cols, cols, table_key, warnings)

        analyses = table.get("statistical_analyses", [])
        if not analyses:
            errors.append(f"{table_key}: type 'csv' requires 'statistical_analyses'")
            continue

        for i, analysis in enumerate(analyses):
            aid = analysis.get("id", f"<missing-id-{i}>")
            all_ids.append(aid)
            print(f"    Analysis: {aid}")

            # ── Required fields (canonical enforcement) ──────────────────────
            # Check the three strictly required fields first (id, type, treatment_condition)
            for req_field, req_desc in REQUIRED_STATS_FIELDS.items():
                val = analysis.get(req_field)
                if val is None or (isinstance(val, str) and not val.strip()):
                    errors.append(
                        f"{config_path} | {aid} | missing required field '{req_field}' "
                        f"({req_desc})"
                    )

            # Check all original required fields (broader set, some overlap)
            for field in REQUIRED_ANALYSIS_FIELDS:
                if field not in analysis:
                    errors.append(f"{aid}: Missing required field '{field}'")

            # Check type value (canonical enum + original set)
            atype = analysis.get("type", "")
            if atype and atype not in VALID_TYPES:
                errors.append(
                    _canonical_field_error(
                        config_path, aid, "type", atype, VALID_TYPES
                    )
                )

            # ── Canonical organism validation ────────────────────────────────
            organism = analysis.get("organism", "")
            if organism and organism not in all_canonical_organisms:
                errors.append(
                    _canonical_organism_error(
                        config_path, aid, organism, all_canonical_organisms
                    )
                )
            elif organism:
                print(f"      organism '{organism}': OK")

            # ── Canonical treatment_organism validation ───────────────────────
            treatment_organism = analysis.get("treatment_organism", "")
            if treatment_organism and treatment_organism not in all_canonical_organisms:
                errors.append(
                    _canonical_organism_error(
                        config_path, f"{aid} treatment_organism", treatment_organism,
                        all_canonical_organisms,
                    )
                )
            elif treatment_organism:
                print(f"      treatment_organism '{treatment_organism}': OK")

            # ── Canonical test_type validation ────────────────────────────────
            test_type = analysis.get("test_type", "")
            if test_type and test_type not in CANONICAL_TEST_TYPES:
                errors.append(
                    _canonical_field_error(
                        config_path, aid, "test_type", test_type, CANONICAL_TEST_TYPES
                    )
                )
            elif test_type:
                print(f"      test_type '{test_type}': OK")

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
                if df is not None and len(df) > 0:
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

            # Validate significance metadata fields
            prefiltered = analysis.get("prefiltered")
            pvalue_threshold = analysis.get("pvalue_threshold")
            logfc_threshold = analysis.get("logfc_threshold")
            pvalue_asterisk = analysis.get("pvalue_asterisk_in_logfc")

            if prefiltered is not None and not isinstance(prefiltered, bool):
                errors.append(f"{aid}: 'prefiltered' must be a boolean, got {type(prefiltered).__name__}")
            if pvalue_threshold is not None:
                if not isinstance(pvalue_threshold, (int, float)) or pvalue_threshold <= 0 or pvalue_threshold > 1:
                    errors.append(f"{aid}: 'pvalue_threshold' must be a number in (0, 1], got {pvalue_threshold}")
            if logfc_threshold is not None:
                if not isinstance(logfc_threshold, (int, float)) or logfc_threshold < 0:
                    errors.append(f"{aid}: 'logfc_threshold' must be a non-negative number, got {logfc_threshold}")
            if prefiltered and (pvalue_threshold is not None or logfc_threshold is not None):
                warnings.append(f"{aid}: 'prefiltered: true' makes 'pvalue_threshold'/'logfc_threshold' redundant (prefiltered takes priority)")
            if prefiltered and pvalue_asterisk:
                warnings.append(f"{aid}: 'prefiltered: true' conflicts with 'pvalue_asterisk_in_logfc: true' (prefiltered takes priority)")

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
    if all_ids:
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
            print(f"  WARNING: {w}")
    if errors:
        print(f"ERRORS ({len(errors)}):")
        for e in errors:
            print(f"  ERROR: {e}")
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
