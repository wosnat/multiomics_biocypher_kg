#!/usr/bin/env python
"""Validate a paperconfig.yaml file for the multiomics knowledge graph.

Usage:
    python scripts/validate_paperconfig.py <path_to_paperconfig.yaml>
    python scripts/validate_paperconfig.py --all

Checks:
    - YAML is parseable
    - Required top-level fields exist (or strain-level resource config with no publication block)
    - PDF file exists
    - CSV files exist and are readable
    - experiments block exists and is well-formed (required fields, canonical vocabulary)
    - Each analysis has 'experiment' reference pointing to a valid experiment key
    - Each analysis has 'timepoint_hours' (or null)
    - Column names (name_col, logfc_col, adjusted_p_value_col) match CSV headers
    - id_columns column names match CSV headers for csv and id_translation entries
    - product_columns column names match CSV headers
    - annotation_gff entries have organism and a valid GFF/GTF extension
    - id_translation entries have organism and id_columns
    - All required analysis fields are present
    - Analysis IDs are unique
    - skip_rows is applied correctly when reading CSV headers
    - Sample data rows are readable and logfc_col contains numeric-like values
    - Canonical vocabulary: organism, treatment_type, test_type, omics_type
"""

import re
import sys
import os
import json
from pathlib import Path

import pandas as pd

# Import shared paperconfig utilities from the main package
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))
from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig,
    load_all_paperconfigs,
    get_publication,
    get_paper_name,
    get_supplementary_materials,
)
from multiomics_kg.vocab.non_de_evidence import (
    COMPARTMENTS,
    EXTENDED_OMICS_TYPES,
    VALUE_KINDS,
    KNOWN_METRIC_TYPES,
    DEFAULT_SKIP_TOKENS,
    VALID_BLANK_POLICIES,
)

# Required fields on each statistical_analyses entry (new experiment-based format)
REQUIRED_ANALYSIS_FIELDS = [
    "id", "experiment", "name_col", "logfc_col",
]

# Strict required fields with descriptions (enforced as errors)
REQUIRED_STATS_FIELDS = {
    "id": "unique identifier within publication",
    "experiment": "reference to a key in the experiments block",
}

# Required fields on each experiment entry (always required)
REQUIRED_EXPERIMENT_FIELDS = [
    "name", "organism", "omics_type",
    "treatment_type", "treatment_condition",
]

# Fields required for experiments with DE analyses, but only recommended for
# cluster-only experiments (profiling studies without a control or statistical test)
DE_REQUIRED_EXPERIMENT_FIELDS = ["test_type", "control_condition"]

# Optional but recommended experiment fields (warn if missing)
RECOMMENDED_EXPERIMENT_FIELDS = [
    "medium", "temperature", "light_condition", "light_intensity",
    "table_scope",
]

VALID_TABLE_SCOPES = {
    "all_detected_genes", "significant_any_timepoint",
    "significant_only", "top_n", "filtered_subset",
}

DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")

# Accepted omics_type values on experiment entries.
# Extended set includes PAIRED_RNASEQ_PROTEOME (Waldbauer 2012 et al.).
VALID_TYPES = set(EXTENDED_OMICS_TYPES)
VALID_ID_TYPES = {
    "locus_tag", "locus_tag_ncbi", "locus_tag_cyanorak",
    "old_locus_tag", "alternative_locus_tag",
    "gene_name", "gene_synonym",
    "protein_id_refseq", "uniprot_accession", "uniprot_entry_name",
    "uniprot_annotation_string",  # NEW — extracts entry_name + GN= via gene_id_utils
    "ncbi_protein_defline",       # NEW — extracts NP_/WP_/CAE/etc. from gi|N|src|ACC|... deflines
    "uniprot_defline",            # NEW — extracts <acc> + <entry_name> from sp|<acc>|<entry_name> deflines
    "jgi_id", "probeset", "rast_id", "annotation_specific", "other",
}

# Column names known to hold free-text narrative (protein descriptions,
# product strings, etc.). Declaring such a column as `id_type: locus_tag`
# causes `build_gene_id_mapping` to tokenise entire sentences into
# `specific_lookup`, polluting the strain's gene_id mapping with thousands
# of junk alt_ids (e.g. single words, punctuation) that all collapse onto
# whichever locus tag appeared in the row. This silently collapses
# hundreds of expression edges after dedup. The correct fix is to extract
# the real identifier into a dedicated column (typically via a `GN=`
# regex in `scripts/build_modified_csv/`) and declare THAT column as a
# real id_type.
FREE_TEXT_COLUMN_NAMES = {
    "Description", "description",
    "Product", "product",
    "Name", "name",
    "Protein",
    "Gene Name",
}

# ── Canonical vocabulary ─────────────────────────────────────────────────────

# Canonical organism names for the 'organism' field in experiments
# and for the 'organism' field in supplementary table entries.
CANONICAL_GENOMIC_ORGANISMS = {
    "Prochlorococcus MED4",
    "Prochlorococcus AS9601",
    "Prochlorococcus MIT9301",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL1A",
    "Prochlorococcus NATL2A",
    "Prochlorococcus MIT9303",
    "Prochlorococcus RSP50",
    "Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)",
    "Synechococcus WH8102",
    "Synechococcus WH7803",
    "Synechococcus CC9311",
    "Synechococcus PCC 7002",
    "Synechococcus elongatus PCC 7942",
    "Synechococcus elongatus UTEX 2973",
    "Synechococcus sp. BL107",
    "Thermosynechococcus vestitus BP-1",
    "Alteromonas macleodii HOT1A3",
    "Alteromonas macleodii EZ55",
    "Alteromonas macleodii MIT1002",
    "Alteromonas macleodii ATCC27126",
    "Alteromonas macleodii AD45",
    "Alteromonas macleodii BS11",
    "Alteromonas macleodii BGP6",
    "Alteromonas (MarRef v6)",
    "Marinobacter (MarRef v6)",
    "Shewanella sp. W3-18-1",
    "Pseudomonas putida KT2440",
    "Ruegeria pomeroyi DSS-3",
}

# Canonical treatment organism names loaded from treatment_organisms.csv
# (genus-level organisms used as coculture partners with no loaded genome).
# These are loaded at runtime from the CSV; this set is the fallback default.
_TREATMENT_ORGANISMS_CSV_DEFAULT = {
    "Phage",
    "Alteromonas",
    "Vibrio parahaemolyticus",
    "Meiothermus ruber",
    "Escherichia coli",
    # Ziegler 2025 organisms (commented out in treatment_organisms.csv until paper is added):
    # "Marinobacter",
    # "Thalassospira",
    # "Pseudohoeflea",
}

# Canonical vocabulary shared by treatment_type and background_factors.
# treatment_type = "what environmental variable is being manipulated"
# background_factors = "what conditions are held constant"
# Same canonical terms; meaning depends on which field they appear in.
CANONICAL_CONDITION_TYPES = {
    # Nutrient / environmental variables
    "nitrogen",         # N-limitation / N-replete
    "phosphorus",       # P-limitation / P-replete
    "iron",             # Fe-limitation / Fe-replete
    "carbon",           # CO2, carbon source, chitosan addition
    "light",            # Light quality, intensity; as bg = continuous light regime
    "darkness",         # Extended dark treatment; as bg = dark regime
    "diel",             # Diel light-dark cycling / circadian
    "temperature",      # Thermal shift / acclimation
    "salt",             # Salinity / osmotic
    # Biotic interactions
    "coculture",        # Co-cultivation with another organism
    "viral",            # Phage infection
    # Chemical / xenobiotic
    "chemical",         # Chemical treatment (e.g., DCMU); as bg = inhibitor present
    "plastic",          # Plastic leachate exposure
    # Growth / baseline
    "growth_phase",     # Growth state / multi-condition comparison
    "mutant",           # Mutant or evolved strain comparison
    # Subcellular fractionation comparisons (vesicle/exoproteome/secretome vs whole cell)
    "compartment",      # Subcellular fraction comparison; specific fraction lives in `compartment` field
    # Background-only (typically not used as treatment_type)
    "axenic",           # Pure culture, no other organisms
}

# Valid cluster_type values for gene_clusters entries.
VALID_CLUSTER_TYPES = {
    "time_course",
    "diel",
    "condition_comparison",
    "expression_bin",
}

# Required fields on gene_clusters supplementary entries.
REQUIRED_CLUSTER_TABLE_FIELDS = [
    "name", "filename", "organism", "gene_id_col", "cluster_col", "cluster_type",
]

OPTIONAL_CLUSTER_TABLE_FIELDS = {
    "cluster_method", "omics_type", "light_condition", "treatment_type",
    "background_factors", "treatment", "experimental_context", "experiments",
    "id_columns", "time_points", "skip_rows", "figure_hint", "extraction_notes",
    "score_col", "p_value_col",
}

# Removed: REQUIRED_CLUSTER_FIELDS, RECOMMENDED_CLUSTER_FIELDS
# Per-cluster data comes from extraction JSON, not paperconfig

VALID_GROWTH_PHASES = {
    "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
    "infected", "recovery", "diel", "darkness", "death", "acute_stress", "unknown",
}


def is_valid_growth_phase(value: str) -> bool:
    """Accepts a canonical enum value or a well-formed `other:<slug>` escape value."""
    if not isinstance(value, str):
        return False
    if value in VALID_GROWTH_PHASES:
        return True
    return value.startswith("other:") and len(value) > len("other:")


# Canonical test_type values for statistical_analyses entries.
CANONICAL_TEST_TYPES = {
    # RNA-seq
    "DESeq2",
    "DESeq",
    "edgeR",
    "Rockhopper",
    # Microarray
    "microarray",
    "microarray_Cyber-T",
    "microarray_LPE",
    "microarray_Goldenspike",
    "SAM",                    # Significance Analysis of Microarrays
    # Proteomics / generic
    "t-test",
    "t-test_Perseus",        # Student's t-test via Perseus
    "ANOVA",
    "iTRAQ_t-test",          # iTRAQ quantification + t-test
    "RPKM_fold_change",      # simple RPKM ratio (no formal test)
    "fold_change",            # fold change without formal test
    "LC-MS/MS",               # label-free proteomics quantification
    "read_coverage",          # raw per-gene sequencing read counts (no formal test)
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


def _load_genome_accessions(project_root: str) -> set:
    """Load assembly accessions (insdc.gcf node IDs) from cyanobacteria_genomes.csv.

    Returns an empty set if the file cannot be read; downstream validation will
    then skip the cross-check rather than report spurious errors.
    """
    csv_path = os.path.join(
        project_root,
        "data", "Prochlorococcus", "genomes", "cyanobacteria_genomes.csv",
    )
    try:
        df = pd.read_csv(csv_path, comment="#")
        return set(df["ncbi_accession"].dropna().str.strip().tolist())
    except Exception:
        return set()


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
        # Free-text column declared as locus_tag — this tokenises whole
        # sentences into gene_id_mapping.specific_lookup. See comment on
        # FREE_TEXT_COLUMN_NAMES above.
        if id_type == "locus_tag" and col_name in FREE_TEXT_COLUMN_NAMES:
            warnings.append(
                f"{table_key}.id_columns | column '{col_name}' declared "
                f"id_type: locus_tag looks like free-text description — "
                f"this can pollute gene_id_mapping.specific_lookup with "
                f"thousands of tokenized alt_ids. Consider extracting a "
                f"locus-tag column upstream in the _modified.csv builder."
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


def _validate_experiments(experiments: dict, config_path: str,
                          all_canonical_organisms: set,
                          errors: list, warnings: list,
                          referenced_exp_keys: set | None = None,
                          genome_accessions: set | None = None) -> None:
    """Validate the experiments block in a publication config."""
    if not experiments:
        # experiments block is optional for cluster-only paperconfigs
        warnings.append(f"{config_path} | missing 'experiments' block in publication (OK for cluster-only configs)")
        return

    if not isinstance(experiments, dict):
        errors.append(f"{config_path} | 'experiments' must be a mapping, got {type(experiments).__name__}")
        return

    for exp_key, exp in experiments.items():
        if not isinstance(exp, dict):
            errors.append(f"{config_path} | experiments.{exp_key} must be a mapping")
            continue

        print(f"\n  [experiment: {exp_key}]")

        # Required fields
        for field in REQUIRED_EXPERIMENT_FIELDS:
            val = exp.get(field)
            if val is None or (isinstance(val, str) and not val.strip()):
                errors.append(
                    f"{config_path} | experiments.{exp_key} | "
                    f"missing required field '{field}'"
                )

        # Fields required for DE experiments, warned for cluster-only
        has_analyses = referenced_exp_keys and exp_key in referenced_exp_keys
        for field in DE_REQUIRED_EXPERIMENT_FIELDS:
            val = exp.get(field)
            if val is None or (isinstance(val, str) and not val.strip()):
                if has_analyses:
                    errors.append(
                        f"{config_path} | experiments.{exp_key} | "
                        f"missing required field '{field}'"
                    )
                else:
                    warnings.append(
                        f"{config_path} | experiments.{exp_key} | "
                        f"missing field '{field}' (optional for cluster-only experiments)"
                    )

        # Canonical organism
        organism = exp.get("organism", "")
        if organism and organism not in all_canonical_organisms:
            errors.append(
                _canonical_organism_error(
                    config_path, f"experiments.{exp_key}", organism,
                    all_canonical_organisms,
                )
            )
        elif organism:
            print(f"    organism '{organism}': OK")

        # Canonical omics_type
        omics_type = exp.get("omics_type", "")
        if omics_type and omics_type not in VALID_TYPES:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "omics_type", omics_type, VALID_TYPES,
                )
            )
        elif omics_type:
            print(f"    omics_type '{omics_type}': OK")

        # Canonical test_type
        test_type = exp.get("test_type", "")
        if test_type and test_type not in CANONICAL_TEST_TYPES:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "test_type", test_type, CANONICAL_TEST_TYPES,
                )
            )
        elif test_type:
            print(f"    test_type '{test_type}': OK")

        # Canonical treatment_type (string or list)
        treatment_type = exp.get("treatment_type", "")
        if isinstance(treatment_type, list):
            for tt in treatment_type:
                if tt not in CANONICAL_CONDITION_TYPES:
                    errors.append(
                        _canonical_field_error(
                            config_path, f"experiments.{exp_key}",
                            "treatment_type", tt, CANONICAL_CONDITION_TYPES,
                        )
                    )
                else:
                    print(f"    treatment_type '{tt}': OK")
        elif treatment_type and treatment_type not in CANONICAL_CONDITION_TYPES:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "treatment_type", treatment_type, CANONICAL_CONDITION_TYPES,
                )
            )
        elif treatment_type:
            print(f"    treatment_type '{treatment_type}': OK")

        # Canonical background_factors (optional list)
        background_factors = exp.get("background_factors", [])
        if isinstance(background_factors, str):
            background_factors = [background_factors]
        for bf in background_factors:
            if bf not in CANONICAL_CONDITION_TYPES:
                errors.append(
                    _canonical_field_error(
                        config_path, f"experiments.{exp_key}",
                        "background_factors", bf, CANONICAL_CONDITION_TYPES,
                    )
                )
            else:
                print(f"    background_factors '{bf}': OK")

        # Canonical treatment_organism (optional, only for coculture experiments)
        treatment_organism = exp.get("treatment_organism", "")
        if treatment_organism and treatment_organism not in all_canonical_organisms:
            errors.append(
                _canonical_organism_error(
                    config_path, f"experiments.{exp_key} treatment_organism",
                    treatment_organism, all_canonical_organisms,
                )
            )
        elif treatment_organism:
            print(f"    treatment_organism '{treatment_organism}': OK")

        # Canonical compartment (optional; default "whole_cell" at ingest time)
        compartment = exp.get("compartment", "")
        if compartment and compartment not in COMPARTMENTS:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "compartment", compartment, COMPARTMENTS,
                )
            )
        elif compartment:
            print(f"    compartment '{compartment}': OK")

        # Validate table_scope enum
        table_scope = exp.get("table_scope", "")
        if table_scope and table_scope not in VALID_TABLE_SCOPES:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "table_scope", table_scope, VALID_TABLE_SCOPES,
                )
            )
        elif table_scope:
            print(f"    table_scope '{table_scope}': OK")

        # Recommended fields (warn if missing)
        for field in RECOMMENDED_EXPERIMENT_FIELDS:
            if field not in exp:
                warnings.append(
                    f"{config_path} | experiments.{exp_key} | "
                    f"missing recommended field '{field}'"
                )

        # Formatting checks
        light_cond = exp.get("light_condition", "")
        if light_cond and "_" in light_cond:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"light_condition '{light_cond}' uses underscores — "
                f"use spaces instead (e.g., 'continuous light')"
            )
        light_int = exp.get("light_intensity", "")
        if light_int and "µ" in str(light_int):
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"light_intensity uses Unicode µ — use ASCII 'umol' instead"
            )

        # Coculture/viral consistency checks
        raw_tt = exp.get("treatment_type", "")
        tt_list = raw_tt if isinstance(raw_tt, list) else ([raw_tt] if raw_tt else [])
        t_org = exp.get("treatment_organism", "")
        if ("coculture" in tt_list or "viral" in tt_list) and not t_org:
            errors.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_type includes 'coculture' or 'viral' but missing treatment_organism"
            )
        if t_org and "coculture" not in tt_list and "viral" not in tt_list:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"has treatment_organism '{t_org}' but treatment_type {tt_list} "
                f"does not include 'coculture' or 'viral'"
            )
        if t_org and str(t_org).strip().lower() in ("phage", "virus", "bacteriophage") and "viral" not in tt_list:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_organism is '{t_org}' but treatment_type {tt_list} "
                f"does not include 'viral' (should it?)"
            )
        if ("coculture" in tt_list or "viral" in tt_list) and "treatment_taxid" not in exp:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_type includes 'coculture' or 'viral' but missing treatment_taxid"
            )

        # treatment_assembly_accession must match an existing OrganismTaxon node id.
        # The omics adapter wraps this value as `insdc.gcf:<accession>` and creates
        # the Tests_coculture_with edge target from it. If the accession is not in
        # cyanobacteria_genomes.csv, neo4j-admin import drops the edge silently
        # (only a line in import.report). See Kratzl 2024 incident, 2026-04-14.
        treatment_acc = exp.get("treatment_assembly_accession", "")
        if treatment_acc and genome_accessions and treatment_acc not in genome_accessions:
            errors.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_assembly_accession '{treatment_acc}' is not present in "
                f"data/Prochlorococcus/genomes/cyanobacteria_genomes.csv | "
                f"the Tests_coculture_with edge will refer to a missing OrganismTaxon "
                f"node and be dropped at import time"
            )
        elif treatment_acc:
            print(f"    treatment_assembly_accession '{treatment_acc}': OK")


def _validate_organism_consistency(
    experiments: dict, supp: dict, config_path: str,
    errors: list, warnings: list,
) -> None:
    """Cross-validate organisms across experiments and supplementary materials.

    Checks:
    1. Every organism in supp materials (id_translation, annotation_gff) appears
       as an experiment organism (primary, not treatment_organism).
    2. Every experiment organism has at least one analysis referencing it.
    3. Every experiment treatment_organism has at least one analysis referencing
       an experiment that uses it.
    """
    # Collect experiment organisms and treatment organisms
    exp_organisms = {}  # organism -> set of experiment keys
    exp_treatment_organisms = {}  # treatment_organism -> set of experiment keys
    for exp_key, exp in experiments.items():
        if not isinstance(exp, dict):
            continue
        org = exp.get("organism", "")
        if org:
            exp_organisms.setdefault(org, set()).add(exp_key)
        t_org = exp.get("treatment_organism", "")
        if t_org:
            exp_treatment_organisms.setdefault(t_org, set()).add(exp_key)

    # Collect organisms from supplementary materials
    supp_organisms = {}  # organism -> list of table keys
    for table_key, table in supp.items():
        if not isinstance(table, dict):
            continue
        table_type = table.get("type", "csv")
        org = table.get("organism", "")
        if org and table_type in ("id_translation", "annotation_gff"):
            supp_organisms.setdefault(org, []).append(table_key)

    # Collect which experiment keys are referenced by analyses
    referenced_exp_keys = set()
    for table_key, table in supp.items():
        if not isinstance(table, dict):
            continue
        for analysis in table.get("statistical_analyses", []):
            if not isinstance(analysis, dict):
                continue
            exp_ref = analysis.get("experiment", "")
            if exp_ref:
                referenced_exp_keys.add(exp_ref)

    print(f"\n  [organism consistency]")
    all_exp_orgs = set(exp_organisms.keys())
    all_exp_treatment_orgs = set(exp_treatment_organisms.keys())
    all_supp_orgs = set(supp_organisms.keys())
    print(f"    experiment organisms: {sorted(all_exp_orgs)}")
    if all_exp_treatment_orgs:
        print(f"    experiment treatment organisms: {sorted(all_exp_treatment_orgs)}")
    if all_supp_orgs:
        print(f"    supplementary material organisms: {sorted(all_supp_orgs)}")

    # Check 1: supp material organisms should be experiment organisms
    for org, table_keys in supp_organisms.items():
        if org not in all_exp_orgs:
            warnings.append(
                f"{config_path} | organism '{org}' in supplementary materials "
                f"({', '.join(table_keys)}) is not an experiment organism "
                f"(experiment organisms: {sorted(all_exp_orgs)})"
            )

    # Check 2: every experiment organism should have at least one analysis
    for org, exp_keys in exp_organisms.items():
        has_analysis = any(ek in referenced_exp_keys for ek in exp_keys)
        if not has_analysis:
            warnings.append(
                f"{config_path} | experiment organism '{org}' "
                f"(experiments: {sorted(exp_keys)}) has no analyses referencing it"
            )

    # Check 3: every experiment treatment_organism should have at least one analysis
    for t_org, exp_keys in exp_treatment_organisms.items():
        has_analysis = any(ek in referenced_exp_keys for ek in exp_keys)
        if not has_analysis:
            warnings.append(
                f"{config_path} | experiment treatment_organism '{t_org}' "
                f"(experiments: {sorted(exp_keys)}) has no analyses referencing it"
            )

    if not (all_supp_orgs - all_exp_orgs):
        print("    supp material organisms match experiment organisms: OK")


def _find_project_root() -> Path:
    """Return project root by locating data/Prochlorococcus/treatment_organisms.csv."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "data" / "Prochlorococcus" / "treatment_organisms.csv").exists():
            return parent
    return Path.cwd()


def _validate_gene_clusters_entry(key, table, config, paperconfig_dir,
                                   all_canonical_organisms, errors, warnings):
    """Validate a type: gene_clusters entry."""

    # Required fields
    for field in REQUIRED_CLUSTER_TABLE_FIELDS:
        if field not in table:
            errors.append(f"  [{key}] Missing required field: {field}")

    # Organism check
    organism = table.get("organism", "")
    if organism:
        print(f"  Organism: {organism}")
        if organism not in all_canonical_organisms:
            warnings.append(f"  [{key}] Organism '{organism}' not in canonical list")

    # cluster_type enum
    ct = table.get("cluster_type", "")
    if ct and ct not in VALID_CLUSTER_TYPES:
        warnings.append(f"  [{key}] cluster_type '{ct}' not in {VALID_CLUSTER_TYPES}")

    # omics_type
    ot = table.get("omics_type", "")
    if ot and ot not in VALID_TYPES:
        warnings.append(f"  [{key}] omics_type '{ot}' not in {VALID_TYPES}")

    # treatment_type must be a list
    tt = table.get("treatment_type")
    if tt is not None:
        if not isinstance(tt, list):
            errors.append(f"  [{key}] treatment_type must be a list, got {type(tt).__name__}")
        else:
            for t in tt:
                if t not in CANONICAL_CONDITION_TYPES:
                    warnings.append(f"  [{key}] treatment_type '{t}' not canonical")

    # Entry key naming
    if key.startswith("cluster_table_"):
        warnings.append(
            f"  [{key}] Entry key should be a meaningful ID "
            f"(e.g., 'med4_kmeans_nstarvation'), not generic"
        )

    # experiments cross-reference
    experiments_ref = table.get("experiments")
    if experiments_ref is None:
        warnings.append(
            f"  [{key}] No 'experiments' field — cluster analysis should link "
            f"to at least one experiment in publication.experiments"
        )
    elif not isinstance(experiments_ref, list):
        errors.append(f"  [{key}] experiments must be a list")
    else:
        if len(experiments_ref) == 0:
            warnings.append(
                f"  [{key}] 'experiments' list is empty — should reference "
                f"at least one experiment key"
            )
        pub_experiments = config.get("publication", {}).get("experiments", {})
        for exp_key in experiments_ref:
            if exp_key not in pub_experiments:
                errors.append(
                    f"  [{key}] experiments references '{exp_key}' not found "
                    f"in publication.experiments"
                )

    # CSV validation
    filename = table.get("filename", "")
    if filename:
        csv_path = Path(filename)
        if not csv_path.is_absolute():
            csv_path = _find_project_root() / csv_path
        if csv_path.exists():
            try:
                sep = table.get("sep", ",")
                skip = table.get("skip_rows", 0)
                df = pd.read_csv(csv_path, sep=sep, skiprows=skip, nrows=5)
                cols = set(df.columns)
                for col_field in ("gene_id_col", "cluster_col", "score_col"):
                    col_name = table.get(col_field)
                    if col_name and col_name not in cols:
                        errors.append(
                            f"  [{key}] {col_field}='{col_name}' not found in "
                            f"CSV columns: {sorted(cols)}"
                        )
            except Exception as e:
                warnings.append(f"  [{key}] Could not read CSV: {e}")
        else:
            warnings.append(f"  [{key}] CSV file not found: {csv_path}")

    # Extraction JSON validation
    if paperconfig_dir:
        json_path = Path(paperconfig_dir) / f"cluster_extraction_{key}.json"
        if json_path.exists():
            try:
                with open(json_path) as f:
                    extraction = json.load(f)
                meta_key = extraction.get("metadata", {}).get("table_key", "")
                if meta_key != key:
                    warnings.append(
                        f"  [{key}] Extraction JSON table_key='{meta_key}' "
                        f"does not match entry key"
                    )
                stage3 = extraction.get("stage3_validation", {})
                for ck, cv in stage3.items():
                    if isinstance(cv, dict) and cv.get("verdict") == "fail":
                        warnings.append(
                            f"  [{key}] Extraction cluster '{ck}' has verdict=fail"
                        )
            except Exception as e:
                warnings.append(f"  [{key}] Could not read extraction JSON: {e}")

    # Warn if clusters: block still present (deprecated)
    if "clusters" in table:
        warnings.append(
            f"  [{key}] 'clusters' block is deprecated — "
            f"per-cluster data comes from extraction JSON"
        )


def _validate_derived_metrics_entry(
    key: str, table: dict, config_path: str,
    all_canonical_organisms: set,
    errors: list, warnings: list,
) -> None:
    """Validate a type: derived_metrics_table supplementary entry.

    Design (option C): per-metric metadata (rankable, has_p_value,
    p_value_threshold, unit, allowed_categories) lives inline on each metric
    entry — it's paper-specific. KNOWN_METRIC_TYPES only pins metric_type →
    value_kind to catch edge-type drift across papers.

    Validation steps per metric:
      1. required top-level fields present (filename, organism, experiment,
         name_col, metrics)
      2. organism is in canonical list
      3. CSV file readable; name_col + value_col exist as columns
      4. metric_type non-empty; value_kind ∈ VALUE_KINDS
      5. if metric_type ∈ KNOWN_METRIC_TYPES: declared value_kind must match;
         if not: warn (novel name)
      6. field_description: required non-empty string (free text) — surfaced
         by downstream MCP tools
      7. per-value_kind required / forbidden field gating:
           numeric     — require value_col, rankable, has_p_value;
                         rankable/has_p_value ∈ {"true","false"};
                         if has_p_value=="true": require p_value_threshold ∈ (0,1]
                           and at least one of p_value_col/adjusted_p_value_col;
                         if has_p_value=="false": those three forbidden;
                         unit optional (free string or omitted);
                         forbidden: true_tokens/false_tokens/skip_tokens/
                                   blank_policy/allowed_categories
           boolean     — require value_col, true_tokens (non-empty list);
                         optional false_tokens / skip_tokens / blank_policy
                         (∈ VALID_BLANK_POLICIES, default "skip");
                         forbidden: unit, rankable, has_p_value, p_value_*,
                                   allowed_categories
                                   (rankable/has_p_value are definitionally
                                   "false" for boolean — adapter sets them
                                   at ingest, paperconfig must not declare);
                         CSV dry-run: hard-error on any unclassified cell value
           categorical — require value_col, allowed_categories (non-empty list);
                         forbidden: unit, rankable, has_p_value,
                                   true_tokens/false_tokens/skip_tokens/
                                   blank_policy, p_value_*
                                   (rankable/has_p_value are definitionally
                                   "false" for categorical — adapter sets them
                                   at ingest, paperconfig must not declare);
                         CSV dry-run: warn on cells outside allowed_categories
                         (adapter hard-errors at ingest in Plan 2)
    """
    # ── Required top-level fields ──
    for req in ("filename", "organism", "experiment", "name_col", "metrics"):
        if req not in table or table.get(req) in (None, "", []):
            errors.append(
                f"{config_path} | {key} | derived_metrics_table missing required "
                f"field '{req}'"
            )

    # ── organism vocab ──
    organism = table.get("organism", "")
    if organism and organism not in all_canonical_organisms:
        errors.append(
            _canonical_organism_error(config_path, key, organism, all_canonical_organisms)
        )

    # ── experiment reference type-check (cross-reference happens elsewhere) ──
    exp_ref = table.get("experiment", "")
    if exp_ref and not isinstance(exp_ref, str):
        errors.append(f"{config_path} | {key} | 'experiment' must be a string key")

    # ── CSV load for dry-runs ──
    filename = table.get("filename", "")
    sep = table.get("sep", ",")
    skip = int(table.get("skip_rows", 0) or 0)
    df = None
    cols: list[str] = []
    if filename and os.path.exists(filename):
        try:
            df = pd.read_csv(filename, sep=sep, skiprows=skip if skip else None, dtype=str)
            cols = list(df.columns)
            print(f"    columns: {cols}")
        except Exception as e:
            errors.append(f"{config_path} | {key} | cannot read CSV: {e}")
    elif filename:
        errors.append(f"{config_path} | {key} | file not found: {filename}")

    # ── name_col presence ──
    name_col = table.get("name_col", "")
    if name_col and cols and name_col not in cols:
        errors.append(f"{config_path} | {key} | name_col '{name_col}' not in CSV columns")

    # ── id_columns validation (reuse existing helper) ──
    id_columns = table.get("id_columns", []) or []
    if cols:
        _validate_id_columns(id_columns, cols, key, errors, warnings)

    # ── metrics ──
    metrics = table.get("metrics") or []
    if not isinstance(metrics, list) or not metrics:
        errors.append(f"{config_path} | {key} | 'metrics' must be a non-empty list")
        return

    for i, m in enumerate(metrics):
        ctx = f"{key}.metrics[{i}]"
        if not isinstance(m, dict):
            errors.append(f"{config_path} | {ctx} | must be a mapping")
            continue

        metric_type = m.get("metric_type", "")
        value_kind = m.get("value_kind", "")
        value_col = m.get("value_col", "")

        # ── metric_type presence ──
        if not metric_type:
            errors.append(f"{config_path} | {ctx} | missing metric_type")
            continue

        # ── value_kind presence + vocabulary ──
        if not value_kind:
            errors.append(f"{config_path} | {ctx} | missing value_kind")
            continue
        if value_kind not in VALUE_KINDS:
            errors.append(
                f"{config_path} | {ctx} | value_kind '{value_kind}' not in "
                f"{sorted(VALUE_KINDS)}"
            )
            continue

        # ── registry cross-check (known metric_type must match value_kind) ──
        if metric_type in KNOWN_METRIC_TYPES:
            expected_kind = KNOWN_METRIC_TYPES[metric_type]
            if value_kind != expected_kind:
                errors.append(
                    f"{config_path} | {ctx} | metric_type '{metric_type}' is "
                    f"registered as value_kind='{expected_kind}' in KNOWN_METRIC_TYPES; "
                    f"paperconfig declares value_kind='{value_kind}'. Use a different "
                    f"metric_type name if this paper measures a different kind."
                )
                continue
        else:
            warnings.append(
                f"{config_path} | {ctx} | metric_type '{metric_type}' is novel "
                f"(not in KNOWN_METRIC_TYPES). Consider adding it to "
                f"multiomics_kg/vocab/non_de_evidence.py if future papers will reuse it."
            )

        # ── value_col presence + CSV column match ──
        if not value_col:
            errors.append(f"{config_path} | {ctx} | missing value_col")
        elif cols and value_col not in cols:
            errors.append(
                f"{config_path} | {ctx} | value_col '{value_col}' not in CSV columns"
            )

        # ── field_description: required non-empty string (free text) ──
        # Human-readable explanation surfaced by downstream MCP tools;
        # without it a DerivedMetric node is opaque.
        fd = m.get("field_description")
        if not isinstance(fd, str) or not fd.strip():
            errors.append(
                f"{config_path} | {ctx} | 'field_description' is required and must be "
                f"a non-empty string (free text; no vocabulary constraint)"
            )

        # ── per-value_kind gating ──
        if value_kind == "numeric":
            # rankable + has_p_value required and must be "true"|"false"
            for flag_field in ("rankable", "has_p_value"):
                if flag_field not in m:
                    errors.append(
                        f"{config_path} | {ctx} | '{flag_field}' is required for "
                        f"value_kind=numeric (\"true\" or \"false\")"
                    )
                elif m[flag_field] not in ("true", "false"):
                    errors.append(
                        f"{config_path} | {ctx} | '{flag_field}' must be \"true\" "
                        f"or \"false\" (string enum), got {m[flag_field]!r}"
                    )
            # p-value field gating
            has_p = m.get("has_p_value")
            threshold = m.get("p_value_threshold")
            pval_col = m.get("p_value_col")
            adj_col = m.get("adjusted_p_value_col")
            if has_p == "true":
                if threshold is None:
                    errors.append(
                        f"{config_path} | {ctx} | p_value_threshold is required "
                        f"when has_p_value='true'"
                    )
                elif not isinstance(threshold, (int, float)) or not (0 < threshold <= 1):
                    errors.append(
                        f"{config_path} | {ctx} | p_value_threshold must be a "
                        f"number in (0, 1], got {threshold!r}"
                    )
                if not pval_col and not adj_col:
                    errors.append(
                        f"{config_path} | {ctx} | at least one of p_value_col / "
                        f"adjusted_p_value_col is required when has_p_value='true'"
                    )
                for col_field in ("p_value_col", "adjusted_p_value_col"):
                    declared = m.get(col_field)
                    if declared and cols and declared not in cols:
                        errors.append(
                            f"{config_path} | {ctx} | {col_field} '{declared}' "
                            f"not in CSV columns"
                        )
            elif has_p == "false":
                for forbid in ("p_value_threshold", "p_value_col", "adjusted_p_value_col"):
                    if forbid in m:
                        errors.append(
                            f"{config_path} | {ctx} | '{forbid}' is forbidden when "
                            f"has_p_value='false'"
                        )
            # Forbidden-field gating
            for forbid in ("true_tokens", "false_tokens", "skip_tokens",
                           "blank_policy", "allowed_categories"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=numeric"
                    )

        elif value_kind == "boolean":
            # true_tokens required
            true_tokens = m.get("true_tokens")
            if not isinstance(true_tokens, list) or not true_tokens:
                errors.append(
                    f"{config_path} | {ctx} | 'true_tokens' is required and must be "
                    f"a non-empty list for value_kind=boolean"
                )
            # Forbidden-field gating. rankable/has_p_value are definitionally
            # "false" for boolean — adapter sets them at ingest, paperconfig
            # must not declare them (forbids redundant / copy-paste noise).
            for forbid in ("unit", "rankable", "has_p_value",
                           "p_value_col", "adjusted_p_value_col",
                           "p_value_threshold", "allowed_categories"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=boolean"
                    )
            # blank_policy vocab
            blank_policy = m.get("blank_policy", "skip")
            if blank_policy not in VALID_BLANK_POLICIES:
                errors.append(
                    f"{config_path} | {ctx} | blank_policy '{blank_policy}' not in "
                    f"{list(VALID_BLANK_POLICIES)}"
                )
            # CSV dry-run: every non-blank cell must be classified
            if df is not None and value_col and value_col in df.columns:
                true_set = set(m.get("true_tokens") or [])
                false_set = set(m.get("false_tokens") or [])
                skip_set = set(m.get("skip_tokens") or list(DEFAULT_SKIP_TOKENS))
                classified = true_set | false_set | skip_set
                distinct = set()
                for v in df[value_col].fillna("").astype(str).str.strip().tolist():
                    if v == "":
                        continue
                    distinct.add(v)
                unknown = distinct - classified
                if unknown:
                    errors.append(
                        f"{config_path} | {ctx} | unclassified token(s) in value_col "
                        f"'{value_col}': {sorted(unknown)}. Add them to true_tokens, "
                        f"false_tokens, or skip_tokens."
                    )

        elif value_kind == "categorical":
            # allowed_categories required
            allowed_cats = m.get("allowed_categories")
            if not isinstance(allowed_cats, list) or not allowed_cats:
                errors.append(
                    f"{config_path} | {ctx} | 'allowed_categories' is required and "
                    f"must be a non-empty list for value_kind=categorical"
                )
            # Forbidden-field gating. rankable/has_p_value are definitionally
            # "false" for categorical — adapter sets them at ingest, paperconfig
            # must not declare them (forbids redundant / copy-paste noise).
            for forbid in ("unit", "rankable", "has_p_value",
                           "true_tokens", "false_tokens", "skip_tokens",
                           "blank_policy", "p_value_col", "adjusted_p_value_col",
                           "p_value_threshold"):
                if forbid in m:
                    errors.append(
                        f"{config_path} | {ctx} | field '{forbid}' is not allowed "
                        f"for value_kind=categorical"
                    )
            # CSV dry-run: warn on cells outside allowed_categories
            if (df is not None and value_col and value_col in df.columns
                    and isinstance(allowed_cats, list) and allowed_cats):
                allowed = set(allowed_cats)
                distinct = set()
                for v in df[value_col].fillna("").astype(str).str.strip().tolist():
                    if v:
                        distinct.add(v)
                unknown = distinct - allowed
                if unknown:
                    warnings.append(
                        f"{config_path} | {ctx} | value_col '{value_col}' contains "
                        f"{len(unknown)} value(s) outside declared allowed_categories: "
                        f"{sorted(unknown)}. Adapter will hard-error on these rows at "
                        f"ingest (Plan 2). Either expand allowed_categories or fix the CSV."
                    )


def validate_paperconfig_content(
    config: dict, config_path: str,
) -> tuple[list[str], list[str]]:
    """Validate a parsed paperconfig dict; returns (errors, warnings).

    Pure function — no stdout prints about results, no SystemExit. Callers
    render the diagnostics they want.

    NOTE: the current implementation still prints progress lines via
    print(); those are debug output, not results, and remain for parity
    with the previous CLI behavior when invoked via validate().
    """
    errors: list[str] = []
    warnings: list[str] = []

    abs_config = os.path.abspath(config_path)
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
    all_canonical_organisms = CANONICAL_GENOMIC_ORGANISMS | treatment_organisms
    genome_accessions = _load_genome_accessions(_project_root or "")

    # Detect strain-level resource config (no 'publication' block)
    is_resource_config = bool(config) and "publication" not in config

    if is_resource_config:
        print("  [strain-level resource config — no publication block]")
        pub = {}
        supp = config.get("supplementary_materials")
        experiments = {}
    else:
        if not config or "publication" not in config:
            errors.append(f"{config_path} | Missing top-level 'publication' key")
            return errors, warnings

        pub = get_publication(config)

        # Validate extraction section (optional)
        extraction = config.get("extraction", {})
        if extraction:
            valid_extraction_keys = {"scope", "additional_pdfs"}
            for ek in extraction:
                if ek not in valid_extraction_keys:
                    warnings.append(f"  extraction: unknown key '{ek}'")
            scope = extraction.get("scope", "analysis")
            if scope not in ("paper", "analysis"):
                errors.append(f"  extraction.scope must be 'paper' or 'analysis', got '{scope}'")
            additional_pdfs = extraction.get("additional_pdfs", [])
            if additional_pdfs:
                if not isinstance(additional_pdfs, list):
                    errors.append("  extraction.additional_pdfs must be a list")
                else:
                    for pdf_path in additional_pdfs:
                        if not Path(pdf_path).exists():
                            warnings.append(f"  extraction.additional_pdfs: file not found: {pdf_path}")

        # --- papername ---
        papername = get_paper_name(config)
        if papername == "unknown":
            errors.append("Missing 'papername'")
        else:
            print(f"  papername: {papername}")

        # --- DOI override ---
        doi_override = pub.get("doi")
        if doi_override is not None:
            if not isinstance(doi_override, str) or not DOI_RE.match(doi_override.strip()):
                errors.append(f"publication.doi '{doi_override}' is not a valid DOI (expected pattern 10.NNNN/...)")
            else:
                print(f"  doi (override): {doi_override.strip()}")

        # --- PDF ---
        pdf = pub.get("papermainpdf", "")
        if not pdf:
            warnings.append("Missing 'papermainpdf' (optional but recommended)")
        elif not os.path.exists(pdf):
            errors.append(f"PDF not found: {pdf}")
        else:
            print(f"  PDF exists: {pdf}")

        # --- Experiments block ---
        experiments = pub.get("experiments", {})
        supp = get_supplementary_materials(config)

        # Pre-compute which experiment keys have DE analyses referencing them
        _referenced_exp_keys = set()
        if supp:
            for _tk, _tv in supp.items():
                if not isinstance(_tv, dict):
                    continue
                for _a in _tv.get("statistical_analyses", []):
                    if isinstance(_a, dict) and _a.get("experiment"):
                        _referenced_exp_keys.add(_a["experiment"])

        _validate_experiments(
            experiments, config_path, all_canonical_organisms,
            errors, warnings, referenced_exp_keys=_referenced_exp_keys,
            genome_accessions=genome_accessions,
        )

    if not supp:
        errors.append(f"{config_path} | Missing 'supplementary_materials'")
        return errors, warnings

    all_ids = []

    for table_key, table in supp.items():
        table_type = table.get("type", "csv")
        print(f"\n  [{table_key}] (type: {table_type})")

        # ── derived_metrics_table (Plan 1 / non-DE-evidence slice) ──────────
        # Dispatched BEFORE the generic filename short-circuits so the helper's
        # detailed required-field validation runs for every DM entry, including
        # entries with missing/bad filename.
        if table_type == "derived_metrics_table":
            _validate_derived_metrics_entry(
                table_key, table, config_path,
                all_canonical_organisms, errors, warnings,
            )
            continue

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

        # ── gene_clusters ───────────────────────────────────────────────────
        if table_type == "gene_clusters":
            paperconfig_dir = os.path.dirname(os.path.abspath(config_path))
            _validate_gene_clusters_entry(
                table_key, table, config, paperconfig_dir,
                all_canonical_organisms, errors, warnings,
            )
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

            # ── Required fields ──────────────────────────────────────────────
            for req_field, req_desc in REQUIRED_STATS_FIELDS.items():
                val = analysis.get(req_field)
                if val is None or (isinstance(val, str) and not val.strip()):
                    errors.append(
                        f"{config_path} | {aid} | missing required field '{req_field}' "
                        f"({req_desc})"
                    )

            for field in REQUIRED_ANALYSIS_FIELDS:
                if field not in analysis:
                    errors.append(f"{aid}: Missing required field '{field}'")

            # ── Stale old-format fields check ─────────────────────────────────
            STALE_FIELDS = {
                "type", "name", "organism", "test_type",
                "control_condition", "treatment_condition",
                "experimental_context",
                "environmental_control_condition_id",
                "environmental_treatment_condition_id",
                "treatment_organism", "treatment_taxid",
                "treatment_assembly_accession",
            }
            stale = STALE_FIELDS & set(analysis.keys())
            if stale:
                warnings.append(
                    f"{aid}: has old-format fields that should be on the experiment: "
                    f"{sorted(stale)}"
                )

            # ── Experiment reference validation ──────────────────────────────
            exp_ref = analysis.get("experiment", "")
            if exp_ref:
                if experiments and exp_ref not in experiments:
                    errors.append(
                        f"{aid}: experiment reference '{exp_ref}' not found in "
                        f"experiments block (valid keys: {list(experiments.keys())})"
                    )
                elif experiments:
                    print(f"      experiment '{exp_ref}': OK")

            # ── timepoint_hours validation ───────────────────────────────────
            if "timepoint_hours" not in analysis:
                warnings.append(
                    f"{aid}: missing 'timepoint_hours' (set to null if not a "
                    "time-course experiment)"
                )
            else:
                tp_hours = analysis.get("timepoint_hours")
                if tp_hours is not None and not isinstance(tp_hours, (int, float)):
                    errors.append(
                        f"{aid}: 'timepoint_hours' must be a number or null, "
                        f"got {type(tp_hours).__name__}"
                    )
                else:
                    print(f"      timepoint_hours: {tp_hours}")

                # Cross-check: if timepoint is set, verify timepoint_hours
                # matches parse_timepoint_hours()
                timepoint = analysis.get("timepoint")
                if timepoint:
                    from multiomics_kg.utils.paperconfig_utils import parse_timepoint_hours
                    expected = parse_timepoint_hours(str(timepoint))
                    if expected is not None and tp_hours is not None:
                        if abs(float(expected) - float(tp_hours)) > 0.01:
                            warnings.append(
                                f"{aid}: timepoint_hours={tp_hours} doesn't match "
                                f"parse_timepoint_hours('{timepoint}')={expected}"
                            )

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

            # Validate fold_change_type
            fold_change_type = analysis.get("fold_change_type")
            if fold_change_type is not None and fold_change_type not in ("log2", "linear"):
                errors.append(f"{aid}: 'fold_change_type' must be 'log2' or 'linear', got '{fold_change_type}'")

            # --- Timepoint / growth_phase backfill checks (2026-04-12) ---
            if "timepoint" not in analysis:
                warnings.append(
                    f"{aid}: missing 'timepoint' (use 'unknown' if paper doesn't state)"
                )

            # NOTE: Use WARNINGS (not errors) for missing growth_phase during rollout.
            # Task 20 will flip this to errors once all 30 paperconfigs are backfilled.
            if "growth_phase" not in analysis:
                # TODO (timepoint-backfill rollout, 2026-04-12): flip to errors.append
                # once all 30 paperconfigs have growth_phase populated. Tracked in
                # docs/superpowers/plans/2026-04-12-timepoint-growth-phase-backfill.md
                warnings.append(
                    f"{aid}: missing required field 'growth_phase' "
                    f"(one of {sorted(VALID_GROWTH_PHASES)} or 'other:<slug>')"
                )
            elif not is_valid_growth_phase(analysis["growth_phase"]):
                errors.append(
                    f"{aid}: invalid growth_phase '{analysis['growth_phase']}' "
                    f"(must be in VALID_GROWTH_PHASES or start with 'other:')"
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

    # --- Organism consistency check ---
    if experiments and supp:
        _validate_organism_consistency(
            experiments, supp, config_path, errors, warnings,
        )

    return errors, warnings


def validate(config_path: str) -> bool:
    """CLI entry: parse YAML, validate, print results, return success bool."""
    try:
        config = load_paperconfig(Path(config_path))
    except Exception as e:
        print(f"FAIL: Cannot parse YAML: {e}")
        return False

    errors, warnings = validate_paperconfig_content(config, config_path)
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
    if len(sys.argv) == 2 and sys.argv[1] == "--all":
        # Validate all paperconfigs from paperconfig_files.txt
        all_configs = load_all_paperconfigs()
        if not all_configs:
            print("No paperconfigs found in paperconfig_files.txt")
            sys.exit(1)
        total = len(all_configs)
        passed = 0
        failed = 0
        for path, _ in all_configs:
            print(f"\n{'='*60}")
            print(f"Validating: {path}\n")
            ok = validate(str(path))
            if ok:
                passed += 1
            else:
                failed += 1
        print(f"\n{'='*60}")
        print(f"Results: {passed}/{total} passed, {failed}/{total} failed")
        sys.exit(0 if failed == 0 else 1)

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_paperconfig.yaml>")
        print(f"       {sys.argv[0]} --all")
        sys.exit(2)

    config_path = sys.argv[1]
    if not os.path.exists(config_path):
        print(f"File not found: {config_path}")
        sys.exit(1)

    print(f"Validating: {config_path}\n")
    ok = validate(config_path)
    sys.exit(0 if ok else 1)
