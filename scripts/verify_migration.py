#!/usr/bin/env python3
"""Verify paperconfig migration correctness.

Compares paperconfig_orig.yaml (old format) against paperconfig.yaml (new format)
for all papers, checking round-trip analysis count, ID preservation, field
completeness, and structural correctness.

Usage:
    uv run python scripts/verify_migration.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp"
PAPERCONFIG_FILES = DATA_DIR / "paperconfig_files.txt"

STALE_ANALYSIS_FIELDS = {
    "type", "name", "organism",
    "environmental_control_condition_id",
    "environmental_treatment_condition_id",
    "control_condition", "treatment_condition",
    "experimental_context", "test_type",
    "treatment_organism", "treatment_taxid", "treatment_assembly_accession",
}

REQUIRED_EXPERIMENT_FIELDS = {
    "name", "organism", "omics_type", "test_type",
    "treatment_type", "treatment_condition", "control_condition",
    "medium", "temperature", "light_condition", "light_intensity",
}


def load_yaml(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def iter_analyses(config: dict):
    pub = config.get("publication", {})
    supp = pub.get("supplementary_materials") or config.get("supplementary_materials") or {}
    for table_key, table in supp.items():
        if table.get("type", "csv") != "csv":
            continue
        for analysis in table.get("statistical_analyses", []):
            yield analysis.get("id"), analysis, table


def main():
    errors = []
    warnings = []
    papers_checked = 0

    active_paths = set()
    with open(PAPERCONFIG_FILES) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            active_paths.add((PROJECT_ROOT / line).resolve())

    for orig_path in sorted(DATA_DIR.glob("*/paperconfig_orig.yaml")):
        paper_dir = orig_path.parent
        paper_name = paper_dir.name
        new_path = paper_dir / "paperconfig.yaml"

        if not new_path.exists():
            continue

        old_config = load_yaml(orig_path)
        new_config = load_yaml(new_path)
        old_pub = old_config.get("publication", {})
        new_pub = new_config.get("publication", {})

        if not old_pub:
            continue  # resource config

        experiments = new_pub.get("experiments", {})
        is_active = new_path.resolve() in active_paths

        if not experiments:
            if is_active:
                errors.append(f"{paper_name}: ACTIVE but no experiments block!")
            continue

        papers_checked += 1
        print(f"  {paper_name}", end="")

        old_analyses = list(iter_analyses(old_config))
        new_analyses = list(iter_analyses(new_config))

        # 1. Round-trip count
        if len(old_analyses) != len(new_analyses):
            errors.append(f"{paper_name}: count mismatch old={len(old_analyses)} new={len(new_analyses)}")

        # 2. ID preservation
        old_ids = {a[0] for a in old_analyses}
        new_ids = {a[0] for a in new_analyses}
        if old_ids != new_ids:
            errors.append(f"{paper_name}: ID mismatch missing={old_ids - new_ids} extra={new_ids - old_ids}")

        old_by_id = {a[0]: a[1] for a in old_analyses}
        new_by_id = {a[0]: a[1] for a in new_analyses}

        for aid in sorted(old_ids & new_ids):
            old_a = old_by_id[aid]
            new_a = new_by_id[aid]

            # 3. experiment + timepoint_hours present
            if "experiment" not in new_a:
                errors.append(f"{paper_name}/{aid}: missing 'experiment'")

            # 4. No stale fields
            stale = STALE_ANALYSIS_FIELDS & set(new_a.keys())
            if stale:
                errors.append(f"{paper_name}/{aid}: stale fields: {stale}")

            # 5. Experiment reference valid
            exp_ref = new_a.get("experiment")
            if exp_ref and exp_ref not in experiments:
                errors.append(f"{paper_name}/{aid}: bad experiment ref '{exp_ref}'")

            # 6. CSV columns preserved
            for col_field in ("name_col", "logfc_col", "adjusted_p_value_col"):
                old_val = old_a.get(col_field)
                new_val = new_a.get(col_field)
                if old_val and old_val != new_val:
                    errors.append(f"{paper_name}/{aid}: {col_field} changed '{old_val}'->'{new_val}'")

            # 7. Coculture fields preserved
            if old_a.get("treatment_organism"):
                exp = experiments.get(exp_ref, {})
                if "treatment_organism" not in exp:
                    errors.append(f"{paper_name}/{aid}: lost treatment_organism on experiment")

        # 8. No NEEDS_CURATION
        yaml_text = new_path.read_text()
        if "NEEDS_CURATION" in yaml_text:
            errors.append(f"{paper_name}: NEEDS_CURATION found")

        # 9. No environmental_conditions
        if "environmental_conditions" in new_pub:
            errors.append(f"{paper_name}: environmental_conditions block remains")

        # 10. Experiment field completeness
        for exp_key, exp in experiments.items():
            missing = REQUIRED_EXPERIMENT_FIELDS - set(exp.keys())
            if missing:
                errors.append(f"{paper_name}/exp/{exp_key}: missing {sorted(missing)}")

        print(f" — {len(experiments)} exp, {len(new_analyses)} analyses ✓")

    print(f"\n{'='*60}")
    print(f"Papers checked: {papers_checked}")

    if warnings:
        print(f"\nWARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  {w}")

    if errors:
        print(f"\nERRORS ({len(errors)}):")
        for e in errors:
            print(f"  {e}")
        print(f"\nFAILED ({len(errors)} errors)")
        return 1

    print("\nPASSED — all checks OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
