"""
Tests that validate all paperconfig.yaml files listed in paperconfig_files.txt.

Reads the file list from data/Prochlorococcus/papers_and_supp/paperconfig_files.txt
and runs the validate_paperconfig.validate() function on each one, ensuring every
config passes validation.

Also contains unit tests for the validator's vocabulary enforcement (canonical
organism names, condition_type, test_type, and required statistical_analyses
fields). These tests use synthetic YAML content via temporary files and do not
depend on the real paperconfig files.
"""

import os
import sys
import textwrap
import pytest
import yaml
from pathlib import Path

# Add the validate script to the import path
VALIDATE_SCRIPT_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, ".claude", "skills", "paperconfig"
)
sys.path.insert(0, os.path.abspath(VALIDATE_SCRIPT_DIR))

from validate_paperconfig import validate, CANONICAL_GENOMIC_ORGANISMS, CANONICAL_CONDITION_TYPES, CANONICAL_TEST_TYPES, REQUIRED_EXPERIMENT_FIELDS

# Project root (one level up from tests/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

PAPERCONFIG_LIST = PROJECT_ROOT / "data" / "Prochlorococcus" / "papers_and_supp" / "paperconfig_files.txt"


def _load_paperconfig_paths() -> list[str]:
    """Read paperconfig file paths from the listing file."""
    assert PAPERCONFIG_LIST.exists(), f"Paperconfig list file not found: {PAPERCONFIG_LIST}"
    paths = []
    with open(PAPERCONFIG_LIST) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(line)
    assert len(paths) > 0, "No paperconfig paths found in listing file"
    return paths


PAPERCONFIG_PATHS = _load_paperconfig_paths()


@pytest.mark.parametrize(
    "config_path",
    PAPERCONFIG_PATHS,
    ids=[Path(p).parent.name for p in PAPERCONFIG_PATHS],
)
def test_paperconfig_validates(config_path: str, monkeypatch):
    """Each paperconfig.yaml listed in paperconfig_files.txt must pass validation."""
    # Run validation from the project root so relative paths in configs resolve
    monkeypatch.chdir(PROJECT_ROOT)

    assert os.path.exists(config_path), f"Config file not found: {config_path}"
    result = validate(config_path)
    assert result is True, f"Validation failed for {config_path}"


# ---------------------------------------------------------------------------
# Helpers for vocabulary unit tests
# ---------------------------------------------------------------------------

def _write_minimal_csv(tmp_path: Path) -> Path:
    """Write a minimal two-column CSV file and return its path."""
    csv_file = tmp_path / "data.csv"
    csv_file.write_text("Gene,log2FC,padj\nPMM0001,1.5,0.01\nPMM0002,-0.8,0.03\n")
    return csv_file


def _write_config(tmp_path: Path, config: dict) -> Path:
    """Serialise *config* to a YAML file inside *tmp_path* and return its path."""
    cfg_file = tmp_path / "paperconfig.yaml"
    cfg_file.write_text(yaml.dump(config))
    return cfg_file


def _make_valid_config(csv_path: Path, overrides: dict = None,
                       experiment_overrides: dict = None) -> dict:
    """Return a minimal valid paperconfig dict that passes all checks.

    The CSV file at *csv_path* must have columns: Gene, log2FC, padj.
    Pass *overrides* as a shallow dict to override keys inside the single
    statistical_analyses entry.
    Pass *experiment_overrides* to override keys inside the experiment entry.
    """
    experiment = {
        "name": "Test coculture experiment",
        "organism": "Prochlorococcus MED4",
        "omics_type": "RNASEQ",
        "test_type": "DESeq2",
        "treatment_type": "coculture",
        "treatment_condition": "Coculture",
        "control_condition": "Axenic",
        "treatment_organism": "Alteromonas",
        "treatment_taxid": "28108",
    }
    if experiment_overrides:
        experiment.update(experiment_overrides)
    analysis = {
        "id": "test_analysis_1",
        "experiment": "exp_coculture",
        "timepoint_hours": None,
        "name_col": "Gene",
        "logfc_col": "log2FC",
        "adjusted_p_value_col": "padj",
    }
    if overrides:
        analysis.update(overrides)
    config = {
        "publication": {
            "papername": "Test Paper 2024",
            "experiments": {
                "exp_coculture": experiment,
            },
            "supplementary_materials": {
                "supp_table_1": {
                    "type": "csv",
                    "filename": str(csv_path),
                    "statistical_analyses": [analysis],
                }
            },
        }
    }
    return config


# ---------------------------------------------------------------------------
# Unit tests — canonical organism validation
# ---------------------------------------------------------------------------

class TestCanonicalOrganism:
    """Validator rejects non-canonical organism names in experiments."""

    def test_non_canonical_organism_in_experiment_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry with an unrecognised organism name must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"organism": "ProchlorococcusXYZ UnknownStrain"},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical organism"

    def test_canonical_organism_in_experiment_is_accepted(self, tmp_path, monkeypatch):
        """All canonical organism names must be accepted without errors."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for organism in sorted(CANONICAL_GENOMIC_ORGANISMS):
            config = _make_valid_config(
                csv,
                experiment_overrides={"organism": organism},
            )
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical organism '{organism}'"

    def test_non_canonical_treatment_organism_is_rejected(self, tmp_path, monkeypatch):
        """A treatment_organism value not in the canonical set must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "treatment_organism": "NonExistentBacterium",
                "treatment_taxid": "99999",
            },
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical treatment_organism"


# ---------------------------------------------------------------------------
# Unit tests — canonical condition_type validation
# ---------------------------------------------------------------------------

class TestCanonicalTreatmentType:
    """Validator rejects non-canonical treatment_type values in experiments."""

    def test_non_canonical_treatment_type_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry with an unknown treatment_type must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_type": "made_up_stress_type"},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical treatment_type"

    def test_canonical_treatment_type_is_accepted(self, tmp_path, monkeypatch):
        """Every canonical treatment_type value must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for ttype in sorted(CANONICAL_CONDITION_TYPES):
            config = _make_valid_config(
                csv,
                experiment_overrides={"treatment_type": ttype},
            )
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical treatment_type '{ttype}'"


# ---------------------------------------------------------------------------
# Unit tests — treatment_type as list
# ---------------------------------------------------------------------------

class TestTreatmentTypeList:
    """Validator accepts treatment_type as a list."""

    def test_treatment_type_list_is_accepted(self, tmp_path, monkeypatch):
        """A list of canonical treatment_type values must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_type": ["coculture", "darkness"]},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass for list treatment_type"

    def test_treatment_type_list_with_invalid_value_rejected(self, tmp_path, monkeypatch):
        """A list containing a non-canonical value must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_type": ["coculture", "fake_stress"]},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for list with non-canonical value"


# ---------------------------------------------------------------------------
# Unit tests — background_factors validation
# ---------------------------------------------------------------------------

class TestBackgroundFactors:
    """Validator accepts and validates background_factors field."""

    def test_background_factors_list_accepted(self, tmp_path, monkeypatch):
        """A list of canonical background_factors values must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "treatment_type": "nitrogen_stress",
                "background_factors": ["axenic", "continuous_light"],
                # Remove coculture-specific fields since treatment_type is not coculture
                "treatment_organism": None,
                "treatment_taxid": None,
            },
        )
        # Remove None-valued keys (treatment_organism, treatment_taxid)
        exp = config["publication"]["experiments"]["exp_coculture"]
        exp.pop("treatment_organism", None)
        exp.pop("treatment_taxid", None)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass for valid background_factors"

    def test_background_factors_with_invalid_value_rejected(self, tmp_path, monkeypatch):
        """A background_factors list with non-canonical value must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "background_factors": ["axenic", "made_up_factor"],
            },
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical background_factors"

    def test_background_factors_optional(self, tmp_path, monkeypatch):
        """Missing background_factors should not cause validation failure."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass without background_factors"


# ---------------------------------------------------------------------------
# Unit tests — canonical test_type validation
# ---------------------------------------------------------------------------

class TestCanonicalTestType:
    """Validator rejects non-canonical test_type values in experiments."""

    def test_non_canonical_test_type_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry with an unknown test_type must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"test_type": "NonExistentMethod"},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical test_type"

    def test_canonical_test_types_are_accepted(self, tmp_path, monkeypatch):
        """Every canonical test_type value must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for ttype in sorted(CANONICAL_TEST_TYPES):
            config = _make_valid_config(
                csv,
                experiment_overrides={"test_type": ttype},
            )
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical test_type '{ttype}'"


# ---------------------------------------------------------------------------
# Unit tests — required fields in statistical_analyses
# ---------------------------------------------------------------------------

class TestRequiredAnalysisFields:
    """Validator enforces id and experiment reference in statistical_analyses."""

    def test_missing_id_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry without 'id' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        # Remove 'id' from the analysis entry
        analysis = config["publication"]["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["id"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'id' is missing"

    def test_missing_experiment_reference_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry without 'experiment' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        analysis = config["publication"]["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["experiment"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'experiment' is missing"

    def test_bad_experiment_reference_is_rejected(self, tmp_path, monkeypatch):
        """An analysis with an experiment reference not in experiments block must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"experiment": "nonexistent_experiment"})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for bad experiment reference"


class TestRequiredExperimentFields:
    """Validator enforces required fields on experiment entries."""

    def test_missing_treatment_condition_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry without 'treatment_condition' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        experiment = config["publication"]["experiments"]["exp_coculture"]
        del experiment["treatment_condition"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'treatment_condition' is missing from experiment"

    def test_empty_treatment_condition_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry with an empty string 'treatment_condition' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_condition": "  "},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when experiment 'treatment_condition' is empty/whitespace"

    def test_missing_organism_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry without 'organism' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        experiment = config["publication"]["experiments"]["exp_coculture"]
        del experiment["organism"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'organism' is missing from experiment"

    def test_missing_omics_type_is_rejected(self, tmp_path, monkeypatch):
        """An experiment entry without 'omics_type' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        experiment = config["publication"]["experiments"]["exp_coculture"]
        del experiment["omics_type"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'omics_type' is missing from experiment"


# ---------------------------------------------------------------------------
# Unit test — fully valid config passes
# ---------------------------------------------------------------------------

class TestValidConfigPasses:
    """A well-formed config with all canonical values must pass validation."""

    def test_valid_config_passes(self, tmp_path, monkeypatch):
        """A minimal but fully correct config should return True."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected a fully valid config to pass validation"

    def test_valid_config_with_stress_experiment_passes(self, tmp_path, monkeypatch):
        """A config with a non-coculture (stress) experiment should pass validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        # Build a config with a phosphorus stress experiment (no treatment_organism)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "name": "Phosphorus depletion experiment",
                "treatment_type": "phosphorus_stress",
                "treatment_condition": "P-depleted medium",
                "control_condition": "P-replete medium",
                # Remove coculture-specific fields
                "treatment_organism": None,
                "treatment_taxid": None,
            },
        )
        # Remove None-valued keys from experiment
        exp = config["publication"]["experiments"]["exp_coculture"]
        exp.pop("treatment_organism", None)
        exp.pop("treatment_taxid", None)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected config with stress experiment to pass validation"

    def test_valid_config_with_timepoint_hours(self, tmp_path, monkeypatch):
        """A config with numeric timepoint_hours on analyses should pass validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"timepoint_hours": 24.0})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected config with numeric timepoint_hours to pass validation"

    def test_valid_config_with_null_timepoint_hours(self, tmp_path, monkeypatch):
        """A config with null timepoint_hours on analyses should pass validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"timepoint_hours": None})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected config with null timepoint_hours to pass validation"

    def test_invalid_timepoint_hours_is_rejected(self, tmp_path, monkeypatch):
        """A config with non-numeric timepoint_hours must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"timepoint_hours": "24h"})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-numeric timepoint_hours"

    def test_missing_experiments_block_warns_but_passes(self, tmp_path, monkeypatch):
        """A publication config without experiments block should warn but pass (cluster-only configs)."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        del config["publication"]["experiments"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass (with warning) when experiments block is missing"
