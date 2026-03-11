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

from validate_paperconfig import validate, CANONICAL_GENOMIC_ORGANISMS, CANONICAL_CONDITION_TYPES, CANONICAL_TEST_TYPES

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


def _make_valid_config(csv_path: Path, overrides: dict = None) -> dict:
    """Return a minimal valid paperconfig dict that passes all checks.

    The CSV file at *csv_path* must have columns: Gene, log2FC, padj.
    Pass *overrides* as a shallow dict to override keys inside the single
    statistical_analyses entry.
    """
    analysis = {
        "id": "test_analysis_1",
        "name": "test analysis",
        "type": "RNASEQ",
        "test_type": "DESeq2",
        "control_condition": "Axenic",
        "treatment_condition": "Coculture",
        "organism": "Prochlorococcus MED4",
        "name_col": "Gene",
        "logfc_col": "log2FC",
        "adjusted_p_value_col": "padj",
        "treatment_organism": "Alteromonas",
        "treatment_taxid": "28108",
    }
    if overrides:
        analysis.update(overrides)
    config = {
        "publication": {
            "papername": "Test Paper 2024",
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
    """Validator rejects non-canonical organism names."""

    def test_non_canonical_organism_in_analysis_is_rejected(self, tmp_path, monkeypatch):
        """A statistical_analyses entry with an unrecognised organism name must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"organism": "ProchlorococcusXYZ UnknownStrain"})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical organism"

    def test_canonical_organism_in_analysis_is_accepted(self, tmp_path, monkeypatch):
        """All canonical organism names must be accepted without errors."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for organism in sorted(CANONICAL_GENOMIC_ORGANISMS):
            config = _make_valid_config(csv, overrides={"organism": organism})
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical organism '{organism}'"

    def test_non_canonical_treatment_organism_is_rejected(self, tmp_path, monkeypatch):
        """A treatment_organism value not in the canonical set must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            overrides={"treatment_organism": "NonExistentBacterium", "treatment_taxid": "99999"},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical treatment_organism"


# ---------------------------------------------------------------------------
# Unit tests — canonical condition_type validation
# ---------------------------------------------------------------------------

class TestCanonicalConditionType:
    """Validator rejects non-canonical condition_type values in environmental_conditions."""

    def _config_with_env_condition(self, csv_path: Path, condition_type: str) -> dict:
        """Build a config that uses an environmental_treatment_condition_id."""
        config = _make_valid_config(csv_path)
        pub = config["publication"]
        pub["environmental_conditions"] = {
            "treatment_cond": {
                "name": "Test condition",
                "condition_type": condition_type,
            }
        }
        # Replace treatment_organism/taxid with env condition reference
        analysis = pub["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["treatment_organism"]
        del analysis["treatment_taxid"]
        analysis["environmental_treatment_condition_id"] = "treatment_cond"
        return config

    def test_non_canonical_condition_type_is_rejected(self, tmp_path, monkeypatch):
        """An environmental_conditions entry with an unknown condition_type must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = self._config_with_env_condition(csv, "made_up_stress_type")
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical condition_type"

    def test_canonical_condition_type_is_accepted(self, tmp_path, monkeypatch):
        """Every canonical condition_type value must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for ctype in sorted(CANONICAL_CONDITION_TYPES):
            config = self._config_with_env_condition(csv, ctype)
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical condition_type '{ctype}'"


# ---------------------------------------------------------------------------
# Unit tests — canonical test_type validation
# ---------------------------------------------------------------------------

class TestCanonicalTestType:
    """Validator rejects non-canonical test_type values in statistical_analyses."""

    def test_non_canonical_test_type_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry with an unknown test_type must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"test_type": "NonExistentMethod"})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical test_type"

    def test_canonical_test_types_are_accepted(self, tmp_path, monkeypatch):
        """Every canonical test_type value must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        for ttype in sorted(CANONICAL_TEST_TYPES):
            config = _make_valid_config(csv, overrides={"test_type": ttype})
            cfg_file = _write_config(tmp_path, config)
            result = validate(str(cfg_file))
            assert result is True, f"Expected validation to pass for canonical test_type '{ttype}'"


# ---------------------------------------------------------------------------
# Unit tests — required fields in statistical_analyses
# ---------------------------------------------------------------------------

class TestRequiredAnalysisFields:
    """Validator enforces id, type, and treatment_condition in statistical_analyses."""

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

    def test_missing_type_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry without 'type' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        analysis = config["publication"]["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["type"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'type' is missing"

    def test_missing_treatment_condition_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry without 'treatment_condition' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        analysis = config["publication"]["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["treatment_condition"]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'treatment_condition' is missing"

    def test_empty_treatment_condition_is_rejected(self, tmp_path, monkeypatch):
        """An analysis entry with an empty string 'treatment_condition' must fail validation."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv, overrides={"treatment_condition": "  "})
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail when 'treatment_condition' is empty/whitespace"


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

    def test_valid_config_with_env_condition_passes(self, tmp_path, monkeypatch):
        """A config using environmental condition references should pass when all values are canonical."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        # Build a config with an environmental condition instead of treatment organism
        config = _make_valid_config(csv)
        pub = config["publication"]
        pub["environmental_conditions"] = {
            "phosphorus_depletion": {
                "name": "Phosphorus depletion",
                "condition_type": "phosphorus_stress",
            }
        }
        analysis = pub["supplementary_materials"]["supp_table_1"]["statistical_analyses"][0]
        del analysis["treatment_organism"]
        del analysis["treatment_taxid"]
        analysis["environmental_treatment_condition_id"] = "phosphorus_depletion"
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected config with canonical env condition to pass validation"
