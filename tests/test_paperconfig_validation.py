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
VALIDATE_SCRIPT_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "scripts")
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
                "treatment_type": "nitrogen",
                "background_factors": ["axenic", "light"],
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
                "treatment_type": "phosphorus",
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


# ---------------------------------------------------------------------------
# Unit tests — treatment_assembly_accession must exist in cyanobacteria_genomes.csv
# ---------------------------------------------------------------------------

class TestTreatmentAssemblyAccession:
    """Validator rejects treatment_assembly_accession not present in the KG genomes registry.

    The omics adapter uses treatment_assembly_accession to construct the
    Tests_coculture_with edge target as `insdc.gcf:<accession>`. If the
    accession is not in cyanobacteria_genomes.csv, neo4j-admin import silently
    drops the edge (only surfaces in import.report). This check makes the
    mismatch a hard error at paperconfig-write time. Regression test for the
    Kratzl 2024 incident (2026-04-14).
    """

    def test_unknown_treatment_assembly_accession_is_rejected(self, tmp_path, monkeypatch):
        """A paperconfig referencing an accession absent from the genomes CSV must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "treatment_organism": "Alteromonas",
                "treatment_taxid": "28108",
                "treatment_assembly_accession": "GCF_999999999.9",
            },
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, (
            "Expected validation to fail when treatment_assembly_accession is not "
            "present in cyanobacteria_genomes.csv"
        )

    def test_known_treatment_assembly_accession_is_accepted(self, tmp_path, monkeypatch):
        """A paperconfig using a real KG accession (MED4 here) must pass."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "treatment_organism": "Prochlorococcus MED4",
                "treatment_taxid": "59919",
                "treatment_assembly_accession": "GCF_000011465.1",
            },
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, (
            "Expected validation to pass when treatment_assembly_accession matches a "
            "row in cyanobacteria_genomes.csv"
        )

    def test_missing_treatment_assembly_accession_is_accepted(self, tmp_path, monkeypatch):
        """The field is optional — absence must not trigger the new error."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)  # no treatment_assembly_accession
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, (
            "Expected validation to pass when treatment_assembly_accession is absent"
        )


# ---------------------------------------------------------------------------
# Unit tests — free-text column declared as locus_tag warning
# ---------------------------------------------------------------------------

class TestFreeTextAsLocusTagWarning:
    """Validator warns when a free-text-looking column is declared id_type: locus_tag.

    This is the bug where ``{column: Description, id_type: locus_tag}`` caused
    build_gene_id_mapping to tokenise entire sentences into specific_lookup,
    polluting a strain's gene_id map with thousands of junk alt_ids.
    """

    def _write_csv_with_description(self, tmp_path: Path) -> Path:
        csv_file = tmp_path / "data.csv"
        csv_file.write_text(
            "Gene,Description,log2FC,padj\n"
            "PMM0001,ATP synthase GN=atpB OS=Prochlorococcus,1.5,0.01\n"
            "PMM0002,Photosystem II subunit GN=psbF OS=Prochlorococcus,-0.8,0.03\n"
        )
        return csv_file

    def test_description_as_locus_tag_emits_warning(self, tmp_path, monkeypatch, capsys):
        """Declaring 'Description' as id_type: locus_tag must emit a warning (still passes)."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = self._write_csv_with_description(tmp_path)
        config = _make_valid_config(csv)
        supp = config["publication"]["supplementary_materials"]["supp_table_1"]
        supp["id_columns"] = [
            {"column": "Description", "id_type": "locus_tag"},
        ]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        # Validator should still pass (warning, not error).
        assert result is True, "Expected validation to pass (warning, not error)"
        out = capsys.readouterr().out
        assert "looks like free-text description" in out, (
            f"Expected free-text warning not found in validator output.\n"
            f"Output:\n{out}"
        )
        assert "Description" in out

    def test_description_as_gene_name_does_not_warn(self, tmp_path, monkeypatch, capsys):
        """Declaring 'Description' as id_type: gene_name (Tier 3) does not trip the warning."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = self._write_csv_with_description(tmp_path)
        config = _make_valid_config(csv)
        supp = config["publication"]["supplementary_materials"]["supp_table_1"]
        supp["id_columns"] = [
            {"column": "Description", "id_type": "gene_name"},
        ]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True
        out = capsys.readouterr().out
        assert "looks like free-text description" not in out, (
            "Warning should not fire for non-locus_tag id_type on free-text column"
        )

    def test_normal_column_as_locus_tag_does_not_warn(self, tmp_path, monkeypatch, capsys):
        """A genuinely id-shaped column declared as locus_tag does not warn."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = self._write_csv_with_description(tmp_path)
        config = _make_valid_config(csv)
        supp = config["publication"]["supplementary_materials"]["supp_table_1"]
        supp["id_columns"] = [
            {"column": "Gene", "id_type": "locus_tag"},
        ]
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True
        out = capsys.readouterr().out
        assert "looks like free-text description" not in out


# ---------------------------------------------------------------------------
# Unit tests — DOI override validation
# ---------------------------------------------------------------------------

class TestDoiOverride:
    """Validator accepts/rejects optional publication.doi field."""

    def test_valid_doi_passes(self, tmp_path, monkeypatch):
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "10.1186/2046-9063-8-7"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is True

    def test_missing_doi_passes(self, tmp_path, monkeypatch):
        """doi is optional — omitting it is fine."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        assert "doi" not in config["publication"]
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is True

    def test_invalid_doi_rejected(self, tmp_path, monkeypatch):
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "not-a-doi"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is False

    def test_doi_url_rejected(self, tmp_path, monkeypatch):
        """Full URL form should be rejected — we want the bare DOI."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        config["publication"]["doi"] = "https://doi.org/10.1186/2046-9063-8-7"
        cfg_file = _write_config(tmp_path, config)
        assert validate(str(cfg_file)) is False


# ---------------------------------------------------------------------------
# validate_paperconfig_content — pure function entry
# ---------------------------------------------------------------------------

def test_validate_paperconfig_content_is_importable():
    """The refactored pure function must be exported from validate_paperconfig."""
    from validate_paperconfig import validate_paperconfig_content

    errors, warnings = validate_paperconfig_content({}, "/tmp/fake.yaml")
    assert isinstance(errors, list)
    assert isinstance(warnings, list)


def test_validate_paperconfig_content_accepts_minimal_valid_config(tmp_path):
    """A minimal well-formed config returns empty errors (warnings OK)."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "paper.pdf"
    pdf.write_bytes(b"%PDF-1.4\n")  # placeholder; validator only checks existence
    config = {
        "publication": {
            "papername": "Test 2026",
            "papermainpdf": str(pdf),
            "experiments": {
                "exp1": {
                    "name": "MED4 N-starvation",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                    "treatment_condition": "N-limited",
                    "control_condition": "N-replete",
                },
            },
            "supplementary_materials": {
                "t1": {
                    "type": "csv",
                    "filename": str(csv_file),
                    "statistical_analyses": [
                        {
                            "id": "an1",
                            "experiment": "exp1",
                            "name_col": "Gene",
                            "logfc_col": "log2FC",
                            "adjusted_p_value_col": "padj",
                            "timepoint": "unknown",
                            "timepoint_hours": None,
                            "growth_phase": "exponential",
                        },
                    ],
                },
            },
        },
    }
    cfg_path = _write_config(tmp_path, config)
    errors, _warnings = validate_paperconfig_content(config, str(cfg_path))
    assert errors == [], f"Expected no errors; got: {errors}"


# ---------------------------------------------------------------------------
# compartment vocabulary check (Task 7)
# ---------------------------------------------------------------------------

def test_validate_rejects_unknown_compartment(tmp_path):
    """Experiment.compartment must be in COMPARTMENTS."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X",
            "papermainpdf": str(pdf),
            "experiments": {
                "exp1": {
                    "name": "e", "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ", "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                    "treatment_condition": "A", "control_condition": "B",
                    "compartment": "nucleus",   # ← not in vocabulary
                },
            },
            "supplementary_materials": {
                "t1": {"type": "csv", "filename": str(csv_file),
                       "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                                 "name_col": "Gene", "logfc_col": "log2FC",
                                                 "adjusted_p_value_col": "padj",
                                                 "timepoint": "unknown", "timepoint_hours": None,
                                                 "growth_phase": "exponential"}]},
            },
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("compartment" in e and "nucleus" in e for e in errors), errors


def test_validate_accepts_whole_cell_compartment(tmp_path):
    """'whole_cell' is the default; explicit 'whole_cell' must be accepted."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp1": {
                "name": "e", "organism": "Prochlorococcus MED4",
                "omics_type": "RNASEQ", "test_type": "DESeq2",
                "treatment_type": ["nitrogen"],
                "treatment_condition": "A", "control_condition": "B",
                "compartment": "whole_cell",
            }},
            "supplementary_materials": {"t1": {
                "type": "csv", "filename": str(csv_file),
                "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                          "name_col": "Gene", "logfc_col": "log2FC",
                                          "adjusted_p_value_col": "padj",
                                          "timepoint": "unknown", "timepoint_hours": None,
                                          "growth_phase": "exponential"}],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert not any("compartment" in e for e in errors), errors


def test_validate_accepts_paired_rnaseq_proteome_omics_type(tmp_path):
    """PAIRED_RNASEQ_PROTEOME is a canonical omics_type after this slice lands."""
    from validate_paperconfig import validate_paperconfig_content

    csv_file = _write_minimal_csv(tmp_path)
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp1": {
                "name": "e", "organism": "Prochlorococcus MED4",
                "omics_type": "PAIRED_RNASEQ_PROTEOME",
                "test_type": "DESeq2", "treatment_type": ["diel"],
                "treatment_condition": "A", "control_condition": "B",
            }},
            "supplementary_materials": {"t1": {
                "type": "csv", "filename": str(csv_file),
                "statistical_analyses": [{"id": "a", "experiment": "exp1",
                                          "name_col": "Gene", "logfc_col": "log2FC",
                                          "adjusted_p_value_col": "padj",
                                          "timepoint": "unknown", "timepoint_hours": None,
                                          "growth_phase": "exponential"}],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert not any("omics_type" in e for e in errors), errors


def test_validate_relaxes_de_fields_for_derived_metrics_only_experiment(tmp_path):
    """Experiments whose only supplementary entries are derived_metrics_table may
    omit control_condition and test_type (warn, don't error)."""
    from validate_paperconfig import validate_paperconfig_content

    # A minimal CSV with a Y/blank boolean column, referenced by the DM entry.
    csv_file = tmp_path / "dm.csv"
    csv_file.write_text("locus_tag,flag_col\nPMN2A_RS00015,Y\nPMN2A_RS00020,\n")
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    config = {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp_dm_only": {
                "name": "e", "organism": "Prochlorococcus NATL2A",
                "omics_type": "RNASEQ",
                "treatment_type": ["darkness"],
                "treatment_condition": "A",
                # ← control_condition + test_type deliberately omitted
            }},
            "supplementary_materials": {"dm1": {
                "type": "derived_metrics_table",
                "filename": str(csv_file),
                "organism": "Prochlorococcus NATL2A",
                "experiment": "exp_dm_only",
                "name_col": "locus_tag",
                "id_columns": [{"column": "locus_tag", "id_type": "locus_tag_ncbi"}],
                "metrics": [{
                    "metric_type": "periodic_in_axenic_LD",
                    "value_kind": "boolean",
                    "value_col": "flag_col",
                    "true_tokens": ["Y"],
                }],
            }},
        },
    }
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert not any("control_condition" in e for e in errors), errors
    assert not any("test_type" in e for e in errors), errors


# ---------------------------------------------------------------------------
# derived_metrics_table dispatch (Task 8) — option C: per-metric metadata
# declared inline on paperconfig; KNOWN_METRIC_TYPES only pins value_kind.
# ---------------------------------------------------------------------------

def _dm_boolean_csv(tmp_path: Path) -> Path:
    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,Periodic in axenic L D cultures\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,\n"
        "PMN2A_RS00025,NA\n"
    )
    return csv


def _dm_wrapper_config(tmp_path: Path, dm_entry: dict) -> dict:
    csv_file = _write_minimal_csv(tmp_path)  # unused but keeps structure
    pdf = tmp_path / "p.pdf"; pdf.write_bytes(b"%PDF-1.4\n")
    return {
        "publication": {
            "papername": "X", "papermainpdf": str(pdf),
            "experiments": {"exp_dm_only": {
                "name": "e", "organism": "Prochlorococcus NATL2A",
                "omics_type": "RNASEQ",
                "treatment_type": ["darkness"],
                "treatment_condition": "ED",
            }},
            "supplementary_materials": {"dm1": dm_entry},
        },
    }


def test_validate_boolean_dm_entry_accepts_shape(tmp_path):
    """Minimal well-formed boolean DM entry — no rankable/has_p_value declared
    (adapter sets them to 'false' at ingest since boolean is definitionally
    non-rankable and has no p-value). field_description is required."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table",
        "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "blank_policy": "skip",
            "field_description": "RAIN periodicity FDR<0.05 in axenic L:D",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert errors == [], errors


def test_validate_dm_requires_field_description(tmp_path):
    """field_description is required on every metric (all value_kinds) —
    free-text human-readable explanation surfaced by downstream MCP tools."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            # ← no field_description
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("field_description" in e for e in errors), errors


def test_validate_dm_rejects_empty_field_description(tmp_path):
    """field_description must be non-empty (whitespace-only is also rejected)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "field_description": "   ",   # whitespace-only
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("field_description" in e for e in errors), errors


def test_validate_boolean_dm_rejects_unknown_token(tmp_path):
    """CSV dry-run finds an unclassified token → hard error."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,flag\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,MAYBE\n"   # ← unknown token
    )
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean", "value_col": "flag",
            "true_tokens": ["Y"], "blank_policy": "skip",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("MAYBE" in e for e in errors), errors


def test_validate_boolean_dm_rejects_missing_true_tokens(tmp_path):
    """true_tokens is required (non-empty list) for value_kind=boolean."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            # ← no true_tokens
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("true_tokens" in e for e in errors), errors


def test_validate_boolean_dm_forbids_rankable_and_has_p_value(tmp_path):
    """rankable and has_p_value are forbidden on boolean entries — adapter sets
    them to "false" at ingest since they're definitionally "false" for boolean.
    Declaring them on the paperconfig (even as "false") is redundant and
    rejected to prevent copy-paste noise."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "rankable": "false",       # ← forbidden even when "false"
            "has_p_value": "false",    # ← forbidden even when "false"
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("rankable" in e for e in errors), errors
    assert any("has_p_value" in e for e in errors), errors


def test_validate_categorical_dm_forbids_rankable_and_has_p_value(tmp_path):
    """Same forbid rule as boolean — adapter sets at ingest."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text("NCBI ID_2,darkness_cluster\nPMN2A_RS00015,a\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            "allowed_categories": ["a"],
            "rankable": "false",       # ← forbidden
            "has_p_value": "false",    # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("rankable" in e for e in errors), errors
    assert any("has_p_value" in e for e in errors), errors


def test_validate_boolean_dm_rejects_forbidden_numeric_fields(tmp_path):
    """unit / p_value_col / p_value_threshold are forbidden on boolean entries."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
            "unit": "h",                  # ← forbidden
            "p_value_threshold": 0.05,    # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("unit" in e for e in errors) and any("p_value_threshold" in e for e in errors), errors


def test_validate_categorical_dm_requires_allowed_categories(tmp_path):
    """value_kind=categorical MUST declare allowed_categories inline (non-empty list)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text("NCBI ID_2,darkness_cluster\nPMN2A_RS00015,darkness_axenic+darkness_coculture\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            # ← no allowed_categories
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("allowed_categories" in e for e in errors), errors


def test_validate_categorical_dm_warns_on_out_of_set_values(tmp_path):
    """CSV dry-run warns if a value_col cell lies outside the paperconfig's
    declared allowed_categories. (Hard error at ingest is the adapter's job
    in Plan 2.)"""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "s5.csv"
    csv.write_text(
        "NCBI ID_2,darkness_cluster\n"
        "PMN2A_RS00015,darkness_axenic+darkness_coculture\n"
        "PMN2A_RS00020,totally_new_category\n"   # ← not in allowed_categories
    )
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "darkness_survival_class",
            "value_kind": "categorical",
            "value_col": "darkness_cluster",
            "allowed_categories": [
                "darkness_axenic+darkness_coculture",
                "darkness_coculture+unique_coculture",
                "darkness_axenic+unique_axenic",
            ],
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert any("totally_new_category" in w for w in warnings), warnings


def test_validate_numeric_dm_requires_value_col_and_rankable(tmp_path):
    """value_kind=numeric requires value_col + rankable + has_p_value."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser_s1.csv"
    csv.write_text("locus_tag,Fourier,Fourier FDR\nPMM0001,0.82,0.01\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "fourier_score",
            "value_kind": "numeric",
            # ← value_col, rankable, has_p_value all missing
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    joined = " | ".join(errors)
    assert "value_col" in joined and "rankable" in joined and "has_p_value" in joined, errors


def test_validate_numeric_dm_has_p_value_true_requires_threshold(tmp_path):
    """When has_p_value='true', p_value_threshold is required and at least one
    of p_value_col / adjusted_p_value_col must be present."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser.csv"
    csv.write_text("locus_tag,Fourier\nPMM0001,0.82\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "fourier_score",
            "value_kind": "numeric",
            "value_col": "Fourier",
            "rankable": "true",
            "has_p_value": "true",
            # ← p_value_threshold and p_value_col/adjusted_p_value_col missing
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("p_value_threshold" in e for e in errors), errors
    assert any("p_value_col" in e or "adjusted_p_value_col" in e for e in errors), errors


def test_validate_numeric_dm_has_p_value_false_forbids_threshold(tmp_path):
    """When has_p_value='false', p_value_threshold is forbidden."""
    from validate_paperconfig import validate_paperconfig_content

    csv = tmp_path / "zinser.csv"
    csv.write_text("locus_tag,Peak\nPMM0001,6.0\n")
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp_dm_only",
        "name_col": "locus_tag",
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag"}],
        "metrics": [{
            "metric_type": "peak_time_h",
            "value_kind": "numeric",
            "value_col": "Peak",
            "unit": "h",
            "rankable": "false",
            "has_p_value": "false",
            "p_value_threshold": 0.05,  # ← forbidden
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    config["publication"]["experiments"]["exp_dm_only"]["organism"] = "Prochlorococcus MED4"
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("p_value_threshold" in e for e in errors), errors


def test_validate_rejects_value_kind_mismatch_against_registry(tmp_path):
    """KNOWN_METRIC_TYPES pins metric_type → value_kind. A paperconfig that
    re-declares a known metric_type under a different value_kind is a hard
    error (prevents silent edge-type changes across papers)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",   # registry: boolean
            "value_kind": "numeric",                    # ← mismatch
            "value_col": "Periodic in axenic L D cultures",
            "rankable": "true",
            "has_p_value": "false",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("periodic_in_axenic_LD" in e and "value_kind" in e for e in errors), errors


def test_validate_novel_metric_type_warns_not_errors(tmp_path):
    """A metric_type absent from KNOWN_METRIC_TYPES is accepted with a warning —
    authors may introduce new names; the registry grows slowly."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "a_brand_new_per_paper_flag",   # not in registry
            "value_kind": "boolean",
            "value_col": "Periodic in axenic L D cultures",
            "true_tokens": ["Y"],
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, warnings = validate_paperconfig_content(config, str(cfg))
    assert not any("a_brand_new_per_paper_flag" in e for e in errors), errors
    assert any("a_brand_new_per_paper_flag" in w and ("novel" in w.lower() or "not in" in w.lower())
               for w in warnings), warnings


def test_validate_rejects_unknown_value_kind(tmp_path):
    """value_kind must be a member of VALUE_KINDS (numeric|boolean|categorical)."""
    from validate_paperconfig import validate_paperconfig_content

    csv = _dm_boolean_csv(tmp_path)
    dm_entry = {
        "type": "derived_metrics_table", "filename": str(csv),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp_dm_only",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "ordinal",  # ← not a VALUE_KINDS member
            "value_col": "Periodic in axenic L D cultures",
        }],
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    assert any("value_kind" in e and "ordinal" in e for e in errors), errors


def test_validate_rejects_missing_top_level_fields(tmp_path):
    """Top-level required fields: filename, organism, experiment, name_col, metrics."""
    from validate_paperconfig import validate_paperconfig_content

    dm_entry = {
        "type": "derived_metrics_table",
        # everything else missing
    }
    config = _dm_wrapper_config(tmp_path, dm_entry)
    cfg = _write_config(tmp_path, config)
    errors, _ = validate_paperconfig_content(config, str(cfg))
    joined = " | ".join(errors)
    for req in ("filename", "organism", "experiment", "name_col", "metrics"):
        assert req in joined, f"'{req}' should be flagged as missing; errors={errors}"
