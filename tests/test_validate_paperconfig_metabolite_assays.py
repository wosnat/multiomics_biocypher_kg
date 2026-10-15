"""Tests for validate_paperconfig _validate_metabolite_assays_entry helper."""
from scripts.validate_paperconfig import _validate_metabolite_assays_entry


def _experiments_with_org(comp="whole_cell", org="Prochlorococcus MIT9303"):
    return {"exp1": {"organism": org, "compartment": comp}}


def _good_entry(filename="data/example.csv"):
    return {
        "type": "metabolite_assays_table",
        "filename": filename,
        "experiment": "exp1",
        "organism": "Prochlorococcus MIT9303",
        "name_col": "compound",
        "assays": [
            {
                "metric_type": "cellular_concentration",
                "value_kind": "numeric",
                "field_description": "test description",
                "sample_columns": [
                    {"condition_label": "control", "replicate_columns": ["c1", "c2"]}
                ],
            }
        ],
    }


def test_validate_metabolite_assays_table_happy_path(tmp_path):
    entry = _good_entry()
    # avoid file-not-found warning by writing the file
    csv = tmp_path / "example.csv"
    csv.write_text("compound,c1,c2\n")
    entry["filename"] = str(csv)
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "entry_key", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert errors == []


def test_validate_metabolite_assays_table_missing_required_field():
    entry = {"type": "metabolite_assays_table", "filename": "x.csv"}
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml", {}, set(), errors, warnings,
    )
    joined = " | ".join(errors)
    for needed in ("experiment", "organism", "name_col", "assays"):
        assert needed in joined, f"missing-field message did not flag {needed}"


def test_validate_metabolite_assays_table_unknown_experiment():
    entry = _good_entry()
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        {"other_exp": {"organism": "x", "compartment": "whole_cell"}},
        {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("exp1" in e and "not found" in e for e in errors)


def test_validate_metabolite_assays_table_compartment_not_in_vocab():
    entry = _good_entry()
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(comp="bogus_compartment"),
        {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("bogus_compartment" in e for e in errors)


def test_validate_metabolite_assays_table_organism_mismatch():
    entry = _good_entry()
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(org="Prochlorococcus MIT9313"),
        {"Prochlorococcus MIT9303", "Prochlorococcus MIT9313"}, errors, warnings,
    )
    assert any("does not match" in e for e in errors)


def test_validate_metabolite_assays_table_boolean_requires_flag_column():
    entry = _good_entry()
    entry["assays"][0]["value_kind"] = "boolean"
    # leave replicate_columns; missing flag_column AND flag_true_value should error
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("flag_column" in e for e in errors)
    assert any("flag_true_value" in e for e in errors)


def test_validate_metabolite_assays_table_duplicate_metric_type():
    entry = _good_entry()
    # Add a second assay with the same metric_type
    entry["assays"].append({
        "metric_type": "cellular_concentration",  # duplicate
        "value_kind": "numeric",
        "field_description": "another",
        "sample_columns": [{"replicate_columns": ["c3"]}],
    })
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("duplicate metric_type" in e for e in errors)


def test_validate_metabolite_assays_table_aggregation_method_enum():
    entry = _good_entry()
    entry["aggregation_method"] = "weird_method"
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("aggregation_method" in e and "weird_method" in e for e in errors)


def test_validate_metabolite_assays_table_numeric_requires_replicate_columns():
    entry = _good_entry()
    entry["assays"][0]["sample_columns"] = [{"condition_label": "x"}]  # no replicate_columns
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("replicate_columns" in e for e in errors)


def test_validate_metabolite_assays_embedded_format_requires_total_replicates(tmp_path):
    """When cell_format=embedded_mean_sd_n, numeric assays must declare
    total_replicates so the adapter can compute n_replicates / detection_status."""
    entry = _good_entry()
    csv = tmp_path / "example.csv"
    csv.write_text("compound,c1\n")
    entry["filename"] = str(csv)
    entry["cell_format"] = "embedded_mean_sd_n"
    entry["assays"][0]["sample_columns"] = [{"condition_label": "x", "replicate_columns": ["c1"]}]
    # total_replicates intentionally omitted — should error
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("total_replicates" in e for e in errors), (
        f"validator should require total_replicates for embedded_mean_sd_n; got: {errors}"
    )


def test_validate_metabolite_assays_embedded_format_with_total_replicates_ok(tmp_path):
    entry = _good_entry()
    csv = tmp_path / "example.csv"
    csv.write_text("compound,c1\n")
    entry["filename"] = str(csv)
    entry["cell_format"] = "embedded_mean_sd_n"
    entry["assays"][0]["sample_columns"] = [{"condition_label": "x", "replicate_columns": ["c1"]}]
    entry["assays"][0]["total_replicates"] = 3
    errors = []
    warnings = []
    _validate_metabolite_assays_entry(
        "e", entry, "/tmp/cfg.yaml",
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert errors == [], f"valid embedded entry should pass; got: {errors}"


# --- alias-file validation (metabolite identity, spec 3.6) ---

def _entry_with_aliases(tmp_path, alias_text):
    csv = tmp_path / "example.csv"
    csv.write_text("compound,c1\n")
    (tmp_path / "aliases.yaml").write_text(alias_text)
    entry = _good_entry(str(csv))
    entry["aliases_file"] = "aliases.yaml"
    return entry


def _run_aliases(tmp_path, alias_text):
    entry = _entry_with_aliases(tmp_path, alias_text)
    errors, warnings = [], []
    _validate_metabolite_assays_entry(
        "e", entry, str(tmp_path / "paperconfig.yaml"),
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    return errors, warnings


def test_alias_file_structured_form_ok(tmp_path):
    errors, warnings = _run_aliases(
        tmp_path,
        "alanine:\n  id: kegg.compound:C00041\n  stereo_assignment: as_reported\n  reason: paper says L\n",
    )
    assert errors == []
    assert not any("legacy" in w for w in warnings)


def test_alias_file_bad_stereo_assignment_errors(tmp_path):
    errors, _ = _run_aliases(
        tmp_path, "alanine:\n  id: kegg.compound:C00041\n  stereo_assignment: maybe\n"
    )
    assert any("stereo_assignment" in e for e in errors)


def test_alias_file_same_name_key_different_ids_errors(tmp_path):
    errors, _ = _run_aliases(
        tmp_path, "Glutathione: kegg.compound:C00051\nglutathione: kegg.compound:C00127\n"
    )
    assert any("name_key" in e for e in errors)


def test_alias_file_formula_case_distinct_ok(tmp_path):
    errors, _ = _run_aliases(
        tmp_path, "CO2: kegg.compound:C00011\nCo2: kegg.compound:C00007\n"
    )
    assert errors == []


def test_alias_file_legacy_string_warns(tmp_path):
    errors, warnings = _run_aliases(tmp_path, "alanine: kegg.compound:C00041\n")
    assert errors == []
    assert any("legacy string form" in w for w in warnings)


def test_alias_file_missing_is_error(tmp_path):
    entry = _entry_with_aliases(tmp_path, "")
    entry["aliases_file"] = "nope.yaml"
    errors, warnings = [], []
    _validate_metabolite_assays_entry(
        "e", entry, str(tmp_path / "paperconfig.yaml"),
        _experiments_with_org(), {"Prochlorococcus MIT9303"}, errors, warnings,
    )
    assert any("aliases_file not found" in e for e in errors)
