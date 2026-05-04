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
