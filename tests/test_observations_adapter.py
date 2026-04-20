"""Tests for ObservationsAdapter (Plan 2)."""
import os
from pathlib import Path

import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.observations_adapter import (
    ObservationsAdapter,
    MultiObservationsAdapter,
    _clean_str,
    _make_derived_metric_id,
    _resolve_csv_path,
    _parse_boolean_cell,
)


def test_clean_str_replaces_single_quote_and_pipe():
    assert _clean_str("it's | great") == "it^s , great"


def test_clean_str_handles_none():
    assert _clean_str(None) == ""


def test_clean_str_passes_non_string_through():
    # ints stringify but don't go through replacement
    assert _clean_str(42) == "42"


def test_make_derived_metric_id_uses_doi_short():
    dm_id = _make_derived_metric_id(
        doi="10.1128/mSystems.00040-18",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:mSystems.00040-18:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_make_derived_metric_id_falls_back_to_paper_slug():
    dm_id = _make_derived_metric_id(
        doi="",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:biller_2018:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_resolve_csv_path_prefers_resolved(tmp_path):
    src = tmp_path / "table.csv"
    resolved = tmp_path / "table_resolved.csv"
    src.write_text("a,b\n1,2\n")
    resolved.write_text("a,b,resolved_locus_tag\n1,2,PMM0001\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == resolved
    assert used_resolved is True


def test_resolve_csv_path_falls_back_to_original(tmp_path):
    src = tmp_path / "table.csv"
    src.write_text("a,b\n1,2\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == src
    assert used_resolved is False


def _bp(value, blank_policy="skip"):
    return _parse_boolean_cell(
        value,
        true_tokens=["Y", "yes"],
        false_tokens=["N", "no"],
        skip_tokens=["NA", "n/a"],
        blank_policy=blank_policy,
    )


def test_boolean_true_token():
    assert _bp("Y") == "true"
    assert _bp("yes") == "true"


def test_boolean_false_token():
    assert _bp("N") == "false"
    assert _bp("no") == "false"


def test_boolean_skip_token_returns_none():
    assert _bp("NA") is None
    assert _bp("n/a") is None


def test_boolean_blank_with_skip_policy():
    assert _bp("", blank_policy="skip") is None
    assert _bp(None, blank_policy="skip") is None
    assert _bp(float("nan"), blank_policy="skip") is None


def test_boolean_blank_with_true_policy():
    assert _bp("", blank_policy="true") == "true"
    assert _bp(None, blank_policy="true") == "true"


def test_boolean_blank_with_false_policy():
    assert _bp("", blank_policy="false") == "false"
    assert _bp(float("nan"), blank_policy="false") == "false"


def test_boolean_unknown_token_raises():
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp("maybe")


def test_boolean_invalid_blank_policy_raises():
    with pytest.raises(ValueError, match="Invalid blank_policy"):
        _parse_boolean_cell("", true_tokens=["Y"], false_tokens=[], skip_tokens=[], blank_policy="whatever")


def test_boolean_whitespace_treated_as_blank():
    # Exercises the strip()-then-empty branch (distinct from None/NaN/"" pre-strip)
    assert _bp("   ") is None
    # With blank_policy="false", whitespace is blank → "false"
    assert _bp("   ", blank_policy="false") == "false"


def test_boolean_pd_na_is_blank():
    # pd.NA (pandas nullable scalar, distinct from float NaN) routes to blank_policy
    assert _bp(pd.NA, blank_policy="skip") is None


def test_boolean_numeric_cells_raise():
    # int cells from CSV reading must NOT be silently coerced — contract per docstring
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp(1)
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp(0)


def test_boolean_case_sensitive():
    # "y" (lowercase) is NOT in true_tokens=["Y", "yes"] — exact-string match, no case folding
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp("y")


def test_denormalized_fields_from_experiment():
    """_denormalized_fields produces the exact 9-field dict DerivedMetric.* needs."""
    # Arrange: a minimal paperconfig + ObservationsAdapter instance
    # Easiest: build the adapter from an in-memory config dict
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "name": "NATL2A extended darkness",
        "organism": "Prochlorococcus NATL2A",
        "omics_type": "RNASEQ",
        "treatment_type": ["darkness"],
        "background_factors": ["axenic", "diel"],
        "treatment_condition": "Extended darkness",
        "light_condition": "continuous darkness",
        "experimental_context": "Axenic NATL2A in Pro99 at 24C",
        "compartment": "whole_cell",
    }
    # Bypass __init__'s file I/O by constructing manually
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1128/mSystems.00040-18"
    fields = adapter._denormalized_fields(exp)

    assert fields == {
        "organism_name": "Prochlorococcus NATL2A",
        "publication_doi": "10.1128/mSystems.00040-18",
        "compartment": "whole_cell",
        "omics_type": "RNASEQ",
        "treatment_type": ["darkness"],
        "background_factors": ["axenic", "diel"],
        "treatment": "Extended darkness",
        "light_condition": "continuous darkness",
        "experimental_context": "Axenic NATL2A in Pro99 at 24C",
    }


def test_denormalized_fields_compartment_defaults():
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "name": "x", "organism": "Prochlorococcus MED4",
        "omics_type": "RNASEQ", "treatment_type": [],
        "background_factors": [], "treatment_condition": "",
        "light_condition": "", "experimental_context": "",
        # no compartment key
    }
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1234/x"
    fields = adapter._denormalized_fields(exp)
    assert fields["compartment"] == "whole_cell"


def test_denormalized_fields_normalizes_scalar_list_fields():
    """paperconfig may accidentally provide a string where list is expected."""
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "organism": "P. MED4", "omics_type": "RNASEQ",
        "treatment_type": "nitrogen",  # scalar, should become ["nitrogen"]
        "background_factors": "axenic",  # scalar → list
        "treatment_condition": "", "light_condition": "",
        "experimental_context": "", "compartment": "whole_cell",
    }
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1234/x"
    fields = adapter._denormalized_fields(exp)
    assert fields["treatment_type"] == ["nitrogen"]
    assert fields["background_factors"] == ["axenic"]


def _write_s4a_like_paperconfig(tmp_path):
    """Minimal in-memory paperconfig mimicking Biller 2018 S4A-axenic shape."""
    csv_path = tmp_path / "s4a_small.csv"
    # 3 rows: 2 Y + 1 blank for axenic_LD; 1 Y + 2 blank for extended_darkness
    csv_path.write_text(
        'NCBI ID_2,"Periodic in axenic, L:D cultures","Periodic in axenic, extended darkness cultures"\n'
        "PMN2A_RS00015,Y,\n"
        "PMN2A_RS00020,Y,Y\n"
        "PMN2A_RS00025,,\n"
    )
    resolved_path = tmp_path / "s4a_small_resolved.csv"
    resolved_path.write_text(
        'NCBI ID_2,"Periodic in axenic, L:D cultures","Periodic in axenic, extended darkness cultures",resolved_locus_tag,resolution_method\n'
        "PMN2A_RS00015,Y,,PMN2A_1328,tier1:NCBI ID_2\n"
        "PMN2A_RS00020,Y,Y,PMN2A_1329,tier1:NCBI ID_2\n"
        "PMN2A_RS00025,,,,unresolved\n"
    )
    config = {
        "publication": {
            "papername": "Biller 2018",
            "doi": "10.1128/mSystems.00040-18",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {
                "axenic_rnaseq": {
                    "name": "NATL2A axenic",
                    "organism": "Prochlorococcus NATL2A",
                    "omics_type": "RNASEQ",
                    "treatment_type": ["darkness"],
                    "background_factors": ["axenic", "diel"],
                    "treatment_condition": "Extended darkness",
                    "light_condition": "continuous darkness",
                    "experimental_context": "Axenic NATL2A in Pro99",
                },
            },
            "supplementary_materials": {
                "s4a_axenic": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus NATL2A",
                    "experiment": "axenic_rnaseq",
                    "name_col": "NCBI ID_2",
                    "metrics": [
                        {
                            "metric_type": "periodic_in_axenic_LD",
                            "value_kind": "boolean",
                            "value_col": "Periodic in axenic, L:D cultures",
                            "true_tokens": ["Y"],
                            "false_tokens": [],
                            "skip_tokens": ["NA", "N/A"],
                            "blank_policy": "skip",
                            "field_description": "RAIN FDR<0.05 axenic L:D",
                        },
                        {
                            "metric_type": "periodic_in_axenic_extended_darkness",
                            "value_kind": "boolean",
                            "value_col": "Periodic in axenic, extended darkness cultures",
                            "true_tokens": ["Y"],
                            "false_tokens": [],
                            "skip_tokens": ["NA", "N/A"],
                            "blank_policy": "skip",
                            "field_description": "RAIN FDR<0.05 axenic extended darkness",
                        },
                    ],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    return str(pc_path)


def test_get_nodes_emits_one_dm_per_boolean_metric(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()

    dm_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "derived_metric"]
    assert len(dm_nodes) == 2

    ids = {nid for nid, _, _ in dm_nodes}
    assert "derived_metric:mSystems.00040-18:s4a_axenic:periodic_in_axenic_LD" in ids
    assert "derived_metric:mSystems.00040-18:s4a_axenic:periodic_in_axenic_extended_darkness" in ids


def test_boolean_dm_node_has_expected_props(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()
    dm = next(
        props for nid, lbl, props in nodes
        if lbl == "derived_metric" and nid.endswith("periodic_in_axenic_LD")
    )
    assert dm["metric_type"] == "periodic_in_axenic_LD"
    assert dm["value_kind"] == "boolean"
    assert dm["rankable"] == "false"
    assert dm["has_p_value"] == "false"
    assert dm["allowed_categories"] == []
    assert dm["unit"] == ""
    assert dm["field_description"] == "RAIN FDR<0.05 axenic L:D"
    # Denormalized from parent Experiment
    assert dm["organism_name"] == "Prochlorococcus NATL2A"
    assert dm["compartment"] == "whole_cell"  # default
    assert dm["omics_type"] == "RNASEQ"
    assert dm["treatment_type"] == ["darkness"]
    assert dm["background_factors"] == ["axenic", "diel"]
    assert dm["treatment"] == "Extended darkness"
    assert dm["light_condition"] == "continuous darkness"
    assert dm["experimental_context"] == "Axenic NATL2A in Pro99"
    assert dm["publication_doi"] == "10.1128/mSystems.00040-18"
    # experiment_id is the raw-doi + exp_key concatenation
    assert dm["experiment_id"] == "10.1128/mSystems.00040-18_axenic_rnaseq"
    assert dm["name"].startswith("periodic_in_axenic_LD")  # default name


def test_categorical_allowed_categories_scalar_string_wraps_to_list(tmp_path):
    """A paperconfig authoring mistake (scalar instead of list) must not
    silently iterate per-character — it should wrap to a 1-element list."""
    csv_path = tmp_path / "cat.csv"
    csv_path.write_text("locus_tag,cluster\nPMM0001,high_to_low\n")
    config = {
        "publication": {
            "papername": "Test", "doi": "10.9999/t",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "custom_cluster",
                        "value_kind": "categorical",
                        "value_col": "cluster",
                        "allowed_categories": "high_to_low",  # scalar, not list
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=str(pc_path))
    nodes = adapter.get_nodes()
    dm_props = next(p for _, lbl, p in nodes if lbl == "derived_metric")
    # Scalar string wrapped to 1-element list — NOT per-character decomposition
    assert dm_props["allowed_categories"] == ["high_to_low"]


def test_numeric_p_value_threshold_empty_string_is_noop(tmp_path, caplog):
    """YAML authors sometimes quote empty values; float('') would crash. Guard."""
    csv_path = tmp_path / "num.csv"
    csv_path.write_text("locus_tag,v\nPMM0001,0.5\n")
    config = {
        "publication": {
            "papername": "Test", "doi": "10.9999/t",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "some_metric",
                        "value_kind": "numeric",
                        "value_col": "v",
                        "rankable": "true",
                        "has_p_value": "true",
                        "p_value_threshold": "",  # blank string — must not crash
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=str(pc_path))
    nodes = adapter.get_nodes()  # must NOT raise
    dm_props = next(p for _, lbl, p in nodes if lbl == "derived_metric")
    assert "p_value_threshold" not in dm_props  # treated as absent


def test_numeric_p_value_threshold_invalid_logs_warning(tmp_path, caplog):
    """A non-numeric p_value_threshold must log a warning and not set the property."""
    import logging
    csv_path = tmp_path / "num.csv"
    csv_path.write_text("locus_tag,v\nPMM0001,0.5\n")
    config = {
        "publication": {
            "papername": "Test", "doi": "10.9999/t",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "some_metric",
                        "value_kind": "numeric",
                        "value_col": "v",
                        "rankable": "true",
                        "has_p_value": "true",
                        "p_value_threshold": "not_a_number",
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=str(pc_path))
    with caplog.at_level(logging.WARNING, logger="multiomics_kg.adapters.observations_adapter"):
        nodes = adapter.get_nodes()
    dm_props = next(p for _, lbl, p in nodes if lbl == "derived_metric")
    assert "p_value_threshold" not in dm_props
    assert any("Invalid p_value_threshold" in rec.message for rec in caplog.records)


def _write_s5_like_paperconfig(tmp_path):
    """Minimal in-memory paperconfig mimicking Biller 2018 S5 categorical shape."""
    csv_path = tmp_path / "s5_small.csv"
    csv_path.write_text(
        "NCBI ID_2,darkness_cluster\n"
        "PMN2A_RS00030,darkness_axenic+darkness_coculture\n"
        "PMN2A_RS00040,darkness_coculture+unique_coculture\n"
        "PMN2A_RS00050,darkness_axenic+unique_axenic\n"
    )
    resolved_path = tmp_path / "s5_small_resolved.csv"
    resolved_path.write_text(
        "NCBI ID_2,darkness_cluster,resolved_locus_tag,resolution_method\n"
        "PMN2A_RS00030,darkness_axenic+darkness_coculture,PMN2A_1330,tier1:NCBI ID_2\n"
        "PMN2A_RS00040,darkness_coculture+unique_coculture,PMN2A_1331,tier1:NCBI ID_2\n"
        "PMN2A_RS00050,darkness_axenic+unique_axenic,PMN2A_1332,tier1:NCBI ID_2\n"
    )
    config = {
        "publication": {
            "papername": "Biller 2018",
            "doi": "10.1128/mSystems.00040-18",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {
                "axenic_rnaseq": {
                    "name": "NATL2A axenic",
                    "organism": "Prochlorococcus NATL2A",
                    "omics_type": "RNASEQ",
                    "treatment_type": ["darkness"],
                    "background_factors": ["axenic", "diel"],
                    "treatment_condition": "Extended darkness",
                    "light_condition": "continuous darkness",
                    "experimental_context": "Axenic NATL2A in Pro99",
                },
            },
            "supplementary_materials": {
                "s5_survival": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus NATL2A",
                    "experiment": "axenic_rnaseq",
                    "name_col": "NCBI ID_2",
                    "metrics": [
                        {
                            "metric_type": "darkness_survival_class",
                            "value_kind": "categorical",
                            "value_col": "darkness_cluster",
                            "allowed_categories": [
                                "darkness_axenic+darkness_coculture",
                                "darkness_coculture+unique_coculture",
                                "darkness_axenic+unique_axenic",
                            ],
                            "field_description": "Transcript presence at 72-144h",
                        },
                    ],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    return str(pc_path)


def test_categorical_dm_node_has_allowed_categories_from_paperconfig(tmp_path):
    pc_path = _write_s5_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()
    dm_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "derived_metric"]
    assert len(dm_nodes) == 1
    _, _, props = dm_nodes[0]
    assert props["value_kind"] == "categorical"
    assert props["rankable"] == "false"
    assert props["has_p_value"] == "false"
    assert props["allowed_categories"] == [
        "darkness_axenic+darkness_coculture",
        "darkness_coculture+unique_coculture",
        "darkness_axenic+unique_axenic",
    ]
    assert props["unit"] == ""
    assert props["metric_type"] == "darkness_survival_class"
