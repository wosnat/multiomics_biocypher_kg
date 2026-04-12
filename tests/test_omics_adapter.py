"""
Tests for OMICS adapter edge property propagation.

Focused tests for specific edge properties that flow from paperconfig analysis
entries through the adapter into Changes_expression_of edge property dicts.
"""

import pytest
import yaml
import pandas as pd

from multiomics_kg.adapters.omics_adapter import OMICSAdapter


def _make_growth_phase_config(tmp_path):
    """Build a minimal paperconfig and DE CSV with growth_phase set."""
    de_csv = tmp_path / "de.csv"
    pd.DataFrame({
        "gene": ["PMM0001", "PMM0002"],
        "log2fc": [1.5, -2.0],
    }).to_csv(de_csv, index=False)

    paperconfig = {
        "publication": {
            "papername": "Test",
            "papermainpdf": str(tmp_path / "dummy.pdf"),
            "experiments": {
                "exp1": {
                    "name": "Test experiment",
                    "organism": "Prochlorococcus MED4",
                    "treatment_condition": "test",
                    "control_condition": "control",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                },
            },
            "supplementary_materials": {
                "tbl": {
                    "type": "csv",
                    "filename": str(de_csv),
                    "statistical_analyses": [{
                        "id": "DE_test",
                        "experiment": "exp1",
                        "timepoint": "24h",
                        "timepoint_hours": 24,
                        "growth_phase": "nutrient_limited",
                        "name_col": "gene",
                        "logfc_col": "log2fc",
                    }],
                },
            },
        },
    }

    config_file = tmp_path / "paperconfig.yaml"
    config_file.write_text(yaml.dump(paperconfig))
    return str(config_file)


def test_growth_phase_flows_to_edge_properties(tmp_path):
    """If an analysis has growth_phase set, the edge carries it."""
    config_file = _make_growth_phase_config(tmp_path)
    adapter = OMICSAdapter(config_file=config_file)
    adapter.extracted_data = {
        "publication": {
            "publication_id": "test_2024",
            "doi": "10.1234/test.2024",
        },
    }

    expression_edges = [
        edge for edge in adapter.get_edges()
        if edge[3] == "changes_expression_of"
    ]

    assert len(expression_edges) == 2, (
        f"Expected 2 expression edges, got {len(expression_edges)}"
    )
    for edge in expression_edges:
        edge_id, source, target, edge_type, props = edge
        assert "growth_phase" in props, (
            f"Edge {edge_id} missing growth_phase property; got props: {props}"
        )
        assert props["growth_phase"] == "nutrient_limited", (
            f"Edge {edge_id} has wrong growth_phase: {props['growth_phase']!r}"
        )


def test_growth_phase_absent_when_not_set(tmp_path):
    """When growth_phase is not in the analysis, the edge omits the property."""
    de_csv = tmp_path / "de.csv"
    pd.DataFrame({
        "gene": ["PMM0001"],
        "log2fc": [1.5],
    }).to_csv(de_csv, index=False)

    paperconfig = {
        "publication": {
            "papername": "Test No GP",
            "papermainpdf": str(tmp_path / "dummy.pdf"),
            "experiments": {
                "exp1": {
                    "name": "Test experiment no gp",
                    "organism": "Prochlorococcus MED4",
                    "treatment_condition": "test",
                    "control_condition": "control",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                },
            },
            "supplementary_materials": {
                "tbl": {
                    "type": "csv",
                    "filename": str(de_csv),
                    "statistical_analyses": [{
                        "id": "DE_no_gp",
                        "experiment": "exp1",
                        "name_col": "gene",
                        "logfc_col": "log2fc",
                    }],
                },
            },
        },
    }

    config_file = tmp_path / "paperconfig_no_gp.yaml"
    config_file.write_text(yaml.dump(paperconfig))

    adapter = OMICSAdapter(config_file=str(config_file))
    adapter.extracted_data = {
        "publication": {
            "publication_id": "test_no_gp",
            "doi": "10.5678/test.2024",
        },
    }

    expression_edges = [
        edge for edge in adapter.get_edges()
        if edge[3] == "changes_expression_of"
    ]

    assert len(expression_edges) == 1
    _, _, _, _, props = expression_edges[0]
    assert "growth_phase" not in props, (
        f"growth_phase should be absent when not set, got: {props.get('growth_phase')!r}"
    )


def test_growth_phase_sanitized(tmp_path):
    """growth_phase with special chars is sanitized via clean_text."""
    de_csv = tmp_path / "de.csv"
    pd.DataFrame({
        "gene": ["PMM0001"],
        "log2fc": [1.5],
    }).to_csv(de_csv, index=False)

    paperconfig = {
        "publication": {
            "papername": "Test Sanitize",
            "papermainpdf": str(tmp_path / "dummy.pdf"),
            "experiments": {
                "exp1": {
                    "name": "Test experiment sanitize",
                    "organism": "Prochlorococcus MED4",
                    "treatment_condition": "test",
                    "control_condition": "control",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                },
            },
            "supplementary_materials": {
                "tbl": {
                    "type": "csv",
                    "filename": str(de_csv),
                    "statistical_analyses": [{
                        "id": "DE_sanitize",
                        "experiment": "exp1",
                        "name_col": "gene",
                        "logfc_col": "log2fc",
                        "growth_phase": "it's|a|test",
                    }],
                },
            },
        },
    }

    config_file = tmp_path / "paperconfig_sanitize.yaml"
    config_file.write_text(yaml.dump(paperconfig))

    adapter = OMICSAdapter(config_file=str(config_file))
    adapter.extracted_data = {
        "publication": {
            "publication_id": "test_sanitize",
            "doi": "10.9999/test.sanitize",
        },
    }

    expression_edges = [
        edge for edge in adapter.get_edges()
        if edge[3] == "changes_expression_of"
    ]

    assert len(expression_edges) == 1
    _, _, _, _, props = expression_edges[0]
    # clean_text strips | → "" and ' → "^"
    assert "growth_phase" in props
    assert "|" not in props["growth_phase"]
    assert "'" not in props["growth_phase"]
