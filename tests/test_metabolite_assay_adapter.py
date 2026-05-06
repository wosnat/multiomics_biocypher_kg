"""Tests for MetaboliteAssayAdapter."""
import yaml
import pandas as pd


# ─── _aggregate_replicates ─────────────────────────────────────────────────────


def test_aggregate_replicates_all_detected():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    mean, sd, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["1.0", "2.0", "3.0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3
    assert n_nz == 3
    assert mean == 2.0
    assert vals == [1.0, 2.0, 3.0]
    assert det == "detected"


def test_aggregate_replicates_all_zero():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, _, det = _aggregate_replicates(
        ["0", "0", "0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3 and n_nz == 0 and det == "not_detected"


def test_aggregate_replicates_sporadic():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, _, det = _aggregate_replicates(
        ["1.5", "0", "0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3 and n_nz == 1 and det == "sporadic"


def test_aggregate_replicates_null_values_become_zero():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["nd", "1.0", "NA"], null_values={"nd", "NA"}, missing_values={""}
    )
    assert n_rep == 3 and n_nz == 1 and det == "sporadic"
    assert vals == [0.0, 1.0, 0.0]


def test_aggregate_replicates_missing_excluded():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["", "1.0", "2.0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 2 and n_nz == 2 and det == "detected"
    assert vals == [1.0, 2.0]


# ─── parse_embedded_mean_sd_n ──────────────────────────────────────────────────


def test_parse_embedded_cell_normal():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    result = parse_embedded_mean_sd_n("0.00054 (8.8e-05), n=2")
    assert result == (0.00054, 8.8e-05, 2)


def test_parse_embedded_cell_nd_returns_zeros():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    assert parse_embedded_mean_sd_n("nd") == (0.0, 0.0, 0)
    assert parse_embedded_mean_sd_n("ND") == (0.0, 0.0, 0)


def test_parse_embedded_cell_na_sd():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    result = parse_embedded_mean_sd_n("0.0016 (NA), n=1")
    assert result[0] == 0.0016
    assert result[1] == 0.0
    assert result[2] == 1


def test_parse_embedded_cell_blank_returns_none():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    assert parse_embedded_mean_sd_n("") is None


# ─── Adapter integration tests ────────────────────────────────────────────────


def _make_minimal_paperconfig(tmp_path, with_resolved=True, value_kind="numeric"):
    """Create a minimal paperconfig with one metabolite_assays_table entry."""
    src = tmp_path / "metab.csv"
    if value_kind == "numeric":
        pd.DataFrame({"compound": ["Glucose"], "c1": [1.0], "c2": [3.0]}).to_csv(src, index=False)
        if with_resolved:
            pd.DataFrame({
                "compound": ["Glucose"],
                "c1": ["1.0"],
                "c2": ["3.0"],
                "metabolite_id": ["kegg.compound:C00031"],
                "resolution_method": ["name_match"],
            }).to_csv(src.with_name("metab_resolved.csv"), index=False)
        sample_columns = [{"condition_label": "control", "replicate_columns": ["c1", "c2"]}]
    else:  # boolean
        pd.DataFrame({"compound": ["GABA", "Glucose"], "intra": ["yes", ""]}).to_csv(src, index=False)
        if with_resolved:
            pd.DataFrame({
                "compound": ["GABA", "Glucose"],
                "intra": ["yes", ""],
                "metabolite_id": ["kegg.compound:C00334", "kegg.compound:C00031"],
                "resolution_method": ["alias_override", "kegg_direct"],
            }).to_csv(src.with_name("metab_resolved.csv"), index=False)
        sample_columns = [{"condition_label": "", "flag_column": "intra", "flag_true_value": "yes"}]

    cfg = {
        "publication": {
            "papername": "Test 2026",
            "doi": "10.1/test",
            "experiments": {
                "exp1": {
                    "name": "test exp",
                    "organism": "Prochlorococcus MIT9303",
                    "compartment": "whole_cell",
                    "omics_type": "METABOLOMICS",
                    "treatment_type": ["nitrogen"],
                    "background_factors": ["axenic"],
                    "treatment_condition": "stress",
                    "experimental_context": "ctx",
                    "light_condition": "continuous",
                },
            },
            "supplementary_materials": {
                "tab_a": {
                    "type": "metabolite_assays_table",
                    "filename": str(src),
                    "experiment": "exp1",
                    "organism": "Prochlorococcus MIT9303",
                    "name_col": "compound",
                    "assays": [{
                        "metric_type": "test_metric",
                        "name": "test assay",
                        "value_kind": value_kind,
                        "unit": "fg/cell" if value_kind == "numeric" else "",
                        "rankable": "true" if value_kind == "numeric" else "false",
                        "field_description": "test desc",
                        "sample_columns": sample_columns,
                    }],
                }
            },
        }
    }
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text(yaml.safe_dump(cfg))
    return pc


def test_adapter_emits_one_node_per_assay(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path)
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    nodes = list(adapter.get_nodes())
    assert len(nodes) == 1
    node_id, label, props = nodes[0]
    assert label == "metabolite_assay"
    assert "tab_a" in node_id
    assert "test_metric" in node_id
    assert props["metric_type"] == "test_metric"
    assert props["value_kind"] == "numeric"
    assert props["compartment"] == "whole_cell"
    assert props["organism_name"] == "Prochlorococcus MIT9303"
    assert props["unit"] == "fg/cell"
    assert props["rankable"] == "true"
    assert props["aggregation_method"] == "mean_across_replicates"
    assert props["treatment_type"] == ["nitrogen"]


def test_adapter_emits_numeric_quantifies_edges(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path)
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    edges = list(adapter.get_edges())
    quant = [e for e in edges if e[3] == "assay_quantifies_metabolite"]
    assert len(quant) == 1
    edge_id, src_id, dst_id, label, props = quant[0]
    assert dst_id == "kegg.compound:C00031"
    assert props["value"] == 2.0
    assert props["n_replicates"] == 2
    assert props["n_non_zero"] == 2
    # Pre-stringified for BioCypher's _write_array_string (which joins without
    # type coercion); neo4j-admin parses them back to floats per the schema's
    # `float[]` declaration.
    assert props["replicate_values"] == ["1.0", "3.0"]
    assert props["detection_status"] == "detected"
    assert props["condition_label"] == "control"


def test_adapter_emits_boolean_flag_edges(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path, value_kind="boolean")
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    edges = [e for e in adapter.get_edges() if e[3] == "assay_flags_metabolite"]
    assert len(edges) == 2
    by_dst = {e[2]: e[4] for e in edges}
    assert by_dst["kegg.compound:C00334"]["flag_value"] == "true"
    assert by_dst["kegg.compound:C00334"]["n_positive"] == 1
    assert by_dst["kegg.compound:C00031"]["flag_value"] == "false"
    assert by_dst["kegg.compound:C00031"]["n_positive"] == 0


def test_adapter_emits_three_binding_edges(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path)
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    # MultiMetaboliteAssayAdapter normally injects this from cyanobacteria_genomes.csv;
    # tests instantiating the single-paper adapter directly must do it explicitly
    # or the organism binding edge is silently dropped (target node ID unknown).
    adapter._organism_lookup = {"Prochlorococcus MIT9303": "insdc.gcf:GCF_000015965.1"}
    labels = [e[3] for e in adapter.get_edges()]
    assert labels.count("publication_has_metabolite_assay") == 1
    assert labels.count("experiment_has_metabolite_assay") == 1
    assert labels.count("metabolite_assay_belongs_to_organism") == 1


def test_adapter_no_resolved_csv_still_emits_binding_edges(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path, with_resolved=False)
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    adapter._organism_lookup = {"Prochlorococcus MIT9303": "insdc.gcf:GCF_000015965.1"}
    edges = list(adapter.get_edges())
    measurement = [e for e in edges if e[3].startswith("assay_")]
    binding = [e for e in edges if not e[3].startswith("assay_")]
    # No measurement edges (no _resolved.csv to look up metabolite_id)
    assert measurement == []
    # But binding edges are still emitted (so MetaboliteAssay isn't an orphan)
    assert len(binding) == 3


# ─── Adapter integration: cell_format=embedded_mean_sd_n ─────────────────────


def _make_embedded_paperconfig(tmp_path, total_replicates=3, drop_undetected=False):
    """Paperconfig using cell_format: embedded_mean_sd_n.

    Three rows exercise the three detection states:
    - Glucose "0.18 (0.049), n=3" → fully detected
    - Aconitic Acid "0.00054 (8.8e-05), n=2" → sporadic
    - Ribose 5 phosphate "nd" → not_detected
    """
    src = tmp_path / "metab.csv"
    pd.DataFrame({
        "compound": ["Glucose", "Aconitic Acid", "Ribose 5 phosphate"],
        "cells": ["0.18 (0.049), n=3", "0.00054 (8.8e-05), n=2", "nd"],
    }).to_csv(src, index=False)
    pd.DataFrame({
        "compound": ["Glucose", "Aconitic Acid", "Ribose 5 phosphate"],
        "cells": ["0.18 (0.049), n=3", "0.00054 (8.8e-05), n=2", "nd"],
        "metabolite_id": [
            "kegg.compound:C00031",
            "kegg.compound:C00417",
            "kegg.compound:C00117",
        ],
        "resolution_method": ["name_match"] * 3,
    }).to_csv(src.with_name("metab_resolved.csv"), index=False)

    cfg = {
        "publication": {
            "papername": "Test 2026 Embedded",
            "doi": "10.1/embedded",
            "experiments": {
                "exp1": {
                    "name": "embedded test exp",
                    "organism": "Prochlorococcus MIT9312",
                    "compartment": "whole_cell",
                    "omics_type": "METABOLOMICS",
                    "treatment_type": ["compartment"],
                    "background_factors": ["axenic", "light"],
                    "treatment_condition": "Whole-cell metabolome",
                    "experimental_context": "ctx",
                    "light_condition": "continuous",
                },
            },
            "supplementary_materials": {
                "tab_a": {
                    "type": "metabolite_assays_table",
                    "filename": str(src),
                    "experiment": "exp1",
                    "organism": "Prochlorococcus MIT9312",
                    "name_col": "compound",
                    "cell_format": "embedded_mean_sd_n",
                    "aggregation_method": "pre_aggregated",
                    "drop_undetected": drop_undetected,
                    "assays": [{
                        "metric_type": "rel_peak_area",
                        "name": "test embedded",
                        "value_kind": "numeric",
                        "unit": "rel area",
                        "rankable": "true",
                        "field_description": "test desc",
                        "total_replicates": total_replicates,
                        "sample_columns": [
                            {
                                "condition_label": "cells",
                                "replicate_columns": ["cells"],
                            }
                        ],
                    }],
                }
            },
        }
    }
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text(yaml.safe_dump(cfg))
    return pc


def test_adapter_embedded_format_three_detection_states(tmp_path):
    """cell_format=embedded_mean_sd_n parses '<mean> (<sd>), n=N' cells.

    n_replicates is taken from paperconfig total_replicates (3 here);
    n_non_zero is the parsed n (number of replicates above LOD).
    Detection status: parsed_n == 0 → not_detected; == total → detected; else sporadic.
    replicate_values stay empty since the embedded format gives only summary stats.
    """
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_embedded_paperconfig(tmp_path)
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    edges = [e for e in adapter.get_edges() if e[3] == "assay_quantifies_metabolite"]
    assert len(edges) == 3
    by_dst = {e[2]: e[4] for e in edges}

    glu = by_dst["kegg.compound:C00031"]
    assert glu["value"] == 0.18
    assert glu["value_sd"] == 0.049
    assert glu["n_replicates"] == 3
    assert glu["n_non_zero"] == 3
    assert glu["detection_status"] == "detected"
    assert glu["replicate_values"] == []

    aco = by_dst["kegg.compound:C00417"]
    assert aco["value"] == 0.00054
    assert aco["value_sd"] == 8.8e-05
    assert aco["n_replicates"] == 3
    assert aco["n_non_zero"] == 2
    assert aco["detection_status"] == "sporadic"

    rib = by_dst["kegg.compound:C00117"]
    assert rib["value"] == 0.0
    assert rib["value_sd"] == 0.0
    assert rib["n_replicates"] == 3
    assert rib["n_non_zero"] == 0
    assert rib["detection_status"] == "not_detected"


def test_multi_adapter_filters_to_metabolomics_papers(tmp_path):
    """MultiMetaboliteAssayAdapter only instantiates per-paper adapters when
    the paperconfig has at least one metabolite_assays_table entry."""
    from multiomics_kg.adapters.metabolite_assay_adapter import MultiMetaboliteAssayAdapter

    pc_with = _make_minimal_paperconfig(tmp_path)

    # Make a paperconfig that does NOT have metabolite_assays_table
    pc_without = tmp_path / "pc_de.yaml"
    pc_without.write_text(yaml.safe_dump({
        "publication": {
            "papername": "DE Paper",
            "doi": "10.1/de",
            "experiments": {"exp1": {"organism": "P", "compartment": "whole_cell"}},
            "supplementary_materials": {
                "tab": {"type": "csv", "filename": "ignored.csv"},
            },
        }
    }))

    multi = MultiMetaboliteAssayAdapter(paperconfig_paths=[pc_with, pc_without])
    assert len(multi.adapters) == 1
    nodes = list(multi.get_nodes())
    assert len(nodes) == 1
