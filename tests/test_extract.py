# tests/test_extract.py
import json
import pytest
from pathlib import Path


def test_build_context_block():
    """Context block includes organism, method, treatment."""
    from multiomics_kg.extraction.cluster.extract import build_context_block

    table = {
        "name": "MIT9313 N-starvation",
        "organism": "Prochlorococcus MIT9313",
        "cluster_method": "K-means (K=7)",
        "cluster_type": "response_pattern",
        "treatment": "N-starvation time course",
        "omics_type": "MICROARRAY",
    }
    result = build_context_block(table)
    assert "Prochlorococcus MIT9313" in result
    assert "K-means (K=7)" in result
    assert "N-starvation time course" in result
    assert "MICROARRAY" in result


def test_format_cluster_summaries():
    """Summary text has all cluster keys."""
    from multiomics_kg.extraction.cluster.extract import format_cluster_summaries

    clusters = {
        "1": {"gene_count": 7, "sample_genes": ["PMT001", "PMT002"]},
        "6": {"gene_count": 81, "sample_genes": ["PMT100"]},
    }
    result = format_cluster_summaries(clusters)
    assert "Cluster 1:" in result
    assert "Cluster 6:" in result
    assert "7 genes" in result
    assert "81 genes" in result


def test_generate_report_stable_order(tmp_path):
    """Report generated twice is identical (no ordering jitter)."""
    from multiomics_kg.extraction.cluster.extract import generate_report
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    metadata = {"paper": "Test 2006"}
    clusters = {
        "2": {"id": "t_2", "name": "Cluster 2", "direction": "down",
              "self_assessment": "medium", "functional_description": "Photo",
              "behavioral_description": "Down", "confidence_notes": ""},
        "1": {"id": "t_1", "name": "Cluster 1", "direction": "up",
              "self_assessment": "high", "functional_description": "Transport",
              "behavioral_description": "Up", "confidence_notes": ""},
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)

    entries = [(tmp_path, "test_entry", {}, {"papername": "Test 2006"})]
    report1 = generate_report(entries)
    report2 = generate_report(entries)
    assert report1 == report2
    assert report1.index("Cluster 1") < report1.index("Cluster 2")
