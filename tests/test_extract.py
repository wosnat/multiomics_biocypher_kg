# tests/test_extract.py
import json
import pytest
from pathlib import Path


def test_build_context_block():
    """Context block includes organism, method, treatment, extraction_notes."""
    from multiomics_kg.extraction.cluster.extract import build_context_block

    table = {
        "name": "MIT9313 N-starvation",
        "organism": "Prochlorococcus MIT9313",
        "cluster_method": "K-means (K=7)",
        "cluster_type": "time_course",
        "treatment": "N-starvation time course",
        "omics_type": "MICROARRAY",
        "extraction_notes": "Paper discusses both MED4 and MIT9313 jointly.",
    }
    result = build_context_block(table)
    assert "Prochlorococcus MIT9313" in result
    assert "K-means (K=7)" in result
    assert "N-starvation time course" in result
    assert "MICROARRAY" in result
    assert "Paper discusses both MED4 and MIT9313 jointly." in result


def test_build_prompt():
    """build_prompt assembles shared rules + type rules."""
    from multiomics_kg.extraction.cluster.extract import build_prompt

    table = {
        "cluster_type": "diel",
        "name": "Test",
        "organism": "Test",
        "cluster_method": "K-means",
        "treatment": "light",
    }
    summaries = {"1": {"gene_count": 50, "sample_genes": ["g1"]}}
    result = build_prompt(table, summaries)
    assert "N/A" in result
    assert "peaks at dawn" in result
    assert "Self-Verification" in result
    assert "EXACTLY 1 clusters" in result


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
        "2": {"id": "t_2", "name": "Cluster 2", "expression_dynamics": "late sustained",
              "self_assessment": "medium", "functional_description": "Photo",
              "temporal_pattern": "Down", "confidence_notes": ""},
        "1": {"id": "t_1", "name": "Cluster 1", "expression_dynamics": "early transient",
              "self_assessment": "high", "functional_description": "Transport",
              "temporal_pattern": "Up", "confidence_notes": ""},
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)

    entries = [(tmp_path, "test_entry", {}, {"papername": "Test 2006"})]
    report1 = generate_report(entries)
    report2 = generate_report(entries)
    assert report1 == report2
    assert report1.index("Cluster 1") < report1.index("Cluster 2")


def test_detect_filler_on_low_confidence():
    """Low-confidence cluster with non-N/A description produces warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "test_low", "name": "C1",
                  "functional_description": "Some vague filler description here",
                  "temporal_pattern": "N/A",
                  "expression_dynamics": "N/A", "self_assessment": "low"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("low confidence" in w for w in warnings)


def test_detect_locus_tags():
    """Locus tag in description produces warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "x", "name": "C1",
                  "functional_description": "Includes gene PMM0042 involved in transport",
                  "temporal_pattern": "Upregulated",
                  "expression_dynamics": "up"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("locus tag" in w and "PMM0042" in w for w in warnings)
