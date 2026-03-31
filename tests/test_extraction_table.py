# tests/test_extraction_table.py
"""Tests for Path: table (no LLM, fast)."""
from pathlib import Path

import pytest

from multiomics_kg.extraction.cluster.table import extract_from_csv, run_table

PROJECT_ROOT = Path(__file__).parent.parent


class TestExtractFromCsv:
    def test_tolonen_med4_basic(self):
        table_config = {
            "filename": "data/Prochlorococcus/papers_and_supp/tolonen 2006/med4_kmeans_clusters.csv",
            "gene_id_col": "gene_id",
            "cluster_col": "cluster",
        }
        result = extract_from_csv(table_config, project_root=PROJECT_ROOT)
        assert "1" in result
        assert result["1"]["gene_count"] == 5
        assert result["6"]["gene_count"] == 124
        assert "PMM0970" in result["1"]["genes"]

    def test_zinser_has_peak_time(self):
        table_config = {
            "filename": "data/Prochlorococcus/papers_and_supp/zinser 2009/diel_clusters.csv",
            "gene_id_col": "gene_id",
            "cluster_col": "cluster",
            "score_col": "membership_score",
        }
        result = extract_from_csv(table_config, project_root=PROJECT_ROOT)
        # Zinser CSV has peak_time column — should be auto-detected
        assert "5" in result
        assert result["5"]["gene_count"] > 0
        assert "peak_time_mean" in result["5"]

    def test_alonso_cluster_keys(self):
        table_config = {
            "filename": "data/Prochlorococcus/papers_and_supp/alonso 2023/soft_clusters.csv",
            "gene_id_col": "locus_tag",
            "cluster_col": "cluster",
            "score_col": "prob_cluster_A",
        }
        result = extract_from_csv(table_config, project_root=PROJECT_ROOT)
        assert "Cluster A" in result
        assert result["Cluster A"]["gene_count"] > 0


class TestRunTable:
    def test_tolonen_with_enrichment_xls(self):
        table_config = {
            "filename": "data/Prochlorococcus/papers_and_supp/tolonen 2006/med4_kmeans_clusters.csv",
            "gene_id_col": "gene_id",
            "cluster_col": "cluster",
        }
        paper_dir = PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp/tolonen 2006"
        cluster_keys = [str(i) for i in range(1, 10)]
        result = run_table(table_config, paper_dir, cluster_keys,
                           organism_hint="MED4", project_root=PROJECT_ROOT)
        # Should have enrichment from XLS
        assert result["1"]["enrichment_category"] == "Transport and binding proteins"
        assert result["7"]["enrichment_category"] == "Translation"
        assert result["7"]["enrichment_pvalue"] < 1e-8
        # Should also have genes from CSV
        assert result["1"]["gene_count"] == 5
