# tests/test_extraction_merge.py
"""Tests for list-based merge (no LLM)."""
from multiomics_kg.extraction.cluster.merge import merge_paths

FIELDS = ["enrichment_category", "genes", "direction", "temporal_pattern"]


class TestMergePaths:
    def test_concatenates_all_sources(self):
        visual = {"1": {"enrichment_category": "Transport", "direction": "up"}}
        semantic = {"1": {"enrichment_category": "Transport and binding", "genes": "urtA, cynA"}}
        table = {"1": {"enrichment_category": "Transport and binding proteins",
                        "enrichment_pvalue": 0.0134, "gene_count": 5}}

        merged = merge_paths(["1"], visual, semantic, table)
        enr = merged["1"]["enrichment_category"]
        assert len(enr) == 3
        # Sorted by confidence: table (very_high) first
        assert enr[0]["source"] == "table"
        assert enr[0]["confidence"] == "very_high"

    def test_single_source_still_list(self):
        visual = {"1": {"direction": "up"}}
        merged = merge_paths(["1"], visual, {}, {})
        assert len(merged["1"]["direction"]) == 1
        assert merged["1"]["direction"][0]["value"] == "up"

    def test_quotes_pooled_by_relevance(self):
        visual = {"1": {"supporting_quotes": [
            {"quote": "from visual", "relevance_score": 0.9}
        ]}}
        semantic = {"1": {"supporting_quotes": [
            {"quote": "from semantic", "relevance_score": 0.7}
        ]}}
        merged = merge_paths(["1"], visual, semantic, {})
        quotes = merged["1"]["supporting_quotes"]
        assert len(quotes) == 2
        assert quotes[0]["relevance_score"] >= quotes[1]["relevance_score"]

    def test_missing_cluster_in_some_paths(self):
        visual = {"1": {"direction": "up"}, "2": {"direction": "down"}}
        semantic = {"1": {"genes": "urtA"}}
        table = {}
        merged = merge_paths(["1", "2"], visual, semantic, table)
        assert "genes" in merged["1"]
        assert "direction" in merged["2"]
        assert len(merged["2"]["direction"]) == 1

    def test_gene_count_from_table(self):
        table = {"1": {"gene_count": 5, "genes": ["PMM0970", "PMM0920"]}}
        merged = merge_paths(["1"], {}, {}, table)
        assert merged["1"]["gene_count"] == 5
