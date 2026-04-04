"""
Tests for ClusterAdapter and MultiClusterAdapter.

Tests verify:
1. get_nodes emits ClusteringAnalysis + GeneCluster nodes
2. get_edges emits 5 edge types (membership, pub→analysis, analysis→cluster, analysis→org, exp→analysis)
3. Unresolved (NaN) genes are skipped
4. Configs without score_col produce edges with empty props
5. Extraction JSON populates descriptions
6. Failed extraction (verdict=fail) skips descriptions
7. Node IDs use entry_key instead of organism suffix
"""
import json
import os
import tempfile

import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.cluster_adapter import (
    ClusterAdapter,
    MultiClusterAdapter,
    _make_cluster_id,
    _make_analysis_id,
    _resolve_csv_path,
    _load_extraction_json,
    _get_extraction_cluster_data,
)


# ─── Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


def _write_cluster_csv(path: str, with_resolved: bool = False, include_nan: bool = False) -> str:
    """Write a cluster CSV with 5 genes across 2 clusters.

    Returns the path actually written (original or _resolved).
    """
    data = {
        "gene_id": ["PMM0001", "PMM0002", "PMM0003", "PMM0004", "PMM0005"],
        "cluster": [1, 1, 2, 2, 2],
        "score": [0.9, 0.8, 0.7, 0.85, 0.6],
    }
    if include_nan:
        data["resolved_locus_tag"] = ["PMM0001", None, "PMM0003", "PMM0004", "PMM0005"]
        data["resolution_method"] = ["tier1:locus_tag", None, "tier1:locus_tag", "tier1:locus_tag", "tier1:locus_tag"]

    df = pd.DataFrame(data)

    if with_resolved:
        stem = os.path.splitext(path)[0]
        resolved_path = f"{stem}_resolved.csv"
        if "resolved_locus_tag" not in df.columns:
            df["resolved_locus_tag"] = df["gene_id"]
            df["resolution_method"] = "tier1:locus_tag"
        df.to_csv(resolved_path, index=False)
        # Also write original so the config filename points to a real file
        df.drop(columns=["resolved_locus_tag", "resolution_method"], errors="ignore").to_csv(path, index=False)
        return resolved_path
    else:
        df.to_csv(path, index=False)
        return path


def _write_paperconfig(
    temp_dir: str,
    csv_path: str,
    doi: str = "10.1234/test.2024",
    score_col: str = "score",
    gene_id_col: str = "gene_id",
) -> str:
    config = {
        "publication": {
            "papername": "Test 2024",
            "doi": doi,
            "supplementary_materials": {
                "cluster_table_1": {
                    "type": "gene_clusters",
                    "name": "MED4 test clusters",
                    "filename": csv_path,
                    "cluster_col": "cluster",
                    "gene_id_col": gene_id_col,
                    "score_col": score_col,
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "cluster_method": "k-means",
                    "cluster_type": "response_pattern",
                    "treatment_type": ["nitrogen"],
                    "treatment": "Nitrogen limitation",
                    "light_condition": "continuous light",
                    "experimental_context": "in Pro99 medium",
                }
            },
        }
    }
    config_path = os.path.join(temp_dir, "paperconfig.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


# ─── Test: ClusteringAnalysis node emission ─────────────────────────────────


def test_get_nodes_emits_clustering_analysis_node(tmp_path):
    """ClusterAdapter should yield one ClusteringAnalysis node per gene_clusters entry."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,1\nPMM0003,2\n")

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "med4_kmeans_test": {
                "type": "gene_clusters",
                "name": "MED4 K-means test clusters",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
                "cluster_method": "K-means (K=2)",
                "omics_type": "MICROARRAY",
                "treatment_type": ["nitrogen_stress"],
                "treatment": "N-starvation",
                "light_condition": "continuous light",
                "experimental_context": "test context",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    nodes = list(adapter.get_nodes())

    analysis_nodes = [(nid, label, props) for nid, label, props in nodes if label == "clustering_analysis"]
    cluster_nodes = [(nid, label, props) for nid, label, props in nodes if label == "gene_cluster"]

    assert len(analysis_nodes) == 1
    assert len(cluster_nodes) == 2

    a_id, a_label, a_props = analysis_nodes[0]
    assert "med4_kmeans_test" in a_id
    assert a_props["name"] == "MED4 K-means test clusters"
    assert a_props["organism_name"] == "Prochlorococcus MED4"
    assert a_props["cluster_count"] == 2
    assert a_props["total_gene_count"] == 3
    assert a_props["cluster_type"] == "response_pattern"


# ─── Test: GeneCluster node ID includes analysis key ────────────────────────


def test_gene_cluster_node_id_includes_analysis_key(tmp_path):
    """GeneCluster node ID should be cluster:{doi}:{analysis_key}:{csv_key}."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,A\nPMM0002,B\n")

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "med4_test_analysis": {
                "type": "gene_clusters",
                "name": "Test analysis",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    nodes = list(adapter.get_nodes())
    cluster_nodes = [(nid, label, props) for nid, label, props in nodes if label == "gene_cluster"]
    cluster_ids = sorted([nid for nid, _, _ in cluster_nodes])
    assert any("med4_test_analysis" in cid and ":A" in cid for cid in cluster_ids)
    assert any("med4_test_analysis" in cid and ":B" in cid for cid in cluster_ids)


# ─── Test: Extraction JSON reading ──────────────────────────────────────────


def test_get_nodes_reads_extraction_json(tmp_path):
    """Adapter should populate descriptions from extraction JSON."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    extraction = {
        "metadata": {"table_key": "test_analysis"},
        "stage2_results": {
            "1": {"id": "med4_up_transport", "name": "MED4 cluster 1 (up, transport)",
                  "functional_description": "Transport genes", "behavioral_description": "Rapid upregulation",
                  "peak_time_hours": None, "period_hours": None},
            "2": {"id": "med4_down_photosynthesis", "name": "MED4 cluster 2 (down, photosynthesis)",
                  "functional_description": "PSI genes", "behavioral_description": "Downregulated at 6h",
                  "peak_time_hours": None, "period_hours": None},
        },
        "stage3_validation": {"1": {"verdict": "pass"}, "2": {"verdict": "pass"}},
    }
    json_path = tmp_path / "cluster_extraction_test_analysis.json"
    json_path.write_text(json.dumps(extraction))

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters",
                "name": "Test analysis",
                "filename": str(csv_path),
                "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id",
                "cluster_col": "cluster",
                "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    nodes = list(adapter.get_nodes())
    cluster_nodes = {props.get("id"): props for _, label, props in nodes if label == "gene_cluster"}

    assert cluster_nodes["med4_up_transport"]["functional_description"] == "Transport genes"
    assert cluster_nodes["med4_up_transport"]["name"] == "MED4 cluster 1 (up, transport)"
    assert cluster_nodes["med4_down_photosynthesis"]["behavioral_description"] == "Downregulated at 6h"


# ─── Test: Failed extraction skipped ────────────────────────────────────────


def test_get_nodes_skips_failed_extraction(tmp_path):
    """Adapter should skip clusters with verdict=fail in extraction JSON."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    extraction = {
        "metadata": {"table_key": "test_analysis"},
        "stage2_results": {
            "1": {"id": "good", "name": "Good", "functional_description": "OK",
                  "behavioral_description": "OK", "peak_time_hours": None, "period_hours": None},
            "2": {"id": "bad", "name": "Bad", "functional_description": "Wrong",
                  "behavioral_description": "Wrong", "peak_time_hours": None, "period_hours": None},
        },
        "stage3_validation": {"1": {"verdict": "pass"}, "2": {"verdict": "fail"}},
    }
    json_path = tmp_path / "cluster_extraction_test_analysis.json"
    json_path.write_text(json.dumps(extraction))

    paperconfig = {
        "publication": {"papername": "Test 2006"},
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters", "name": "Test",
                "filename": str(csv_path), "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id", "cluster_col": "cluster", "cluster_type": "response_pattern",
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    nodes = list(adapter.get_nodes())
    cluster_nodes = {props.get("id", ""): props for _, label, props in nodes if label == "gene_cluster"}

    assert cluster_nodes["good"]["functional_description"] == "OK"
    # Cluster 2 verdict=fail → empty descriptions
    failed = [p for p in cluster_nodes.values() if p.get("id") in ("bad", "")]
    assert len(failed) == 1
    assert failed[0]["functional_description"] == ""


# ─── Test: New edge types ───────────────────────────────────────────────────


def test_get_edges_new_structure(tmp_path):
    """Adapter should emit 5 edge types including ClusteringAnalysis edges."""
    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMM0001,1\nPMM0002,2\n")

    paperconfig = {
        "publication": {
            "papername": "Test 2006",
            "doi": "10.1234/test.2006",
            "experiments": {"exp1": {"name": "test experiment"}},
        },
        "supplementary_materials": {
            "test_analysis": {
                "type": "gene_clusters", "name": "Test",
                "filename": str(csv_path), "organism": "Prochlorococcus MED4",
                "gene_id_col": "gene_id", "cluster_col": "cluster",
                "cluster_type": "response_pattern", "experiments": ["exp1"],
            }
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    edges = list(adapter.get_edges())
    edge_types = set(e[3] for e in edges)

    assert "gene_in_gene_cluster" in edge_types
    assert "publication_has_clustering_analysis" in edge_types
    assert "clustering_analysis_has_gene_cluster" in edge_types
    assert "clusteringanalysis_belongs_to_organism" in edge_types
    assert "experiment_has_clustering_analysis" in edge_types
    assert "publication_has_gene_cluster" not in edge_types
    assert "genecluster_belongs_to_organism" not in edge_types

    # 2 membership + 1 pub→analysis + 2 analysis→cluster + 1 analysis→org + 1 exp→analysis = 7
    assert len(edges) == 7


# ─── Test: get_nodes emits correct GeneCluster nodes (updated) ──────────────


def test_get_nodes_emits_cluster_nodes(temp_dir):
    """Two clusters → ClusteringAnalysis + 2 GeneCluster nodes with correct properties."""
    csv_path = os.path.join(temp_dir, "clusters.csv")
    _write_cluster_csv(csv_path, with_resolved=False)

    config_path = _write_paperconfig(temp_dir, csv_path)
    adapter = ClusterAdapter(config_file=config_path)
    nodes = adapter.get_nodes()

    # 1 ClusteringAnalysis + 2 GeneCluster
    analysis_nodes = [n for n in nodes if n[1] == "clustering_analysis"]
    cluster_nodes = [n for n in nodes if n[1] == "gene_cluster"]
    assert len(analysis_nodes) == 1
    assert len(cluster_nodes) == 2

    # Check analysis node
    a_id, a_label, a_props = analysis_nodes[0]
    assert "cluster_table_1" in a_id
    assert a_props["cluster_count"] == 2
    assert a_props["total_gene_count"] == 5
    assert a_props["organism_name"] == "Prochlorococcus MED4"

    # Check cluster IDs use entry_key
    ids = {n[0] for n in cluster_nodes}
    doi_short = "test.2024"
    expected_id_1 = f"cluster:{doi_short}:cluster_table_1:1"
    expected_id_2 = f"cluster:{doi_short}:cluster_table_1:2"
    assert expected_id_1 in ids, f"{expected_id_1} not in {ids}"
    assert expected_id_2 in ids, f"{expected_id_2} not in {ids}"

    # Check properties on cluster 1
    node_1 = next(n for n in cluster_nodes if n[0] == expected_id_1)
    props = node_1[2]
    assert props["organism_name"] == "Prochlorococcus MED4"
    assert props["omics_type"] == "RNASEQ"
    assert props["member_count"] == 2  # PMM0001 and PMM0002
    assert props["treatment_type"] == ["nitrogen"]

    # Check properties on cluster 2
    node_2 = next(n for n in cluster_nodes if n[0] == expected_id_2)
    props2 = node_2[2]
    assert props2["member_count"] == 3  # PMM0003, PMM0004, PMM0005


# ─── Test: get_edges emits membership edges with scores ─────────────────────


def test_get_edges_membership_with_scores(temp_dir):
    """5 gene membership edges total, ncbigene: prefix, scores as floats."""
    csv_path = os.path.join(temp_dir, "clusters.csv")
    _write_cluster_csv(csv_path, with_resolved=True)

    config_path = _write_paperconfig(temp_dir, csv_path)
    adapter = ClusterAdapter(config_file=config_path)
    edges = adapter.get_edges()

    membership_edges = [e for e in edges if e[3] == "gene_in_gene_cluster"]
    assert len(membership_edges) == 5, f"Expected 5 membership edges, got {len(membership_edges)}"

    # All gene target IDs must be prefixed with ncbigene:
    for edge_id, src, tgt, label, props in membership_edges:
        assert tgt.startswith("ncbigene:"), f"Gene ID '{tgt}' missing ncbigene: prefix"
        assert "membership_score" in props, f"Expected membership_score in edge props: {props}"
        assert isinstance(props["membership_score"], float)

    # Check pub→analysis edges exist (not pub→cluster)
    pub_edges = [e for e in edges if e[3] == "publication_has_clustering_analysis"]
    assert len(pub_edges) == 1, f"Expected 1 pub→analysis edge, got {len(pub_edges)}"
    assert pub_edges[0][1] == "doi:10.1234/test.2024"

    # Check analysis→cluster edges
    analysis_cluster_edges = [e for e in edges if e[3] == "clustering_analysis_has_gene_cluster"]
    assert len(analysis_cluster_edges) == 2


# ─── Test: Unresolved (NaN) genes are skipped ──────────────────────────────


def test_nan_locus_tags_skipped(temp_dir):
    """NaN locus_tag values produce no membership edges."""
    csv_path = os.path.join(temp_dir, "clusters.csv")

    # Build a resolved CSV manually with one NaN locus_tag
    df = pd.DataFrame({
        "gene_id": ["PMM0001", "PMM0002", "PMM0003", "PMM0004", "PMM0005"],
        "cluster": [1, 1, 2, 2, 2],
        "score": [0.9, 0.8, 0.7, 0.85, 0.6],
        "resolved_locus_tag": ["PMM0001", None, "PMM0003", "PMM0004", "PMM0005"],
        "resolution_method": [
            "tier1:locus_tag", None, "tier1:locus_tag", "tier1:locus_tag", "tier1:locus_tag"
        ],
    })
    # Write original (required by config filename)
    df.drop(columns=["resolved_locus_tag", "resolution_method"]).to_csv(csv_path, index=False)
    # Write resolved alongside
    resolved_path = csv_path.replace(".csv", "_resolved.csv")
    df.to_csv(resolved_path, index=False)

    config_path = _write_paperconfig(temp_dir, csv_path)
    adapter = ClusterAdapter(config_file=config_path)
    edges = adapter.get_edges()

    membership_edges = [e for e in edges if e[3] == "gene_in_gene_cluster"]
    # NaN in cluster 1 → only 1 gene, cluster 2 → 3 genes = 4 total
    assert len(membership_edges) == 4, (
        f"Expected 4 membership edges (NaN skipped), got {len(membership_edges)}"
    )

    # Ensure no edge target is NaN or empty
    for edge_id, src, tgt, label, props in membership_edges:
        assert tgt != "ncbigene:", f"Target should not be 'ncbigene:' (empty locus_tag)"
        assert tgt.startswith("ncbigene:")


# ─── Test: No score_col means edges have empty props ────────────────────────


def test_no_score_col_produces_empty_edge_props(temp_dir):
    """Without score_col the membership edge props dict should be empty."""
    csv_path = os.path.join(temp_dir, "clusters_noscore.csv")

    # Write a resolved CSV (no score column)
    df = pd.DataFrame({
        "gene_id": ["PMM0001", "PMM0002"],
        "cluster": [1, 1],
        "resolved_locus_tag": ["PMM0001", "PMM0002"],
        "resolution_method": ["tier1:locus_tag", "tier1:locus_tag"],
    })
    df.drop(columns=["resolved_locus_tag", "resolution_method"]).to_csv(csv_path, index=False)
    resolved_path = csv_path.replace(".csv", "_resolved.csv")
    df.to_csv(resolved_path, index=False)

    config = {
        "publication": {
            "papername": "Test No Score 2024",
            "doi": "10.1234/noscore.2024",
            "supplementary_materials": {
                "cluster_table_noscore": {
                    "type": "gene_clusters",
                    "filename": csv_path,
                    "cluster_col": "cluster",
                    "gene_id_col": "gene_id",
                    # No score_col
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "cluster_method": "hierarchical",
                    "cluster_type": "response_pattern",
                    "treatment_type": [],
                }
            },
        }
    }
    config_path = os.path.join(temp_dir, "paperconfig_noscore.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    adapter = ClusterAdapter(config_file=config_path)
    edges = adapter.get_edges()

    membership_edges = [e for e in edges if e[3] == "gene_in_gene_cluster"]
    assert len(membership_edges) == 2

    for edge_id, src, tgt, label, props in membership_edges:
        assert props == {}, f"Expected empty props dict, got {props}"


# ─── Test: _make_cluster_id with entry_key ──────────────────────────────────


def test_make_cluster_id_no_doi():
    """Without DOI, falls back to paper_name slug + entry_key."""
    cid = _make_cluster_id("", "Test Paper 2024", "med4_analysis", "cluster_A")
    assert cid == "cluster:test_paper_2024:med4_analysis:cluster_A"


def test_make_cluster_id_with_doi():
    """With DOI, uses last path component + entry_key."""
    cid = _make_cluster_id("10.1234/my.paper", "Test Paper 2024", "med4_analysis", "k1")
    assert cid == "cluster:my.paper:med4_analysis:k1"


def test_make_analysis_id_with_doi():
    """Analysis ID uses DOI short form + entry_key."""
    aid = _make_analysis_id("10.1234/my.paper", "Test Paper 2024", "med4_analysis")
    assert aid == "clustering_analysis:my.paper:med4_analysis"


def test_make_analysis_id_no_doi():
    """Analysis ID falls back to paper_name slug."""
    aid = _make_analysis_id("", "Test Paper 2024", "med4_analysis")
    assert aid == "clustering_analysis:test_paper_2024:med4_analysis"


# ─── Test: _resolve_csv_path prefers _resolved.csv ─────────────────────────


def test_resolve_csv_path_prefers_resolved(temp_dir):
    original = os.path.join(temp_dir, "data.csv")
    resolved = os.path.join(temp_dir, "data_resolved.csv")

    # Only original exists
    open(original, "w").close()
    path, is_resolved = _resolve_csv_path(original)
    assert str(path) == original
    assert not is_resolved

    # Now create resolved — should be preferred
    open(resolved, "w").close()
    path, is_resolved = _resolve_csv_path(original)
    assert str(path) == resolved
    assert is_resolved


# ─── Test: MultiClusterAdapter skips configs without gene_clusters ──────────


def test_get_nodes_resolved_csv_with_skip_rows(tmp_path):
    """get_nodes should ignore skip_rows when using a resolved CSV.

    Regression test: resolved CSVs are clean single-header files produced by
    resolve_paper_ids.py.  When skip_rows was applied to them, the real header
    was skipped and the cluster_col could not be found.
    """
    # Write original CSV with 3 junk header rows before the real header
    orig_csv = tmp_path / "table_s1.csv"
    orig_csv.write_text(
        "junk1\n"
        "junk2\n"
        "junk3\n"
        "gene_id,cluster,score\n"
        "PMM0001,5,0.67\n"
        "PMM0002,6,0.61\n"
    )
    # Write resolved CSV — clean single-header file
    resolved_csv = tmp_path / "table_s1_resolved.csv"
    resolved_csv.write_text(
        "gene_id,cluster,score,resolved_locus_tag,resolution_method\n"
        "PMM0001,5.0,0.67,PMM0001,locus_tag:gene_id\n"
        "PMM0002,6.0,0.61,PMM0002,locus_tag:gene_id\n"
    )

    paperconfig = {
        "publication": {
            "papername": "SkipRows 2024",
            "doi": "10.1234/skiprows.2024",
            "supplementary_materials": {
                "med4_diel_clusters": {
                    "type": "gene_clusters",
                    "name": "MED4 diel clusters",
                    "filename": str(orig_csv),
                    "skip_rows": 3,
                    "cluster_col": "cluster",
                    "gene_id_col": "gene_id",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "MICROARRAY",
                    "cluster_method": "k-means",
                    "cluster_type": "response_pattern",
                    "treatment_type": ["diel"],
                }
            },
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(paperconfig))

    adapter = ClusterAdapter(str(config_path))
    nodes = adapter.get_nodes()

    # Should produce 1 ClusteringAnalysis + 2 GeneCluster nodes
    analysis_nodes = [n for n in nodes if n[1] == "clustering_analysis"]
    cluster_nodes = [n for n in nodes if n[1] == "gene_cluster"]
    assert len(analysis_nodes) == 1
    assert len(cluster_nodes) == 2
    assert analysis_nodes[0][2]["total_gene_count"] == 2


# ─── Test: MultiClusterAdapter skips configs without gene_clusters ──────────


def test_multi_cluster_adapter_skips_non_cluster_configs(temp_dir):
    """MultiClusterAdapter only loads paperconfigs that have gene_clusters entries."""
    csv_path = os.path.join(temp_dir, "clusters.csv")
    _write_cluster_csv(csv_path, with_resolved=False)

    # Config WITH clusters
    config_with = _write_paperconfig(temp_dir, csv_path, doi="10.1/with")
    # Config WITHOUT clusters
    config_without = os.path.join(temp_dir, "paperconfig_no_clusters.yaml")
    with open(config_without, "w") as f:
        yaml.dump(
            {
                "publication": {
                    "papername": "No Clusters Paper",
                    "doi": "10.1/without",
                    "supplementary_materials": {
                        "supp_1": {
                            "type": "csv",
                            "filename": csv_path,
                            "statistical_analyses": [],
                        }
                    },
                }
            },
            f,
        )

    list_file = os.path.join(temp_dir, "paperconfig_files.txt")
    with open(list_file, "w") as f:
        f.write(config_with + "\n")
        f.write(config_without + "\n")

    adapter = MultiClusterAdapter(config_list_file=list_file)
    assert len(adapter.adapters) == 1
    assert adapter.adapters[0].doi == "10.1/with"


# ─── Test: New versioned extraction cache structure ───────────────────────────


def test_load_extraction_from_cache_dir(tmp_path):
    """_load_extraction_json reads stage files from .extraction_cache/{entry}/current/."""
    from multiomics_kg.extraction.cluster.run_manager import RunManager

    paperconfig_dir = tmp_path
    entry_key = "med4_kmeans_test"

    # Set up cache directory structure
    # RunManager expects cache_root = the .extraction_cache dir
    cache_dir = paperconfig_dir / ".extraction_cache"
    cache_dir.mkdir()
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.create_run()  # also sets current symlink

    # Write stage files
    stage2_data = {
        "1": {"id": "c1", "name": "Cluster 1", "functional_description": "Transport"},
        "2": {"id": "c2", "name": "Cluster 2", "functional_description": "Photosynthesis"},
    }
    stage3_data = {"1": {"verdict": "pass"}, "2": {"verdict": "fail"}}
    rm.write_stage(run_dir, 2, stage2_data)
    rm.write_stage(run_dir, 3, stage3_data)
    rm.finalize_run(run_dir)

    result = _load_extraction_json(paperconfig_dir, entry_key)
    assert "stage2_results" in result
    assert "stage3_validation" in result
    assert result["stage2_results"]["1"]["functional_description"] == "Transport"
    assert result["stage3_validation"]["2"]["verdict"] == "fail"
    # stage1 and stage4 not written, so not present
    assert "stage1_merged" not in result
    assert "stage4_review" not in result


def test_adapter_respects_review_reject():
    """stage4 reject overrides stage3 pass — cluster data should be empty."""
    extraction = {
        "stage2_results": {
            "1": {"id": "c1", "name": "Cluster 1", "functional_description": "Transport"},
        },
        "stage3_validation": {"1": {"verdict": "pass"}},
        "stage4_review": {"1": {"status": "reject", "reason": "incorrect extraction"}},
    }
    result = _get_extraction_cluster_data(extraction, "1")
    assert result == {}


def test_adapter_uses_edited_fields():
    """stage4 edit with edited_fields should override stage2 values."""
    extraction = {
        "stage2_results": {
            "1": {
                "id": "c1",
                "name": "Cluster 1",
                "functional_description": "Transport genes",
                "behavioral_description": "Upregulated early",
            },
        },
        "stage3_validation": {"1": {"verdict": "pass"}},
        "stage4_review": {
            "1": {
                "status": "edit",
                "edited_fields": {
                    "functional_description": "Nitrogen transport genes",
                },
            },
        },
    }
    result = _get_extraction_cluster_data(extraction, "1")
    assert result["functional_description"] == "Nitrogen transport genes"
    # Non-edited fields should retain stage2 values
    assert result["behavioral_description"] == "Upregulated early"
    assert result["name"] == "Cluster 1"
