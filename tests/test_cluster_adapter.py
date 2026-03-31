"""
Tests for ClusterAdapter and MultiClusterAdapter.

Tests verify:
1. get_nodes emits correct GeneCluster nodes (count, IDs, properties)
2. get_edges emits membership edges with scores and correct gene IDs
3. Unresolved (NaN) genes are skipped
4. Configs without score_col produce edges with empty props
"""
import os
import tempfile

import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.cluster_adapter import (
    ClusterAdapter,
    MultiClusterAdapter,
    _make_cluster_id,
    _resolve_csv_path,
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
        data["locus_tag"] = ["PMM0001", None, "PMM0003", "PMM0004", "PMM0005"]
        data["resolution_method"] = ["tier1:locus_tag", None, "tier1:locus_tag", "tier1:locus_tag", "tier1:locus_tag"]

    df = pd.DataFrame(data)

    if with_resolved:
        stem = os.path.splitext(path)[0]
        resolved_path = f"{stem}_resolved.csv"
        if "locus_tag" not in df.columns:
            df["locus_tag"] = df["gene_id"]
            df["resolution_method"] = "tier1:locus_tag"
        df.to_csv(resolved_path, index=False)
        # Also write original so the config filename points to a real file
        df.drop(columns=["locus_tag", "resolution_method"], errors="ignore").to_csv(path, index=False)
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
                    "filename": csv_path,
                    "cluster_col": "cluster",
                    "gene_id_col": gene_id_col,
                    "score_col": score_col,
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "cluster_method": "k-means",
                    "treatment_type": ["nitrogen"],
                    "treatment": "Nitrogen limitation",
                    "light_condition": "continuous light",
                    "experimental_context": "in Pro99 medium",
                    "clusters": {
                        "1": {
                            "name": "Early response",
                            "cluster_type": "up",
                            "functional_description": "Photosynthesis genes",
                            "behavioral_description": "Peaks at 4h",
                            "peak_time_hours": 4.0,
                        },
                        "2": {
                            "name": "Late response",
                            "cluster_type": "down",
                            "functional_description": "Nitrogen assimilation",
                            "behavioral_description": "Sustained suppression",
                            "period_hours": 24.0,
                        },
                    },
                }
            },
        }
    }
    config_path = os.path.join(temp_dir, "paperconfig.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


# ─── Test 1: get_nodes emits correct GeneCluster nodes ───────────────────────


def test_get_nodes_emits_cluster_nodes(temp_dir):
    """Two clusters → two GeneCluster nodes with correct properties."""
    csv_path = os.path.join(temp_dir, "clusters.csv")
    _write_cluster_csv(csv_path, with_resolved=False)

    config_path = _write_paperconfig(temp_dir, csv_path)
    adapter = ClusterAdapter(config_file=config_path)
    nodes = adapter.get_nodes()

    assert len(nodes) == 2, f"Expected 2 nodes, got {len(nodes)}"

    ids = {n[0] for n in nodes}
    doi_short = "test.2024"
    expected_id_1 = f"cluster:{doi_short}:med4:1"
    expected_id_2 = f"cluster:{doi_short}:med4:2"
    assert expected_id_1 in ids, f"{expected_id_1} not in {ids}"
    assert expected_id_2 in ids, f"{expected_id_2} not in {ids}"

    # Check labels
    for node_id, label, props in nodes:
        assert label == "gene_cluster"

    # Check properties on cluster 1
    node_1 = next(n for n in nodes if n[0] == expected_id_1)
    props = node_1[2]
    assert props["name"] == "Early response"
    assert props["organism_name"] == "Prochlorococcus MED4"
    assert props["omics_type"] == "RNASEQ"
    assert props["member_count"] == 2  # PMM0001 and PMM0002
    assert props["treatment_type"] == ["nitrogen"]
    assert props["cluster_type"] == "up"
    assert props["peak_time_hours"] == 4.0

    # Check properties on cluster 2
    node_2 = next(n for n in nodes if n[0] == expected_id_2)
    props2 = node_2[2]
    assert props2["member_count"] == 3  # PMM0003, PMM0004, PMM0005
    assert props2["period_hours"] == 24.0


# ─── Test 2: get_edges emits membership edges with scores ────────────────────


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

    # Check publication edges exist
    pub_edges = [e for e in edges if e[3] == "publication_has_gene_cluster"]
    assert len(pub_edges) == 2, f"Expected 2 pub→cluster edges, got {len(pub_edges)}"
    for e in pub_edges:
        assert e[1] == "doi:10.1234/test.2024"


# ─── Test 3: Unresolved (NaN) genes are skipped ──────────────────────────────


def test_nan_locus_tags_skipped(temp_dir):
    """NaN locus_tag values produce no membership edges."""
    csv_path = os.path.join(temp_dir, "clusters.csv")

    # Build a resolved CSV manually with one NaN locus_tag
    df = pd.DataFrame({
        "gene_id": ["PMM0001", "PMM0002", "PMM0003", "PMM0004", "PMM0005"],
        "cluster": [1, 1, 2, 2, 2],
        "score": [0.9, 0.8, 0.7, 0.85, 0.6],
        "locus_tag": ["PMM0001", None, "PMM0003", "PMM0004", "PMM0005"],
        "resolution_method": [
            "tier1:locus_tag", None, "tier1:locus_tag", "tier1:locus_tag", "tier1:locus_tag"
        ],
    })
    # Write original (required by config filename)
    df.drop(columns=["locus_tag", "resolution_method"]).to_csv(csv_path, index=False)
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


# ─── Test 4: No score_col means edges have empty props ───────────────────────


def test_no_score_col_produces_empty_edge_props(temp_dir):
    """Without score_col the membership edge props dict should be empty."""
    csv_path = os.path.join(temp_dir, "clusters_noscore.csv")

    # Write a resolved CSV (no score column)
    df = pd.DataFrame({
        "gene_id": ["PMM0001", "PMM0002"],
        "cluster": [1, 1],
        "locus_tag": ["PMM0001", "PMM0002"],
        "resolution_method": ["tier1:locus_tag", "tier1:locus_tag"],
    })
    df.drop(columns=["locus_tag", "resolution_method"]).to_csv(csv_path, index=False)
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
                    "treatment_type": [],
                    "clusters": {
                        "1": {"name": "Cluster A"},
                    },
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


# ─── Test 5: _make_cluster_id fallback when no DOI ───────────────────────────


def test_make_cluster_id_no_doi():
    """Without DOI, falls back to paper_name slug + organism suffix."""
    cid = _make_cluster_id("", "Test Paper 2024", "Prochlorococcus MED4", "cluster_A")
    assert cid == "cluster:test_paper_2024:med4:cluster_A"


def test_make_cluster_id_with_doi():
    """With DOI, uses last path component + organism suffix."""
    cid = _make_cluster_id("10.1234/my.paper", "Test Paper 2024", "Prochlorococcus MED4", "k1")
    assert cid == "cluster:my.paper:med4:k1"


# ─── Test 6: _resolve_csv_path prefers _resolved.csv ────────────────────────


def test_resolve_csv_path_prefers_resolved(temp_dir):
    original = os.path.join(temp_dir, "data.csv")
    resolved = os.path.join(temp_dir, "data_resolved.csv")

    # Only original exists
    open(original, "w").close()
    assert str(_resolve_csv_path(original)) == original

    # Now create resolved — should be preferred
    open(resolved, "w").close()
    assert str(_resolve_csv_path(original)) == resolved


# ─── Test 7: MultiClusterAdapter skips configs without gene_clusters ─────────


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
