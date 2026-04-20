"""
KG validity tests for ClusteringAnalysis nodes and related edges.

Validates:
- ClusteringAnalysis nodes exist with correct properties
- GeneCluster nodes exist with correct structure
- Edge structure: Publication → ClusteringAnalysis → GeneCluster → Gene
- Old edge types are removed
- Member counts and cluster counts are consistent
- Extraction-derived properties are populated where expected
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# ClusteringAnalysis node presence and counts
# ---------------------------------------------------------------------------

def test_clustering_analysis_nodes_exist(run_query):
    """ClusteringAnalysis nodes must exist in the graph."""
    result = run_query("MATCH (ca:ClusteringAnalysis) RETURN count(ca) AS cnt")
    assert result[0]["cnt"] >= 12, "Expected at least 12 ClusteringAnalysis nodes (across remaining papers after Biller 2018 retrofit removal)"


def test_gene_cluster_nodes_exist(run_query):
    """GeneCluster nodes must exist in the graph."""
    result = run_query("MATCH (gc:GeneCluster) RETURN count(gc) AS cnt")
    assert result[0]["cnt"] >= 90, "Expected at least 90 GeneCluster nodes (across remaining papers after Biller 2018 retrofit removal)"


# ---------------------------------------------------------------------------
# ClusteringAnalysis required properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prop", [
    "name",
    "organism_name",
    "cluster_method",
    "cluster_type",
    "cluster_count",
    "total_gene_count",
    "omics_type",
])
def test_clustering_analysis_required_properties(run_query, prop):
    """Every ClusteringAnalysis node must have key properties populated."""
    result = run_query(
        f"MATCH (ca:ClusteringAnalysis) WHERE ca.{prop} IS NULL RETURN count(ca) AS cnt"
    )
    assert result[0]["cnt"] == 0, f"Found ClusteringAnalysis nodes with null {prop}"


def test_clustering_analysis_cluster_type_values(run_query):
    """cluster_type must be from the valid enum."""
    valid_types = {
        "time_course", "diel", "condition_comparison", "classification",
    }
    result = run_query(
        "MATCH (ca:ClusteringAnalysis) RETURN DISTINCT ca.cluster_type AS ct"
    )
    for row in result:
        assert row["ct"] in valid_types, f"Invalid cluster_type: {row['ct']}"


# ---------------------------------------------------------------------------
# GeneCluster required properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prop", [
    "name",
    "organism_name",
    "member_count",
])
def test_gene_cluster_required_properties(run_query, prop):
    """Every GeneCluster node must have key properties populated."""
    result = run_query(
        f"MATCH (gc:GeneCluster) WHERE gc.{prop} IS NULL RETURN count(gc) AS cnt"
    )
    assert result[0]["cnt"] == 0, f"Found GeneCluster nodes with null {prop}"


def test_gene_cluster_no_source_paper(run_query):
    """source_paper property should not exist on GeneCluster (removed in redesign)."""
    result = run_query(
        "MATCH (gc:GeneCluster) WHERE gc.source_paper IS NOT NULL RETURN count(gc) AS cnt"
    )
    assert result[0]["cnt"] == 0, "Found GeneCluster nodes still carrying source_paper"


# ---------------------------------------------------------------------------
# Edge structure: new edges present
# ---------------------------------------------------------------------------

def test_publication_has_clustering_analysis_edges(run_query):
    """Publication → ClusteringAnalysis edges must exist."""
    result = run_query(
        "MATCH (:Publication)-[r:PublicationHasClusteringAnalysis]->(:ClusteringAnalysis) "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] >= 12, "Expected at least 12 Publication → ClusteringAnalysis edges"


def test_clustering_analysis_has_gene_cluster_edges(run_query):
    """ClusteringAnalysis → GeneCluster edges must exist."""
    result = run_query(
        "MATCH (:ClusteringAnalysis)-[r:ClusteringAnalysisHasGeneCluster]->(:GeneCluster) "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] >= 90, "Expected at least 90 ClusteringAnalysis → GeneCluster edges"


def test_clustering_analysis_belongs_to_organism_edges(run_query):
    """ClusteringAnalysis → OrganismTaxon edges must exist."""
    result = run_query(
        "MATCH (:ClusteringAnalysis)-[r:ClusteringanalysisBelongsToOrganism]->(:OrganismTaxon) "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] >= 12, "Expected at least 12 ClusteringAnalysis → OrganismTaxon edges"


def test_experiment_has_clustering_analysis_edges(run_query):
    """Experiment → ClusteringAnalysis edges must exist for linked experiments."""
    result = run_query(
        "MATCH (:Experiment)-[r:ExperimentHasClusteringAnalysis]->(:ClusteringAnalysis) "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] >= 10, "Expected at least 10 Experiment → ClusteringAnalysis edges"


def test_gene_in_gene_cluster_edges(run_query):
    """GeneCluster → Gene membership edges must exist."""
    result = run_query(
        "MATCH (:GeneCluster)-[r:Gene_in_gene_cluster]->(:Gene) RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] >= 2000, "Expected at least 2000 gene membership edges"


# ---------------------------------------------------------------------------
# Old edges removed
# ---------------------------------------------------------------------------

def test_no_old_publication_has_gene_cluster(run_query):
    """Old Publication_has_gene_cluster edges must not exist."""
    result = run_query(
        "MATCH ()-[r:Publication_has_gene_cluster]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, "Old Publication_has_gene_cluster edges still exist"


def test_no_old_genecluster_belongs_to_organism(run_query):
    """Old Genecluster_belongs_to_organism edges must not exist."""
    result = run_query(
        "MATCH ()-[r:Genecluster_belongs_to_organism]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, "Old Genecluster_belongs_to_organism edges still exist"


# ---------------------------------------------------------------------------
# Consistency checks
# ---------------------------------------------------------------------------

def test_cluster_count_matches_actual(run_query):
    """ClusteringAnalysis.cluster_count must match actual GeneCluster children."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)
        OPTIONAL MATCH (ca)-[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
        WITH ca, ca.cluster_count AS declared, count(gc) AS actual
        WHERE declared <> actual
        RETURN ca.name, declared, actual
    """)
    assert len(result) == 0, f"cluster_count mismatch: {result}"


def test_member_count_matches_actual(run_query):
    """GeneCluster.member_count must match actual Gene_in_gene_cluster edges."""
    result = run_query("""
        MATCH (gc:GeneCluster)
        OPTIONAL MATCH (gc)-[:Gene_in_gene_cluster]->(g:Gene)
        WITH gc, gc.member_count AS declared, count(g) AS actual
        WHERE declared <> actual
        RETURN gc.name, gc.organism_name, declared, actual
    """)
    assert len(result) == 0, f"member_count mismatch: {result}"


def test_total_gene_count_gte_member_sum(run_query):
    """ClusteringAnalysis.total_gene_count (from CSV) must be >= sum of member_count (actual edges).

    total_gene_count counts CSV rows; member_count counts resolved Gene edges.
    Some CSV gene IDs may not resolve, so total_gene_count >= sum(member_count).
    """
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
        WITH ca, ca.total_gene_count AS declared, sum(gc.member_count) AS actual
        WHERE declared < actual
        RETURN ca.name, declared, actual
    """)
    assert len(result) == 0, f"total_gene_count < sum(member_count), unexpected: {result}"


def test_every_gene_cluster_has_analysis_parent(run_query):
    """Every GeneCluster must be connected to exactly one ClusteringAnalysis."""
    result = run_query("""
        MATCH (gc:GeneCluster)
        OPTIONAL MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(gc)
        WITH gc, count(ca) AS parent_count
        WHERE parent_count <> 1
        RETURN gc.name, gc.organism_name, parent_count
    """)
    assert len(result) == 0, f"GeneCluster without exactly 1 ClusteringAnalysis parent: {result}"


def test_every_analysis_has_publication_parent(run_query):
    """Every ClusteringAnalysis must be connected to exactly one Publication."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)
        OPTIONAL MATCH (p:Publication)-[:PublicationHasClusteringAnalysis]->(ca)
        WITH ca, count(p) AS parent_count
        WHERE parent_count <> 1
        RETURN ca.name, parent_count
    """)
    assert len(result) == 0, f"ClusteringAnalysis without exactly 1 Publication parent: {result}"


def test_every_analysis_has_organism(run_query):
    """Every ClusteringAnalysis must be connected to exactly one OrganismTaxon."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)
        OPTIONAL MATCH (ca)-[:ClusteringanalysisBelongsToOrganism]->(o:OrganismTaxon)
        WITH ca, count(o) AS org_count
        WHERE org_count <> 1
        RETURN ca.name, org_count
    """)
    assert len(result) == 0, f"ClusteringAnalysis without exactly 1 organism: {result}"


def test_organism_name_consistency(run_query):
    """ClusteringAnalysis.organism_name must match the linked OrganismTaxon.preferred_name."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)-[:ClusteringanalysisBelongsToOrganism]->(o:OrganismTaxon)
        WHERE ca.organism_name <> o.preferred_name
        RETURN ca.name, ca.organism_name, o.preferred_name
    """)
    assert len(result) == 0, f"organism_name mismatch with OrganismTaxon: {result}"


# ---------------------------------------------------------------------------
# Denormalized property consistency
# ---------------------------------------------------------------------------

def test_gene_cluster_denormalized_organism(run_query):
    """GeneCluster.organism_name must match parent ClusteringAnalysis.organism_name."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
        WHERE gc.organism_name <> ca.organism_name
        RETURN gc.name, gc.organism_name, ca.organism_name
    """)
    assert len(result) == 0, f"Denormalized organism_name mismatch: {result}"


# ---------------------------------------------------------------------------
# Spot checks: Tolonen 2006
# ---------------------------------------------------------------------------

def test_tolonen_med4_analysis(run_query):
    """Tolonen 2006 MED4 clustering analysis should have 9 clusters, 410 genes."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)
        WHERE ca.name = 'MED4 K-means N-starvation clusters'
        RETURN ca.cluster_count AS clusters, ca.total_gene_count AS genes,
               ca.cluster_method AS method, ca.organism_name AS org
    """)
    assert len(result) == 1
    row = result[0]
    assert row["clusters"] == 9
    assert row["genes"] == 410
    assert row["method"] == "K-means (K=9)"
    assert row["org"] == "Prochlorococcus MED4"


def test_tolonen_mit9313_analysis(run_query):
    """Tolonen 2006 MIT9313 clustering analysis should have 7 clusters, 559 genes."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis)
        WHERE ca.name = 'MIT9313 K-means N-starvation clusters'
        RETURN ca.cluster_count AS clusters, ca.total_gene_count AS genes,
               ca.cluster_method AS method
    """)
    assert len(result) == 1
    row = result[0]
    assert row["clusters"] == 7
    assert row["genes"] == 559
    assert row["method"] == "K-means (K=7)"


def test_tolonen_med4_cluster6_largest(run_query):
    """MED4 cluster 6 (translation) should be the largest with 124 genes."""
    result = run_query("""
        MATCH (ca:ClusteringAnalysis {name: 'MED4 K-means N-starvation clusters'})
              -[:ClusteringAnalysisHasGeneCluster]->(gc:GeneCluster)
        RETURN gc.name, gc.member_count
        ORDER BY gc.member_count DESC LIMIT 1
    """)
    assert len(result) == 1
    assert result[0]["gc.member_count"] == 124


# ---------------------------------------------------------------------------
# Extraction-derived properties on GeneCluster
# ---------------------------------------------------------------------------

def test_some_clusters_have_functional_description(run_query):
    """Some GeneCluster nodes should have extraction-derived functional_description."""
    result = run_query("""
        MATCH (gc:GeneCluster)
        WHERE gc.functional_description IS NOT NULL AND gc.functional_description <> 'N/A'
        RETURN count(gc) AS cnt
    """)
    assert result[0]["cnt"] >= 30, "Expected at least 30 clusters with functional_description"


def test_some_clusters_have_temporal_pattern(run_query):
    """Some GeneCluster nodes should have extraction-derived temporal_pattern."""
    result = run_query("""
        MATCH (gc:GeneCluster)
        WHERE gc.temporal_pattern IS NOT NULL AND gc.temporal_pattern <> 'N/A'
        RETURN count(gc) AS cnt
    """)
    assert result[0]["cnt"] >= 30, "Expected at least 30 clusters with temporal_pattern"


def test_some_clusters_have_expression_dynamics(run_query):
    """Some GeneCluster nodes should have extraction-derived expression_dynamics."""
    result = run_query("""
        MATCH (gc:GeneCluster)
        WHERE gc.expression_dynamics IS NOT NULL AND gc.expression_dynamics <> 'N/A'
        RETURN count(gc) AS cnt
    """)
    assert result[0]["cnt"] >= 30, "Expected at least 30 clusters with expression_dynamics"


# ---------------------------------------------------------------------------
# Gene routing signals for clusters (post-import computed)
# ---------------------------------------------------------------------------

def test_gene_cluster_membership_count(run_query):
    """Genes in clusters should have cluster_membership_count > 0."""
    result = run_query("""
        MATCH (g:Gene) WHERE g.cluster_membership_count > 0
        RETURN count(g) AS cnt
    """)
    assert result[0]["cnt"] >= 5000, "Expected at least 5000 genes with cluster_membership_count > 0"


def test_gene_cluster_types_populated(run_query):
    """Genes in clusters should have cluster_types routing signal."""
    result = run_query("""
        MATCH (g:Gene) WHERE g.cluster_types IS NOT NULL AND size(g.cluster_types) > 0
        RETURN count(g) AS cnt
    """)
    assert result[0]["cnt"] >= 5000, "Expected at least 5000 genes with cluster_types populated"
