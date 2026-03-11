"""
Post-import Cypher validation tests.

The post-import.cypher script runs after neo4j-admin import to:
1. Denormalize function_description from Protein onto Gene (uniform functional
   text access for all organisms without a multi-hop traversal)
2. Denormalize go_biological_processes from Protein onto Gene (where Gene
   doesn't already have it from Cyanorak)
3. Create bidirectional Gene_is_homolog_of_gene edges between genes sharing
   a Cyanorak_cluster (both A→B and B→A, with distance property)
4. Propagate expression ortholog edges: if X changes expression of gene A,
   and A is homolog of B, then X changes_expression_of_ortholog B
   (Condition_changes_expression_of_ortholog and Coculture_changes_expression_of_ortholog)

These tests verify that post-import ran correctly and produced valid results.

Note: test_gene_go_biological_processes_set was removed (Feb 2026).
The Cypher denormalization step for go_biological_processes was never implemented
in post-import.cypher (only function_description was added in commit 0a9a52c), and
the go_biological_processes field was subsequently removed from the gene schema in
the schema refactor (commit 9fdd673). If go_biological_processes denormalization is
needed in future, add it to post-import.cypher and the gene schema before re-adding
the test.
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Denormalized Gene.function_description (copied from Protein)
# ---------------------------------------------------------------------------

def test_gene_function_description_set(run_query):
    """
    After post-import, genes linked to a Protein with function_description
    should have Gene.function_description populated.
    Allow up to 10% ungapped (some proteins may have null function_description
    in UniProt; the threshold catches a complete post-import failure).
    """
    result = run_query("""
        MATCH (p:Protein)-[:Gene_encodes_protein]->(g:Gene)
        WHERE p.function_description IS NOT NULL
        WITH count(g) AS total,
             count(CASE WHEN g.function_description IS NOT NULL THEN 1 END) AS filled
        RETURN total, filled, total - filled AS missing
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Protein→Gene pairs with function_description found")
    missing_fraction = row["missing"] / row["total"]
    assert missing_fraction < 0.10, (
        f"{row['missing']} / {row['total']} genes ({missing_fraction:.1%}) are missing "
        f"function_description despite their linked Protein having it. "
        f"Did the denormalization step in post-import.cypher run?"
    )


def test_alteromonas_genes_have_function_description(run_query):
    """
    Alteromonas genes have no Cyanorak annotations, so function_description
    must come entirely from the post-import denormalization step.
    At least 50% of Alteromonas genes linked to a reviewed Protein should
    have function_description set.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.genus = 'Alteromonas'
        WITH count(g) AS total,
             count(CASE WHEN g.function_description IS NOT NULL THEN 1 END) AS filled
        RETURN total, filled
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Alteromonas Gene nodes found")
    filled_fraction = row["filled"] / row["total"]
    assert filled_fraction > 0.50, (
        f"Only {row['filled']} / {row['total']} Alteromonas genes ({filled_fraction:.1%}) "
        f"have function_description. Denormalization from Protein may not have run."
    )


# ---------------------------------------------------------------------------
# Gene_is_homolog_of_gene edges
# ---------------------------------------------------------------------------

def test_homolog_edges_exist(run_query):
    """Gene_is_homolog_of_gene edges must exist after post-import."""
    result = run_query(
        "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No Gene_is_homolog_of_gene edges found. "
        "Did post-import.cypher run successfully?"
    )


def test_homolog_edge_count_reasonable(run_query):
    """
    With 12 strains and thousands of shared clusters, expect millions of
    homolog edge pairs. A low count indicates post-import failure.
    """
    result = run_query(
        "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 100_000, (
        f"Only {cnt} Gene_is_homolog_of_gene edges found; "
        f"expected > 100,000 for 12 strains with shared Cyanorak clusters"
    )


def test_homolog_edges_are_bidirectional(run_query):
    """
    post-import.cypher creates both A→B and B→A homolog edges so LLM agents
    can traverse from any direction. For a sample of 200 edges, verify that
    the reverse edge always exists.
    """
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_is_homolog_of_gene]->(g2:Gene)
        WITH g1, g2 LIMIT 200
        WHERE NOT (g2)-[:Gene_is_homolog_of_gene]->(g1)
        RETURN count(*) AS missing_reverse
    """)
    missing = result[0]["missing_reverse"]
    assert missing == 0, (
        f"{missing} forward homolog edge(s) have no corresponding reverse edge. "
        f"post-import.cypher bidirectional MERGE may have failed."
    )


def test_homolog_edges_have_distance_property(run_query):
    """Every Gene_is_homolog_of_gene edge must carry a distance property."""
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.distance IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Gene_is_homolog_of_gene edges are missing the distance property"
    )


def test_homolog_edges_have_cluster_id(run_query):
    """Every homolog edge must reference the Cyanorak cluster that caused it."""
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.cluster_id IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Gene_is_homolog_of_gene edges are missing cluster_id"
    )


def test_homolog_edges_source_is_known(run_query):
    """
    All Gene_is_homolog_of_gene edges must have a known source value.
    Valid sources: cyanorak_cluster (from post-import.cypher) and eggNOG
    orthologous-group sources (eggnog_alteromonadaceae_og, eggnog_bacteria_cog_og).
    """
    valid_sources = {
        "cyanorak_cluster",
        "eggnog_alteromonadaceae_og",
        "eggnog_bacteria_cog_og",
    }
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.source IS NULL OR NOT r.source IN $valid
        RETURN count(r) AS bad_rows,
               collect(DISTINCT r.source)[..5] AS bad_sources
    """, valid=list(valid_sources))
    bad = result[0]["bad_rows"]
    assert bad == 0, (
        f"{bad} Gene_is_homolog_of_gene edges have unexpected source: "
        f"{result[0]['bad_sources']}"
    )


def test_homolog_distance_valid_values(run_query):
    """
    distance values are produced by a CASE statement in post-import.cypher
    and by the eggNOG homolog adapter.
    All values must start with one of the expected prefixes.
    """
    valid_prefixes = [
        "same strain:",
        "same clade:",
        "same species:",
        "same genus:",
        "same family:",
        "same order:",
        "cross order",
        "cross phylum",
    ]
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.distance IS NOT NULL
        RETURN collect(DISTINCT r.distance) AS all_distances
    """)
    if not result or not result[0]["all_distances"]:
        pytest.skip("No distance values to check")

    all_distances = result[0]["all_distances"]
    bad = [
        d for d in all_distances
        if not any(d.startswith(p) for p in valid_prefixes)
    ]
    assert not bad, (
        f"Unexpected distance values found: {bad}\n"
        f"Valid prefixes: {valid_prefixes}"
    )


def test_no_self_homolog_edges(run_query):
    """A gene should not be its own homolog (same node → same node)."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_is_homolog_of_gene]->(g)
        RETURN count(g) AS self_loops
    """)
    self_loops = result[0]["self_loops"]
    assert self_loops == 0, (
        f"{self_loops} Gene node(s) have self-referential Gene_is_homolog_of_gene edges"
    )


# ---------------------------------------------------------------------------
# eggNOG-sourced homolog edges
# ---------------------------------------------------------------------------

def test_eggnog_alteromonadaceae_homolog_edges_exist(run_query):
    """
    eggnog_alteromonadaceae_og homolog edges must exist, connecting the three
    Alteromonas strains (MIT1002, EZ55, HOT1A3) via shared family-level OGs.
    """
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.source = 'eggnog_alteromonadaceae_og'
        RETURN count(r) AS cnt
    """)
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No eggnog_alteromonadaceae_og Gene_is_homolog_of_gene edges found"
    )


def test_eggnog_bacteria_cog_homolog_edges_exist(run_query):
    """
    eggnog_bacteria_cog_og homolog edges must exist, connecting genes across
    all organisms via conserved bacterial COG orthologous groups.
    """
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.source = 'eggnog_bacteria_cog_og'
        RETURN count(r) AS cnt
    """)
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No eggnog_bacteria_cog_og Gene_is_homolog_of_gene edges found"
    )


def test_eggnog_alteromonadaceae_homologs_connect_alteromonas(run_query):
    """
    Alteromonadaceae-level OG edges should connect genes that belong to
    Alteromonas organisms (same family). At least one such cross-strain pair
    must exist.
    """
    result = run_query("""
        MATCH (g1:Gene)-[r:Gene_is_homolog_of_gene]->(g2:Gene)
        WHERE r.source = 'eggnog_alteromonadaceae_og'
          AND r.distance = 'same family: Alteromonadaceae'
        MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
        MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
        WHERE o1.strain_name <> o2.strain_name
        RETURN count(r) AS cnt
    """)
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No cross-strain Alteromonadaceae homolog edges found. "
        "Expected edges connecting MIT1002, EZ55, and HOT1A3."
    )


def test_eggnog_bacteria_cog_connects_cross_phylum(run_query):
    """
    Bacteria COG OG edges must connect genes from different genera
    (e.g., Alteromonas ↔ Prochlorococcus). All these edges should carry
    distance='cross phylum'.
    """
    result = run_query("""
        MATCH (g1:Gene)-[r:Gene_is_homolog_of_gene]->(g2:Gene)
        WHERE r.source = 'eggnog_bacteria_cog_og'
        MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
        MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
        WHERE o1.genus <> o2.genus
        RETURN count(r) AS cnt
    """)
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No cross-genus eggnog_bacteria_cog_og homolog edges found. "
        "Expected edges connecting Alteromonas genes to Prochlorococcus/Synechococcus genes."
    )


def test_eggnog_homologs_have_cluster_id(run_query):
    """
    All eggNOG-sourced homolog edges must carry cluster_id (the OG identifier
    like 'COG0168' or '46764@72275,Alteromonadaceae').
    """
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.source IN ['eggnog_alteromonadaceae_og', 'eggnog_bacteria_cog_og']
          AND r.cluster_id IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} eggNOG-sourced Gene_is_homolog_of_gene edges are missing cluster_id"
    )


def test_eggnog_homologs_are_bidirectional(run_query):
    """
    eggNOG homolog edges must be bidirectional (A→B implies B→A),
    consistent with Cyanorak-derived homolog edges.
    """
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_is_homolog_of_gene]->(g2:Gene)
        WHERE EXISTS {
            MATCH ()-[r:Gene_is_homolog_of_gene]->()
            WHERE r.source IN ['eggnog_alteromonadaceae_og', 'eggnog_bacteria_cog_og']
              AND startNode(r) = g1 AND endNode(r) = g2
        }
        WITH g1, g2 LIMIT 500
        WHERE NOT (g2)-[:Gene_is_homolog_of_gene]->(g1)
        RETURN count(*) AS missing_reverse
    """)
    missing = result[0]["missing_reverse"]
    assert missing == 0, (
        f"{missing} eggNOG homolog edge(s) have no corresponding reverse edge"
    )


# ---------------------------------------------------------------------------
# Expression ortholog edges (Condition_changes_expression_of_ortholog and
# Coculture_changes_expression_of_ortholog)
# ---------------------------------------------------------------------------

def test_expression_homolog_edges_exist(run_query):
    """Expression ortholog edges must exist after post-import."""
    result = run_query(
        "MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No Condition_changes_expression_of_ortholog or Coculture_changes_expression_of_ortholog edges found. "
        "Did post-import.cypher propagation step run?"
    )


def test_expression_homolog_has_original_gene(run_query):
    """Every expression ortholog edge must record the original_gene ID."""
    result = run_query("""
        MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->()
        WHERE r.original_gene IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} expression ortholog edges are missing original_gene"
    )


def test_expression_homolog_has_homology_metadata(run_query):
    """
    Expression ortholog edges must carry homology_source,
    homology_cluster_id, and distance (copied from the homolog edge).
    """
    result = run_query("""
        MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->()
        WHERE r.homology_source IS NULL
           OR r.homology_cluster_id IS NULL
           OR r.distance IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} expression ortholog edges are missing homology metadata "
        f"(homology_source / homology_cluster_id / distance)"
    )


def test_expression_homolog_inherits_log2fc(run_query):
    """
    Expression ortholog edges must carry log2_fold_change (always present on the
    original edge). Note: adjusted_p_value may be null if the original study
    did not report it — that is expected and not tested here.
    """
    result = run_query("""
        MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->()
        WHERE r.log2_fold_change IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} expression ortholog edges are missing log2_fold_change"
    )


def test_expression_homolog_count_vs_direct(run_query):
    """
    The number of expression ortholog edges should be substantially larger
    than the number of direct expression edges (homology amplifies coverage).
    If ortholog_count < direct_count, propagation likely failed.
    """
    result = run_query("""
        MATCH ()-[:Condition_changes_expression_of|Coculture_changes_expression_of]->()
        WITH count(*) AS direct
        MATCH ()-[:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->()
        WITH direct, count(*) AS homolog
        RETURN direct, homolog
    """)
    row = result[0]
    assert row["homolog"] > row["direct"], (
        f"Expression ortholog edge count ({row['homolog']}) is not greater "
        f"than direct expression edge count ({row['direct']}). "
        f"Ortholog propagation may have failed or graph has very few shared clusters."
    )


def test_no_cross_phylum_coculture_ortholog_edges(run_query):
    """
    Coculture_changes_expression_of_ortholog edges must not span cross-phylum
    distances. Coculture experiments are intra-genus (Prochlorococcus/Alteromonas),
    so cross-phylum propagation would be biologically nonsensical.
    """
    result = run_query("""
        MATCH ()-[e:Coculture_changes_expression_of_ortholog]->()
        WHERE e.distance = 'cross phylum'
        RETURN count(e) AS cnt
    """)
    cnt = result[0]["cnt"]
    assert cnt == 0, (
        f"{cnt} Coculture_changes_expression_of_ortholog edges have distance='cross phylum'. "
        f"Cross-phylum coculture ortholog edges should not exist."
    )
