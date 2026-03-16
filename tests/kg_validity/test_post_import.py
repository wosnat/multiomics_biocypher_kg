"""
Post-import validation tests.

The post-import.sh script runs after neo4j-admin import to create indexes.
Homolog edges and expression propagation have been replaced by OrthologGroup
nodes and Gene_in_ortholog_group edges (adapter-created, query-time joins).

These tests verify:
1. OrthologGroup nodes exist with correct properties
2. Gene_in_ortholog_group membership edges exist
3. No old homolog/ortholog edges remain
4. No old Cyanorak_cluster nodes remain
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# OrthologGroup nodes
# ---------------------------------------------------------------------------

def test_ortholog_group_nodes_exist(run_query):
    """OrthologGroup nodes must exist."""
    result = run_query("MATCH (og:OrthologGroup) RETURN count(og) AS cnt")
    assert result[0]["cnt"] > 5000, (
        f"Only {result[0]['cnt']} OrthologGroup nodes found; expected > 5000"
    )


def test_ortholog_group_sources(run_query):
    """OrthologGroup source must be 'cyanorak' or 'eggnog'."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        WHERE NOT og.source IN ['cyanorak', 'eggnog']
        RETURN count(og) AS bad, collect(DISTINCT og.source)[..5] AS bad_sources
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} OrthologGroup nodes have unexpected source: {result[0]['bad_sources']}"
    )


def test_ortholog_group_per_source_counts(run_query):
    """Each source must have a minimum number of nodes (catches silent zero)."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        RETURN og.source AS source, count(og) AS cnt
        ORDER BY cnt DESC
    """)
    counts = {row["source"]: row["cnt"] for row in result}
    assert counts.get("cyanorak", 0) > 5000, (
        f"Cyanorak OrthologGroup count {counts.get('cyanorak', 0)} < 5000"
    )
    assert counts.get("eggnog", 0) > 2000, (
        f"EggNOG OrthologGroup count {counts.get('eggnog', 0)} < 2000"
    )


def test_ortholog_group_has_taxonomic_level(run_query):
    """All OrthologGroup nodes must have a non-null taxonomic_level."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        WHERE og.taxonomic_level IS NULL
        RETURN count(og) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} OrthologGroup nodes are missing taxonomic_level"
    )


# ---------------------------------------------------------------------------
# Gene_in_ortholog_group edges
# ---------------------------------------------------------------------------

def test_gene_in_ortholog_group_edges_exist(run_query):
    """Membership edges must exist."""
    result = run_query(
        "MATCH ()-[e:Gene_in_ortholog_group]->() RETURN count(e) AS cnt"
    )
    assert result[0]["cnt"] > 50000, (
        f"Only {result[0]['cnt']} Gene_in_ortholog_group edges; expected > 50000"
    )


def test_gene_in_ortholog_group_per_source_counts(run_query):
    """Each source must have a minimum number of membership edges."""
    result = run_query("""
        MATCH ()-[e:Gene_in_ortholog_group]->(og:OrthologGroup)
        RETURN og.source AS source, count(e) AS cnt
        ORDER BY cnt DESC
    """)
    counts = {row["source"]: row["cnt"] for row in result}
    assert counts.get("cyanorak", 0) > 18000, (
        f"Cyanorak membership edges {counts.get('cyanorak', 0)} < 18000"
    )
    assert counts.get("eggnog", 0) > 24000, (
        f"EggNOG membership edges {counts.get('eggnog', 0)} < 24000"
    )


# ---------------------------------------------------------------------------
# Cross-strain bridging via OrthologGroup
# ---------------------------------------------------------------------------

def test_cyanorak_og_pro_genes_share_group(run_query):
    """Two Pro genes from different strains with same Cyanorak cluster → shared OG."""
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'cyanorak'})<-[:Gene_in_ortholog_group]-(g2:Gene)
        WHERE g1 <> g2
        MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
        MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
        WHERE o1.strain_name <> o2.strain_name
          AND o1.genus = 'Prochlorococcus'
          AND o2.genus = 'Prochlorococcus'
        RETURN count(*) AS cnt LIMIT 1
    """)
    assert result[0]["cnt"] > 0, (
        "No cross-strain Prochlorococcus gene pairs sharing a Cyanorak OrthologGroup"
    )


def test_eggnog_lowest_alt_genes_share_group(run_query):
    """Two Alt genes from different strains sharing an Alteromonadaceae-level OG."""
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'eggnog', taxonomic_level: 'Alteromonadaceae'})<-[:Gene_in_ortholog_group]-(g2:Gene)
        WHERE g1 <> g2
        MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
        MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
        WHERE o1.strain_name <> o2.strain_name
          AND o1.genus = 'Alteromonas'
          AND o2.genus = 'Alteromonas'
        RETURN count(*) AS cnt LIMIT 1
    """)
    assert result[0]["cnt"] > 0, (
        "No cross-strain Alteromonas gene pairs sharing an Alteromonadaceae OrthologGroup"
    )


def test_cross_phylum_cog_bridging(run_query):
    """Pro gene + Alt gene sharing a Bacteria-level COG → same OG."""
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'eggnog', taxonomic_level: 'Bacteria'})<-[:Gene_in_ortholog_group]-(g2:Gene)
        WHERE g1 <> g2
        MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
        MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
        WHERE o1.genus <> o2.genus
        RETURN count(*) AS cnt LIMIT 1
    """)
    assert result[0]["cnt"] > 0, (
        "No cross-genus gene pairs sharing a Bacteria-level COG OrthologGroup"
    )


# ---------------------------------------------------------------------------
# Absence of old edge/node types
# ---------------------------------------------------------------------------

def test_no_old_homolog_edges(run_query):
    """Gene_is_homolog_of_gene edges must not exist (replaced by OrthologGroup)."""
    result = run_query(
        "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} Gene_is_homolog_of_gene edges still exist"
    )


def test_no_old_ortholog_expression_edges(run_query):
    """*_expression_of_ortholog edges must not exist (replaced by query-time joins)."""
    result = run_query(
        "MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->() "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} *_expression_of_ortholog edges still exist"
    )


def test_no_cyanorak_cluster_nodes(run_query):
    """Cyanorak_cluster nodes must not exist (replaced by OrthologGroup)."""
    result = run_query("MATCH (c:Cyanorak_cluster) RETURN count(c) AS cnt")
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} Cyanorak_cluster nodes still exist"
    )


def test_no_gene_in_cyanorak_cluster_edges(run_query):
    """Gene_in_cyanorak_cluster edges must not exist."""
    result = run_query(
        "MATCH ()-[r:Gene_in_cyanorak_cluster]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} Gene_in_cyanorak_cluster edges still exist"
    )


# ---------------------------------------------------------------------------
# OrthologGroup enrichment properties
# ---------------------------------------------------------------------------

def test_ortholog_group_has_member_count(run_query):
    """All OrthologGroup nodes must have member_count > 0."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        WHERE og.member_count IS NULL OR og.member_count < 1
        RETURN count(og) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} OrthologGroup nodes have null or zero member_count"
    )


def test_ortholog_group_has_organism_count(run_query):
    """All OrthologGroup nodes must have organism_count > 0."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        WHERE og.organism_count IS NULL OR og.organism_count < 1
        RETURN count(og) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} OrthologGroup nodes have null or zero organism_count"
    )


def test_ortholog_group_member_count_matches_edges(run_query):
    """member_count should match actual Gene_in_ortholog_group edge count."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        OPTIONAL MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og)
        WITH og, og.member_count AS declared, count(g) AS actual
        WHERE declared <> actual
        RETURN count(og) AS mismatched,
               collect(og.name)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"{result[0]['mismatched']} OrthologGroup nodes have member_count != edge count. "
        f"Examples: {result[0]['examples']}"
    )


def test_ortholog_group_cross_genus_flag_consistent(run_query):
    """has_cross_genus_members should match actual genera from member genes."""
    result = run_query("""
        MATCH (og:OrthologGroup)<-[:Gene_in_ortholog_group]-(g:Gene)
        WITH og, collect(DISTINCT split(g.organism_strain, ' ')[0]) AS actual_genera
        WHERE (og.has_cross_genus_members = 'cross_genus' AND size(actual_genera) < 2)
           OR (og.has_cross_genus_members = 'single_genus' AND size(actual_genera) > 1)
        RETURN count(og) AS inconsistent,
               collect(og.name)[..5] AS examples
    """)
    assert result[0]["inconsistent"] == 0, (
        f"{result[0]['inconsistent']} OrthologGroup nodes have inconsistent has_cross_genus_members. "
        f"Examples: {result[0]['examples']}"
    )


def test_ortholog_group_consensus_product_coverage(run_query):
    """Majority of OrthologGroup nodes should have a consensus_product."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        WITH count(og) AS total,
             count(CASE WHEN og.consensus_product IS NOT NULL THEN 1 END) AS with_product
        RETURN total, with_product,
               toFloat(with_product) / total AS coverage
    """)
    coverage = result[0]["coverage"]
    assert coverage > 0.5, (
        f"Only {coverage:.1%} of OrthologGroup nodes have consensus_product; expected > 50%"
    )


def test_specificity_rank_index_exists(run_query):
    """The specificity_rank index should exist."""
    result = run_query("SHOW INDEXES YIELD name WHERE name = 'ortholog_group_rank_idx' RETURN count(*) AS cnt")
    assert result[0]["cnt"] == 1, "ortholog_group_rank_idx index not found"


# ---------------------------------------------------------------------------
# Denormalized Gene.function_description (copied from Protein) — kept from
# original post-import tests, still validated
# ---------------------------------------------------------------------------

def test_gene_function_description_set(run_query):
    """
    After post-import, genes linked to a Protein with function_description
    should have Gene.function_description populated.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_encodes_protein]->(p:Protein)
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
        f"function_description despite their linked Protein having it."
    )
