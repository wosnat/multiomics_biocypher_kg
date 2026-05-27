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
        WITH og, collect(DISTINCT split(g.organism_name, ' ')[0]) AS actual_genera
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


def test_ortholog_group_fulltext_index_exists(run_query):
    """The orthologGroupFullText index should exist."""
    result = run_query(
        "SHOW INDEXES YIELD name, type WHERE name = 'orthologGroupFullText' RETURN count(*) AS cnt"
    )
    assert result[0]["cnt"] == 1, "orthologGroupFullText index not found"


def test_ortholog_group_fulltext_search_returns_results(run_query):
    """Fulltext search should return results for 'photosynthesis'."""
    result = run_query("""
        CALL db.index.fulltext.queryNodes('orthologGroupFullText', 'photosynthesis')
        YIELD node, score
        RETURN count(node) AS matches
    """)
    assert result[0]["matches"] > 0, (
        "orthologGroupFullText search for 'photosynthesis' returned no results"
    )


def test_ortholog_group_description_coverage(run_query):
    """EggNOG OrthologGroup nodes should have description populated (majority)."""
    result = run_query("""
        MATCH (og:OrthologGroup {source: 'eggnog'})
        RETURN count(og) AS total,
               count(og.description) AS with_description,
               toFloat(count(og.description)) / count(og) AS coverage
    """)
    coverage = result[0]["coverage"]
    assert coverage > 0.5, (
        f"Only {coverage:.1%} of eggNOG OrthologGroup nodes have description; expected > 50%"
    )


def test_ortholog_group_functional_description_coverage(run_query):
    """Some OrthologGroup nodes should have functional_description."""
    result = run_query("""
        MATCH (og:OrthologGroup)
        RETURN count(og) AS total,
               count(og.functional_description) AS with_func_desc
    """)
    assert result[0]["with_func_desc"] > 100, (
        f"Only {result[0]['with_func_desc']} OrthologGroup nodes have "
        f"functional_description; expected > 100"
    )


def test_og_in_cog_category_edges_exist(run_query):
    """Og_in_cog_category edges should exist."""
    result = run_query(
        "MATCH ()-[r:Og_in_cog_category]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] > 100, (
        f"Only {result[0]['cnt']} Og_in_cog_category edges; expected > 100"
    )


def test_og_has_cyanorak_role_edges_exist(run_query):
    """Og_has_cyanorak_role edges should exist."""
    result = run_query(
        "MATCH ()-[r:Og_has_cyanorak_role]->() RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] > 100, (
        f"Only {result[0]['cnt']} Og_has_cyanorak_role edges; expected > 100"
    )


def test_og_in_cog_category_targets_valid(run_query):
    """Og_in_cog_category edges must point to CogFunctionalCategory nodes."""
    result = run_query("""
        MATCH (og:OrthologGroup)-[:Og_in_cog_category]->(c)
        WHERE NOT c:CogFunctionalCategory
        RETURN count(*) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Og_in_cog_category edges point to non-CogFunctionalCategory nodes"
    )


def test_og_has_cyanorak_role_targets_valid(run_query):
    """Og_has_cyanorak_role edges must point to CyanorakRole nodes."""
    result = run_query("""
        MATCH (og:OrthologGroup)-[:Og_has_cyanorak_role]->(r)
        WHERE NOT r:CyanorakRole
        RETURN count(*) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Og_has_cyanorak_role edges point to non-CyanorakRole nodes"
    )


# ---------------------------------------------------------------------------
# Denormalized Gene.function_description (copied from Protein) — kept from
# original post-import tests, still validated
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Gene routing signals (pre-computed by post-import Cypher)
# ---------------------------------------------------------------------------

def test_annotation_types_populated(run_query):
    """All genes should have annotation_types (list, possibly empty)."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.annotation_types IS NULL
        RETURN count(g) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} genes missing annotation_types property"
    )


def test_annotation_types_values_valid(run_query):
    """annotation_types should only contain known ontology type labels."""
    result = run_query("""
        MATCH (g:Gene)
        UNWIND g.annotation_types AS t
        WITH DISTINCT t
        WHERE NOT t IN ['go_bp', 'go_mf', 'go_cc', 'pfam', 'cog_category',
                        'kegg', 'brite', 'ec', 'cyanorak_role', 'tigr_role',
                        'tcdb', 'cazy']
        RETURN collect(t) AS bad
    """)
    assert result[0]["bad"] == [], (
        f"Unexpected annotation_types values: {result[0]['bad']}"
    )


def test_annotation_types_spot_check(run_query):
    """A gene with known GO edges should have go_bp or go_mf in annotation_types."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_involved_in_biological_process]->()
        WITH g LIMIT 1
        RETURN g.locus_tag AS lt, g.annotation_types AS types
    """)
    if not result:
        pytest.skip("No Gene_involved_in_biological_process edges found")
    assert 'go_bp' in result[0]["types"], (
        f"Gene {result[0]['lt']} has GO BP edges but annotation_types = {result[0]['types']}"
    )


def test_expression_edge_count_populated(run_query):
    """All genes should have expression_edge_count (int, possibly 0)."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.expression_edge_count IS NULL
        RETURN count(g) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} genes missing expression_edge_count"
    )


def test_expression_edge_count_accurate(run_query):
    """expression_edge_count should match actual edge count for a sample gene."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.expression_edge_count > 0
        WITH g LIMIT 5
        OPTIONAL MATCH (g)<-[e:Changes_expression_of]-()
        WITH g, g.expression_edge_count AS declared, count(e) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag) AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"expression_edge_count mismatch for: {result[0]['examples']}"
    )


def test_significant_counts_lte_total(run_query):
    """significant_up_count + significant_down_count must be <= expression_edge_count."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE (g.significant_up_count + g.significant_down_count) > g.expression_edge_count
        RETURN count(g) AS bad, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} genes have sig_up + sig_down > total: {result[0]['examples']}"
    )


def test_closest_ortholog_group_size_spot_check(run_query):
    """Genes with OG edges should have closest_ortholog_group_size > 0."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_in_ortholog_group]->()
        WITH DISTINCT g LIMIT 10
        WHERE g.closest_ortholog_group_size IS NULL OR g.closest_ortholog_group_size < 1
        RETURN count(g) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"Genes with OG edges but null/zero closest_ortholog_group_size"
    )


def test_closest_ortholog_genera_is_list(run_query):
    """closest_ortholog_genera should be a list where present."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.closest_ortholog_group_size IS NOT NULL
          AND g.closest_ortholog_group_size > 0
          AND g.closest_ortholog_genera IS NULL
        RETURN count(g) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} genes have OG size but null genera"
    )


# ---------------------------------------------------------------------------
# Denormalized Gene.function_description (copied from Protein)
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


# ---------------------------------------------------------------------------
# Publication computed properties (pre-computed by post-import Cypher)
# ---------------------------------------------------------------------------

def test_publication_fulltext_index_exists(run_query):
    """publicationFullText index must exist."""
    result = run_query(
        "SHOW INDEXES YIELD name WHERE name = 'publicationFullText' RETURN count(*) AS cnt"
    )
    assert result[0]["cnt"] == 1, "publicationFullText index not found"


def test_publication_experiment_count_populated(run_query):
    """All Publications should have experiment_count (not null)."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE p.experiment_count IS NULL
        RETURN count(p) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} publications missing experiment_count"
    )


def test_publication_experiment_count_accurate(run_query):
    """experiment_count should match actual Has_experiment edge count."""
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
        WITH p, p.experiment_count AS declared, count(e) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"experiment_count mismatch for: {result[0]['examples']}"
    )


def test_publication_treatment_types_no_nulls(run_query):
    """treatment_types list should contain no null entries."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE ANY(x IN p.treatment_types WHERE x IS NULL)
        RETURN count(p) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} publications have null entries in treatment_types"
    )


def test_publication_background_factors_no_nulls(run_query):
    """background_factors list should contain no null entries."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE p.background_factors IS NOT NULL
          AND ANY(x IN p.background_factors WHERE x IS NULL)
        RETURN count(p) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} publications have null entries in background_factors"
    )


def test_publication_omics_types_no_nulls(run_query):
    """omics_types list should contain no null entries."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE ANY(x IN p.omics_types WHERE x IS NULL)
        RETURN count(p) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} publications have null entries in omics_types"
    )


def test_publication_organisms_from_experiments(run_query):
    """Publications with experiments should have non-empty organisms list."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE p.experiment_count > 0
          AND (p.organisms IS NULL OR size(p.organisms) = 0)
        RETURN count(p) AS bad, collect(p.doi)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} publications with experiments have empty organisms: {result[0]['examples']}"
    )


# ── Experiment summary properties ──


def test_experiment_gene_count_not_null(run_query):
    """All Experiment nodes should have gene_count set (not null)."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.gene_count IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} experiments missing gene_count"
    )


def test_experiment_gene_count_accurate(run_query):
    """gene_count should match actual Changes_expression_of edge count."""
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->(g:Gene)
        WITH e, e.gene_count AS declared, count(r) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"gene_count mismatch for: {result[0]['examples']}"
    )


def test_experiment_distinct_gene_count_not_null(run_query):
    """All Experiment nodes should have distinct_gene_count set (not null)."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.distinct_gene_count IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} experiments missing distinct_gene_count"
    )


def test_experiment_distinct_gene_count_accurate(run_query):
    """distinct_gene_count should equal count(DISTINCT gene) over Changes_expression_of edges."""
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->(g:Gene)
        WITH e, e.distinct_gene_count AS declared, count(DISTINCT g) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"distinct_gene_count mismatch for: {result[0]['examples']}"
    )


def test_experiment_distinct_gene_count_le_gene_count(run_query):
    """Invariant: distinct_gene_count <= gene_count always (gene_count is cumulative across TPs)."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.distinct_gene_count > e.gene_count
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiments violate distinct_gene_count <= gene_count: {result[0]['examples']}"
    )


def test_experiment_distinct_gene_count_equals_gene_count_for_single_tp(run_query):
    """For non-time-course or single-TP experiments, distinct_gene_count == gene_count.

    With a single timepoint, every edge points to a distinct gene (the adapter
    emits one row per gene per timepoint), so cumulative == distinct.
    """
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.time_point_count <= 1
          AND e.distinct_gene_count <> e.gene_count
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} single-TP experiments have distinct_gene_count != gene_count: {result[0]['examples']}"
    )


def test_experiment_significant_up_count_accurate(run_query):
    """significant_up_count should match live aggregation."""
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->(g:Gene)
        WHERE r.expression_status = 'significant_up'
        WITH e, e.significant_up_count AS declared, count(r) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"significant_up_count mismatch for: {result[0]['examples']}"
    )


def test_experiment_significant_down_count_accurate(run_query):
    """significant_down_count should match live aggregation."""
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->(g:Gene)
        WHERE r.expression_status = 'significant_down'
        WITH e, e.significant_down_count AS declared, count(r) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"significant_down_count mismatch for: {result[0]['examples']}"
    )


def test_experiment_time_point_count_accurate(run_query):
    """time_point_count should match distinct time points from edges."""
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->(g:Gene)
        WITH e, count(DISTINCT r.time_point_order) AS actual
        WHERE e.time_point_count <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"time_point_count mismatch for: {result[0]['examples']}"
    )


def test_experiment_no_empty_string_coculture_partner(run_query):
    """No Experiment should have coculture_partner = '' (empty string)."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.coculture_partner = ''
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiments have empty-string coculture_partner"
    )


def test_experiment_coculture_partner_null_when_no_partner(run_query):
    """Experiments without a treatment organism should have null coculture_partner.

    Both coculture and viral experiments have a treatment organism (partner).
    All other treatment types should have null coculture_partner.
    """
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT 'coculture' IN e.treatment_type
          AND NOT 'viral' IN e.treatment_type
          AND e.coculture_partner IS NOT NULL
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} non-partner experiments have coculture_partner set: {result[0]['examples']}"
    )


def test_experiment_parallel_arrays_aligned(run_query):
    """All time_point arrays should have same length = time_point_count."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE e.time_point_count <> size(e.time_point_labels)
           OR e.time_point_count <> size(e.time_point_orders)
           OR e.time_point_count <> size(e.time_point_hours)
           OR e.time_point_count <> size(e.time_point_totals)
           OR e.time_point_count <> size(e.time_point_significant_up)
           OR e.time_point_count <> size(e.time_point_significant_down)
           OR e.time_point_count <> size(e.time_point_growth_phases)
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiments have misaligned time_point arrays: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Expression edge computed properties
# ---------------------------------------------------------------------------

def test_expression_status_populated(run_query):
    """Every Changes_expression_of edge must have expression_status set."""
    result = run_query("""
        MATCH ()-[r:Changes_expression_of]->()
        WHERE r.expression_status IS NULL
        RETURN count(r) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} Changes_expression_of edges missing expression_status"
    )


def test_expression_status_values_valid(run_query):
    """expression_status values must be in the known set."""
    result = run_query("""
        MATCH ()-[r:Changes_expression_of]->()
        WITH DISTINCT r.expression_status AS status
        WHERE NOT status IN ['significant_up', 'significant_down', 'not_significant']
        RETURN collect(status) AS bad
    """)
    assert result[0]["bad"] == [], (
        f"Unexpected expression_status values: {result[0]['bad']}"
    )


def test_rank_by_effect_contiguous(run_query):
    """rank_by_effect should be 1..N contiguous per (experiment, time_point_order)."""
    result = run_query("""
        MATCH (e:Experiment)-[r:Changes_expression_of]->()
        WITH e.id AS eid, r.time_point_order AS tp, collect(r.rank_by_effect) AS ranks
        WITH eid, tp, ranks,
             apoc.coll.sort(ranks) AS sorted,
             size(ranks) AS n
        WHERE sorted <> range(1, n)
        RETURN count(*) AS bad, collect([eid, tp])[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} (experiment, timepoint) groups have non-contiguous rank_by_effect: "
        f"{result[0]['examples']}"
    )


def test_rank_up_only_on_significant_up(run_query):
    """rank_up must be non-null iff expression_status = 'significant_up'."""
    result = run_query("""
        MATCH ()-[r:Changes_expression_of]->()
        WITH
          count(CASE WHEN r.rank_up IS NOT NULL AND r.expression_status <> 'significant_up'
                     THEN 1 END) AS rank_up_on_wrong_status,
          count(CASE WHEN r.rank_up IS NULL AND r.expression_status = 'significant_up'
                     THEN 1 END) AS sig_up_missing_rank
        RETURN rank_up_on_wrong_status, sig_up_missing_rank
    """)
    row = result[0]
    assert row["rank_up_on_wrong_status"] == 0, (
        f"{row['rank_up_on_wrong_status']} edges have rank_up set but status != significant_up"
    )
    assert row["sig_up_missing_rank"] == 0, (
        f"{row['sig_up_missing_rank']} significant_up edges are missing rank_up"
    )


def test_rank_down_only_on_significant_down(run_query):
    """rank_down must be non-null iff expression_status = 'significant_down'."""
    result = run_query("""
        MATCH ()-[r:Changes_expression_of]->()
        WITH
          count(CASE WHEN r.rank_down IS NOT NULL AND r.expression_status <> 'significant_down'
                     THEN 1 END) AS rank_down_on_wrong_status,
          count(CASE WHEN r.rank_down IS NULL AND r.expression_status = 'significant_down'
                     THEN 1 END) AS sig_down_missing_rank
        RETURN rank_down_on_wrong_status, sig_down_missing_rank
    """)
    row = result[0]
    assert row["rank_down_on_wrong_status"] == 0, (
        f"{row['rank_down_on_wrong_status']} edges have rank_down set but status != significant_down"
    )
    assert row["sig_down_missing_rank"] == 0, (
        f"{row['sig_down_missing_rank']} significant_down edges are missing rank_down"
    )


# ---------------------------------------------------------------------------
# Gene cluster routing signals
# ---------------------------------------------------------------------------

def test_cluster_membership_count_populated(run_query):
    """All Gene nodes must have non-null cluster_membership_count."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.cluster_membership_count IS NULL
        RETURN count(g) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} genes missing cluster_membership_count"
    )


def test_cluster_types_is_list(run_query):
    """All Gene nodes must have cluster_types as a list (possibly empty)."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.cluster_types IS NULL
        RETURN count(g) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} genes missing cluster_types"
    )


# ---------------------------------------------------------------------------
# BriteCategory computed properties
# ---------------------------------------------------------------------------

def test_britecategory_computed_populated(run_query):
    """All BriteCategory nodes must have non-null member_ko_count, gene_count, organism_count."""
    result = run_query("""
        MATCH (b:BriteCategory)
        WHERE b.member_ko_count IS NULL
           OR b.gene_count IS NULL
           OR b.organism_count IS NULL
        RETURN count(b) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} BriteCategory nodes missing a computed property"
    )


# ---------------------------------------------------------------------------
# Index presence — hardcoded list of what post-import.sh should create.
# If a CREATE INDEX line is dropped during refactoring, this test fails.
# ---------------------------------------------------------------------------

EXPECTED_INDEXES = {
    # scalar
    "gene_locus_tag_idx",
    "gene_name_idx",
    "gene_organism_name_idx",
    "ortholog_group_id_idx",
    "ortholog_group_name_idx",
    "ortholog_group_level_idx",
    "ortholog_group_rank_idx",
    "pfam_name_idx",
    "pfam_clan_name_idx",
    "brite_category_tree_idx",
    "brite_category_level_idx",
    "brite_category_name_idx",
    # TCDB / CAZy scalar
    "tcdb_family_level_idx",
    "tcdb_family_level_kind_idx",
    "tcdb_family_tcdb_id_idx",
    "tcdb_family_tc_class_id_idx",
    "cazy_family_level_idx",
    "cazy_family_level_kind_idx",
    "cazy_family_cazy_id_idx",
    "experiment_id_idx",
    "experiment_organism_idx",
    "experiment_treatment_type_idx",
    "experiment_background_factors_idx",
    "experiment_omics_type_idx",
    "organism_type_idx",
    "clustering_analysis_organism_idx",
    "clustering_analysis_method_idx",
    "clustering_analysis_type_idx",
    # DerivedMetric scalar (Plan 3)
    "derived_metric_metric_type_idx",
    "derived_metric_value_kind_idx",
    "derived_metric_compartment_idx",
    "derived_metric_omics_type_idx",
    "derived_metric_treatment_type_idx",
    "derived_metric_organism_idx",
    "derived_metric_experiment_idx",
    # Experiment.compartment scalar (Plan 2 adapter-emitted; Plan 3 indexed)
    "experiment_compartment_idx",
    # full-text
    "geneFullText",
    "biologicalProcessFullText",
    "molecularFunctionFullText",
    "cellularComponentFullText",
    "ecNumberFullText",
    "keggFullText",
    "cogCategoryFullText",
    "cyanorakRoleFullText",
    "tigrRoleFullText",
    "orthologGroupFullText",
    "pfamFullText",
    "pfamClanFullText",
    "briteCategoryFullText",
    "tcdbFamilyFullText",
    "cazyFamilyFullText",
    "publicationFullText",
    "experimentFullText",
    "clusteringAnalysisFullText",
    "geneClusterFullText",
    # DerivedMetric full-text (Plan 3)
    "derivedMetricFullText",
}


def test_expected_indexes_present_and_online(run_query):
    """All indexes listed in post-import.sh must exist and be ONLINE."""
    result = run_query(
        "SHOW INDEXES YIELD name, type, state WHERE type <> 'LOOKUP' "
        "RETURN name, state"
    )
    actual = {row["name"]: row["state"] for row in result}
    missing = EXPECTED_INDEXES - set(actual.keys())
    not_online = {
        name: state for name, state in actual.items()
        if name in EXPECTED_INDEXES and state != "ONLINE"
    }
    assert not missing, f"Missing indexes: {sorted(missing)}"
    assert not not_online, f"Indexes not ONLINE: {not_online}"


# ---------------------------------------------------------------------------
# TCDB / CAZy / Metabolite UNION rollups (Commit 7)
# ---------------------------------------------------------------------------


@pytest.mark.kg
def test_tcdb_family_has_gene_count_and_organism_count(run_query):
    """Every TcdbFamily has gene_count + organism_count populated (default 0)."""
    rows = run_query(
        "MATCH (t:TcdbFamily) WHERE t.gene_count IS NULL OR t.organism_count IS NULL "
        "RETURN count(t) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_tcdb_family_member_count_populated_on_classes(run_query):
    rows = run_query(
        "MATCH (t:TcdbFamily {level_kind: 'tc_class'}) RETURN t.member_count AS mc"
    )
    assert all(r["mc"] is not None and r["mc"] >= 0 for r in rows)


@pytest.mark.kg
def test_tcdb_family_metabolite_count_subtree_aggregates(run_query):
    """At least one tc_class node has metabolite_count > 0 (some descendant has substrates)."""
    rows = run_query(
        "MATCH (cls:TcdbFamily {level_kind: 'tc_class'}) WHERE cls.metabolite_count > 0 "
        "RETURN count(cls) AS n"
    )
    assert rows[0]["n"] > 0


@pytest.mark.kg
def test_tcdb_family_tc_class_id_self_on_class_nodes(run_query):
    """Every tc_class node has tc_class_id == its own id."""
    rows = run_query(
        "MATCH (t:TcdbFamily {level_kind: 'tc_class'}) "
        "WHERE t.tc_class_id <> t.id RETURN count(t) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_tcdb_family_tc_class_id_points_to_class_for_descendants(run_query):
    """Non-class TcdbFamily nodes have tc_class_id pointing at a real tc_class."""
    rows = run_query("""
        MATCH (t:TcdbFamily) WHERE t.level_kind <> 'tc_class'
        OPTIONAL MATCH (cls:TcdbFamily {id: t.tc_class_id, level_kind: 'tc_class'})
        WITH count(t) AS total, count(cls) AS resolved
        RETURN total, resolved
    """)
    assert rows[0]["total"] == rows[0]["resolved"]


@pytest.mark.kg
def test_cazy_family_has_gene_count_and_organism_count(run_query):
    rows = run_query(
        "MATCH (c:CazyFamily) WHERE c.gene_count IS NULL OR c.organism_count IS NULL "
        "RETURN count(c) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_gene_tcdb_family_count_populated(run_query):
    """Every Gene has tcdb_family_count and cazy_family_count populated (default 0)."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.tcdb_family_count IS NULL OR g.cazy_family_count IS NULL "
        "RETURN count(g) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_gene_annotation_types_includes_tcdb_when_edges_present(run_query):
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_tcdb_family]->() RETURN g.annotation_types AS at LIMIT 5"
    )
    assert all("tcdb" in r["at"] for r in rows)


@pytest.mark.kg
def test_gene_annotation_types_includes_cazy_when_edges_present(run_query):
    rows = run_query(
        "MATCH (g:Gene)-[:Gene_has_cazy_family]->() RETURN g.annotation_types AS at LIMIT 5"
    )
    assert all("cazy" in r["at"] for r in rows)


@pytest.mark.kg
def test_gene_metabolite_count_populated_for_genes_with_chemistry(run_query):
    """Every Gene whose chemistry edges actually reach a Metabolite has
    metabolite_count > 0.

    Looser than "has any chemistry edge": some genes annotate to KEGG glycan-only
    reactions (Reaction.compounds = []) or to TCDB sub-specifications without
    substrate annotations (e.g. 2.A.7.11.2 has no substrates in TCDB), so their
    metabolite_count is legitimately 0. We only assert population for the genes
    whose UNION-2hop actually yields a metabolite.
    """
    rows = run_query("""
        MATCH (g:Gene)
        WHERE EXISTS { (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(:Metabolite) }
           OR EXISTS {
                (g)-[:Gene_has_tcdb_family]->(:TcdbFamily)
                  <-[:Tcdb_family_is_a_tcdb_family*0..]-(:TcdbFamily {level_kind: 'tc_specificity'})
                  -[:Tcdb_family_transports_metabolite]->(:Metabolite)
              }
        RETURN count(CASE WHEN g.metabolite_count IS NULL OR g.metabolite_count = 0 THEN 1 END) AS n
    """)
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_metabolite_transporter_count_populated(run_query):
    rows = run_query(
        "MATCH (m:Metabolite) WHERE m.transporter_count IS NULL RETURN count(m) AS n"
    )
    assert rows[0]["n"] == 0


@pytest.mark.kg
def test_organism_has_metabolite_includes_transport_path(run_query):
    """At least some Organism_has_metabolite edges exist that come ONLY from transport
    (organism has no Reaction touching that metabolite, but does have a TCDB
    family with substrate edge to it)."""
    rows = run_query("""
        MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
        WHERE 'transport' IN m.evidence_sources
          AND NOT 'metabolism' IN m.evidence_sources
        RETURN count(*) AS n
    """)
    assert rows[0]["n"] > 0


@pytest.mark.kg
def test_metabolomics_evidence_source_appears_when_papers_integrated(run_query):
    """After Capovilla/Kujawinski integration, at least some Metabolite nodes
    must carry 'metabolomics' in evidence_sources. Skipped pre-integration."""
    rows = run_query(
        "MATCH (m:Metabolite) WHERE 'metabolomics' IN m.evidence_sources RETURN count(m) AS n"
    )
    n = rows[0]["n"]
    if n == 0:
        pytest.skip("No metabolomics-source metabolites yet (Phase 2 papers not integrated)")
    assert n > 0


# ---------------------------------------------------------------------------
# gene_neighbors composite index (docs/kg-specs/kg-spec-gene-neighbors.md)
#
# A composite RANGE index Gene(organism_name, contig, start) backs the
# gene_neighbors genomic-window query: equality prefix (organism_name, contig)
# + ordered range suffix start. Without it the planner scans an organism's full
# gene set and Top-sorts in memory; with it the seek folds in the contig
# equality and the start range, touching only the bounded window.
# ---------------------------------------------------------------------------

GENE_NEIGHBORS_INDEX = "gene_org_contig_start_idx"


@pytest.mark.kg
def test_gene_neighbors_composite_index_online(run_query):
    """The composite RANGE index backing gene_neighbors exists, is ONLINE, and
    covers (organism_name, contig, start) in that order."""
    rows = run_query(
        "SHOW INDEXES YIELD name, type, state, properties "
        f"WHERE name = '{GENE_NEIGHBORS_INDEX}' RETURN type, state, properties"
    )
    assert len(rows) == 1, f"{GENE_NEIGHBORS_INDEX} not found"
    assert rows[0]["state"] == "ONLINE", f"{GENE_NEIGHBORS_INDEX} state = {rows[0]['state']}"
    assert rows[0]["type"] == "RANGE", f"{GENE_NEIGHBORS_INDEX} type = {rows[0]['type']}"
    assert rows[0]["properties"] == ["organism_name", "contig", "start"], (
        f"{GENE_NEIGHBORS_INDEX} properties = {rows[0]['properties']}"
    )


def _plan_details(plan):
    """Concatenate every operator's `Details` string from a Neo4j plan tree.

    The query plan is a nested dict ({operatorType, args, children, ...}); the
    index a seek uses is described by its property signature in args['Details']
    (e.g. "RANGE INDEX x:Gene(organism_name, contig, start) WHERE ...") rather
    than by index name, so we match on the signature.
    """
    parts = [str((plan.get("args") or {}).get("Details", ""))]
    for child in plan.get("children", []) or []:
        parts.append(_plan_details(child))
    return " ".join(parts)


@pytest.mark.kg
def test_gene_neighbors_query_uses_composite_index(neo4j_driver, run_query):
    """The bounded-window neighbor subquery must plan against the composite index
    (contig equality + start range folded into the seek), not a full-organism
    scan + Filter. Discriminator: the "before" plan referenced only
    Gene(organism_name)."""
    anchor = run_query(
        "MATCH (g:Gene) WHERE g.contig IS NOT NULL AND g.start IS NOT NULL "
        "RETURN g.locus_tag AS lt LIMIT 1"
    )
    assert anchor, "no positioned Gene found to anchor the plan check"
    lt = anchor[0]["lt"]
    cypher = (
        "EXPLAIN "
        "MATCH (a:Gene {locus_tag: $lt}) "
        "MATCH (d:Gene) WHERE d.organism_name = a.organism_name "
        "AND d.contig = a.contig AND d.start > a.start "
        "WITH d ORDER BY d.start ASC LIMIT 5 RETURN collect(d) AS ds"
    )
    with neo4j_driver.session() as session:
        plan = session.run(cypher, lt=lt).consume().plan
    details = _plan_details(plan)
    assert "organism_name, contig, start" in details, (
        "neighbor query plan does not seek the composite index "
        f"Gene(organism_name, contig, start); plan details:\n{details}"
    )
