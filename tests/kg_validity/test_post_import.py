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
                        'kegg', 'ec', 'cyanorak_role', 'tigr_role']
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
        WHERE NOT e.treatment_type IN ['coculture', 'viral']
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
           OR e.time_point_count <> size(e.time_point_totals)
           OR e.time_point_count <> size(e.time_point_significant_up)
           OR e.time_point_count <> size(e.time_point_significant_down)
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiments have misaligned time_point arrays: {result[0]['examples']}"
    )
