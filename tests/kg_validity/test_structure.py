"""
Structural integrity tests for the multi-omics knowledge graph.

Validates:
- All expected node types are present and populated
- Minimum counts (sanity thresholds, not exact — data evolves)
- Orphan detection: nodes that are missing required relationships
- Key properties are populated on critical node types
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Node type presence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", [
    "Gene",
    "Protein",
    "OrganismTaxon",
    "Publication",
    "EnvironmentalCondition",
    "Cyanorak_cluster",
])
def test_node_type_exists(run_query, label):
    """Every expected node type must have at least one node in the graph."""
    result = run_query(f"MATCH (n:{label}) RETURN count(n) AS cnt")
    assert result[0]["cnt"] > 0, f"No nodes found with label :{label}"


# ---------------------------------------------------------------------------
# Minimum node counts (thresholds, not exact)
# ---------------------------------------------------------------------------

def test_gene_count_minimum(run_query):
    """MED4 alone has ~1900 genes; 12 strains combined should exceed 5000."""
    result = run_query("MATCH (g:Gene) RETURN count(g) AS cnt")
    assert result[0]["cnt"] > 5000, (
        f"Only {result[0]['cnt']} Gene nodes found; expected > 5000 for all strains"
    )


def test_organism_count(run_query):
    """
    At least 12 OrganismTaxon nodes expected (7 Prochlorococcus + Synechococcus +
    Alteromonas strains + treatment organisms from papers).
    """
    result = run_query("MATCH (o:OrganismTaxon) RETURN count(o) AS cnt")
    assert result[0]["cnt"] >= 12, (
        f"Only {result[0]['cnt']} OrganismTaxon nodes; expected >= 12"
    )


def test_publication_count(run_query):
    """At least one publication must be present."""
    result = run_query("MATCH (p:Publication) RETURN count(p) AS cnt")
    assert result[0]["cnt"] >= 1, "No Publication nodes found"


def test_protein_count_minimum(run_query):
    """Proteins from UniProt for all strains should exceed 5000."""
    result = run_query("MATCH (p:Protein) RETURN count(p) AS cnt")
    assert result[0]["cnt"] > 5000, (
        f"Only {result[0]['cnt']} Protein nodes; expected > 5000"
    )


def test_cyanorak_cluster_count(run_query):
    """Cyanorak clusters group homologous genes; expect at least 1000 clusters."""
    result = run_query("MATCH (c:Cyanorak_cluster) RETURN count(c) AS cnt")
    assert result[0]["cnt"] > 1000, (
        f"Only {result[0]['cnt']} Cyanorak_cluster nodes; expected > 1000"
    )


# ---------------------------------------------------------------------------
# Orphan detection
# ---------------------------------------------------------------------------

def test_no_orphan_genes(run_query):
    """Every Gene must be linked to an OrganismTaxon via Gene_belongs_to_organism."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE NOT (g)-[:Gene_belongs_to_organism]->(:OrganismTaxon)
        RETURN count(g) AS orphans
    """)
    orphans = result[0]["orphans"]
    assert orphans == 0, (
        f"{orphans} Gene node(s) have no Gene_belongs_to_organism edge"
    )


def test_no_orphan_proteins(run_query):
    """Every Protein must be linked to an OrganismTaxon via Protein_belongs_to_organism."""
    result = run_query("""
        MATCH (p:Protein)
        WHERE NOT (p)-[:Protein_belongs_to_organism]->(:OrganismTaxon)
        RETURN count(p) AS orphans
    """)
    orphans = result[0]["orphans"]
    assert orphans == 0, (
        f"{orphans} Protein node(s) have no Protein_belongs_to_organism edge"
    )


def test_gene_encodes_protein_edges_exist(run_query):
    """Gene_encodes_protein edges must exist (created by UniProt adapter via RefSeq join)."""
    result = run_query(
        "MATCH ()-[r:Gene_encodes_protein]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 5000, (
        f"Only {cnt} Gene_encodes_protein edges found; expected > 5000 "
        f"for all strains"
    )


def test_gene_encodes_protein_links_correct_types(run_query):
    """Gene_encodes_protein must connect Protein -> Gene (source is protein, target is gene)."""
    result = run_query("""
        MATCH (src)-[r:Gene_encodes_protein]->(tgt)
        WHERE NOT src:Protein OR NOT tgt:Gene
        RETURN count(r) AS bad
    """)
    bad = result[0]["bad"]
    assert bad == 0, (
        f"{bad} Gene_encodes_protein edges connect wrong node types "
        f"(expected Protein -> Gene)"
    )


def test_no_orphan_proteins_without_gene(run_query):
    """
    Most proteins should be linked to a gene via Gene_encodes_protein.
    Allow up to 15% unlinked (some UniProt entries may lack RefSeq cross-refs).
    """
    result = run_query("""
        MATCH (p:Protein)
        OPTIONAL MATCH (p)-[:Gene_encodes_protein]->(g:Gene)
        WITH count(p) AS total, count(g) AS linked
        RETURN total, linked, total - linked AS unlinked
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Protein nodes found")
    unlinked_fraction = row["unlinked"] / row["total"]
    assert unlinked_fraction < 0.15, (
        f"{row['unlinked']} / {row['total']} proteins ({unlinked_fraction:.1%}) "
        f"have no Gene_encodes_protein edge; threshold is < 15%"
    )


def test_gene_encodes_protein_no_duplicates(run_query):
    """
    Each (Protein, Gene) pair should have at most one Gene_encodes_protein edge.
    Duplicates would indicate a bug in the RefSeq join logic.
    """
    result = run_query("""
        MATCH (p:Protein)-[r:Gene_encodes_protein]->(g:Gene)
        WITH p, g, count(r) AS edge_count
        WHERE edge_count > 1
        RETURN count(*) AS duplicates
    """)
    duplicates = result[0]["duplicates"]
    assert duplicates == 0, (
        f"{duplicates} (Protein, Gene) pair(s) have duplicate Gene_encodes_protein edges"
    )


def test_prochlorococcus_genes_in_cyanorak_clusters(run_query):
    """
    Prochlorococcus genes should mostly belong to a Cyanorak cluster.
    Allow up to 20% missing (some genes are absent from Cyanorak but present
    in NCBI; Alteromonas/Synechococcus genes are not covered by Cyanorak at all).
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.genus = 'Prochlorococcus'
        OPTIONAL MATCH (g)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster)
        WITH count(g) AS total, count(c) AS in_cluster
        RETURN total, in_cluster, total - in_cluster AS missing
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No Prochlorococcus Gene nodes found")
    missing_fraction = row["missing"] / row["total"]
    assert missing_fraction < 0.20, (
        f"{row['missing']} / {row['total']} Prochlorococcus genes ({missing_fraction:.1%}) "
        f"are not in any Cyanorak_cluster; threshold is < 20%"
    )


# ---------------------------------------------------------------------------
# Property presence on critical node types
# ---------------------------------------------------------------------------

def test_gene_locus_tag_present(run_query):
    """Every Gene must have a locus_tag (primary human-readable identifier)."""
    result = run_query("""
        MATCH (g:Gene) WHERE g.locus_tag IS NULL
        RETURN count(g) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} Gene node(s) are missing the locus_tag property"


def test_protein_name_present(run_query):
    """Every Protein must have a protein_name."""
    result = run_query("""
        MATCH (p:Protein) WHERE p.protein_name IS NULL
        RETURN count(p) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} Protein node(s) are missing the protein_name property"


def test_organism_name_present(run_query):
    """Every genomic OrganismTaxon (with a strain_name) must have organism_name."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IS NOT NULL AND o.organism_name IS NULL
        RETURN count(o) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} OrganismTaxon node(s) are missing organism_name"


def test_organism_taxon_id_present(run_query):
    """Every genomic OrganismTaxon (with a strain_name) must have ncbi_taxon_id."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IS NOT NULL AND o.ncbi_taxon_id IS NULL
        RETURN count(o) AS missing, collect(o.strain_name) AS strains
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} OrganismTaxon node(s) are missing ncbi_taxon_id: "
        f"{result[0]['strains']}"
    )


def test_publication_doi_or_pmid_present(run_query):
    """Every Publication should have at least a DOI or PMID."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE p.doi IS NULL AND p.pmid IS NULL
        RETURN count(p) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Publication node(s) have neither DOI nor PMID"
    )
