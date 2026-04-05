"""
OrganismTaxon validation tests.

Validates:
- No garbage taxonomy nodes (wrong NCBI taxon mapping)
- Species populated on Alteromonas strains
- All organisms have bacterial lineage (except Phage)
- Correct organism count
- Precomputed stats (gene_count, publication_count, etc.) are consistent
- Publication.organisms values align with OrganismTaxon.preferred_name
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Garbage node detection
# ---------------------------------------------------------------------------

def test_no_garbage_taxonomy_nodes(run_query):
    """No OrganismTaxon nodes should map to non-marine-bacteria taxa."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.preferred_name IN ['Pseudohoeflea', 'Marinobacter', 'Thalassospira']
        RETURN count(o) AS cnt
    """)
    assert result[0]["cnt"] == 0, (
        "Garbage taxonomy nodes still present (wrong NCBI taxon IDs)"
    )


def test_no_non_bacterial_organisms(run_query):
    """All organisms with superkingdom must be Bacteria (except Phage)."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.superkingdom IS NOT NULL
          AND o.superkingdom <> 'Bacteria'
          AND o.preferred_name <> 'Phage'
        RETURN o.preferred_name AS name, o.superkingdom AS sk, o.ncbi_taxon_id AS taxid
    """)
    assert len(result) == 0, (
        f"Non-bacterial organisms found: {result}"
    )


def test_organism_count(run_query):
    """28 OrganismTaxon nodes: 23 genome strains + 5 treatment organisms."""
    result = run_query("MATCH (o:OrganismTaxon) RETURN count(o) AS cnt")
    assert result[0]["cnt"] == 28, (
        f"Expected 28 OrganismTaxon nodes, got {result[0]['cnt']}"
    )


# ---------------------------------------------------------------------------
# Species on Alteromonas strains
# ---------------------------------------------------------------------------

def test_alteromonas_strains_have_species(run_query):
    """All Alteromonas strain-level nodes must have species = 'Alteromonas macleodii'."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.genus = 'Alteromonas' AND o.strain_name IS NOT NULL
        RETURN o.preferred_name AS name, o.species AS species
    """)
    assert len(result) == 3, f"Expected 3 Alteromonas strains, got {len(result)}"
    for r in result:
        assert r["species"] == "Alteromonas macleodii", (
            f"{r['name']} has species={r['species']!r}, expected 'Alteromonas macleodii'"
        )


# ---------------------------------------------------------------------------
# Precomputed stats consistency
# ---------------------------------------------------------------------------

def test_gene_count_matches_live(run_query):
    """Precomputed gene_count must match actual Gene_belongs_to_organism edge count."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o)
        WITH o, count(g) AS live_count
        WHERE o.gene_count <> live_count
        RETURN o.preferred_name AS name, o.gene_count AS stored, live_count
    """)
    assert len(result) == 0, (
        f"gene_count mismatch on: {result}"
    )


def test_publication_count_matches_live(run_query):
    """Precomputed publication_count must match actual Publication.organisms join."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (p:Publication)
          WHERE ANY(org IN p.organisms WHERE org = o.preferred_name)
        WITH o, count(DISTINCT p) AS live_count
        WHERE o.publication_count <> live_count
        RETURN o.preferred_name AS name, o.publication_count AS stored, live_count
    """)
    assert len(result) == 0, (
        f"publication_count mismatch on: {result}"
    )


def test_experiment_count_matches_live(run_query):
    """Precomputed experiment_count must match sum of publication experiment_counts."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (p:Publication)
          WHERE ANY(org IN p.organisms WHERE org = o.preferred_name)
        WITH o, CASE WHEN count(p) > 0 THEN sum(p.experiment_count) ELSE 0 END AS live_count
        WHERE o.experiment_count <> live_count
        RETURN o.preferred_name AS name, o.experiment_count AS stored, live_count
    """)
    assert len(result) == 0, (
        f"experiment_count mismatch on: {result}"
    )


def test_no_null_precomputed_lists(run_query):
    """treatment_types and omics_types must never be null (empty list OK)."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.treatment_types IS NULL OR o.omics_types IS NULL
        RETURN o.preferred_name AS name
    """)
    assert len(result) == 0, (
        f"Null precomputed lists on: {[r['name'] for r in result]}"
    )


def test_gene_count_not_null(run_query):
    """gene_count must be set on all OrganismTaxon nodes."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.gene_count IS NULL
        RETURN o.preferred_name AS name
    """)
    assert len(result) == 0, (
        f"gene_count is null on: {[r['name'] for r in result]}"
    )


# ---------------------------------------------------------------------------
# Publication.organisms alignment
# ---------------------------------------------------------------------------

def test_publication_organisms_match_organism_taxon(run_query):
    """Every value in Publication.organisms must match an OrganismTaxon.preferred_name."""
    result = run_query("""
        MATCH (p:Publication)
        UNWIND p.organisms AS pub_org
        WITH DISTINCT pub_org
        WHERE NOT EXISTS { MATCH (o:OrganismTaxon) WHERE o.preferred_name = pub_org }
        RETURN pub_org
    """)
    assert len(result) == 0, (
        f"Publication.organisms values with no matching OrganismTaxon: "
        f"{[r['pub_org'] for r in result]}"
    )


# ---------------------------------------------------------------------------
# Spot checks on specific organisms
# ---------------------------------------------------------------------------

def test_med4_has_most_publications(run_query):
    """MED4 is the most-studied strain; should have the most publications."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IS NOT NULL
        RETURN o.preferred_name AS name, o.publication_count AS pubs
        ORDER BY pubs DESC LIMIT 1
    """)
    assert "MED4" in result[0]["name"], (
        f"Expected MED4 to have most publications, got {result[0]['name']} "
        f"with {result[0]['pubs']}"
    )


def test_phage_has_zero_genes(run_query):
    """Phage is a treatment organism only — no genes loaded."""
    result = run_query("""
        MATCH (o:OrganismTaxon {preferred_name: 'Phage'})
        RETURN o.gene_count AS gc
    """)
    assert result[0]["gc"] == 0, (
        f"Phage should have 0 genes, got {result[0]['gc']}"
    )
