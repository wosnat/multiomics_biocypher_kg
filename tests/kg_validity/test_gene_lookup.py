"""
Gene lookup infrastructure tests for the multi-omics knowledge graph.

Validates the schema/index changes needed for MCP get_gene and find_gene services:
- Scalar indexes on locus_tag, gene_name, organism_strain
- Full-text index geneFullText
- New computed fields: organism_strain, gene_summary, all_identifiers
- go_term_descriptions stored as list (str[])
- get_gene and find_gene query patterns work correctly
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Index existence
# ---------------------------------------------------------------------------

def test_scalar_indexes_exist(run_query):
    """Scalar indexes on locus_tag, gene_name, organism_strain must exist."""
    rows = run_query("SHOW INDEXES YIELD name, labelsOrTypes, properties RETURN name, labelsOrTypes, properties")
    index_keys = {(tuple(r["labelsOrTypes"]), tuple(r["properties"])) for r in rows if r["labelsOrTypes"] is not None}
    expected = [
        (("Gene",), ("locus_tag",)),
        (("Gene",), ("gene_name",)),
        (("Gene",), ("organism_strain",)),
    ]
    for label_tuple, prop_tuple in expected:
        assert (label_tuple, prop_tuple) in index_keys, (
            f"Missing scalar index on Gene({prop_tuple[0]})"
        )


def test_fulltext_index_exists(run_query):
    """Full-text index geneFullText must exist with expected properties."""
    rows = run_query(
        "SHOW FULLTEXT INDEXES YIELD name, labelsOrTypes, properties "
        "WHERE name = 'geneFullText' "
        "RETURN name, labelsOrTypes, properties"
    )
    assert len(rows) == 1, "geneFullText index not found"
    props = set(rows[0]["properties"])
    expected_props = {
        "gene_summary", "gene_synonyms", "product_cyanorak",
        "alternate_functional_descriptions", "go_term_descriptions",
        "pfam_names", "pfam_descriptions", "eggnog_og_descriptions",
    }
    assert expected_props == props, f"Index properties mismatch: missing={expected_props - props}, extra={props - expected_props}"


# ---------------------------------------------------------------------------
# New field population
# ---------------------------------------------------------------------------

def test_organism_strain_populated_on_all_genes(run_query):
    """Every Gene node must have a non-null organism_strain."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.organism_strain IS NULL RETURN count(g) AS n"
    )
    assert rows[0]["n"] == 0, "Some Gene nodes are missing organism_strain"


def test_organism_strain_values(run_query):
    """organism_strain values should include known strains."""
    rows = run_query(
        "MATCH (g:Gene) RETURN DISTINCT g.organism_strain AS os ORDER BY os"
    )
    values = {r["os"] for r in rows}
    # Spot-check a few expected values
    for expected in ["Prochlorococcus MED4", "Synechococcus CC9311"]:
        assert expected in values, f"Expected organism_strain '{expected}' not found"


def test_gene_summary_populated(run_query):
    """gene_summary should be populated on well-annotated genes (quality >= 2)."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.gene_summary IS NOT NULL "
        "RETURN count(g) AS n"
    )
    # All genes should have a summary (even if just the locus_tag)
    total = run_query("MATCH (g:Gene) RETURN count(g) AS n")[0]["n"]
    assert rows[0]["n"] > total * 0.95, (
        f"Only {rows[0]['n']}/{total} genes have gene_summary — expected >95%"
    )


def test_gene_summary_format(run_query):
    """gene_summary for a known gene should contain gene_name and product."""
    rows = run_query(
        "MATCH (g:Gene {locus_tag: 'PMM0001'}) RETURN g.gene_summary AS s"
    )
    assert len(rows) == 1
    summary = rows[0]["s"]
    assert "dnaN" in summary, "gene_summary should contain gene_name"
    assert "::" in summary, "gene_summary should use :: separator"


def test_all_identifiers_populated(run_query):
    """all_identifiers should be a non-empty array on most genes."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.all_identifiers IS NOT NULL "
        "AND size(g.all_identifiers) > 0 "
        "RETURN count(g) AS n"
    )
    total = run_query("MATCH (g:Gene) RETURN count(g) AS n")[0]["n"]
    assert rows[0]["n"] > total * 0.90, (
        f"Only {rows[0]['n']}/{total} genes have all_identifiers — expected >90%"
    )


def test_all_identifiers_contains_known_ids(run_query):
    """PMM0001 all_identifiers should include its RefSeq protein ID and Cyanorak locus tag."""
    rows = run_query(
        "MATCH (g:Gene {locus_tag: 'PMM0001'}) "
        "RETURN g.all_identifiers AS ids"
    )
    assert len(rows) == 1
    ids = rows[0]["ids"]
    assert isinstance(ids, list), "all_identifiers should be a list"
    # Should contain at least the Cyanorak locus tag
    assert any("CK_" in i for i in ids), f"Expected Cyanorak ID in all_identifiers: {ids}"


def test_go_term_descriptions_is_list(run_query):
    """go_term_descriptions should be stored as a list (str[]), not a string."""
    rows = run_query(
        "MATCH (g:Gene {locus_tag: 'PMM0001'}) "
        "RETURN g.go_term_descriptions AS gtd"
    )
    assert len(rows) == 1
    gtd = rows[0]["gtd"]
    assert gtd is not None, "PMM0001 should have go_term_descriptions"
    assert isinstance(gtd, list), f"go_term_descriptions should be list, got {type(gtd).__name__}"


# ---------------------------------------------------------------------------
# get_gene query pattern
# ---------------------------------------------------------------------------

def test_get_gene_by_locus_tag(run_query):
    """get_gene: exact match on locus_tag (scalar index)."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.locus_tag = $id "
        "RETURN g.locus_tag AS lt, g.organism_strain AS os",
        id="PMM0001",
    )
    assert len(rows) == 1
    assert rows[0]["lt"] == "PMM0001"
    assert rows[0]["os"] is not None


def test_get_gene_by_gene_name(run_query):
    """get_gene: exact match on gene_name (scalar index)."""
    rows = run_query(
        "MATCH (g:Gene) WHERE g.gene_name = $id "
        "RETURN g.locus_tag AS lt, g.organism_strain AS os",
        id="dnaN",
    )
    # dnaN exists in multiple strains
    assert len(rows) >= 1
    locus_tags = {r["lt"] for r in rows}
    assert "PMM0001" in locus_tags


def test_get_gene_by_all_identifiers(run_query):
    """get_gene: match via all_identifiers array scan."""
    # Look up a Cyanorak ID that should be in all_identifiers
    rows = run_query(
        "MATCH (g:Gene) WHERE $id IN g.all_identifiers "
        "RETURN g.locus_tag AS lt",
        id="CK_Pro_MED4_00001",
    )
    assert len(rows) == 1
    assert rows[0]["lt"] == "PMM0001"


def test_get_gene_organism_filter(run_query):
    """get_gene: organism filter narrows results to one strain."""
    rows = run_query(
        "MATCH (g:Gene) "
        "WHERE g.gene_name = $id AND g.organism_strain = $organism "
        "RETURN g.locus_tag AS lt",
        id="dnaN",
        organism="Prochlorococcus MED4",
    )
    assert len(rows) == 1
    assert rows[0]["lt"] == "PMM0001"


def test_get_gene_full_query(run_query):
    """get_gene: combined query pattern from the plan."""
    rows = run_query(
        "MATCH (g:Gene) "
        "WHERE (g.locus_tag = $id OR g.gene_name = $id OR $id IN g.all_identifiers) "
        "  AND ($organism IS NULL OR g.organism_strain = $organism) "
        "RETURN g.locus_tag AS lt, g.organism_strain AS os "
        "LIMIT 5",
        id="PMM0001",
        organism=None,
    )
    assert len(rows) == 1
    assert rows[0]["lt"] == "PMM0001"


# ---------------------------------------------------------------------------
# find_gene query pattern
# ---------------------------------------------------------------------------

def test_find_gene_known_query(run_query):
    """find_gene: full-text search for 'photosystem' should return results."""
    rows = run_query(
        "CALL db.index.fulltext.queryNodes('geneFullText', $search_text) "
        "YIELD node AS g, score "
        "RETURN g.locus_tag AS lt, score "
        "ORDER BY score DESC LIMIT 10",
        search_text="photosystem",
    )
    assert len(rows) > 0, "Full-text search for 'photosystem' returned no results"


def test_find_gene_organism_filter(run_query):
    """find_gene: organism filter narrows full-text results."""
    rows = run_query(
        "CALL db.index.fulltext.queryNodes('geneFullText', $search_text) "
        "YIELD node AS g, score "
        "WHERE g.organism_strain = $organism "
        "RETURN g.locus_tag AS lt, g.organism_strain AS os "
        "ORDER BY score DESC LIMIT 10",
        search_text="dnaN",
        organism="Prochlorococcus MED4",
    )
    assert len(rows) >= 1
    # All results should be MED4
    for r in rows:
        assert r["os"] == "Prochlorococcus MED4"


def test_find_gene_quality_filter(run_query):
    """find_gene: min_quality filter excludes low-quality genes."""
    rows = run_query(
        "CALL db.index.fulltext.queryNodes('geneFullText', $search_text) "
        "YIELD node AS g, score "
        "WHERE g.annotation_quality >= $min_quality "
        "RETURN g.locus_tag AS lt, g.annotation_quality AS q "
        "ORDER BY score DESC LIMIT 10",
        search_text="hypothetical",
        min_quality=2,
    )
    for r in rows:
        assert r["q"] >= 2


def test_find_gene_scores_descending(run_query):
    """find_gene: results should be ordered by score descending."""
    rows = run_query(
        "CALL db.index.fulltext.queryNodes('geneFullText', $search_text) "
        "YIELD node AS g, score "
        "RETURN score "
        "ORDER BY score DESC LIMIT 10",
        search_text="iron transport",
    )
    scores = [r["score"] for r in rows]
    assert scores == sorted(scores, reverse=True), "Scores should be in descending order"
