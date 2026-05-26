"""KG validity: DataSource nodes (F2)."""

import pytest

pytestmark = pytest.mark.kg


def test_five_data_source_nodes(run_query):
    """Deployment has exactly 5 DataSource nodes (incl. psortb, added Phase-2).

    BioCypher prepends the `data_source:` CURIE prefix per schema_config.
    """
    result = run_query("MATCH (d:DataSource) RETURN d.id AS id ORDER BY d.id")
    ids = [row["id"] for row in result]
    assert ids == [
        "data_source:cyanorak",
        "data_source:eggnog",
        "data_source:ncbi",
        "data_source:psortb",
        "data_source:uniprot",
    ]


def test_info_types_non_empty(run_query):
    """Every DataSource must have at least one info_type listed."""
    result = run_query("""
        MATCH (d:DataSource)
        WHERE size(d.info_types) = 0
        RETURN d.id AS id
    """)
    assert result == [], f"DataSources with empty info_types: {result}"


def test_eggnog_provenance_is_tool_run(run_query):
    result = run_query("""
        MATCH (d:DataSource {id: 'data_source:eggnog'}) RETURN d.provenance AS p
    """)
    assert result and result[0]["p"] == "tool_run"


def test_psortb_provenance_is_tool_run(run_query):
    result = run_query("""
        MATCH (d:DataSource {id: 'data_source:psortb'}) RETURN d.provenance AS p
    """)
    assert result and result[0]["p"] == "tool_run"


def test_cyanorak_organism_restricted(run_query):
    result = run_query("""
        MATCH (d:DataSource {id: 'data_source:cyanorak'})
        RETURN d.scope AS scope, size(d.applies_to_organisms) AS n
    """)
    assert result and result[0]["scope"] == "organism_restricted"
    assert result[0]["n"] > 0


def test_every_gene_has_ncbi_in_contributing_sources(run_query):
    """The ncbi source is universal — every Gene must list it."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE NOT 'ncbi' IN g.contributing_sources
        RETURN count(*) AS missing
    """)
    assert result[0]["missing"] == 0
