"""Live checks: every DataSource carries citations; Schema_info carries the
build framework's (config/data_source_citations.yaml)."""

import pytest

pytestmark = pytest.mark.kg

EXPECTED_DATA_SOURCES = 17


def test_schema_info_build_citations(run_query):
    rows = run_query("MATCH (s:Schema_info {id: 'schema_info'}) RETURN s.build_citations AS c")
    assert rows and rows[0]["c"], "Schema_info.build_citations missing or empty"
    assert any("BioCypher" in c for c in rows[0]["c"])


def test_every_data_source_has_citations(run_query):
    rows = run_query(
        "MATCH (d:DataSource) RETURN d.id AS id, d.citations AS c, d.dois AS dois, d.url AS url"
    )
    assert len(rows) == EXPECTED_DATA_SOURCES, sorted(r["id"] for r in rows)
    for r in rows:
        assert r["c"], f"{r['id']}: no citations"
        assert len(r["dois"] or []) == len(r["c"]), f"{r['id']}: dois not parallel to citations"
        assert r["url"], f"{r['id']}: no url"
