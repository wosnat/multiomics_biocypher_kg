"""KG validity tests for TCDB and CAZy ontologies."""
from __future__ import annotations

import pytest


@pytest.mark.kg
def test_tcdb_family_node_count_in_range(run_query):
    """Pruned subhierarchy walks above + below gene-annotated TCDB IDs.

    Local step 6 run produces ~12,881 kept IDs from 535 gene-annotated seeds
    after the hierarchy was widened to seed from acc2tcid.tsv (May 2026).
    The upper bound is intentionally generous to absorb growth as more strains
    or annotations land.
    """
    n = run_query("MATCH (t:TcdbFamily) RETURN count(t) AS n")[0]["n"]
    assert 100 <= n <= 25000, f"TcdbFamily count {n} outside expected 100-25000"


@pytest.mark.kg
def test_cazy_family_node_count_in_range(run_query):
    """Observed-only — small ontology."""
    n = run_query("MATCH (c:CazyFamily) RETURN count(c) AS n")[0]["n"]
    assert 5 <= n <= 100, f"CazyFamily count {n} outside expected 5-100"


@pytest.mark.kg
def test_gene_has_tcdb_family_edge_count(run_query):
    n = run_query("MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n")[0]["n"]
    assert 500 <= n <= 15000, f"Gene_has_tcdb_family count {n} outside 500-15000"


@pytest.mark.kg
def test_gene_has_cazy_family_edge_count(run_query):
    n = run_query("MATCH ()-[r:Gene_has_cazy_family]->() RETURN count(r) AS n")[0]["n"]
    assert 100 <= n <= 1500, f"Gene_has_cazy_family count {n} outside 100-1500"


@pytest.mark.kg
def test_tcdb_family_transports_metabolite_edge_count(run_query):
    """Substrate edges are rolled up from tc_specificity leaves to every ancestor
    so the upper bound covers ~5x the leaf-only count."""
    n = run_query(
        "MATCH ()-[r:Tcdb_family_transports_metabolite]->() RETURN count(r) AS n"
    )[0]["n"]
    assert 1000 <= n <= 100000, f"Tcdb_family_transports_metabolite count {n} outside 1000-100000"


@pytest.mark.kg
def test_tcdb_family_substrate_edges_present_at_every_level(run_query):
    """Substrate edges must exist at every TcdbFamily level (rolled up from
    tc_specificity leaves to ancestors). Levels with no descendant having a
    substrate may legitimately be empty, so this just asserts >0 levels overall."""
    rows = run_query("""
        MATCH (t:TcdbFamily)-[:Tcdb_family_transports_metabolite]->()
        RETURN DISTINCT t.level_kind AS lk
    """)
    levels = {r["lk"] for r in rows}
    assert "tc_specificity" in levels, "No leaf-level substrate edges (rollup base case missing)"
    # Expect at least one ancestor level present (tc_class, tc_subclass, tc_family,
    # or tc_subfamily) since the pruned hierarchy has multiple levels.
    ancestor_levels = levels - {"tc_specificity"}
    assert ancestor_levels, "No ancestor-level substrate edges — rollup not applied?"


@pytest.mark.kg
def test_tcdb_family_levels_are_valid(run_query):
    rows = run_query(
        "MATCH (t:TcdbFamily) RETURN DISTINCT t.level_kind AS lk, t.level AS lv"
    )
    expected = {
        ("tc_class", 0), ("tc_subclass", 1), ("tc_family", 2),
        ("tc_subfamily", 3), ("tc_specificity", 4),
    }
    actual = {(r["lk"], r["lv"]) for r in rows}
    assert actual.issubset(expected), f"Unexpected level_kind/level pairs: {actual - expected}"


@pytest.mark.kg
def test_cazy_family_levels_are_valid(run_query):
    rows = run_query(
        "MATCH (c:CazyFamily) RETURN DISTINCT c.level_kind AS lk, c.level AS lv"
    )
    expected = {("cazy_class", 0), ("cazy_family", 1), ("cazy_subfamily", 2)}
    actual = {(r["lk"], r["lv"]) for r in rows}
    assert actual.issubset(expected), f"Unexpected level_kind/level pairs: {actual - expected}"


@pytest.mark.kg
def test_tcdb_family_no_orphans_below_class(run_query):
    """Every non-class TcdbFamily must have a parent edge."""
    n_bad = run_query("""
        MATCH (t:TcdbFamily) WHERE t.level_kind <> 'tc_class'
          AND NOT EXISTS { (t)-[:Tcdb_family_is_a_tcdb_family]->() }
        RETURN count(t) AS n
    """)[0]["n"]
    assert n_bad == 0, f"{n_bad} non-class TcdbFamily nodes missing parent edge"


@pytest.mark.kg
def test_cazy_family_no_orphans_below_class(run_query):
    """Every non-class CazyFamily must have a parent edge."""
    n_bad = run_query("""
        MATCH (c:CazyFamily) WHERE c.level_kind <> 'cazy_class'
          AND NOT EXISTS { (c)-[:Cazy_family_is_a_cazy_family]->() }
        RETURN count(c) AS n
    """)[0]["n"]
    assert n_bad == 0, f"{n_bad} non-class CazyFamily nodes missing parent edge"


@pytest.mark.kg
def test_tcdb_full_text_index_exists(run_query):
    rows = run_query("SHOW INDEXES YIELD name WHERE name = 'tcdbFamilyFullText' RETURN name")
    assert len(rows) == 1, "tcdbFamilyFullText index missing"


@pytest.mark.kg
def test_cazy_full_text_index_exists(run_query):
    rows = run_query("SHOW INDEXES YIELD name WHERE name = 'cazyFamilyFullText' RETURN name")
    assert len(rows) == 1, "cazyFamilyFullText index missing"


@pytest.mark.kg
def test_transport_only_metabolites_exist(run_query):
    rows = run_query("""
        MATCH (m:Metabolite)
        WHERE 'transport' IN m.evidence_sources
          AND NOT 'metabolism' IN m.evidence_sources
        RETURN count(m) AS n
    """)
    assert rows[0]["n"] >= 100, f"Expected at least 100 transport-only Metabolites, found {rows[0]['n']}"


@pytest.mark.kg
def test_metabolites_have_evidence_sources(run_query):
    """Every Metabolite has at least one evidence source."""
    n_bad = run_query("""
        MATCH (m:Metabolite)
        WHERE m.evidence_sources IS NULL OR size(m.evidence_sources) = 0
        RETURN count(m) AS n
    """)[0]["n"]
    assert n_bad == 0, f"{n_bad} Metabolite nodes lack evidence_sources"
