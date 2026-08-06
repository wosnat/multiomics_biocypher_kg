"""KG validity: InterProScan InterPro-entry ontology (Phase-2 integration).

Validates the deployed graph — node presence/properties, the scored
Gene_has_interpro_entry edge, the is-a hierarchy, the Pfam bridge, post-import
rollups, the annotation_types fold-in, and no-orphan/no-dangling invariants.
See docs/kg-changes/interproscan-extension.md.
"""

import pytest

pytestmark = pytest.mark.kg

VALID_TYPES = {
    "FAMILY", "DOMAIN", "HOMOLOGOUS_SUPERFAMILY", "REPEAT",
    "CONSERVED_SITE", "ACTIVE_SITE", "BINDING_SITE", "PTM", "",
}


# ── nodes ─────────────────────────────────────────────────────────────────────

def test_interpro_entry_nodes_exist(run_query):
    n = run_query("MATCH (e:InterproEntry) RETURN count(e) AS n")[0]["n"]
    assert n > 5000, f"expected thousands of InterproEntry nodes, got {n}"


def test_interpro_entry_required_properties(run_query):
    """Every node has interpro_id + a valid type + a non-negative level."""
    bad = run_query("""
        MATCH (e:InterproEntry)
        WHERE e.interpro_id IS NULL OR NOT e.interpro_id STARTS WITH 'IPR'
           OR e.level IS NULL OR e.level < 0
           OR NOT e.interpro_type IN $types
        RETURN e.id AS id, e.interpro_type AS t LIMIT 5
    """, types=sorted(VALID_TYPES))
    assert bad == [], f"InterproEntry nodes with bad props: {bad}"


def test_interpro_node_id_uses_colon_prefix(run_query):
    """interpro is a registered bioregistry prefix → colon form."""
    row = run_query("MATCH (e:InterproEntry) RETURN e.id AS id LIMIT 1")[0]
    assert row["id"].startswith("interpro:IPR")


# ── Gene_has_interpro_entry edge (scored) ─────────────────────────────────────

def test_gene_has_interpro_entry_edges_exist(run_query):
    n = run_query("MATCH ()-[r:Gene_has_interpro_entry]->() RETURN count(r) AS n")[0]["n"]
    assert n > 100000, f"expected >100K gene→entry edges, got {n}"


def test_gene_interpro_edge_properties_sane(run_query):
    """match_count>=1, libraries non-empty, start<=end, evalue>=0 when present."""
    bad = run_query("""
        MATCH ()-[r:Gene_has_interpro_entry]->()
        WHERE r.match_count < 1
           OR r.libraries IS NULL OR size(r.libraries) = 0
           OR (r.start IS NOT NULL AND r.end IS NOT NULL AND r.start > r.end)
           OR (r.evalue IS NOT NULL AND r.evalue < 0)
        RETURN r.match_count AS mc, r.start AS s, r.end AS e LIMIT 5
    """)
    assert bad == [], f"edges with bad evidence props: {bad}"


def test_no_orphan_gene_interpro_edges(run_query):
    """Every Gene_has_interpro_entry endpoint resolves (no dangling)."""
    n = run_query("""
        MATCH (g)-[r:Gene_has_interpro_entry]->(e)
        WHERE NOT g:Gene OR NOT e:InterproEntry
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


# ── hierarchy ─────────────────────────────────────────────────────────────────

def test_hierarchy_edges_exist_and_resolve(run_query):
    n = run_query("MATCH (:InterproEntry)-[r:Interpro_entry_is_a_interpro_entry]->(:InterproEntry) RETURN count(r) AS n")[0]["n"]
    assert n > 500, f"expected hierarchy edges, got {n}"


def test_hierarchy_child_deeper_than_parent(run_query):
    """is-a child.level == parent.level + 1 (ParentChildTree depth)."""
    bad = run_query("""
        MATCH (c:InterproEntry)-[:Interpro_entry_is_a_interpro_entry]->(p:InterproEntry)
        WHERE c.level <> p.level + 1
        RETURN c.id AS child, c.level AS cl, p.level AS pl LIMIT 5
    """)
    assert bad == [], f"hierarchy level mismatches: {bad}"


# ── Pfam bridge ───────────────────────────────────────────────────────────────

def test_pfam_bridge_edges_resolve(run_query):
    """Pfam_in_interpro_entry: Pfam → InterproEntry, no dangling."""
    total = run_query("MATCH ()-[r:Pfam_in_interpro_entry]->() RETURN count(r) AS n")[0]["n"]
    assert total > 100, f"expected Pfam bridge edges, got {total}"
    bad = run_query("""
        MATCH (p)-[r:Pfam_in_interpro_entry]->(e)
        WHERE NOT p:Pfam OR NOT e:InterproEntry
        RETURN count(r) AS n
    """)[0]["n"]
    assert bad == 0


# ── post-import rollups ───────────────────────────────────────────────────────

def test_rollups_computed(run_query):
    """gene_count / organism_count / is_promiscuous set on every node post-import."""
    n = run_query("""
        MATCH (e:InterproEntry)
        WHERE e.gene_count IS NULL OR e.organism_count IS NULL OR e.is_promiscuous IS NULL
        RETURN count(e) AS n
    """)[0]["n"]
    assert n == 0, f"{n} InterproEntry nodes missing post-import rollups"


def test_is_promiscuous_matches_threshold(run_query):
    """is_promiscuous == (gene_count >= 1000)."""
    bad = run_query("""
        MATCH (e:InterproEntry)
        WHERE e.is_promiscuous <> (coalesce(e.gene_count,0) >= 1000)
        RETURN count(e) AS n
    """)[0]["n"]
    assert bad == 0


def test_gene_count_direct_matches_edges(run_query):
    """gene_count is DIRECT: equals distinct genes with an edge to the entry."""
    row = run_query("""
        MATCH (e:InterproEntry) WHERE e.gene_count > 0
        WITH e LIMIT 1
        MATCH (g:Gene)-[:Gene_has_interpro_entry]->(e)
        RETURN e.gene_count AS stored, count(DISTINCT g) AS actual
    """)[0]
    assert row["stored"] == row["actual"]


# ── functional fold-in (annotation_types only, NOT quality) ───────────────────

def test_annotation_types_includes_interpro(run_query):
    """Genes with an interpro edge carry 'interpro' in annotation_types."""
    n = run_query("""
        MATCH (g:Gene)-[:Gene_has_interpro_entry]->()
        WHERE NOT 'interpro' IN g.annotation_types
        RETURN count(DISTINCT g) AS n
    """)[0]["n"]
    assert n == 0


def test_interpro_not_in_informative_annotation_types(run_query):
    """Deferred by design: interpro must NOT appear in informative_annotation_types."""
    n = run_query("""
        MATCH (g:Gene)
        WHERE 'interpro' IN coalesce(g.informative_annotation_types, [])
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, "interpro leaked into informative_annotation_types (should be deferred)"


def test_interpro_entry_count_set(run_query):
    row = run_query("""
        MATCH (g:Gene)-[:Gene_has_interpro_entry]->() WITH g LIMIT 1
        RETURN g.interpro_entry_count AS c
    """)[0]
    assert row["c"] is not None and row["c"] >= 1
