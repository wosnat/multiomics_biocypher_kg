"""
KG-validity tests for publication "discusses" literature-index edges.

Validates the Publication_discusses_gene / Publication_discusses_kegg_pathway
edges produced by publication_topics_adapter:
- edges exist and connect the correct node types (no dangling/wrong endpoints)
- prominence is a valid enum; evidence is present
- pathway targets are pathway-level KeggTerm nodes
- post-import rollups are computed

Auto-skips when Neo4j is unreachable (run against a built Docker graph).
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Edge presence + endpoint correctness (no dangling)
# ---------------------------------------------------------------------------


def test_discuss_gene_edges_exist(run_query):
    cnt = run_query(
        "MATCH ()-[r:Publication_discusses_gene]->() RETURN count(r) AS cnt"
    )[0]["cnt"]
    assert cnt > 100, f"Only {cnt} Publication_discusses_gene edges; expected >100"


def test_discuss_pathway_edges_exist(run_query):
    cnt = run_query(
        "MATCH ()-[r:Publication_discusses_kegg_pathway]->() RETURN count(r) AS cnt"
    )[0]["cnt"]
    assert cnt > 20, f"Only {cnt} Publication_discusses_kegg_pathway edges; expected >20"


def test_discuss_gene_endpoints_correct(run_query):
    """Source must be Publication, target must be Gene (no dangling/wrong-type)."""
    bad = run_query("""
        MATCH (src)-[r:Publication_discusses_gene]->(tgt)
        WHERE NOT src:Publication OR NOT tgt:Gene
        RETURN count(r) AS bad
    """)[0]["bad"]
    assert bad == 0, f"{bad} Publication_discusses_gene edges have wrong endpoints"


def test_discuss_pathway_endpoints_correct(run_query):
    """Source Publication, target a pathway-level KeggTerm."""
    bad = run_query("""
        MATCH (src)-[r:Publication_discusses_kegg_pathway]->(tgt)
        WHERE NOT src:Publication OR NOT tgt:KeggTerm
        RETURN count(r) AS bad
    """)[0]["bad"]
    assert bad == 0, f"{bad} Publication_discusses_kegg_pathway edges have wrong endpoints"

    non_pathway = run_query("""
        MATCH (:Publication)-[:Publication_discusses_kegg_pathway]->(k:KeggTerm)
        WHERE k.level_kind <> 'pathway'
        RETURN count(*) AS cnt
    """)[0]["cnt"]
    assert non_pathway == 0, f"{non_pathway} pathway edges target non-pathway KeggTerm nodes"


# ---------------------------------------------------------------------------
# Edge properties
# ---------------------------------------------------------------------------


def test_discuss_prominence_enum(run_query):
    bad = run_query("""
        MATCH ()-[r:Publication_discusses_gene|Publication_discusses_kegg_pathway]->()
        WHERE NOT r.prominence IN ['central', 'peripheral']
        RETURN count(r) AS bad
    """)[0]["bad"]
    assert bad == 0, f"{bad} discuss-edges have prominence outside {{central, peripheral}}"


def test_discuss_evidence_present(run_query):
    """Every discuss-edge carries an evidence property (string, may be empty)."""
    missing = run_query("""
        MATCH ()-[r:Publication_discusses_gene|Publication_discusses_kegg_pathway]->()
        WHERE r.evidence IS NULL
        RETURN count(r) AS missing
    """)[0]["missing"]
    assert missing == 0, f"{missing} discuss-edges missing the evidence property"


# ---------------------------------------------------------------------------
# Post-import rollups
# ---------------------------------------------------------------------------


def test_publication_discuss_rollups(run_query):
    """Publications that discuss genes have discussed_gene_count > 0 and consistent."""
    rows = run_query("""
        MATCH (p:Publication)
        WHERE EXISTS { (p)-[:Publication_discusses_gene]->() }
        WITH p, count { (p)-[:Publication_discusses_gene]->() } AS actual
        RETURN count(*) AS pubs,
               sum(CASE WHEN p.discussed_gene_count = actual THEN 1 ELSE 0 END) AS ok
    """)[0]
    assert rows["pubs"] > 0, "No publications with discuss-gene edges found"
    assert rows["ok"] == rows["pubs"], (
        f"discussed_gene_count mismatch on {rows['pubs'] - rows['ok']} publications"
    )


def test_gene_discussed_in_publication_count(run_query):
    """At least some genes carry a positive discussed_in_publication_count."""
    cnt = run_query("""
        MATCH (g:Gene)
        WHERE coalesce(g.discussed_in_publication_count, 0) > 0
        RETURN count(g) AS cnt
    """)[0]["cnt"]
    assert cnt > 50, f"Only {cnt} genes with discussed_in_publication_count > 0"
