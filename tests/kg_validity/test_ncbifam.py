"""KG validity: NCBIfam ontology (InterProScan Phase-2, Task 19 of the InterPro
multi-ontology redesign).

Layer-3 spot checks from
docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md
§6 — verified against committed calls.json / reference data before the
rebuild; this file asserts the same facts survive into the deployed graph.

See also tests/kg_validity/test_interpro.py (sibling InterPro-entry ontology).
"""

import pytest

pytestmark = pytest.mark.kg


# ── node counts / properties ────────────────────────────────────────────────

def test_ncbifam_family_count_in_range(run_query):
    """Spec §6 predicts ~4,957 NcbifamFamily nodes (observed-only, flat)."""
    n = run_query("MATCH (n:NcbifamFamily) RETURN count(n) AS n")[0]["n"]
    assert 4900 <= n <= 5050, f"expected NcbifamFamily count in [4900, 5050], got {n}"


def test_ncbifam_node_id_uses_underscore_prefix(run_query):
    """`ncbifam` is not a bioregistry prefix -> underscore fallback (psortb/signalp precedent)."""
    row = run_query("MATCH (n:NcbifamFamily) RETURN n.id AS id LIMIT 1")[0]
    assert row["id"].startswith("ncbifam_")


def test_tigr00198_node_is_equivalog_katg(run_query):
    """TIGR00198 (catalase-peroxidase katG) is a curated 'equivalog' family
    with a non-null description."""
    row = run_query(
        "MATCH (n:NcbifamFamily {ncbifam_id: 'TIGR00198'}) "
        "RETURN n.family_type AS family_type, n.description AS description, "
        "n.gene_symbol AS gene_symbol"
    )
    assert row, "TIGR00198 node not found"
    assert row[0]["family_type"] == "equivalog"
    assert row[0]["description"] is not None and row[0]["description"] != ""


# ── Gene_has_ncbifam_family edges ───────────────────────────────────────────

def test_gene_has_ncbifam_family_edges_exist(run_query):
    n = run_query("MATCH ()-[r:Gene_has_ncbifam_family]->() RETURN count(r) AS n")[0]["n"]
    assert n > 40000, f"expected tens of thousands of Gene_has_ncbifam_family edges, got {n}"


def test_no_orphan_gene_ncbifam_edges(run_query):
    n = run_query("""
        MATCH (g)-[r:Gene_has_ncbifam_family]->(n)
        WHERE NOT g:Gene OR NOT n:NcbifamFamily
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


# ── katG (Black Queen Hypothesis) at the NCBIfam layer ──────────────────────

def test_katg_present_ez55_and_mit1002(run_query):
    """TIGR00198 (katG) is called on Alteromonas EZ55 + MIT1002 genes."""
    rows = run_query("""
        MATCH (g:Gene)-[:Gene_has_ncbifam_family]->(:NcbifamFamily {ncbifam_id: 'TIGR00198'})
        MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        RETURN DISTINCT o.strain_name AS strain
    """)
    strains = {r["strain"] for r in rows}
    assert "EZ55" in strains, f"expected EZ55 in katG strains, got {strains}"
    assert "MIT1002" in strains, f"expected MIT1002 in katG strains, got {strains}"


def test_katg_absent_in_prochlorococcus(run_query):
    """Extends the existing Black Queen katG absence test to the NCBIfam layer."""
    n = run_query("""
        MATCH (g:Gene)-[:Gene_has_ncbifam_family]->(:NcbifamFamily {ncbifam_id: 'TIGR00198'})
        MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.genus = 'Prochlorococcus'
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"expected zero Prochlorococcus katG (TIGR00198) genes, got {n}"


# ── domain-inferred GO edge via the InterPro Layer-B path ──────────────────

def test_kt2440_mfs_domain_inferred_go(run_query):
    """KT2440 WP_010953880.1 (hypothetical protein + MFS domain, IPR020846)
    gets GO:0022857 (transmembrane transporter activity) with 'interpro' in
    sources — the DOMAIN-typed entry contributes domain_inferred-strength
    evidence via the Layer-B router.

    The design spec's Layer-3 spot check (§6) predicted the aggregate
    `evidence` field would read 'domain_inferred', on the assumption interpro
    was the sole contributor. In the actual merge, uniprot independently
    annotates this same (gene, GO) pair with 'curated' evidence
    (gene_annotations_merged.json: go_terms_source['GO:0022857'] =
    ['interpro', 'uniprot'], go_terms_evidence['GO:0022857'] = 'curated') —
    confirmed against cache/data/Pseudomonas/genomes/KT2440/. Since one edge
    carries the UNION of sources and `evidence` is the strongest label across
    them (curated > signature > family_inferred > domain_inferred, see
    CLAUDE.md), the edge legitimately reads 'curated'. The router fact this
    check exists to verify — interpro contributing to the edge — still holds
    via `sources`."""
    rows = run_query("""
        MATCH (g:Gene {protein_id: 'WP_010953880.1'})-[r:Gene_enables_molecular_function]->(t:MolecularFunction {id: 'go:0022857'})
        RETURN r.evidence AS evidence, r.sources AS sources
    """)
    assert rows, "expected a Gene_enables_molecular_function edge to GO:0022857 for WP_010953880.1"
    assert rows[0]["evidence"] in ("domain_inferred", "curated")
    assert "interpro" in (rows[0]["sources"] or [])


# ── informativeness (F1.1, ncbifam property-valued rule) ────────────────────

def test_hypoth_equivalog_flagged_uninformative(run_query):
    n = run_query("""
        MATCH (n:NcbifamFamily {family_type: 'hypoth_equivalog'})
        WHERE n.is_uninformative = 'true'
        RETURN count(n) AS n
    """)[0]["n"]
    assert n >= 1, "expected at least one hypoth_equivalog NcbifamFamily flagged is_uninformative"


def test_duf_named_interpro_entry_flagged(run_query):
    """DUF/Domain-of-unknown-function-named InterproEntry nodes are flagged
    is_uninformative (interpro_entry name_patterns rule, F1.1)."""
    n = run_query("""
        MATCH (e:InterproEntry)
        WHERE e.name STARTS WITH 'Domain of unknown function'
          AND e.is_uninformative = 'true'
        RETURN count(e) AS n
    """)[0]["n"]
    assert n >= 1, "expected at least one DUF-named InterproEntry flagged is_uninformative"
