"""KG-SYNC-005 (2026-08-27) — annotation-trust surface acceptance criteria.

Driver: multiomics_explorer/docs/kg-specs/2026-08-27-annotation-trust-kg-asks.md
(ONT-001…015). Plan: plans/kg_sync_005_annotation_trust.md.
"""
import pytest

pytestmark = pytest.mark.kg

GENE_ONTOLOGY_EDGES = [
    "Gene_involved_in_biological_process", "Gene_enables_molecular_function",
    "Gene_located_in_cellular_component", "Gene_catalyzes_ec_number",
    "Gene_has_pfam", "Gene_has_cazy_family", "Gene_has_tcdb_family",
    "Gene_has_kegg_ko", "Gene_in_cog_category", "Gene_has_cyanorak_role",
    "Gene_has_tigr_role", "Gene_has_interpro_entry", "Gene_has_ncbifam_family",
    "Gene_has_merops_family",
]
LADDER = {"curated", "signature", "homology", "family_inferred", "domain_inferred"}

HIERARCHICAL = [  # label, hierarchy rel(s) — subtree gene_count + direct_gene_count
    "BiologicalProcess", "MolecularFunction", "CellularComponent", "EcNumber",
    "KeggTerm", "CyanorakRole", "InterproEntry", "TcdbFamily", "CazyFamily",
    "MeropsFamily", "TigrRole",
]
FLAT = ["Pfam", "PfamClan", "CogFunctionalCategory", "NcbifamFamily",
        "SubcellularLocalization", "SignalPeptideType", "BriteCategory"]


@pytest.mark.parametrize("edge", GENE_ONTOLOGY_EDGES)
def test_sources_and_evidence_on_every_edge(run_query, edge):
    """AC1 (ONT-007/008): 100% of every gene→ontology edge type carries both."""
    row = run_query(
        f"MATCH ()-[r:{edge}]->() RETURN count(r) AS n, count(r.sources) AS s, "
        f"count(r.evidence) AS e, collect(DISTINCT r.evidence) AS values"
    )[0]
    assert row["n"] > 0
    assert row["s"] == row["n"], f"{edge}: sources on {row['s']}/{row['n']}"
    assert row["e"] == row["n"], f"{edge}: evidence on {row['e']}/{row['n']}"
    assert set(row["values"]) <= LADDER, row["values"]


def test_evidence_rungs_by_edge(run_query):
    for edge, expected in {
        "Gene_has_kegg_ko": {"family_inferred"},
        "Gene_in_cog_category": {"family_inferred"},
        "Gene_has_cyanorak_role": {"curated"},
        "Gene_has_tigr_role": {"curated", "family_inferred"},
        "Gene_has_interpro_entry": {"signature"},
        "Gene_has_ncbifam_family": {"signature"},
        "Gene_has_merops_family": {"homology"},
        "Gene_has_tcdb_family": {"homology", "family_inferred"},
    }.items():
        vals = run_query(f"MATCH ()-[r:{edge}]->() RETURN collect(DISTINCT r.evidence) AS v")[0]["v"]
        assert set(vals) == expected, (edge, vals)


def test_sources_join_data_source_nodes(run_query):
    """R2: every sources value on every edge type resolves to a DataSource node."""
    for edge in GENE_ONTOLOGY_EDGES:
        rows = run_query(
            f"MATCH ()-[r:{edge}]->() UNWIND r.sources AS s WITH DISTINCT s "
            f"OPTIONAL MATCH (d:DataSource {{id: 'data_source:' + s}}) "
            f"RETURN s, d IS NOT NULL AS ok"
        )
        bad = [r["s"] for r in rows if not r["ok"]]
        assert not bad, (edge, bad)


def test_tcdb_evidence_matches_sources(run_query):
    rows = run_query(
        "MATCH ()-[r:Gene_has_tcdb_family]->() "
        "RETURN 'tcdb_diamond' IN r.sources AS diamond, r.evidence AS ev, count(*) AS n"
    )
    for r in rows:
        assert r["ev"] == ("homology" if r["diamond"] else "family_inferred"), r


def test_merops_evidence_score(run_query):
    """AC2 (ONT-006): 2-signal score on every MEROPS edge; round(score*2) = fired."""
    row = run_query(
        "MATCH ()-[r:Gene_has_merops_family]->() "
        "RETURN count(r) AS n, count(r.evidence_score) AS scored, "
        "collect(DISTINCT r.evidence_score) AS vals, "
        "sum(CASE WHEN round(r.evidence_score * 2) = "
        "  (CASE WHEN r.tier <= 2 THEN 1 ELSE 0 END) + "
        "  (CASE WHEN r.pfam_support = 'corroborated' THEN 1 ELSE 0 END) THEN 1 ELSE 0 END) AS consistent"
    )[0]
    assert row["scored"] == row["n"] > 0
    assert set(row["vals"]) <= {0.0, 0.5, 1.0}, row["vals"]
    assert row["consistent"] == row["n"]


def test_merops_evidence_score_max_sparse(run_query):
    row = run_query(
        "MATCH (g:Gene) WITH g, EXISTS { (g)-[:Gene_has_merops_family]->() } AS has "
        "RETURN sum(CASE WHEN has AND g.merops_evidence_score_max IS NULL THEN 1 ELSE 0 END) AS missing, "
        "sum(CASE WHEN NOT has AND g.merops_evidence_score_max IS NOT NULL THEN 1 ELSE 0 END) AS spurious"
    )[0]
    assert row["missing"] == 0 and row["spurious"] == 0


def test_attachment_depth_materialized(run_query):
    """AC3 (ONT-010): property on 100% of edges and equal to the inline predicate."""
    row = run_query(
        "MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily) "
        "WITH r, EXISTS { MATCH (g)-[:Gene_has_tcdb_family]->(d:TcdbFamily) "
        "  WHERE (d)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t) } AS sup "
        "RETURN count(r) AS n, count(r.attachment_depth) AS set_, "
        "sum(CASE WHEN r.attachment_depth = (CASE WHEN sup THEN 'superseded' ELSE 'most_specific' END) "
        "    THEN 1 ELSE 0 END) AS agree, collect(DISTINCT r.attachment_depth) AS vals"
    )[0]
    assert row["set_"] == row["n"] == row["agree"] > 0
    assert set(row["vals"]) == {"most_specific", "superseded"}


def test_attachment_depth_med4_spot_check(run_query):
    row = run_query(
        "MATCH (g:Gene {organism_name: 'Prochlorococcus MED4'})-[r:Gene_has_tcdb_family]->() "
        "RETURN count(r) AS rows, sum(CASE WHEN r.attachment_depth = 'superseded' THEN 1 ELSE 0 END) AS sup"
    )[0]
    assert row["rows"] == 670 and row["sup"] == 73


def test_transport_counts_read_most_specific_only(run_query):
    """The gene scalar and the metabolite scalar are two projections of one set."""
    a = run_query(
        "MATCH (g:Gene)-[:Gene_has_tcdb_family {attachment_depth: 'most_specific'}]->(:TcdbFamily)"
        "-[:Tcdb_family_transports_metabolite]->(m:Metabolite) "
        "RETURN count(DISTINCT [g, m]) AS pairs"
    )[0]["pairs"]
    b = run_query("MATCH (g:Gene) RETURN sum(coalesce(g.transported_metabolite_count, 0)) AS s")[0]["s"]
    c = run_query("MATCH (m:Metabolite) RETURN sum(coalesce(m.transporter_gene_count, 0)) AS s")[0]["s"]
    assert a == b == c


@pytest.mark.parametrize("label", HIERARCHICAL)
def test_hierarchical_counts(run_query, label):
    """AC4 (ONT-009/015): subtree gene_count >= direct_gene_count, both present."""
    row = run_query(
        f"MATCH (n:{label}) RETURN count(n) AS n, count(n.gene_count) AS gc, "
        f"count(n.direct_gene_count) AS dgc, count(n.organism_count) AS oc, "
        f"sum(CASE WHEN n.gene_count < n.direct_gene_count THEN 1 ELSE 0 END) AS bad, "
        f"sum(CASE WHEN n.gene_count > n.direct_gene_count THEN 1 ELSE 0 END) AS rolled"
    )[0]
    assert row["gc"] == row["dgc"] == row["oc"] == row["n"] > 0, row
    assert row["bad"] == 0
    assert row["rolled"] > 0, f"{label}: subtree never exceeds direct — rollup not applied?"


@pytest.mark.parametrize("label", FLAT)
def test_flat_counts(run_query, label):
    row = run_query(
        f"MATCH (n:{label}) RETURN count(n) AS n, count(n.gene_count) AS gc, count(n.organism_count) AS oc"
    )[0]
    assert row["gc"] == row["oc"] == row["n"] > 0, row


def test_interpro_direct_count_is_the_old_direct_semantics(run_query):
    row = run_query(
        "MATCH (e:InterproEntry) OPTIONAL MATCH (e)<-[:Gene_has_interpro_entry]-(g:Gene) "
        "WITH e, count(DISTINCT g) AS live RETURN sum(CASE WHEN live <> e.direct_gene_count THEN 1 ELSE 0 END) AS bad"
    )[0]
    assert row["bad"] == 0


def test_merops_peptidase_organism_count(run_query):
    row = run_query(
        "MATCH (m:MeropsFamily) RETURN count(m) AS n, count(m.peptidase_organism_count) AS c, "
        "sum(CASE WHEN m.peptidase_organism_count > m.organism_count THEN 1 ELSE 0 END) AS bad"
    )[0]
    assert row["c"] == row["n"] and row["bad"] == 0


def test_ncbifam_bit_score_and_retired(run_query):
    """AC5 (ONT-004/011)."""
    e = run_query(
        "MATCH ()-[r:Gene_has_ncbifam_family]->() RETURN count(r) AS n, count(r.bit_score) AS bs, "
        "count(r.score) AS old, count(r.match_count) AS mc"
    )[0]
    assert e["bs"] == e["n"] > 0 and e["old"] == 0 and e["mc"] == 0
    n = run_query(
        "MATCH (f:NcbifamFamily) RETURN sum(CASE WHEN f.family_type IS NULL THEN 1 ELSE 0 END) AS nulls, "
        "sum(CASE WHEN f.family_type = 'retired' THEN 1 ELSE 0 END) AS retired"
    )[0]
    assert n["nulls"] == 0 and n["retired"] > 0


def test_merops_family_class_renamed(run_query):
    row = run_query(
        "MATCH (m:MeropsFamily) RETURN count(m) AS n, count(m.family_class) AS fc, count(m.family_type) AS old, "
        "collect(DISTINCT m.family_class) AS vals"
    )[0]
    assert row["fc"] == row["n"] and row["old"] == 0
    assert set(row["vals"]) == {"peptidase", "inhibitor"}


EXPLORER_VOCAB_PAIRS = [  # §7.3 of the asks doc — every (applies_to, property) the explorer asserts
    *[(e, p) for e in GENE_ONTOLOGY_EDGES for p in ("sources", "evidence")],
    *[(e, "evidence_score") for e in GENE_ONTOLOGY_EDGES if e not in (
        "Gene_has_kegg_ko", "Gene_in_cog_category", "Gene_has_cyanorak_role", "Gene_has_tigr_role",
        "Gene_has_interpro_entry", "Gene_has_ncbifam_family")],
    *[(e, p) for e in ("Gene_has_tcdb_family", "Gene_has_merops_family")
      for p in ("tier", "confidence_score", "identity", "qcov", "evalue", "consensus_n", "pfam_support")],
    *[("Gene_has_tcdb_family", p) for p in ("source_agreement", "go_support", "attachment_depth")],
    *[("Gene_has_merops_family", p) for p in ("call_class", "best_hit_kind")],
    *[("Gene_has_interpro_entry", p) for p in ("libraries", "evalue_library", "evalue", "match_count")],
    *[("Gene_has_ncbifam_family", p) for p in ("evalue", "bit_score")],
    ("Gene_has_subcellular_localization", "score"),
    *[("Gene_has_signal_peptide_type", p) for p in ("probability", "cleavage_site", "cleavage_probability")],
    ("InterproEntry", "interpro_type"), ("NcbifamFamily", "family_type"), ("MeropsFamily", "family_class"),
    ("MeropsFamily", "catalytic_type"), ("MeropsFamily", "level_kind"), ("TcdbFamily", "level_kind"),
    ("CazyFamily", "level_kind"), ("BriteCategory", "tree"), ("Gene", "merops_classes"),
]


def test_explorer_vocab_pairs_have_nodes(run_query):
    """AC6 (ONT-014): §7.3 list fully covered by ControlledVocabulary nodes."""
    ids = {r["id"] for r in run_query("MATCH (v:ControlledVocabulary) RETURN v.id AS id")}
    missing = [f"{a}.{p}" for a, p in EXPLORER_VOCAB_PAIRS if f"{a}.{p}" not in ids]
    assert not missing, missing


# ── DOC-001 / DOC-006 (explorer docs review, 2026-08-29) ─────────────────────

MERGED_ANNOTATION_EDGES = [
    "Gene_involved_in_biological_process", "Gene_enables_molecular_function",
    "Gene_located_in_cellular_component", "Gene_catalyzes_ec_number",
    "Gene_has_pfam", "Gene_has_cazy_family",
]


@pytest.mark.parametrize("edge", MERGED_ANNOTATION_EDGES)
def test_eggnog_only_edges_are_family_inferred(run_query, edge):
    """DOC-001: eggNOG orthology transfer is `family_inferred` on every edge
    type, never `curated` — the ladder rung must not depend on which adapter
    emitted the edge. Conversely a curated source (ncbi/cyanorak/uniprot)
    always reads `curated`."""
    rows = run_query(
        f"MATCH ()-[r:{edge}]->() WHERE r.sources = ['eggnog'] "
        f"RETURN r.evidence AS e, count(*) AS n"
    )
    assert rows, f"{edge}: no eggNOG-only edges?"
    assert {r["e"] for r in rows} == {"family_inferred"}, rows
    bad = run_query(
        f"MATCH ()-[r:{edge}]->() WHERE any(s IN r.sources WHERE s IN "
        f"['ncbi','cyanorak','uniprot']) AND r.evidence <> 'curated' "
        f"RETURN count(*) AS n"
    )[0]["n"]
    assert bad == 0, f"{edge}: {bad} curated-backed edges not 'curated'"


def test_pfam_and_cazy_agree_on_the_eggnog_interproscan_pair(run_query):
    """The same source pair used to read `signature` on Pfam and `curated` on
    CAZy (a key-alignment artefact in the step-2 merge). Now: Pfam keeps the
    direct-hit rung; CAZy/EC/GO take eggNOG's family-level floor."""
    pf = run_query("MATCH ()-[r:Gene_has_pfam]->() WHERE r.sources = ['eggnog','interproscan'] "
                   "RETURN collect(DISTINCT r.evidence) AS v")[0]["v"]
    assert set(pf) == {"signature"}, pf
    for edge in ("Gene_has_cazy_family", "Gene_catalyzes_ec_number",
                 "Gene_involved_in_biological_process"):
        v = run_query(f"MATCH ()-[r:{edge}]->() WHERE r.sources = ['eggnog','interproscan'] "
                      "RETURN collect(DISTINCT r.evidence) AS v")[0]["v"]
        assert set(v) <= {"family_inferred"}, (edge, v)


def test_kegg_direct_gene_count_only_on_kos(run_query):
    """DOC-006: genes attach to KOs only, so direct_gene_count is omitted (not
    stored as 0) on pathway / subcategory / category KeggTerm nodes —
    BriteCategory / PfamClan precedent."""
    rows = run_query("MATCH (t:KeggTerm) RETURN t.level_kind AS k, count(*) AS n, "
                     "count(t.direct_gene_count) AS dgc")
    by = {r["k"]: r for r in rows}
    assert by["ko"]["dgc"] == by["ko"]["n"]
    for k in ("pathway", "subcategory", "category"):
        assert by[k]["dgc"] == 0, (k, by[k])
