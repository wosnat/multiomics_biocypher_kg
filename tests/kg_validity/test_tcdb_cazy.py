"""KG validity tests for the TCDB and CAZy ontologies.

TCDB part follows the TCDB cleanup (project 3a):
docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md — §9 "Acceptance" is the
checklist the tests below encode. Magnitudes were measured by running the adapter
on the 2026-10-23 caches (9,389 TcdbFamily nodes; 75,823 gene edges); bands are
used instead of exact values except where the spec demands exactness.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
TCDB_PRUNED = REPO / "cache/data/tcdb/tcdb_pruned.json"
TCDB_HIERARCHY = REPO / "cache/data/tcdb/tcdb_hierarchy.json"

GATE_PASSING = ["direct", "annotated"]
TCDB_CALLS = {"eggnog", "diamond_tier", "diamond_best_hit"}
TCDB_GATES = {"direct", "annotated", "below_gate"}
LEVEL_KINDS = ["tc_class", "tc_subclass", "tc_family", "tc_subfamily", "tc_system"]



def _unnamed_systems_in_caches() -> int:
    """Kept tc_system ids whose committed hierarchy name is missing or equals the id —
    the same count step 6 logs ("134 kept TCDB systems have no name", 2026-10-23).
    Those render as the bare TCID by design (unnamed upstream)."""
    for p in (TCDB_PRUNED, TCDB_HIERARCHY):
        if not p.exists():
            pytest.skip(f"{p.relative_to(REPO)} missing; run prepare_data step 6")
    kept = set(json.loads(TCDB_PRUNED.read_text())["kept_tcdb_ids"])
    hier = json.loads(TCDB_HIERARCHY.read_text())
    return sum(1 for t in kept
               if hier.get(t, {}).get("level_kind") == "tc_system"
               and (not hier[t].get("name") or hier[t]["name"] == t))

# The tcdb_gate rule recomputed from the edge's own call kind / rank / tier and the
# target's level (spec §5.3). Kept as Cypher fragments so several tests share one
# definition that is independent of the stored tcdb_gate value.
_PASSES = ("(r.tcdb_call = 'eggnog'"
           " OR (r.tcdb_call = 'diamond_tier' AND coalesce(r.tier <= 2, false))"
           " OR (r.tcdb_call = 'diamond_best_hit' AND coalesce(r.rank = 1 AND r.tier <= 2, false)))")
_DIRECT_ARM = ("(r.tcdb_call = 'eggnog'"
               " OR (r.tcdb_call = 'diamond_tier' AND coalesce(r.tier = 1, false))"
               " OR (r.tcdb_call = 'diamond_best_hit' AND coalesce(r.rank = 1 AND r.tier <= 2, false)))")
_IS_DIRECT = f"(t.level_kind = 'tc_system' AND {_DIRECT_ARM})"


# ============================================================================
# Node / edge counts
# ============================================================================


@pytest.mark.kg
def test_tcdb_family_node_count_in_range(run_query):
    """Kept set = gene-anchored seeds (eggNOG + diamond tier calls + best-hit systems)
    + ancestors + every substrate-bearing system under a kept family (§5.4).
    Measured 9,389 on the 2026-10-23 caches (was 1,515 before the cleanup)."""
    n = run_query("MATCH (t:TcdbFamily) RETURN count(t) AS n")[0]["n"]
    assert 8500 <= n <= 10500, f"TcdbFamily count {n} outside expected 8,500-10,500"


@pytest.mark.kg
def test_range_only_systems_exist_with_zero_gene_count(run_query):
    """Range-only systems (§5.4 step 3) carry substrates but no gate-passing gene."""
    n = run_query("""
        MATCH (t:TcdbFamily {level_kind: 'tc_system'})
        WHERE t.gene_count = 0 AND t.tcdb_substrate_count > 0
        RETURN count(t) AS n""")[0]["n"]
    assert n > 1000, f"only {n} substrate-bearing systems with gene_count = 0"


@pytest.mark.kg
def test_cazy_family_node_count_in_range(run_query):
    """Observed-only — small ontology."""
    n = run_query("MATCH (c:CazyFamily) RETURN count(c) AS n")[0]["n"]
    assert 5 <= n <= 100, f"CazyFamily count {n} outside expected 5-100"


@pytest.mark.kg
def test_gene_has_tcdb_family_edge_count(run_query):
    """One edge per (gene, TC node, tcdb_call) — measured 75,823."""
    n = run_query("MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n")[0]["n"]
    assert 68000 <= n <= 84000, f"Gene_has_tcdb_family count {n} outside 68,000-84,000"


@pytest.mark.kg
def test_gene_edge_counts_by_tcdb_call_and_gate(run_query):
    """Measured: eggnog 17,069 / diamond_tier 41,242 / diamond_best_hit 17,512;
    gate direct 11,158 / annotated 19,866 / below_gate 44,799. The exact per-call
    match against the merge inputs lives in test_tcdb_best_hit.py."""
    calls = {r["k"]: r["n"] for r in run_query(
        "MATCH ()-[r:Gene_has_tcdb_family]->() RETURN r.tcdb_call AS k, count(*) AS n")}
    gates = {r["k"]: r["n"] for r in run_query(
        "MATCH ()-[r:Gene_has_tcdb_family]->() RETURN r.tcdb_gate AS k, count(*) AS n")}
    assert set(calls) == TCDB_CALLS, calls
    assert set(gates) == TCDB_GATES, gates
    bands = {"eggnog": (15000, 20000), "diamond_tier": (37000, 46000),
             "diamond_best_hit": (15500, 19500),
             "direct": (10000, 12500), "annotated": (18000, 22000), "below_gate": (40000, 50000)}
    for k, n in {**calls, **gates}.items():
        lo, hi = bands[k]
        assert lo <= n <= hi, f"{k}: {n} outside {lo}-{hi}"


@pytest.mark.kg
def test_gene_has_cazy_family_edge_count(run_query):
    """Upper bound raised from 2000 by the InterPro two-layer integration
    (2026-08-10): Layer B adds ~642 CAZy ids from FAMILY/DOMAIN entry xrefs."""
    n = run_query("MATCH ()-[r:Gene_has_cazy_family]->() RETURN count(r) AS n")[0]["n"]
    assert 100 <= n <= 5000, f"Gene_has_cazy_family count {n} outside 100-5000"


@pytest.mark.kg
def test_tcdb_substrate_edge_counts(run_query):
    """Measured 10,778 metabolite + 367 role edges, all from tc_system nodes."""
    met = run_query("MATCH ()-[r:Tcdb_family_transports_metabolite]->() RETURN count(r) AS n")[0]["n"]
    role = run_query("MATCH ()-[r:Tcdb_family_transports_chemical_role]->() RETURN count(r) AS n")[0]["n"]
    assert 9500 <= met <= 12500, f"Tcdb_family_transports_metabolite {met} outside 9,500-12,500"
    assert 300 <= role <= 450, f"Tcdb_family_transports_chemical_role {role} outside 300-450"


# ============================================================================
# §9: substrate edges leave systems only
# ============================================================================


@pytest.mark.kg
def test_substrate_edges_leave_systems_only(run_query):
    row = run_query("""
        MATCH (t:TcdbFamily)-[r:Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role]->()
        RETURN count(r) AS total,
               sum(CASE WHEN t.level_kind <> 'tc_system' THEN 1 ELSE 0 END) AS off_system""")[0]
    assert row["total"] > 0
    assert row["off_system"] == 0, f"{row['off_system']} substrate/role edges leave a non-system node"


@pytest.mark.kg
def test_substrate_edges_carry_no_properties(run_query):
    """substrate_depth is gone; BioCypher writes only `id` on relationships."""
    rows = run_query("""
        MATCH ()-[r:Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role]->()
        UNWIND keys(r) AS k WITH DISTINCT k WHERE k <> 'id' RETURN collect(k) AS ks""")
    assert rows[0]["ks"] == [], f"substrate edges carry properties: {rows[0]['ks']}"


# ============================================================================
# §9: gene edge shape
# ============================================================================


@pytest.mark.kg
def test_gene_edge_shape_by_tcdb_call(run_query):
    """One tcdb_call and a one-element sources per edge; each kind carries only its
    own properties; exactly one rank = 1 per gene with best hits; tcdb_gate equals
    the rule recomputed in Cypher."""
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        RETURN count(r) AS n,
          sum(CASE WHEN r.tcdb_call IS NULL OR NOT r.tcdb_call IN $calls THEN 1 ELSE 0 END) AS bad_call,
          sum(CASE WHEN r.sources IS NULL OR size(r.sources) <> 1 THEN 1 ELSE 0 END) AS bad_sources_len,
          sum(CASE WHEN r.sources <> CASE r.tcdb_call WHEN 'eggnog' THEN ['eggnog']
                                                      ELSE ['tcdb_diamond'] END THEN 1 ELSE 0 END) AS bad_sources,
          sum(CASE WHEN r.evidence <> CASE r.tcdb_call WHEN 'eggnog' THEN 'family_inferred'
                                                       ELSE 'homology' END THEN 1 ELSE 0 END) AS bad_evidence,
          sum(CASE WHEN r.tcdb_call = 'eggnog' AND (r.tier IS NOT NULL OR r.identity IS NOT NULL
                     OR r.qcov IS NOT NULL OR r.scov IS NOT NULL OR r.evalue IS NOT NULL
                     OR r.rank IS NOT NULL OR r.bitscore IS NOT NULL OR r.bitscore_ratio IS NOT NULL)
                   THEN 1 ELSE 0 END) AS eggnog_leak,
          sum(CASE WHEN r.tcdb_call <> 'eggnog' AND (r.tier IS NULL OR NOT r.tier IN [1, 2, 3])
                   THEN 1 ELSE 0 END) AS diamond_no_tier,
          sum(CASE WHEN r.tcdb_call = 'diamond_tier' AND (r.rank IS NOT NULL OR r.bitscore IS NOT NULL
                     OR r.bitscore_ratio IS NOT NULL) THEN 1 ELSE 0 END) AS tier_leak,
          sum(CASE WHEN r.tcdb_call = 'diamond_best_hit' AND (r.rank IS NULL OR r.bitscore IS NULL
                     OR r.bitscore_ratio IS NULL OR r.identity IS NULL) THEN 1 ELSE 0 END) AS bh_missing,
          sum(CASE WHEN r.tcdb_call <> 'eggnog' AND (r.identity IS NULL OR r.qcov IS NULL
                     OR r.scov IS NULL OR r.evalue IS NULL) THEN 1 ELSE 0 END) AS diamond_stats_missing,
          sum(CASE WHEN r.tcdb_gate IS NULL OR NOT r.tcdb_gate IN $gates THEN 1 ELSE 0 END) AS bad_gate
    """, calls=sorted(TCDB_CALLS), gates=sorted(TCDB_GATES))[0]
    assert row["n"] > 0
    bad = {k: v for k, v in row.items() if k != "n" and v}
    assert not bad, bad

    # Best-hit statistics: rank 1..3, bitscore_ratio 1.0 at rank 1 and >= 0.9 for
    # near-ties, tier = identity band (>=70 -> 1, >=40 -> 2, else 3).
    bh = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->()
        RETURN sum(CASE WHEN r.rank < 1 OR r.rank > 3 THEN 1 ELSE 0 END) AS bad_rank,
               sum(CASE WHEN r.rank = 1 AND abs(r.bitscore_ratio - 1.0) > 1e-9 THEN 1 ELSE 0 END) AS bad_top,
               sum(CASE WHEN r.rank > 1 AND (r.bitscore_ratio < 0.9 - 1e-6 OR r.bitscore_ratio > 1.0 + 1e-9)
                        THEN 1 ELSE 0 END) AS bad_ratio,
               sum(CASE WHEN r.tier <> CASE WHEN r.identity >= 70 THEN 1
                                           WHEN r.identity >= 40 THEN 2 ELSE 3 END
                        THEN 1 ELSE 0 END) AS bad_band""")[0]
    assert not any(bh.values()), bh

    n_rank1 = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->()
        WITH g, sum(CASE WHEN r.rank = 1 THEN 1 ELSE 0 END) AS k WHERE k <> 1
        RETURN count(g) AS n""")[0]["n"]
    assert n_rank1 == 0, f"{n_rank1} genes with best hits do not have exactly one rank = 1 edge"
    dup_rank = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->()
        WITH g, r.rank AS rk, count(*) AS c WHERE c > 1
        RETURN count(DISTINCT g) AS n""")[0]["n"]
    assert dup_rank == 0, f"{dup_rank} genes with two best-hit edges of the same rank"

    mismatch = run_query(f"""
        MATCH (:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
        WITH r, CASE WHEN NOT {_PASSES} THEN 'below_gate'
                     WHEN {_IS_DIRECT} THEN 'direct'
                     ELSE 'annotated' END AS expected
        WHERE r.tcdb_gate <> expected
        RETURN count(r) AS n, collect([r.tcdb_call, r.tcdb_gate, expected])[..5] AS ex""")[0]
    assert mismatch["n"] == 0, f"{mismatch['n']} edges whose tcdb_gate disagrees with the rule: {mismatch['ex']}"


@pytest.mark.kg
def test_one_edge_per_gene_tc_and_call(run_query):
    """Edge ids are {locus_tag}-has_tcdb-{tc_id}-{tcdb_call}: at most one edge per
    (gene, TC node, call kind); up to three kinds may share a (gene, node)."""
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
        WITH g, t, r.tcdb_call AS k, count(r) AS c WHERE c > 1
        RETURN count(*) AS n""")[0]["n"]
    assert n == 0, f"{n} (gene, TC, tcdb_call) triples have parallel edges"
    shared = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
        WITH g, t, count(r) AS c WHERE c > 1 RETURN count(*) AS n""")[0]["n"]
    assert shared > 0, "no (gene, node) pair carries edges from two call kinds — sources merged again?"


@pytest.mark.kg
def test_sources_join_data_source_nodes(run_query):
    rows = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->() UNWIND r.sources AS s WITH DISTINCT s
        OPTIONAL MATCH (d:DataSource {id: 'data_source:' + s})
        RETURN s, d IS NOT NULL AS ok""")
    assert {r["s"] for r in rows} == {"eggnog", "tcdb_diamond"}
    assert all(r["ok"] for r in rows), rows


# ============================================================================
# §9 / §7: annotation gate and the transport arm
# ============================================================================


@pytest.mark.kg
def test_annotation_types_tcdb_uses_the_gate(run_query):
    """'tcdb' in annotation_types (and 'transporter' in informative_annotation_types)
    iff the gene has a gate-passing edge (tcdb_gate IN direct/annotated). Both
    directions."""
    rows = run_query("""
        MATCH (g:Gene)
        WITH g, EXISTS { (g)-[r:Gene_has_tcdb_family]->() WHERE r.tcdb_gate IN $gp } AS passes
        RETURN sum(CASE WHEN passes AND NOT 'tcdb' IN coalesce(g.annotation_types, []) THEN 1 ELSE 0 END) AS under,
               sum(CASE WHEN NOT passes AND 'tcdb' IN coalesce(g.annotation_types, []) THEN 1 ELSE 0 END) AS over,
               sum(CASE WHEN passes <> ('transporter' IN coalesce(g.informative_annotation_types, []))
                        THEN 1 ELSE 0 END) AS informative_mismatch,
               sum(CASE WHEN passes THEN 1 ELSE 0 END) AS annotated""", gp=GATE_PASSING)[0]
    assert rows["under"] == 0 and rows["over"] == 0 and rows["informative_mismatch"] == 0, rows
    # 14,745 -> ~15,800 with the best-hit arm (spec §7).
    assert 14500 <= rows["annotated"] <= 17500, rows


@pytest.mark.kg
def test_best_hit_arm_adds_annotated_genes(run_query):
    """The best-hit arm (rank 1, tier <= 2) reaches genes the old gate missed
    (~1,055 genes with no eggNOG edge and no diamond_tier edge with tier <= 2)."""
    n = run_query("""
        MATCH (g:Gene) WHERE 'tcdb' IN g.annotation_types
          AND NOT EXISTS { (g)-[r:Gene_has_tcdb_family]->()
                           WHERE r.tcdb_call = 'eggnog'
                              OR (r.tcdb_call = 'diamond_tier' AND r.tier <= 2) }
        RETURN count(g) AS n""")[0]["n"]
    assert 800 <= n <= 1400, f"{n} genes annotated through the best-hit arm only (expected ~1,055)"


@pytest.mark.kg
def test_transport_arm_nests_in_annotation(run_query):
    n = run_query("""
        MATCH (g:Gene) WHERE g.transported_metabolite_count > 0
          AND NOT 'tcdb' IN coalesce(g.annotation_types, [])
        RETURN count(g) AS n""")[0]["n"]
    assert n == 0, f"{n} genes have a transport arm without 'tcdb' in annotation_types"


@pytest.mark.kg
def test_transported_metabolite_count_is_direct_systems_only(run_query):
    """Recomputed from the direct rule (not the stored tcdb_gate): distinct
    Metabolites of the tc_system nodes the gene points at directly."""
    row = run_query(f"""
        MATCH (g:Gene) WHERE EXISTS {{ (g)-[:Gene_has_tcdb_family]->() }}
        OPTIONAL MATCH (g)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
          WHERE {_IS_DIRECT}
        WITH g, count(DISTINCT m) AS expected
        RETURN sum(CASE WHEN coalesce(g.transported_metabolite_count, -1) <> expected THEN 1 ELSE 0 END) AS bad,
               sum(CASE WHEN expected > 0 THEN 1 ELSE 0 END) AS with_arm""")[0]
    assert row["bad"] == 0, f"{row['bad']} genes whose transported_metabolite_count is not the direct-system count"
    # Measured 6,995 on the rebuilt graph (2026-10-23; was ~30K before the cleanup).
    # Fewer than the 8,312 genes with a direct edge, because some direct systems carry
    # no curated substrate in TCDB. Band ~±15%.
    assert 5950 <= row["with_arm"] <= 8050, row


@pytest.mark.kg
def test_gene_catalyzed_metabolite_count_is_catalysis_only(run_query):
    """catalyzed_metabolite_count must match the catalysis arm exactly, with no
    transport contribution folded in."""
    n = run_query("""
        MATCH (g:Gene)
        OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m:Metabolite)
        WITH g, count(DISTINCT m) AS expected
        WHERE coalesce(g.catalyzed_metabolite_count, 0) <> expected
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"{n} genes whose catalyzed_metabolite_count disagrees with the catalysis arm"


@pytest.mark.kg
def test_gene_and_metabolite_transport_projections_agree(run_query):
    """Gene.transported_metabolite_count and Metabolite.transporter_gene_count are
    two projections of ONE direct (gene, metabolite) set, so their totals match."""
    from_genes = run_query(
        "MATCH (g:Gene) RETURN sum(coalesce(g.transported_metabolite_count, 0)) AS n")[0]["n"]
    from_metabolites = run_query(
        "MATCH (m:Metabolite) RETURN sum(coalesce(m.transporter_gene_count, 0)) AS n")[0]["n"]
    assert from_genes == from_metabolites > 0, (from_genes, from_metabolites)


@pytest.mark.kg
def test_metabolite_transporter_counts_are_direct_anchored(run_query):
    """transporter_count = distinct systems carrying m that are the target of >= 1
    direct edge; transporter_gene_count = distinct genes with such an edge.
    Range-only systems never count."""
    row = run_query(f"""
        MATCH (m:Metabolite)
        OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(t:TcdbFamily)<-[r:Gene_has_tcdb_family]-(g:Gene)
          WHERE {_IS_DIRECT}
        WITH m, count(DISTINCT t) AS sys, count(DISTINCT g) AS genes
        RETURN sum(CASE WHEN coalesce(m.transporter_count, -1) <> sys THEN 1 ELSE 0 END) AS bad_sys,
               sum(CASE WHEN coalesce(m.transporter_gene_count, -1) <> genes THEN 1 ELSE 0 END) AS bad_genes,
               sum(CASE WHEN sys > 0 THEN 1 ELSE 0 END) AS reached""")[0]
    assert row["bad_sys"] == 0 and row["bad_genes"] == 0, row
    assert row["reached"] > 0


@pytest.mark.kg
def test_transport_tag_does_not_imply_gene_reach(run_query):
    """'transport' IN evidence_sources = TCDB lists it as a substrate of a kept
    system, not "a gene here reaches it" — range-only systems make both states exist."""
    row = run_query("""
        MATCH (m:Metabolite) WHERE 'transport' IN m.evidence_sources
        RETURN sum(CASE WHEN m.transporter_gene_count = 0 THEN 1 ELSE 0 END) AS unreached,
               sum(CASE WHEN m.transporter_gene_count > 0 THEN 1 ELSE 0 END) AS reached""")[0]
    assert row["unreached"] > 0 and row["reached"] > 0, row


@pytest.mark.kg
def test_organism_metabolite_counts_are_split_by_evidence_arm(run_query):
    n = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (o)-[r:Organism_has_metabolite]->(m:Metabolite)
        WITH o,
             count(DISTINCT CASE WHEN 'metabolism' IN r.evidence_sources THEN m END) AS cat,
             count(DISTINCT CASE WHEN 'transport'  IN r.evidence_sources THEN m END) AS tr
        WHERE coalesce(o.catalyzed_metabolite_count, 0) <> cat
           OR coalesce(o.transported_metabolite_count, 0) <> tr
        RETURN count(o) AS n
    """)[0]["n"]
    assert n == 0, f"{n} organisms whose metabolite counts disagree with their edge evidence_sources"


@pytest.mark.kg
def test_organism_transport_arm_uses_direct_edges(run_query):
    """The transport arm of Organism_has_metabolite = the direct (gene, system) set."""
    row = run_query(f"""
        MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
              -[:Tcdb_family_transports_metabolite]->(m:Metabolite)
        WHERE {_IS_DIRECT}
        WITH DISTINCT o, m
        OPTIONAL MATCH (o)-[e:Organism_has_metabolite]->(m)
        RETURN count(*) AS pairs,
               sum(CASE WHEN e IS NULL OR NOT 'transport' IN e.evidence_sources THEN 1 ELSE 0 END) AS missing""")[0]
    assert row["pairs"] > 0 and row["missing"] == 0, row
    total = run_query("""
        MATCH ()-[e:Organism_has_metabolite]->() WHERE 'transport' IN e.evidence_sources
        RETURN count(e) AS n""")[0]["n"]
    assert total == row["pairs"], f"{total} transport-tagged organism edges vs {row['pairs']} direct pairs"


# ============================================================================
# §5.1: TcdbFamily computed properties
# ============================================================================


@pytest.mark.kg
def test_tcdb_substrate_count_dense_on_every_level(run_query):
    rows = run_query("""
        MATCH (t:TcdbFamily)
        RETURN t.level_kind AS lk, count(t) AS n, count(t.tcdb_substrate_count) AS set_,
               sum(CASE WHEN t.tcdb_substrate_count > 0 THEN 1 ELSE 0 END) AS pos""")
    assert {r["lk"] for r in rows} == set(LEVEL_KINDS), rows
    for r in rows:
        assert r["set_"] == r["n"], f"{r['lk']}: tcdb_substrate_count on {r['set_']}/{r['n']}"
        assert r["pos"] > 0, f"{r['lk']}: no node with tcdb_substrate_count > 0"


@pytest.mark.kg
def test_tcdb_substrate_count_equals_walk(run_query):
    """On EVERY TcdbFamily node (~9.4K): distinct Metabolite + ChemicalRole targets of
    the tc_system nodes in the node's *0.. subtree (a system: its own)."""
    rows = run_query("""
        MATCH (t:TcdbFamily)
        OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(s:TcdbFamily {level_kind: 'tc_system'})
                       -[:Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role]->(x)
        WITH t, count(DISTINCT x) AS walk
        RETURN t.id AS id, t.level_kind AS lk, t.tcdb_substrate_count AS prop, walk""")
    assert {r["lk"] for r in rows} == set(LEVEL_KINDS)
    bad = [r for r in rows if r["prop"] != r["walk"]]
    assert not bad, bad[:5]


@pytest.mark.kg
def test_tcdb_gene_count_uses_gate(run_query):
    """gene_count (subtree) / direct_gene_count (self) / organism_count count only
    gate-passing edges, recomputed from the rule on EVERY TcdbFamily node."""
    rows = run_query(f"""
        MATCH (n:TcdbFamily)
        OPTIONAL MATCH (n)<-[:Tcdb_family_is_a_tcdb_family*0..]-(t:TcdbFamily)<-[r:Gene_has_tcdb_family]-(g:Gene)
          WHERE {_PASSES}
        WITH n, count(DISTINCT g) AS gc, count(DISTINCT CASE WHEN t = n THEN g END) AS dgc,
             count(DISTINCT g.organism_name) AS oc
        RETURN n.id AS id, n.gene_count AS gene_count, gc, n.direct_gene_count AS direct_gene_count, dgc,
               n.organism_count AS organism_count, oc""")
    bad = [r for r in rows if (r["gene_count"], r["direct_gene_count"], r["organism_count"])
           != (r["gc"], r["dgc"], r["oc"])]
    assert not bad, bad[:5]
    # A node reached only by below_gate edges reads 0.
    n = run_query("""
        MATCH (t:TcdbFamily)<-[:Gene_has_tcdb_family {tcdb_gate: 'below_gate'}]-()
        WHERE NOT EXISTS { (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(:TcdbFamily)
                           <-[r:Gene_has_tcdb_family]-() WHERE r.tcdb_gate IN $gp }
        WITH DISTINCT t
        RETURN count(t) AS n, sum(CASE WHEN t.gene_count <> 0 THEN 1 ELSE 0 END) AS bad""", gp=GATE_PASSING)[0]
    assert n["n"] > 0 and n["bad"] == 0, n


@pytest.mark.kg
def test_tcdb_family_has_member_count_and_class_id(run_query):
    n = run_query("""
        MATCH (t:TcdbFamily)
        OPTIONAL MATCH (c:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(t)
        WITH t, count(c) AS mc
        WHERE t.member_count IS NULL OR t.member_count <> mc OR t.tc_class_id IS NULL
        RETURN count(t) AS n""")[0]["n"]
    assert n == 0, f"{n} TcdbFamily nodes with a wrong member_count or no tc_class_id"


@pytest.mark.kg
def test_tcdb_assignment_count_rule(run_query):
    """Gene.tcdb_assignment_count (spec §5.3): the gate-passing targets, minus any target
    that is an ancestor of another target, counted distinct. Recomputed with TCID
    prefixes (independent of the hierarchy edges): a diamond_tier edge on 3.A.1.4 +
    a rank-1 best hit on 3.A.1.4.4 = 1; + eggNOG on 3.A.1.4.5 = 2; eggNOG and a best
    hit on one system = 1; near-ties and tier-3 calls never count."""
    row = run_query("""
        MATCH (g:Gene)
        OPTIONAL MATCH (g)-[r:Gene_has_tcdb_family]->(t:TcdbFamily) WHERE r.tcdb_gate IN $gp
        WITH g, collect(DISTINCT t.tcdb_id) AS ts
        WITH g, size(ts) AS targets,
             size([x IN ts WHERE NOT any(y IN ts WHERE y STARTS WITH x + '.')]) AS expected
        RETURN sum(CASE WHEN coalesce(g.tcdb_assignment_count, -1) <> expected THEN 1 ELSE 0 END) AS bad,
               sum(CASE WHEN targets > expected THEN 1 ELSE 0 END) AS collapsed,
               sum(CASE WHEN expected > 1 THEN 1 ELSE 0 END) AS multi""", gp=GATE_PASSING)[0]
    assert row["bad"] == 0, f"{row['bad']} genes whose tcdb_assignment_count breaks the §5.3 rule"
    assert row["collapsed"] > 0, "no gene where an ancestor target was dropped — rule not exercised"
    assert row["multi"] > 0, "no gene with two sibling assignments — rule not exercised"


@pytest.mark.kg
def test_tcdb_assignment_count_ignores_below_gate_only_genes(run_query):
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
        WITH g, collect(r.tcdb_gate) AS gates WHERE all(x IN gates WHERE x = 'below_gate')
        RETURN count(g) AS n, sum(CASE WHEN coalesce(g.tcdb_assignment_count, -1) <> 0 THEN 1 ELSE 0 END) AS bad""")[0]
    assert n["n"] > 0 and n["bad"] == 0, n


@pytest.mark.kg
def test_gene_tcdb_family_count_renamed(run_query):
    """Gene.tcdb_family_count was renamed tcdb_assignment_count (same meaning)."""
    assert not run_query("MATCH (g:Gene) WHERE g.tcdb_family_count IS NOT NULL "
                         "RETURN g.locus_tag LIMIT 3"), "old name remains"


# ============================================================================
# Direction + names
# ============================================================================


@pytest.mark.kg
def test_direction_dense_and_basis_sparse(run_query):
    """direction dense on every level; direction_basis present iff direction is
    stated; class/subclass names are mechanism classes, so those levels are unknown.
    Measured uptake 716 / export 696 / unknown 7,958 (conservative rule + overrides, 2026-10-24)."""
    rows = run_query("""
        MATCH (t:TcdbFamily)
        RETURN t.direction AS d, t.direction_basis AS b, t.level_kind AS lk, count(*) AS n""")
    by_dir: dict = {}
    for r in rows:
        assert r["d"] in {"uptake", "export", "unknown"}, r
        assert (r["b"] is None) == (r["d"] == "unknown"), r
        assert r["b"] in {None, "abc_numbering", "tcdb_name"}, r
        if r["lk"] in ("tc_class", "tc_subclass"):
            assert r["d"] == "unknown", r
        by_dir[r["d"]] = by_dir.get(r["d"], 0) + r["n"]
    bands = {"uptake": (600, 850), "export": (580, 820), "unknown": (7400, 8700)}
    for d, (lo, hi) in bands.items():
        assert lo <= by_dir.get(d, 0) <= hi, (d, by_dir)


@pytest.mark.kg
def test_direction_inherited_downward(run_query):
    """A node's direction comes from itself or an ancestor, so a stated parent
    direction is never contradicted by a child whose own name states nothing —
    weakest safe check: a child of an ABC-numbered node either inherits it or is
    unknown (its name says export, so numbering and name conflict)."""
    parents = run_query("""
        MATCH (:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(p:TcdbFamily {direction_basis: 'abc_numbering'})
        RETURN count(DISTINCT p) AS n""")[0]["n"]
    assert parents >= 1, "no abc_numbering node with children — check not exercised"
    n = run_query("""
        MATCH (c:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(p:TcdbFamily {direction_basis: 'abc_numbering'})
        WHERE c.direction <> 'unknown'
          AND (c.direction <> p.direction OR c.direction_basis <> 'abc_numbering')
        RETURN count(c) AS n""")[0]["n"]
    assert n == 0, f"{n} children of ABC-numbered nodes contradict the direction"


@pytest.mark.kg
def test_systems_are_named(run_query):
    row = run_query("""
        MATCH (t:TcdbFamily)
        RETURN sum(CASE WHEN t.name IS NULL OR t.name = '' THEN 1 ELSE 0 END) AS missing,
               sum(CASE WHEN t.level_kind = 'tc_system' AND t.name = t.tcdb_id THEN 1 ELSE 0 END) AS unnamed,
               sum(CASE WHEN t.level_kind = 'tc_system' THEN 1 ELSE 0 END) AS systems""")[0]
    assert row["missing"] == 0
    ceiling = _unnamed_systems_in_caches()
    assert row["unnamed"] <= ceiling, (
        f"{row['unnamed']} kept systems named by their TCID (caches leave {ceiling} unnamed)")
    assert row["systems"] > 5000, row


# ============================================================================
# MED4 spot checks (spec §9)
# ============================================================================


@pytest.mark.kg
def test_med4_spot_checks(run_query):
    # urtA: rank-1 best hit to the urea ABC system, which carries urea itself.
    rows = run_query("""
        MATCH (g:Gene {locus_tag: 'PMM0970', organism_name: 'Prochlorococcus MED4'})
              -[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->(t:TcdbFamily)
        WHERE r.rank = 1
        RETURN t.id AS tc, [(t)-[:Tcdb_family_transports_metabolite]->(m) | m.id] AS subs""")
    assert [r["tc"] for r in rows] == ["tcdb:3.A.1.4.4"], rows
    assert "kegg.compound:C00086" in rows[0]["subs"], rows   # urea

    # pstS best hit.
    rows = run_query("""
        MATCH (g:Gene {locus_tag: 'PMM0710', organism_name: 'Prochlorococcus MED4'})
              -[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->(t:TcdbFamily)
        WHERE r.rank = 1 RETURN t.id AS tc""")
    assert [r["tc"] for r in rows] == ["tcdb:3.A.1.7.6"], rows

    # ABC superfamily: no own substrate edges, direction unset at its own level.
    row = run_query("""
        MATCH (t:TcdbFamily {id: 'tcdb:3.A.1'})
        RETURN t.direction AS d, t.direction_basis AS b,
               COUNT { (t)-[:Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role]->() } AS own,
               t.tcdb_substrate_count AS subs""")[0]
    assert row["own"] == 0 and row["d"] == "unknown" and row["b"] is None, row
    assert row["subs"] > 100, row   # the walk still reaches TCDB's full ABC range

    row = run_query("MATCH (t:TcdbFamily {id: 'tcdb:3.A.1.4'}) "
                    "RETURN t.direction AS d, t.direction_basis AS b")[0]
    assert (row["d"], row["b"]) == ("uptake", "abc_numbering"), row


# ============================================================================
# Removed properties (spec §6)
# ============================================================================


@pytest.mark.kg
def test_removed_tcdb_properties_absent(run_query):
    node_props = ["substrate_breadth", "metabolite_count", "is_promiscuous"]
    for p in node_props:
        n = run_query(f"MATCH (t:TcdbFamily) WHERE t.`{p}` IS NOT NULL RETURN count(t) AS n")[0]["n"]
        assert n == 0, f"TcdbFamily.{p} still on {n} nodes"
    edge_props = ["attachment_depth", "evidence_score", "source_agreement", "pfam_support",
                  "go_support", "confidence_score", "consensus_n", "best_hit_tcids",
                  "best_hit_names", "best_hit_identity", "best_hit_substrate_tcids",
                  "best_hit_substrate_ids", "best_hit_substrate_names"]
    for p in edge_props:
        n = run_query(f"MATCH ()-[r:Gene_has_tcdb_family]->() WHERE r.`{p}` IS NOT NULL "
                      f"RETURN count(r) AS n")[0]["n"]
        assert n == 0, f"Gene_has_tcdb_family.{p} still on {n} edges"
    for p in ["transport_substrate_resolution", "tcdb_evidence_score_max"]:
        n = run_query(f"MATCH (g:Gene) WHERE g.`{p}` IS NOT NULL RETURN count(g) AS n")[0]["n"]
        assert n == 0, f"Gene.{p} still on {n} genes"
    n = run_query("MATCH (t:TcdbFamily {level_kind: 'tc_specificity'}) RETURN count(t) AS n")[0]["n"]
    assert n == 0, "level_kind 'tc_specificity' should be renamed 'tc_system'"


# ============================================================================
# Structure (TCDB + CAZy)
# ============================================================================


@pytest.mark.kg
def test_tcdb_family_levels_are_valid(run_query):
    rows = run_query(
        "MATCH (t:TcdbFamily) RETURN DISTINCT t.level_kind AS lk, t.level AS lv"
    )
    expected = {
        ("tc_class", 0), ("tc_subclass", 1), ("tc_family", 2),
        ("tc_subfamily", 3), ("tc_system", 4),
    }
    actual = {(r["lk"], r["lv"]) for r in rows}
    assert actual == expected, f"level_kind/level pairs: {actual ^ expected}"


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
    """914 in the 2026-10-23 kegg_data.json (65 transport-only nodes were removed)."""
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
