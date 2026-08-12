"""KG validity tests for TCDB and CAZy ontologies."""
from __future__ import annotations

import pytest


@pytest.mark.kg
def test_tcdb_family_node_count_in_range(run_query):
    """Pruned to gene-annotated TCDB IDs + their ANCESTORS only.

    The downward arm was removed 2026-08-06 (it left 94.5% of nodes gene-less);
    diamond then joined eggNOG as a second source, taking seeds 540 -> 1,397 and
    kept IDs 704 -> ~1,515. Bounds stay generous to absorb further growth.
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
    """Band widened 2026-08-06: diamond joined eggNOG as a second evidence source,
    taking the edge count from ~16.8K to ~53.8K."""
    n = run_query("MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n")[0]["n"]
    assert 30000 <= n <= 100000, f"Gene_has_tcdb_family count {n} outside 30000-100000"


@pytest.mark.kg
def test_gene_has_tcdb_family_sources_populated(run_query):
    """Every edge carries source provenance from {eggnog, diamond}."""
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        RETURN count(r) AS total,
               sum(CASE WHEN r.sources IS NULL OR size(r.sources) = 0 THEN 1 ELSE 0 END) AS no_source,
               sum(CASE WHEN any(s IN r.sources WHERE NOT s IN ['eggnog','diamond']) THEN 1 ELSE 0 END) AS bad_source
    """)[0]
    assert row["no_source"] == 0, f"{row['no_source']} edges without source provenance"
    assert row["bad_source"] == 0, f"{row['bad_source']} edges with an unknown source label"


@pytest.mark.kg
def test_diamond_evidence_is_sparse_and_well_formed(run_query):
    """Diamond evidence rides ONLY on edges diamond called, and every such edge
    has a tier in 1-3. eggNOG-only edges must carry none of it."""
    egg = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->() WHERE r.sources = ['eggnog']
        RETURN sum(CASE WHEN r.tier IS NOT NULL THEN 1 ELSE 0 END) AS leaked_tier,
               sum(CASE WHEN r.confidence_score IS NOT NULL THEN 1 ELSE 0 END) AS leaked_conf
    """)[0]
    assert egg["leaked_tier"] == 0 and egg["leaked_conf"] == 0

    dia = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->() WHERE 'diamond' IN r.sources
        RETURN count(r) AS n, sum(CASE WHEN r.tier IS NULL THEN 1 ELSE 0 END) AS missing,
               min(r.tier) AS lo, max(r.tier) AS hi
    """)[0]
    assert dia["n"] > 0
    assert dia["missing"] == 0, f"{dia['missing']} diamond edges without a tier"
    assert dia["lo"] >= 1 and dia["hi"] <= 3


@pytest.mark.kg
def test_one_edge_per_gene_and_tc(run_query):
    """A TC called by BOTH sources is one edge with two `sources` entries — never
    two parallel edges. That collapse is the corroboration signal."""
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
        WITH g, t, count(r) AS c WHERE c > 1
        RETURN count(*) AS n
    """)[0]["n"]
    assert n == 0, f"{n} (gene, TC) pairs have parallel edges"


@pytest.mark.kg
def test_annotation_types_tcdb_respects_the_tier_gate(run_query):
    """'tcdb' counts toward annotation_types only for eggNOG-sourced or tier<=2
    edges, so ~22K conservative tier-3 calls stay findable without inflating
    annotation_quality. Asserted in BOTH directions."""
    over = run_query("""
        MATCH (g:Gene) WHERE 'tcdb' IN g.annotation_types
          AND NOT EXISTS { MATCH (g)-[r:Gene_has_tcdb_family]->()
                           WHERE 'eggnog' IN r.sources OR r.tier <= 2 }
        RETURN count(g) AS n
    """)[0]["n"]
    assert over == 0, f"{over} genes counted as tcdb-annotated on tier-3-only evidence"

    under = run_query("""
        MATCH (g:Gene) WHERE EXISTS { MATCH (g)-[r:Gene_has_tcdb_family]->()
                                      WHERE 'eggnog' IN r.sources OR r.tier <= 2 }
          AND NOT 'tcdb' IN g.annotation_types
        RETURN count(g) AS n
    """)[0]["n"]
    assert under == 0, f"{under} genes with qualifying evidence missing from annotation_types"


@pytest.mark.kg
def test_tcdb_family_count_is_not_tier_gated(run_query):
    """Routing counts cover ALL edges — only the quality buckets are gated."""
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
        WITH g, count(r) AS actual WHERE g.tcdb_family_count <> actual
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"{n} genes whose tcdb_family_count disagrees with their edge count"


@pytest.mark.kg
def test_gene_has_cazy_family_edge_count(run_query):
    """Upper bound raised from 2000 by the InterPro two-layer integration
    (2026-08-10): Layer B adds ~642 CAZy ids from FAMILY/DOMAIN entry xrefs."""
    n = run_query("MATCH ()-[r:Gene_has_cazy_family]->() RETURN count(r) AS n")[0]["n"]
    assert 100 <= n <= 5000, f"Gene_has_cazy_family count {n} outside 100-5000"


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


# ============================================================================
# TCDB evidence score (post-import, advisory)
# ============================================================================


@pytest.mark.kg
def test_tcdb_evidence_score_populated_and_in_range(run_query):
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        RETURN count(r) AS total,
               sum(CASE WHEN r.tcdb_evidence_score IS NULL THEN 1 ELSE 0 END) AS missing,
               min(r.tcdb_evidence_score) AS lo, max(r.tcdb_evidence_score) AS hi
    """)[0]
    assert row["missing"] == 0, f"{row['missing']} edges without an evidence score"
    assert row["lo"] >= 0 and row["hi"] <= 5


@pytest.mark.kg
def test_tcdb_evidence_score_equals_sum_of_its_components(run_query):
    """The score must be reconstructable from the stored components + sources/tier.

    This is the guard that keeps it honest: if the total ever drifts from its
    parts, the number has become an opaque verdict — the failure mode that got
    `filter_action` deleted.
    """
    n = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WITH r, (CASE WHEN 'eggnog' IN r.sources THEN 1 ELSE 0 END
               + CASE WHEN r.agrees_across_sources THEN 1 ELSE 0 END
               + CASE WHEN coalesce(r.tier <= 2, false) THEN 1 ELSE 0 END
               + CASE WHEN r.pfam_corroborated THEN 1 ELSE 0 END
               + CASE WHEN r.go_corroborated THEN 1 ELSE 0 END) AS recomputed
        WHERE r.tcdb_evidence_score <> recomputed
        RETURN count(*) AS n
    """)[0]["n"]
    assert n == 0, f"{n} edges whose score disagrees with its components"


@pytest.mark.kg
def test_tcdb_evidence_components_are_booleans_not_null(run_query):
    n = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WHERE r.agrees_across_sources IS NULL
           OR r.pfam_corroborated IS NULL
           OR r.go_corroborated IS NULL
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0, f"{n} edges with a null evidence component"


@pytest.mark.kg
def test_exact_two_source_edges_always_agree(run_query):
    """size(sources)=2 is the strictest form of agreement, so those edges must
    always have agrees_across_sources = true."""
    n = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WHERE size(r.sources) = 2 AND NOT r.agrees_across_sources
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0, f"{n} two-source edges not flagged as agreeing"


@pytest.mark.kg
def test_hierarchical_agreement_dominates_exact(run_query):
    """Agreement is mostly at DIFFERENT hierarchy depths — eggNOG names
    subfamilies, diamond's tier-3 truncation names the parent family. If this
    ever inverts, the agreement definition has silently narrowed to exact-node
    and the score has lost ~83% of its corroboration signal."""
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WHERE r.agrees_across_sources
        RETURN count(r) AS agreeing,
               sum(CASE WHEN size(r.sources) = 2 THEN 1 ELSE 0 END) AS exact
    """)[0]
    assert row["agreeing"] > 0
    assert row["exact"] < row["agreeing"] / 2, (
        f"only {row['exact']}/{row['agreeing']} agreements are hierarchical — "
        "agreement may have narrowed to exact-node"
    )


@pytest.mark.kg
def test_evidence_score_spans_its_full_range(run_query):
    """Every score level is populated — a degenerate distribution (everything 0,
    or everything 5) would mean the score separates nothing."""
    rows = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        RETURN r.tcdb_evidence_score AS score, count(*) AS n ORDER BY score
    """)
    seen = {r["score"] for r in rows}
    assert seen == {0, 1, 2, 3, 4, 5}, f"score levels present: {sorted(seen)}"
    assert all(r["n"] > 100 for r in rows), "a score level is near-empty"


@pytest.mark.kg
def test_is_promiscuous_is_level_gated(run_query):
    """Only tc_family and deeper can be promiscuous.

    Substrate/member counts scale mechanically with level (the step-6 rollup puts
    every descendant's substrates on each ancestor), so a class-level flag was
    vacuous: "Channels and Pores transports many things" is what a class IS. The
    previous absolute-only rule fired on 5 of 7 tc_class and 7 of 34 tc_subclass.
    """
    n = run_query("""
        MATCH (t:TcdbFamily) WHERE t.is_promiscuous AND coalesce(t.level, 0) < 2
        RETURN count(t) AS n
    """)[0]["n"]
    assert n == 0, f"{n} class/subclass nodes flagged promiscuous"


@pytest.mark.kg
def test_is_promiscuous_matches_its_thresholds(run_query):
    n = run_query("""
        MATCH (t:TcdbFamily)
        WITH t, (coalesce(t.level,0) >= 2 AND coalesce(t.metabolite_count,0) >= 50) AS expected
        WHERE coalesce(t.is_promiscuous, false) <> expected
        RETURN count(t) AS n
    """)[0]["n"]
    assert n == 0, f"{n} nodes whose is_promiscuous disagrees with the thresholds"


@pytest.mark.kg
def test_is_promiscuous_flags_the_canonical_broad_families(run_query):
    """ABC (3.A.1) and MFS (2.A.1) are the textbook multi-substrate transporter
    superfamilies — if the rule stops flagging them it has drifted."""
    rows = run_query("""
        MATCH (t:TcdbFamily) WHERE t.tcdb_id IN ['3.A.1', '2.A.1']
        RETURN t.tcdb_id AS tc, t.is_promiscuous AS p
    """)
    assert {r["tc"] for r in rows} == {"3.A.1", "2.A.1"}
    assert all(r["p"] for r in rows), "a canonical broad family is not flagged"


@pytest.mark.kg
def test_gene_best_evidence_score_is_sparse_and_correct(run_query):
    """Set iff the gene has a TCDB edge, and equal to the max edge score.

    Sparse by design: a gene with NO TCDB evidence must stay distinguishable from
    one whose evidence is weak. Writing 0 for both would collapse "no transporter
    signal" into "a poor one".
    """
    row = run_query("""
        MATCH (g:Gene)
        RETURN count(CASE WHEN g.tcdb_best_evidence_score IS NULL
                           AND EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN 1 END) AS missing,
               count(CASE WHEN g.tcdb_best_evidence_score IS NOT NULL
                           AND NOT EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN 1 END) AS spurious
    """)[0]
    assert row["missing"] == 0, f"{row['missing']} annotated genes without a rollup"
    assert row["spurious"] == 0, f"{row['spurious']} unannotated genes with a rollup"

    mismatch = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
        WITH g, max(r.tcdb_evidence_score) AS best
        WHERE g.tcdb_best_evidence_score <> best
        RETURN count(g) AS n
    """)[0]["n"]
    assert mismatch == 0, f"{mismatch} genes whose rollup <> max(edge score)"


@pytest.mark.kg
def test_is_promiscuous_means_substrate_breadth_only(run_query):
    """Every flagged family must actually have many substrates.

    Guards a regression shipped and reverted on 2026-08-07: a `gene_count >= 500`
    arm flagged large-but-substrate-poor families (e.g. 9.B.34 KPSH, which has
    ZERO substrates) as "promiscuous" — the opposite of what the term means. The
    flag is substrate breadth ONLY; use `t.gene_count` directly for family size.
    """
    n = run_query("""
        MATCH (t:TcdbFamily)
        WHERE t.is_promiscuous AND coalesce(t.metabolite_count, 0) < 50
        RETURN count(t) AS n
    """)[0]["n"]
    assert n == 0, f"{n} families flagged promiscuous without substrate breadth"


# ─────────────────────────────────────────────────────────────────────────────
# Substrate rollup depth + the arm split (TCDB rollup fix)
#
# The step-6 rollup materialises every descendant's substrates onto each
# ancestor. Before this fix that inflation flowed unlabelled into three scalars:
# Metabolite.transporter_count read 0 for 83% of transported metabolites, and
# Gene/Metabolite/Organism metabolite counts unioned a p90=11 catalysis arm with
# a p90=554 transport arm under one name.
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.kg
def test_substrate_edges_carry_a_valid_depth_marker(run_query):
    rows = run_query("""
        MATCH ()-[r:Tcdb_family_transports_metabolite]->()
        RETURN r.substrate_depth AS d, count(*) AS n
    """)
    seen = {r["d"]: r["n"] for r in rows}
    assert set(seen) == {"deepest", "ancestor"}, f"unexpected substrate_depth values: {seen}"
    assert seen["deepest"] > 0 and seen["ancestor"] > 0


@pytest.mark.kg
def test_substrate_depth_agrees_with_the_hierarchy(run_query):
    """'deepest' must mean exactly: no kept child of this node carries the same
    substrate. Guards the adapter's build-time shortcut (it checks DIRECT
    children only, valid because the kept set is ancestor-closed) against the
    graph-side definition."""
    n = run_query("""
        MATCH (t:TcdbFamily)-[r:Tcdb_family_transports_metabolite]->(m:Metabolite)
        WITH t, r, m, EXISTS {
            MATCH (c:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(t)
            WHERE (c)-[:Tcdb_family_transports_metabolite]->(m)
        } AS has_child_with_substrate
        WHERE (r.substrate_depth = 'deepest') <> (NOT has_child_with_substrate)
        RETURN count(*) AS n
    """)[0]["n"]
    assert n == 0, f"{n} substrate edges whose substrate_depth contradicts the hierarchy"


@pytest.mark.kg
def test_every_transported_metabolite_has_a_nonzero_transporter_count(run_query):
    """The regression this fix exists for: with the old tc_specificity-only
    filter, 1,218 of 1,462 transported metabolites read transporter_count = 0."""
    n = run_query("""
        MATCH (m:Metabolite)<-[:Tcdb_family_transports_metabolite]-()
        WITH DISTINCT m WHERE coalesce(m.transporter_count, 0) = 0
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"{n} transported metabolites still have transporter_count = 0"


@pytest.mark.kg
def test_transporter_count_never_double_counts_an_ancestor(run_query):
    """transporter_count counts 'deepest' edges only, so it must never exceed the
    number of distinct TcdbFamily nodes with any substrate edge to the metabolite."""
    n = run_query("""
        MATCH (m:Metabolite)<-[:Tcdb_family_transports_metabolite]-(t:TcdbFamily)
        WITH m, count(DISTINCT t) AS all_levels
        WHERE coalesce(m.transporter_count, 0) > all_levels
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0, f"{n} metabolites whose transporter_count exceeds their substrate-edge sources"


@pytest.mark.kg
def test_gene_metabolite_count_is_catalysis_only(run_query):
    """BREAKING split: metabolite_count must match the catalysis arm exactly, with
    no transport contribution folded in."""
    n = run_query("""
        MATCH (g:Gene)
        OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m:Metabolite)
        WITH g, count(DISTINCT m) AS expected
        WHERE coalesce(g.metabolite_count, 0) <> expected
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"{n} genes whose metabolite_count disagrees with the catalysis arm"


@pytest.mark.kg
def test_transported_metabolite_count_uses_deepest_attachments_only(run_query):
    n = run_query("""
        MATCH (g:Gene)
        OPTIONAL MATCH (g)-[:Gene_has_tcdb_family]->(t:TcdbFamily)
        WHERE NOT EXISTS {
            MATCH (g)-[:Gene_has_tcdb_family]->(d:TcdbFamily)
            WHERE (d)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t)
        }
        OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
        WITH g, count(DISTINCT m) AS expected
        WHERE coalesce(g.transported_metabolite_count, 0) <> expected
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"{n} genes whose transported_metabolite_count is not the deepest-attachment count"


@pytest.mark.kg
def test_transport_substrate_resolution_vocabulary_and_sparsity(run_query):
    rows = run_query("""
        MATCH (g:Gene) WHERE g.transport_substrate_resolution IS NOT NULL
        RETURN g.transport_substrate_resolution AS v, count(*) AS n
    """)
    seen = {r["v"]: r["n"] for r in rows}
    assert set(seen) <= {"resolved", "family_inferred"}, f"unexpected vocabulary: {seen}"
    assert seen.get("resolved", 0) > 0 and seen.get("family_inferred", 0) > 0

    # Sparse by design: set only on genes that actually carry a TCDB edge, so
    # "no transporter evidence" stays distinguishable from a weak substrate claim.
    n = run_query("""
        MATCH (g:Gene)
        WHERE g.transport_substrate_resolution IS NOT NULL
          AND NOT (g)-[:Gene_has_tcdb_family]->()
        RETURN count(g) AS n
    """)[0]["n"]
    assert n == 0, f"{n} genes carry a resolution verdict without any TCDB edge"


@pytest.mark.kg
def test_resolution_is_not_a_restatement_of_tier(run_query):
    """Substrate resolution and evidence tier are orthogonal — a tier gate would
    discard the narrow 2.A.x carriers that are 'resolved' but tier-3-only. If this
    ever hits 0, the flag has silently collapsed into the tier gate."""
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
        WHERE g.transport_substrate_resolution = 'resolved'
        WITH g, collect(coalesce(r.tier, 9)) AS tiers, collect(r.sources) AS srcs
        WHERE none(t IN tiers WHERE t <= 2)
          AND none(s IN srcs WHERE 'eggnog' IN s)
        RETURN count(g) AS n
    """)[0]["n"]
    assert n > 0, "no 'resolved' tier-3-only genes — resolution has collapsed into the tier gate"


@pytest.mark.kg
def test_gene_and_metabolite_transport_projections_agree(run_query):
    """Gene.transported_metabolite_count and Metabolite.transporter_gene_count are
    two projections of ONE (gene, metabolite) set, so their totals must match."""
    from_genes = run_query(
        "MATCH (g:Gene) RETURN sum(coalesce(g.transported_metabolite_count, 0)) AS n"
    )[0]["n"]
    from_metabolites = run_query(
        "MATCH (m:Metabolite) RETURN sum(coalesce(m.transporter_gene_count, 0)) AS n"
    )[0]["n"]
    assert from_genes == from_metabolites, (
        f"gene-side total {from_genes} != metabolite-side total {from_metabolites}")


@pytest.mark.kg
def test_organism_metabolite_counts_are_split_by_evidence_arm(run_query):
    n = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (o)-[r:Organism_has_metabolite]->(m:Metabolite)
        WITH o,
             count(DISTINCT CASE WHEN 'metabolism' IN r.evidence_sources THEN m END) AS cat,
             count(DISTINCT CASE WHEN 'transport'  IN r.evidence_sources THEN m END) AS tr
        WHERE coalesce(o.metabolite_count, 0) <> cat
           OR coalesce(o.transported_metabolite_count, 0) <> tr
        RETURN count(o) AS n
    """)[0]["n"]
    assert n == 0, f"{n} organisms whose metabolite counts disagree with their edge evidence_sources"
