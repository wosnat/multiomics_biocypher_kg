"""
DE table-quality properties on the live graph (spec 2026-11-18-de-table-quality-design.md).

Checks the Experiment caveat / scale stamps, the six per-time-point statistics
computed post-import, and that those statistics agree with the pure-Python
`time_point_stats` recomputed from the live edges.
"""

import pytest

from multiomics_kg.utils.de_calls import time_point_stats

pytestmark = pytest.mark.kg

ISMEJ = "10.1038/ismej.2016.70_coculture_alteromonas_hot1a3_med4_rnaseq"

STAT_ARRAYS = [
    "time_point_negative_fractions",
    "time_point_median_log2fc",
    "time_point_max_abs_log2fc",
    "time_point_large_fc_not_significant",
    "time_point_extreme_fc_count",
    "time_point_low_signal_count",
]


def test_de_experiments_have_dense_caveats_and_scale(run_query):
    rows = run_query("""
        MATCH (e:Experiment)
        WHERE EXISTS { (e)-[:Changes_expression_of]->() }
          AND (e.data_caveats IS NULL OR e.fold_change_provenance IS NULL
               OR NOT e.fold_change_type IN ['log2', 'linear', 'signed_linear'])
        RETURN e.id AS id, e.data_caveats AS caveats, e.fold_change_type AS fct
        LIMIT 10
    """)
    assert rows == []


# time_point_low_signal_count is present only where level columns were declared.
DENSE_STAT_ARRAYS = [a for a in STAT_ARRAYS if a != "time_point_low_signal_count"]


def test_statistics_arrays_parallel(run_query):
    conds = " OR ".join(
        f"e.{a} IS NULL OR size(e.{a}) <> size(e.time_point_orders)"
        for a in DENSE_STAT_ARRAYS)
    rows = run_query(f"""
        MATCH (e:Experiment)
        WHERE EXISTS {{ (e)-[:Changes_expression_of]->() }}
          AND ({conds})
        RETURN e.id AS id LIMIT 10
    """)
    assert rows == []


def _live_vs_stored(run_query, exp_id):
    exp = run_query(
        "MATCH (e:Experiment {id: $id}) RETURN e", id=exp_id)[0]["e"]
    edges = run_query("""
        MATCH (:Experiment {id: $id})-[r:Changes_expression_of]->()
        RETURN r.time_point_order AS tp, r.log2_fold_change AS fc,
               r.expression_status AS st, r.low_signal_reason AS lsr
    """, id=exp_id)
    by_tp = {}
    for r in edges:
        by_tp.setdefault(r["tp"], []).append((r["fc"], r["st"], r["lsr"]))
    orders = list(exp["time_point_orders"])
    assert edges, f"{exp_id}: no Changes_expression_of edges"
    assert orders, f"{exp_id}: empty time_point_orders"
    assert sorted(by_tp) == sorted(orders)
    # arrays are ordered by time_point_order
    for i, order in enumerate(orders):
        s = time_point_stats(by_tp[order])
        assert exp["time_point_negative_fractions"][i] == pytest.approx(s.negative_fraction, abs=1e-9)
        assert exp["time_point_median_log2fc"][i] == pytest.approx(s.median_log2fc, abs=1e-9)
        assert exp["time_point_max_abs_log2fc"][i] == pytest.approx(s.max_abs_log2fc, abs=1e-9)
        assert exp["time_point_large_fc_not_significant"][i] == s.large_fc_not_significant
        assert exp["time_point_extreme_fc_count"][i] == s.extreme_fc_count
        if exp.get("low_signal_basis") is None:
            assert exp.get("time_point_low_signal_count") is None
        else:
            assert exp["time_point_low_signal_count"][i] == s.low_signal_count


def _sample_ids(run_query):
    ids = [ISMEJ]
    multi = run_query("""
        MATCH (e:Experiment) WHERE e.time_point_count > 3
        RETURN e.id AS id ORDER BY e.id LIMIT 1""")
    sig_only = run_query("""
        MATCH (e:Experiment) WHERE e.table_scope = 'significant_only'
        RETURN e.id AS id ORDER BY e.id LIMIT 1""")
    for rows in (multi, sig_only):
        assert rows, "no sample experiment found"
        ids.append(rows[0]["id"])
    return ids


def test_statistics_match_python(run_query):
    for exp_id in _sample_ids(run_query):
        _live_vs_stored(run_query, exp_id)


def test_curated_caveats_live(run_query):
    rows = run_query("""
        MATCH (e:Experiment {id: $id})
        RETURN e.data_caveats AS caveats, e.data_caveats_reviewed_on AS reviewed
    """, id=ISMEJ)
    assert rows, f"{ISMEJ} missing"
    assert "fold_change_inflated_at_low_counts" in rows[0]["caveats"]
    assert rows[0]["reviewed"]


def test_low_signal_count_is_edge_rollup(run_query):
    rows = run_query("""
        MATCH (e:Experiment)
        WHERE EXISTS { (e)-[:Changes_expression_of]->() }
        OPTIONAL MATCH (e)-[r:Changes_expression_of]->()
          WHERE r.low_signal_reason IS NOT NULL
        WITH e, count(r) AS n_edges
        WHERE e.low_signal_basis IS NOT NULL
          AND (e.time_point_low_signal_count IS NULL
               OR size(e.time_point_low_signal_count) <> size(e.time_point_orders)
               OR reduce(s = 0, x IN e.time_point_low_signal_count | s + x) <> n_edges)
        RETURN e.id AS id LIMIT 10
    """)
    assert rows == []


def test_low_signal_count_absent_when_not_assessed(run_query):
    """No low_signal_basis = not assessed: the count is absent, never [0, ...]."""
    rows = run_query("""
        MATCH (e:Experiment)
        WHERE e.low_signal_basis IS NULL AND e.time_point_low_signal_count IS NOT NULL
        RETURN e.id AS id LIMIT 10
    """)
    assert rows == []


def test_low_signal_assessed_somewhere(run_query):
    rows = run_query("""
        MATCH (e:Experiment) WHERE e.low_signal_basis IS NOT NULL
        RETURN count(e) AS n
    """)
    assert rows[0]["n"] > 0


def test_statistics_absent_without_de(run_query):
    """The six statistics are absent (not []) on experiments with no DE edges."""
    conds = " OR ".join(f"e.{a} IS NOT NULL" for a in STAT_ARRAYS)
    rows = run_query(f"""
        MATCH (e:Experiment)
        WHERE NOT EXISTS {{ (e)-[:Changes_expression_of]->() }} AND ({conds})
        RETURN e.id AS id LIMIT 10
    """)
    assert rows == []
