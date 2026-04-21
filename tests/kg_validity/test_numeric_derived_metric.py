"""
KG validity tests for numeric DerivedMetric post-import computations:
rank_by_metric, metric_percentile, metric_bucket, significant.

These assertions require numeric DMs in the graph. The synthetic paperconfig
fixture at tests/fixtures/non_de/ provides deterministic expected values:
100 rows × 3 numeric metrics (fourier_score, peak_time_h, peak_fit_r_squared).

Fixture metric config (tests/fixtures/non_de/synthetic_paperconfig.yaml):
  fourier_score:      rankable=true,  has_p_value=true,  p_value_threshold=0.05
  peak_time_h:        rankable=false, has_p_value=false
  peak_fit_r_squared: rankable=true,  has_p_value=false

Expected distributions (100 rows each):
  Rankable metrics: rank 1..100 contiguous; percentile 100.0 → 0.0
  Significance (fourier_score only): all 100 edges have `significant` set
    (values depend on per-row adjusted_p_value vs 0.05 threshold)
"""

import pytest

pytestmark = pytest.mark.kg


FIXTURE_DOI = "10.9999/synthetic-numeric-dm"


def _fixture_loaded(run_query):
    """Return True if the synthetic fixture publication is in the graph."""
    result = run_query(
        "MATCH (p:Publication {doi: $doi}) RETURN count(p) AS cnt",
        doi=FIXTURE_DOI,
    )
    return result[0]["cnt"] > 0


@pytest.fixture(scope="module", autouse=True)
def skip_if_no_fixture(run_query):
    if not _fixture_loaded(run_query):
        pytest.skip(
            f"Synthetic numeric-DM fixture ({FIXTURE_DOI}) not in graph — "
            f"wire tests/fixtures/non_de/paperconfig_files.txt into create_knowledge_graph.py "
            f"and rebuild Docker"
        )


# ---------------------------------------------------------------------------
# Fixture presence sanity
# ---------------------------------------------------------------------------

def test_fixture_dm_nodes_present(run_query):
    result = run_query("""
        MATCH (p:Publication {doi: $doi})-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        RETURN dm.metric_type AS metric_type, dm.rankable AS rankable,
               dm.has_p_value AS has_p_value
        ORDER BY metric_type
    """, doi=FIXTURE_DOI)
    metric_types = {r["metric_type"] for r in result}
    assert metric_types == {"fourier_score", "peak_time_h", "peak_fit_r_squared"}


def test_fixture_quantifies_edge_count(run_query):
    """300 edges total: 100 rows × 3 numeric metrics."""
    result = run_query("""
        MATCH (p:Publication {doi: $doi})-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
          -[r:Derived_metric_quantifies_gene]->()
        RETURN count(r) AS cnt
    """, doi=FIXTURE_DOI)
    assert result[0]["cnt"] == 300


# ---------------------------------------------------------------------------
# Rank contiguity on rankable DMs
# ---------------------------------------------------------------------------

def test_rank_by_metric_only_on_rankable(run_query):
    """rank_by_metric non-null iff parent DM.rankable='true' (quantifies edges only)."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WITH
          count(CASE WHEN r.rank_by_metric IS NOT NULL AND dm.rankable <> 'true'
                     THEN 1 END) AS rank_on_non_rankable,
          count(CASE WHEN r.rank_by_metric IS NULL AND dm.rankable = 'true'
                     THEN 1 END) AS rankable_missing_rank
        RETURN rank_on_non_rankable, rankable_missing_rank
    """)
    row = result[0]
    assert row["rank_on_non_rankable"] == 0, (
        f"{row['rank_on_non_rankable']} edges have rank_by_metric set on non-rankable DM"
    )
    assert row["rankable_missing_rank"] == 0, (
        f"{row['rankable_missing_rank']} rankable-DM edges are missing rank_by_metric"
    )


def test_rank_contiguous_per_dm(run_query):
    """rank_by_metric should be 1..N contiguous per DerivedMetric."""
    result = run_query("""
        MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WITH dm.id AS dm_id, collect(r.rank_by_metric) AS ranks,
             size(collect(r)) AS n
        WITH dm_id, ranks, n, apoc.coll.sort(ranks) AS sorted
        WHERE sorted <> range(1, n)
        RETURN dm_id, n, sorted
    """)
    assert len(result) == 0, f"Non-contiguous ranks: {result}"


# ---------------------------------------------------------------------------
# Percentile and bucket consistency
# ---------------------------------------------------------------------------

def test_percentile_in_range(run_query):
    """metric_percentile ∈ [0, 100] on every quantifies edge where it's set."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_percentile IS NOT NULL
          AND (r.metric_percentile < 0.0 OR r.metric_percentile > 100.0)
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_bucket_matches_pinned_thresholds(run_query):
    """metric_bucket follows pinned thresholds: top_decile>=90, top_quartile>=75<90, mid>=25<75, low<25."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_bucket IS NOT NULL
        WITH r, r.metric_percentile AS pct, r.metric_bucket AS bucket
        WITH r, pct, bucket,
             CASE
               WHEN pct >= 90.0 THEN 'top_decile'
               WHEN pct >= 75.0 THEN 'top_quartile'
               WHEN pct >= 25.0 THEN 'mid'
               ELSE 'low'
             END AS expected
        WHERE bucket <> expected
        RETURN count(r) AS bad, collect([pct, bucket, expected])[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have bucket != pinned-threshold match: {result[0]['examples']}"
    )


def test_rank_1_is_top_decile(run_query):
    """The top-ranked gene per rankable DM is always top_decile."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene {rank_by_metric: 1}]->()
        WHERE r.metric_bucket <> 'top_decile'
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_highest_value_has_rank_1(run_query):
    """For each rankable DM, the edge with highest value has rank_by_metric=1."""
    result = run_query("""
        MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WITH dm, max(r.value) AS max_val
        MATCH (dm)-[r2:Derived_metric_quantifies_gene]->()
        WHERE r2.value = max_val
        WITH dm, r2
        WHERE r2.rank_by_metric <> 1
        RETURN dm.id AS dm_id, r2.rank_by_metric AS rank, r2.value AS value
    """)
    assert len(result) == 0, f"Max-value edge has rank != 1: {result}"


# ---------------------------------------------------------------------------
# Significance gating
# ---------------------------------------------------------------------------

def test_significant_only_on_has_p_value(run_query):
    """significant non-null only when parent DM.has_p_value='true' AND p_value_threshold IS NOT NULL
    AND edge.adjusted_p_value IS NOT NULL."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND (dm.has_p_value <> 'true'
               OR dm.p_value_threshold IS NULL
               OR r.adjusted_p_value IS NULL)
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_significant_computed_correctly(run_query):
    """significant='true' iff adjusted_p_value < threshold."""
    result = run_query("""
        MATCH (dm:DerivedMetric {has_p_value: 'true'})-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND dm.p_value_threshold IS NOT NULL
          AND r.adjusted_p_value IS NOT NULL
        WITH r, dm,
             CASE WHEN r.adjusted_p_value < dm.p_value_threshold THEN 'true' ELSE 'false' END AS expected
        WHERE r.significant <> expected
        RETURN count(r) AS bad, collect([r.adjusted_p_value, dm.p_value_threshold, r.significant, expected])[..3] AS examples
    """)
    assert result[0]["bad"] == 0, f"{result[0]['bad']}: {result[0]['examples']}"


def test_significant_enum(run_query):
    """When non-null, significant ∈ {'true','false'}."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
          AND NOT r.significant IN ['true', 'false']
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_fourier_score_has_significance(run_query):
    """Synthetic fixture fourier_score (has_p_value='true', threshold=0.05) has 100 edges with significant set."""
    result = run_query("""
        MATCH (dm:DerivedMetric {metric_type: 'fourier_score'})
          -[r:Derived_metric_quantifies_gene]->()
        RETURN count(r) AS total, count(r.significant) AS with_sig
    """)
    assert result[0]["total"] == 100
    assert result[0]["with_sig"] == 100, (
        f"Expected all 100 fourier_score edges to have significant set; got {result[0]['with_sig']}"
    )


def test_peak_time_has_no_significance(run_query):
    """Synthetic fixture peak_time_h (has_p_value='false') should have no significant set."""
    result = run_query("""
        MATCH (dm:DerivedMetric {metric_type: 'peak_time_h'})
          -[r:Derived_metric_quantifies_gene]->()
        WHERE r.significant IS NOT NULL
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0
