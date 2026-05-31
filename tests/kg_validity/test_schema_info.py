"""
Release-metadata tests for the `Schema_info` node.

`Schema_info` is created by BioCypher (`bc.write_schema_info(as_node=True)`); the
post-import scripts (`scripts/post-import.sh` Group 4 + `scripts/post-import.cypher`)
stamp it with release identity and live counts on every build:

  - Identity: version, built_at, git_sha[_short], git_branch, git_dirty,
    mcp_min_version, release_notes_url
  - Counts (computed, not hardcoded — track data drift):
    paper_count, experiment_count, gene_count, organism_count,
    expression_edge_count

Dev builds stamp `version="0.0.0-dev"` / `git_sha="unknown"` (KG_* env unset →
bash :- defaults → coalesce-guarded params). Release builds stamp the tagged
version + the real git identity.

See `plans/alpha_release.md` §2.1.
"""

import pytest


pytestmark = pytest.mark.kg


# Property → expected python type. release_notes_url MAY be empty (operator-set);
# everything else must be non-empty / >= 0.
STRING_PROPS = [
    "version",
    "built_at",
    "git_sha",
    "git_sha_short",
    "git_branch",
    "git_dirty",
    "mcp_min_version",
    "release_notes_url",
]

INT_PROPS = [
    "paper_count",
    "experiment_count",
    "gene_count",
    "organism_count",
    "expression_edge_count",
]

# (Schema_info property, Cypher COUNT subquery) — must agree.
COUNT_PAIRS = [
    ("paper_count",           "COUNT { (:Publication) }"),
    ("experiment_count",      "COUNT { (:Experiment) }"),
    ("gene_count",            "COUNT { (:Gene) }"),
    ("organism_count",        "COUNT { (:OrganismTaxon) }"),
    ("expression_edge_count", "COUNT { ()-[:Changes_expression_of]->() }"),
]


@pytest.fixture(scope="module")
def schema_info(run_query):
    rows = run_query(
        "MATCH (s:Schema_info {id: 'schema_info'}) RETURN properties(s) AS props"
    )
    assert len(rows) == 1, (
        f"Expected exactly 1 Schema_info node with id='schema_info'; got {len(rows)}. "
        "BioCypher's write_schema_info(as_node=True) should produce it during build."
    )
    return rows[0]["props"]


# ---------------------------------------------------------------------------
# Release-identity properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("key", STRING_PROPS)
def test_release_identity_property_present_and_string(key, schema_info):
    """Every release-identity property must be present and string-typed."""
    assert key in schema_info, (
        f"Schema_info is missing release property '{key}'. "
        "Post-import Group 4 did not run, or did not stamp this property."
    )
    assert isinstance(schema_info[key], str), (
        f"Schema_info.{key} = {schema_info[key]!r} (type {type(schema_info[key]).__name__}); "
        "expected str."
    )


@pytest.mark.parametrize("key", [k for k in STRING_PROPS if k != "release_notes_url"])
def test_release_identity_property_nonempty(key, schema_info):
    """
    Non-URL identity props must be non-empty.

    Specifically guards against the empty-string regression where
    `${KG_RELEASE_VERSION:-}` + `coalesce($version, '0.0.0-dev')` would have
    left version=''. The fix is bash `${VAR:-default}` (the `:-` form fires on
    unset OR empty), so dev builds reach this with the default string.
    """
    val = schema_info[key]
    assert val != "", (
        f"Schema_info.{key} is empty. "
        "Check the bash default in post-import.sh Group 4 (${{VAR:-default}}) — "
        "coalesce() alone does NOT catch empty strings."
    )


def test_version_format(schema_info):
    """version is either the dev sentinel or a release-scheme string."""
    import re
    version = schema_info["version"]
    pattern = r"^(?:0\.0\.0-dev|\d+\.\d+\.\d+(?:-(?:alpha|beta|rc)\.\d+)?)$"
    assert re.match(pattern, version), (
        f"Schema_info.version = {version!r}; expected '0.0.0-dev' (dev build) "
        f"or 'X.Y.Z[-(alpha|beta|rc).N]' (release build)."
    )


def test_built_at_is_iso8601(schema_info):
    """built_at must be an ISO 8601 string (toString(datetime()))."""
    from datetime import datetime
    built_at = schema_info["built_at"]
    # Neo4j toString(datetime()) emits e.g. "2026-05-25T11:52:44.149Z"
    try:
        datetime.fromisoformat(built_at.replace("Z", "+00:00"))
    except ValueError as exc:
        pytest.fail(f"Schema_info.built_at = {built_at!r} is not ISO 8601: {exc}")


def test_git_dirty_is_boolean_string(schema_info):
    """git_dirty follows the project convention of bool-as-string."""
    val = schema_info["git_dirty"]
    assert val in {"true", "false", "unknown"}, (
        f"Schema_info.git_dirty = {val!r}; expected 'true' | 'false' | 'unknown' "
        "(project convention: booleans serialized as strings)."
    )


# ---------------------------------------------------------------------------
# Count properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("key", INT_PROPS)
def test_count_property_present_and_int(key, schema_info):
    assert key in schema_info, f"Schema_info is missing count property '{key}'."
    val = schema_info[key]
    assert isinstance(val, int), (
        f"Schema_info.{key} = {val!r} (type {type(val).__name__}); expected int."
    )
    assert val >= 0, f"Schema_info.{key} = {val} (negative)."


@pytest.mark.parametrize("prop, count_expr", COUNT_PAIRS)
def test_count_matches_live_graph(prop, count_expr, schema_info, run_query):
    """
    Counts on Schema_info must match the live `COUNT { ... }` of the graph.

    Catches stale stats: post-import block didn't run, or data drifted after
    stamping (e.g. a new adapter wrote nodes after Group 4).
    """
    live = run_query(f"RETURN {count_expr} AS cnt")[0]["cnt"]
    stamped = schema_info[prop]
    assert stamped == live, (
        f"Schema_info.{prop} = {stamped} but {count_expr} = {live}. "
        "Re-run post-process so the stamp matches the data."
    )
