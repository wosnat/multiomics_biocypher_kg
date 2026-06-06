"""
Unit tests for `.claude/skills/release-kg/release_kg.py` helpers.

Covers the deterministic, network-free helpers:
- VERSION_RE
- cut_changelog (idempotency, missing [Unreleased], empty/with-entries body)
- extract_changelog_fragment
- build_metadata shape
- _parse_plain_row

The integration paths (subprocess to git/docker/gh/cypher-shell) are covered by
`release_kg.py --dry-run`, not unit tests.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


# ─── Module loader (skill lives outside the package tree) ────────────────────
REPO_ROOT = Path(__file__).resolve().parents[1]
RKG_PATH = REPO_ROOT / ".claude" / "skills" / "release-kg" / "release_kg.py"


@pytest.fixture(scope="module")
def rkg():
    spec = importlib.util.spec_from_file_location("release_kg", RKG_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["release_kg"] = mod
    spec.loader.exec_module(mod)
    return mod


# ─── VERSION_RE ──────────────────────────────────────────────────────────────
@pytest.mark.parametrize("v", [
    "0.1.0",
    "1.0.0",
    "10.20.30",
    "0.1.0-alpha.1",
    "0.1.0-beta.12",
    "1.2.3-rc.99",
])
def test_version_regex_accepts(v, rkg):
    assert rkg.VERSION_RE.match(v), f"should accept {v!r}"


@pytest.mark.parametrize("v", [
    "",
    "0.1",
    "0.1.0.0",
    "0.1.0-dev",          # dev is the post-import default; not a valid release tag
    "0.1.0-alpha",        # missing .N
    "0.1.0-alpha.1+sha",  # no build metadata in our scheme
    "v0.1.0",
    "0.1.0 ",
    "0.0.0-pre.1",        # not in {alpha, beta, rc}
])
def test_version_regex_rejects(v, rkg):
    assert not rkg.VERSION_RE.match(v), f"should reject {v!r}"


# ─── cut_changelog ───────────────────────────────────────────────────────────
SAMPLE_WITH_ENTRIES = """# Changelog

## [Unreleased]

### Added
- Item A
- Item B

### Fixed
- Bug X

## [0.0.1] - 2026-01-01

### Added
- Genesis.
"""


SAMPLE_EMPTY_UNRELEASED = """# Changelog

## [Unreleased]

### Added

### Changed

### Fixed

## [0.0.1] - 2026-01-01

- Genesis.
"""


SAMPLE_NO_UNRELEASED = """# Changelog

## [0.0.1] - 2026-01-01

- Genesis.
"""


def test_cut_with_entries_moves_body(rkg):
    new, performed = rkg.cut_changelog(SAMPLE_WITH_ENTRIES, "0.1.0-alpha.1", "2026-05-25")
    assert performed
    # New empty [Unreleased] appears above the version section
    assert "## [Unreleased]" in new
    assert "## [0.1.0-alpha.1] - 2026-05-25" in new
    # The old entries moved into the version section
    assert "- Item A" in new.split("## [0.1.0-alpha.1]")[1]
    assert "- Bug X" in new.split("## [0.1.0-alpha.1]")[1]
    # The old prior version stays intact
    assert "## [0.0.1] - 2026-01-01" in new
    # Order: new Unreleased FIRST, then new version, then old version
    iu = new.index("## [Unreleased]")
    iv = new.index("## [0.1.0-alpha.1]")
    io = new.index("## [0.0.1]")
    assert iu < iv < io


def test_cut_empty_unreleased_gets_placeholder(rkg):
    new, performed = rkg.cut_changelog(SAMPLE_EMPTY_UNRELEASED, "0.1.0-alpha.1", "2026-05-25")
    assert performed
    assert "## [0.1.0-alpha.1] - 2026-05-25" in new
    after_version = new.split("## [0.1.0-alpha.1]")[1].split("## [0.0.1]")[0]
    assert "_No entries — placeholder._" in after_version


def test_cut_idempotent(rkg):
    """Re-cutting the same version is a no-op."""
    new1, performed1 = rkg.cut_changelog(SAMPLE_WITH_ENTRIES, "0.1.0-alpha.1", "2026-05-25")
    new2, performed2 = rkg.cut_changelog(new1, "0.1.0-alpha.1", "2026-05-26")  # different date
    assert performed1 is True
    assert performed2 is False
    assert new1 == new2  # content unchanged on second call


def test_cut_missing_unreleased_raises(rkg):
    with pytest.raises(ValueError, match="\\[Unreleased\\]"):
        rkg.cut_changelog(SAMPLE_NO_UNRELEASED, "0.1.0-alpha.1", "2026-05-25")


def test_cut_preserves_header_lines(rkg):
    new, _ = rkg.cut_changelog(SAMPLE_WITH_ENTRIES, "0.1.0-alpha.1", "2026-05-25")
    assert new.startswith("# Changelog\n")


# ─── extract_changelog_fragment ──────────────────────────────────────────────
def test_extract_fragment(rkg, tmp_path: Path):
    cut, _ = rkg.cut_changelog(SAMPLE_WITH_ENTRIES, "0.1.0-alpha.1", "2026-05-25")
    p = tmp_path / "CHANGELOG.md"
    p.write_text(cut)
    fragment = rkg.extract_changelog_fragment(p, "0.1.0-alpha.1")
    assert fragment.startswith("## [0.1.0-alpha.1] - 2026-05-25")
    assert "- Item A" in fragment
    # Should NOT bleed into the next version's section
    assert "## [0.0.1]" not in fragment
    assert "Genesis." not in fragment


# ─── build_metadata ──────────────────────────────────────────────────────────
def test_build_metadata_shape(rkg):
    ctx = rkg.Context(
        version="0.1.0-alpha.1", target="staging", mcp_min="0.1.0",
        allow_dirty=False, draft=False, dry_run=True, resume=False,
        skip_kg_tests=False,
        git_sha="a" * 40, git_sha_short="aaaaaaa", git_branch="main", git_dirty=False,
        schema_info={
            "version": "0.1.0-alpha.1",
            "papers": 43, "experiments": 197, "genes": 99871,
            "organisms": 37, "expr_edges": 232758,
            "built_at": "2026-05-25T11:52:44.149Z",
        },
    )
    meta = rkg.build_metadata(ctx)
    assert meta["version"] == "0.1.0-alpha.1"
    assert meta["tag"] == "kg-0.1.0-alpha.1"
    assert meta["counts"]["genes"] == 99871
    assert meta["counts"]["expression_edges"] == 232758
    assert meta["built_at"] == "2026-05-25T11:52:44.149Z"
    # stamped_at must be present and ISO-shaped
    assert "T" in meta["stamped_at"]
    # target carried through
    assert meta["target"] == "staging"
    # per_publication_edges is always present (empty dict when no per-pub data)
    assert meta["per_publication_edges"] == {}


def test_build_metadata_includes_per_publication_edges(rkg):
    ctx = rkg.Context(
        version="0.1.0-alpha.1", target="staging", mcp_min="0.1.0",
        allow_dirty=False, draft=False, dry_run=True, resume=False,
        skip_kg_tests=False,
        git_sha="a" * 40, git_sha_short="aaaaaaa", git_branch="main", git_dirty=False,
        schema_info={"version": "0.1.0-alpha.1", "built_at": "2026-05-25T11:52:44.149Z"},
        per_publication_edges={"10.b/zzz": 100, "10.a/aaa": 50},
    )
    meta = rkg.build_metadata(ctx)
    # Sorted on the way out so diffs against future releases are stable
    assert list(meta["per_publication_edges"].keys()) == ["10.a/aaa", "10.b/zzz"]
    assert meta["per_publication_edges"]["10.a/aaa"] == 50
    assert meta["per_publication_edges"]["10.b/zzz"] == 100


# ─── _parse_plain_row ────────────────────────────────────────────────────────
def test_parse_plain_row(rkg):
    sample = (
        "version, git_sha, papers, experiments, genes, organisms, expr_edges, mcp_min, built_at\n"
        '"0.1.0-alpha.1", "abc123", 43, 197, 99871, 37, 232758, "0.1.0", "2026-05-25T11:52:44.149Z"'
    )
    row = rkg._parse_plain_row(sample)
    assert row["version"] == "0.1.0-alpha.1"
    assert row["git_sha"] == "abc123"
    assert row["papers"] == 43
    assert row["experiments"] == 197
    assert row["genes"] == 99871
    assert row["organisms"] == 37
    assert row["expr_edges"] == 232758
    assert row["mcp_min"] == "0.1.0"
    assert row["built_at"] == "2026-05-25T11:52:44.149Z"


# ─── _load_default_mcp_min ───────────────────────────────────────────────────
def test_load_default_mcp_min_reads_from_pyproject(rkg, tmp_path):
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text('[tool.release-kg]\nmcp_min_version = "0.7.3"\n')
    assert rkg._load_default_mcp_min(pyproject) == "0.7.3"


def test_load_default_mcp_min_repo_pyproject_value_matches_default(rkg):
    # Sanity: the value the script imports at module-load time must match what
    # the real pyproject.toml carries. Catches drift where someone edits the
    # config but the import-time DEFAULT_MCP_MIN was captured under old data.
    repo_pyproject = REPO_ROOT / "pyproject.toml"
    assert rkg._load_default_mcp_min(repo_pyproject) == rkg.DEFAULT_MCP_MIN


def test_load_default_mcp_min_falls_back_when_section_missing(rkg, tmp_path):
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text('[tool.other]\nfoo = "bar"\n')
    assert rkg._load_default_mcp_min(pyproject) == "0.1.0"


def test_load_default_mcp_min_falls_back_when_file_missing(rkg, tmp_path):
    assert rkg._load_default_mcp_min(tmp_path / "does-not-exist.toml") == "0.1.0"


def test_load_default_mcp_min_falls_back_on_malformed_toml(rkg, tmp_path):
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text("this is not [valid TOML at all\n=== oh no")
    assert rkg._load_default_mcp_min(pyproject) == "0.1.0"


# ─── parse_args ──────────────────────────────────────────────────────────────
def test_parse_args_defaults(rkg):
    ctx = rkg.parse_args(["0.1.0-alpha.1"])
    assert ctx.version == "0.1.0-alpha.1"
    assert ctx.target == "staging"
    # mcp_min default comes from pyproject [tool.release-kg].mcp_min_version
    # via DEFAULT_MCP_MIN; assert it tracks rather than hardcoding the value
    # so future config bumps don't break this test.
    assert ctx.mcp_min == rkg.DEFAULT_MCP_MIN
    assert ctx.allow_dirty is False
    assert ctx.draft is False
    assert ctx.dry_run is False
    assert ctx.resume is False
    assert ctx.skip_kg_tests is False


def test_parse_args_all_flags(rkg):
    ctx = rkg.parse_args([
        "0.1.0-alpha.1", "--target", "local",
        "--mcp-min", "0.2.0",
        "--allow-dirty", "--draft", "--dry-run", "--resume", "--skip-kg-tests",
    ])
    assert ctx.target == "local"
    assert ctx.mcp_min == "0.2.0"
    assert ctx.allow_dirty
    assert ctx.draft
    assert ctx.dry_run
    assert ctx.resume
    assert ctx.skip_kg_tests


# ─── _parse_per_publication_edges ────────────────────────────────────────────
def test_parse_per_publication_edges_basic(rkg):
    sample = (
        "doi, edges\n"
        '"10.1234/foo", 4200\n'
        '"10.5678/bar", 4000'
    )
    out = rkg._parse_per_publication_edges(sample)
    assert out == {"10.1234/foo": 4200, "10.5678/bar": 4000}


def test_parse_per_publication_edges_empty_graph(rkg):
    # Header-only output (no rows) — empty graph, not an error
    assert rkg._parse_per_publication_edges("doi, edges\n") == {}
    # Truly empty stdout — also not an error
    assert rkg._parse_per_publication_edges("") == {}


def test_parse_per_publication_edges_rsplits_on_last_comma(rkg):
    # DOI with embedded comma (rare but legal). rsplit on the last comma
    # keeps the DOI intact.
    sample = 'doi, edges\n"10.1234/foo,bar", 1234'
    out = rkg._parse_per_publication_edges(sample)
    assert out == {"10.1234/foo,bar": 1234}


def test_parse_per_publication_edges_skips_malformed(rkg):
    sample = (
        "doi, edges\n"
        '"10.1234/good", 100\n'
        "malformed_row_no_comma\n"
        '"10.5678/also_good", 200'
    )
    out = rkg._parse_per_publication_edges(sample)
    assert out == {"10.1234/good": 100, "10.5678/also_good": 200}


# ─── render_diff_block ───────────────────────────────────────────────────────
def _make_metadata(tag, papers=43, experiments=197, genes=99871, organisms=37,
                   expression_edges=232758, per_pub=None):
    return {
        "tag": tag,
        "counts": {
            "papers": papers, "experiments": experiments, "genes": genes,
            "organisms": organisms, "expression_edges": expression_edges,
        },
        "per_publication_edges": per_pub or {},
    }


def test_render_diff_block_no_changes_returns_empty(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1")
    cur = _make_metadata("kg-0.1.0-alpha.2")
    assert rkg.render_diff_block(prior, cur) == ""


def test_render_diff_block_headline_change_only(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1", papers=43, expression_edges=232758)
    cur = _make_metadata("kg-0.1.0-alpha.2", papers=45, expression_edges=240958)
    block = rkg.render_diff_block(prior, cur)
    assert "What changed since kg-0.1.0-alpha.1" in block
    assert "Headline counts" in block
    assert "+2" in block  # papers delta
    assert "+8,200" in block  # expression_edges delta with thousands separator
    # No per-publication changes → that subsection is omitted
    assert "Per-publication" not in block


def test_render_diff_block_negative_delta_uses_minus(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1", genes=100000)
    cur = _make_metadata("kg-0.1.0-alpha.2", genes=95000)
    block = rkg.render_diff_block(prior, cur)
    assert "-5,000" in block  # plain minus (Python int formatting), not en-dash


def test_render_diff_block_new_publication(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1",
                            per_pub={"10.a/old": 1000})
    cur = _make_metadata("kg-0.1.0-alpha.2",
                          per_pub={"10.a/old": 1000, "10.b/new": 500})
    block = rkg.render_diff_block(prior, cur)
    assert "New publications:" in block
    assert "10.b/new" in block
    assert "500" in block
    assert "Changed:" not in block
    assert "Removed" not in block


def test_render_diff_block_removed_publication(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1",
                            per_pub={"10.a/kept": 100, "10.b/gone": 999})
    cur = _make_metadata("kg-0.1.0-alpha.2",
                          per_pub={"10.a/kept": 100})
    block = rkg.render_diff_block(prior, cur)
    assert "Removed publications:" in block
    assert "10.b/gone" in block
    assert "999" in block


def test_render_diff_block_changed_publication(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1",
                            per_pub={"10.a/kept": 100, "10.b/grew": 800})
    cur = _make_metadata("kg-0.1.0-alpha.2",
                          per_pub={"10.a/kept": 100, "10.b/grew": 1000})
    block = rkg.render_diff_block(prior, cur)
    assert "Changed:" in block
    assert "10.b/grew" in block
    assert "800" in block and "1,000" in block
    assert "+200" in block


def test_render_diff_block_full_diff(rkg):
    prior = _make_metadata("kg-0.1.0-alpha.1", papers=43,
                            per_pub={"10.a/old": 100, "10.b/changed": 200})
    cur = _make_metadata("kg-0.1.0-alpha.2", papers=44,
                          per_pub={"10.a/old": 100, "10.b/changed": 250,
                                   "10.c/added": 75})
    block = rkg.render_diff_block(prior, cur)
    assert "Headline counts" in block
    assert "Per-publication" in block
    assert "New publications:" in block
    assert "Changed:" in block
    # No removed → that subsection is omitted
    assert "Removed publications:" not in block


def test_parse_args_rejects_bad_target(rkg):
    with pytest.raises(SystemExit):
        rkg.parse_args(["0.1.0-alpha.1", "--target", "nope"])
