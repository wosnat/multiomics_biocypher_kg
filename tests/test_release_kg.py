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


# ─── parse_args ──────────────────────────────────────────────────────────────
def test_parse_args_defaults(rkg):
    ctx = rkg.parse_args(["0.1.0-alpha.1"])
    assert ctx.version == "0.1.0-alpha.1"
    assert ctx.target == "staging"
    assert ctx.mcp_min == "0.1.0"
    assert ctx.allow_dirty is False
    assert ctx.draft is False
    assert ctx.dry_run is False
    assert ctx.resume is False


def test_parse_args_all_flags(rkg):
    ctx = rkg.parse_args([
        "0.1.0-alpha.1", "--target", "local",
        "--mcp-min", "0.2.0",
        "--allow-dirty", "--draft", "--dry-run", "--resume",
    ])
    assert ctx.target == "local"
    assert ctx.mcp_min == "0.2.0"
    assert ctx.allow_dirty
    assert ctx.draft
    assert ctx.dry_run
    assert ctx.resume


def test_parse_args_rejects_bad_target(rkg):
    with pytest.raises(SystemExit):
        rkg.parse_args(["0.1.0-alpha.1", "--target", "nope"])
