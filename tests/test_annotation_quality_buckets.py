"""Fast static gate for the annotation_quality/annotation_state source-bucket
list (F1.2/F1.3, see CLAUDE.md "Source bucket maintenance" and
docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md).

The bucket list is explicitly enumerated in three places inside the
post-import Cypher: a human-readable `SOURCE_BUCKETS` comment block, the
`has_<bucket>` EXISTS lines, and the `informative_count` CASE-sum. All three
must agree, and `scripts/post-import.cypher` (the reference copy used
outside Docker) must stay byte-identical in this logic to the authoritative
`scripts/post-import.sh` (see CLAUDE.md's Docker Pipeline Stages / step 3).

This test parses both files statically — no Neo4j required — so it belongs
in the fast suite (`pytest -m "not slow and not kg"`) and acts as the
"bucket-count test" the maintenance procedure promises bumping whenever a
bucket is added or removed.
"""
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
CYPHER_FILE = REPO_ROOT / "scripts" / "post-import.cypher"
SH_FILE = REPO_ROOT / "scripts" / "post-import.sh"

EXPECTED_BUCKET_COUNT = 9
EXPECTED_BUCKETS = {
    "go",
    "kegg",
    "pfam",
    "ec",
    "role",
    "reaction",
    "transporter",
    "cazy",
    "ncbifam",
}

SOURCE_BUCKETS_RE = re.compile(
    r"SOURCE_BUCKETS:start.*?live \((\d+)\):\s*(.+?)\n.*?SOURCE_BUCKETS:end",
    re.DOTALL,
)

# The informative_count sum: a run of `CASE WHEN has_<bucket> THEN 1 ELSE 0 END`
# terms combined with `+`, immediately followed by `) AS informative_count,`.
INFORMATIVE_COUNT_RE = re.compile(
    r"\(CASE WHEN has_(\w+) THEN 1 ELSE 0 END\s*"
    r"((?:\+\s*CASE WHEN has_\w+ THEN 1 ELSE 0 END\s*)*)"
    r"\)\s*AS informative_count",
)
HAS_TERM_RE = re.compile(r"CASE WHEN has_(\w+) THEN 1 ELSE 0 END")


def _parse_source_buckets_comment(text: str) -> list[str]:
    m = SOURCE_BUCKETS_RE.search(text)
    assert m, "SOURCE_BUCKETS:start/end comment block not found"
    declared_count = int(m.group(1))
    buckets = [b.strip() for b in m.group(2).split(",")]
    assert len(buckets) == declared_count, (
        f"SOURCE_BUCKETS comment says live ({declared_count}) but lists "
        f"{len(buckets)}: {buckets}"
    )
    return buckets


def _parse_informative_count_buckets(text: str) -> list[str]:
    m = INFORMATIVE_COUNT_RE.search(text)
    assert m, "informative_count CASE-sum block not found"
    full_block = m.group(0)
    return HAS_TERM_RE.findall(full_block)


@pytest.mark.parametrize("path", [CYPHER_FILE, SH_FILE], ids=["cypher", "sh"])
def test_source_bucket_comment_matches_expected(path):
    text = path.read_text()
    buckets = _parse_source_buckets_comment(text)
    assert len(buckets) == EXPECTED_BUCKET_COUNT, (
        f"{path.name}: expected {EXPECTED_BUCKET_COUNT} source buckets, "
        f"found {len(buckets)}: {buckets}"
    )
    assert set(buckets) == EXPECTED_BUCKETS, (
        f"{path.name}: bucket set mismatch. "
        f"missing={EXPECTED_BUCKETS - set(buckets)} "
        f"extra={set(buckets) - EXPECTED_BUCKETS}"
    )


@pytest.mark.parametrize("path", [CYPHER_FILE, SH_FILE], ids=["cypher", "sh"])
def test_informative_count_sum_matches_expected(path):
    """The has_<bucket> terms actually summed into informative_count must
    match the declared SOURCE_BUCKETS list 1:1 (order-independent) — this is
    what makes the comment block a true description of the logic, not just
    documentation that can drift."""
    text = path.read_text()
    declared = _parse_source_buckets_comment(text)
    summed = _parse_informative_count_buckets(text)
    assert len(summed) == EXPECTED_BUCKET_COUNT, (
        f"{path.name}: informative_count sums {len(summed)} has_<bucket> "
        f"terms, expected {EXPECTED_BUCKET_COUNT}: {summed}"
    )
    assert set(summed) == set(declared), (
        f"{path.name}: informative_count buckets {sorted(summed)} != "
        f"SOURCE_BUCKETS comment buckets {sorted(declared)}"
    )
    assert set(summed) == EXPECTED_BUCKETS


def test_cypher_and_sh_agree_on_bucket_logic():
    """scripts/post-import.cypher is a reference copy of scripts/post-import.sh
    kept in sync for non-Docker use (see CLAUDE.md Docker Pipeline Stages,
    step 3). Both must declare and sum the exact same bucket set."""
    cypher_text = CYPHER_FILE.read_text()
    sh_text = SH_FILE.read_text()

    cypher_declared = set(_parse_source_buckets_comment(cypher_text))
    sh_declared = set(_parse_source_buckets_comment(sh_text))
    assert cypher_declared == sh_declared, (
        f"SOURCE_BUCKETS comment differs between post-import.cypher and "
        f"post-import.sh: {cypher_declared} vs {sh_declared}"
    )

    cypher_summed = set(_parse_informative_count_buckets(cypher_text))
    sh_summed = set(_parse_informative_count_buckets(sh_text))
    assert cypher_summed == sh_summed, (
        f"informative_count bucket sum differs between post-import.cypher "
        f"and post-import.sh: {cypher_summed} vs {sh_summed}"
    )
