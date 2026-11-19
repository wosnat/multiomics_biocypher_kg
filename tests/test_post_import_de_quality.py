"""Static check: post-import scripts compute the DE table-quality arrays."""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FILES = [ROOT / "scripts" / "post-import.sh", ROOT / "scripts" / "post-import.cypher"]
PROPS = [
    "time_point_negative_fractions",
    "time_point_median_log2fc",
    "time_point_max_abs_log2fc",
    "time_point_large_fc_not_significant",
    "time_point_extreme_fc_count",
    "time_point_low_signal_count",
]


def _code_lines(path):
    keys = PROPS + ["data_caveats"]
    out = set()
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("//"):
            continue
        if any(k in s for k in keys):
            out.add(s)
    return out


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_post_import_sets_de_quality_arrays(path):
    text = path.read_text()
    for prop in PROPS:
        assert prop in text, prop
    assert "data_caveats = coalesce(e.data_caveats, [])" in text


def test_post_import_files_agree_on_de_quality_lines():
    assert _code_lines(FILES[0]) == _code_lines(FILES[1])
