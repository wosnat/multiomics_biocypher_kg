import json
from pathlib import Path

import pytest

from multiomics_kg.utils.go_plus import (
    build_go_transport_cache, go_plus_release, load_target_filter, parse_transport_axioms)

FIX = Path(__file__).parent / "fixtures" / "go_plus_mini.json"
EXPECTED = {"GO:0015110": ["CHEBI:29195"], "GO:0015704": ["CHEBI:29195"],
            "GO:0005274": ["CHEBI:15378", "CHEBI:15676"], "GO:0015891": ["CHEBI:26672"]}


@pytest.fixture
def mini():
    return json.loads(FIX.read_text())["graphs"][0]


def _setup(tmp_path, release, raw=True):
    d = tmp_path / "go_terms"
    (d / "raw").mkdir(parents=True)
    if raw:
        (d / "raw" / "go-plus.json").write_text(FIX.read_text())
    (d / "go_release.txt").write_text(release + "\n")
    return tmp_path


def test_parse_transport_axioms(mini):
    assert parse_transport_axioms(mini) == EXPECTED


def test_release(mini):
    assert go_plus_release(mini) == "2026-07-26"


def test_target_filter_contents():
    f = load_target_filter(Path("config/go_transport_target_filter.yaml"))
    assert set(f) == {f"CHEBI:{n}" for n in (24870, 24867, 23906, 23905, 25213, 33515, 51143,
                                              33694, 24431, 24834, 15378, 29101)}
    # protein left the filter 2026-10-25: it is the substrate of protein translocation
    assert "CHEBI:36080" not in f
    assert "CHEBI:29103" not in f and "CHEBI:17996" not in f
    assert f["CHEBI:15378"] == "hydron"


def test_build_cache_from_raw(tmp_path):
    root = _setup(tmp_path, "2026-07-26")
    c = build_go_transport_cache(root)
    assert c["release"] == "2026-07-26" and c["axioms"] == EXPECTED and len(c["md5"]) == 32
    on_disk = json.loads((root / "go_terms" / "go_transport_substrates.json").read_text())
    assert on_disk == c


def test_build_cache_release_mismatch_raises(tmp_path):
    root = _setup(tmp_path, "2026-03-01")
    with pytest.raises(RuntimeError):
        build_go_transport_cache(root)


def test_build_cache_offline_reuses_committed(tmp_path):
    root = _setup(tmp_path, "2026-07-26", raw=False)
    committed = {"release": "2026-07-26", "md5": "x", "axioms": EXPECTED}
    (root / "go_terms" / "go_transport_substrates.json").write_text(json.dumps(committed))
    assert build_go_transport_cache(root) == committed


def test_build_cache_offline_stale_committed_raises(tmp_path):
    root = _setup(tmp_path, "2026-07-26", raw=False)
    (root / "go_terms" / "go_transport_substrates.json").write_text(
        json.dumps({"release": "old", "md5": "x", "axioms": {}}))
    with pytest.raises(RuntimeError, match="refresh_go_release"):
        build_go_transport_cache(root)
