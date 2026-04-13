"""Unit tests for brite_utils.py — all HTTP calls mocked."""
import json
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from multiomics_kg.utils.brite_utils import (
    BRITE_TREES,
    compute_level_kind,
    download_brite_tree,
    load_brite_trees,
)


# ---------------------------------------------------------------------------
# compute_level_kind
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("depth,expected", [
    (0, "brite_class"),
    (1, "brite_subclass"),
    (2, "brite_family"),
    (3, "brite_subfamily"),
])
def test_compute_level_kind(depth, expected):
    assert compute_level_kind(depth) == expected


def test_compute_level_kind_too_deep():
    with pytest.raises(ValueError, match="4"):
        compute_level_kind(4)


# ---------------------------------------------------------------------------
# download_brite_tree
# ---------------------------------------------------------------------------

SAMPLE_BRITE_JSON = {
    "name": "ko02000",
    "children": [
        {
            "name": "1. ABC Transporters",
            "children": [
                {
                    "name": "Phosphate transport",
                    "children": [
                        {"name": "K02036  pstB; phosphate transport"},
                    ],
                }
            ],
        }
    ],
}


def test_download_brite_tree_fetches_correct_url(tmp_path):
    """Verifies correct KEGG REST URL is used."""
    with patch(
        "multiomics_kg.utils.brite_utils._fetch_json",
        return_value=SAMPLE_BRITE_JSON,
    ) as mock_fetch:
        result = download_brite_tree("ko02000", tmp_path, cache=False)

    mock_fetch.assert_called_once_with(
        "https://rest.kegg.jp/get/br:ko02000/json"
    )
    assert result == SAMPLE_BRITE_JSON


def test_download_brite_tree_writes_cache(tmp_path):
    """First fetch writes a cache file at the expected path."""
    with patch(
        "multiomics_kg.utils.brite_utils._fetch_json",
        return_value=SAMPLE_BRITE_JSON,
    ):
        download_brite_tree("ko02000", tmp_path, cache=True)

    cache_file = tmp_path / "kegg" / "brite_ko02000.json"
    assert cache_file.exists()
    with open(cache_file) as fh:
        assert json.load(fh) == SAMPLE_BRITE_JSON


def test_download_brite_tree_cache_hit_skips_http(tmp_path):
    """Second call with cache=True reads cache file, makes no HTTP call."""
    cache_dir = tmp_path / "kegg"
    cache_dir.mkdir()
    (cache_dir / "brite_ko02000.json").write_text(json.dumps(SAMPLE_BRITE_JSON))

    with patch("multiomics_kg.utils.brite_utils._fetch_json") as mock_fetch:
        result = download_brite_tree("ko02000", tmp_path, cache=True)

    mock_fetch.assert_not_called()
    assert result == SAMPLE_BRITE_JSON


def test_download_brite_tree_cache_false_forces_refetch(tmp_path):
    """cache=False re-fetches even when cache file exists."""
    cache_dir = tmp_path / "kegg"
    cache_dir.mkdir()
    stale = {"name": "ko02000", "children": []}
    (cache_dir / "brite_ko02000.json").write_text(json.dumps(stale))

    with patch(
        "multiomics_kg.utils.brite_utils._fetch_json",
        return_value=SAMPLE_BRITE_JSON,
    ) as mock_fetch:
        result = download_brite_tree("ko02000", tmp_path, cache=False)

    mock_fetch.assert_called_once()
    assert result == SAMPLE_BRITE_JSON


def test_download_brite_tree_http_error_raises(tmp_path):
    """Non-200 response propagates as HTTPError."""
    with patch(
        "multiomics_kg.utils.brite_utils._fetch_json",
        side_effect=requests.HTTPError("404 Not Found"),
    ):
        with pytest.raises(requests.HTTPError):
            download_brite_tree("ko99999", tmp_path, cache=False)


# ---------------------------------------------------------------------------
# load_brite_trees
# ---------------------------------------------------------------------------

def test_load_brite_trees_returns_keyed_dict(tmp_path):
    """Returns {tree_code: parsed_json} for all requested trees."""
    fixtures = {
        "ko02000": SAMPLE_BRITE_JSON,
        "ko01002": {"name": "ko01002", "children": []},
    }
    with patch("multiomics_kg.utils.brite_utils.download_brite_tree") as mock_dl:
        mock_dl.side_effect = lambda code, root, **kwargs: fixtures[code]
        result = load_brite_trees(tmp_path, ["ko02000", "ko01002"])

    assert set(result.keys()) == {"ko02000", "ko01002"}
    assert result["ko02000"] == SAMPLE_BRITE_JSON


def test_load_brite_trees_rate_limit_sleep(tmp_path):
    """A sleep is inserted between sequential tree fetches (n-1 times for n trees)."""
    with patch("multiomics_kg.utils.brite_utils.download_brite_tree",
               return_value={}):
        with patch("multiomics_kg.utils.brite_utils.time.sleep") as mock_sleep:
            load_brite_trees(tmp_path, ["ko02000", "ko01002", "ko03000"])

    assert mock_sleep.call_count == 2  # n-1 sleeps for n=3 trees


def test_load_brite_trees_single_tree_no_sleep(tmp_path):
    """A single tree fetch requires no sleep."""
    with patch("multiomics_kg.utils.brite_utils.download_brite_tree",
               return_value={}):
        with patch("multiomics_kg.utils.brite_utils.time.sleep") as mock_sleep:
            load_brite_trees(tmp_path, ["ko02000"])

    mock_sleep.assert_not_called()


# ---------------------------------------------------------------------------
# BRITE_TREES constant
# ---------------------------------------------------------------------------

def test_brite_trees_has_12_entries():
    assert len(BRITE_TREES) == 12


def test_brite_trees_keys_are_ko_prefixed():
    for code in BRITE_TREES:
        assert code.startswith("ko"), f"Expected ko-prefixed code, got {code!r}"


def test_brite_trees_values_are_snake_case_strings():
    for code, name in BRITE_TREES.items():
        assert isinstance(name, str)
        assert " " not in name, f"Tree name {name!r} contains spaces (expected snake_case)"
