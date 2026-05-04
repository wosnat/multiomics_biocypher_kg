"""Tests for sub-step 6 download module.

HTTP is mocked via monkeypatch; tests assert that the right URLs are
requested and that bytes are written to the right cache paths.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from multiomics_kg.download import download_metabolism_reference as dmr


class _FakeResponse:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code
        self._chunks = [content]

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int = 8192):
        return iter(self._chunks)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _patch_requests_get(monkeypatch, response_map: dict[str, _FakeResponse]):
    """Patch requests.get to look up responses by URL."""
    def fake_get(url, *args, **kwargs):
        if url not in response_map:
            raise AssertionError(f"unexpected URL: {url}")
        return response_map[url]
    monkeypatch.setattr(dmr.requests, "get", fake_get)


# ---------------------------------------------------------------------------
# Structure tests
# ---------------------------------------------------------------------------

def test_sources_by_group_keys():
    assert set(dmr.SOURCES_BY_GROUP.keys()) == {"mnx", "tcdb"}
    assert len(dmr.MNX_SOURCES) == 4
    assert len(dmr.TCDB_SOURCES) == 4


def test_group_dir_layout(tmp_path, monkeypatch):
    """MNX honors MNX_DATA_DIR; TCDB lives under <cache_root>/tcdb/raw/."""
    monkeypatch.setenv("MNX_DATA_DIR", str(tmp_path / "shared_mnx"))
    assert dmr._group_dir("mnx", tmp_path) == tmp_path / "shared_mnx"
    assert dmr._group_dir("tcdb", tmp_path) == tmp_path / "tcdb" / "raw"
    with pytest.raises(ValueError, match="Unknown source group"):
        dmr._group_dir("bogus", tmp_path)


# ---------------------------------------------------------------------------
# download_all filtering tests
# ---------------------------------------------------------------------------

def test_download_all_defaults_to_all_groups(tmp_path, monkeypatch):
    """No --sources arg → downloads both groups."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=None)
    assert len(called) == 8  # 4 MNX + 4 TCDB


def test_download_all_filters_to_mnx(tmp_path, monkeypatch):
    """--sources mnx → only MNX URLs."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=["mnx"])
    assert len(called) == 4
    assert all("metanetx.org" in u for u in called)


def test_download_all_filters_to_tcdb(tmp_path, monkeypatch):
    """--sources tcdb → only TCDB URLs."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=["tcdb"])
    assert len(called) == 4
    assert all("tcdb.org" in u for u in called)


def test_download_all_invalid_source_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown source group"):
        dmr.download_all(cache_root=tmp_path, sources=["bogus"])


def test_download_all_empty_sources_downloads_nothing(tmp_path, monkeypatch):
    """sources=[] is distinct from sources=None: empty list = download nothing."""
    called: list[str] = []
    monkeypatch.setattr(dmr, "download_one",
                        lambda url, dest, force: (called.append(url), True)[1])
    dmr.download_all(cache_root=tmp_path, force=False, sources=[])
    assert called == []  # no downloads triggered


# ---------------------------------------------------------------------------
# Full download + caching behaviour (HTTP-mocked)
# ---------------------------------------------------------------------------

def _expected_path(group: str, filename: str, cache_root: Path, mnx_dir: Path) -> Path:
    if group == "mnx":
        return mnx_dir / filename
    return cache_root / "tcdb" / "raw" / filename


def test_download_writes_all_files(monkeypatch, tmp_path):
    """One download per source; each writes to its group-specific destination."""
    mnx_dir = tmp_path / "shared_mnx"
    monkeypatch.setenv("MNX_DATA_DIR", str(mnx_dir))

    all_sources = [
        (group, key, url, filename)
        for group, sources in dmr.SOURCES_BY_GROUP.items()
        for key, (url, filename) in sources.items()
    ]
    fake_responses = {
        url: _FakeResponse(f"BODY-{key}".encode())
        for _group, key, url, _filename in all_sources
    }
    _patch_requests_get(monkeypatch, fake_responses)

    dmr.download_all(cache_root=tmp_path, force=True)

    for group, key, _url, filename in all_sources:
        out = _expected_path(group, filename, tmp_path, mnx_dir)
        assert out.exists(), f"missing output for {key}: {out}"
        assert out.read_bytes() == f"BODY-{key}".encode()


def test_skip_when_cached_unless_force(monkeypatch, tmp_path):
    """Existing files are skipped when force=False; re-downloaded when force=True."""
    mnx_dir = tmp_path / "shared_mnx"
    monkeypatch.setenv("MNX_DATA_DIR", str(mnx_dir))

    # Pre-create the MNX chem_prop cache file with stale content
    target = mnx_dir / dmr.MNX_SOURCES["mnx_chem_prop"][1]
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b"STALE")

    # Mock all URLs to return fresh content
    all_urls = [
        url
        for sources in dmr.SOURCES_BY_GROUP.values()
        for url, _filename in sources.values()
    ]
    fake_responses = {url: _FakeResponse(b"FRESH") for url in all_urls}
    _patch_requests_get(monkeypatch, fake_responses)

    # Without --force: stale content preserved
    dmr.download_all(cache_root=tmp_path, force=False)
    assert target.read_bytes() == b"STALE"

    # With --force: stale content overwritten
    dmr.download_all(cache_root=tmp_path, force=True)
    assert target.read_bytes() == b"FRESH"
