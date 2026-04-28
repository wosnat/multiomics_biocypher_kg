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


def test_sources_table_has_six_entries():
    """Six expected sources: 4 MNX TSVs + 2 TCDB TSVs.

    CAZy is intentionally not downloaded — its hierarchy is bootstrapped from
    observed eggNOG `CAZy` columns + mechanical ID parsing in Phase 1.1B.
    """
    assert len(dmr.SOURCES) == 6
    assert {"mnx_chem_prop", "mnx_chem_xref",
            "mnx_reac_prop", "mnx_reac_xref",
            "tcdb_families", "tcdb_substrates"} == set(dmr.SOURCES)


def test_download_writes_all_files(monkeypatch, tmp_path):
    """One download per source; each writes to the configured relative path."""
    fake_responses = {
        url: _FakeResponse(f"BODY-{key}".encode())
        for key, (url, _rel) in dmr.SOURCES.items()
    }
    _patch_requests_get(monkeypatch, fake_responses)

    dmr.download_all(cache_root=tmp_path, force=True)

    for key, (_url, rel) in dmr.SOURCES.items():
        out = tmp_path / rel
        assert out.exists(), f"missing output for {key}: {out}"
        assert out.read_bytes() == f"BODY-{key}".encode()


def test_skip_when_cached_unless_force(monkeypatch, tmp_path):
    """Existing files are skipped when force=False; re-downloaded when force=True."""
    # Pre-create the MNX chem_prop cache file with stale content
    target = tmp_path / dmr.SOURCES["mnx_chem_prop"][1]
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b"STALE")

    # Mock all URLs to return fresh content
    fake_responses = {
        url: _FakeResponse(b"FRESH")
        for url, _rel in dmr.SOURCES.values()
    }
    _patch_requests_get(monkeypatch, fake_responses)

    # Without --force: stale content preserved
    dmr.download_all(cache_root=tmp_path, force=False)
    assert target.read_bytes() == b"STALE"

    # With --force: stale content overwritten
    dmr.download_all(cache_root=tmp_path, force=True)
    assert target.read_bytes() == b"FRESH"
