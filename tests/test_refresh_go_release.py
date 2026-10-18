"""Offline tests for refresh_go_release (fake downloader)."""
import json
import pytest

from multiomics_kg.download import refresh_go_release as m
from multiomics_kg.utils.go_utils import read_go_release

OBO = """format-version: 1.2
data-version: releases/{rel}

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
alt_id: GO:0000099
"""


def _plus(rel):
    return json.dumps({"graphs": [{"meta": {"basicPropertyValues": [
        {"pred": "http://www.w3.org/2002/07/owl#versionInfo", "val": rel}]}}]})


def _fake(monkeypatch, basic_rel, plus_rel):
    calls = []

    def fake(url, dest):
        calls.append(url)
        dest.write_text(OBO.format(rel=basic_rel) if url == m.GO_BASIC_URL else _plus(plus_rel))

    monkeypatch.setattr(m, "_download", fake)
    return calls


def test_refresh_writes_all_four_files_and_release(tmp_path, monkeypatch):
    _fake(monkeypatch, "2026-07-26", "2026-07-26")
    assert m.refresh_go_release(tmp_path) == "2026-07-26"
    g = tmp_path / "go_terms"
    for f in ("raw/go-basic.obo", "raw/go-plus.json", "go_namespace_cache.json", "go_release.txt"):
        assert (g / f).exists()
    assert read_go_release(tmp_path) == "2026-07-26"
    assert json.loads((g / "go_alt_ids.json").read_text()) == {"GO:0000099": "GO:0008150"}
    assert "GO:0008150" in json.loads((g / "go_namespace_cache.json").read_text())


def test_release_mismatch_writes_nothing(tmp_path, monkeypatch):
    _fake(monkeypatch, "2026-07-26", "2026-08-10")
    with pytest.raises(m.GoReleaseMismatch):
        m.refresh_go_release(tmp_path)
    assert not [p for p in tmp_path.rglob("*") if p.is_file()]


def test_no_force_no_download(tmp_path, monkeypatch):
    _fake(monkeypatch, "2026-07-26", "2026-07-26")
    m.refresh_go_release(tmp_path)
    calls = _fake(monkeypatch, "2026-09-01", "2026-09-01")
    assert m.refresh_go_release(tmp_path) == "2026-07-26"
    assert calls == []
    m.refresh_go_release(tmp_path, force=True)
    assert len(calls) == 2
