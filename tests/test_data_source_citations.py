"""Unit tests for config/data_source_citations.yaml and its loader.

The yaml is the single source of per-DataSource metadata (name, description,
citations, dois, url, license, version) plus the build framework's citations
stamped on Schema_info. See docs/superpowers/specs/2026-11-18-data-source-citations-design.md.
"""

import re
from pathlib import Path

import pytest
import yaml

from multiomics_kg.adapters.data_source_adapter import data_source_ids
from multiomics_kg.utils import data_source_citations as dsc
from multiomics_kg.utils.data_source_citations import (
    FORBIDDEN_CHARS,
    load_citations,
    resolve_version,
)

DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


@pytest.fixture(scope="module")
def citations():
    return load_citations()


def _all_strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from _all_strings(k)
            yield from _all_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _all_strings(v)


def test_yaml_and_data_sources_match(citations):
    with open("config/gene_annotations_config.yaml") as f:
        expected = data_source_ids(yaml.safe_load(f))
    got = set(citations["sources"])
    assert got == expected, (
        f"missing yaml entries: {sorted(expected - got)}; "
        f"yaml entries with no DataSource: {sorted(got - expected)}"
    )


def test_every_entry_complete(citations):
    for sid, e in citations["sources"].items():
        assert e.get("name"), sid
        assert e.get("description"), sid
        assert e.get("url"), sid
        cites = e.get("citations")
        assert isinstance(cites, list) and cites and all(isinstance(c, str) and c for c in cites), sid
        dois = e.get("dois")
        assert isinstance(dois, list) and len(dois) == len(cites), sid
        for d in dois:
            assert d == "" or DOI_RE.match(d), f"{sid}: bad doi {d!r}"


def test_curated_text_has_no_forbidden_chars(citations):
    for s in _all_strings(citations):
        bad = [c for c in FORBIDDEN_CHARS if c in s]
        assert not bad, f"forbidden {bad} in {s!r}"


def test_interpro_version_matches_runner_pin(citations):
    src = Path(".claude/skills/interproscan-run/run_interproscan.py").read_text()
    pin = re.search(r'INTERPRO_RELEASE = "([^"]+)"', src).group(1)
    assert citations["sources"]["interproscan"]["version"] == pin


def test_version_from_missing_file_is_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(dsc, "GO_RELEASE_PATH", tmp_path / "missing.txt")
    monkeypatch.setattr(dsc, "CHEBI_CACHE_PATH", tmp_path / "missing.json")
    assert resolve_version({"version_from": "go_release"}) == ""
    assert resolve_version({"version_from": "chebi_cache"}) == ""


def test_unknown_version_from_raises():
    with pytest.raises(ValueError):
        resolve_version({"version_from": "nope"})


def test_literal_version_and_default():
    assert resolve_version({"version": "12.5"}) == "12.5"
    assert resolve_version({}) == ""


def test_go_release_reader():
    p = Path("cache/data/go_terms/go_release.txt")
    if not p.exists():
        pytest.skip("go_release.txt not present")
    assert resolve_version({"version_from": "go_release"}) == p.read_text().strip() != ""


def test_build_citations_nonempty(citations):
    build = citations["build_citations"]
    assert isinstance(build, list) and build and all(isinstance(c, str) and c for c in build)
    assert any("BioCypher" in c for c in build)
    # DOI carried inside the text: Schema_info has no parallel dois list
    assert all(re.search(r"doi: 10\.\d{4,9}/\S+$", c) for c in build)


def test_write_build_citations_sidecar(tmp_path):
    import json
    from multiomics_kg.utils.data_source_citations import write_build_citations
    path = write_build_citations(tmp_path)
    assert path == tmp_path / "build_citations.json"
    assert json.loads(path.read_text()) == load_citations()["build_citations"]
