"""Unit tests for the 5 new KEGG endpoints introduced in Spec 1.2."""
from __future__ import annotations

import textwrap

from multiomics_kg.utils import kegg_utils


REACTION_LIST_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpyruvate kinase reaction
    rn:R00010\tphosphofructokinase reaction
    invalid_line
    rn:R12345\t
""")

COMPOUND_LIST_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tD-glucose; alpha-D-glucopyranose
    cpd:C00002\tATP; adenosine 5'-triphosphate
    bad_prefix:Cxxxx\tjunk
""")

LINK_CR_FIXTURE = textwrap.dedent("""\
    rn:R00200\tcpd:C00074
    rn:R00200\tcpd:C00008
    rn:R00010\tcpd:C00031
""")

LINK_PR_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpath:rn00010
    rn:R00200\tpath:rn00710
    rn:R00010\tpath:rn00010
""")

LINK_PC_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tpath:map00010
    cpd:C00031\tpath:map00500
""")


def test_parse_reaction_names():
    out = kegg_utils._parse_reaction_names(REACTION_LIST_FIXTURE)
    assert out == {
        "R00200": "pyruvate kinase reaction",
        "R00010": "phosphofructokinase reaction",
        "R12345": "",
    }


def test_parse_compound_names():
    out = kegg_utils._parse_compound_names(COMPOUND_LIST_FIXTURE)
    # Multi-name field: take first synonym (before first '; ')
    assert out["C00031"] == "D-glucose"
    assert out["C00002"] == "ATP"
    assert "Cxxxx" not in out


def test_parse_reaction_to_compounds():
    out = kegg_utils._parse_reaction_to_compounds(LINK_CR_FIXTURE)
    assert sorted(out["R00200"]) == ["C00008", "C00074"]
    assert out["R00010"] == ["C00031"]


def test_parse_reaction_to_pathways_strips_rn_prefix():
    """KEGG `/link/pathway/reaction` returns rn-prefixed pathway IDs.
    We normalize them to the ko-prefixed form used elsewhere in the KG.
    """
    out = kegg_utils._parse_reaction_to_pathways(LINK_PR_FIXTURE)
    assert sorted(out["R00200"]) == ["ko00010", "ko00710"]
    assert out["R00010"] == ["ko00010"]


def test_parse_compound_to_pathways_strips_map_prefix():
    out = kegg_utils._parse_compound_to_pathways(LINK_PC_FIXTURE)
    # map-prefixed pathways → ko-prefixed for consistency with KeggTerm node IDs
    assert sorted(out["C00031"]) == ["ko00010", "ko00500"]


def test_download_kegg_data_includes_metabolism_keys(tmp_path, monkeypatch):
    """download_kegg_data should populate the 5 new metabolism keys."""

    # Stub the network: return synthetic fixtures for every URL we know about.
    def fake_fetch_text(url: str) -> str:
        if url.endswith("/list/ko"):
            return "ko:K02338\tDNA polymerase III\n"
        if url.endswith("/link/pathway/ko"):
            return "ko:K02338\tpath:ko03030\n"
        if url.endswith("/list/pathway/ko"):
            return "path:ko03030\tDNA replication\n"
        if url.endswith("/list/reaction"):
            return REACTION_LIST_FIXTURE
        if url.endswith("/list/compound"):
            return COMPOUND_LIST_FIXTURE
        if url.endswith("/link/compound/reaction"):
            return LINK_CR_FIXTURE
        if url.endswith("/link/pathway/reaction"):
            return LINK_PR_FIXTURE
        if url.endswith("/link/pathway/compound"):
            return LINK_PC_FIXTURE
        raise AssertionError(f"unexpected URL {url}")

    def fake_fetch_json(url: str) -> dict:
        # Minimal BRITE skeleton — children empty, parsers handle that gracefully
        return {"children": []}

    monkeypatch.setattr(kegg_utils, "_fetch_text", fake_fetch_text)
    monkeypatch.setattr(kegg_utils, "_fetch_json", fake_fetch_json)

    data = kegg_utils.download_kegg_data(tmp_path, force=True)

    # Pre-existing keys still present
    for key in ("ko_names", "pathway_names", "ko_to_pathways"):
        assert key in data

    # New Spec 1.2 keys
    assert data["reaction_names"] == {
        "R00200": "pyruvate kinase reaction",
        "R00010": "phosphofructokinase reaction",
        "R12345": "",
    }
    assert data["compound_names"]["C00031"] == "D-glucose"
    assert sorted(data["reaction_to_compounds"]["R00200"]) == ["C00008", "C00074"]
    assert sorted(data["reaction_to_pathways"]["R00200"]) == ["ko00010", "ko00710"]
    assert sorted(data["compound_to_pathways"]["C00031"]) == ["ko00010", "ko00500"]

    # Cache file written and round-trips
    cache_file = tmp_path / "kegg" / "kegg_data.json"
    assert cache_file.exists()
    import json
    reloaded = json.loads(cache_file.read_text())
    assert reloaded["reaction_names"] == data["reaction_names"]


def test_download_kegg_raw_populates_raw_dir(tmp_path, monkeypatch):
    """download_kegg_raw should populate cache/data/kegg/raw/ with 8 .txt + 1 .json."""
    text_calls: list[str] = []
    json_calls: list[str] = []

    def fake_text(url):
        text_calls.append(url)
        return "stub"

    def fake_json(url):
        json_calls.append(url)
        return {"children": []}

    monkeypatch.setattr(kegg_utils, "_fetch_text", fake_text)
    monkeypatch.setattr(kegg_utils, "_fetch_json", fake_json)

    kegg_utils.download_kegg_raw(tmp_path, force=False)
    raw_dir = tmp_path / "kegg" / "raw"
    txt_files = sorted(p.name for p in raw_dir.glob("*.txt"))
    json_files = sorted(p.name for p in raw_dir.glob("*.json"))
    assert len(txt_files) == 8
    assert len(json_files) == 1
    assert "br_ko00001.json" in json_files
    # Anchor: the patched fakes were actually invoked (not a real network call)
    assert text_calls, "expected _fetch_text to be invoked on first pass"
    assert json_calls, "expected _fetch_json to be invoked on first pass"

    # Second call without force: no new fetches (raw cache is read from disk)
    text_calls.clear()
    json_calls.clear()
    kegg_utils.download_kegg_raw(tmp_path, force=False)
    assert text_calls == []
    assert json_calls == []

    # Force: re-fetches everything
    kegg_utils.download_kegg_raw(tmp_path, force=True)
    assert len(text_calls) == 8
    assert len(json_calls) == 1


def test_download_kegg_data_uses_raw_cache_on_second_call(tmp_path, monkeypatch):
    """Second call with force=False must not hit the network — raw cache is reused."""
    fetch_count = {"text": 0, "json": 0}

    def counting_fetch_text(url):
        fetch_count["text"] += 1
        if url.endswith("/list/ko"):
            return "ko:K02338\tDNA polymerase III\n"
        if url.endswith("/link/pathway/ko"):
            return "ko:K02338\tpath:ko03030\n"
        if url.endswith("/list/pathway/ko"):
            return "path:ko03030\tDNA replication\n"
        if url.endswith("/list/reaction"):
            return REACTION_LIST_FIXTURE
        if url.endswith("/list/compound"):
            return COMPOUND_LIST_FIXTURE
        if url.endswith("/link/compound/reaction"):
            return LINK_CR_FIXTURE
        if url.endswith("/link/pathway/reaction"):
            return LINK_PR_FIXTURE
        if url.endswith("/link/pathway/compound"):
            return LINK_PC_FIXTURE
        raise AssertionError(f"unexpected URL {url}")

    def counting_fetch_json(url):
        fetch_count["json"] += 1
        return {"children": []}

    monkeypatch.setattr(kegg_utils, "_fetch_text", counting_fetch_text)
    monkeypatch.setattr(kegg_utils, "_fetch_json", counting_fetch_json)

    # First call: 8 text fetches + 1 JSON fetch
    kegg_utils.download_kegg_data(tmp_path, force=True)
    first_text = fetch_count["text"]
    first_json = fetch_count["json"]
    assert first_text == 8
    assert first_json == 1

    # Delete kegg_data.json so download_kegg_data has to re-parse
    (tmp_path / "kegg" / "kegg_data.json").unlink()

    # Second call: raw cache exists, no network, but parses from disk
    kegg_utils.download_kegg_data(tmp_path, force=False)
    assert fetch_count["text"] == first_text  # no new text fetches
    assert fetch_count["json"] == first_json  # no new JSON fetches

    # Third call with force=True: re-fetches everything
    kegg_utils.download_kegg_data(tmp_path, force=True)
    assert fetch_count["text"] == first_text + 8
    assert fetch_count["json"] == first_json + 1
