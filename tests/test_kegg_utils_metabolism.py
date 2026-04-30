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
    cpd:C00074\trn:R00200
    cpd:C00008\trn:R00200
    cpd:C00031\trn:R00010
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
