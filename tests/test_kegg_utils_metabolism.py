"""Unit tests for the 5 new KEGG endpoints introduced in Spec 1.2."""
from __future__ import annotations

import json
import textwrap

import pytest

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


def test_load_kegg_data_returns_pruned_dict(tmp_path):
    """load_kegg_data should round-trip a synthetic nested kegg_data.json."""
    pruned = {
        "kos": {"K02338": {"name": "dnaN", "pathways": ["ko03030"]}},
        "pathways": {"ko03030": {"name": "DNA replication", "subcategory": "09124"}},
        "subcategories": {"09124": {"name": "Replication and repair", "category": "09120"}},
        "categories": {"09120": {"name": "Genetic Information Processing"}},
        "reactions": {},
        "compounds": {},
    }
    kegg_dir = tmp_path / "kegg"
    kegg_dir.mkdir()
    (kegg_dir / "kegg_data.json").write_text(json.dumps(pruned), encoding="utf-8")

    result = kegg_utils.load_kegg_data(tmp_path)

    assert result == pruned
    assert result["kos"]["K02338"]["name"] == "dnaN"
    assert result["pathways"]["ko03030"]["subcategory"] == "09124"


def test_load_kegg_data_missing_file_raises(tmp_path):
    """load_kegg_data should raise FileNotFoundError with a helpful message."""
    with pytest.raises(FileNotFoundError) as exc_info:
        kegg_utils.load_kegg_data(tmp_path)
    msg = str(exc_info.value)
    assert "kegg_data.json" in msg
    assert "prepare_data.sh" in msg
    assert "--steps 6" in msg


def test_load_kegg_data_corrupted_file_raises(tmp_path):
    """Corrupted (non-JSON) file produces a clear RuntimeError pointing at remediation."""
    (tmp_path / "kegg").mkdir()
    (tmp_path / "kegg" / "kegg_data.json").write_text("{not valid json")
    with pytest.raises(RuntimeError) as exc_info:
        kegg_utils.load_kegg_data(tmp_path)
    msg = str(exc_info.value)
    assert "corrupted" in msg
    assert "prepare_data.sh" in msg
    assert "--steps 6" in msg
    assert "--force" in msg
