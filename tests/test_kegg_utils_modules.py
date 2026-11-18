"""KEGG module raw parsers and the step-6 module blocks (kegg_utils)."""
import json

import pytest

from multiomics_kg.utils import kegg_utils as ku

GET_TEXT = """\
ENTRY       M00031            Pathway   Module
NAME        Lysine biosynthesis, mediated by LysW, 2-aminoadipate => lysine
DEFINITION  K05826
            K05827 K05828
ORTHOLOGY   K05826  lysine biosynthesis protein LysW
CLASS       Pathway modules; Amino acid metabolism; Lysine metabolism
///
ENTRY       M00625            Signature Module
NAME        Methicillin resistance
DEFINITION  K02547 K02545
///
"""

INFO_TEXT = """\
kegg\tKEGG (Kyoto Encyclopedia of Genes and Genomes)
\tpathway        588  2026/11/18
\tbrite          205  2026/11/17
\tmodule         583  2026/11/09
\tko          28,592  2026/11/18
\tgenes   68,234,194  2026/11/18
\t  ag         4,498  2026/11/17
\tcompound    19,692  2026/11/18
\treaction    12,548  2026/11/18
\tenzyme       6,972  2026/11/16
\tKanehisa Laboratories
"""

BRITE = {"name": "ko00002", "children": [
    {"name": "Pathway modules", "children": [
        {"name": "Amino acid metabolism", "children": [
            {"name": "Lysine metabolism", "children": [
                {"name": "M00031  Lysine biosynthesis, mediated by LysW [PATH:map00300 map01100]"}]}]}]},
    {"name": "Signature modules", "children": [
        {"name": "Drug resistance", "children": [
            {"name": "Beta-Lactam resistance", "children": [
                {"name": "M00625  Methicillin resistance [PATH:map01501]"}]}]}]}]}


def test_parse_module_entries_multiline_definition_and_type():
    e = ku.parse_module_entries(GET_TEXT)
    assert set(e) == {"M00031", "M00625"}
    assert e["M00031"] == {"name": "Lysine biosynthesis, mediated by LysW, 2-aminoadipate => lysine",
                           "module_type": "pathway", "definition": "K05826\nK05827 K05828"}
    assert e["M00625"]["module_type"] == "signature"


def test_parse_module_entries_fails_loud_on_format_drift():
    with pytest.raises(ValueError):
        ku.parse_module_entries("SOMETHING   else entirely\n///\n")


def test_parse_module_classes():
    c = ku.parse_module_classes(BRITE)
    assert c["M00031"] == ["Pathway modules", "Amino acid metabolism", "Lysine metabolism"]
    assert c["M00625"][0] == "Signature modules"
    with pytest.raises(ValueError):
        ku.parse_module_classes({"name": "ko00002", "children": [{"name": "x", "children": []}]})


def test_parse_kegg_release():
    assert ku.parse_kegg_release(INFO_TEXT) == (
        "pathway 2026-11-18, brite 2026-11-17, module 2026-11-09, ko 2026-11-18, "
        "compound 2026-11-18, reaction 2026-11-18, enzyme 2026-11-16")
    with pytest.raises(ValueError):
        ku.parse_kegg_release("kegg\tKEGG\n\tpathway 588 2026/11/18\n")


def _write_module_raws(raw):
    raw.mkdir(parents=True, exist_ok=True)
    (raw / "list_module.txt").write_text("M00031\tLysine\nM00625\tMethicillin\n")
    (raw / "get_modules.txt").write_text(GET_TEXT)
    (raw / "link_pathway_module.txt").write_text(
        "md:M00031\tpath:map00300\nmd:M00031\tpath:map01100\nmd:M00625\tpath:map01501\n")
    (raw / "link_ko_enzyme.txt").write_text(
        "ec:2.3.1.1\tko:K05826\nec:2.3.1.1\tko:K99999\nec:9.9.9.9\tko:K99998\n")
    (raw / "br_ko00002.json").write_text(json.dumps(BRITE))
    (raw / "info_kegg.txt").write_text(INFO_TEXT)


def test_build_module_blocks(tmp_path):
    raw = tmp_path / "kegg" / "raw"
    _write_module_raws(raw)
    b = ku.build_module_blocks(raw)
    assert set(b) == {"modules", "module_ec_to_kos", "kegg_release"}
    m = b["modules"]["M00031"]
    assert m["class_path"] == ["Pathway modules", "Amino acid metabolism", "Lysine metabolism"]
    assert m["pathways"] == ["ko00300", "ko01100"]                  # map… → ko…, sorted
    assert m["definition"] == "K05826\nK05827 K05828"
    assert b["module_ec_to_kos"] == {"2.3.1.1": ["K05826"]}          # KOs outside every definition dropped
    assert b["kegg_release"].startswith("pathway 2026-11-18")


def test_build_module_blocks_fails_on_unparseable_definition(tmp_path):
    raw = tmp_path / "kegg" / "raw"
    _write_module_raws(raw)
    (raw / "get_modules.txt").write_text(GET_TEXT.replace("K02547 K02545", "K02547 ?? K02545"))
    with pytest.raises(ValueError):
        ku.build_module_blocks(raw)


def test_module_raw_files_named():
    assert set(ku.MODULE_RAW_FILES) == {"list_module.txt", "link_pathway_module.txt",
                                        "link_ko_enzyme.txt", "br_ko00002.json",
                                        "get_modules.txt", "info_kegg.txt"}


def test_download_module_raw_batches_get_by_ten(tmp_path, monkeypatch):
    ids = [f"M{n:05d}" for n in range(1, 24)]
    calls = []

    def fake_text(url):
        calls.append(url)
        if url.endswith("/list/module"):
            return "".join(f"{m}\tname\n" for m in ids)
        if "/get/" in url:
            return "".join(f"ENTRY       {m}            Pathway   Module\n///\n"
                           for m in url.rsplit("/", 1)[1].split("+"))
        return "x\n"
    monkeypatch.setattr(ku, "_fetch_text", fake_text)
    monkeypatch.setattr(ku, "_fetch_json", lambda url: {"name": "ko00002"})
    raw = tmp_path / "raw"
    ku.download_module_raw(raw, force=False)
    gets = [u for u in calls if "/get/M" in u]
    assert len(gets) == 3 and gets[0].endswith("+".join(ids[:10]))
    assert (raw / "get_modules.txt").read_text().count("///") == 23
    assert all((raw / f).exists() for f in ku.MODULE_RAW_FILES)
    calls.clear()
    ku.download_module_raw(raw, force=False)                          # cached: no network
    assert calls == []
