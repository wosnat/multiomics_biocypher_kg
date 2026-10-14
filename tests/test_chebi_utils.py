from pathlib import Path

from multiomics_kg.utils.chebi_utils import (
    build_chebi_cache, index_from_cache, parse_chebi_obo, ChebiIndex)

FIX = Path(__file__).parent / "fixtures" / "chebi_mini.obo"


def _index():
    release, terms = parse_chebi_obo(FIX)
    return ChebiIndex(terms, release)


def test_parse_release_and_fields():
    release, terms = parse_chebi_obo(FIX)
    assert release == "256"
    t = terms["CHEBI:48107"]
    assert t["kegg"] == ["C00244"] and t["charge"] == 0
    assert terms["CHEBI:17632"]["conjugates"] == ["CHEBI:48107"]
    assert terms["CHEBI:16977"]["enantiomers"] == ["CHEBI:15570"]
    assert terms["CHEBI:17522"]["generalized_formula"] == "C2H6O2(CH2O)n"


def test_secondary_maps_to_primary_and_is_recorded():
    ix = _index()
    assert ix.primary("chebi:10004") == "CHEBI:33229"
    assert ix.primary("10004") == "CHEBI:33229"
    assert ix.primary("CHEBI:999999") is None
    assert "CHEBI:33229" in ix.touched and "CHEBI:999999" in ix.missing
    assert ix.secondary_ids("CHEBI:33229") == ["CHEBI:10004"]


def test_kegg_terms_and_roles():
    ix = _index()
    assert ix.kegg_terms("C00244") == ["CHEBI:48107"]
    assert ix.is_role("CHEBI:38155") and ix.is_role("chebi:10004")
    assert not ix.is_role("CHEBI:16977")
    assert ix.ancestors("CHEBI:38155") >= {"CHEBI:26672", "CHEBI:24432", "CHEBI:50906"}


def test_class_and_carbohydrate():
    ix = _index()
    assert ix.is_class("CHEBI:17522")          # alditol: no inchikey, generic formula
    assert not ix.is_class("CHEBI:16977")
    assert ix.is_carbohydrate("CHEBI:15903")


def test_cache_coverage_and_round_trip():
    ix = _index()
    ix.primary("chebi:10004")                   # touched → kept with ancestors
    cache = build_chebi_cache(ix)
    assert cache["release"] == "256"
    kept = set(cache["terms"])
    assert {"CHEBI:48107", "CHEBI:17632"} <= kept      # KEGG xref + its conjugate
    assert {"CHEBI:33229", "CHEBI:24432", "CHEBI:50906"} <= kept
    assert cache["terms"]["CHEBI:33229"]["is_role"] is True
    again = index_from_cache(cache)
    assert again.primary("10004") == "CHEBI:33229" and again.release == "256"
