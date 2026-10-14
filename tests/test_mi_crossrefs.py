import pytest

from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.crossrefs import (
    build_identity_fields, is_isotopologue_name, pick_chebi_id)

CHEBI = make_chebi(
    CHEBI_17632=term("nitrate", charge=-1, alt=["CHEBI:14654", "CHEBI:44487"], kegg=["C00244"],
                     inchikey="NHNBFGGVMKEFGY-UHFFFAOYSA-N"),
    CHEBI_48107=term("nitric acid", charge=0, kegg=["C00244"], inchikey="GRYLNZFGIOXLOG-UHFFFAOYSA-N"),
    CHEBI_25545=term("nitrates"),
    CHEBI_31000=term("succinic acid-d4", inchikey="KDYFGRWQOYBRFD-NMQOAUCRSA-N"),
    CHEBI_15741=term("succinic acid", charge=0, kegg=["C00042"], inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-N"),
    CHEBI_26806=term("succinate(2-)", charge=-2, inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-L"),
    CHEBI_30031=term("succinate(1-)", charge=-1, inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-M"))


def test_chebi_id_prefers_kegg_xref_charge0():
    assert pick_chebi_id("kegg.compound:C00244", cpd("Nitrate", chebi_id="14654"), CHEBI) == "48107"


def test_isotopologue_never_chebi_id():
    assert pick_chebi_id("kegg.compound:C00042", cpd("Succinate", chebi_id="31000"), CHEBI) == "15741"


def test_mnx_fallback_ranks_aliases_and_skips_isotopologues_and_classes():
    e = cpd("Succinate-ish")
    e["mnx_chebi_aliases"] = ["31000", "25545", "30031", "26806"]
    assert pick_chebi_id("kegg.compound:C09999", e, CHEBI) == "26806"


@pytest.mark.parametrize("name", ["succinic acid-d4", "UREA C 14", "calcium-43(2+)"])
def test_isotopologue_names(name):
    assert is_isotopologue_name(name)


@pytest.mark.parametrize("name", ["1H-indole-3-acetate", "Vitamin B12", "(2E)-hexenal", "(3S)-citramalate"])
def test_not_isotopologue_names(name):
    assert not is_isotopologue_name(name)


def test_identity_fields_nitrate():
    table = {"kegg.compound:C00244": cpd("Nitrate", chebi_id="14654", mnxm_id="MNXM732399",
                                         hmdb_id="HMDB02878", inchikey="GRYLNZFGIOXLOG-UHFFFAOYSA-N")}
    old = {"kegg.compound:C00244": dict(table["kegg.compound:C00244"]),
           "chebi:14654": cpd("nitrate", chebi_id="14654", mnxm_id="MNXM732398")}
    ident = IdentityResult(
        canonical_of={"kegg.compound:C00244": "kegg.compound:C00244", "chebi:14654": "kegg.compound:C00244"},
        members={"kegg.compound:C00244": ["chebi:14654", "kegg.compound:C00244"]},
        basis={"kegg.compound:C00244": ["chebi_conjugate"]},
        member_basis={"chebi:14654": "chebi_conjugate"})
    index = build_identity_fields(table, ident, {}, CHEBI, old)
    e = table["kegg.compound:C00244"]
    assert e["chebi_id"] == "48107" and e["chebi_ids"][0] == "48107"
    assert set(e["chebi_ids"]) == {"48107", "17632"}
    assert e["mnxm_ids"] == ["MNXM732399", "MNXM732398"]
    assert e["hmdb_id"] == "HMDB0002878" and e["identity_basis"] == ["chebi_conjugate"]
    assert {"kegg.compound:C00244", "chebi:14654", "chebi:17632", "chebi:44487", "chebi:48107",
            "mnx:MNXM732398", "mnx:MNXM732399", "hmdb:HMDB0002878",
            "inchikey:GRYLNZFGIOXLOG-UHFFFAOYSA-N"} <= set(e["all_identifiers"])
    assert all(":" in x for x in e["all_identifiers"])
    assert index["chebi:14654"] == "kegg.compound:C00244"


def test_rekeyed_node_has_basis_and_old_id():
    table = {"kegg.compound:C00244": cpd("Nitrate")}
    ident = IdentityResult({"kegg.compound:C00244": "kegg.compound:C00244"}, {}, {}, {})
    build_identity_fields(table, ident, {"chebi:17632": "kegg.compound:C00244"}, CHEBI,
                          {"kegg.compound:C00244": dict(table["kegg.compound:C00244"])})
    e = table["kegg.compound:C00244"]
    assert e["identity_basis"] == ["chebi_kegg_xref"] and "chebi:17632" in e["all_identifiers"]


def test_chebi_keyed_node_chebi_id_equals_local_part():
    table = {"chebi:17632": cpd("nitrate", chebi_id="44487")}
    ident = IdentityResult({"chebi:17632": "chebi:17632"}, {}, {}, {})
    build_identity_fields(table, ident, {}, CHEBI, {"chebi:17632": dict(table["chebi:17632"])})
    assert table["chebi:17632"]["chebi_id"] == "17632"
    assert table["chebi:17632"]["name"] == "nitrate"
    assert table["chebi:17632"].get("kegg_compound_ids") is None


@pytest.mark.parametrize("name", ["pyridin-2(1H)-one", "furan-2(5H)-one", "4(3H)-quinazolinone",
                                  "naphthalen-1(2H)-one", "25-hydroxyvitamin-D3"])
def test_more_not_isotopologue_names(name):
    assert not is_isotopologue_name(name)


def test_chebi_id_from_merged_member_pool():
    table = {"kegg.compound:C09998": cpd("x")}
    old = {"kegg.compound:C09998": dict(table["kegg.compound:C09998"]), "kegg.compound:C00244": cpd("y")}
    ident = IdentityResult({"kegg.compound:C09998": "kegg.compound:C09998",
                            "kegg.compound:C00244": "kegg.compound:C09998"},
                           {"kegg.compound:C09998": ["kegg.compound:C00244", "kegg.compound:C09998"]},
                           {"kegg.compound:C09998": ["name_key"]}, {})
    build_identity_fields(table, ident, {}, CHEBI, old)
    e = table["kegg.compound:C09998"]
    assert e["chebi_id"] == "48107" and e["chebi_ids"][0] == "48107"


def test_no_usable_term_both_none():
    table = {"kegg.compound:C09997": cpd("z")}
    ident = IdentityResult({"kegg.compound:C09997": "kegg.compound:C09997"}, {}, {}, {})
    build_identity_fields(table, ident, {}, CHEBI, {"kegg.compound:C09997": dict(table["kegg.compound:C09997"])})
    e = table["kegg.compound:C09997"]
    assert e["chebi_id"] is None and e["chebi_ids"] is None


def test_hmdb_from_member():
    table = {"kegg.compound:C09996": cpd("z")}
    old = {"kegg.compound:C09996": dict(table["kegg.compound:C09996"]),
           "chebi:25545": cpd("w", hmdb_id="HMDB02878")}
    ident = IdentityResult({"kegg.compound:C09996": "kegg.compound:C09996", "chebi:25545": "kegg.compound:C09996"},
                           {"kegg.compound:C09996": ["chebi:25545", "kegg.compound:C09996"]}, {}, {})
    build_identity_fields(table, ident, {}, CHEBI, old)
    e = table["kegg.compound:C09996"]
    assert e["hmdb_id"] == "HMDB0002878" and "hmdb:HMDB0002878" in e["all_identifiers"]


def test_secondary_mnx_alias_maps_to_primary():
    e = cpd("q")
    e["mnx_chebi_aliases"] = ["14654"]
    assert pick_chebi_id("kegg.compound:C09995", e, make_chebi(
        CHEBI_17632=term("nitrate", charge=-1, alt=["CHEBI:14654"], inchikey="NHNBFGGVMKEFGY-UHFFFAOYSA-N"))) == "17632"
