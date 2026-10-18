import pytest

from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult, load_identity_overrides
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
                                  "naphthalen-1(2H)-one", "25-hydroxyvitamin-D3",
                                  "(25S)-cholest-5-ene", "(12S)-HETE", "(20S)-x", "Vitamin C 12"])
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


SHARED = make_chebi(
    CHEBI_28061=term("alpha-D-galactose", charge=0, alt=["CHEBI:10231"], kegg=["C00984", "C01582"],
                     inchikey="WQZGKKKJIJFFOK-PHYPRBDBSA-N"),
    CHEBI_28260=term("galactose", charge=0, kegg=["C01582"], inchikey="GZCGUPFRVQAUEE-UHFFFAOYSA-N"),
    CHEBI_4139=term("D-galactopyranose", charge=0, kegg=["C00124"], inchikey="WQZGKKKJIJFFOK-SVZMEOIVSA-N"))


def _two_kegg_table():
    table = {"kegg.compound:C00984": cpd("alpha-D-Galactose"), "kegg.compound:C01582": cpd("Galactose")}
    ident = IdentityResult(canonical_of={n: n for n in table}, members={}, basis={}, member_basis={})
    return table, ident


def test_shared_chebi_terms_reports_term_on_two_separate_kegg_nodes():
    from multiomics_kg.utils.metabolite_identity.crossrefs import shared_chebi_terms
    table, ident = _two_kegg_table()
    # C00984 xrefs only alpha-D-galactose; C01582 also xrefs galactose: C00984 owns the term
    assert shared_chebi_terms(table, ident, SHARED) == {"CHEBI:28061": ["kegg.compound:C01582"]}


def test_shared_term_not_shared_when_both_kegg_ids_on_one_canonical_node():
    from multiomics_kg.utils.metabolite_identity.crossrefs import shared_chebi_terms
    table = {"kegg.compound:C00984": cpd("alpha-D-Galactose")}
    ident = IdentityResult(
        canonical_of={"kegg.compound:C00984": "kegg.compound:C00984", "kegg.compound:C01582": "kegg.compound:C00984"},
        members={"kegg.compound:C00984": ["kegg.compound:C00984", "kegg.compound:C01582"]},
        basis={}, member_basis={})
    assert shared_chebi_terms(table, ident, SHARED) == {}


def test_shared_term_goes_to_its_single_owner_only():
    from multiomics_kg.utils.metabolite_identity.crossrefs import shared_chebi_terms
    table, ident = _two_kegg_table()
    shared = shared_chebi_terms(table, ident, SHARED)
    old = {n: dict(e) for n, e in table.items()}
    index = build_identity_fields(table, ident, {}, SHARED, old, shared=shared)
    a, b = table["kegg.compound:C00984"], table["kegg.compound:C01582"]
    assert a["chebi_id"] == "28061" and a["chebi_ids"] == ["28061"]
    assert b["chebi_id"] == "28260" and b["chebi_ids"] == ["28260"]
    assert "chebi:10231" in a["all_identifiers"]
    assert "chebi:28061" not in b["all_identifiers"] and "chebi:10231" not in b["all_identifiers"]
    assert index["chebi:28061"] == "kegg.compound:C00984"
    assert pick_chebi_id("kegg.compound:C01582", cpd("Galactose"), SHARED, shared) == "28260"
    assert pick_chebi_id("kegg.compound:C00984", cpd("x", chebi_id="10231"), SHARED, shared) == "28061"


TIED = make_chebi(
    CHEBI_17594=term("hydroquinone", kegg=["C00530", "C15603"], inchikey="QIGBRXMKCJKVMJ-UHFFFAOYSA-N"),
    CHEBI_17169=term("benzaldehyde", kegg=["C00193", "C00261"], inchikey="HUMNYLRZRPPJDN-UHFFFAOYSA-N"),
    CHEBI_33855=term("arenecarbaldehyde", kegg=["C00193"]))


def test_shared_term_without_single_owner_withheld_from_all():
    from multiomics_kg.utils.metabolite_identity.crossrefs import shared_chebi_terms
    table = {f"kegg.compound:{k}": cpd(k) for k in ("C00530", "C15603", "C00193", "C00261")}
    ident = IdentityResult(canonical_of={n: n for n in table}, members={}, basis={}, member_basis={})
    shared = shared_chebi_terms(table, ident, TIED)
    # both hydroquinone ids xref only hydroquinone: a tie, so neither publishes it
    assert shared["CHEBI:17594"] == ["kegg.compound:C00530", "kegg.compound:C15603"]
    # the class node C00193 also xrefs a class term, so benzaldehyde goes to C00261
    assert shared["CHEBI:17169"] == ["kegg.compound:C00193"]


@pytest.mark.parametrize("name", ["[13C]glucose", "(15N)urea", "[2H]water", "(18O)water", "N 15 ammonium"])
def test_more_isotopologue_names(name):
    assert is_isotopologue_name(name)


# ── one chebi_id pick for stages G and H (classes only as a last resort) ──────

CLS = make_chebi(
    CHEBI_15889=term("sterol", kegg=["C00370"]),                       # class (no InChIKey)
    CHEBI_1=term("thing", charge=0, kegg=["C09001"], inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N"),
    CHEBI_2=term("things", kegg=["C09001"]))                           # class on the same KEGG id


def test_class_term_rekeyed_onto_concrete_kegg_node_is_generic_and_scrubbed():
    from multiomics_kg.utils.metabolite_identity.crossrefs import pick_chebi_id, node_olds
    from multiomics_kg.utils.metabolite_identity.structure import apply_structure
    node = "kegg.compound:C00370"
    table = {node: cpd("Sterol", formula="C27H46O", mass=386.7, inchikey="HVYWMOMLDIMFJA-DPAQBDIFSA-N",
                       smiles="CC")}
    ident = IdentityResult(canonical_of={node: node}, members={}, basis={}, member_basis={})
    rekeys = {"chebi:15889": node}
    old = {n: dict(e) for n, e in table.items()}
    olds = node_olds(node, ident, rekeys)
    cid = pick_chebi_id(node, table[node], CLS, olds=olds, old_entries=old)
    assert cid == "15889"
    apply_structure(table, CLS, {node: f"CHEBI:{cid}"})
    e = table[node]
    assert e["structure_kind"] == "generic"
    assert e["inchikey"] is None and e["smiles"] is None and e["formula"] is None and e["mass"] is None
    assert e["elements"] == ["C", "H", "O"]
    build_identity_fields(table, ident, rekeys, CLS, old)
    assert e["chebi_id"] == "15889" and e["chebi_ids"] == ["15889"]


def test_entity_term_beats_class_term_on_one_kegg_node():
    from multiomics_kg.utils.metabolite_identity.crossrefs import pick_chebi_id
    node = "kegg.compound:C09001"
    table = {node: cpd("Thing", formula="C2H6", inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N")}
    ident = IdentityResult(canonical_of={node: node}, members={}, basis={}, member_basis={})
    rekeys = {"chebi:2": node}                                         # the class re-keyed in too
    assert pick_chebi_id(node, table[node], CLS) == "1"
    build_identity_fields(table, ident, rekeys, CLS, {n: dict(e) for n, e in table.items()})
    assert table[node]["chebi_id"] == "1" and table[node]["chebi_ids"] == ["1"]


def test_rekeyed_id_follows_its_target_into_the_canonical_node():
    from multiomics_kg.utils.metabolite_identity.crossrefs import node_olds
    a, b = "kegg.compound:C00193", "kegg.compound:C00261"
    table = {b: cpd("Benzaldehyde-ish", formula="C7H6O")}
    ident = IdentityResult(canonical_of={a: b, b: b}, members={b: [a, b]},
                           basis={b: ["mnx_same_entity"]}, member_basis={a: "mnx_same_entity"})
    rekeys = {"chebi:33855": a}
    assert "chebi:33855" in node_olds(b, ident, rekeys)
    index = build_identity_fields(table, ident, rekeys, make_chebi(),
                                  {a: cpd("A"), b: dict(table[b])})
    assert "chebi:33855" in table[b]["all_identifiers"] and index["chebi:33855"] == b
    assert table[b]["identity_basis"] == ["chebi_kegg_xref", "mnx_same_entity"]


# Final review items 2 + 7: isotope variants (structural) and ChEBI role terms never reach
# chebi_id / chebi_ids / all_identifiers.
ISO_ROLE = make_chebi(
    CHEBI_24431=term("chemical entity"),
    CHEBI_33250=term("atom", is_a=["CHEBI:24431"]),
    CHEBI_26833=term("sulfur atom", charge=0, kegg=["C00087"], is_a=["CHEBI:33250"],
                     inchikey="NINIDFKCEFEMDL-UHFFFAOYSA-N"),
    CHEBI_37979=term("sulfur-32 atom", charge=0, kegg=["C00087"], is_a=["CHEBI:26833"],
                     inchikey="NINIDFKCEFEMDL-AKLPVKDBSA-N"),
    CHEBI_139358=term("isotopically modified compound", is_a=["CHEBI:24431"]),
    CHEBI_32285=term("((13)C)urea", charge=0, is_a=["CHEBI:139358"], inchikey="XSQUKJJJFZCRTK-OUBTZVSYSA-N"),
    CHEBI_16199=term("urea", charge=0, kegg=["C00086"], is_a=["CHEBI:24431"],
                     inchikey="XSQUKJJJFZCRTK-UHFFFAOYSA-N"),
    CHEBI_50906=term("role"),
    CHEBI_51086=term("chemical role", is_a=["CHEBI:50906"]),
    CHEBI_15339=term("acceptor", kegg=["C00028"], is_a=["CHEBI:51086"]))


def _fields(node, entry, olds=None):
    table = {node: entry}
    members = {node: sorted(olds)} if olds else {}
    ident = IdentityResult({n: node for n in (olds or [node])} | {node: node}, members, {}, {})
    build_identity_fields(table, ident, {}, ISO_ROLE, {node: dict(entry), **{o: cpd(o) for o in olds or []}})
    return table[node]


def test_atom_isotope_xref_not_published():
    e = _fields("kegg.compound:C00087", cpd("Sulfur"))
    assert e["chebi_id"] == "26833" and e["chebi_ids"] == ["26833"]
    assert "chebi:37979" not in e["all_identifiers"]


def test_isotopically_modified_alias_and_member_not_published():
    e = cpd("Urea")
    e["mnx_chebi_aliases"] = ["32285"]
    out = _fields("kegg.compound:C00086", e)
    assert out["chebi_ids"] == ["16199"] and "chebi:32285" not in out["all_identifiers"]


def test_role_term_never_chebi_id():
    e = _fields("kegg.compound:C00028", cpd("Acceptor"))
    assert e["chebi_id"] is None and e["chebi_ids"] is None
    assert not any(x.startswith("chebi:") for x in e["all_identifiers"])


# ── xref: overrides (explorer ask 5, 2026-10-15) ───────────────────────────────────────────────

def _xref_table():
    def node(n, hmdb=None, chebi=()):
        ids = {f"kegg.compound:{n}"} | {f"chebi:{c}" for c in chebi}
        if hmdb:
            ids.add(f"hmdb:{hmdb}")
        return {"hmdb_id": hmdb, "chebi_id": chebi[0] if chebi else None,
                "chebi_ids": list(chebi) or None, "all_identifiers": sorted(ids)}
    table = {"kegg.compound:C00221": node("C00221", "HMDB0000122", ["15903"]),
             "kegg.compound:C00031": node("C00031", "HMDB0304632", ["4167"])}
    index = {i: n for n, e in table.items() for i in e["all_identifiers"]}
    return table, index


GLC = make_chebi(CHEBI_17634=term("D-glucose", alt=["CHEBI:12965", "CHEBI:20999"]),
                 CHEBI_4167=term("D-glucopyranose", kegg=["C00031"]),
                 CHEBI_15903=term("beta-D-glucose", kegg=["C00221"]))


def test_xref_override_moves_hmdb_to_generic():
    from multiomics_kg.utils.metabolite_identity.canonical import XrefOverride
    from multiomics_kg.utils.metabolite_identity.crossrefs import apply_xref_overrides
    table, index = _xref_table()
    applied, stale = apply_xref_overrides(table, index, [XrefOverride(
        "hmdb:HMDB0000122", "kegg.compound:C00031", "kegg.compound:C00221")], GLC)
    assert len(applied) == 1 and not stale
    b, g = table["kegg.compound:C00221"], table["kegg.compound:C00031"]
    assert b["hmdb_id"] is None and "hmdb:HMDB0000122" not in b["all_identifiers"]
    assert g["hmdb_id"] == "HMDB0000122"
    assert {"hmdb:HMDB0000122", "hmdb:HMDB0304632"} <= set(g["all_identifiers"])
    assert index["hmdb:HMDB0000122"] == "kegg.compound:C00031"


def test_xref_override_adds_chebi_with_secondary_ids():
    from multiomics_kg.utils.metabolite_identity.canonical import XrefOverride
    from multiomics_kg.utils.metabolite_identity.crossrefs import apply_xref_overrides
    table, index = _xref_table()
    applied, stale = apply_xref_overrides(
        table, index, [XrefOverride("chebi:17634", "kegg.compound:C00031")], GLC)
    g = table["kegg.compound:C00031"]
    assert applied and not stale
    assert g["chebi_id"] == "4167" and g["chebi_ids"] == ["4167", "17634"]
    assert {"chebi:17634", "chebi:12965", "chebi:20999"} <= set(g["all_identifiers"])
    assert index["chebi:20999"] == "kegg.compound:C00031"


def test_xref_override_stale_when_identifier_elsewhere():
    from multiomics_kg.utils.metabolite_identity.canonical import XrefOverride
    from multiomics_kg.utils.metabolite_identity.crossrefs import apply_xref_overrides
    table, index = _xref_table()
    before = {n: dict(e) for n, e in table.items()}
    # wrong `from`, an add of an id another node holds, and an unwritten target
    applied, stale = apply_xref_overrides(table, index, [
        XrefOverride("hmdb:HMDB0000122", "kegg.compound:C00031", "kegg.compound:C99999"),
        XrefOverride("chebi:15903", "kegg.compound:C00031"),
        XrefOverride("hmdb:HMDB0000122", "kegg.compound:C99998", "kegg.compound:C00221")], GLC)
    assert not applied and len(stale) == 3
    assert table == before


def test_load_xref_override_validates(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("xref:\n  - {id: 'mnx:MNXM1', to: kegg.compound:C1, reason: x}\n")
    with pytest.raises(ValueError):
        load_identity_overrides(p)
    p.write_text("xref:\n  - {id: 'hmdb:HMDB0000122', to: kegg.compound:C00031, "
                 "from: kegg.compound:C00221, reason: x}\n")
    ov = load_identity_overrides(p)
    assert ov.xrefs[0].frm == "kegg.compound:C00221"


# ── narrowed D (explorer ask 8): a chebi: node's structure / ChEBI ids are its own term's ─────

PIP = make_chebi(
    CHEBI_17964=term("pipecolic acid", inchikey="HXEACLLIILLPRG-UHFFFAOYSA-N"),
    CHEBI_41582=term("D-pipecolic acid", is_a=["CHEBI:17964"], alt=["CHEBI:41581"],
                     inchikey="HXEACLLIILLPRG-RXMQYKEDSA-N"),
    CHEBI_29223=term("bromate", inchikey="SXDBWCPKPHAZSM-UHFFFAOYSA-M", conj=["CHEBI:29222"]),
    CHEBI_29222=term("bromic acid", inchikey="SXDBWCPKPHAZSM-UHFFFAOYSA-N", conj=["CHEBI:29223"]),
    CHEBI_30867=term("dihydroorotate", inchikey="UFIVEPVSAGBUSI-UHFFFAOYSA-M", conj=["CHEBI:30865"]),
    CHEBI_30865=term("dihydroorotic acid", inchikey="UFIVEPVSAGBUSI-UHFFFAOYSA-N", conj=["CHEBI:30867"]),
    CHEBI_30866=term("(R)-dihydroorotic acid", is_a=["CHEBI:30865"], inchikey="UFIVEPVSAGBUSI-UWTATZPHSA-N"))


def _chebi_node(n, ik, chebi_ids, extra=()):
    ids = {n, f"inchikey:{ik}"} | {f"chebi:{c}" for c in chebi_ids} | set(extra)
    return {"inchikey": ik, "smiles": "X", "chebi_id": n.split(":")[1], "chebi_ids": list(chebi_ids),
            "all_identifiers": sorted(ids), "formula": "F", "structure_kind": "defined"}


def test_chebi_node_structure_moves_stereo_mismatch_keeps_variants():
    from multiomics_kg.utils.metabolite_identity.crossrefs import chebi_node_structure
    table = {"chebi:17964": _chebi_node("chebi:17964", "HXEACLLIILLPRG-RXMQYKEDSA-M", ["17964", "41582"],
                                        ["chebi:41581"]),
             "chebi:29223": _chebi_node("chebi:29223", "SXDBWCPKPHAZSM-OTHERREPSA-N", ["29223", "29222"]),
             "chebi:30867": _chebi_node("chebi:30867", "UFIVEPVSAGBUSI-UWTATZPHSA-N", ["30867", "30865", "30866"])}
    index = {i: n for n, e in table.items() for i in e["all_identifiers"]}
    rows = chebi_node_structure(table, index, PIP)
    p, b, d = table["chebi:17964"], table["chebi:29223"], table["chebi:30867"]
    # stereo mismatch: id, its secondary and the MNX InChIKey move out; own structure comes in
    assert p["chebi_ids"] == ["17964"] and p["inchikey"] == "HXEACLLIILLPRG-UHFFFAOYSA-N" and p["smiles"] is None
    assert p["removed_stereo_xrefs"] == ["chebi:41581", "chebi:41582", "inchikey:HXEACLLIILLPRG-RXMQYKEDSA-M"]
    assert "chebi:41582" not in p["all_identifiers"] and "chebi:41582" not in index
    assert index["inchikey:HXEACLLIILLPRG-UHFFFAOYSA-N"] == "chebi:17964"
    # representation variant: own InChIKey published, the MNX one stays resolvable, nothing removed
    assert b["inchikey"] == "SXDBWCPKPHAZSM-UHFFFAOYSA-M" and "removed_stereo_xrefs" not in b
    assert "inchikey:SXDBWCPKPHAZSM-OTHERREPSA-N" in b["all_identifiers"] and b["chebi_ids"] == ["29223", "29222"]
    # (R)-acid sits below the anion's conjugate acid: removed; the acid itself stays
    assert d["chebi_ids"] == ["30867", "30865"]
    assert set(d["removed_stereo_xrefs"]) == {"chebi:30866", "inchikey:UFIVEPVSAGBUSI-UWTATZPHSA-N"}
    assert len(rows) == 3
    for e in (p, b, d):
        assert e["formula"] == "F" and e["structure_kind"] == "defined"
