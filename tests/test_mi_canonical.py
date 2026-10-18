from tests.metabolite_identity_helpers import cpd, make_chebi, term
import pytest

from multiomics_kg.utils.metabolite_identity.canonical import (
    CanonicalPin, IdentityOverrides, MergeLink, apply_mnx_split, canonicalise, collect_merge_links,
    load_identity_overrides, rewrite_to_canonical, split_shared_mnx)

CHEBI = make_chebi(
    CHEBI_17632=term("nitrate", charge=-1, alt=["CHEBI:14654"], conj=["CHEBI:48107"]),
    CHEBI_48107=term("nitric acid", charge=0, kegg=["C00244"], conj=["CHEBI:17632"]),
    CHEBI_15377=term("water", charge=0, kegg=["C00001"], conj=["CHEBI:16234"]),
    CHEBI_16234=term("hydroxide", charge=-1, kegg=["C01328"], conj=["CHEBI:15377"]),
    CHEBI_17234=term("glucose", kegg=["C00293"]),
    CHEBI_28260=term("galactose", kegg=["C01582"]),
)
NO_OV = IdentityOverrides([], set())


def _canon(table, links, **kw):
    kw.setdefault("overrides", NO_OV)
    kw.setdefault("reaction_refs", {})
    kw.setdefault("table_pins", [])
    return canonicalise(table, links, chebi=CHEBI, **kw)


def test_nitrate_merges_by_conjugate_onto_kegg():
    table = {"kegg.compound:C00244": cpd("Nitrate"),
             "chebi:14654": cpd("nitrate", chebi_id="14654", ev=("transport",))}
    links = collect_merge_links(table, CHEBI)
    assert MergeLink("chebi:14654", "kegg.compound:C00244", "chebi_conjugate") in links
    res = _canon(table, links)
    assert res.canonical_of["chebi:14654"] == "kegg.compound:C00244"
    assert res.basis["kegg.compound:C00244"] == ["chebi_conjugate"]


def test_mnx_same_entity_winner_by_reaction_count():
    table = {"kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    res = _canon(table, collect_merge_links(table, CHEBI),
                 reaction_refs={"kegg.compound:C00288": 40, "kegg.compound:C01353": 3})
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C00288"
    assert res.basis["kegg.compound:C00288"] == ["mnx_same_entity"]


def test_chebi_only_kegg_link_rejected_water_hydroxide():
    table = {"kegg.compound:C00001": cpd("H2O"), "kegg.compound:C01328": cpd("HO-")}
    res = _canon(table, collect_merge_links(table, CHEBI))
    assert res.canonical_of["kegg.compound:C01328"] == "kegg.compound:C01328"
    assert [(l.a, l.b) for l in res.rejected] == [("kegg.compound:C00001", "kegg.compound:C01328")]


def test_manual_merge_overrides_guard_and_keep_separate_blocks_mnx():
    table = {"kegg.compound:C00001": cpd("H2O", mnxm_id="M1"),
             "kegg.compound:C01328": cpd("HO-"),
             "kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    ov = IdentityOverrides(
        merges=[MergeLink("kegg.compound:C00001", "kegg.compound:C01328", "manual")],
        keep_separate={frozenset({"kegg.compound:C00288", "kegg.compound:C01353"})})
    res = _canon(table, collect_merge_links(table, CHEBI), overrides=ov,
                 reaction_refs={"kegg.compound:C00001": 50, "kegg.compound:C01328": 4})
    assert res.canonical_of["kegg.compound:C01328"] == "kegg.compound:C00001"
    assert "manual" in res.basis["kegg.compound:C00001"]
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C01353"


def test_stale_identity_override_reported_and_skipped():
    ov = IdentityOverrides(merges=[MergeLink("kegg.compound:C99998", "kegg.compound:C99999", "manual")],
                           keep_separate=set())
    res = _canon({"kegg.compound:C00001": cpd("H2O")}, [], overrides=ov)
    assert res.stale_overrides and res.stale_overrides[0]["kind"] == "merge"


def test_assay_collision_merge_refused():
    table = {"kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    pins = [("p|t", "bicarbonate", "kegg.compound:C00288"),
            ("p|t", "carbonic acid", "kegg.compound:C01353")]
    res = _canon(table, collect_merge_links(table, CHEBI), table_pins=pins)
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C01353"
    assert res.refused_collisions


def test_glucose_galactose_stay_apart():
    table = {"kegg.compound:C00293": cpd("Glucose", inchikey="WQZGKKKJIJFFOK-UHFFFAOYSA-N"),
             "kegg.compound:C01582": cpd("Galactose", inchikey="WQZGKKKJIJFFOK-UHFFFAOYSA-N")}
    res = _canon(table, collect_merge_links(table, CHEBI))
    assert not res.members


def test_chebi_group_prefers_charge0_then_lowest():
    chebi = make_chebi(CHEBI_100=term("x acid", charge=0, conj=["CHEBI:50"]),
                       CHEBI_50=term("x anion", charge=-1, conj=["CHEBI:100"]))
    table = {"chebi:100": cpd("x acid"), "chebi:50": cpd("x anion")}
    res = canonicalise(table, collect_merge_links(table, chebi), overrides=NO_OV,
                       reaction_refs={}, chebi=chebi, table_pins=[])
    assert res.canonical_of["chebi:50"] == "chebi:100"


def test_rewrite_unions_lists_and_dedups_reactions():
    table = {"kegg.compound:C00288": cpd("HCO3-", pathways=["ko00910"], ev=("metabolism",)),
             "kegg.compound:C01353": cpd("Carbonic acid", pathways=["ko00720"], ev=("transport",))}
    canon = {"kegg.compound:C00288": "kegg.compound:C00288",
             "kegg.compound:C01353": "kegg.compound:C00288"}
    rxns = {"R1": {"compounds": ["C00288", "C01353", "C00001"]}}
    subs = {"2.A.1": ["kegg.compound:C01353", "kegg.compound:C00288"]}
    new = rewrite_to_canonical(table, rxns, subs, canon)
    assert list(new) == ["kegg.compound:C00288"]
    assert new["kegg.compound:C00288"]["pathways"] == ["ko00720", "ko00910"]
    assert new["kegg.compound:C00288"]["evidence_sources"] == ["metabolism", "transport"]
    assert rxns["R1"]["compounds"] == ["C00001", "C00288"]
    assert subs["2.A.1"] == ["kegg.compound:C00288"]


def test_load_overrides_contradiction_fails(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("merge:\n  - {a: kegg.compound:C1, b: kegg.compound:C2, reason: x}\n"
                 "keep_separate:\n  - {a: kegg.compound:C2, b: kegg.compound:C1, reason: y}\n")
    with pytest.raises(ValueError):
        load_identity_overrides(p)


# ── explorer ask 2026-10-15: MNX stereo guard, canonical margin, pins, MNX field split ─────────

STEREO = make_chebi(
    # 3-phosphoglycerate: stereo-unspecified (C00597) vs D (C00197)
    CHEBI_17050=term("3-phosphoglyceric acid", kegg=["C00597"], inchikey="OSJPPGNTCRNQQC-UHFFFAOYSA-N"),
    CHEBI_17794=term("3-phospho-D-glyceric acid", kegg=["C00197"], is_a=["CHEBI:17050"],
                     inchikey="OSJPPGNTCRNQQC-UWTATZPHSA-N"),
    # fructose 6-phosphate: generic (C00085) vs beta anomer (C05345)
    CHEBI_61553=term("D-fructofuranose 6-phosphate", kegg=["C00085"], inchikey="BGWGXPAPYGQALX-VRPWFDPXSA-N"),
    CHEBI_16084=term("beta-D-fructofuranose 6-phosphate", kegg=["C05345"], is_a=["CHEBI:61553"],
                     inchikey="BGWGXPAPYGQALX-ARQDHWQXSA-N"),
    # caffeate: KEGG lists generic + trans as two ids of one compound (geometric: exempt)
    CHEBI_16433=term("trans-caffeic acid", kegg=["C01197", "C01481"], inchikey="QAIPRVGONGVQAS-DUXPYHPUSA-N"),
    CHEBI_36281=term("caffeic acid", kegg=["C01481"], inchikey="QAIPRVGONGVQAS-UHFFFAOYSA-N"),
    # conjugate pair, same stereo layer
    CHEBI_17996=term("chloride", charge=-1, kegg=["C00698"], inchikey="VEXZGXHMUGYJMC-UHFFFAOYSA-M"),
    CHEBI_17883=term("hydrogen chloride", charge=0, kegg=["C01327"], inchikey="VEXZGXHMUGYJMC-UHFFFAOYSA-N"),
)


def _mnx_pair(a, b, mnx="MNXM1", ik=None):
    return {f"kegg.compound:{a}": cpd(a, mnxm_id=mnx, inchikey=ik, smiles="X", hmdb_id="HMDB0000001"),
            f"kegg.compound:{b}": cpd(b, mnxm_id=mnx, inchikey=ik, smiles="X", hmdb_id="HMDB0000001")}


def _canon_s(table, **kw):
    kw.setdefault("overrides", NO_OV)
    kw.setdefault("table_pins", [])
    return canonicalise(table, collect_merge_links(table, STEREO), chebi=STEREO, **kw)


@pytest.mark.parametrize("a,b", [("C00197", "C00597"), ("C05345", "C00085")])
def test_mnx_link_between_stereo_forms_refused(a, b):
    res = _canon_s(_mnx_pair(a, b), reaction_refs={f"kegg.compound:{a}": 10, f"kegg.compound:{b}": 9})
    assert not res.members
    assert res.stereo_refused and res.stereo_refused[0]["link"] == sorted(
        [f"kegg.compound:{a}", f"kegg.compound:{b}"])


@pytest.mark.parametrize("a,b", [("C01197", "C01481"), ("C00698", "C01327")])
def test_mnx_link_geometric_or_conjugate_pair_kept(a, b):
    res = _canon_s(_mnx_pair(a, b), reaction_refs={f"kegg.compound:{a}": 10})
    assert res.canonical_of[f"kegg.compound:{b}"] == f"kegg.compound:{a}"
    assert not res.stereo_refused


def test_manual_merge_of_stereo_forms_exempt_from_guard():
    ov = IdentityOverrides([MergeLink("kegg.compound:C00085", "kegg.compound:C05345", "manual")], set())
    res = _canon_s(_mnx_pair("C05345", "C00085"), overrides=ov,
                   reaction_refs={"kegg.compound:C00085": 30})
    assert res.canonical_of["kegg.compound:C05345"] == "kegg.compound:C00085"


def test_canonical_margin_below_2x_is_unresolved():
    res = _canon_s(_mnx_pair("C00698", "C01327"),
                   reaction_refs={"kegg.compound:C00698": 10, "kegg.compound:C01327": 9})
    assert res.canonical_unresolved == [{
        "members": ["kegg.compound:C00698", "kegg.compound:C01327"],
        "gene_reachable_reactions": {"kegg.compound:C00698": 10, "kegg.compound:C01327": 9},
        "all_kegg_reactions": {"kegg.compound:C00698": 0, "kegg.compound:C01327": 0}}]


def test_canonical_margin_one_vs_zero_passes():
    res = _canon_s(_mnx_pair("C00698", "C01327"), reaction_refs={"kegg.compound:C01327": 1})
    assert res.canonical_of["kegg.compound:C00698"] == "kegg.compound:C01327"
    assert not res.canonical_unresolved


def test_canonical_fallback_all_kegg_refs_only_when_zero_vs_zero():
    res = _canon_s(_mnx_pair("C00698", "C01327"), reaction_refs={},
                   all_reaction_refs={"kegg.compound:C00698": 47, "kegg.compound:C01327": 6})
    assert res.canonical_of["kegg.compound:C01327"] == "kegg.compound:C00698"
    # with any gene-reachable signal, the fallback never runs
    res = _canon_s(_mnx_pair("C00698", "C01327"),
                   reaction_refs={"kegg.compound:C00698": 3, "kegg.compound:C01327": 2},
                   all_reaction_refs={"kegg.compound:C00698": 47, "kegg.compound:C01327": 6})
    assert res.canonical_unresolved


def test_canonical_pin_overrides_margin():
    ov = IdentityOverrides([], set(), [CanonicalPin("kegg.compound:C00698",
                                                    ("kegg.compound:C00698", "kegg.compound:C01327"))])
    res = _canon_s(_mnx_pair("C00698", "C01327"), overrides=ov,
                   reaction_refs={"kegg.compound:C00698": 1, "kegg.compound:C01327": 13})
    assert res.canonical_of["kegg.compound:C01327"] == "kegg.compound:C00698"
    assert res.canonical_pinned == ["kegg.compound:C00698"]


def test_canonical_pin_on_group_that_no_longer_forms_is_stale():
    ov = IdentityOverrides([], set(), [CanonicalPin("kegg.compound:C00085",
                                                    ("kegg.compound:C00085", "kegg.compound:C05345"))])
    res = _canon_s(_mnx_pair("C05345", "C00085"), overrides=ov, reaction_refs={})
    assert res.stale_overrides == [{"kind": "canonical", "id": "kegg.compound:C00085",
                                    "group": ["kegg.compound:C00085", "kegg.compound:C05345"]}]


def test_load_canonical_pin_must_be_in_its_group(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("canonical:\n  - {id: kegg.compound:C1, group: [kegg.compound:C2, kegg.compound:C3]}\n")
    with pytest.raises(ValueError):
        load_identity_overrides(p)


def test_split_shared_mnx_keeps_fields_on_structure_match():
    table = _mnx_pair("C05345", "C00085", ik="BGWGXPAPYGQALX-ARQDHWQXSA-L")
    res = _canon_s(table, reaction_refs={"kegg.compound:C00085": 99})
    strip, report = split_shared_mnx(table, res.canonical_of, STEREO,
                                     {"kegg.compound:C00085": 99})
    # MNX's InChIKey is the beta anomer's, so beta keeps the MNX fields despite fewer reactions
    assert strip == {"kegg.compound:C00085": "BGWGXPAPYGQALX-VRPWFDPXSA-N"}
    assert report[0]["kept_on"] == "kegg.compound:C05345"
    assert report[0]["basis"] == "structure_match"
    apply_mnx_split(table, strip)
    e = table["kegg.compound:C00085"]
    assert (e["mnxm_id"], e["smiles"], e["hmdb_id"]) == (None, None, None)
    assert e["inchikey"] == "BGWGXPAPYGQALX-VRPWFDPXSA-N"
    assert table["kegg.compound:C05345"]["mnxm_id"] == "MNXM1"
