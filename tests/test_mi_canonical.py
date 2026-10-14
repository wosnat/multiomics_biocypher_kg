from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.canonical import (
    IdentityOverrides, MergeLink, canonicalise, collect_merge_links, load_identity_overrides,
    rewrite_to_canonical)

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
    kw.setdefault("pathway_refs", {})
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
    res = _canon(table, collect_merge_links(table, CHEBI), overrides=ov)
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
                       reaction_refs={}, pathway_refs={}, chebi=chebi, table_pins=[])
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
    import pytest
    with pytest.raises(ValueError):
        load_identity_overrides(p)
