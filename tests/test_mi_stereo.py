from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.stereo import (
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)

ALA, GLC = dict(gformula="C3H7NO2", charge=0), dict(gformula="C6H12O6", charge=0)
CHEBI = make_chebi(
    CHEBI_16449=term("alanine", **ALA),
    CHEBI_16977=term("L-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:15570"], **ALA),
    CHEBI_15570=term("D-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:16977"], **ALA),
    CHEBI_16414=term("L-valine", enan=["CHEBI:27266"]),
    CHEBI_27266=term("D-valine", enan=["CHEBI:16414"]),
    CHEBI_16646=term("carbohydrate"),
    CHEBI_4167=term("D-glucopyranose", is_a=["CHEBI:16646"], **GLC),
    CHEBI_15903=term("beta-D-glucose", is_a=["CHEBI:4167"], **GLC),
    CHEBI_17522=term("alditol"),
    CHEBI_17754=term("glycerol", is_a=["CHEBI:17522"]),
    CHEBI_100=term("citronellol"),
    CHEBI_101=term("beta-citronellol", enan=["CHEBI:102"]),
    CHEBI_102=term("(-)-citronellol", enan=["CHEBI:101"]))

def _d(name, ik):
    e = cpd(name, inchikey=ik)
    e["structure_kind"] = "defined"
    return e

T = {"kegg.compound:C01401": _d("Alanine", "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"),
     "kegg.compound:C00041": _d("L-Alanine", "QNAYBMKLOCPYGJ-REOHCLBHSA-N"),
     "kegg.compound:C00133": _d("D-Alanine", "QNAYBMKLOCPYGJ-UWTATZPHSA-N"),
     "kegg.compound:C00183": _d("L-Valine", "KZSNJWFQEVHDMF-BYPYZUCNSA-N"),
     "kegg.compound:C06417": _d("D-Valine", "KZSNJWFQEVHDMF-SCSAIBSYSA-N"),
     "kegg.compound:C00031": _d("D-Glucose", "WQZGKKKJIJFFOK-GASJEMHNSA-N"),
     "kegg.compound:C00221": _d("beta-D-Glucose", "WQZGKKKJIJFFOK-VFUOTHLCSA-N"),
     "chebi:17522": dict(cpd("alditol"), structure_kind="generic"),
     "kegg.compound:C00116": _d("Glycerol", "PEDCQBHIVMGVHV-UHFFFAOYSA-N"),
     "chebi:101": _d("beta-Citronellol", "QMVPMAAFGQKVCJ-JTQLQIEISA-N"),
     "chebi:102": _d("(-)-Citronellol", "QMVPMAAFGQKVCJ-SNVBAGLBSA-N")}
TERMS = {"kegg.compound:C01401": {"CHEBI:16449"}, "kegg.compound:C00041": {"CHEBI:16977"},
         "kegg.compound:C00133": {"CHEBI:15570"}, "kegg.compound:C00183": {"CHEBI:16414"},
         "kegg.compound:C06417": {"CHEBI:27266"}, "kegg.compound:C00031": {"CHEBI:4167"},
         "kegg.compound:C00221": {"CHEBI:15903"}, "chebi:17522": {"CHEBI:17522"},
         "kegg.compound:C00116": {"CHEBI:17754"}, "chebi:101": {"CHEBI:101"}, "chebi:102": {"CHEBI:102"}}


def test_links_need_isa_or_enantiomer_and_same_first_block():
    s, e, _ = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00041", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00133", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00041", "kegg.compound:C00133") in e
    assert not any("chebi:17522" in p for p in s)          # class scrubbed: no InChIKey


def test_enantiomers_without_unspecified_node_stored_once():
    _, e, _ = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00183", "kegg.compound:C06417") in e
    assert ("kegg.compound:C06417", "kegg.compound:C00183") not in e


def test_families_and_unspecified():
    s, e, _ = stereo_links(T, TERMS, CHEBI)
    fam = stereo_families(set(T), s, e)
    assert fam["kegg.compound:C00041"] == fam["kegg.compound:C00133"] == fam["kegg.compound:C01401"]
    assert "kegg.compound:C01401" in unspecified_members(s)


def test_anomer_guard():
    assert is_anomer_pair("kegg.compound:C00031", "kegg.compound:C00221", T, TERMS, CHEBI)
    assert not is_anomer_pair("chebi:101", "chebi:102", T, TERMS, CHEBI)   # beta-Citronellol: not a sugar
    assert not is_anomer_pair("kegg.compound:C00041", "kegg.compound:C00133", T, TERMS, CHEBI)


def test_finalize_drops_pairs_with_non_survivors():
    s, e, _ = stereo_links(T, TERMS, CHEBI)
    table = {k: dict(v) for k, v in T.items()}
    survivors = set(table) - {"kegg.compound:C00133"}
    counts = finalize_stereo(table, s, e, survivors)
    assert "enantiomer_of" not in table["kegg.compound:C00041"]
    assert table["kegg.compound:C00041"]["stereo_form_of"] == ["kegg.compound:C01401"]
    assert table["kegg.compound:C00183"]["enantiomer_of"] == ["kegg.compound:C06417"]
    assert counts["enantiomer_of"] >= 1


def test_geometric_isomers_are_never_stereo_linked():
    aco = dict(gformula="C6H3O6", charge=-3)
    chebi = make_chebi(CHEBI_22211=term("aconitate", **aco),
                       CHEBI_16383=term("cis-aconitate", is_a=["CHEBI:22211"], **aco),
                       CHEBI_16449=term("alanine", **ALA),
                       CHEBI_16977=term("L-alanine", is_a=["CHEBI:16449"], **ALA))
    t = {"kegg.compound:C02341": _d("Aconitate", "GTZCVFVGUGFEME-UHFFFAOYSA-N"),
         "kegg.compound:C00417": _d("cis-Aconitate", "GTZCVFVGUGFEME-IWQZZHSRSA-N"),
         "kegg.compound:C01401": _d("Alanine", "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"),
         "kegg.compound:C00041": _d("L-Alanine", "QNAYBMKLOCPYGJ-REOHCLBHSA-N")}
    terms = {"kegg.compound:C02341": {"CHEBI:22211"}, "kegg.compound:C00417": {"CHEBI:16383"},
             "kegg.compound:C01401": {"CHEBI:16449"}, "kegg.compound:C00041": {"CHEBI:16977"}}
    s, e, _ = stereo_links(t, terms, chebi)
    assert not any("kegg.compound:C00417" in p for p in s | e)
    assert ("kegg.compound:C00041", "kegg.compound:C01401") in s


def test_strip_geometric_tokens():
    from multiomics_kg.utils.metabolite_identity.names import strip_geometric_tokens
    for name, want in [("cis-Aconitate", "aconitate"), ("trans-Aconitate", "aconitate"),
                       ("(2E)-hexenal", "hexenal"), ("(2Z,4E)-dienoate", "dienoate"),
                       ("all-trans-Retinal", "retinal"), ("L-Alanine", "l-alanine"),
                       ("(3S)-citramalate", "(3s)-citramalate")]:
        assert strip_geometric_tokens(name) == want


# ── explorer ask 8 (2026-10-16/18): ChEBI path + formula / charge rule ─────────────────────────

SUG = make_chebi(
    CHEBI_37684=term("mannose", gformula="C6H12O6", charge=0),                 # no structure
    CHEBI_4208=term("D-mannose", is_a=["CHEBI:37684"], gformula="C6H12O6", charge=0),
    CHEBI_9=term("D-mannopyranose", is_a=["CHEBI:4208"], gformula="C6H12O6", charge=0),
    CHEBI_10=term("beta-D-mannopyranose", is_a=["CHEBI:9"], gformula="C6H12O6", charge=0),
    CHEBI_20=term("hexose", gformula=None),                                    # class
    CHEBI_21=term("D-aldohexose x", is_a=["CHEBI:20"], gformula="C6H12O6", charge=0),
    CHEBI_30=term("mannose 6-phosphate", gformula="C6H13O9P", charge=0),
    CHEBI_31=term("D-mannose 6-phosphate", is_a=["CHEBI:30"], gformula="C6H13O9P", charge=0),
    CHEBI_32=term("D-mannose 6-phosphate(2-)", is_a=["CHEBI:31"], gformula="C6H11O9P", charge=-2,
                  conj=["CHEBI:31"]))


def _sug(name, ik=None, kind="defined"):
    e = cpd(name, inchikey=ik)
    e["structure_kind"] = kind
    return e


def test_two_step_path_without_structure_on_the_generic():
    t = {"chebi:37684": _sug("mannose"), "kegg.compound:C00159": _sug("D-Mannose", "WQZG-QTV")}
    s, _, _ = stereo_links(t, {"chebi:37684": {"CHEBI:37684"}, "kegg.compound:C00159": {"CHEBI:9"}}, SUG)
    assert s == {("kegg.compound:C00159", "chebi:37684")}


def test_path_stops_at_formula_or_charge_change_and_at_classes():
    # D-mannose 6-phosphate(2-) climbs to D-mannose 6-phosphate only if charge matched: it doesn't
    t = {"a": _sug("M6P"), "b": _sug("D-M6P(2-)"), "c": _sug("hexose", kind="generic"), "d": _sug("x")}
    s, _, _ = stereo_links(t, {"a": {"CHEBI:31"}, "b": {"CHEBI:32"}, "c": {"CHEBI:20"}, "d": {"CHEBI:21"}}, SUG)
    assert ("b", "a") not in s
    assert not any(p[1] == "c" for p in s)


def test_conjugate_members_count_on_either_side():
    # the generic node publishes both the acid and the anion; the path is written against the acid
    t = {"g": _sug("mannose 6-phosphate"), "s": _sug("D-mannose 6-phosphate")}
    s, _, _ = stereo_links(t, {"g": {"CHEBI:30"}, "s": {"CHEBI:31", "CHEBI:32"}}, SUG)
    assert s == {("s", "g")}


def test_transitive_reduction_curated_blocked_and_mutual():
    t = {"m": _sug("mannose"), "dm": _sug("D-Mannose"), "b": _sug("beta-D-Mannose"),
         "x": _sug("x"), "y": _sug("y"), "h": _sug("D-Hexose")}
    terms = {"m": {"CHEBI:37684"}, "dm": {"CHEBI:9"}, "b": {"CHEBI:10"}, "x": set(), "y": set(),
             "h": {"CHEBI:21"}}
    s, _, rep = stereo_links(t, terms, SUG, curated=[("x", "m"), ("y", "x"), ("x", "y")],
                             blocked=frozenset({"h"}))
    assert ("b", "dm") in s and ("dm", "m") in s and ("b", "m") not in s     # reduced
    assert ("x", "m") in s                                                  # curated
    assert rep["mutual_dropped"] == [["x", "y"]]
    assert not any("h" in p for p in s)                                     # blocked


def test_enantiomers_are_not_stereo_forms():
    s, e, _ = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00041", "kegg.compound:C00133") in e
    assert not ({("kegg.compound:C00041", "kegg.compound:C00133"),
                 ("kegg.compound:C00133", "kegg.compound:C00041")} & s)


def test_formula_less_intermediate_is_passed_through_never_a_result():
    chebi = make_chebi(CHEBI_26546=term("rhamnose", gformula="C6H12O5", charge=0),
                       CHEBI_62345=term("L-rhamnose", is_a=["CHEBI:26546"]),                 # no formula
                       CHEBI_62346=term("L-rhamnopyranose", is_a=["CHEBI:62345"], gformula="C6H12O5", charge=0))
    t = {"g": _sug("rhamnose"), "lr": _sug("L-Rhamnose-only"), "s": _sug("L-Rhamnose")}
    s, _, _ = stereo_links(t, {"g": {"CHEBI:26546"}, "lr": {"CHEBI:62345"}, "s": {"CHEBI:62346"}}, chebi)
    assert ("s", "g") in s and not any(p[1] == "lr" for p in s)


def test_curated_direction_wins_over_computed_reverse():
    t = {"m": _sug("mannose"), "dm": _sug("D-Mannose")}
    s, _, rep = stereo_links(t, {"m": {"CHEBI:37684"}, "dm": {"CHEBI:9"}}, SUG, curated=[("m", "dm")])
    assert s == {("m", "dm")} and rep["mutual_dropped"] == []
