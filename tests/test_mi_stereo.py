from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.stereo import (
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)

CHEBI = make_chebi(
    CHEBI_16449=term("alanine"),
    CHEBI_16977=term("L-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:15570"]),
    CHEBI_15570=term("D-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:16977"]),
    CHEBI_16414=term("L-valine", enan=["CHEBI:27266"]),
    CHEBI_27266=term("D-valine", enan=["CHEBI:16414"]),
    CHEBI_16646=term("carbohydrate"),
    CHEBI_4167=term("D-glucopyranose", is_a=["CHEBI:16646"]),
    CHEBI_15903=term("beta-D-glucose", is_a=["CHEBI:4167"]),
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
    s, e = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00041", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00133", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00041", "kegg.compound:C00133") in e
    assert not any("chebi:17522" in p for p in s)          # class scrubbed: no InChIKey


def test_enantiomers_without_unspecified_node_stored_once():
    _, e = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00183", "kegg.compound:C06417") in e
    assert ("kegg.compound:C06417", "kegg.compound:C00183") not in e


def test_families_and_unspecified():
    s, e = stereo_links(T, TERMS, CHEBI)
    fam = stereo_families(set(T), s, e)
    assert fam["kegg.compound:C00041"] == fam["kegg.compound:C00133"] == fam["kegg.compound:C01401"]
    assert "kegg.compound:C01401" in unspecified_members(s)


def test_anomer_guard():
    assert is_anomer_pair("kegg.compound:C00031", "kegg.compound:C00221", T, TERMS, CHEBI)
    assert not is_anomer_pair("chebi:101", "chebi:102", T, TERMS, CHEBI)   # beta-Citronellol: not a sugar
    assert not is_anomer_pair("kegg.compound:C00041", "kegg.compound:C00133", T, TERMS, CHEBI)


def test_finalize_drops_pairs_with_non_survivors():
    s, e = stereo_links(T, TERMS, CHEBI)
    table = {k: dict(v) for k, v in T.items()}
    survivors = set(table) - {"kegg.compound:C00133"}
    counts = finalize_stereo(table, s, e, survivors)
    assert "enantiomer_of" not in table["kegg.compound:C00041"]
    assert table["kegg.compound:C00041"]["stereo_form_of"] == ["kegg.compound:C01401"]
    assert table["kegg.compound:C00183"]["enantiomer_of"] == ["kegg.compound:C06417"]
    assert counts["enantiomer_of"] >= 1


def test_geometric_isomers_are_never_stereo_linked():
    chebi = make_chebi(CHEBI_22211=term("aconitate"),
                       CHEBI_16383=term("cis-aconitate", is_a=["CHEBI:22211"]),
                       CHEBI_16449=term("alanine"),
                       CHEBI_16977=term("L-alanine", is_a=["CHEBI:16449"]))
    t = {"kegg.compound:C02341": _d("Aconitate", "GTZCVFVGUGFEME-UHFFFAOYSA-N"),
         "kegg.compound:C00417": _d("cis-Aconitate", "GTZCVFVGUGFEME-IWQZZHSRSA-N"),
         "kegg.compound:C01401": _d("Alanine", "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"),
         "kegg.compound:C00041": _d("L-Alanine", "QNAYBMKLOCPYGJ-REOHCLBHSA-N")}
    terms = {"kegg.compound:C02341": {"CHEBI:22211"}, "kegg.compound:C00417": {"CHEBI:16383"},
             "kegg.compound:C01401": {"CHEBI:16449"}, "kegg.compound:C00041": {"CHEBI:16977"}}
    s, e = stereo_links(t, terms, chebi)
    assert not any("kegg.compound:C00417" in p for p in s | e)
    assert ("kegg.compound:C00041", "kegg.compound:C01401") in s


def test_strip_geometric_tokens():
    from multiomics_kg.utils.metabolite_identity.names import strip_geometric_tokens
    for name, want in [("cis-Aconitate", "aconitate"), ("trans-Aconitate", "aconitate"),
                       ("(2E)-hexenal", "hexenal"), ("(2Z,4E)-dienoate", "dienoate"),
                       ("all-trans-Retinal", "retinal"), ("L-Alanine", "l-alanine"),
                       ("(3S)-citramalate", "(3s)-citramalate")]:
        assert strip_geometric_tokens(name) == want
