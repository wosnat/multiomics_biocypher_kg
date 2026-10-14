import pytest

from multiomics_kg.utils.metabolite_identity import names as N


@pytest.mark.parametrize("a,b", [("Glutamic acid", "glutamic acid"), ("NAD", "nad"),
                                 ("  L-Alanine ", "l-alanine"), ("GABA", "gaba")])
def test_name_key_folds(a, b):
    assert N.name_key(a) == N.name_key(b)


@pytest.mark.parametrize("a,b", [("Co", "CO"), ("NO", "No"), ("CO2", "Co2"), ("HF", "Hf")])
def test_name_key_keeps_formula_case(a, b):
    assert N.name_key(a) != N.name_key(b)


def test_formula_like():
    assert N.is_formula_like("NH4+") and N.is_formula_like("Fe3+") and N.is_formula_like("CO2")
    assert N.is_formula_like("AMP")          # Am + P: harmless false positive, not folded
    assert not N.is_formula_like("NAD") and not N.is_formula_like("DMSP")


@pytest.mark.parametrize("name", [
    "N-Acetyl-D-glucosamine", "S-Adenosyl-L-methionine", "UDP-D-glucose", "Glycyl-L-leucine",
    "(2R,3S)-isocitrate", "(S)-malate", "L-Alanine", "beta-D-glucose", "(R,S)-lactate",
    "(+)-catechin", "(−)-citronellol", "(-)-citronellol", "meso-2,6-diaminopimelate",
    "DL-Methionine", "alpha-L-fucose"])
def test_stereo_qualified(name):
    assert N.is_stereo_qualified(name)


@pytest.mark.parametrize("name", [
    "beta-Citronellol", "alpha-Ketoglutarate", "S-Adenosylmethionine", "Alanine",
    "Glutamate", "R-lipoate", "trans-Aconitate", "cis-Aconitate", "(E)-cinnamate"])
def test_stereo_unqualified(name):
    assert not N.is_stereo_qualified(name)


def test_strip_anomer_tokens():
    assert N.strip_anomer_tokens("beta-D-Glucose") == N.strip_anomer_tokens("D-Glucose")
    assert (N.strip_anomer_tokens("N-Acetyl-beta-D-glucosamine")
            == N.strip_anomer_tokens("N-Acetyl-D-glucosamine"))


@pytest.mark.parametrize("f,generic,polymer", [
    ("C6H12O6", False, False), ("Rb", False, False), ("Ru", False, False), ("ReO4", False, False),
    ("C5H8O2R", True, False), ("C5H8O2Rn", True, False), ("C3H3O6R", True, False),
    ("H2A", True, False), ("RX", True, False), ("C10H17N3O6S*", True, False),
    ("(C6H10O5)n", False, True), ("C12H22O11(C6H10O5)n", False, True), ("(C4H6O2)n", False, True),
    ("C10H17O8PR2(C5H8O5PR)n", True, True)])
def test_formula_classes(f, generic, polymer):
    assert N.is_generic_formula(f) is generic
    assert N.is_polymer_formula(f) is polymer


@pytest.mark.parametrize("f,els", [
    ("C6H12O6", ["C", "H", "O"]), ("Rb", ["Rb"]), ("C5H8O2Rn", ["C", "H", "O"]),
    ("C3H3O6R", ["C", "H", "O"]), ("(C6H10O5)n", ["C", "H", "O"]), ("H2A", ["H"]),
    ("Fe+3", ["Fe"]), (None, []), ("", [])])
def test_parse_elements_allow_list(f, els):
    assert N.parse_elements(f) == els


def test_normalize_hmdb():
    assert N.normalize_hmdb("HMDB02878") == "HMDB0002878"
    assert N.normalize_hmdb("HMDB0002878") == "HMDB0002878"
    assert N.normalize_hmdb("") is None and N.normalize_hmdb(None) is None
