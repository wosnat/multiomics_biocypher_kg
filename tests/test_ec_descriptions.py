"""Expasy enzyme.dat DE lines: pypath keeps only the last one; we join them."""
from multiomics_kg.adapters.ec_adapter import enzyme_dat_descriptions, repair_truncated_descriptions

DAT = """\
ID   1.7.99.4
DE   Transferred entry: 1.7.1.1, 1.7.1.2, 1.7.1.3, 1.7.5.1, 1.7.7.2 and
DE   1.9.6.1.
//
ID   1.13.11.90
DE   [1-hydroxy-2-(trimethylamino)ethyl]phosphonate dioxygenase (glycine-
DE   betaine-forming).
//
ID   1.1.1.1
DE   alcohol dehydrogenase.
//
"""


def test_descriptions_join_continuation_lines():
    d = enzyme_dat_descriptions(DAT)
    assert d["1.7.99.4"] == "Transferred entry: 1.7.1.1, 1.7.1.2, 1.7.1.3, 1.7.5.1, 1.7.7.2 and 1.9.6.1"
    # Expasy wraps at a hyphen: no space is inserted.
    assert d["1.13.11.90"] == "[1-hydroxy-2-(trimethylamino)ethyl]phosphonate dioxygenase (glycine-betaine-forming)"
    assert d["1.1.1.1"] == "alcohol dehydrogenase"


def test_repair_only_touches_truncated_entries_and_prunes_moved_ones():
    enzymes = {
        "1.7.99.4": {"de": "1.9.6.1"},                     # truncated, really transferred
        "1.13.11.90": {"de": "betaine-forming)"},          # truncated name
        "1.1.1.1": {"de": "alcohol dehydrogenase (old)"},  # differs upstream, not a truncation
    }
    ec_dict = {"1.-.-.-": {"1.7.-.-": {"1.7.99.-": {"entries": ["1.7.99.4", "1.7.99.5"]}}}}
    repaired = repair_truncated_descriptions(enzymes, ec_dict, enzyme_dat_descriptions(DAT))
    assert repaired == ["1.13.11.90", "1.7.99.4"]
    assert enzymes["1.7.99.4"]["de"].startswith("Transferred entry:")
    assert enzymes["1.1.1.1"]["de"] == "alcohol dehydrogenase (old)"
    assert ec_dict["1.-.-.-"]["1.7.-.-"]["1.7.99.-"]["entries"] == ["1.7.99.5"]
