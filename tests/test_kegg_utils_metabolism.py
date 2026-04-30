"""Unit tests for the 5 new KEGG endpoints introduced in Spec 1.2."""
from __future__ import annotations

import textwrap

from multiomics_kg.utils import kegg_utils


REACTION_LIST_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpyruvate kinase reaction
    rn:R00010\tphosphofructokinase reaction
    invalid_line
    rn:R12345\t
""")

COMPOUND_LIST_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tD-glucose; alpha-D-glucopyranose
    cpd:C00002\tATP; adenosine 5'-triphosphate
    bad_prefix:Cxxxx\tjunk
""")

LINK_CR_FIXTURE = textwrap.dedent("""\
    cpd:C00074\trn:R00200
    cpd:C00008\trn:R00200
    cpd:C00031\trn:R00010
""")

LINK_PR_FIXTURE = textwrap.dedent("""\
    rn:R00200\tpath:rn00010
    rn:R00200\tpath:rn00710
    rn:R00010\tpath:rn00010
""")

LINK_PC_FIXTURE = textwrap.dedent("""\
    cpd:C00031\tpath:map00010
    cpd:C00031\tpath:map00500
""")


def test_parse_reaction_names():
    out = kegg_utils._parse_reaction_names(REACTION_LIST_FIXTURE)
    assert out == {
        "R00200": "pyruvate kinase reaction",
        "R00010": "phosphofructokinase reaction",
        "R12345": "",
    }


def test_parse_compound_names():
    out = kegg_utils._parse_compound_names(COMPOUND_LIST_FIXTURE)
    # Multi-name field: take first synonym (before first '; ')
    assert out["C00031"] == "D-glucose"
    assert out["C00002"] == "ATP"
    assert "Cxxxx" not in out


def test_parse_reaction_to_compounds():
    out = kegg_utils._parse_reaction_to_compounds(LINK_CR_FIXTURE)
    assert sorted(out["R00200"]) == ["C00008", "C00074"]
    assert out["R00010"] == ["C00031"]


def test_parse_reaction_to_pathways_strips_rn_prefix():
    """KEGG `/link/pathway/reaction` returns rn-prefixed pathway IDs.
    We normalize them to the ko-prefixed form used elsewhere in the KG.
    """
    out = kegg_utils._parse_reaction_to_pathways(LINK_PR_FIXTURE)
    assert sorted(out["R00200"]) == ["ko00010", "ko00710"]
    assert out["R00010"] == ["ko00010"]


def test_parse_compound_to_pathways_strips_map_prefix():
    out = kegg_utils._parse_compound_to_pathways(LINK_PC_FIXTURE)
    # map-prefixed pathways → ko-prefixed for consistency with KeggTerm node IDs
    assert sorted(out["C00031"]) == ["ko00010", "ko00500"]
