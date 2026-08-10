"""Unit tests for the pure InterPro reference parser (prepare_data step 9 core)."""

from __future__ import annotations

from multiomics_kg.utils.interpro_reference import (
    build_reference,
    normalize_type,
    parse_entry_db_xrefs,
    parse_entry_list,
    parse_interpro2go,
    parse_parent_child_tree,
    parse_pathway_xrefs,
)

ENTRY_LIST = "\n".join(
    [
        "ENTRY_AC\tENTRY_TYPE\tENTRY_NAME",
        "IPR000001\tDomain\tKringle",
        "IPR000002\tFamily\tCdc20/Fizzy",
        "IPR000003\tHomologous_superfamily\tRetinoid X receptor",
        "IPR000004\tConserved_site\tSome conserved site",
        "IPR000685\tFamily\tRuBisCO large subunit",
        "malformed line without tabs",
    ]
)

# Depth = leading "--" pairs. IPR000002 is a child of IPR000001; IPR000004 a
# grandchild. IPR000685 is a separate root. IPR000003 is not in the tree.
PARENT_CHILD = "\n".join(
    [
        "IPR000001::Kringle::",
        "--IPR000002::Cdc20/Fizzy::",
        "----IPR000004::Some conserved site::",
        "IPR000685::RuBisCO large subunit::",
    ]
)


# Real interpro2go shape: a `!` comment header, then one line per (entry, GO).
# IPR000685 carries two terms; the malformed tail line must be ignored.
INTERPRO2GO = "\n".join(
    [
        "!date: 2026/05/24 10:10:21",
        "!Mapping of InterPro entries to GO",
        "InterPro:IPR000001 Kringle > GO:protein binding ; GO:0005515",
        "InterPro:IPR000685 RuBisCO large subunit > GO:magnesium ion binding ; GO:0000287",
        "InterPro:IPR000685 RuBisCO large subunit > GO:RuBisCO activity ; GO:0016984",
        "InterPro:IPR000685 RuBisCO large subunit > GO:duplicate ; GO:0016984",
        "InterPro:IPR999999 Retired entry > GO:nothing ; GO:0000001",
        "InterPro:IPR000002 no semicolon or GO id here",
    ]
)

# Real interpro.xml shape. IPR000001 carries a Reactome xref (species-expanded)
# plus MetaCyc and a CAZy family; IPR000685 MetaCyc + an EC number. The PUBMED
# xref, the stray MetaCyc/EC xrefs outside any <interpro> element, and the
# member-list PFAM signature must never be picked up as entry xrefs.
INTERPRO_XML = """<?xml version="1.0"?>
<interprodb>
<db_xref db="METACYC" dbkey="PWY-HEADER-LEAK"/>
<db_xref db="EC" dbkey="9.9.9.9"/>
<interpro id="IPR000001" short_name="Kringle" type="Domain">
  <db_xref db="PUBMED" dbkey="12345678"/>
  <db_xref db="REACTOME" dbkey="R-BTA-6798695"/>
  <db_xref db="METACYC" dbkey="PWY-1042"/>
  <db_xref db="CAZY" dbkey="CBM5"/>
  <member_list>
    <db_xref db="PFAM" dbkey="PF00051"/>
  </member_list>
</interpro>
<interpro id="IPR000685" short_name="RuBisCO" type="Family">
  <db_xref db="METACYC" dbkey="PWY-5723"/>
  <db_xref db="METACYC" dbkey="PWY-5532"/>
  <db_xref db="EC" dbkey="4.1.1.39"/>
</interpro>
</interprodb>
"""


def test_normalize_type_matches_callsjson_upper():
    assert normalize_type("Homologous_superfamily") == "HOMOLOGOUS_SUPERFAMILY"
    assert normalize_type("Conserved_site") == "CONSERVED_SITE"
    assert normalize_type("Family") == "FAMILY"
    assert normalize_type(None) is None


def test_parse_entry_list_skips_header_and_malformed():
    entries = parse_entry_list(ENTRY_LIST)
    assert set(entries) == {"IPR000001", "IPR000002", "IPR000003", "IPR000004", "IPR000685"}
    assert entries["IPR000001"] == {"type": "DOMAIN", "name": "Kringle"}
    assert entries["IPR000003"]["type"] == "HOMOLOGOUS_SUPERFAMILY"


def test_parse_parent_child_tree_depth_and_parent():
    tree = parse_parent_child_tree(PARENT_CHILD)
    assert tree["IPR000001"] == {"parent": None, "level": 0}
    assert tree["IPR000002"] == {"parent": "IPR000001", "level": 1}
    assert tree["IPR000004"] == {"parent": "IPR000002", "level": 2}
    assert tree["IPR000685"] == {"parent": None, "level": 0}


def test_build_reference_combines_and_defaults():
    ref = build_reference(ENTRY_LIST, PARENT_CHILD)
    # name/type from entry.list; parent/level from tree
    assert ref["IPR000002"] == {
        "name": "Cdc20/Fizzy",
        "type": "FAMILY",
        "parent": "IPR000001",
        "level": 1,
    }
    # entry not in the tree → parentless root at level 0
    assert ref["IPR000003"]["parent"] is None
    assert ref["IPR000003"]["level"] == 0
    assert ref["IPR000003"]["type"] == "HOMOLOGOUS_SUPERFAMILY"


def test_build_reference_hierarchy_never_dangles():
    """Every entry's parent (when set) is itself a key in the reference."""
    ref = build_reference(ENTRY_LIST, PARENT_CHILD)
    for acc, meta in ref.items():
        if meta["parent"] is not None:
            assert meta["parent"] in ref, f"{acc} parent {meta['parent']} missing"


# --------------------------------------------------------------------------
# interpro2go → go_terms
# --------------------------------------------------------------------------

def test_parse_interpro2go_dedupes_and_skips_comments():
    go = parse_interpro2go(INTERPRO2GO)
    assert go["IPR000001"] == ["GO:0005515"]
    # two distinct terms; the duplicated GO:0016984 line collapses
    assert go["IPR000685"] == ["GO:0000287", "GO:0016984"]
    # `!` comment lines and the malformed line contribute nothing
    assert "IPR000002" not in go


# --------------------------------------------------------------------------
# interpro.xml → pathways
# --------------------------------------------------------------------------

def test_parse_pathway_xrefs_metacyc_only_by_default():
    pw = parse_pathway_xrefs(INTERPRO_XML.splitlines())
    assert pw["IPR000001"] == ["MetaCyc:PWY-1042"]  # Reactome excluded by default
    assert pw["IPR000685"] == ["MetaCyc:PWY-5532", "MetaCyc:PWY-5723"]


def test_parse_pathway_xrefs_uses_interproscan_db_casing():
    """`MetaCyc:PWY-1042`, matching the DB:id form in calls.json — not `METACYC:`."""
    pw = parse_pathway_xrefs(INTERPRO_XML.splitlines(), include_dbs=("METACYC", "REACTOME"))
    assert pw["IPR000001"] == ["MetaCyc:PWY-1042", "Reactome:R-BTA-6798695"]


def test_parse_pathway_xrefs_ignores_non_pathway_and_out_of_entry_xrefs():
    pw = parse_pathway_xrefs(INTERPRO_XML.splitlines())
    # the PUBMED xref inside IPR000001 is not a pathway
    assert not any("12345678" in p for pws in pw.values() for p in pws)
    # the xref before the first <interpro> element is attributed to nothing
    assert not any("HEADER-LEAK" in p for pws in pw.values() for p in pws)


# --------------------------------------------------------------------------
# build_reference with xrefs
# --------------------------------------------------------------------------

def test_build_reference_attaches_go_and_pathways():
    ref = build_reference(
        ENTRY_LIST, PARENT_CHILD,
        go_map=parse_interpro2go(INTERPRO2GO),
        pathway_map=parse_pathway_xrefs(INTERPRO_XML.splitlines()),
    )
    assert ref["IPR000685"]["go_terms"] == ["GO:0000287", "GO:0016984"]
    assert ref["IPR000685"]["pathways"] == ["MetaCyc:PWY-5532", "MetaCyc:PWY-5723"]
    assert ref["IPR000001"]["go_terms"] == ["GO:0005515"]


def test_build_reference_xref_fields_are_sparse():
    """Entries without GO/pathways omit the keys entirely (keeps the JSON small)."""
    ref = build_reference(
        ENTRY_LIST, PARENT_CHILD,
        go_map=parse_interpro2go(INTERPRO2GO),
        pathway_map=parse_pathway_xrefs(INTERPRO_XML.splitlines()),
    )
    assert "go_terms" not in ref["IPR000003"]
    assert "pathways" not in ref["IPR000003"]
    # IPR000001 has GO but no... it does have MetaCyc; IPR000002 has neither
    assert "go_terms" not in ref["IPR000002"]
    assert "pathways" not in ref["IPR000002"]


def test_build_reference_drops_xrefs_for_unknown_accessions():
    """A GO mapping for an accession absent from entry.list must not create a node."""
    ref = build_reference(ENTRY_LIST, PARENT_CHILD, go_map=parse_interpro2go(INTERPRO2GO))
    assert "IPR999999" not in ref


def test_build_reference_without_xref_maps_is_unchanged():
    """Backward compatibility: the two-arg call still produces the original shape."""
    ref = build_reference(ENTRY_LIST, PARENT_CHILD)
    assert ref["IPR000685"] == {
        "name": "RuBisCO large subunit",
        "type": "FAMILY",
        "parent": None,
        "level": 0,
    }


# --------------------------------------------------------------------------
# interpro.xml → EC / CAZy (generic db_xref extractor)
# --------------------------------------------------------------------------

def _ec_cazy_maps():
    raw = parse_entry_db_xrefs(INTERPRO_XML.splitlines(), include_dbs=("EC", "CAZY"))
    ec_map = {acc: dbs["EC"] for acc, dbs in raw.items() if "EC" in dbs}
    cazy_map = {acc: dbs["CAZY"] for acc, dbs in raw.items() if "CAZY" in dbs}
    return ec_map, cazy_map


def test_parse_entry_db_xrefs_returns_raw_keys_grouped_by_db():
    raw = parse_entry_db_xrefs(INTERPRO_XML.splitlines(), include_dbs=("EC", "CAZY"))
    assert raw["IPR000685"]["EC"] == ["4.1.1.39"]  # raw, no label prefix
    assert raw["IPR000001"]["CAZY"] == ["CBM5"]
    # IPR000001 has no EC; IPR000685 has no CAZy
    assert "EC" not in raw.get("IPR000001", {})
    assert "CAZY" not in raw.get("IPR000685", {})


def test_parse_entry_db_xrefs_ignores_out_of_entry_and_member_list():
    raw = parse_entry_db_xrefs(
        INTERPRO_XML.splitlines(), include_dbs=("EC", "CAZY", "METACYC", "PFAM")
    )
    flat = {key for dbs in raw.values() for keys in dbs.values() for key in keys}
    assert "9.9.9.9" not in flat          # EC before the first <interpro> element
    assert "PWY-HEADER-LEAK" not in flat  # MetaCyc header leak
    # PFAM inside <member_list> IS attributed to the open entry (it is inside it);
    # it is simply never in EC/CAZy's include set at build time. Confirm the guard
    # is about entry boundaries, not element nesting:
    assert raw["IPR000001"].get("PFAM") == ["PF00051"]


def test_build_reference_attaches_ec_and_cazy_sparse():
    ec_map, cazy_map = _ec_cazy_maps()
    ref = build_reference(ENTRY_LIST, PARENT_CHILD, ec_map=ec_map, cazy_map=cazy_map)
    assert ref["IPR000685"]["ec_numbers"] == ["4.1.1.39"]
    assert ref["IPR000001"]["cazy_ids"] == ["CBM5"]
    # sparse: absent where the entry has none
    assert "ec_numbers" not in ref["IPR000001"]
    assert "cazy_ids" not in ref["IPR000685"]
    assert "ec_numbers" not in ref["IPR000002"]


def test_build_reference_ec_numbers_stored_raw_unnormalized():
    """EC is stored verbatim; normalize_ec / bare-3-level is the consumer's job."""
    ref = build_reference(ENTRY_LIST, PARENT_CHILD, ec_map={"IPR000685": ["3.4.21"]})
    assert ref["IPR000685"]["ec_numbers"] == ["3.4.21"]
