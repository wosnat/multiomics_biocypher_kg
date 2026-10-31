"""Tests for utils/tcdb_utils.py — TCDB hierarchy accessor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.utils import tcdb_utils as tu


TINY_TCDB = {
    "1":         {"name": "Channels and Pores", "level": 0, "level_kind": "tc_class", "parent": None},
    "1.A":       {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
    "1.A.1":     {"name": "VIC Family", "level": 2, "level_kind": "tc_family", "parent": "1.A", "abbreviation": "VIC"},
    "1.A.1.1":   {"name": "", "level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.1.1": {"name": "", "level": 4, "level_kind": "tc_system", "parent": "1.A.1.1",
                  "substrate_classes": ["potassium(1+)"]},
}


@pytest.fixture(autouse=True)
def patch_tcdb(monkeypatch, tmp_path):
    """Point load_tcdb at a tiny synthetic JSON; reset the module cache."""
    p = tmp_path / "tcdb_hierarchy.json"
    p.write_text(json.dumps(TINY_TCDB))
    monkeypatch.setattr(tu, "DEFAULT_PATH", p)
    monkeypatch.setattr(tu, "_CACHE", None)
    yield


def test_load_tcdb_returns_dict():
    h = tu.load_tcdb()
    assert isinstance(h, dict)
    assert h["1"]["level"] == 0


def test_load_tcdb_caches():
    """Second call returns the same object without re-reading the file."""
    a = tu.load_tcdb()
    b = tu.load_tcdb()
    assert a is b


def test_is_valid_tcdb_present():
    assert tu.is_valid_tcdb("1.A.1.1.1") is True
    assert tu.is_valid_tcdb("1.A") is True
    assert tu.is_valid_tcdb("1") is True


def test_is_valid_tcdb_absent():
    assert tu.is_valid_tcdb("99.X.99") is False
    assert tu.is_valid_tcdb("") is False


def test_tcdb_ancestors_full_chain():
    """Specificity-level TCID → list of all ancestors (root → parent)."""
    assert tu.tcdb_ancestors("1.A.1.1.1") == ["1", "1.A", "1.A.1", "1.A.1.1"]


def test_tcdb_ancestors_partial():
    assert tu.tcdb_ancestors("1.A.1") == ["1", "1.A"]
    assert tu.tcdb_ancestors("1") == []


def test_tcdb_ancestors_unknown_returns_empty():
    assert tu.tcdb_ancestors("99.X.99") == []


# ============================================================================
# Cross-ontology bridges (Pfam / GO -> TcdbFamily)
# ============================================================================

from multiomics_kg.utils.tcdb_utils import build_tcdb_bridges

# 1.A.1.5.2 is kept (gene-annotated); 1.A.1.9.1 is not, and rolls up to 1.A.1.
_H = {
    "1": {"level": 0, "level_kind": "tc_class", "parent": None},
    "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
    "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
    "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.5.2": {"level": 4, "level_kind": "tc_system", "parent": "1.A.1.5"},
    "1.A.1.9": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
    "1.A.1.9.1": {"level": 4, "level_kind": "tc_system", "parent": "1.A.1.9"},
    "2": {"level": 0, "level_kind": "tc_class", "parent": None},
    "2.A": {"level": 1, "level_kind": "tc_subclass", "parent": "2"},
    "2.A.9.1.1": {"level": 4, "level_kind": "tc_system", "parent": "2.A"},
}
_KEPT = {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2", "2", "2.A"}


def test_bridge_attaches_directly_when_tcid_survives():
    b = build_tcdb_bridges({"PF00001": {"1.A.1.5.2"}}, _H, _KEPT)
    assert b == {"1.A.1.5.2": {"PF00001": ["1.A.1.5.2"]}}


def test_bridge_rolls_up_to_nearest_surviving_ancestor():
    """1.A.1.9.1 has no node; its nearest kept ancestor is the family 1.A.1."""
    b = build_tcdb_bridges({"PF00002": {"1.A.1.9.1"}}, _H, _KEPT)
    assert b == {"1.A.1": {"PF00002": ["1.A.1.9.1"]}}


def test_bridge_drops_targets_shallower_than_family():
    """2.A.9.1.1 rolls up only as far as the subclass 2.A -> uninformative, dropped.

    An edge saying "this domain relates to Electrochemical Potential-driven
    Transporters" is true but carries no usable signal.
    """
    assert build_tcdb_bridges({"PF00003": {"2.A.9.1.1"}}, _H, _KEPT) == {}


def test_bridge_records_all_curated_tcids_after_rollup():
    """Roll-up must not lose provenance: every original TCID is retained."""
    b = build_tcdb_bridges({"PF00004": {"1.A.1.9.1", "1.A.1.5.2"}}, _H, _KEPT)
    assert b["1.A.1"]["PF00004"] == ["1.A.1.9.1"]
    assert b["1.A.1.5.2"]["PF00004"] == ["1.A.1.5.2"]


def test_bridge_curated_tcids_are_sorted_and_deduped():
    b = build_tcdb_bridges({"PF00005": {"1.A.1.9.1"}}, _H, _KEPT | {"1.A.1.9"})
    assert b["1.A.1.9"]["PF00005"] == ["1.A.1.9.1"]


def test_bridge_skips_tcids_absent_from_the_hierarchy():
    assert build_tcdb_bridges({"PF00006": {"9.Z.99.9.9"}}, _H, _KEPT) == {}


def test_bridge_handles_go_ids_identically():
    """The builder is xref-agnostic — GO ids roll up the same way as Pfam."""
    b = build_tcdb_bridges({"GO:0005215": {"1.A.1.9.1"}}, _H, _KEPT)
    assert b == {"1.A.1": {"GO:0005215": ["1.A.1.9.1"]}}


def test_bridge_empty_input_returns_empty():
    assert build_tcdb_bridges({}, _H, _KEPT) == {}


# ── tcdb_direction (spec §5.1; conservative rule after the 2026-10-24 audit) ─

def _dir_h(*nodes):
    """Hierarchy from (id, name[, description]) tuples; levels/parents from the id."""
    kinds = ["tc_class", "tc_subclass", "tc_family", "tc_subfamily", "tc_system"]
    h = {}
    for node in nodes:
        tc, name, desc = (node + ("",))[:3]
        parts = tc.split(".")
        h[tc] = {"name": name, "description": desc, "level": len(parts) - 1,
                 "level_kind": kinds[len(parts) - 1],
                 "parent": ".".join(parts[:-1]) or None}
    return h


_ABC_H = _dir_h(("3", "Primary Active Transporters"), ("3.A", "P-P-bond hydrolysis"),
            ("3.A.1", "The ATP-binding Cassette (ABC) Superfamily"),
            ("3.A.1.4", "The Hydrophobic Amino Acid Uptake Transporter (HAAT) Family"),
            ("3.A.1.4.4", "The high-affinity urea porter"),
            ("3.A.1.4.5", ""),
            ("3.A.1.4.6", "ABC transporter of unknown function"),
            ("3.A.1.4.7", "Putative drug exporter"),
            ("3.A.1.4.8", "General amino acid porter", "General porter (catalyzes both uptake and efflux)."),
            ("3.A.1.4.9", "Peptide uptake and export system"),
            ("3.A.1.22", "The Nickel Uptake Transporter (NiT) Family"),
            ("3.A.1.22.2", "The ABC exporter, CbiMNQO."),
            ("3.A.1.122", "The Macrolide Exporter (MacB) Family"),
            ("3.A.1.122.1", "Macrolide efflux pump, MacAB"),
            ("3.A.1.122.2", ""),
            ("3.A.1.203", "The Peroxysomal Fatty Acyl CoA Transporter Family"),
            ("3.A.1.203.4", "The BacA porter (selective for the uptake of bleomycin)"))


@pytest.mark.parametrize("tc, expected", [
    ("3.A.1.4", ("uptake", "abc_numbering")),       # ABC uptake subfamily, by number
    ("3.A.1.4.4", ("uptake", "abc_numbering")),     # ... and its systems
    ("3.A.1.4.5", ("uptake", "abc_numbering")),     # unnamed system: number, not name
    ("3.A.1.4.6", ("unknown", None)),               # own name: function unknown
    ("3.A.1.4.7", ("unknown", None)),               # tentative export claim still blocks
    ("3.A.1.4.8", ("unknown", None)),               # description: both directions
    ("3.A.1.4.9", ("unknown", None)),               # name states both
    ("3.A.1.22", ("uptake", "abc_numbering")),
    ("3.A.1.22.2", ("unknown", None)),              # number says uptake, name says export
    ("3.A.1", ("unknown", None)),                   # superfamily: numbering starts below
    ("3.A.1.122", ("unknown", None)),               # 3.A.1.100+: family names never decide
    ("3.A.1.122.1", ("export", "tcdb_name")),       # ... only a system's own name does
    ("3.A.1.122.2", ("unknown", None)),             # unnamed: nothing inherited
    ("3.A.1.203.4", ("uptake", "tcdb_name")),       # BacA: importer numbered among exporters
    ("3", ("unknown", None)), ("3.A", ("unknown", None)),
    ("9.Z.9.9.9", ("unknown", None)),               # not in the hierarchy
])
def test_tcdb_direction_abc(tc, expected):
    assert tu.tcdb_direction(tc, _ABC_H, {}) == expected


def _sys(name, desc=""):
    return _dir_h(("2", ""), ("2.A", "Porters"), ("2.A.7", "The Drug Efflux (DE) Family"),
              ("2.A.7.1", "The Sugar Efflux subfamily"), ("2.A.7.1.1", name, desc))


@pytest.mark.parametrize("name, desc, expected", [
    ("Sugar efflux porter, SetA", "", ("export", "tcdb_name")),
    ("Glucose uptake porter, GlcP", "", ("uptake", "tcdb_name")),       # own name beats family
    ("Ammonium uptake permease, AmtB (may need AmtB to concentrate)", "",
     ("uptake", "tcdb_name")),                                           # hedge AFTER the word
    ("Drug Exporter-1 protein", "", ("export", "tcdb_name")),            # hyphenated agent noun
    ("Nucleotide importers", "", ("uptake", "tcdb_name")),               # plural
    ("The importance of being transported", "", ("unknown", None)),     # whole word only
    ("Uptake and export of things", "", ("unknown", None)),             # both
    ("", "", ("unknown", None)),                                        # family name never inherited
    ("Glucose porter of 400 aas", "", ("unknown", None)),
    ("Probable sugar:H+ efflux porter", "", ("unknown", None)),         # hedge across a colon
    ("Putative heavy metal (Me 2+) efflux porter", "", ("unknown", None)),
    ("Lipid efflux homologue of 593 aas", "", ("unknown", None)),       # resemblance after the word
    ("It is not involved in the efflux of drugs", "", ("unknown", None)),  # negation
    ("Heme efflux porter?", "", ("unknown", None)),                     # literal ?
    ("Rv0205 (may function with a heme uptake transporter)", "", ("unknown", None)),
    ("Zn uptake (into synaptic vesicles)", "", ("unknown", None)),      # organelle after the word
    ("Ferrous iron exporter, MbfA, of the erythrin-vacuolar family", "", ("unknown", None)),
    ("Protein import into the ER", "", ("unknown", None)),
    ("Sugar efflux porter", "Catalyzes both uptake and efflux of sugars.", ("unknown", None)),
    ("Sugar efflux porter", "Its role is under debate.", ("unknown", None)),
    ("Sugar efflux porter", "Sugar uptake porter, SetA. Unrelated detail.", ("unknown", None)),
    ("Sugar efflux porter", "Sugar efflux porter, SetA. Assay used dye uptake.",
     ("export", "tcdb_name")),                                          # later sentence ignored
    ("Sugar efflux porter", "SetA. Its partner, of unknown function, binds it.",
     ("export", "tcdb_name")),                                          # unknown function in desc
])
def test_tcdb_direction_system_name(name, desc, expected):
    assert tu.tcdb_direction("2.A.7.1.1", _sys(name, desc), {}) == expected


def test_family_names_never_decide():
    h = _sys("Sugar efflux porter, SetA")
    for tc in ("2", "2.A", "2.A.7", "2.A.7.1"):
        assert tu.tcdb_direction(tc, h, {}) == ("unknown", None), tc


def test_direction_overrides_subtree():
    h = _sys("Sugar efflux porter, SetA")
    ov = {"2.A.7.1": {"scope": "subtree", "reason": "t"}}
    assert tu.tcdb_direction("2.A.7.1.1", h, ov) == ("unknown", None)
    abc = {"3.A.1.4": {"scope": "subtree", "reason": "t"}}
    assert tu.tcdb_direction("3.A.1.4.4", _ABC_H, abc) == ("unknown", None)
    assert tu.tcdb_direction("3.A.1.4.99", _ABC_H, abc) == ("unknown", None)  # not in hierarchy


def test_direction_overrides_file_valid():
    ov = tu.load_direction_overrides()
    assert len(ov) > 50
    hierarchy = json.loads(Path("cache/data/tcdb/tcdb_hierarchy.json").read_text())
    missing = sorted(k for k in ov if k not in hierarchy)
    assert not missing, f"override ids not in the TCDB hierarchy: {missing}"
    for k, v in ov.items():
        assert v["scope"] == "subtree" and v["reason"], k


def test_direction_overrides_rejects_bad_entry(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("2.A.7:\n  scope: node\n  reason: x\n")
    with pytest.raises(ValueError):
        tu.load_direction_overrides(p)


def test_direction_real_hierarchy_examples():
    hierarchy = json.loads(Path("cache/data/tcdb/tcdb_hierarchy.json").read_text())
    d = lambda t: tu.tcdb_direction(t, hierarchy)
    assert d("3.A.1.4") == ("uptake", "abc_numbering")
    assert d("3.A.1.4.4") == ("uptake", "abc_numbering")         # urea ABC (MED4 urtA)
    assert d("3.A.1.27.3") == ("unknown", None)                   # Mla (override)
    assert d("2.A.103.1.1") == ("unknown", None)                  # FtsW (override)
    assert d("3.A.1.3.8") == ("unknown", None)                    # Aap: both directions
    assert d("3.A.1.5.24") == ("unknown", None)                   # YejABEF: exporter in full name
    assert d("2.A.6.1.4") == ("export", "tcdb_name")              # CusABCF
    assert d("2.A.6") == ("unknown", None)                        # RND family: names never decide
    assert d("2.A.86.1.7") == ("unknown", None)                   # old documented false positive


# ── Retired TCDB ids: moved-to map (post-3a P2) ──────────────────────────────

_MOVED_PAGE = (
    "<div>The following TC entries have been moved</br>"
    "9.B.99.1.1 has been moved to <a href='result.php?tc=1.a.17.4.1'>1.a.17.4.1</a></br>"
    "9.B.99.1.1 has been moved to <a href='result.php?tc=1.B.93.1.1'>1.B.93.1.1</a></br>"
    "9.B.99.1.2 has been moved to <a href='result.php?tc=1.B.93.1.2'>1.B.93.1.2</a></br>"
    "</div>"
)


def _h(*ids):
    return {t: {"name": "", "parent": t.rsplit(".", 1)[0] if "." in t else None} for t in ids}


def test_parse_moved_tcids_real_line_format_and_lower_case_target():
    from multiomics_kg.utils.tcdb_utils import parse_moved_tcids
    moves = parse_moved_tcids([_MOVED_PAGE, "<p>nothing here</p>"])
    assert moves == {"9.B.99.1.1": {"1.A.17.4.1", "1.B.93.1.1"},
                     "9.B.99.1.2": {"1.B.93.1.2"}}


def test_resolve_moved_aliases_system_and_coarse_all_agree():
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("1", "1.B", "1.B.166", "1.B.166.1", "1.B.166.1.1", "1.B.166.1.2",
           "1.B.101", "1.B.101.1", "1.B.101.1.1", "1.B.101.1.2")
    raw = {"1.B.166.1.1": {"1.B.101.1.1"}, "1.B.166.1.2": {"1.B.101.1.2"}}
    aliases, ambiguous = resolve_moved_aliases(raw, h)
    assert aliases == {"1.B.166.1.1": "1.B.101.1.1", "1.B.166.1.2": "1.B.101.1.2",
                       "1.B.166.1": "1.B.101.1", "1.B.166": "1.B.101"}
    assert ambiguous == []


def test_resolve_moved_aliases_disagreeing_targets_give_no_coarse_alias():
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("9", "9.A", "9.A.40", "9.A.40.1", "9.A.40.1.1", "9.A.40.1.5", "9.A.40.3",
           "9.A.40.3.5", "1", "1.C", "1.C.126", "1.C.126.1", "1.C.126.1.1")
    raw = {"9.A.40.1.1": {"1.C.126.1.1"}, "9.A.40.1.5": {"9.A.40.3.5"}}
    aliases, ambiguous = resolve_moved_aliases(raw, h)
    assert aliases == {"9.A.40.1.1": "1.C.126.1.1", "9.A.40.1.5": "9.A.40.3.5"}
    assert ambiguous == []


def test_resolve_moved_aliases_coarse_needs_every_system_moved():
    """A subfamily with a system that was NOT moved keeps no coarse alias."""
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("1", "1.B", "1.B.166", "1.B.166.1", "1.B.166.1.1", "1.B.166.1.7",
           "1.B.101", "1.B.101.1", "1.B.101.1.1")
    aliases, _ = resolve_moved_aliases({"1.B.166.1.1": {"1.B.101.1.1"}}, h)
    assert aliases == {"1.B.166.1.1": "1.B.101.1.1"}


def test_resolve_moved_aliases_ambiguous_and_override_wins():
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("9", "9.B", "9.B.99", "9.B.99.1", "9.B.99.1.1", "9.B.99.1.2",
           "1", "1.A", "1.A.17", "1.A.17.4", "1.A.17.4.1",
           "1.B", "1.B.93", "1.B.93.1", "1.B.93.1.1", "1.B.93.1.2")
    raw = {"9.B.99.1.1": {"1.A.17.4.1", "1.B.93.1.1"}, "9.B.99.1.2": {"1.B.93.1.2"}}
    aliases, ambiguous = resolve_moved_aliases(raw, h)
    assert aliases == {"9.B.99.1.2": "1.B.93.1.2"}       # 1.1 ambiguous, no coarse
    assert ambiguous == ["9.B.99.1.1"]
    aliases, _ = resolve_moved_aliases(raw, h, overrides={"9.B.99.1": "1.B.93.1"})
    # the override also decides its ambiguous descendants, by prefix substitution
    assert aliases == {"9.B.99.1.1": "1.B.93.1.1", "9.B.99.1.2": "1.B.93.1.2",
                       "9.B.99.1": "1.B.93.1"}


def test_override_ambiguous_descendant_falls_back_to_override_target():
    """Substituted target missing from the hierarchy -> the override target itself."""
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("9", "9.B", "9.B.99", "9.B.99.1", "9.B.99.1.3", "1", "1.A", "1.A.17",
           "1.A.17.4", "1.A.17.4.2", "1.B", "1.B.93", "1.B.93.1")
    raw = {"9.B.99.1.3": {"1.A.17.4.2", "1.B.93.1.3"}}
    aliases, ambiguous = resolve_moved_aliases(raw, h, overrides={"9.B.99.1": "1.B.93.1"})
    assert aliases == {"9.B.99.1.3": "1.B.93.1", "9.B.99.1": "1.B.93.1"}
    assert ambiguous == ["9.B.99.1.3"]


def test_still_current_old_ids_are_not_aliased():
    """A parsed old id TCDB still lists (raw TSVs) is never aliased."""
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases, still_current_moves
    h = _h("1", "1.B", "1.B.166", "1.B.166.1", "1.B.166.1.1", "1.B.166.1.2",
           "1.B.101", "1.B.101.1", "1.B.101.1.1", "1.B.101.1.2")
    raw = {"1.B.166.1.1": {"1.B.101.1.1"}, "1.B.166.1.2": {"1.B.101.1.2"}}
    current = {"1.B.166", "1.B.166.1", "1.B.166.1.2"}
    assert still_current_moves(raw, current) == ["1.B.166.1.2"]
    aliases, _ = resolve_moved_aliases(raw, h, current=current)
    # 1.B.166.1.2 is current -> no alias; coarse ids then have an unmoved system -> none
    assert aliases == {"1.B.166.1.1": "1.B.101.1.1"}
    assert resolve_moved_aliases(raw, h, current=None)[0]["1.B.166"] == "1.B.101"


def test_current_tcids_from_raw_includes_subfamily_prefixes():
    from multiomics_kg.utils.tcdb_utils import current_tcids
    cur = current_tcids(families={"1.B.101"}, systems={"1.B.101.1.1", "2.A.1.3.5"})
    assert cur == {"1.B.101", "1.B.101.1", "1.B.101.1.1", "2.A.1.3", "2.A.1.3.5"}


def test_resolve_moved_aliases_target_outside_hierarchy_is_skipped():
    from multiomics_kg.utils.tcdb_utils import resolve_moved_aliases
    h = _h("1", "1.B", "1.B.166", "1.B.166.1", "1.B.166.1.1")
    aliases, _ = resolve_moved_aliases({"1.B.166.1.1": {"1.B.101.1.1"}}, h)
    assert aliases == {}


def test_load_moved_overrides_validates(tmp_path):
    import pytest
    from multiomics_kg.utils.tcdb_utils import load_moved_overrides
    h = _h("9", "9.B", "9.B.99", "9.B.99.1", "1", "1.B", "1.B.93", "1.B.93.1")
    p = tmp_path / "o.yaml"
    p.write_text('9.B.99.1:\n  target: 1.B.93.1\n  reason: "MipA genes"\n')
    assert load_moved_overrides(p, h) == {"9.B.99.1": "1.B.93.1"}
    p.write_text('9.B.99.1:\n  target: 1.B.93.1\n  reason: "x"\n  extra: 1\n')
    with pytest.raises(ValueError):
        load_moved_overrides(p, h)
    p.write_text('9.B.99.1:\n  target: 1.B.93.9\n  reason: "x"\n')
    with pytest.raises(ValueError):
        load_moved_overrides(p, h)
    p.write_text('9.B.99.1:\n  target: 1.B.93.1\n  reason: ""\n')
    with pytest.raises(ValueError):
        load_moved_overrides(p, h)
    assert load_moved_overrides(tmp_path / "missing.yaml", h) == {}


def test_committed_moved_overrides_load_against_committed_hierarchy():
    import json
    from pathlib import Path
    from multiomics_kg.utils.tcdb_utils import load_moved_overrides
    root = Path(__file__).resolve().parents[1]
    h = json.loads((root / "cache/data/tcdb/tcdb_hierarchy.json").read_text())
    assert load_moved_overrides(root / "config/tcdb_moved_overrides.yaml", h) == \
        {"9.B.99.1": "1.B.93.1"}


@pytest.mark.parametrize("name, expected", [
    ("Type VI secretion system, TssA - TssM.", ("export", "tcdb_name")),       # secretion is export
    ("General secretory pathway (Sec-SRP) complex.", ("export", "tcdb_name")),
    ("K+ channel, Kir4.2, involved in bile secretion.", ("unknown", None)),     # role phrase
    ("TRPC1 Ca2+ channel controls salivary gland fluid secretion", ("unknown", None)),
    ("DedA family member involved in selenite uptake", ("unknown", None)),
    ("Golgi/secretory granule Zn2+ uptake permease", ("unknown", None)),         # organelle
])
def test_tcdb_direction_secretion_and_role_phrases(name, expected):
    assert tu.tcdb_direction("2.A.7.1.1", _sys(name), {}) == expected


def test_tcdb_display_name_names_unnamed_members_after_nearest_named_ancestor():
    h = {
        "3": {"name": "Primary Active Transporters", "level_kind": "tc_class", "parent": None},
        "3.A": {"name": "", "level_kind": "tc_subclass", "parent": "3"},
        "3.A.1": {"name": "The ABC Superfamily", "level_kind": "tc_family", "parent": "3.A"},
        "3.A.1.17": {"name": "", "level_kind": "tc_subfamily", "parent": "3.A.1"},
        "3.A.1.17.1": {"name": "", "level_kind": "tc_system", "parent": "3.A.1.17"},
        "3.A.1.17.2": {"name": "Taurine porter", "level_kind": "tc_system", "parent": "3.A.1.17"},
    }
    assert tu.tcdb_display_name("3.A.1.17.2", h) == "Taurine porter"
    assert tu.tcdb_display_name("3.A.1.17", h) == "The ABC Superfamily — subfamily 3.A.1.17"
    # skips the unnamed subfamily, lands on the family
    assert tu.tcdb_display_name("3.A.1.17.1", h) == "The ABC Superfamily — system 3.A.1.17.1"
    # only subfamilies / systems borrow a name
    assert tu.tcdb_display_name("3.A", h) == "3.A"
    assert tu.tcdb_display_name("9.Z.9", h) == "9.Z.9"
