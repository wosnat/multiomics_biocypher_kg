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


# ── tcdb_direction (spec §5.1) ───────────────────────────────────────────────

_DIR_H = {
    "3":           {"name": "Primary Active Transporters", "level": 0, "parent": None},
    "3.A":         {"name": "P-P-bond-hydrolysis-driven transporters", "level": 1, "parent": "3"},
    "3.A.1":       {"name": "The ATP-binding Cassette (ABC) Superfamily", "level": 2, "parent": "3.A"},
    "3.A.1.4":     {"name": "The Hydrophobic Amino Acid Uptake Transporter (HAAT) Family",
                    "level": 3, "parent": "3.A.1"},
    "3.A.1.4.4":   {"name": "The high-affinity urea porter", "level": 4, "parent": "3.A.1.4"},
    "3.A.1.122":   {"name": "The Macrolide Exporter (MacB) Family", "level": 3, "parent": "3.A.1"},
    "3.A.1.122.2": {"name": "", "level": 4, "parent": "3.A.1.122"},
    "3.A.1.22":    {"name": "The Nickel Uptake Transporter (NiT) Family", "level": 3, "parent": "3.A.1"},
    "3.A.1.22.2":  {"name": "The ABC exporter, CbiMNQO.", "level": 4, "parent": "3.A.1.22"},
    "3.A.1.27":    {"name": "The HCH Family", "level": 3, "parent": "3.A.1"},
    "3.A.1.27.3":  {"name": "MlaABCDEF", "level": 4, "parent": "3.A.1.27"},
    "3.A.1.120":   {"name": "The (Putative) Drug Resistance ATPase-1 Family", "level": 3, "parent": "3.A.1"},
    "3.A.1.120.6": {"name": "The Uup protein", "level": 4, "parent": "3.A.1.120"},
    "3.A.1.203":   {"name": "The Peroxysomal Fatty Acyl CoA Transporter Family", "level": 3,
                    "parent": "3.A.1"},
    "3.A.1.203.4": {"name": "The BacA porter (selective for the uptake of bleomycin)", "level": 4,
                    "parent": "3.A.1.203"},
    "3.A.1.203.11": {"name": "Putative fatty acid exporter; homodimer.", "level": 4,
                     "parent": "3.A.1.203"},
    "2":           {"name": "Electrochemical Potential-driven Transporters (export uptake)",
                    "level": 0, "parent": None},
    "2.A":         {"name": "Porters (efflux)", "level": 1, "parent": "2"},
    "2.A.6":       {"name": "The Resistance-Nodulation-Cell Division (RND) Efflux Family",
                    "level": 2, "parent": "2.A"},
    "2.A.6.1":     {"name": "", "level": 3, "parent": "2.A.6"},
    "2.A.6.1.1":   {"name": "", "level": 4, "parent": "2.A.6.1"},
    "2.A.9":       {"name": "Uptake and export of things", "level": 2, "parent": "2.A"},
    "2.A.9.1.1":   {"name": "", "level": 4, "parent": "2.A.9"},
    "2.A.10":      {"name": "The Importer family", "level": 2, "parent": "2.A"},
    "2.A.11":      {"name": "The Capsular Polysaccharide Exporter (CPSE) Family",
                    "level": 2, "parent": "2.A"},
    "2.A.12":      {"name": "Drug Exporter-1 Family", "level": 2, "parent": "2.A"},
    "2.A.13":      {"name": "Nucleotide importers", "level": 2, "parent": "2.A"},
    "2.A.14":      {"name": "The importance of being transported", "level": 2, "parent": "2.A"},
    "2.A.15":      {"name": "Sugar efflux family", "level": 2, "parent": "2.A"},
    "2.A.15.1":    {"name": "Uptake and export subfamily", "level": 3, "parent": "2.A.15"},
    "2.A.15.1.1":  {"name": "", "level": 4, "parent": "2.A.15.1"},
    "2.A.6.1.2":   {"name": "Putative arabinose efflux porter, AraJ.", "level": 4, "parent": "2.A.6.1"},
    "2.A.6.1.3":   {"name": "Cu+ efflux pump, CusABCF (may pump ions to the medium)", "level": 4,
                    "parent": "2.A.6.1"},
    "2.A.6.1.4":   {"name": "Rv0205 (may function with a heme uptake transporter)", "level": 4,
                    "parent": "2.A.6.1"},
    "2.A.6.1.5":   {"name": "Organic acid uptake porter, DctA of 8 - 10 putative TMSs.", "level": 4,
                    "parent": "2.A.6.1"},
    "2.A.6.1.6":   {"name": "Ferrous iron exporter, MbfA, of the erythrin-vacuolar iron family", "level": 4,
                    "parent": "2.A.6.1"},
    "2.A.20":      {"name": "The Peroxisomal Protein Importer (PPI) Family", "level": 2, "parent": "2.A"},
    "2.A.20.1.1":  {"name": "", "level": 4, "parent": "2.A.20"},
    "2.A.21":      {"name": "The Uncharacterized DMT-1 (possible Mg2+ uptake porter) Family", "level": 2,
                    "parent": "2.A"},
    "2.A.50":      {"name": "The Glycerol Uptake (GUP) or Membrane-bound Acyl Transferase Family",
                    "level": 2, "parent": "2.A"},
    "2.A.50.1.1":  {"name": "", "level": 4, "parent": "2.A.50"},
}


@pytest.mark.parametrize("tc, expected", [
    ("3.A.1.4.4", ("uptake", "abc_numbering")),
    ("3.A.1.4", ("uptake", "abc_numbering")),
    ("3.A.1.122.2", ("export", "tcdb_name")),  # 3.A.1.100+: export only when a name says so
    ("3.A.1.122", ("export", "tcdb_name")),
    ("3.A.1.120.6", ("unknown", None)),        # ABC-F, no direction word: numbering is not enough
    ("3.A.1.203.4", ("uptake", "tcdb_name")),  # BacA: an importer numbered among exporters
    ("3.A.1.203.11", ("unknown", None)),       # hedged ("Putative ... exporter")
    ("3.A.1.22", ("uptake", "abc_numbering")),
    ("3.A.1.22.2", ("unknown", None)),         # numbering says uptake, name says export
    ("3.A.1.27", ("unknown", None)),           # Mla: _NO_DIRECTION overrides numbering
    ("3.A.1.27.3", ("unknown", None)),
    ("3.A.1", ("unknown", None)),           # superfamily itself: rule 1 starts at subfamilies
    ("2.A.6", ("export", "tcdb_name")),
    ("2.A.6.1.1", ("export", "tcdb_name")),  # inherits the family's name
    ("2.A.9.1.1", ("unknown", None)),        # both word groups in the only usable name
    ("2.A.10", ("uptake", "tcdb_name")),     # agent noun "Importer"
    ("2.A.11", ("export", "tcdb_name")),     # "Exporter"
    ("2.A.12", ("export", "tcdb_name")),     # hyphenated "Exporter-1"
    ("2.A.13", ("uptake", "tcdb_name")),     # plural "importers"
    ("2.A.14", ("unknown", None)),           # whole word: "importance" is not "import"
    ("2.A.15.1.1", ("unknown", None)),       # a name stating both directions stops the walk
    ("2.A", ("unknown", None)),              # subclass names are never used
    ("2", ("unknown", None)),
    ("3", ("unknown", None)),
    ("9.Z.9", ("unknown", None)),            # not in the hierarchy
    ("2.A.6.1.2", ("unknown", None)),        # hedge in the clause: stops, no parent fallback
    ("2.A.6.1.3", ("export", "tcdb_name")),  # "may" after the direction word's clause
    ("2.A.6.1.4", ("unknown", None)),        # "(may function with a heme uptake ...)"
    ("2.A.6.1.5", ("uptake", "tcdb_name")),  # "putative TMSs" is a later clause
    ("2.A.6.1.6", ("unknown", None)),        # organelle word anywhere in the name
    ("2.A.20.1.1", ("unknown", None)),       # organelle import is not uptake
    ("2.A.21", ("unknown", None)),           # "(possible Mg2+ uptake ...)"
    ("2.A.50.1.1", ("unknown", None)),       # _NO_DIRECTION node stops the walk
])
def test_tcdb_direction(tc, expected):
    # The rule alone: the curated overrides file is tested separately below.
    overrides = {k: {"scope": "subtree", "reason": "t"} for k in ("3.A.1.27", "2.A.50")}
    assert tu.tcdb_direction(tc, _DIR_H, overrides) == expected


def test_tcdb_direction_nearest_name_wins():
    h = {"2": {"name": "", "level": 0, "parent": None},
         "2.A": {"name": "", "level": 1, "parent": "2"},
         "2.A.7": {"name": "Drug efflux family", "level": 2, "parent": "2.A"},
         "2.A.7.1": {"name": "Sugar UPTAKE subfamily", "level": 3, "parent": "2.A.7"},
         "2.A.7.1.1": {"name": "", "level": 4, "parent": "2.A.7.1"}}
    assert tu.tcdb_direction("2.A.7.1.1", h, {}) == ("uptake", "tcdb_name")
    assert tu.tcdb_direction("2.A.7", h, {}) == ("export", "tcdb_name")


def _fam(family_name, **systems):
    """Hierarchy: 2 > 2.A > 2.A.7 (family_name) > 2.A.7.1 (unnamed) > 2.A.7.1.<n>."""
    h = {"2": {"name": "", "level": 0, "parent": None},
         "2.A": {"name": "", "level": 1, "parent": "2"},
         "2.A.7": {"name": family_name, "level": 2, "parent": "2.A"},
         "2.A.7.1": {"name": "", "level": 3, "parent": "2.A.7"}}
    for n, (name, desc) in systems.items():
        h[f"2.A.7.1.{n[1:]}"] = {"name": name, "description": desc, "level": 4,
                                 "parent": "2.A.7.1"}
    return h


@pytest.mark.parametrize("name, desc, expected", [
    ("Sugar efflux porter, SetA", "", ("export", "tcdb_name")),
    ("Probable sugar:H+ efflux porter", "", ("unknown", None)),        # hedge across a colon
    ("Putative heavy metal (Me 2+) efflux porter", "", ("unknown", None)),
    ("Lipid efflux homologue of 593 aas", "", ("unknown", None)),       # resemblance after the word
    ("It is not involved in the efflux of drugs", "", ("unknown", None)),  # negation
    ("Heme efflux porter?", "", ("unknown", None)),                     # literal ?
    ("Sugar efflux porter", "Catalyzes both uptake and efflux of sugars.", ("unknown", None)),
    ("Sugar efflux porter", "Its role is under debate.", ("unknown", None)),
    ("Sugar efflux porter", "Sugar uptake porter, SetA. Unrelated detail.", ("unknown", None)),
    ("Sugar efflux porter", "Sugar efflux porter, SetA. Assay used dye uptake.",
     ("export", "tcdb_name")),                                          # later sentence ignored
    ("Sugar efflux porter", "SetA. Its partner, of unknown function, binds it.",
     ("export", "tcdb_name")),                                          # unknown function in desc
    ("", "", ("export", "tcdb_name")),                                  # empty name inherits
    ("Glucose porter of 400 aas", "", ("export", "tcdb_name")),         # plain name inherits
    ("Uncharacterized MFS protein", "", ("unknown", None)),             # hedged name does not
    ("CorC homologue, YfjD", "", ("unknown", None)),
    ("DUF2899 protein", "", ("unknown", None)),
    ("Auxin efflux carrier family member", "", ("unknown", None)),
    ("Mitochondrial ABC protein, ABCB8", "", ("unknown", None)),        # organellar name
    ("Zn uptake (into synaptic vesicles)", "", ("unknown", None)),
    ("Protein import into the ER", "", ("unknown", None)),
])
def test_tcdb_direction_weak_evidence(name, desc, expected):
    h = _fam("The Drug Efflux (DE) Family", s1=(name, desc))
    assert tu.tcdb_direction("2.A.7.1.1", h, {}) == expected


def test_abc_numbering_escapes():
    h = {"3": {"name": "", "level": 0, "parent": None},
         "3.A": {"name": "", "level": 1, "parent": "3"},
         "3.A.1": {"name": "ABC", "level": 2, "parent": "3.A"},
         "3.A.1.5": {"name": "The Peptide Uptake Family", "level": 3, "parent": "3.A.1"},
         "3.A.1.5.1": {"name": "Dipeptide porter, DppABCDF", "level": 4, "parent": "3.A.1.5"},
         "3.A.1.5.2": {"name": "ABC transporter of unknown function", "level": 4,
                       "parent": "3.A.1.5"},
         "3.A.1.5.3": {"name": "Putative drug exporter", "level": 4, "parent": "3.A.1.5"},
         "3.A.1.5.4": {"name": "General amino acid porter", "level": 4, "parent": "3.A.1.5",
                       "description": "General porter (catalyzes both uptake and efflux)."},
         "3.A.1.5.5": {"name": "Peptide uptake and export system", "level": 4,
                       "parent": "3.A.1.5"}}
    d = lambda t: tu.tcdb_direction(t, h, {})
    assert d("3.A.1.5.1") == ("uptake", "abc_numbering")
    assert d("3.A.1.5.2") == ("unknown", None)    # unknown function in its own name
    assert d("3.A.1.5.3") == ("unknown", None)    # a tentative export claim still blocks
    assert d("3.A.1.5.4") == ("unknown", None)    # description: both directions
    assert d("3.A.1.5.5") == ("unknown", None)    # name states both


def test_direction_overrides_scopes():
    h = _fam("The Drug Efflux (DE) Family", s1=("Sugar efflux porter", ""), s2=("", ""),
             s3=("Glucose porter", ""))
    assert tu.tcdb_direction("2.A.7.1.1", h, {"2.A.7.1": {"scope": "subtree", "reason": "t"}})         == ("unknown", None)
    node = {"2.A.7": {"scope": "node", "reason": "t"}}
    assert tu.tcdb_direction("2.A.7", h, node) == ("unknown", None)
    assert tu.tcdb_direction("2.A.7.1.1", h, node) == ("export", "tcdb_name")  # own name
    assert tu.tcdb_direction("2.A.7.1.2", h, node) == ("unknown", None)        # no inheritance
    assert tu.tcdb_direction("2.A.7.1.3", h, node) == ("unknown", None)


def test_direction_overrides_file_valid():
    ov = tu.load_direction_overrides()
    assert len(ov) > 50
    hierarchy = json.loads(Path("cache/data/tcdb/tcdb_hierarchy.json").read_text())
    missing = sorted(k for k in ov if k not in hierarchy)
    assert not missing, f"override ids not in the TCDB hierarchy: {missing}"
    for k, v in ov.items():
        assert v["scope"] in {"node", "subtree"} and v["reason"], k


def test_direction_overrides_rejects_bad_entry(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("2.A.7:\n  scope: everything\n  reason: x\n")
    with pytest.raises(ValueError):
        tu.load_direction_overrides(p)


def test_direction_overrides_real_file_examples():
    hierarchy = json.loads(Path("cache/data/tcdb/tcdb_hierarchy.json").read_text())
    d = lambda t: tu.tcdb_direction(t, hierarchy)
    assert d("3.A.1.4.4") == ("uptake", "abc_numbering")       # urea ABC (MED4 urtA)
    assert d("3.A.1.27.3") == ("unknown", None)                 # Mla
    assert d("2.A.103.1.1") == ("unknown", None)                # FtsW
    assert d("2.A.7.3.55")[0] == "unknown"                      # DME family name not inherited
    assert d("2.A.6.1.4") == ("export", "tcdb_name")            # CusABCF


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
