"""Paper-metabolite resolution (spec 2026-10-13 §3.6, §3.6a, §6.5, §6.6, §6.8)."""
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.utils.metabolite_identity.paper import (
    Override, PaperRow, PaperTable, StereoContext, dominant_form, load_paper_overrides,
    resolve_paper_rows)

ALA, DALA, ALA_U = "kegg.compound:C00041", "kegg.compound:C00133", "kegg.compound:C01401"
GLU, DGLU, GLU_U = "kegg.compound:C00025", "kegg.compound:C00217", "kegg.compound:C00302"
GLCNAC, BGLCNAC = "kegg.compound:C00140", "kegg.compound:C03878"
LYS, DLYS = "kegg.compound:C00047", "kegg.compound:C00739"
GSH, GSSG = "kegg.compound:C00051", "kegg.compound:C00127"
FAM = {
    ALA: frozenset({ALA, DALA, ALA_U}), DALA: frozenset({ALA, DALA, ALA_U}), ALA_U: frozenset({ALA, DALA, ALA_U}),
    GLU: frozenset({GLU, DGLU, GLU_U}), DGLU: frozenset({GLU, DGLU, GLU_U}), GLU_U: frozenset({GLU, DGLU, GLU_U}),
    GLCNAC: frozenset({GLCNAC, BGLCNAC}), BGLCNAC: frozenset({GLCNAC, BGLCNAC}),
    LYS: frozenset({LYS, DLYS}), DLYS: frozenset({LYS, DLYS}),
    GSH: frozenset({GSH}), GSSG: frozenset({GSSG})}
COUNTS = {ALA: 520, DALA: 103, ALA_U: 0, GLU: 1905, DGLU: 101, GLU_U: 69, GLCNAC: 65, BGLCNAC: 0,
          LYS: 300, DLYS: 10}


def _ctx(counts=COUNTS):
    return StereoContext(families=FAM, unspecified={ALA_U, GLU_U, GLCNAC},
                         is_anomer=lambda a, b: {a, b} == {GLCNAC, BGLCNAC}, catalyst_counts=counts)


def _table(*rows, paper="p", table="t"):
    return PaperTable(paper, table, None, 0,
                      [PaperRow(paper, table, i, n, d, "kegg_direct" if d else None) for i, (n, d) in enumerate(rows)])


IDX = {x: x for x in FAM} | {"kegg.compound:C99999": ALA}


def _run(tables, *, cands=None, overrides=None, lock=None, ctx=None):
    return resolve_paper_rows(tables, overrides=overrides or {}, name_candidates=cands or {},
                              id_index=IDX, ctx=ctx or _ctx(), lock=lock or {"names": {}, "direct_ids": {}},
                              today="2026-10-15")


def test_unqualified_name_goes_to_dominant_form_and_freezes():
    r = _run([_table(("Alanine", None))], cands={"alanine": {DALA, ALA}})
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method, d.stereo_assignment) == (ALA, "dominant_stereo_form", "dominant_form_assumed")
    assert r.lock["names"]["alanine"]["id"] == ALA
    assert r.lock["names"]["alanine"]["counts"] == {"winner": 520, "runner_up": 103}


def test_qualified_name_stays():
    r = _run([_table(("D-Alanine", None))], cands={"d-alanine": {DALA}})
    assert r.decisions[("p", "t", 0)].metabolite_id == DALA
    assert r.decisions[("p", "t", 0)].stereo_assignment == "as_reported"


def test_exact_match_to_unspecified_node_goes_through_rule():
    r = _run([_table(("Glutamate", None))], cands={"glutamate": {GLU_U}})
    assert r.decisions[("p", "t", 0)].metabolite_id == GLU


def test_unspecified_node_can_win():
    counts = dict(COUNTS, **{GLU_U: 9000})
    r = _run([_table(("Glutamate", None))], cands={"glutamate": {GLU_U}}, ctx=_ctx(counts))
    d = r.decisions[("p", "t", 0)]
    assert d.metabolite_id == GLU_U and d.stereo_assignment == "dominant_form_assumed"


def test_margin_below_2x_is_ambiguous():
    counts = dict(COUNTS, **{DALA: 300})
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, ctx=_ctx(counts))
    assert r.ambiguous and ("p", "t", 0) not in r.decisions


def test_single_specified_match_is_assumed():
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}})
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method, d.stereo_assignment) == (LYS, "name_match", "dominant_form_assumed")


def test_anomer_only_family_keeps_target_as_reported():
    r = _run([_table(("N-Acetylglucosamine", GLCNAC))])
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.stereo_assignment) == (GLCNAC, "as_reported")


def test_direct_id_canonicalised():
    r = _run([_table(("Whatever", "kegg.compound:C99999"))])
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA


def test_stereo_unspecified_direct_id_goes_through_rule():
    r = _run([_table(("Alanine", ALA_U))])
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.stereo_assignment) == (ALA, "dominant_form_assumed")
    assert r.lock["direct_ids"][ALA_U]["id"] == ALA


def test_override_beats_direct_id():
    ov = {"p": {"glutathione": Override(GSH, "as_reported", "GSH, reduced", "Glutathione")}}
    r = _run([_table(("Glutathione", GSSG))], overrides=ov)
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method) == (GSH, "alias_override")


def test_lock_wins_and_drift_reported_not_applied():
    lock = {"names": {"alanine": {"id": DALA, "frozen_on": "2026-01-01",
                                  "counts": {"winner": 1, "runner_up": 0}, "runner_up_id": ALA}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock)
    assert r.decisions[("p", "t", 0)].metabolite_id == DALA
    assert r.drift and r.drift[0]["key"] == "alanine"


def test_lock_target_remapped_through_aliases():
    lock = {"names": {"alanine": {"id": "kegg.compound:C99999", "frozen_on": "x",
                                  "counts": {"winner": 1, "runner_up": 0}, "runner_up_id": None}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock)
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA and r.lock_remapped


def test_stale_lock_entry_reported_not_deleted():
    lock = {"names": {"sarcosine": {"id": ALA, "frozen_on": "x", "counts": {}, "runner_up_id": None}},
            "direct_ids": {}}
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}}, lock=lock)
    assert r.stale_lock and "sarcosine" in r.lock["names"]


def test_unmatched_is_unresolved_not_fatal():
    r = _run([_table(("s-(1-2-dicarboxyethyl)glutathione pos", None))])
    assert r.unresolved and not r.ambiguous and ("p", "t", 0) not in r.decisions


def test_distinct_names_same_target_collide_unless_both_overridden():
    t = _table(("Glutathione", GSSG), ("Glutathione Disulfide", GSSG))
    r = _run([t])
    assert r.collisions
    ov = {"p": {"glutathione": Override(GSSG, "as_reported", "", "Glutathione"),
                "glutathione disulfide": Override(GSSG, "as_reported", "", "Glutathione Disulfide")}}
    assert not _run([t], overrides=ov).collisions


def test_case_variants_are_one_key_not_a_collision():
    r = _run([_table(("Glutamic acid", None), ("glutamic acid", None))], cands={"glutamic acid": {GLU, DGLU, GLU_U}})
    assert not r.collisions and r.decisions[("p", "t", 1)].metabolite_id == GLU


def test_two_papers_map_one_name_differently():
    ov = {"a": {"glutathione": Override(GSH, "as_reported", "", "Glutathione")},
          "b": {"glutathione": Override(GSSG, "as_reported", "", "Glutathione")}}
    r = _run([_table(("Glutathione", None), paper="a"), _table(("Glutathione", None), paper="b")], overrides=ov)
    assert r.decisions[("a", "t", 0)].metabolite_id == GSH and r.decisions[("b", "t", 0)].metabolite_id == GSSG


def test_legacy_override_on_stereo_target_needs_assignment():
    ov = {"p": {"alanine": Override(ALA, None, "", "Alanine")}}
    r = _run([_table(("Alanine", None))], overrides=ov)
    assert r.errors and "stereo_assignment" in r.errors[0]


def test_stale_alias_override_reported():
    ov = {"p": {"cystine": Override(GSSG, "as_reported", "", "Cystine")}}
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}}, overrides=ov)
    assert r.stale_overrides == [{"paper": "p", "name_key": "cystine", "id": GSSG}]


def test_unresolvable_override_target_fails():
    ov = {"p": {"lysine": Override("kegg.compound:C00000", "as_reported", "", "Lysine")}}
    r = _run([_table(("Lysine", None))], overrides=ov)
    assert r.errors


def test_dominant_form_tie_refused():
    counts = dict(COUNTS, **{DALA: 520})
    winner, info = dominant_form({ALA, DALA}, _ctx(counts))
    assert winner is None and info["counts"] == {"winner": 520, "runner_up": 520}


def test_load_overrides_union_and_key_collision(tmp_path):
    (tmp_path / "a.yaml").write_text('Glutathione: {id: "kegg.compound:C00051", stereo_assignment: as_reported, reason: GSH}\n'
                                     'CO2: "kegg.compound:C00011"\nCo2: "kegg.compound:C00175"\n')
    cfg = {"publication": {"supplementary_materials": {
        "t1": {"type": "metabolite_assays_table", "aliases_file": "a.yaml"},
        "t2": {"type": "metabolite_assays_table"}}}}
    ov = load_paper_overrides(tmp_path / "paperconfig.yaml", cfg)
    assert ov["glutathione"].id == "kegg.compound:C00051"
    assert ov["CO2"].id != ov["Co2"].id
    (tmp_path / "a.yaml").write_text('Glutathione: "kegg.compound:C00051"\nglutathione: "kegg.compound:C00127"\n')
    with pytest.raises(ValueError):
        load_paper_overrides(tmp_path / "paperconfig.yaml", cfg)


def test_overrides_apply_to_every_table_of_the_paper():
    """Review focus 1: Kujawinski references its alias file from 1 of 7 tables."""
    ov = {"p": {"glutathione": Override(GSH, "as_reported", "", "Glutathione")}}
    r = _run([_table(("Glutathione", GSSG), table="kegg_table"), _table(("Glutathione", None), table="flags")],
             overrides=ov)
    assert r.decisions[("p", "kegg_table", 0)].metabolite_id == GSH
    assert r.decisions[("p", "flags", 0)].metabolite_id == GSH


def _stub(tmp_path: Path) -> sqlite3.Connection:
    """Tiny MNX-schema resolver (copied from test_build_kegg_metabolism_xrefs_metabolomics) + CO rows."""
    conn = sqlite3.connect(str(tmp_path / "stub.db"))
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT, mass REAL, inchikey TEXT);
        CREATE TABLE compound_aliases (mnxm_id TEXT, source TEXT, value TEXT);
        CREATE TABLE compound_names (mnxm_id TEXT, name_normalized TEXT);
        INSERT INTO compounds VALUES ('MNXM41', 'D-Glucose', 'C6H12O6', 180.06, 'WQZGKKKJIJFFOK-GASJEMHNSA-N');
        INSERT INTO compounds VALUES ('MNXM58', 'Pyruvate', 'C3H3O3', 87.03, 'LCTONWCANYUPML-UHFFFAOYSA-M');
        INSERT INTO compound_aliases VALUES ('MNXM41', 'kegg.compound', 'C00031');
        INSERT INTO compound_aliases VALUES ('MNXM58', 'kegg.compound', 'C00022');
        INSERT INTO compound_names VALUES ('MNXM41', 'glucose');
        INSERT INTO compound_names VALUES ('MNXM58', 'pyruvate');
        INSERT INTO compound_names VALUES ('MNXM13', 'co');
        INSERT INTO compound_aliases VALUES ('MNXM13', 'name', 'CO');
    """)
    conn.commit()
    return conn


def test_resolver_skips_lowercased_step_for_formula_names(tmp_path):
    from multiomics_kg.utils.metabolite_utils import resolve_metabolite_candidates
    conn = _stub(tmp_path)          # helper defined above in this file
    assert resolve_metabolite_candidates("CO", conn) == ["MNXM13"]      # exact alias hit
    assert resolve_metabolite_candidates("Co", conn) == []              # no alias, name step skipped
    assert resolve_metabolite_candidates("glucose", conn) == ["MNXM41"]


def test_lock_entry_with_vanished_target_is_stale_kept_and_not_applied():
    """A lock target that no longer exists is stale_lock_entry, kept in the lock, never used."""
    lock = {"names": {"alanine": {"id": "kegg.compound:C00000", "frozen_on": "x",
                                  "counts": {"winner": 1, "runner_up": 0}, "runner_up_id": None}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock)
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA
    assert [s["key"] for s in r.stale_lock] == ["alanine"]
    assert r.lock["names"]["alanine"]["id"] == "kegg.compound:C00000"


def test_dominant_form_empty_candidates_refused():
    assert dominant_form(set(), _ctx())[0] is None


def test_harvest_paper_tables(tmp_path):
    from multiomics_kg.utils.metabolite_identity.paper import harvest_paper_tables
    d = tmp_path / "data" / "P 2026"
    d.mkdir(parents=True)
    (d / "m.csv").write_text("junk line\ncompound,KEGG\nAlanine,C00041\nGlutathione,0\n,C00001\n")
    (d / "metabolite_aliases.yaml").write_text('Glutathione: {id: "kegg.compound:C00051", '
                                               'stereo_assignment: as_reported, reason: GSH}\n')
    cfg = {"publication": {"supplementary_materials": {"tab": {
        "type": "metabolite_assays_table", "filename": "data/P 2026/m.csv", "skip_rows": 1,
        "id_col": "KEGG", "id_type": "kegg.compound", "aliases_file": "metabolite_aliases.yaml"}}}}
    tables, ovs = harvest_paper_tables([(d / "paperconfig.yaml", cfg)], tmp_path)
    assert [t.key for t in tables] == ["data/P 2026/paperconfig.yaml|tab"]
    rows = tables[0].rows
    assert [(r.name, r.direct_id, r.id_method) for r in rows] == [
        ("Alanine", "kegg.compound:C00041", "kegg_direct"), ("Glutathione", None, None)]
    assert ovs["data/P 2026/paperconfig.yaml"]["glutathione"].id == "kegg.compound:C00051"


def test_lock_drift_reported_when_margin_falls_below_2x():
    lock = {"names": {"alanine": {"id": ALA, "frozen_on": "x",
                                  "counts": {"winner": 520, "runner_up": 103}, "runner_up_id": DALA}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock,
             ctx=_ctx(dict(COUNTS, **{DALA: 300})))
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA
    assert r.drift and r.drift[0]["recomputed"] is None
    assert not _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock).drift


# ---- Fix round 1 ----------------------------------------------------------------------------

def _two_runs(tables, cands, ctx=None):
    r1 = _run(tables, cands=cands, ctx=ctx)
    r2 = _run(tables, cands=cands, ctx=ctx, lock=r1.lock)
    return r1, r2


@pytest.mark.parametrize("tables,cands,ctx", [
    ([_table(("Glutamate", None))], {"glutamate": {GLU_U}}, _ctx(dict(COUNTS, **{GLU_U: 9000}))),
    ([_table(("Alanine", None))], {"alanine": {ALA, DALA}}, None),
    ([_table(("Alanine", ALA_U))], {}, None),
], ids=["unspecified_node_wins", "alanine_family", "unspecified_direct_id"])
def test_decisions_identical_across_builds_with_lock(tables, cands, ctx):
    r1, r2 = _two_runs(tables, cands, ctx)
    assert r1.decisions and r1.decisions == r2.decisions
    assert not r2.drift and r1.lock == r2.lock


def test_name_matching_both_anomers_keeps_anomer_unspecified_node():
    r = _run([_table(("N-Acetylglucosamine", None))], cands={"n-acetylglucosamine": {GLCNAC, BGLCNAC}},
             ctx=_ctx(dict(COUNTS, **{BGLCNAC: 200})))
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.stereo_assignment) == (GLCNAC, "as_reported")
    assert not r.ambiguous and not r.lock["names"]


def test_dominant_form_never_picks_anomer_of_candidate():
    winner, _ = dominant_form({GLCNAC, BGLCNAC}, _ctx(dict(COUNTS, **{BGLCNAC: 200})))
    assert winner == GLCNAC


def _alias_cfg():
    return {"publication": {"supplementary_materials": {
        "t1": {"type": "metabolite_assays_table", "aliases_file": "a.yaml"}}}}


@pytest.mark.parametrize("text,needle", [
    ('NO: "kegg.compound:C00533"\n', "quote the key"),
    ('- Glutathione\n', "mapping"),
    ('Glutathione: {stereo_assignment: as_reported}\n', "no 'id'"),
], ids=["yaml11_bool_key", "non_mapping", "dict_without_id"])
def test_load_overrides_rejects_bad_alias_files(tmp_path, text, needle):
    (tmp_path / "a.yaml").write_text(text)
    with pytest.raises(ValueError, match=needle):
        load_paper_overrides(tmp_path / "paperconfig.yaml", _alias_cfg())


def test_missing_aliases_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="a.yaml"):
        load_paper_overrides(tmp_path / "paperconfig.yaml", _alias_cfg())


@pytest.mark.parametrize("extra", [{"name_col": "Name"}, {"id_col": "KEGG_ID", "id_type": "kegg.compound"}])
def test_harvest_rejects_missing_columns(tmp_path, extra):
    from multiomics_kg.utils.metabolite_identity.paper import harvest_paper_tables
    (tmp_path / "m.csv").write_text("compound,KEGG\nAlanine,C00041\n")
    cfg = {"publication": {"supplementary_materials": {"tab": dict(
        {"type": "metabolite_assays_table", "filename": "m.csv"}, **extra)}}}
    with pytest.raises(ValueError, match="not in"):
        harvest_paper_tables([(tmp_path / "paperconfig.yaml", cfg)], tmp_path)


def test_name_lookup_excludes_overridden_names():
    ov = {"p": {"glutathione": Override(GSH, "as_reported", "", "Glutathione")}}
    r = _run([_table(("Glutathione", None), ("Lysine", None))],
             cands={"glutathione": {GSSG}, "lysine": {LYS}}, overrides=ov)
    assert r.name_lookup == {"lysine": LYS}
