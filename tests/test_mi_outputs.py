import json
from pathlib import Path

import pandas as pd

from tests.metabolite_identity_helpers import cpd
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.outputs import (
    remap_rows, repoint_reason, validate_build, write_all_atomic, write_resolved_tables)
from multiomics_kg.utils.metabolite_identity.paper import Decision, PaperRow, PaperTable, ResolutionResult


def _entry(ids, **kw):
    e = cpd("x", **kw)
    e["all_identifiers"] = ids
    return e


def test_validate_flags_duplicate_identifier_across_roles():
    table = {"chebi:1": _entry(["chebi:1", "chebi:10004"])}
    roles = {"chebi.role:33229": {"parents": [], "all_identifiers": ["chebi.role:33229", "chebi:10004"]}}
    errs = validate_build(table=table, roles=roles, subtree_substrates={}, subtree_roles={},
                          class_trees={}, resolution=ResolutionResult(), chebi_missing_needed=set())
    assert any("chebi:10004" in e for e in errs)


def test_validate_flags_dangling_and_missing_assignment():
    table = {"kegg.compound:C1": dict(_entry(["kegg.compound:C1"]), stereo_form_of=["kegg.compound:C9"])}
    res = ResolutionResult(decisions={("p", "t", 0): Decision("kegg.compound:C1", "name_match", None)})
    errs = validate_build(table=table, roles={}, subtree_substrates={"2.A.1": ["kegg.compound:C7"]},
                          subtree_roles={}, class_trees={}, resolution=res, chebi_missing_needed=set())
    assert any("C9" in e for e in errs) and any("C7" in e for e in errs)
    assert any("stereo_assignment" in e for e in errs)


def test_write_all_atomic(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    write_all_atomic({a: "1", b: "2"})
    assert a.read_text() == "1" and not list(tmp_path.glob("*.tmp"))


def test_write_resolved_tables(tmp_path):
    csv = tmp_path / "m.csv"
    pd.DataFrame({"compound": ["Alanine", "Mystery"], "v": ["1", "2"]}).to_csv(csv, index=False)
    t = PaperTable("p", "t", csv, 0, [PaperRow("p", "t", 0, "Alanine", None, None),
                                      PaperRow("p", "t", 1, "Mystery", None, None)])
    out = write_resolved_tables([t], {("p", "t", 0): Decision("kegg.compound:C00041", "dominant_stereo_form",
                                                               "dominant_form_assumed")})
    text = out[tmp_path / "m_resolved.csv"]
    df = pd.read_csv(__import__("io").StringIO(text), dtype=str, keep_default_na=False)
    assert list(df.columns)[-3:] == ["metabolite_id", "resolution_method", "stereo_assignment"]
    assert df.iloc[1].tolist()[-3:] == ["", "unresolved", ""]
    rep = json.loads(out[tmp_path / "m_resolution_report.json"])
    assert rep["resolved"] == 1 and rep["unresolved"] == 1


def test_repoint_reason_priority():
    ident = IdentityResult({"kegg.compound:C00133": "kegg.compound:C00133"}, {}, {}, {})
    row = PaperRow("p", "t", 0, "Alanine", None, None)
    r = repoint_reason(old_id="kegg.compound:C00133",
                       new=Decision("kegg.compound:C00041", "dominant_stereo_form", "dominant_form_assumed"),
                       row=row, identity=ident, rekeys={}, dropped={}, case_fold_keys=set())
    assert r == "dominant_stereo_form"
    assert repoint_reason(old_id="kegg.compound:C00041", new=Decision("kegg.compound:C00041", "name_match", "x"),
                          row=row, identity=ident, rekeys={}, dropped={}, case_fold_keys=set()) is None


def test_remap_rows_bases():
    table = {"kegg.compound:C00244": _entry(["kegg.compound:C00244", "chebi:14654"]),
             "chebi:17315": _entry(["chebi:17315", "chebi:12961"])}
    ident = IdentityResult({"chebi:14654": "kegg.compound:C00244", "kegg.compound:C00244": "kegg.compound:C00244"},
                           {"kegg.compound:C00244": ["chebi:14654", "kegg.compound:C00244"]},
                           {"kegg.compound:C00244": ["chebi_conjugate"]}, {"chebi:14654": "chebi_conjugate"})
    rows, unaccounted = remap_rows(
        old_ids={"kegg.compound:C00244", "chebi:14654", "chebi:10004", "mnx:MNXM1565340", "chebi:12961", "chebi:999"},
        table=table, identity=ident, rekeys={}, role_of={"chebi:10004": "chebi.role:33229"},
        dropped={"mnx:MNXM1565340": "kegg.compound:C00062"},
        normalize=lambda x: {"chebi:12961": "chebi:17315"}.get(x, x))
    by = {r["old_id"]: r for r in rows}
    assert by["chebi:14654"]["basis"] == "chebi_conjugate"
    assert by["chebi:10004"] == {"old_id": "chebi:10004", "new_id": "", "basis": "removed_role", "role_id": "chebi.role:33229"}
    assert by["mnx:MNXM1565340"]["basis"] == "removed_placeholder"
    assert by["chebi:12961"]["basis"] == "chebi_primary_id"
    assert unaccounted == ["chebi:999"]


import io  # noqa: E402

import pytest  # noqa: E402

from multiomics_kg.utils.metabolite_identity.outputs import build_identity_report, build_mapping_v2  # noqa: E402


def _v(**kw):
    args = dict(table={}, roles={}, subtree_substrates={}, subtree_roles={}, class_trees={},
                resolution=ResolutionResult(), chebi_missing_needed=set())
    args.update(kw)
    return validate_build(**args)


def test_validate_resolution_problems():
    res = ResolutionResult(ambiguous=["foo"], collisions=["bar"], errors=["boom"])
    errs = _v(resolution=res)
    assert any("foo" in e for e in errs) and any("bar" in e for e in errs) and "boom" in errs


def test_validate_duplicate_inchikey():
    t = {"a": dict(_entry(["a"]), inchikey="K"), "b": dict(_entry(["b"]), inchikey="K")}
    assert any("InChIKey K" in e for e in _v(table=t))


def test_validate_class_tree_dangling():
    tree = {"nodes": [{"id": "n1", "parent": "nX"}],
            "members": [{"metabolite_id": "zz", "class_id": "n1"}]}
    errs = _v(class_trees={"br1": tree})
    assert any("nX" in e for e in errs) and any("zz" in e for e in errs)


def test_validate_roles_and_chebi():
    roles = {"r1": {"parents": ["r9"], "all_identifiers": []}}
    errs = _v(roles=roles, subtree_roles={"2.A": ["r7"]}, chebi_missing_needed={"chebi:5"})
    assert any("r9" in e for e in errs) and any("r7" in e for e in errs) and any("chebi:5" in e for e in errs)


def test_write_all_atomic_failure(tmp_path, monkeypatch):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    orig = Path.write_text
    calls = []

    def bad(self, *args, **kw):
        calls.append(self)
        if len(calls) == 2:
            raise OSError("disk")
        return orig(self, *args, **kw)
    monkeypatch.setattr(Path, "write_text", bad)
    with pytest.raises(OSError):
        write_all_atomic({a: "1", b: "2"})
    assert not list(tmp_path.glob("*.tmp")) and not a.exists()


def test_write_resolved_tables_skip_rows(tmp_path):
    csv = tmp_path / "m.csv"
    csv.write_text("preamble one\npreamble two\ncompound,v\nAlanine,1\nGlycine,2\n")
    t = PaperTable("p", "t", csv, 2, [PaperRow("p", "t", 0, "Alanine", None, None),
                                      PaperRow("p", "t", 1, "Glycine", None, None)])
    out = write_resolved_tables([t], {("p", "t", 0): Decision("kegg.compound:C00041", "name_match", "as_reported")})
    text = out[tmp_path / "m_resolved.csv"]
    assert "preamble" not in text
    df = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
    assert df.iloc[0]["compound"] == "Alanine" and df.iloc[0]["metabolite_id"] == "kegg.compound:C00041"


def test_build_mapping_v2():
    row_ok, row_bad = PaperRow("p", "t", 0, "Alanine", None, None), PaperRow("p", "t", 1, "Mystery", None, None)
    t = PaperTable("p", "t", Path("x.csv"), 0, [row_ok, row_bad])
    res = ResolutionResult(decisions={("p", "t", 0): Decision("kegg.compound:C1", "name_match", "as_reported"),
                                      ("p", "t", 1): Decision(None, "unresolved", None)},
                           unresolved=[{"paper": "p", "table": "t", "name": "Mystery"}])
    m = build_mapping_v2(id_aliases={}, resolution=res, overrides={}, tables=[t],
                         compounds_summary={}, per_paper_stats={})
    assert m["schema_version"] == 2 and "specific_lookup" not in m and "multi_lookup" not in m
    for k in ("id_aliases", "name_lookup", "paper_overrides", "decisions", "unresolved", "conflicts"):
        assert k in m
    assert list(m["decisions"]["p"]["t"].values()) == [
        {"metabolite_id": "kegg.compound:C1", "resolution_method": "name_match", "stereo_assignment": "as_reported"}]
    assert m["unresolved"]["p"]["t"] == ["Mystery"]


def test_build_identity_report_sets():
    assert build_identity_report(a={"z", "b"}, b=1) == {"a": ["b", "z"], "b": 1}


def test_repoint_reasons():
    ident = IdentityResult({}, {}, {}, {})
    row = PaperRow("p", "t", 0, "Alanine", None, None)

    def rr(old, method, **kw):
        a = dict(rekeys={}, dropped={}, case_fold_keys=set())
        a.update(kw)
        return repoint_reason(old_id=old, new=Decision("kegg.compound:C2", method, "x"),
                              row=row, identity=ident, **a)
    assert rr("kegg.compound:C1", "name_match", case_fold_keys={row.name_key}) == "case_fold_convergence"
    assert rr("kegg.compound:C1", "kegg_direct") == "direct_id_canonicalised"
    assert rr("mnx:MNXM1", "name_match") == "placeholder_recurated"
    assert rr("chebi:5", "name_match", dropped={"chebi:5": "kegg.compound:C2"}) == "remeasured_elsewhere"
    # old target gone or still written, no identity link to the new node: the name re-resolved
    assert rr("kegg.compound:C1", "name_match") == "name_re_resolved"
    assert rr("kegg.compound:C1", "name_match", written={"kegg.compound:C1", "kegg.compound:C2"}) == "name_re_resolved"
    # merged only when identity says the old id merged into the new node (or was re-keyed into it)
    merged = IdentityResult({"kegg.compound:C1": "kegg.compound:C2"}, {}, {}, {"kegg.compound:C1": "mnx_same_entity"})
    assert repoint_reason(old_id="kegg.compound:C1", new=Decision("kegg.compound:C2", "name_match", "x"),
                          row=row, identity=merged, rekeys={}, dropped={}, case_fold_keys=set()) == "merged"
    assert rr("chebi:9", "name_match", rekeys={"chebi:9": "kegg.compound:C2"}) == "merged"


def test_remap_rows_rekey_and_remeasured():
    table = {"kegg.compound:C1": _entry(["kegg.compound:C1"])}
    ident = IdentityResult({"kegg.compound:C1": "kegg.compound:C1"}, {}, {}, {})
    rows, un = remap_rows(old_ids={"chebi:7", "chebi:8"}, table=table, identity=ident,
                          rekeys={"chebi:7": "kegg.compound:C1"}, role_of={},
                          dropped={"chebi:8": "kegg.compound:C1"}, normalize=lambda x: x)
    by = {r["old_id"]: r for r in rows}
    assert by["chebi:7"]["basis"] == "rekey_chebi_kegg_xref" and by["chebi:7"]["new_id"] == "kegg.compound:C1"
    assert by["chebi:8"]["basis"] == "remeasured_elsewhere" and un == []


def test_harvest_id_col_with_prefix(tmp_path):
    """Ported from step 7: a prefixed id cell keeps its prefix, a bare one gets id_type,
    and placeholder cells ("0", "NA") are no id."""
    from multiomics_kg.utils.metabolite_identity.paper import harvest_paper_tables
    (tmp_path / "m.csv").write_text("compound,KEGG\nX,kegg.compound:C00031\nY,C00031\nZ,0\nW,NA\n")
    cfg = {"publication": {"supplementary_materials": {"tab": {
        "type": "metabolite_assays_table", "filename": "m.csv", "name_col": "compound",
        "id_col": "KEGG", "id_type": "kegg.compound"}}}}
    tables, _ = harvest_paper_tables([(tmp_path / "paperconfig.yaml", cfg)], tmp_path)
    assert [(r.name, r.direct_id) for r in tables[0].rows] == [
        ("X", "kegg.compound:C00031"), ("Y", "kegg.compound:C00031"), ("Z", None), ("W", None)]


def _remap_fallback(old_ids, *, id_index, previous_entries=None, classes=()):
    table = {"kegg.compound:C1": _entry(["kegg.compound:C1"]), "chebi:20": _entry(["chebi:20"]),
             "chebi:30": _entry(["chebi:30"])}
    ident = IdentityResult({n: n for n in table}, {}, {}, {})
    return remap_rows(old_ids=set(old_ids), table=table, identity=ident, rekeys={}, role_of={},
                      dropped={}, normalize=lambda x: {"chebi:3": "chebi:4"}.get(x, x),
                      id_index=id_index, previous_entries=previous_entries or {},
                      is_isotopologue=lambda n: "-43" in (n or ""), is_class=lambda c: c in classes)


def test_remap_fallback_identifier_on_kegg_node_is_rekey():
    rows, un = _remap_fallback({"chebi:1"}, id_index={"chebi:1": "kegg.compound:C1"})
    assert rows == [{"old_id": "chebi:1", "new_id": "kegg.compound:C1", "basis": "rekey_chebi_kegg_xref",
                     "role_id": ""}] and un == []


def test_remap_fallback_identifier_via_normalised_id_on_chebi_node():
    rows, un = _remap_fallback({"chebi:3"}, id_index={"chebi:4": "chebi:20"})
    assert rows[0]["new_id"] == "chebi:20" and rows[0]["basis"] == "chebi_primary_id" and un == []


@pytest.mark.parametrize("prev,classes,basis", [
    ({"name": "calcium-43 atom", "mnxm_id": "MNXM9"}, (), "isotopologue_rejected"),
    ({"name": "mycolate", "mnxm_id": "MNXM9"}, ("chebi:5",), "class_term_rejected"),
    ({"name": "acetyltaurine(1-)", "mnxm_id": "MNXM9"}, (), "mnx_same_entity"),
])
def test_remap_fallback_via_previous_mnx_id(prev, classes, basis):
    rows, un = _remap_fallback({"chebi:5"}, id_index={"mnx:MNXM9": "chebi:30"},
                               previous_entries={"chebi:5": prev}, classes=classes)
    assert rows == [{"old_id": "chebi:5", "new_id": "chebi:30", "basis": basis, "role_id": ""}] and un == []


def test_remap_fallback_no_longer_reachable():
    rows, un = _remap_fallback({"chebi:6"}, id_index={"mnx:MNXM1": "chebi:30"},
                               previous_entries={"chebi:6": {"name": "x", "mnxm_id": "MNXM2"}})
    assert rows == [{"old_id": "chebi:6", "new_id": "", "basis": "no_longer_reachable", "role_id": ""}]
    assert un == []


def test_remap_fallback_ignores_index_target_not_written():
    rows, un = _remap_fallback({"chebi:1"}, id_index={"chebi:1": "kegg.compound:C999"})
    assert rows[0]["basis"] == "no_longer_reachable" and un == []


def test_metabolite_identity_tsvs(tmp_path):
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(
        "mi_tsvs", Path(__file__).resolve().parents[1] / "scripts" / "metabolite_identity_tsvs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    report = {"remap": [{"old_id": "chebi:2", "new_id": "kegg.compound:C1", "basis": "manual", "role_id": ""},
                        {"old_id": "chebi:1", "new_id": "", "basis": "removed_role", "role_id": "chebi.role:9"}],
              "repoints": [{"paper": "p", "table": "t", "paper_name": "Alanine", "old_id": "kegg.compound:C00133",
                            "new_id": "kegg.compound:C00041", "stereo_assignment": "dominant_form_assumed",
                            "reason": "dominant_stereo_form"}]}
    files = mod.write_tsvs(report, tmp_path / "out")
    assert (tmp_path / "out" / "metabolite_id_remap.tsv").read_text() == (
        "old_id\tnew_id\tbasis\trole_id\n"
        "chebi:1\t\tremoved_role\tchebi.role:9\n"
        "chebi:2\tkegg.compound:C1\tmanual\t\n")
    assert files[tmp_path / "out" / "metabolite_assay_repoints.tsv"] == (
        "paper\ttable\tpaper_name\told_id\tnew_id\tstereo_assignment\treason\n"
        "p\tt\tAlanine\tkegg.compound:C00133\tkegg.compound:C00041\tdominant_form_assumed\tdominant_stereo_form\n")


# Final review item 3: an old id whose TCDB substrate string now resolves to a written node
# remaps there (old string → primary map rebuilt with chebi=None), with a ChEBI-term basis.
from tests.metabolite_identity_helpers import make_chebi, term  # noqa: E402
from multiomics_kg.utils.metabolite_identity.crossrefs import is_isotopologue_name  # noqa: E402

_SUB_CHEBI = make_chebi(
    CHEBI_24431=term("chemical entity"),
    CHEBI_33250=term("atom", is_a=["CHEBI:24431"]),
    CHEBI_22984=term("calcium atom", kegg=["C00076"], is_a=["CHEBI:33250"]),
    CHEBI_176566=term("calcium-43 atom", is_a=["CHEBI:22984"]),
    CHEBI_176629=term("cadmium-111(2+)", is_a=["CHEBI:24431"]),
    CHEBI_3556=term("cercosporin", kegg=["C10309"], conj=["CHEBI:146017"], inchikey="X-Y-N"),
    CHEBI_146017=term("cercosporin(2-)", conj=["CHEBI:3556"], inchikey="X-Y-L"),
    CHEBI_15444=term("(1->4)-alpha-D-glucan", alt=["CHEBI:10809"], inchikey="G-H-N"),
    CHEBI_60524=term("mycolates"))


def _substrate_remap(old, string_new, node_entry, *, id_index=None):
    node = "kegg.compound:C9"
    table = {node: node_entry}
    ident = IdentityResult({node: node}, {}, {}, {})
    s = f"CHEBI:{old.split(':')[1]};x"
    return remap_rows(old_ids={old}, table=table, identity=ident, rekeys={}, role_of={}, dropped={},
                      normalize=lambda x: x, id_index=id_index or {node: node}, previous_entries={},
                      is_isotopologue=is_isotopologue_name, is_class=_SUB_CHEBI.is_class,
                      old_string_to_primary={s: old}, string_to_primary={s: string_new},
                      chebi=_SUB_CHEBI)


@pytest.mark.parametrize("old,node_chebi,basis", [
    ("chebi:176566", "22984", "isotopologue_rejected"),     # structural (atom + mass number)
    ("chebi:176629", None, "isotopologue_rejected"),        # by name (charged isotope)
    ("chebi:60524", None, "class_term_rejected"),
    ("chebi:146017", "3556", "chebi_conjugate"),
    ("chebi:10809", "15444", "chebi_primary_id"),
])
def test_remap_via_re_resolved_substrate_string(old, node_chebi, basis):
    e = _entry(["kegg.compound:C9"])
    e["chebi_id"] = node_chebi
    rows, un = _substrate_remap(old, "kegg.compound:C9", e)
    assert rows == [{"old_id": old, "new_id": "kegg.compound:C9", "basis": basis, "role_id": ""}] and un == []


def test_remap_substrate_string_through_id_index():
    e = _entry(["kegg.compound:C9"])
    rows, _ = _substrate_remap("chebi:146017", "kegg.compound:C8", dict(e, chebi_id="3556"),
                               id_index={"kegg.compound:C8": "kegg.compound:C9"})
    assert rows[0]["new_id"] == "kegg.compound:C9" and rows[0]["basis"] == "chebi_conjugate"


def test_remap_substrate_string_unresolved_stays_unreachable():
    rows, _ = _substrate_remap("chebi:146017", "kegg.compound:C7", _entry(["kegg.compound:C9"]))
    assert rows[0]["basis"] == "no_longer_reachable" and rows[0]["new_id"] == ""


def test_resolution_report_paths_relative_to_project_root(tmp_path, monkeypatch):
    """Final review item 5: report paths are repo-relative (raw path outside the repo) and the
    JSON ends with a newline."""
    import multiomics_kg.utils.metabolite_identity.outputs as outputs_mod
    (tmp_path / "data").mkdir()
    csv = tmp_path / "data" / "m.csv"
    pd.DataFrame({"compound": ["Alanine"]}).to_csv(csv, index=False)
    t = PaperTable("p", "t", csv, 0, [PaperRow("p", "t", 0, "Alanine", None, None)])
    monkeypatch.setattr(outputs_mod, "PROJECT_ROOT", tmp_path)
    text = write_resolved_tables([t], {})[tmp_path / "data" / "m_resolution_report.json"]
    rep = json.loads(text)
    assert rep["source_csv"] == "data/m.csv" and rep["resolved_csv"] == "data/m_resolved.csv"
    assert text.endswith("}\n")
    monkeypatch.setattr(outputs_mod, "PROJECT_ROOT", tmp_path / "elsewhere")
    rep = json.loads(write_resolved_tables([t], {})[tmp_path / "data" / "m_resolution_report.json"])
    assert rep["source_csv"] == str(csv)
