import json

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
