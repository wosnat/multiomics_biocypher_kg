"""Tests for the tigr_roles.json sub-builder of step 9 (offline: raw files are
written into a tmp cache dir; no network)."""
import json

import pytest

from multiomics_kg.download import build_ncbifam_reference as b


ROLE_NAMES = (
    "role_id:\t132\tmainrole:\tDNA metabolism\n"
    "role_id:\t132\tsub1role:\tDNA replication, recombination, and repair\n"
    "role_id:\t719\tsub1role:\tOrphan\n"
)
ROLE_LINK = "TIGR00001\t132\nTIGR00003\t719\n"


@pytest.fixture
def cache(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "TIGR_ROLE_NAMES").write_text(ROLE_NAMES)
    (raw / "TIGRFAMS_ROLE_LINK").write_text(ROLE_LINK)
    monkeypatch.setattr(b, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(b, "RAW_DIR", raw)
    monkeypatch.setattr(b, "TIGR_ROLE_NAMES_RAW", raw / "TIGR_ROLE_NAMES")
    monkeypatch.setattr(b, "TIGRFAMS_ROLE_LINK_RAW", raw / "TIGRFAMS_ROLE_LINK")
    monkeypatch.setattr(b, "TIGR_ROLES_JSON", tmp_path / "tigr_roles.json")
    return tmp_path


def test_build_tigr_roles_writes_expected_shape(cache):
    out = b.build_tigr_roles(force=True)
    assert out["release"].startswith("TIGRFAMs 15.0")
    assert out["roles"] == {"132": {"mainrole": "DNA metabolism",
                                    "sub1role": "DNA replication, recombination, and repair"}}
    assert out["family_role"] == {"TIGR00001": "132"}
    on_disk = json.loads((cache / "tigr_roles.json").read_text())
    assert on_disk == out


def test_build_tigr_roles_reuses_existing_without_force(cache):
    (cache / "tigr_roles.json").write_text(json.dumps({"release": "x", "roles": {}, "family_role": {}}))
    assert b.build_tigr_roles(force=False) == {"release": "x", "roles": {}, "family_role": {}}


def test_build_tigr_roles_download_failure_reuses_committed(cache, monkeypatch):
    (cache / "tigr_roles.json").write_text(json.dumps({"release": "kept", "roles": {}, "family_role": {}}))
    def boom(url, dest):
        raise RuntimeError("ftp down")
    monkeypatch.setattr(b, "_download", boom)
    out = b.build_tigr_roles(force=True, refetch_raw=True)
    assert out["release"] == "kept"
