"""Unit tests for multiomics_kg.utils.tool_calls_io (shared per-strain I/O).

Pins the file-naming convention (`.limited_<N>.` infix, indented + sorted JSON)
and the cross-strain iteration so a future schema drift breaks loudly.
"""

import json

import pytest

from multiomics_kg.utils import tool_calls_io as tcio


def test_calls_path_naming(tmp_path):
    p = tcio.calls_path(tmp_path, "interproscan", "MED4")
    assert p == tmp_path / "interproscan" / "MED4.interproscan.calls.json"


def test_calls_path_limited_infix(tmp_path):
    p = tcio.calls_path(tmp_path, "interproscan", "MED4", limited=100)
    assert p.name == "MED4.interproscan.limited_100.calls.json"


def test_skill_summary_path_naming(tmp_path):
    p = tcio.skill_summary_path(tmp_path, "interproscan", "MED4")
    assert p.name == "MED4.interproscan.skill_summary.json"
    p2 = tcio.skill_summary_path(tmp_path, "interproscan", "MED4", limited=50)
    assert p2.name == "MED4.interproscan.limited_50.skill_summary.json"


def test_calls_round_trip(tmp_path):
    calls = {"WP_2": {"match_count": 0}, "WP_1": {"match_count": 3}}
    tcio.save_calls(tmp_path, "interproscan", "MED4", calls)
    assert tcio.load_calls(tmp_path, "interproscan", "MED4") == calls


def test_calls_written_sorted_and_indented(tmp_path):
    tcio.save_calls(tmp_path, "interproscan", "MED4", {"WP_2": 2, "WP_1": 1})
    text = tcio.calls_path(tmp_path, "interproscan", "MED4").read_text()
    assert text.index('"WP_1"') < text.index('"WP_2"')  # sort_keys
    assert "\n  " in text  # indented


def test_summary_round_trip(tmp_path):
    summary = {"strain": "MED4", "calls_made": 5}
    tcio.save_skill_summary(tmp_path, "interproscan", "MED4", summary)
    assert tcio.load_skill_summary(tmp_path, "interproscan", "MED4") == summary


def test_limited_artifacts_are_separate_files(tmp_path):
    tcio.save_calls(tmp_path, "interproscan", "MED4", {"a": 1})
    tcio.save_calls(tmp_path, "interproscan", "MED4", {"b": 2}, limited=10)
    assert tcio.load_calls(tmp_path, "interproscan", "MED4") == {"a": 1}
    limited = tcio.calls_path(tmp_path, "interproscan", "MED4", limited=10)
    assert json.loads(limited.read_text()) == {"b": 2}


def test_iter_strain_calls(tmp_path, monkeypatch):
    # Two fake strains; only one has a calls.json for the tool.
    dir_a = tmp_path / "A"
    dir_b = tmp_path / "B"
    tcio.save_calls(dir_a, "interproscan", "A", {"WP_A": {"match_count": 1}})
    rows = [
        {"strain_name": "A", "data_dir": str(dir_a)},
        {"strain_name": "B", "data_dir": str(dir_b)},
    ]
    monkeypatch.setattr(tcio, "load_genome_rows", lambda strains=None: rows)

    got = list(tcio.iter_strain_calls("interproscan"))
    assert [s for s, _, _ in got] == ["A"]  # B skipped (no calls.json)
    assert got[0][2] == {"WP_A": {"match_count": 1}}


def test_iter_strain_calls_status_yields_none_for_missing(tmp_path, monkeypatch):
    dir_a = tmp_path / "A"
    dir_b = tmp_path / "B"
    tcio.save_calls(dir_a, "interproscan", "A", {"WP_A": 1})
    rows = [
        {"strain_name": "A", "data_dir": str(dir_a)},
        {"strain_name": "B", "data_dir": str(dir_b)},
    ]
    monkeypatch.setattr(tcio, "load_genome_rows", lambda strains=None: rows)

    status = {s: calls for s, _, calls in tcio.iter_strain_calls_status("interproscan")}
    assert status["A"] == {"WP_A": 1}
    assert status["B"] is None
