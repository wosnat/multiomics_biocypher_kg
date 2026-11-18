"""assembly_properties: NCBI assembly report -> OrganismTaxon properties."""
import copy
import json
from pathlib import Path

import pytest

from multiomics_kg.utils.assembly_report import (
    ASSEMBLY_REPORT_FILENAME,
    assembly_properties,
    load_assembly_report,
)

FIX = Path(__file__).parent / "fixtures" / "assembly_reports"
ACC = {"MED4": "GCF_000011465.1", "MIT9301": "GCF_000015965.1",
       "WH7803": "GCF_000063505.1", "CB0101": "GCF_000179235.1"}


def _load(name):
    return json.loads((FIX / f"{name}.json").read_text())


def _props(name):
    return assembly_properties(_load(name), ACC[name])


def test_med4():
    p = _props("MED4")
    assert p["assembly_level"] == "complete_genome"
    assert p["assembly_status"] == "current"
    assert p["checkm_completeness"] == 99.45
    assert p["checkm_contamination"] == 0.0
    assert p["checkm_marker_set_rank"] == "species"
    assert p["assembly_report_date"] == "2026-11-18"
    assert p["genome_length"] == 1657990 and p["contig_count"] == 1
    assert "assembly_current_accession" not in p
    assert "assembly_suppression_reason" not in p


def test_mit9301_int_completeness():
    p = _props("MIT9301")
    assert p["checkm_completeness"] == 100.0
    assert isinstance(p["checkm_completeness"], float)
    assert isinstance(p["gc_percent"], float)


def test_wh7803_suppressed_no_checkm():
    p = _props("WH7803")
    assert p["assembly_status"] == "suppressed"
    assert p["assembly_suppression_reason"] == "missing strain identifier"
    assert not [k for k in p if k.startswith("checkm_")]


def test_cb0101_previous():
    p = _props("CB0101")
    assert p["assembly_level"] == "contig"
    assert p["assembly_status"] == "previous"
    assert p["assembly_current_accession"] == "GCF_000179235.2"


def test_reason_without_clause():
    w = copy.deepcopy(_load("WH7803"))
    w["report"]["assembly_info"]["suppression_reason"] = (
        "This record was removed by RefSeq staff. Please contact "
        "info@ncbi.nlm.nih.gov for further details.")
    assert assembly_properties(w, ACC["WH7803"])["assembly_suppression_reason"] == \
        "removed by RefSeq staff"


def test_unknown_status_raises():
    w = copy.deepcopy(_load("MED4"))
    w["report"]["assembly_info"]["assembly_status"] = "replaced"
    with pytest.raises(ValueError):
        assembly_properties(w, ACC["MED4"])


def test_unknown_level_raises():
    w = copy.deepcopy(_load("MED4"))
    w["report"]["assembly_info"]["assembly_level"] = "Plasmid Only"
    with pytest.raises(ValueError):
        assembly_properties(w, ACC["MED4"])


def test_load_missing_returns_none(tmp_path):
    assert load_assembly_report(str(tmp_path)) is None


def test_load_reads_the_strain_file(tmp_path):
    (tmp_path / ASSEMBLY_REPORT_FILENAME).write_text((FIX / "MED4.json").read_text())
    assert load_assembly_report(str(tmp_path))["report"]["accession"] == ACC["MED4"]


@pytest.mark.parametrize("name", sorted(ACC))
def test_no_empty_strings(name):
    assert all(v not in ("", None) for v in _props(name).values())
