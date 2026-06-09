"""Tests for publication-topics extraction utilities (deterministic helpers)."""
import json

import pytest

from multiomics_kg.extraction.topics.extraction_utils import (
    LEGACY_TOPICS_FILENAME,
    get_paper_strains,
    has_topics,
    load_topics,
    migrate_legacy,
    resolve_pdf_path,
    save_topics,
    topics_path,
)


# ── get_paper_strains ──


def test_get_paper_strains_dedup_and_order():
    """Strains are the order-preserving, de-duplicated union of experiment
    organisms plus any treatment_organism."""
    config = {
        "publication": {
            "experiments": {
                "exp1": {"organism": "Prochlorococcus MED4"},
                "exp2": {
                    "organism": "Prochlorococcus MIT9313",
                    "treatment_organism": "Alteromonas macleodii HOT1A3",
                },
                "exp3": {"organism": "Prochlorococcus MED4"},  # duplicate
            }
        }
    }
    assert get_paper_strains(config) == [
        "Prochlorococcus MED4",
        "Prochlorococcus MIT9313",
        "Alteromonas macleodii HOT1A3",
    ]


def test_get_paper_strains_empty():
    assert get_paper_strains({}) == []
    assert get_paper_strains({"publication": {}}) == []
    assert get_paper_strains({"publication": {"experiments": {}}}) == []


def test_get_paper_strains_skips_non_dict_and_blank():
    config = {
        "publication": {
            "experiments": {
                "exp1": "not-a-dict",
                "exp2": {"organism": ""},          # blank ignored
                "exp3": {"organism": "Prochlorococcus MED4"},
            }
        }
    }
    assert get_paper_strains(config) == ["Prochlorococcus MED4"]


# ── resolve_pdf_path ──


def test_resolve_pdf_path_relative_to_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pdf = tmp_path / "paper.pdf"
    pdf.write_bytes(b"%PDF-1.4 fake")
    config = {"publication": {"papermainpdf": "paper.pdf"}}
    assert resolve_pdf_path(config, tmp_path) == pdf


def test_resolve_pdf_path_missing_returns_none(tmp_path):
    config = {"publication": {"papermainpdf": "does_not_exist.pdf"}}
    assert resolve_pdf_path(config, tmp_path) is None


def test_resolve_pdf_path_falls_back_to_dir_glob(tmp_path):
    pdf = tmp_path / "some_paper.pdf"
    pdf.write_bytes(b"%PDF-1.4 fake")
    config = {"publication": {}}  # no papermainpdf
    assert resolve_pdf_path(config, tmp_path) == pdf


# ── save/load/has topics ──


def test_save_and_load_topics_roundtrip(tmp_path):
    genes = [{"gene_name": "ntcA", "identifiers": ["PMM0552"], "strain": "Prochlorococcus MED4",
              "prominence": "central", "surface_form": "NtcA", "evidence": "q"}]
    pathways = [{"surface_form": "nitrogen metabolism", "kegg_pathway_id": "ko00910",
                 "prominence": "central", "evidence": "q"}]
    metadata = {"paper": "Test 2026", "doi": "10.x/y", "uncaptured_identifiers": "some"}

    assert not has_topics(tmp_path)
    out = save_topics(tmp_path, metadata, genes, pathways)
    # Written into the publication_topics/ subfolder.
    assert out == topics_path(tmp_path)
    assert out.parent.name == "publication_topics"
    assert has_topics(tmp_path)

    data = load_topics(tmp_path)
    assert data["metadata"]["paper"] == "Test 2026"
    assert data["genes"] == genes
    assert data["pathways"] == pathways
    assert json.loads(out.read_text(encoding="utf-8"))["genes"][0]["gene_name"] == "ntcA"


def test_load_topics_missing_returns_empty(tmp_path):
    assert load_topics(tmp_path) == {}


def test_legacy_flat_file_is_recognized_and_migrated(tmp_path):
    """A pre-subfolder publication_topics.json is read, then moved into the subfolder."""
    legacy = tmp_path / LEGACY_TOPICS_FILENAME
    legacy.write_text(json.dumps({"metadata": {"paper": "Old"}, "genes": [], "pathways": []}),
                      encoding="utf-8")
    # Recognized + readable from the legacy location.
    assert has_topics(tmp_path)
    assert load_topics(tmp_path)["metadata"]["paper"] == "Old"
    # Migration moves it into the subfolder and removes the flat file.
    assert migrate_legacy(tmp_path) is True
    assert topics_path(tmp_path).exists()
    assert not legacy.exists()
    assert migrate_legacy(tmp_path) is False  # idempotent
