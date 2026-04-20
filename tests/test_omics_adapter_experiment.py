"""Tests for Experiment.compartment emission (Plan 2 Task 1)."""
import pytest
import yaml

from multiomics_kg.adapters.omics_adapter import OMICSAdapter


@pytest.fixture
def stub_pdf_extractor(monkeypatch):
    """PDFPublicationExtractor hits disk/network; stub it out."""
    def _fake_extract(self, pdf_path):
        return {"publication": {"title": "stub", "doi": None}}
    from multiomics_kg.adapters.pdf_publication_extraction import PDFPublicationExtractor
    monkeypatch.setattr(PDFPublicationExtractor, "extract_from_pdf", _fake_extract)


def _write_paperconfig(tmp_path, compartment=None):
    exp = {
        "name": "Test exp",
        "organism": "Prochlorococcus MED4",
        "treatment_condition": "N-limit",
        "control_condition": "replete",
        "omics_type": "RNASEQ",
        "test_type": "DESeq2",
        "treatment_type": ["nitrogen"],
        "background_factors": [],
    }
    if compartment is not None:
        exp["compartment"] = compartment
    config = {
        "publication": {
            "papername": "Test 2024",
            "doi": "10.1234/test.2024",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"exp_a": exp},
            "supplementary_materials": {},
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(config))
    (tmp_path / "fake.pdf").write_bytes(b"")
    return str(config_path)


def test_experiment_compartment_defaults_to_whole_cell(tmp_path, stub_pdf_extractor):
    config_path = _write_paperconfig(tmp_path)  # no compartment in paperconfig
    adapter = OMICSAdapter(config_file=config_path)
    adapter.download_data()
    nodes = adapter.get_nodes()
    exp_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "experiment"]
    assert len(exp_nodes) == 1
    _, _, props = exp_nodes[0]
    assert props["compartment"] == "whole_cell"


def test_experiment_compartment_honours_paperconfig(tmp_path, stub_pdf_extractor):
    config_path = _write_paperconfig(tmp_path, compartment="vesicle")
    adapter = OMICSAdapter(config_file=config_path)
    adapter.download_data()
    nodes = adapter.get_nodes()
    exp_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "experiment"]
    _, _, props = exp_nodes[0]
    assert props["compartment"] == "vesicle"
