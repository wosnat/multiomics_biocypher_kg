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


def _write_paperconfig(tmp_path, compartment=None, **exp_overrides):
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
    exp.update(exp_overrides)
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


# ─── Coculture partners: list property + one role-tagged edge per partner ───

HOT1A3 = {"organism": "Alteromonas macleodii HOT1A3", "taxid": 28108,
          "assembly_accession": "GCF_001578515.1"}
MARINO = {"organism": "Marinobacter (MarRef v6)", "assembly_accession": "GCF_000166295.1"}


def _build(tmp_path, **exp_overrides):
    adapter = OMICSAdapter(config_file=_write_paperconfig(tmp_path, **exp_overrides))
    adapter.download_data()
    exp_props = next(p for _, lbl, p in adapter.get_nodes() if lbl == "experiment")
    partner_edges = [
        (eid, src, tgt, props) for eid, src, tgt, lbl, props in adapter.get_edges()
        if lbl == "tests_coculture_with"
    ]
    return exp_props, partner_edges


class TestCoculturePartners:
    def test_treatment_partner_edge_is_role_treatment_and_list_has_its_name(self, tmp_path, stub_pdf_extractor):
        props, edges = _build(
            tmp_path, treatment_type=["coculture"], background_factors=["light"],
            treatment_organism=HOT1A3["organism"], treatment_taxid=HOT1A3["taxid"],
            treatment_assembly_accession=HOT1A3["assembly_accession"],
        )
        assert props["coculture_partner"] == ["Alteromonas macleodii HOT1A3"]
        assert len(edges) == 1
        _, _, target, eprops = edges[0]
        assert target == "insdc.gcf:GCF_001578515.1"
        assert eprops["partner_role"] == "treatment"

    def test_phage_partner_is_role_treatment(self, tmp_path, stub_pdf_extractor):
        props, edges = _build(
            tmp_path, treatment_type=["viral"], background_factors=["light"],
            treatment_organism="Phage", treatment_taxid=10239,
        )
        assert props["coculture_partner"] == ["Phage"]
        assert edges[0][2] == "ncbitaxon:10239"
        assert edges[0][3]["partner_role"] == "treatment"

    def test_treatment_organism_outside_coculture_or_viral_is_role_background(self, tmp_path, stub_pdf_extractor):
        """Biller 2016 shape: within-coculture time contrast, partner is the constant."""
        props, edges = _build(
            tmp_path, treatment_type=["growth_phase"], background_factors=["light", "coculture"],
            treatment_organism="Prochlorococcus NATL2A", treatment_taxid=59920,
            treatment_assembly_accession="GCF_000012465.1",
        )
        assert props["coculture_partner"] == ["Prochlorococcus NATL2A"]
        assert edges[0][3]["partner_role"] == "background"

    def test_background_organisms_give_background_edges(self, tmp_path, stub_pdf_extractor):
        props, edges = _build(
            tmp_path, treatment_type=["nitrogen"], background_factors=["coculture", "light"],
            background_organisms=[HOT1A3],
        )
        assert props["coculture_partner"] == ["Alteromonas macleodii HOT1A3"]
        assert [(e[2], e[3]["partner_role"]) for e in edges] == [
            ("insdc.gcf:GCF_001578515.1", "background")
        ]

    def test_two_background_partners_sorted_list_and_distinct_edge_ids(self, tmp_path, stub_pdf_extractor):
        """Moreno 2023 shape: a non-axenic culture with two identified heterotrophs."""
        props, edges = _build(
            tmp_path, treatment_type=["carbon"], background_factors=["coculture", "light"],
            background_organisms=[MARINO, HOT1A3],
        )
        assert props["coculture_partner"] == ["Alteromonas macleodii HOT1A3", "Marinobacter (MarRef v6)"]
        assert len(edges) == 2
        assert len({e[0] for e in edges}) == 2, "edge ids must be distinct per partner"
        assert {e[2] for e in edges} == {"insdc.gcf:GCF_000166295.1", "insdc.gcf:GCF_001578515.1"}
        assert {e[3]["partner_role"] for e in edges} == {"background"}

    def test_treatment_and_background_partners_combine(self, tmp_path, stub_pdf_extractor):
        props, edges = _build(
            tmp_path, treatment_type=["coculture"], background_factors=["coculture", "light"],
            treatment_organism="Prochlorococcus MED4", treatment_taxid=59919,
            treatment_assembly_accession="GCF_000011465.1",
            background_organisms=[MARINO],
        )
        assert props["coculture_partner"] == ["Marinobacter (MarRef v6)", "Prochlorococcus MED4"]
        roles = {e[2]: e[3]["partner_role"] for e in edges}
        assert roles == {"insdc.gcf:GCF_000011465.1": "treatment",
                         "insdc.gcf:GCF_000166295.1": "background"}

    def test_axenic_experiment_has_empty_list_and_no_edges(self, tmp_path, stub_pdf_extractor):
        props, edges = _build(tmp_path, background_factors=["axenic"])
        assert props["coculture_partner"] == []
        assert edges == []

    def test_one_organism_under_two_spellings_is_one_treatment_edge(self, tmp_path, stub_pdf_extractor):
        """Partners dedupe by OrganismTaxon id, not by display name."""
        props, edges = _build(
            tmp_path, treatment_type=["coculture"], background_factors=["coculture", "light"],
            treatment_organism=HOT1A3["organism"],
            treatment_assembly_accession=HOT1A3["assembly_accession"],
            background_organisms=[{"organism": "Alteromonas HOT1A3",
                                   "assembly_accession": HOT1A3["assembly_accession"]}],
        )
        assert [(e[2], e[3]["partner_role"]) for e in edges] == [
            ("insdc.gcf:GCF_001578515.1", "treatment")
        ]
        assert props["coculture_partner"] == ["Alteromonas macleodii HOT1A3"]


def test_empty_supplementary_materials_key_keeps_structural_edges(tmp_path, stub_pdf_extractor):
    """`supplementary_materials:` left empty reads as None — the experiment still gets its
    Has_experiment and partner edges."""
    config_path = _write_paperconfig(
        tmp_path, treatment_type=["coculture"], background_factors=["light"],
        treatment_organism=HOT1A3["organism"],
        treatment_assembly_accession=HOT1A3["assembly_accession"])
    cfg = yaml.safe_load(open(config_path))
    cfg["publication"]["supplementary_materials"] = None
    open(config_path, "w").write(yaml.dump(cfg))
    adapter = OMICSAdapter(config_file=config_path)
    adapter.download_data()
    labels = sorted(lbl for _, _, _, lbl, _ in adapter.get_edges())
    assert labels == ["has_experiment", "tests_coculture_with"]
