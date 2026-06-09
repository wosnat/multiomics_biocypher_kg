"""Tests for the publication discuss-edges adapter (Stage 3)."""
import json

import pytest

from multiomics_kg.adapters.publication_topics_adapter import PublicationTopicsAdapter
from multiomics_kg.extraction.topics.extraction_utils import resolved_path


def _write_paperconfig(paper_dir, doi="10.1234/test"):
    pc = paper_dir / "paperconfig.yaml"
    pc.write_text(
        "publication:\n"
        "  papername: Test 2026\n"
        f"  doi: \"{doi}\"\n"
        "  papermainpdf: data/fake.pdf\n",
        encoding="utf-8",
    )
    return pc


def _write_resolved(paper_dir, genes, pathways):
    p = resolved_path(paper_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"metadata": {}, "genes": genes, "pathways": pathways}),
                 encoding="utf-8")


def test_emits_gene_and_pathway_edges(tmp_path):
    pc = _write_paperconfig(tmp_path)
    _write_resolved(
        tmp_path,
        genes=[
            {"locus_tag": "PMM1157", "prominence": "peripheral", "evidence": "ev1"},
            {"locus_tag": None, "prominence": "central", "evidence": "skip"},  # unresolved
        ],
        pathways=[
            {"pathway_id": "ko00010", "prominence": "central", "evidence": "pev"},
            {"pathway_id": None, "prominence": "peripheral", "evidence": "skip"},  # unresolved
        ],
    )
    edges = list(PublicationTopicsAdapter(str(pc)).get_edges())

    gene_edges = [e for e in edges if e[3] == "publication_discusses_gene"]
    path_edges = [e for e in edges if e[3] == "publication_discusses_kegg_pathway"]
    assert len(gene_edges) == 1 and len(path_edges) == 1  # unresolved skipped

    eid, src, tgt, label, props = gene_edges[0]
    assert src == "doi:10.1234/test"
    assert tgt == "ncbigene:PMM1157"
    assert props == {"prominence": "peripheral", "evidence": "ev1"}

    _, _, ptgt, _, pprops = path_edges[0]
    assert ptgt == "kegg.pathway:ko00010"
    assert pprops["prominence"] == "central"


def test_dedup_central_wins(tmp_path):
    pc = _write_paperconfig(tmp_path)
    _write_resolved(
        tmp_path,
        genes=[
            {"locus_tag": "PMM0552", "prominence": "peripheral", "evidence": "a"},
            {"locus_tag": "PMM0552", "prominence": "central", "evidence": "b"},  # same gene
        ],
        pathways=[],
    )
    edges = list(PublicationTopicsAdapter(str(pc)).get_edges())
    gene_edges = [e for e in edges if e[3] == "publication_discusses_gene"]
    assert len(gene_edges) == 1
    assert gene_edges[0][4]["prominence"] == "central"  # central-wins


def test_string_sanitization(tmp_path):
    pc = _write_paperconfig(tmp_path)
    _write_resolved(
        tmp_path,
        genes=[{"locus_tag": "PMM0001", "prominence": "central",
                "evidence": "it's a |risky| quote"}],
        pathways=[],
    )
    edges = list(PublicationTopicsAdapter(str(pc)).get_edges())
    ev = edges[0][4]["evidence"]
    assert "'" not in ev and "|" not in ev
    assert ev == "it^s a risky quote"


def test_no_doi_skips_edges(tmp_path):
    """A paperconfig with no doi and no PDF cache yields no edges (best effort)."""
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("publication:\n  papername: NoDoi 2026\n  papermainpdf: data/none.pdf\n",
                  encoding="utf-8")
    _write_resolved(tmp_path,
                    genes=[{"locus_tag": "PMM0001", "prominence": "central", "evidence": "x"}],
                    pathways=[])
    edges = list(PublicationTopicsAdapter(str(pc)).get_edges())
    assert edges == []


def test_missing_resolved_file_yields_no_edges(tmp_path):
    pc = _write_paperconfig(tmp_path)
    edges = list(PublicationTopicsAdapter(str(pc)).get_edges())
    assert edges == []
