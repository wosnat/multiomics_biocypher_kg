"""Unit test for _compute_contributing_sources in build_gene_annotations.

F2 ship 2.4 derives Gene.contributing_sources from existing source-provenance
fields in gene_annotations_merged.json.
"""

import pytest
from multiomics_kg.download.build_gene_annotations import _compute_contributing_sources


def test_ncbi_only_gene_returns_ncbi():
    """A gene with no source-tagged fields beyond NCBI gets ['ncbi']."""
    gene = {"locus_tag": "TX50_RS09500"}
    assert _compute_contributing_sources(gene) == ["ncbi"]


def test_gene_with_eggnog_match():
    """Eggnog hit (seed_ortholog non-null) → eggnog in sources."""
    gene = {
        "locus_tag": "TX50_RS09501",
        "seed_ortholog": "59919.PMM0001",
        "seed_ortholog_evalue": 1e-100,
    }
    sources = _compute_contributing_sources(gene)
    assert "ncbi" in sources
    assert "eggnog" in sources


def test_gene_with_uniprot_match():
    """UniProt accession non-null → uniprot in sources."""
    gene = {
        "locus_tag": "PMM0001",
        "uniprot_accession": "Q7VFD7",
    }
    sources = _compute_contributing_sources(gene)
    assert "uniprot" in sources


def test_gene_with_cyanorak_locus():
    """Cyanorak locus_tag non-null → cyanorak in sources."""
    gene = {
        "locus_tag": "PMM0001",
        "locus_tag_cyanorak": "PMM0001",
    }
    sources = _compute_contributing_sources(gene)
    assert "cyanorak" in sources


def test_full_house_gene():
    """A well-annotated gene shows all 4 sources."""
    gene = {
        "locus_tag": "PMM0001",
        "locus_tag_cyanorak": "PMM0001",
        "uniprot_accession": "Q7VFD7",
        "seed_ortholog": "59919.PMM0001",
        "eggnog_ogs": ["COG0001@1|root"],
    }
    sources = _compute_contributing_sources(gene)
    assert sources == ["cyanorak", "eggnog", "ncbi", "uniprot"]


def test_returns_sorted_list():
    """Output must be a sorted list (deterministic)."""
    gene = {"locus_tag": "X", "uniprot_accession": "Z", "seed_ortholog": "1.A"}
    sources = _compute_contributing_sources(gene)
    assert sources == sorted(sources)


# ─── per-token *_source maps (Phase 1b union attribution) ─────────────────────

from multiomics_kg.download.build_gene_annotations import _has_source_label


def test_has_source_label_reads_per_token_map():
    """A union field's per-token map {token: [sources]} is recognised."""
    gene = {
        "locus_tag": "X",
        "go_terms": ["GO:0003677"],
        "go_terms_source": {"GO:0003677": ["uniprot", "interproscan"]},
    }
    assert _has_source_label(gene, "uniprot")
    assert _has_source_label(gene, "interproscan")
    assert not _has_source_label(gene, "eggnog")


def test_has_source_label_still_reads_scalar_fields():
    """Scalar *_source fields (single-resolver) keep working alongside maps."""
    gene = {"locus_tag": "X", "product_source": "cyanorak"}
    assert _has_source_label(gene, "cyanorak")
    assert not _has_source_label(gene, "ncbi")


def test_contributing_sources_sees_union_cyanorak_go():
    """A Cyanorak-only GO token surfaces cyanorak via the per-token map."""
    gene = {
        "locus_tag": "X",
        "go_terms": ["GO:0009523"],
        "go_terms_source": {"GO:0009523": ["cyanorak"]},
    }
    assert "cyanorak" in _compute_contributing_sources(gene)
