"""Verify collect_legacy_protein_accessions skips genome-derived sources
(genomic_gca/genomic_gcf/cds_from_genomic.fna) so we don't waste NCBI IPG
fetches on accessions whose locus_tag is already locally known via GFF
coordinates.

Regression: prior to this filter, ~30K INSDC accessions per pipeline run
were fetched, bloating the IPG cache to 100+ MB and adding ~30 min of
network latency.
"""
from __future__ import annotations

from multiomics_kg.download.build_gene_id_mapping import (
    collect_legacy_protein_accessions,
)
from multiomics_kg.download.gene_id_graph import GeneIdGraph


def test_skips_genomic_gca_source():
    """Accessions extracted from genomic_gca.gff are NOT enriched."""
    graph = GeneIdGraph()
    rows = [
        ([("CAE18549.1", "protein_id_refseq")], "genomic_gca/MED4"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == []


def test_skips_genomic_gcf_source():
    graph = GeneIdGraph()
    rows = [
        ([("WP_011131728.1", "protein_id_refseq")], "genomic_gcf/MED4"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == []


def test_skips_cds_from_genomic_fna_source():
    graph = GeneIdGraph()
    rows = [
        ([("NP_892211.1", "protein_id_refseq")], "cds_from_genomic.fna"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == []


def test_includes_paperconfig_sources():
    """Accessions from real paperconfigs ARE enriched."""
    graph = GeneIdGraph()
    rows = [
        ([("NP_892211.1", "protein_id_refseq")], "Biller 2014/s2_med4_vesicle_proteome"),
        ([("NP_894023.1", "protein_id_refseq")], "Biller 2014/s3_mit9313_vesicle_proteome"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == [
        "NP_892211.1", "NP_894023.1",
    ]


def test_skips_when_already_in_graph():
    """Already-known accessions (in specific_lookup or multi_lookup) are skipped."""
    graph = GeneIdGraph()
    graph.specific_lookup["NP_892211.1"] = "PMM0XXX"
    rows = [
        ([("NP_892211.1", "protein_id_refseq")], "Biller 2014/s2_med4_vesicle_proteome"),
        ([("NP_892218.1", "protein_id_refseq")], "Biller 2014/s2_med4_vesicle_proteome"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == ["NP_892218.1"]


def test_skips_non_protein_id_refseq_tokens():
    """Only protein_id_refseq tokens trigger IPG enrichment, not others."""
    graph = GeneIdGraph()
    rows = [
        ([("Q31L36", "uniprot_accession"), ("NP_892211.1", "protein_id_refseq")],
         "Some Paper/table"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == ["NP_892211.1"]


def test_skips_non_legacy_format_tokens():
    """Plain locus tags incorrectly typed as protein_id_refseq are skipped."""
    graph = GeneIdGraph()
    rows = [
        ([("PMM0001", "protein_id_refseq")], "Some Paper/table"),
        ([("TX50_RS00020", "protein_id_refseq")], "Some Paper/table"),
    ]
    assert collect_legacy_protein_accessions(rows, graph) == []


def test_mixed_sources_only_paperconfig_kept():
    """Genome-derived rows in the same list as paperconfig rows: only
    paperconfig accessions are enriched."""
    graph = GeneIdGraph()
    rows = [
        ([("CAE18549.1", "protein_id_refseq")], "genomic_gca/MED4"),       # skip
        ([("WP_011131728.1", "protein_id_refseq")], "genomic_gcf/MED4"),    # skip
        ([("NP_111.1", "protein_id_refseq")], "cds_from_genomic.fna"),      # skip
        ([("NP_892211.1", "protein_id_refseq")], "Biller 2014/s2_med4"),   # keep
    ]
    assert collect_legacy_protein_accessions(rows, graph) == ["NP_892211.1"]
