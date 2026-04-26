"""Unit tests for multiomics_kg.utils.ncbi_protein_xref.

Network-free: exercises the parser and order-based splitter against the
canonical IPG TSV format (sample taken from a live efetch call for
NP_892211.1 + NP_892218.1, two MED4 proteins from Biller 2014 Table S2).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from multiomics_kg.utils.ncbi_protein_xref import (
    IpgEntry,
    _parse_ipg_tsv,
    _split_returned_ipg_by_accession,
    fetch_ipg_xrefs,
    find_assembly_match,
    looks_like_legacy_protein_acc,
)

# Two-group response captured from a live efetch?db=protein&rettype=ipg call.
SAMPLE_TSV = (
    "Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\t"
    "Protein Name\tOrganism\tStrain\tAssembly\n"
    "1441914\tRefSeq\tNC_005072.1\t95155\t96303\t-\tWP_011131728.1\t"
    "trypsin-like peptidase domain-containing protein\t"
    "Prochlorococcus marinus subsp. pastoris str. CCMP1986\tMED4\t"
    "GCF_000011465.1\n"
    "1441914\tINSDC\tBX548174.1\t95155\t96303\t-\tCAE18549.1\t"
    "possible serine protease\t"
    "Prochlorococcus marinus subsp. pastoris str. CCMP1986\tMED4\t"
    "GCA_000011465.1\n"
    "1608989\tRefSeq\tNC_005072.1\t99535\t100965\t-\tWP_011131735.1\t"
    "TolC family protein\t"
    "Prochlorococcus marinus subsp. pastoris str. CCMP1986\tMED4\t"
    "GCF_000011465.1\n"
    "1608989\tINSDC\tBX548174.1\t99535\t100965\t-\tCAE18556.1\t"
    "possible RND family outer membrane efflux protein\t"
    "Prochlorococcus marinus subsp. pastoris str. CCMP1986\tMED4\t"
    "GCA_000011465.1\n"
)


def test_parse_ipg_tsv_extracts_all_rows():
    entries = _parse_ipg_tsv(SAMPLE_TSV)
    assert len(entries) == 4
    assert entries[0] == IpgEntry(
        source_db="RefSeq",
        nucleotide_accession="NC_005072.1",
        start=95155,
        stop=96303,
        strand="-",
        protein_accession="WP_011131728.1",
        protein_name="trypsin-like peptidase domain-containing protein",
        organism="Prochlorococcus marinus subsp. pastoris str. CCMP1986",
        strain="MED4",
        assembly="GCF_000011465.1",
    )


def test_parse_ipg_tsv_handles_empty_or_header_only():
    assert _parse_ipg_tsv("") == []
    assert _parse_ipg_tsv("Id\tSource\n") == []


def test_split_attributes_groups_by_input_order():
    requested = ["NP_892211.1", "NP_892218.1"]
    bodies = _split_returned_ipg_by_accession(requested, SAMPLE_TSV)
    assert set(bodies) == set(requested)
    # First group (1441914) → first input
    assert "1441914" in bodies["NP_892211.1"]
    assert "1608989" not in bodies["NP_892211.1"]
    # Second group (1608989) → second input
    assert "1608989" in bodies["NP_892218.1"]
    assert "1441914" not in bodies["NP_892218.1"]


def test_split_handles_more_inputs_than_groups():
    """A requested accession that produces no IPG group still gets a
    cacheable header-only body."""
    requested = ["NP_892211.1", "NP_892218.1", "NP_INVALID.0"]
    bodies = _split_returned_ipg_by_accession(requested, SAMPLE_TSV)
    assert len(bodies) == 3
    # The third input has no group → header-only
    assert bodies["NP_INVALID.0"].startswith("Id\tSource")
    assert "1441914" not in bodies["NP_INVALID.0"]
    assert "1608989" not in bodies["NP_INVALID.0"]


def test_find_assembly_match_prefers_refseq():
    entries = _parse_ipg_tsv(SAMPLE_TSV)
    # Group 1 has both RefSeq (GCF) and INSDC (GCA) for MED4 — RefSeq wins.
    med4_group = entries[:2]
    hit = find_assembly_match(med4_group, "GCF_000011465.1")
    assert hit is not None
    assert hit.protein_accession == "WP_011131728.1"
    assert hit.source_db == "RefSeq"


def test_find_assembly_match_strips_version():
    entries = _parse_ipg_tsv(SAMPLE_TSV)
    hit = find_assembly_match(entries[:2], "GCF_000011465.99")
    assert hit is not None
    assert hit.protein_accession == "WP_011131728.1"


def test_find_assembly_match_falls_back_to_insdc():
    """When only INSDC (GCA) is present for the target, return that."""
    entries = _parse_ipg_tsv(SAMPLE_TSV)
    insdc_only = [e for e in entries[:2] if e.source_db == "INSDC"]
    hit = find_assembly_match(insdc_only, "GCA_000011465.1")
    assert hit is not None
    assert hit.protein_accession == "CAE18549.1"


def test_find_assembly_match_returns_none_when_no_match():
    entries = _parse_ipg_tsv(SAMPLE_TSV)
    assert find_assembly_match(entries[:2], "GCF_999999999.1") is None


@pytest.mark.parametrize("acc,expected", [
    ("NP_892211.1", True),
    ("WP_011131728.1", True),
    ("YP_001234.1", True),
    ("XP_009876.2", True),
    ("CAE18549.1", True),
    ("KGF88150.1", True),
    ("PMM0001", False),         # locus tag — not a protein accession
    ("TX50_RS00020", False),    # NCBI RS-format locus tag
    ("Q7V3R7", False),          # UniProt accession
    ("", False),
    ("NP_892211", True),        # version-less is also acceptable
])
def test_looks_like_legacy_protein_acc(acc, expected):
    assert looks_like_legacy_protein_acc(acc) is expected


def test_fetch_ipg_xrefs_uses_cache(tmp_path: Path):
    """When all accessions are cached, no network call is made."""
    cache_dir = tmp_path / "ipg"
    cache_dir.mkdir()
    # Pre-populate cache
    (cache_dir / "NP_892211.1.tsv").write_text(SAMPLE_TSV)

    # If this hits the network, it'd raise (network blocked or wrong url),
    # but with cache it never tries.
    result = fetch_ipg_xrefs(["NP_892211.1"], cache_dir=cache_dir)
    assert "NP_892211.1" in result
    # The cached file has 4 rows total (both groups), and that's what we
    # parse — find_assembly_match would still narrow it to MED4 GCF.
    assert len(result["NP_892211.1"]) == 4
