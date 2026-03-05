"""Unit tests for scripts/map_img_to_ncbi_proteins.py.

Tests the pure functions used for cross-assembly protein sequence bridging.
Diamond-dependent tests are skipped if diamond is not installed.
"""

import csv
import shutil
import textwrap
from pathlib import Path

import pytest

# Import functions from the script
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from map_img_to_ncbi_proteins import (
    load_gene_mapping,
    map_by_exact_match,
    map_by_subsequence,
    parse_fasta,
    parse_img_gff_gene_ids,
    remap_img_headers,
)

# Conditionally import map_by_diamond (needs diamond binary for full test)
from map_img_to_ncbi_proteins import map_by_diamond


# --- parse_fasta ---


class TestParseFasta:
    def test_single_sequence(self, tmp_path):
        faa = tmp_path / "test.faa"
        faa.write_text(">prot1 some description\nMKTLL\nAVGK\n")
        result = parse_fasta(faa)
        assert result == {"prot1": "MKTLLAVGK"}

    def test_multiple_sequences(self, tmp_path):
        faa = tmp_path / "test.faa"
        faa.write_text(">prot1\nMKTLL\n>prot2\nAVGKM\nLLAA\n>prot3\nM\n")
        result = parse_fasta(faa)
        assert result == {"prot1": "MKTLL", "prot2": "AVGKMLLAA", "prot3": "M"}

    def test_header_takes_first_token(self, tmp_path):
        faa = tmp_path / "test.faa"
        faa.write_text(">AEZ55_0001 product=hypothetical\nMKTLL\n")
        result = parse_fasta(faa)
        assert "AEZ55_0001" in result

    def test_empty_file(self, tmp_path):
        faa = tmp_path / "test.faa"
        faa.write_text("")
        result = parse_fasta(faa)
        assert result == {}

    def test_img_oid_headers(self, tmp_path):
        """IMG FASTA uses numeric OIDs as first token, gene ID as second."""
        faa = tmp_path / "test.faa"
        faa.write_text(">2785712097 altEZ55_00001 product=DnaA\nMKTLL\n")
        result = parse_fasta(faa)
        assert "2785712097" in result  # Takes first token


# --- load_gene_mapping ---


class TestLoadGeneMapping:
    def test_basic_mapping(self, tmp_path):
        csv_path = tmp_path / "gene_mapping.csv"
        csv_path.write_text(
            "locus_tag,protein_id,gene_name\n"
            "EZ55_00001,WP_123456.1,dnaA\n"
            "EZ55_00002,WP_789012.1,dnaN\n"
        )
        result = load_gene_mapping(csv_path)
        assert result == {"WP_123456.1": "EZ55_00001", "WP_789012.1": "EZ55_00002"}

    def test_skips_empty_protein_id(self, tmp_path):
        csv_path = tmp_path / "gene_mapping.csv"
        csv_path.write_text(
            "locus_tag,protein_id\n" "EZ55_00001,WP_123456.1\n" "EZ55_00002,\n"
        )
        result = load_gene_mapping(csv_path)
        assert len(result) == 1

    def test_skips_empty_locus_tag(self, tmp_path):
        csv_path = tmp_path / "gene_mapping.csv"
        csv_path.write_text(
            "locus_tag,protein_id\n" ",WP_123456.1\n" "EZ55_00002,WP_789012.1\n"
        )
        result = load_gene_mapping(csv_path)
        assert len(result) == 1
        assert result["WP_789012.1"] == "EZ55_00002"


# --- map_by_exact_match ---


class TestMapByExactMatch:
    def test_unique_match(self):
        img_seqs = {"AEZ55_0001": "MKTLLAVGK", "AEZ55_0002": "MDNAKRFYQ"}
        ncbi_seqs = {"WP_111.1": "MKTLLAVGK", "WP_222.1": "MDNAKRFYQ"}
        protein_to_lt = {"WP_111.1": "EZ55_00001", "WP_222.1": "EZ55_00002"}

        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {"AEZ55_0001": "EZ55_00001", "AEZ55_0002": "EZ55_00002"}
        assert len(unmatched) == 0
        assert used_ncbi == {"WP_111.1", "WP_222.1"}

    def test_no_match(self):
        img_seqs = {"AEZ55_0001": "MKTLLAVGK"}
        ncbi_seqs = {"WP_111.1": "COMPLETELYDIFFERENT"}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {}
        assert unmatched == {"AEZ55_0001"}

    def test_ambiguous_paralogs_single_locus_tag(self):
        """Two NCBI proteins with identical sequence but one locus_tag → resolves."""
        img_seqs = {"AEZ55_0001": "MKTLLAVGK"}
        ncbi_seqs = {"WP_111.1": "MKTLLAVGK", "WP_222.1": "MKTLLAVGK"}
        # Only one has a locus_tag in gene_mapping
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {"AEZ55_0001": "EZ55_00001"}

    def test_ambiguous_paralogs_multiple_locus_tags(self):
        """Two NCBI proteins with identical sequence AND different locus_tags → ambiguous."""
        img_seqs = {"AEZ55_0001": "MKTLLAVGK"}
        ncbi_seqs = {"WP_111.1": "MKTLLAVGK", "WP_222.1": "MKTLLAVGK"}
        protein_to_lt = {"WP_111.1": "EZ55_00001", "WP_222.1": "EZ55_00002"}

        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {}
        assert "AEZ55_0001" in unmatched

    def test_protein_not_in_gene_mapping(self):
        """NCBI protein exists but has no locus_tag in gene_mapping."""
        img_seqs = {"AEZ55_0001": "MKTLLAVGK"}
        ncbi_seqs = {"WP_111.1": "MKTLLAVGK"}
        protein_to_lt = {}  # No mapping

        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {}
        assert "AEZ55_0001" in unmatched


# --- map_by_subsequence ---


class TestMapBySubsequence:
    def test_shorter_img_contained_in_ncbi(self):
        """IMG protein is shorter (truncated start) but fully contained in NCBI."""
        base = "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDN"  # 32 aa base
        ncbi_seq = base + "AR"  # 34 aa
        img_seq = base + "A"  # 33 aa — missing last 1 aa (33/34 = 97% overlap)
        img_seqs = {"AEZ55_0001": img_seq}
        ncbi_seqs = {"WP_111.1": ncbi_seq}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        result = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert result == {"AEZ55_0001": "EZ55_00001"}

    def test_ncbi_shorter_contained_in_img(self):
        """NCBI protein is shorter but fully contained in IMG (different stop codon)."""
        base = "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDN"  # 32 aa base
        img_seq = base + "AR"  # 34 aa
        ncbi_seq = base + "A"  # 33 aa (33/34 = 97% overlap)
        img_seqs = {"AEZ55_0001": img_seq}
        ncbi_seqs = {"WP_111.1": ncbi_seq}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        result = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert result == {"AEZ55_0001": "EZ55_00001"}

    def test_below_overlap_threshold(self):
        """Subsequence too short relative to the longer sequence → no match."""
        img_seq = "MKTLL"  # 5 aa
        ncbi_seq = "MKTLLAVGKRFYQNDAMKTLL"  # 21 aa — 5/21 = 24% overlap
        img_seqs = {"AEZ55_0001": img_seq}
        ncbi_seqs = {"WP_111.1": ncbi_seq}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        result = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert result == {}

    def test_skips_very_short_sequences(self):
        """Sequences shorter than 30 aa are skipped."""
        img_seqs = {"AEZ55_0001": "MKTLL"}  # 5 aa
        ncbi_seqs = {"WP_111.1": "MKTLLAVGK"}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        result = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert result == {}

    def test_not_a_subsequence(self):
        """Different sequences, no containment."""
        img_seqs = {"AEZ55_0001": "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDNDIFFERENT"}
        ncbi_seqs = {"WP_111.1": "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDNSOMETHING"}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        result = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert result == {}


# --- map_by_diamond ---


class TestFragmentDeduplication:
    """Test the fragment dedup logic directly (no Diamond needed)."""

    def test_longest_fragment_wins(self):
        """When multiple fragments hit the same locus_tag, keep longest."""
        from collections import defaultdict

        # Simulate Diamond hits: 3 fragments hitting same gene
        all_hits = {
            "AEZ55_0001": ("EZ55_00001", 200),  # longest
            "AEZ55_0002": ("EZ55_00001", 150),  # shorter
            "AEZ55_0003": ("EZ55_00001", 80),  # shortest
        }

        lt_to_fragments = defaultdict(list)
        for img_id, (lt, qlen) in all_hits.items():
            lt_to_fragments[lt].append((img_id, qlen))

        extra_mapping = {}
        discarded = []
        for lt, fragments in lt_to_fragments.items():
            fragments.sort(key=lambda x: x[1], reverse=True)
            winner_id, winner_len = fragments[0]
            extra_mapping[winner_id] = lt
            for loser_id, loser_len in fragments[1:]:
                discarded.append((loser_id, lt, loser_len))

        assert extra_mapping == {"AEZ55_0001": "EZ55_00001"}
        assert len(discarded) == 2
        assert discarded[0][0] == "AEZ55_0002"
        assert discarded[1][0] == "AEZ55_0003"

    def test_independent_genes_not_deduped(self):
        """Fragments hitting different locus_tags are kept independently."""
        from collections import defaultdict

        all_hits = {
            "AEZ55_0001": ("EZ55_00001", 200),
            "AEZ55_0002": ("EZ55_00002", 150),  # different gene
        }

        lt_to_fragments = defaultdict(list)
        for img_id, (lt, qlen) in all_hits.items():
            lt_to_fragments[lt].append((img_id, qlen))

        extra_mapping = {}
        discarded = []
        for lt, fragments in lt_to_fragments.items():
            fragments.sort(key=lambda x: x[1], reverse=True)
            winner_id, winner_len = fragments[0]
            extra_mapping[winner_id] = lt
            for loser_id, loser_len in fragments[1:]:
                discarded.append((loser_id, lt, loser_len))

        assert len(extra_mapping) == 2
        assert discarded == []


@pytest.mark.skipif(
    not shutil.which("diamond"), reason="diamond not installed"
)
class TestMapByDiamond:
    def test_near_identical_match(self):
        """Two proteins with a few substitutions should match via Diamond."""
        # Create a realistic-length protein (100 aa)
        base_seq = "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDNARRFYQKISPLWCHGTVEDNARR" * 2
        img_seq = base_seq  # identical
        ncbi_seq = base_seq[:5] + "X" + base_seq[6:]  # 1 substitution

        img_seqs = {"AEZ55_0001": img_seq}
        ncbi_seqs = {"WP_111.1": ncbi_seq}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        mapping, discarded = map_by_diamond(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert mapping == {"AEZ55_0001": "EZ55_00001"}
        assert discarded == []

    def test_no_match_below_identity(self):
        """Completely different proteins should not match."""
        # Two unrelated proteins of realistic length
        img_seqs = {"AEZ55_0001": "MKTLLAVGKRFYQNDAMEK" * 5}
        ncbi_seqs = {"WP_111.1": "WCHGTVEDNARRPSQDTVA" * 5}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        mapping, discarded = map_by_diamond(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert mapping == {}

    def test_empty_query_returns_empty(self):
        """When all query proteins are too short, returns empty without crashing."""
        img_seqs = {"AEZ55_0001": "MKTLL"}  # 5 aa — too short for query
        ncbi_seqs = {"WP_111.1": "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDN" * 3}
        protein_to_lt = {"WP_111.1": "EZ55_00001"}

        mapping, discarded = map_by_diamond(
            img_seqs, ncbi_seqs, protein_to_lt, {"AEZ55_0001"}, set()
        )
        assert mapping == {}
        assert discarded == []


# --- remap_img_headers ---


class TestRemapImgHeaders:
    def test_already_aez55_headers(self):
        """If headers already start with AEZ55_, return unchanged."""
        seqs = {"AEZ55_0001": "MKTLL", "AEZ55_0002": "AVGKM"}
        result = remap_img_headers(seqs, None)
        assert result == seqs

    def test_remap_via_gff(self, tmp_path):
        """Numeric OID headers remapped using GFF gene_id attribute."""
        gff = tmp_path / "test.gff"
        gff.write_text(
            "scf1\tgenbank_gmh_script\tgene\t1\t500\t.\t+\t.\tID=12345; gene_id AEZ55_0001\n"
            "scf1\tgenbank_gmh_script\tgene\t600\t1000\t.\t-\t.\tID=12346; gene_id AEZ55_0002\n"
        )
        seqs = {"12345": "MKTLL", "12346": "AVGKM"}
        result = remap_img_headers(seqs, gff)
        assert result == {"AEZ55_0001": "MKTLL", "AEZ55_0002": "AVGKM"}

    def test_no_gff_warns(self):
        """If no GFF and headers are not AEZ55_, returns unchanged with warning."""
        seqs = {"12345": "MKTLL"}
        result = remap_img_headers(seqs, None)
        assert result == seqs  # Returned unchanged


# --- parse_img_gff_gene_ids ---


class TestParseImgGffGeneIds:
    def test_extracts_gene_ids(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(
            "# GFF comment\n"
            "scf1\tgenbank_gmh_script\tgene\t1\t500\t.\t+\t.\tgene_id AEZ55_0001\n"
            "scf1\tgenbank_gmh_script\tgene\t600\t1000\t.\t-\t.\tgene_id AEZ55_0002\n"
        )
        result = parse_img_gff_gene_ids(gff)
        assert result == {"AEZ55_0001": "AEZ55_0001", "AEZ55_0002": "AEZ55_0002"}

    def test_skips_comments_and_short_lines(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text("# comment\nshort line\n")
        result = parse_img_gff_gene_ids(gff)
        assert result == {}


# --- Integration-style test (all phases without Diamond) ---


class TestEndToEndWithoutDiamond:
    def test_exact_plus_subsequence_pipeline(self):
        """Test phases 1+2 together on a small dataset."""
        # Sequences must be >= 30 aa for subsequence matching
        exact_seq = "MKTLLAVGKRFYQNDAMEKISPLWCHGTVEDNARR"  # 35 aa
        ncbi_longer = "AVGKMLLAARFYQNDAMEKISPLWCHGTVEDNARR"  # 35 aa
        img_subseq = "AVGKMLLAARFYQNDAMEKISPLWCHGTVEDNAR"  # 34 aa, subseq of ncbi_longer (34/35 = 97%)

        img_seqs = {
            "AEZ55_0001": exact_seq,  # exact match
            "AEZ55_0002": img_subseq,  # subsequence of NCBI
            "AEZ55_0003": "COMPLETELYDIFFERENTPROTEINSEQUENCEXXXXXXXXXXXXXXXX",  # no match
        }
        ncbi_seqs = {
            "WP_111.1": exact_seq,
            "WP_222.1": ncbi_longer,
        }
        protein_to_lt = {"WP_111.1": "EZ55_00001", "WP_222.1": "EZ55_00002"}

        # Phase 1
        mapping, unmatched, used_ncbi = map_by_exact_match(
            img_seqs, ncbi_seqs, protein_to_lt
        )
        assert mapping == {"AEZ55_0001": "EZ55_00001"}
        assert "AEZ55_0002" in unmatched
        assert "AEZ55_0003" in unmatched

        # Phase 2
        extra = map_by_subsequence(
            img_seqs, ncbi_seqs, protein_to_lt, unmatched, used_ncbi
        )
        mapping.update(extra)
        assert mapping.get("AEZ55_0002") == "EZ55_00002"
        assert "AEZ55_0003" not in mapping


# --- Output CSV format test ---


class TestOutputFormat:
    def test_csv_has_correct_columns(self, tmp_path):
        """Verify the output CSV format matches what paperconfig expects."""
        output = tmp_path / "id_translation.csv"
        mapping = {"AEZ55_0001": "EZ55_00001", "AEZ55_0002": "EZ55_00002"}

        with open(output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["aez55_id", "locus_tag"])
            for img_id in sorted(mapping.keys()):
                writer.writerow([img_id, mapping[img_id]])

        # Read back and verify
        with open(output) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert reader.fieldnames == ["aez55_id", "locus_tag"]
        assert len(rows) == 2
        assert rows[0]["aez55_id"] == "AEZ55_0001"
        assert rows[0]["locus_tag"] == "EZ55_00001"
