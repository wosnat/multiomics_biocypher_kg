"""Unit tests for multiomics_kg/download/build_gene_mapping.py.

Coverage
--------
- _get_cynaorak_ID
- _get_cyanorak_id_map_from_gbk
- _get_ec_numbers_from_gbff
- ncbi_merge_cds_and_gene_entries
- load_gff_from_ncbi_only
- load_gff_from_ncbi_and_cyanorak
- build_gene_mapping
- step5_gene_mapping (from download_genome_data)
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.build_gene_mapping import (
    _get_cynaorak_ID,
    _get_cyanorak_id_map_from_gbk,
    _get_ec_numbers_from_gbff,
    build_gene_mapping,
    load_gff_from_ncbi_and_cyanorak,
    load_gff_from_ncbi_only,
    ncbi_merge_cds_and_gene_entries,
)


# ─── module-level fixture data ────────────────────────────────────────────────

# Minimal Cyanorak GBK with two CDS entries (for SeqIO.read, single record)
CYANORAK_GBK_CONTENT = """\
LOCUS       PRO_MED4          30 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Prochlorococcus marinus str. MED4.
ACCESSION   CYKMED4
VERSION     CYKMED4.1
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="Prochlorococcus marinus str. MED4"
     CDS             1..10
                     /locus_tag="PMM0001"
                     /product="DNA polymerase III, beta subunit"
                     /note="cyanorak ORF Id:CK_Pro_MED4_00001"
                     /note="cyanorak cluster number:CK_00000364"
     CDS             11..20
                     /locus_tag="PMM0002"
                     /product="hypothetical protein"
                     /note="cyanorak ORF Id:CK_Pro_MED4_00002"
ORIGIN
        1 atgcatgcat atgcatgcat atgcatgcat
//
"""

# Minimal NCBI GBFF with three CDS entries (for SeqIO.parse, multi-record OK)
NCBI_GBFF_CONTENT = """\
LOCUS       NC_005072         30 bp    DNA     circular BCT 01-JAN-2024
DEFINITION  Prochlorococcus marinus str. MED4, complete genome.
ACCESSION   NC_005072
VERSION     NC_005072.1
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="Prochlorococcus marinus str. MED4"
                     /db_xref="taxon:59919"
     CDS             1..10
                     /locus_tag="TX50_RS00020"
                     /old_locus_tag="PMM0001"
                     /EC_number="2.7.7.7"
                     /product="DNA polymerase III subunit beta"
     CDS             11..20
                     /locus_tag="TX50_RS00025"
                     /old_locus_tag="PMM0002"
                     /product="hypothetical protein"
     CDS             21..30
                     /locus_tag="TX50_RS00030"
                     /EC_number="6.3.5.3"
                     /EC_number="3.6.1.4"
                     /product="phosphoribosylformylglycinamidine synthase"
ORIGIN
        1 atgcatgcat atgcatgcat atgcatgcat
//
"""


# ─── helpers ─────────────────────────────────────────────────────────────────


def _make_ncbi_gff_df(
    locus_tag_ncbi: str = "TX50_RS00020",
    old_locus_tag: str = "PMM0001",
    gene_name: str = "dnaN",
    protein_id: str = "WP_011131639.1",
    include_old_locus_tag: bool = True,
    start: int = 174,
    end: int = 1331,
    strand: str = "+",
) -> pd.DataFrame:
    """Return a minimal DataFrame mimicking gffpandas output for one NCBI gene+CDS pair.

    Both gene and CDS rows share the same set of columns (NaN where unused),
    which is exactly what gffpandas produces after attributes_to_columns().

    When ``include_old_locus_tag=False``, the ``old_locus_tag`` column is
    completely absent (as in Alteromonas genomes whose GFF has no old_locus_tag
    attribute, so gffpandas never creates that column).
    """
    common: dict = dict(
        Name=gene_name,
        gene=gene_name,
        locus_tag=locus_tag_ncbi,
        source="RefSeq",
        start=start,
        end=end,
        strand=strand,
        Note=None,
        exception=None,
        inference="COORDINATES: similar to AA sequence:RefSeq:WP_002806737.1",
        product="DNA polymerase III subunit beta",
        protein_id=protein_id,
        Ontology_term="GO:0006260,GO:0003677",
        go_component="DNA polymerase III complex|0009360||IEA",
        go_function="DNA binding|0003677||IEA",
        go_process="DNA replication|0006260||IEA",
        gene_synonym=None,
        Parent=None,
    )
    if include_old_locus_tag:
        common["old_locus_tag"] = old_locus_tag
    gene_row = {"type": "gene", "ID": f"gene-{locus_tag_ncbi}", **common}
    cds_row = {
        "type": "CDS",
        "ID": f"cds-{protein_id}",
        **dict(common, Parent=f"gene-{locus_tag_ncbi}"),
    }
    return pd.DataFrame([gene_row, cds_row])


def _make_cyan_gff_df(
    cyanorak_id: str = "CK_Pro_MED4_00001",
    locus_tag: str = "PMM0001",
    start: int = 174,
    end: int = 1331,
    strand: str = "+",
) -> pd.DataFrame:
    """Return a minimal Cyanorak GFF DataFrame with one CDS row."""
    return pd.DataFrame([{
        "type": "CDS",
        "ID": cyanorak_id,
        "Name": "dnaN",
        "start": start,
        "end": end,
        "strand": strand,
        "locus_tag": locus_tag,
        "product": "DNA polymerase III, beta subunit",
        "cluster_number": "CK_00000364",
        "Ontology_term": "GO:0006260",
    }])


# ─── _get_cynaorak_ID ─────────────────────────────────────────────────────────

class TestGetCynaorakID:
    def _mock_rec(self, notes: list[str], rec_id: str = "GENE1") -> MagicMock:
        rec = MagicMock()
        rec.qualifiers = {"note": notes}
        rec.id = rec_id
        return rec

    def test_extracts_cyanorak_id(self):
        rec = self._mock_rec(
            ["cyanorak ORF Id:CK_Pro_MED4_00001", "cyanorak cluster number:CK_00000364"]
        )
        assert _get_cynaorak_ID(rec) == "CK_Pro_MED4_00001"

    def test_strips_whitespace_from_id(self):
        rec = self._mock_rec(["cyanorak ORF Id: CK_Pro_MED4_00001"])
        assert _get_cynaorak_ID(rec) == "CK_Pro_MED4_00001"

    def test_warns_on_multiple_ids(self, capsys):
        rec = self._mock_rec(
            ["cyanorak ORF Id:CK_Pro_MED4_00001", "cyanorak ORF Id:CK_Pro_MED4_99999"]
        )
        result = _get_cynaorak_ID(rec)
        assert result == "CK_Pro_MED4_00001"
        assert "Warning" in capsys.readouterr().out

    def test_returns_first_id_when_multiple(self):
        rec = self._mock_rec(
            ["cyanorak ORF Id:FIRST", "cyanorak ORF Id:SECOND"]
        )
        assert _get_cynaorak_ID(rec) == "FIRST"


# ─── _get_cyanorak_id_map_from_gbk ───────────────────────────────────────────

class TestGetCyanorakIdMapFromGbk:
    def test_maps_cyanorak_id_to_locus_tag(self, tmp_path):
        gbk = tmp_path / "Pro_MED4.gbk"
        gbk.write_text(CYANORAK_GBK_CONTENT)
        result = _get_cyanorak_id_map_from_gbk(str(gbk))
        assert result["CK_Pro_MED4_00001"] == "PMM0001"
        assert result["CK_Pro_MED4_00002"] == "PMM0002"

    def test_returns_only_cds_entries(self, tmp_path):
        gbk = tmp_path / "Pro_MED4.gbk"
        gbk.write_text(CYANORAK_GBK_CONTENT)
        result = _get_cyanorak_id_map_from_gbk(str(gbk))
        # Only CDS features contribute; source feature is excluded
        assert len(result) == 2


# ─── _get_ec_numbers_from_gbff ───────────────────────────────────────────────

class TestGetEcNumbersFromGbff:
    def test_extracts_single_ec_number(self, tmp_path):
        gbff = tmp_path / "genomic.gbff"
        gbff.write_text(NCBI_GBFF_CONTENT)
        result = _get_ec_numbers_from_gbff(str(gbff))
        assert result["TX50_RS00020"] == ["2.7.7.7"]

    def test_cds_without_ec_not_included(self, tmp_path):
        gbff = tmp_path / "genomic.gbff"
        gbff.write_text(NCBI_GBFF_CONTENT)
        result = _get_ec_numbers_from_gbff(str(gbff))
        # TX50_RS00025 has no EC_number
        assert "TX50_RS00025" not in result

    def test_accumulates_multiple_ec_numbers(self, tmp_path):
        gbff = tmp_path / "genomic.gbff"
        gbff.write_text(NCBI_GBFF_CONTENT)
        result = _get_ec_numbers_from_gbff(str(gbff))
        # TX50_RS00030 has two EC_number qualifiers
        assert set(result["TX50_RS00030"]) == {"6.3.5.3", "3.6.1.4"}

    def test_returns_empty_dict_for_no_cds(self, tmp_path):
        content = """\
LOCUS       EMPTY            10 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Empty.
ACCESSION   EMPTY
VERSION     EMPTY.1
FEATURES             Location/Qualifiers
     source          1..10
                     /organism="Test"
ORIGIN
        1 atgcatgcat
//
"""
        gbff = tmp_path / "genomic.gbff"
        gbff.write_text(content)
        result = _get_ec_numbers_from_gbff(str(gbff))
        assert result == {}


# ─── ncbi_merge_cds_and_gene_entries ─────────────────────────────────────────

class TestNcbiMergeCdsAndGeneEntries:
    def test_produces_one_row_per_gene(self):
        df = _make_ncbi_gff_df()
        result = ncbi_merge_cds_and_gene_entries(df)
        assert len(result) == 1

    def test_renames_locus_tag_cds_to_locus_tag_ncbi(self):
        df = _make_ncbi_gff_df()
        result = ncbi_merge_cds_and_gene_entries(df)
        assert "locus_tag_ncbi" in result.columns
        assert result["locus_tag_ncbi"].iloc[0] == "TX50_RS00020"

    def test_renames_old_locus_tag_gene_to_locus_tag(self):
        df = _make_ncbi_gff_df()
        result = ncbi_merge_cds_and_gene_entries(df)
        # old_locus_tag_gene → locus_tag
        assert "locus_tag" in result.columns
        assert result["locus_tag"].iloc[0] == "PMM0001"

    def test_url_encoded_old_locus_tag_is_decoded(self):
        # PMT0003%2CPMT_0003 (URL-encoded comma) → "PMT0003,PMT_0003"
        df = _make_ncbi_gff_df(old_locus_tag="PMT0003%2CPMT_0003")
        result = ncbi_merge_cds_and_gene_entries(df)
        assert result["locus_tag"].iloc[0] == "PMT0003,PMT_0003"

    def test_old_locus_tags_column_stores_decoded_value(self):
        df = _make_ncbi_gff_df(old_locus_tag="PMT0003%2CPMT_0003")
        result = ncbi_merge_cds_and_gene_entries(df)
        assert "old_locus_tags" in result.columns
        assert result["old_locus_tags"].iloc[0] == "PMT0003,PMT_0003"

    def test_product_comes_from_cds_row(self):
        df = _make_ncbi_gff_df()
        result = ncbi_merge_cds_and_gene_entries(df)
        assert result["product"].iloc[0] == "DNA polymerase III subunit beta"

    def test_missing_old_locus_tag_column_handled_gracefully(self):
        # Alteromonas genomes lack old_locus_tag in their GFF
        df = _make_ncbi_gff_df(include_old_locus_tag=False)
        result = ncbi_merge_cds_and_gene_entries(df)
        # locus_tag column is absent (old_locus_tag_gene not in columns)
        assert "locus_tag" not in result.columns
        # locus_tag_ncbi is present
        assert "locus_tag_ncbi" in result.columns

    def test_two_genes_produce_two_rows(self):
        df1 = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS00020", old_locus_tag="PMM0001",
            gene_name="dnaN", protein_id="WP_011131639.1",
        )
        df2 = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS00025", old_locus_tag="PMM0002",
            gene_name="", protein_id="WP_011131640.1",
        )
        combined = pd.concat([df1, df2], ignore_index=True)
        result = ncbi_merge_cds_and_gene_entries(combined)
        assert len(result) == 2
        assert set(result["locus_tag_ncbi"]) == {"TX50_RS00020", "TX50_RS00025"}


# ─── load_gff_from_ncbi_only ─────────────────────────────────────────────────

_LOAD_GFF = "multiomics_kg.download.build_gene_mapping.load_gff"


class TestLoadGffFromNcbiOnly:
    def test_ncbi_only_uses_old_locus_tag_as_locus_tag(self):
        ncbi_df = _make_ncbi_gff_df()
        with patch(_LOAD_GFF, return_value=ncbi_df):
            result = load_gff_from_ncbi_only("ncbi.gff")
        assert result["locus_tag"].iloc[0] == "PMM0001"

    def test_falls_back_to_locus_tag_ncbi_when_old_locus_tag_absent(self):
        ncbi_df = _make_ncbi_gff_df(include_old_locus_tag=False)
        with patch(_LOAD_GFF, return_value=ncbi_df):
            result = load_gff_from_ncbi_only("ncbi.gff")
        assert result["locus_tag"].iloc[0] == "TX50_RS00020"

    def test_drops_rows_with_no_locus_tag(self):
        # Create two genes, one with old_locus_tag and one without (no fallback)
        df1 = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        df2 = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS00025", old_locus_tag=None,
            gene_name="", protein_id="WP_011131640.1",
        )
        # Remove locus_tag_ncbi from df2 to force NaN locus_tag
        df2["locus_tag_ncbi"] = None
        combined = pd.concat([df1, df2], ignore_index=True)
        with patch(_LOAD_GFF, return_value=combined):
            result = load_gff_from_ncbi_only("ncbi.gff")
        # Row with no resolvable locus_tag must be dropped
        assert all(result["locus_tag"].notna())


# ─── load_gff_from_ncbi_and_cyanorak ─────────────────────────────────────────

_GET_CYAN_MAP = "multiomics_kg.download.build_gene_mapping._get_cyanorak_id_map_from_gbk"


class TestLoadGffFromNcbiAndCyanorak:
    def _run(self, ncbi_df, cyan_df, cyan_id_map):
        with (
            patch(_LOAD_GFF, side_effect=[cyan_df, ncbi_df]),
            patch(_GET_CYAN_MAP, return_value=cyan_id_map),
        ):
            return load_gff_from_ncbi_and_cyanorak("ncbi.gff", "cyan.gff", "cyan.gbk")

    def test_ncbi_gene_matched_to_cyanorak_carries_data_from_both(self):
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        cyan_df = _make_cyan_gff_df(cyanorak_id="CK_Pro_MED4_00001", locus_tag="PMM0001")
        result = self._run(ncbi_df, cyan_df, {"CK_Pro_MED4_00001": "PMM0001"})
        # Should have a row for PMM0001 with data from both sources
        assert len(result) == 1
        row = result.iloc[0]
        # NCBI data: locus_tag_ncbi is present
        assert row["locus_tag_ncbi"] == "TX50_RS00020"
        # Cyanorak data: cluster_number column is present
        assert "cluster_number" in result.columns
        assert row["cluster_number"] == "CK_00000364"

    def test_ncbi_gene_without_cyanorak_match_appears_in_output(self):
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        # Cyanorak has a different gene, so NCBI gene won't match
        cyan_df = _make_cyan_gff_df(cyanorak_id="CK_Pro_MED4_00099", locus_tag="PMM0099")
        result = self._run(ncbi_df, cyan_df, {"CK_Pro_MED4_00099": "PMM0099"})
        # PMM0001 (NCBI-only) must be in output
        assert any(result["locus_tag"] == "PMM0001")

    def test_cyanorak_only_gene_not_dropped(self):
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        # Cyanorak has its own extra gene
        cyan_extra = pd.DataFrame([{
            "type": "CDS", "ID": "CK_Pro_MED4_00099",
            "Name": "hypothetical", "start": 5000, "end": 5300, "strand": "+",
            "locus_tag": "PMM0099", "product": "hypothetical protein",
            "cluster_number": None, "Ontology_term": None,
        }])
        cyan_combined = pd.concat(
            [_make_cyan_gff_df("CK_Pro_MED4_00001", "PMM0001"), cyan_extra],
            ignore_index=True,
        )
        result = self._run(
            ncbi_df, cyan_combined,
            {"CK_Pro_MED4_00001": "PMM0001", "CK_Pro_MED4_00099": "PMM0099"},
        )
        # Both genes should appear
        locus_tags = set(result["locus_tag"].dropna())
        assert "PMM0001" in locus_tags
        assert "PMM0099" in locus_tags

    def test_ncbi_gene_with_two_old_locus_tags_produces_one_output_row(self):
        # URL-encoded comma: "PMM0001%2CPMT0001" → "PMM0001,PMT0001"
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001%2CPMT0001")
        cyan_df = _make_cyan_gff_df(cyanorak_id="CK_Pro_MED4_00001", locus_tag="PMM0001")
        result = self._run(ncbi_df, cyan_df, {"CK_Pro_MED4_00001": "PMM0001"})
        # Should be deduplicated to one row per locus_tag_ncbi
        assert result["locus_tag_ncbi"].nunique() == 1
        assert len(result[result["locus_tag_ncbi"] == "TX50_RS00020"]) == 1

    def test_final_rename_applied(self):
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        cyan_df = _make_cyan_gff_df()
        result = self._run(ncbi_df, cyan_df, {"CK_Pro_MED4_00001": "PMM0001"})
        # _FINAL_MERGE_RENAME maps Name_ncbi → gene_names
        assert "gene_names" in result.columns
        # start_ncbi → start
        assert "start" in result.columns

    # --- Position fallback tests ---

    def test_position_fallback_merges_unmatched_pair(self):
        """MIT9313-style: NCBI has PMT_0107, Cyanorak has PMT0107 at same coords."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="AKG35_RS00545",
            old_locus_tag="PMT_0107",
            start=122194, end=123066, strand="-",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MIT9313_00107",
            locus_tag="PMT0107",
            start=122194, end=123066, strand="-",
        )
        result = self._run(
            ncbi_df, cyan_df,
            {"CK_Pro_MIT9313_00107": "PMT0107"},
        )
        # Should merge into one row with data from both sources
        assert len(result) == 1
        row = result.iloc[0]
        assert row["locus_tag_ncbi"] == "AKG35_RS00545"
        assert row["locus_tag_cyanoak"] == "CK_Pro_MIT9313_00107"
        assert row["cluster_number"] == "CK_00000364"
        # Position merge note should be set
        assert "position_merge_note" in result.columns
        assert pd.notna(row["position_merge_note"])
        assert "position_merge" in str(row["position_merge_note"])

    def test_position_fallback_skips_different_strand(self):
        """Same coords but different strand → no merge."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=122194, end=123066, strand="-",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=122194, end=123066, strand="+",  # different strand
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        # Should remain as two separate rows (no fallback merge)
        assert len(result) == 2

    def test_position_fallback_skips_low_overlap(self):
        """Overlap < 90% → no merge."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=800, end=2000, strand="+",  # only 200bp overlap / 1200bp max = 16%
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_skips_large_coord_diff(self):
        """End diff > 3bp → no merge even with high overlap."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=100, end=1010, strand="+",  # 10bp end diff, ~99% overlap
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_skips_conflict(self):
        """Two Cyanorak entries at same position as one NCBI entry → skip both."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_2613", start=1128809, end=1129006, strand="-",
        )
        cyan1 = _make_cyan_gff_df(
            cyanorak_id="CK_02281", locus_tag="PMT2281",
            start=1128809, end=1129006, strand="-",
        )
        cyan2 = _make_cyan_gff_df(
            cyanorak_id="CK_02283", locus_tag="PMT2283",
            start=1128809, end=1129006, strand="-",
        )
        cyan_combined = pd.concat([cyan1, cyan2], ignore_index=True)
        result = self._run(
            ncbi_df, cyan_combined,
            {"CK_02281": "PMT2281", "CK_02283": "PMT2283"},
        )
        # All three should remain separate (conflict = no merge)
        assert len(result) == 3
        # NCBI entry should NOT have Cyanorak data
        ncbi_row = result[result["locus_tag_ncbi"] == "TX50_RS00020"].iloc[0]
        assert pd.isna(ncbi_row["locus_tag_cyanoak"])


# ─── build_gene_mapping ───────────────────────────────────────────────────────

_NCBI_ONLY = "multiomics_kg.download.build_gene_mapping.load_gff_from_ncbi_only"
_NCBI_CYAN = "multiomics_kg.download.build_gene_mapping.load_gff_from_ncbi_and_cyanorak"
_GET_EC    = "multiomics_kg.download.build_gene_mapping._get_ec_numbers_from_gbff"


class TestBuildGeneMapping:
    def _minimal_df(self):
        return pd.DataFrame([{
            "locus_tag": "PMM0001",
            "locus_tag_ncbi": "TX50_RS00020",
            "gene_names": "dnaN",
            "product": "DNA polymerase III subunit beta",
        }])

    def test_ncbi_only_when_no_cyanorak_files(self):
        df = self._minimal_df()
        with patch(_NCBI_ONLY, return_value=df) as mock_ncbi, \
             patch(_NCBI_CYAN) as mock_cyan:
            result = build_gene_mapping(ncbi_gff_file="ncbi.gff")
        mock_ncbi.assert_called_once_with("ncbi.gff")
        mock_cyan.assert_not_called()

    def test_ncbi_and_cyanorak_when_both_files_provided(self, tmp_path):
        cyan_gff = tmp_path / "Pro_MED4.gff"
        cyan_gbk = tmp_path / "Pro_MED4.gbk"
        cyan_gff.write_text("##gff-version 3\n")
        cyan_gbk.write_text(CYANORAK_GBK_CONTENT)
        df = self._minimal_df()
        with patch(_NCBI_CYAN, return_value=df) as mock_cyan, \
             patch(_NCBI_ONLY) as mock_ncbi:
            result = build_gene_mapping(
                ncbi_gff_file="ncbi.gff",
                cyan_gff_file=str(cyan_gff),
                cyan_gbk_file=str(cyan_gbk),
            )
        mock_cyan.assert_called_once_with("ncbi.gff", str(cyan_gff), str(cyan_gbk))
        mock_ncbi.assert_not_called()

    def test_ec_numbers_merged_when_gbff_exists(self, tmp_path):
        gbff = tmp_path / "genomic.gbff"
        gbff.write_text(NCBI_GBFF_CONTENT)
        df = pd.DataFrame([{
            "locus_tag": "PMM0001",
            "locus_tag_ncbi": "TX50_RS00020",
            "gene_names": "dnaN",
        }])
        with patch(_NCBI_ONLY, return_value=df):
            result = build_gene_mapping(
                ncbi_gff_file="ncbi.gff",
                ncbi_gbff_file=str(gbff),
            )
        assert "ec_numbers" in result.columns
        assert result.loc[result["locus_tag_ncbi"] == "TX50_RS00020", "ec_numbers"].iloc[0] == "2.7.7.7"

    def test_no_ec_column_when_gbff_absent(self):
        df = self._minimal_df()
        with patch(_NCBI_ONLY, return_value=df):
            result = build_gene_mapping(
                ncbi_gff_file="ncbi.gff",
                ncbi_gbff_file=None,
            )
        assert "ec_numbers" not in result.columns

    def test_no_ec_column_when_gbff_path_not_on_disk(self, tmp_path):
        df = self._minimal_df()
        with patch(_NCBI_ONLY, return_value=df):
            result = build_gene_mapping(
                ncbi_gff_file="ncbi.gff",
                ncbi_gbff_file=str(tmp_path / "nonexistent.gbff"),
            )
        assert "ec_numbers" not in result.columns


# ─── step5_gene_mapping ───────────────────────────────────────────────────────

_BUILD_GENE_MAPPING = "multiomics_kg.download.build_gene_mapping.build_gene_mapping"
_STEP5_PROJECT_ROOT = "multiomics_kg.download.download_genome_data.PROJECT_ROOT"


class TestStep5GeneMapping:
    """Tests for step5_gene_mapping from download_genome_data."""

    @pytest.fixture(autouse=True)
    def _import_step5(self):
        from multiomics_kg.download.download_genome_data import step5_gene_mapping
        self.step5_gene_mapping = step5_gene_mapping

    def _genome(self, tmp_path, strain="MED4", cyan_org="Pro_MED4"):
        """Helper: genome dict using absolute data_dir so PROJECT_ROOT is bypassed."""
        d = tmp_path / "genomes" / strain
        d.mkdir(parents=True)
        return {"strain_name": strain, "data_dir": str(d), "cyanorak_organism": cyan_org}

    def test_skips_when_gene_mapping_exists_and_no_force(self, tmp_path, capsys):
        genome = self._genome(tmp_path)
        (tmp_path / "genomes" / "MED4" / "gene_mapping.csv").write_text("locus_tag\nPMM0001\n")
        with patch(_BUILD_GENE_MAPPING) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        mock_build.assert_not_called()

    def test_skip_status_is_reported_for_existing_csv(self, tmp_path, capsys):
        genome = self._genome(tmp_path)
        (tmp_path / "genomes" / "MED4" / "gene_mapping.csv").write_text("locus_tag\nPMM0001\n")
        self.step5_gene_mapping([genome], force=False)
        assert "SKIP_EXISTS" in capsys.readouterr().out

    def test_skips_when_ncbi_gff_missing(self, tmp_path, capsys):
        genome = self._genome(tmp_path)
        # No genomic.gff → should report SKIP_NO_GFF without calling build_gene_mapping
        with patch(_BUILD_GENE_MAPPING) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        mock_build.assert_not_called()
        assert "SKIP_NO_GFF" in capsys.readouterr().out

    def test_ncbi_gbff_passed_as_none_when_absent(self, tmp_path):
        genome = self._genome(tmp_path)
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")
        # No genomic.gbff → ncbi_gbff should be None

        mock_df = pd.DataFrame([{"locus_tag": "PMM0001"}])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        _, kwargs = mock_build.call_args
        assert kwargs["ncbi_gbff_file"] is None

    def test_cyanorak_files_passed_when_present(self, tmp_path):
        genome = self._genome(tmp_path, cyan_org="Pro_MED4")
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")
        cyan_dir = d / "cyanorak"
        cyan_dir.mkdir()
        (cyan_dir / "Pro_MED4.gff").write_text("##gff-version 3\n")
        (cyan_dir / "Pro_MED4.gbk").write_text(CYANORAK_GBK_CONTENT)

        mock_df = pd.DataFrame([{"locus_tag": "PMM0001"}])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        _, kwargs = mock_build.call_args
        assert kwargs["cyan_gff_file"] is not None
        assert kwargs["cyan_gbk_file"] is not None

    def test_cyanorak_files_none_when_not_on_disk(self, tmp_path):
        genome = self._genome(tmp_path, cyan_org="Pro_MED4")
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")
        # cyanorak/ dir exists but GFF/GBK files are absent

        mock_df = pd.DataFrame([{"locus_tag": "PMM0001"}])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        _, kwargs = mock_build.call_args
        assert kwargs["cyan_gff_file"] is None
        assert kwargs["cyan_gbk_file"] is None

    def test_cyanorak_files_none_when_no_cyanorak_organism(self, tmp_path):
        genome = self._genome(tmp_path, cyan_org="")  # no Cyanorak organism
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")

        mock_df = pd.DataFrame([{"locus_tag": "PMM0001"}])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df) as mock_build:
            self.step5_gene_mapping([genome], force=False)
        _, kwargs = mock_build.call_args
        assert kwargs["cyan_gff_file"] is None
        assert kwargs["cyan_gbk_file"] is None

    def test_exception_from_build_gene_mapping_reports_failed(self, tmp_path, capsys):
        genome = self._genome(tmp_path)
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")

        with patch(_BUILD_GENE_MAPPING, side_effect=RuntimeError("GFF parse error")):
            self.step5_gene_mapping([genome], force=False)
        assert "FAILED" in capsys.readouterr().out

    def test_force_overwrites_existing_gene_mapping(self, tmp_path):
        genome = self._genome(tmp_path)
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")
        old_csv = d / "gene_mapping.csv"
        old_csv.write_text("locus_tag\nPMM9999\n")

        mock_df = pd.DataFrame([{"locus_tag": "PMM0001"}])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df):
            self.step5_gene_mapping([genome], force=True)
        # File should be overwritten with new content
        assert "PMM0001" in old_csv.read_text()

    def test_successful_run_reports_ok_with_gene_count(self, tmp_path, capsys):
        genome = self._genome(tmp_path)
        d = tmp_path / "genomes" / "MED4"
        (d / "genomic.gff").write_text("##gff-version 3\n")

        mock_df = pd.DataFrame([
            {"locus_tag": "PMM0001"},
            {"locus_tag": "PMM0002"},
        ])
        with patch(_BUILD_GENE_MAPPING, return_value=mock_df):
            self.step5_gene_mapping([genome], force=False)
        out = capsys.readouterr().out
        assert "OK" in out
        assert "2" in out
