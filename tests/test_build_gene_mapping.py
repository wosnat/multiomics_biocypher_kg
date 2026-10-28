"""Unit tests for multiomics_kg/download/build_gene_mapping.py.

Coverage
--------
- _get_cyanorak_id
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
    _get_cyanorak_id,
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
        seq_id="NC_005072.1",
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


# ─── _get_cyanorak_id ─────────────────────────────────────────────────────────

class TestGetCyanorakId:
    def _mock_rec(self, qualifiers: dict | None = None) -> MagicMock:
        rec = MagicMock()
        rec.qualifiers = qualifiers if qualifiers is not None else {}
        return rec

    def test_extracts_cyanorak_id(self):
        rec = self._mock_rec(
            {"note": ["cyanorak ORF Id:CK_Pro_MED4_00001", "cyanorak cluster number:CK_00000364"],
             "locus_tag": ["PMM0001"]}
        )
        assert _get_cyanorak_id(rec) == "CK_Pro_MED4_00001"

    def test_strips_whitespace_from_id(self):
        rec = self._mock_rec({"note": ["cyanorak ORF Id: CK_Pro_MED4_00001"], "locus_tag": ["PMM0001"]})
        assert _get_cyanorak_id(rec) == "CK_Pro_MED4_00001"

    def test_returns_none_when_no_note_qualifier(self):
        rec = self._mock_rec({"locus_tag": ["PMM0001"]})
        assert _get_cyanorak_id(rec) is None

    def test_returns_none_when_no_matching_note(self):
        rec = self._mock_rec({"note": ["some other note"], "locus_tag": ["PMM0001"]})
        assert _get_cyanorak_id(rec) is None

    def test_warns_on_multiple_ids(self):
        rec = self._mock_rec(
            {"note": ["cyanorak ORF Id:CK_Pro_MED4_00001", "cyanorak ORF Id:CK_Pro_MED4_99999"],
             "locus_tag": ["PMM0001"]}
        )
        result = _get_cyanorak_id(rec)
        assert result == "CK_Pro_MED4_00001"

    def test_returns_first_id_when_multiple(self):
        rec = self._mock_rec(
            {"note": ["cyanorak ORF Id:FIRST", "cyanorak ORF Id:SECOND"], "locus_tag": ["PMM0001"]}
        )
        assert _get_cyanorak_id(rec) == "FIRST"


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

    def test_locus_tag_only_on_gene_row_old_gff_format(self):
        # Older NCBI annotwriter GFFs (HNLC1, CB0101) put locus_tag on the gene
        # row only; the CDS row has none. locus_tag_ncbi must still be set.
        df = _make_ncbi_gff_df()
        df.loc[df["type"] == "CDS", "locus_tag"] = None
        result = ncbi_merge_cds_and_gene_entries(df)
        assert result["locus_tag_ncbi"].iloc[0] == "TX50_RS00020"

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

    def test_seqid_preserved_as_first_column(self):
        df = _make_ncbi_gff_df()
        result = ncbi_merge_cds_and_gene_entries(df)
        assert "seqid" in result.columns
        assert result.columns[0] == "seqid"
        assert result["seqid"].iloc[0] == "NC_005072.1"

    def test_seqid_absent_when_gff_lacks_seq_id(self):
        # Some minimal/test GFF DataFrames may not have seq_id — graceful skip
        df = _make_ncbi_gff_df()
        df = df.drop(columns=["seq_id"])
        result = ncbi_merge_cds_and_gene_entries(df)
        assert "seqid" not in result.columns


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

    def _run_anchored(self, ncbi_df, cyan_df, cyan_id_map):
        """Like _run, plus one locus_tag-paired anchor gene at offset 0 far away.

        The position fallback measures the NCBI→Cyanorak offset on paired genes
        and never assumes one (assuming 0 paired genes by coincidence when the
        frames differ — KORDI-49/52/100, 2026-10-26), so a fixture testing it
        needs a paired gene to measure from. The anchor is dropped from the
        result.
        """
        seqid = ncbi_df["seq_id"].dropna().iloc[0]
        anchor_n = _make_ncbi_gff_df(locus_tag_ncbi="RS_ANCHOR", old_locus_tag="ANCHOR",
                                     protein_id="WP_anchor", start=900_001, end=900_300)
        anchor_n["seq_id"] = seqid
        anchor_c = _make_cyan_gff_df(cyanorak_id="CK_ANCHOR", locus_tag="ANCHOR",
                                     start=900_001, end=900_300)
        result = self._run(pd.concat([ncbi_df, anchor_n], ignore_index=True),
                           pd.concat([cyan_df, anchor_c], ignore_index=True),
                           {**cyan_id_map, "CK_ANCHOR": "ANCHOR"})
        return result[result["locus_tag"] != "ANCHOR"].reset_index(drop=True)

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
        result = self._run_anchored(
            ncbi_df, cyan_df,
            {"CK_Pro_MIT9313_00107": "PMT0107"},
        )
        # Should merge into one row with data from both sources
        assert len(result) == 1
        row = result.iloc[0]
        assert row["locus_tag_ncbi"] == "AKG35_RS00545"
        assert row["locus_tag_cyanorak"] == "CK_Pro_MIT9313_00107"
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

    def test_position_fallback_skips_different_stop_codon_far(self):
        """Stop codons 1000 bp apart (and out of frame) → no merge."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=800, end=2000, strand="+",  # 3' end 1000 bp past N's, not in frame
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_skips_different_stop_codon_near(self):
        """Stop codons 10 bp apart (not a multiple of 3) → no merge."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=100, end=1010, strand="+",  # 3' end 10 bp past N's, out of frame
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_merges_start_codon_recall_plus_strand(self):
        """MIT9313 PMT0040/PMT_0040: same stop codon, NCBI re-called the start
        93 bp downstream (reciprocal overlap 0.76). Same CDS -> must merge."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="AKG35_RS00210", old_locus_tag="PMT_0040",
            start=41255, end=41551, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MIT9313_00040", locus_tag="PMT0040",
            start=41162, end=41551, strand="+",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_Pro_MIT9313_00040": "PMT0040"})
        assert len(result) == 1
        row = result.iloc[0]
        assert row["locus_tag_cyanorak"] == "CK_Pro_MIT9313_00040"
        assert "PMT0040" in str(row["old_locus_tags"])
        assert "position_merge:PMT0040" in str(row["position_merge_note"])

    def test_position_fallback_merges_start_codon_recall_minus_strand(self):
        """MIT9313 PMT0110/PMT_0110 (- strand): the stop codon is the genomic
        *start*; the start-codon re-call moves the genomic *end* by 24 bp."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="AKG35_RS00560", old_locus_tag="PMT_0110",
            start=125340, end=126095, strand="-",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MIT9313_00110", locus_tag="PMT0110",
            start=125340, end=126119, strand="-",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_Pro_MIT9313_00110": "PMT0110"})
        assert len(result) == 1
        assert result.iloc[0]["locus_tag_cyanorak"] == "CK_Pro_MIT9313_00110"

    def test_position_fallback_merges_cyanorak_stop_codon_excluded_plus(self):
        """MED4 PMM1858a/TX50_RS09595: Cyanorak coordinates exclude the stop
        codon (3' end 3 bp short, in frame). Same ORF -> must merge."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS09595", old_locus_tag="PMM1858",
            start=627532, end=627696, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MED4_01858a", locus_tag="PMM1858a",
            start=627532, end=627693, strand="+",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_Pro_MED4_01858a": "PMM1858a"})
        assert len(result) == 1
        assert result.iloc[0]["locus_tag_cyanorak"] == "CK_Pro_MED4_01858a"

    def test_position_fallback_merges_cyanorak_stop_codon_excluded_minus(self):
        """MED4 PMM1719/TX50_RS09730 (- strand): stop codon at the genomic
        start; Cyanorak's start is 3 bp inside NCBI's, in frame."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS09730", old_locus_tag="PMM1719x",
            start=965317, end=965487, strand="-",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MED4_01719", locus_tag="PMM1719",
            start=965320, end=965487, strand="-",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_Pro_MED4_01719": "PMM1719"})
        assert len(result) == 1

    def test_position_fallback_skips_out_of_frame_3prime_end(self):
        """3' ends differ by 4 bp: different reading frame -> different gene."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=100, end=996, strand="+",
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_skips_shared_start_different_stop(self):
        """Same start coordinate but a different stop codon is a different
        reading frame / different gene -> no merge, however large the overlap."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=100, end=1003, strand="+",
        )
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_merges_nested_in_frame_orf(self):
        """A Cyanorak ORF ending in frame well inside the NCBI span (600 bp
        short of N's stop) is the same reading frame with no stop codon in
        between → merged. Pins the deliberate absence of a distance bound."""
        ncbi_df = _make_ncbi_gff_df(
            old_locus_tag="PMT_0107", start=100, end=1000, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_00107", locus_tag="PMT0107",
            start=100, end=400, strand="+",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 1
        assert result.iloc[0]["position_merge_note"] == "position_merge:PMT0107→PMT_0107"

    def test_position_fallback_note_names_ncbi_locus_tag_without_old_locus_tag(self):
        """NCBI gene with no old_locus_tag: the note must name locus_tag_ncbi,
        not the still-unfilled NaN locus_tag."""
        ncbi_df = _make_ncbi_gff_df(
            locus_tag_ncbi="TX50_RS09915", include_old_locus_tag=False,
            start=595314, end=595439, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_Pro_MED4_00624", locus_tag="PMM0624",
            start=595314, end=595439, strand="+",
        )
        result = self._run_anchored(ncbi_df, cyan_df, {"CK_Pro_MED4_00624": "PMM0624"})
        assert len(result) == 1
        assert result.iloc[0]["position_merge_note"] == "position_merge:PMM0624→TX50_RS09915"
        assert result.iloc[0]["locus_tag"] == "TX50_RS09915"

    def test_position_fallback_never_assumes_an_offset(self):
        """No locus_tag pair → the frame offset is unknown → no position merge,
        even when coordinates coincide. Assuming offset 0 here paired ~330-400
        unrelated genes per strain on KORDI-49/52/100 and CC9616 (2026-10-26)."""
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMT_0107", start=100, end=1000, strand="+")
        cyan_df = _make_cyan_gff_df(cyanorak_id="CK_00107", locus_tag="PMT0107",
                                    start=100, end=1000, strand="+")
        result = self._run(ncbi_df, cyan_df, {"CK_00107": "PMT0107"})
        assert len(result) == 2

    def test_position_fallback_skips_reverse_conflict(self):
        """One Cyanorak call in frame inside two nested NCBI genes → neither merged.
        Outer 1–4587 and inner 4432–4584 (+): a Cyanorak ORF ending at 4584 is
        in frame with both stops (diff 3 and 0) and inside both intervals."""
        outer = _make_ncbi_gff_df(
            locus_tag_ncbi="RS_OUTER", old_locus_tag="OUTERx", protein_id="WP_o",
            start=1, end=4587, strand="+",
        )
        inner = _make_ncbi_gff_df(
            locus_tag_ncbi="RS_INNER", old_locus_tag="INNERx", protein_id="WP_i",
            start=4432, end=4584, strand="+",
        )
        cyan_df = _make_cyan_gff_df(
            cyanorak_id="CK_X", locus_tag="CY_X", start=4000, end=4584, strand="+",
        )
        ncbi_df = pd.concat([outer, inner], ignore_index=True)
        result = self._run(ncbi_df, cyan_df, {"CK_X": "CY_X"})
        assert len(result) == 3
        ncbi_rows = result[result["locus_tag_ncbi"].notna()]
        assert len(ncbi_rows) == 2 and ncbi_rows["locus_tag_cyanorak"].isna().all()
        assert (result["locus_tag"] == "CY_X").sum() == 1  # stays Cyanorak-only

    def test_position_fallback_is_contig_aware(self):
        """Multi-contig assembly (PAC1/SB shape): Cyanorak stores the genome as
        ONE concatenated record, so NCBI contig-relative coordinates only
        compare after adding the contig's offset. The offset is derived from
        genes already matched by locus_tag on that contig. A Cyanorak gene at
        the raw (un-shifted) coordinates of a contig-2 NCBI gene must NOT
        merge; one at the shifted coordinates must."""
        # contig 2, offset 1000: an anchor gene matched by locus_tag ...
        anchor = _make_ncbi_gff_df(
            locus_tag_ncbi="EV03_RS00100", old_locus_tag="EV03_0100",
            protein_id="WP_1", start=100, end=400, strand="+",
        )
        anchor["seq_id"] = "contig2"
        anchor_cy = _make_cyan_gff_df(
            cyanorak_id="CK_0100", locus_tag="EV03_0100",
            start=1100, end=1400, strand="+",
        )
        # ... and two unmatched NCBI genes on contig 2
        n1 = _make_ncbi_gff_df(
            locus_tag_ncbi="EV03_RS00200", old_locus_tag="EV03_0200x",
            protein_id="WP_2", start=2000, end=2300, strand="+",
        )
        n1["seq_id"] = "contig2"
        n2 = _make_ncbi_gff_df(
            locus_tag_ncbi="EV03_RS00300", old_locus_tag="EV03_0300x",
            protein_id="WP_3", start=5000, end=5300, strand="+",
        )
        n2["seq_id"] = "contig2"
        # Cyanorak: one at n1's RAW coords (a coincidence from another contig),
        # one at n2's SHIFTED coords (the real counterpart)
        cy_raw = _make_cyan_gff_df(
            cyanorak_id="CK_0200", locus_tag="EV03_0200",
            start=2000, end=2300, strand="+",
        )
        cy_shift = _make_cyan_gff_df(
            cyanorak_id="CK_0300", locus_tag="EV03_0300",
            start=6000, end=6300, strand="+",
        )
        ncbi_df = pd.concat([anchor, n1, n2], ignore_index=True)
        cyan_df = pd.concat([anchor_cy, cy_raw, cy_shift], ignore_index=True)
        result = self._run(
            ncbi_df, cyan_df,
            {"CK_0100": "EV03_0100", "CK_0200": "EV03_0200", "CK_0300": "EV03_0300"},
        )
        by_ncbi = result.set_index("locus_tag_ncbi")
        assert by_ncbi.loc["EV03_RS00300", "locus_tag_cyanorak"] == "CK_0300"
        assert pd.isna(by_ncbi.loc["EV03_RS00200", "locus_tag_cyanorak"])
        # the raw-coordinate Cyanorak gene survives as a Cyanorak-only row
        assert (result["locus_tag"] == "EV03_0200").sum() == 1

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
        assert pd.isna(ncbi_row["locus_tag_cyanorak"])


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


# ─── Cyanorak-only coordinate fallback ───────────────────────────────────────

from multiomics_kg.download.build_gene_mapping import (  # noqa: E402
    _cyanorak_tiles,
    _ncbi_contig_lengths,
)


def _ncbi_genes(*specs: dict) -> pd.DataFrame:
    """Concatenate several ``_make_ncbi_gff_df`` genes; each spec may set ``seq_id``."""
    frames = []
    for spec in specs:
        seq_id = spec.pop("seq_id", "NC_005072.1")
        df = _make_ncbi_gff_df(**spec)
        df["seq_id"] = seq_id
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _cyan_genes(*specs: dict) -> pd.DataFrame:
    return pd.concat([_make_cyan_gff_df(**s) for s in specs], ignore_index=True)


def _write_ncbi_gff_header(path, lengths: dict[str, int]) -> str:
    lines = ["##gff-version 3"] + [
        f"##sequence-region {sid} 1 {length}" for sid, length in lengths.items()
    ]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


class TestCyanorakOnlyCoordinates:
    def _run(self, ncbi_gff_path, ncbi_df, cyan_df, cyan_id_map):
        with (
            patch(_LOAD_GFF, side_effect=[cyan_df, ncbi_df]),
            patch(_GET_CYAN_MAP, return_value=cyan_id_map),
        ):
            return load_gff_from_ncbi_and_cyanorak(ncbi_gff_path, "cyan.gff", "cyan.gbk")

    def test_ncbi_contig_lengths_parses_sequence_region_directives(self, tmp_path):
        gff = _write_ncbi_gff_header(
            tmp_path / "genomic.gff",
            {"NZ_JNAX01000001.1": 32549, "NZ_JNAX01000002.1": 30478},
        )
        assert _ncbi_contig_lengths(gff) == {
            "NZ_JNAX01000001.1": 32549,
            "NZ_JNAX01000002.1": 30478,
        }

    def test_ncbi_contig_lengths_missing_file_is_empty(self, tmp_path):
        assert _ncbi_contig_lengths(str(tmp_path / "nope.gff")) == {}

    def test_cyanorak_tiles_are_one_based_sorted_by_offset(self):
        # Cyanorak record = contig c1 (50,000 bp) + 50 bp spacer + contig c2 (30,000 bp)
        tiles = _cyanorak_tiles({"c2": 50050, "c1": 0}, {"c1": 50000, "c2": 30000})
        assert tiles == [(1, 50000, "c1"), (50051, 80050, "c2")]

    def test_cyanorak_tiles_unknown_length_single_contig_is_open_ended(self):
        tiles = _cyanorak_tiles({"NC_005072.1": 0}, {})
        assert tiles == [(1, float("inf"), "NC_005072.1")]

    def test_closed_genome_cyanorak_only_gene_gets_ncbi_contig_and_coordinates(self, tmp_path):
        gff = _write_ncbi_gff_header(tmp_path / "genomic.gff", {"NC_005072.1": 10000})
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001", start=174, end=1331)
        cyan_df = _cyan_genes(
            dict(cyanorak_id="CK_01", locus_tag="PMM0001", start=174, end=1331),
            dict(cyanorak_id="CK_99", locus_tag="PMM0099", start=5000, end=5300, strand="-"),
        )
        result = self._run(gff, ncbi_df, cyan_df, {"CK_01": "PMM0001", "CK_99": "PMM0099"})
        row = result.set_index("locus_tag").loc["PMM0099"]
        assert row["seqid"] == "NC_005072.1"
        assert row["start"] == 5000
        assert row["end"] == 5300
        assert row["strand"] == "-"
        assert row["coordinate_source"] == "cyanorak"
        assert result.set_index("locus_tag").loc["PMM0001", "coordinate_source"] == "ncbi"

    def test_draft_genome_cyanorak_only_gene_lands_on_second_contig(self, tmp_path):
        # Two 1,000 bp contigs; Cyanorak concatenates them with a 50 bp spacer,
        # so contig 2 sits at offset 1,050 in Cyanorak's coordinate space.
        gff = _write_ncbi_gff_header(tmp_path / "genomic.gff", {"c1": 1000, "c2": 1000})
        ncbi_df = _ncbi_genes(
            dict(seq_id="c1", locus_tag_ncbi="RS_1", old_locus_tag="PMM0001", start=100, end=400),
            dict(seq_id="c2", locus_tag_ncbi="RS_2", old_locus_tag="PMM0002", start=100, end=400),
        )
        cyan_df = _cyan_genes(
            dict(cyanorak_id="CK_01", locus_tag="PMM0001", start=100, end=400),
            dict(cyanorak_id="CK_02", locus_tag="PMM0002", start=1150, end=1450),
            dict(cyanorak_id="CK_03", locus_tag="PMM0003", start=1600, end=1900, strand="-"),
        )
        result = self._run(
            gff, ncbi_df, cyan_df,
            {"CK_01": "PMM0001", "CK_02": "PMM0002", "CK_03": "PMM0003"},
        )
        row = result.set_index("locus_tag").loc["PMM0003"]
        assert row["seqid"] == "c2"
        assert row["start"] == 550
        assert row["end"] == 850
        assert row["strand"] == "-"
        assert row["coordinate_source"] == "cyanorak"

    def test_cyanorak_only_gene_crossing_contig_boundary_is_left_unplaced(self, tmp_path):
        gff = _write_ncbi_gff_header(tmp_path / "genomic.gff", {"c1": 1000, "c2": 1000})
        ncbi_df = _ncbi_genes(
            dict(seq_id="c1", locus_tag_ncbi="RS_1", old_locus_tag="PMM0001", start=100, end=400),
            dict(seq_id="c2", locus_tag_ncbi="RS_2", old_locus_tag="PMM0002", start=100, end=400),
        )
        cyan_df = _cyan_genes(
            dict(cyanorak_id="CK_01", locus_tag="PMM0001", start=100, end=400),
            dict(cyanorak_id="CK_02", locus_tag="PMM0002", start=1150, end=1450),
            dict(cyanorak_id="CK_03", locus_tag="PMM0003", start=900, end=1200),  # spans the spacer
        )
        result = self._run(
            gff, ncbi_df, cyan_df,
            {"CK_01": "PMM0001", "CK_02": "PMM0002", "CK_03": "PMM0003"},
        )
        row = result.set_index("locus_tag").loc["PMM0003"]
        assert pd.isna(row["seqid"])
        assert pd.isna(row["start"])
        assert pd.isna(row["coordinate_source"])
        # the Cyanorak coordinates themselves are untouched
        assert row["start_cyanorak"] == 900

    # Neighbour-based placement (2026-10-28): Cyanorak records laid out unlike NCBI.

    @staticmethod
    def _paired(*genes):
        """Paired rows: (ncbi_start, ncbi_end, ncbi_strand, cyan_start, cyan_end, cyan_strand, seqid)."""
        return bgm._paired_relations(pd.DataFrame(
            genes, columns=["start_ncbi", "end_ncbi", "strand_ncbi",
                            "start_cyanorak", "end_cyanorak", "strand_cyanorak", "seqid"]))

    @staticmethod
    def _orphan(start, end, strand="+"):
        return pd.DataFrame([{"locus_tag": "CK_X", "start_cyanorak": start,
                              "end_cyanorak": end, "strand_cyanorak": strand}])

    def test_rotated_record_places_gene_by_its_own_segment(self):
        # Cyanorak starts the circle 600 bp later: genes before NCBI 600 sit at +400 in
        # Cyanorak, genes after it at -600. One tile offset would misplace one segment.
        paired = self._paired(
            (700, 800, "+", 100, 200, "+", "c1"), (900, 950, "+", 300, 350, "+", "c1"),
            (100, 200, "+", 500, 600, "+", "c1"), (300, 350, "+", 700, 750, "+", "c1"))
        cyan = self._orphan(400, 450)         # between pairs 2 and 3: a layout break
        assert bgm._fill_cyanorak_only_coordinates(cyan, [(1, 1000, "c1")], paired, {"c1": 1000}) == 0
        cyan = self._orphan(620, 680)         # between pairs 3 and 4: the +400 segment
        assert bgm._fill_cyanorak_only_coordinates(cyan, [(1, 1000, "c1")], paired, {"c1": 1000}) == 1
        assert (cyan.loc[0, "start_ncbi"], cyan.loc[0, "end_ncbi"]) == (220, 280)
        assert cyan.loc[0, "strand_ncbi"] == "+"

    def test_reverse_complemented_record_reflects_and_flips_strand(self):
        # Cyanorak = reverse complement of a 1,000 bp contig: cyan = 1001 - ncbi.
        paired = self._paired(
            (801, 900, "-", 101, 200, "+", "c1"), (601, 700, "+", 301, 400, "-", "c1"))
        cyan = self._orphan(221, 280, "+")
        assert bgm._fill_cyanorak_only_coordinates(cyan, [], paired, {"c1": 1000}) == 1
        assert (cyan.loc[0, "start_ncbi"], cyan.loc[0, "end_ncbi"]) == (721, 780)
        assert cyan.loc[0, "strand_ncbi"] == "-"

    def test_end_gene_needs_its_single_flank_to_agree_with_the_tile(self):
        paired = self._paired((100, 200, "+", 500, 600, "+", "c1"))   # shift 400
        cyan = self._orphan(700, 760)
        assert bgm._fill_cyanorak_only_coordinates(cyan, [(1, 1000, "c1")], paired) == 0
        cyan = self._orphan(700, 760)
        assert bgm._fill_cyanorak_only_coordinates(cyan, [(401, 1400, "c1")], paired) == 1
        assert cyan.loc[0, "start_ncbi"] == 300

    def test_placement_outside_the_contig_is_refused(self):
        paired = self._paired((100, 200, "+", 500, 600, "+", "c1"),
                              (300, 400, "+", 700, 800, "+", "c1"))
        cyan = self._orphan(650, 690)        # maps to 250-290, past the end of a 280 bp contig
        assert bgm._fill_cyanorak_only_coordinates(cyan, [], paired, {"c1": 280}) == 0

    def test_ncbi_only_mode_marks_every_gene_ncbi(self):
        ncbi_df = _make_ncbi_gff_df(old_locus_tag="PMM0001")
        with patch(_LOAD_GFF, return_value=ncbi_df):
            result = load_gff_from_ncbi_only("ncbi.gff")
        assert (result["coordinate_source"] == "ncbi").all()


# ─── protein-sequence fallback + merge QC (Ford 2021 batch, 2026-10-25) ───────

import numpy as np  # noqa: E402
from multiomics_kg.download import build_gene_mapping as bgm  # noqa: E402

_PROTS = {
    "WP_A": "M" + "ACDEFGHIKLMNPQRSTVWY" * 4,
    "WP_B": "M" + "WYVTSRQPNMLKIHGFEDCA" * 4,
    "WP_C": "M" + "KLMNPQRSTVACDEFGHIWY" * 4,
    "WP_P": "M" + "PQRSTVWYACDEFGHIKLMN" * 4,   # paralog-like, matches a far-away call
}


def _seq_frames(offset=100_000):
    """Three NCBI genes whose Cyanorak calls sit `offset` bp downstream (a
    re-laid-out record) plus one isolated candidate at an inconsistent offset."""
    ncbi = pd.DataFrame({
        "ID": [None, None, None, None],
        "protein_id": ["WP_A", "WP_B", "WP_C", "WP_P"],
        "locus_tag_ncbi": ["N1", "N2", "N3", "N4"],
        "locus_tag": [None, None, None, None],
        "seqid": ["c1"] * 4,
        "start_ncbi": [1000.0, 3000.0, 5000.0, 7000.0],
        "end_ncbi": [1243.0, 3243.0, 5243.0, 7243.0],
        "strand_ncbi": ["+"] * 4,
        "start_cyanorak": [np.nan] * 4, "end_cyanorak": [np.nan] * 4,
        "strand_cyanorak": [None] * 4, "product_cyanorak": [None] * 4,
    })
    cyan = pd.DataFrame({
        "locus_tag": ["C1", "C2", "C3", "C4"],
        "ID": ["CK_1", "CK_2", "CK_3", "CK_4"],
        "start_cyanorak": [1000.0 + offset, 3000.0 + offset, 5000.0 + offset, 400_000.0],
        "end_cyanorak": [1243.0 + offset, 3243.0 + offset, 5243.0 + offset, 400_243.0],
        "strand_cyanorak": ["+"] * 4,
        "product_cyanorak": ["a", "b", "c", "p"],
    })
    translations = {"C1": _PROTS["WP_A"], "C2": _PROTS["WP_B"], "C3": _PROTS["WP_C"], "C4": _PROTS["WP_P"]}
    return ncbi, cyan, translations


class TestSequenceFallback:
    def test_clean_translation(self):
        assert bgm._clean_translation("vkla*") == "MKLA"
        assert bgm._clean_translation("MKLA#") == "MKLA"
        assert bgm._clean_translation(None) is None
        assert bgm._clean_translation("MKL+AQ*") == "MKLXAQ"

    def test_syntenic_filter_keeps_consistent_run_rejects_isolated(self):
        ncbi, cyan, _ = _seq_frames()
        cands = {i: {"cyanorak_locus_tag": f"C{i + 1}"} for i in range(4)}
        keep = bgm._syntenic_candidates(ncbi, cyan, cands)
        assert keep == {0, 1, 2}   # C4 sits at its own offset: a paralog, not this locus

    def test_syntenic_filter_accepts_reverse_complemented_record(self):
        ncbi, cyan, _ = _seq_frames()
        L = 1_000_000.0
        # Cyanorak record reverse-complemented: c = L - n on the opposite strand
        cyan["start_cyanorak"], cyan["end_cyanorak"] = L - ncbi["end_ncbi"], L - ncbi["start_ncbi"]
        cyan["strand_cyanorak"] = "-"
        cands = {i: {"cyanorak_locus_tag": f"C{i + 1}"} for i in range(4)}
        assert bgm._syntenic_candidates(ncbi, cyan, cands) == {0, 1, 2, 3}

    def test_below_threshold_does_not_run(self, tmp_path):
        ncbi, cyan, tr = _seq_frames()
        n, consumed, qc = bgm._sequence_fallback_merge(ncbi, cyan, str(tmp_path / "x.faa"), tr)
        assert (n, consumed, qc["ran"]) == (0, set(), False)

    def test_end_to_end_merges_1to1_syntenic_only(self, tmp_path, monkeypatch):
        monkeypatch.setattr(bgm, "_SEQ_FALLBACK_MIN_UNPAIRED", 1)
        ncbi, cyan, tr = _seq_frames()
        faa = tmp_path / "protein.faa"
        faa.write_text("".join(f">{k}\n{v}\n" for k, v in _PROTS.items()))
        n, consumed, qc = bgm._sequence_fallback_merge(ncbi, cyan, str(faa), tr)
        assert n == 3 and consumed == {"C1", "C2", "C3"}
        assert qc["rejected_non_syntenic"] == 1
        assert qc["merged_by_method"] == {"exact": 3}
        assert ncbi.loc[0, "ID"] == "CK_1"
        assert ncbi.loc[0, "position_merge_note"] == "sequence_merge_exact:C1→N1"
        assert "C1" in str(ncbi.loc[0, "old_locus_tags"])
        assert pd.isna(ncbi.loc[3, "ID"])

    def test_merge_qc_flags_duplication(self):
        df = pd.DataFrame({"locus_tag": [f"g{i}" for i in range(180)]})
        qc = bgm._merge_qc(df, 100, 100, {}, {"ran": False})
        assert qc["duplication_ratio"] == 1.8
        assert any("duplication_ratio" in f for f in qc["flags"])
        ok = bgm._merge_qc(df.iloc[:105], 100, 100, {}, {"ran": False})
        assert ok["flags"] == []
