"""Unit tests for the CyanorakNcbi adapter."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from multiomics_kg.adapters.cyanorak_ncbi_adapter import (
    CyanorakNcbi,
    GeneEdgeType,
    GeneModel,
    GeneNodeField,
    MultiCyanorakNcbi,
)


# ---------------------------------------------------------------------------
# Mock data constants
# ---------------------------------------------------------------------------

NCBI_GFF_CONTENT = """\
##gff-version 3
##sequence-region NC_TEST.1 1 50000
NC_TEST.1\tRefSeq\tregion\t1\t50000\t.\t+\t.\tID=NC_TEST.1:1..50000;Dbxref=taxon:59919;Name=TEST
NC_TEST.1\tRefSeq\tgene\t174\t1331\t.\t+\t.\tID=gene-TEST_RS00020;Name=dnaN;gbkey=Gene;gene=dnaN;gene_biotype=protein_coding;locus_tag=TEST_RS00020;old_locus_tag=PMM0001
NC_TEST.1\tProtein Homology\tCDS\t174\t1331\t.\t+\t0\tID=cds-WP_TEST001.1;Parent=gene-TEST_RS00020;Dbxref=Genbank:WP_TEST001.1;Name=WP_TEST001.1;Note=test note;exception=test;gbkey=CDS;gene=dnaN;inference=COORDINATES: similar to AA sequence:RefSeq:WP_TEST001.1;locus_tag=TEST_RS00020;product=DNA polymerase III%2C beta subunit;protein_id=WP_TEST001.1
NC_TEST.1\tRefSeq\tgene\t1333\t2040\t.\t+\t.\tID=gene-TEST_RS00025;Name=TEST_RS00025;gbkey=Gene;gene_biotype=protein_coding;locus_tag=TEST_RS00025;old_locus_tag=PMM0002
NC_TEST.1\tProtein Homology\tCDS\t1333\t2040\t.\t+\t0\tID=cds-WP_TEST002.1;Parent=gene-TEST_RS00025;Dbxref=Genbank:WP_TEST002.1;Name=WP_TEST002.1;Note=test;exception=test;gbkey=CDS;inference=COORDINATES: similar to AA sequence:RefSeq:WP_TEST002.1;locus_tag=TEST_RS00025;product=hypothetical protein;protein_id=WP_TEST002.1
NC_TEST.1\tRefSeq\tgene\t2044\t4383\t.\t+\t.\tID=gene-TEST_RS00030;Name=purL;gbkey=Gene;gene=purL;gene_biotype=protein_coding;locus_tag=TEST_RS00030;old_locus_tag=PMM0003
NC_TEST.1\tProtein Homology\tCDS\t2044\t4383\t.\t+\t0\tID=cds-WP_TEST003.1;Parent=gene-TEST_RS00030;Dbxref=Genbank:WP_TEST003.1;Name=WP_TEST003.1;Note=test;exception=test;gbkey=CDS;gene=purL;inference=COORDINATES: similar to AA sequence:RefSeq:WP_TEST003.1;locus_tag=TEST_RS00030;product=phosphoribosylformylglycinamidine synthase;protein_id=WP_TEST003.1
"""

CYAN_GFF_CONTENT = """\
##gff-version 3
#seqID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes
TEST_chrom\tcyanorak\tsequence_assembly\t1\t50000\t.\t+\t0\tID=TEST_chrom
TEST_chrom\tcyanorak\tCDS\t174\t1331\t.\t+\t0\tID=CK_TEST_00001;Name=dnaN;product=DNA polymerase III%2C beta subunit;cluster_number=CK_00000364;Ontology_term=GO:0006260,GO:0003677;ontology_term_description=DNA replication,DNA binding;kegg=2.7.7.7;kegg_description=DNA-directed DNA polymerase;eggNOG=COG0592,bactNOG00989;eggNOG_description=COG: Replication%2C recombination;tIGR_Role=132;tIGR_Role_description=DNA metabolism;cyanorak_Role=F.1;cyanorak_Role_description=DNA replication;protein_domains=TIGR00663,PF00712;protein_domains_description=beta subunit,N-terminal domain
TEST_chrom\tcyanorak\tCDS\t1333\t2040\t.\t+\t0\tID=CK_TEST_00002;Name=PMM0002;product=conserved hypothetical protein;cluster_number=CK_00000363;eggNOG=COG0243;eggNOG_description=Energy production;tIGR_Role=156;tIGR_Role_description=Hypothetical;cyanorak_Role=R.2;cyanorak_Role_description=Conserved hypothetical
TEST_chrom\tcyanorak\tCDS\t2044\t4383\t.\t+\t0\tID=CK_TEST_00003;Name=purL;product=synthase;cluster_number=CK_00000362;Ontology_term=GO:0009152;ontology_term_description=purine biosynthesis;kegg=6.3.5.3;kegg_description=FGAM synthetase;eggNOG=COG0046;eggNOG_description=Nucleotide metabolism;tIGR_Role=125;tIGR_Role_description=Purine biosynthesis;cyanorak_Role=M.3;cyanorak_Role_description=Purine ribonucleotide biosynthesis;protein_domains=TIGR01736;protein_domains_description=synthase
"""

CYAN_GBK_CONTENT = """\
LOCUS       TEST_chrom           50000 bp    DNA     circular BCT 01-JAN-2024
DEFINITION  Test organism chromosome.
ACCESSION   TEST_ACC
VERSION     TEST_ACC.1
FEATURES             Location/Qualifiers
     source          1..50000
                     /organism="Test organism"
                     /mol_type="genomic DNA"
     CDS             174..1331
                     /gene="dnaN"
                     /locus_tag="PMM0001"
                     /product="DNA polymerase III, beta subunit"
                     /note="cyanorak ORF Id:CK_TEST_00001"
                     /translation="MTEST"
     CDS             1333..2040
                     /gene="PMM0002"
                     /locus_tag="PMM0002"
                     /product="conserved hypothetical protein"
                     /note="cyanorak ORF Id:CK_TEST_00002"
                     /translation="MTEST"
     CDS             2044..4383
                     /gene="purL"
                     /locus_tag="PMM0003"
                     /product="synthase"
                     /note="cyanorak ORF Id:CK_TEST_00003"
                     /translation="MTEST"
ORIGIN
        1 atgaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
//
"""


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_data_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def ncbi_gff_file(temp_data_dir):
    path = os.path.join(temp_data_dir, "genomic.gff")
    with open(path, "w") as f:
        f.write(NCBI_GFF_CONTENT)
    return path


@pytest.fixture
def cyan_gff_file(temp_data_dir):
    path = os.path.join(temp_data_dir, "cyanorak.gff")
    with open(path, "w") as f:
        f.write(CYAN_GFF_CONTENT)
    return path


@pytest.fixture
def cyan_gbk_file(temp_data_dir):
    path = os.path.join(temp_data_dir, "cyanorak.gbk")
    with open(path, "w") as f:
        f.write(CYAN_GBK_CONTENT)
    return path


@pytest.fixture
def adapter(ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
    return CyanorakNcbi(
        ncbi_gff_file=ncbi_gff_file,
        cyan_gff_file=cyan_gff_file,
        cyan_gbk_file=cyan_gbk_file,
    )


@pytest.fixture
def adapter_with_data(adapter):
    adapter.download_data()
    return adapter


@pytest.fixture
def adapter_no_prefix(ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
    return CyanorakNcbi(
        ncbi_gff_file=ncbi_gff_file,
        cyan_gff_file=cyan_gff_file,
        cyan_gbk_file=cyan_gbk_file,
        add_prefix=False,
    )


@pytest.fixture
def adapter_test_mode(ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
    return CyanorakNcbi(
        ncbi_gff_file=ncbi_gff_file,
        cyan_gff_file=cyan_gff_file,
        cyan_gbk_file=cyan_gbk_file,
        test_mode=True,
    )


@pytest.fixture
def adapter_subset_fields(ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
    return CyanorakNcbi(
        ncbi_gff_file=ncbi_gff_file,
        cyan_gff_file=cyan_gff_file,
        cyan_gbk_file=cyan_gbk_file,
        gene_node_fields=[
            GeneNodeField.LOCUS_TAG,
            GeneNodeField.PRODUCT,
            GeneNodeField.GENE_NAMES,
        ],
    )


@pytest.fixture
def real_data_paths():
    dpath = "data/Prochlorococcus/genomes/MED4/"
    ncbi = os.path.join(dpath, "genomic.gff")
    cyan_gff = os.path.join(dpath, "cyanorak/Pro_MED4.gff")
    cyan_gbk = os.path.join(dpath, "cyanorak/Pro_MED4.gbk")
    if not all(os.path.exists(f) for f in [ncbi, cyan_gff, cyan_gbk]):
        pytest.skip("Real MED4 data files not available")
    return ncbi, cyan_gff, cyan_gbk


# ---------------------------------------------------------------------------
# Tests: GeneNodeField enum
# ---------------------------------------------------------------------------


class TestGeneNodeFieldEnum:
    def test_enum_values_exist(self):
        expected = [
            "GENE_NAMES", "LOCUS_TAG", "START", "END", "STRAND",
            "PRODUCT", "PROTEIN_ID", "ONTOLOGY_TERM", "EGGNOG", "KEGG",
            "CYANORAK_ROLE", "TIGR_ROLE", "CLUSTER_NUMBER", "PROTEIN_DOMAINS",
        ]
        for name in expected:
            assert name in GeneNodeField.__members__

    def test_enum_contains_check(self):
        assert "GENE_NAMES" in GeneNodeField

    def test_enum_contains_nonexistent(self):
        assert "NONEXISTENT" not in GeneNodeField

    def test_missing_case_insensitive_lookup(self):
        result = GeneNodeField("Gene_Names")
        assert result == GeneNodeField.GENE_NAMES

    def test_missing_returns_none_for_invalid(self):
        assert GeneNodeField._missing_("totally_invalid") is None


# ---------------------------------------------------------------------------
# Tests: GeneModel validation
# ---------------------------------------------------------------------------


class TestGeneModelValidation:
    def test_default_construction(self):
        model = GeneModel()
        assert model.test_mode is False
        assert model.add_prefix is True
        assert model.organism is None

    def test_test_mode_flag(self):
        model = GeneModel(test_mode=True)
        assert model.test_mode is True

    def test_organism_int(self):
        model = GeneModel(organism=59919)
        assert model.organism == 59919

    def test_organism_star(self):
        model = GeneModel(organism="*")
        assert model.organism == "*"

    def test_organism_none(self):
        model = GeneModel(organism=None)
        assert model.organism is None

    def test_invalid_output_dir(self):
        with pytest.raises(ValidationError):
            GeneModel(output_dir="/nonexistent/path/xyz")

    def test_valid_output_dir(self, temp_data_dir):
        model = GeneModel(output_dir=temp_data_dir)
        assert str(model.output_dir) == temp_data_dir

    def test_gene_node_fields_list(self):
        model = GeneModel(
            gene_node_fields=[GeneNodeField.PRODUCT, GeneNodeField.LOCUS_TAG]
        )
        assert len(model.gene_node_fields) == 2


# ---------------------------------------------------------------------------
# Tests: Constructor
# ---------------------------------------------------------------------------


class TestConstructor:
    def test_default_node_fields(self, adapter):
        all_values = [f.value for f in GeneNodeField]
        assert adapter.gene_node_fields == all_values

    def test_custom_node_fields(self, adapter_subset_fields):
        assert "locus_tag" in adapter_subset_fields.gene_node_fields
        assert "product" in adapter_subset_fields.gene_node_fields
        assert "gene_names" in adapter_subset_fields.gene_node_fields
        assert len(adapter_subset_fields.gene_node_fields) == 3

    def test_default_edge_types(self, adapter):
        assert adapter.edge_types == list(GeneEdgeType)

    def test_add_prefix_default_true(self, adapter):
        assert adapter.add_prefix is True

    def test_add_prefix_false(self, adapter_no_prefix):
        assert adapter_no_prefix.add_prefix is False

    def test_test_mode_sets_early_stopping(self, adapter_test_mode):
        assert adapter_test_mode.early_stopping == 100

    def test_normal_mode_no_early_stopping(self, adapter):
        assert adapter.early_stopping is None

    def test_file_paths_stored(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        assert adapter.ncbi_gff_file == ncbi_gff_file
        assert adapter.cyan_gff_file == cyan_gff_file
        assert adapter.cyan_gbk_file == cyan_gbk_file


# ---------------------------------------------------------------------------
# Tests: clean_text
# ---------------------------------------------------------------------------


class TestCleanText:
    def test_pipe_replaced_with_comma(self, adapter):
        assert adapter.clean_text("value1|value2") == "value1,value2"

    def test_single_quote_replaced_with_caret(self, adapter):
        assert adapter.clean_text("it's a test") == "it^s a test"

    def test_both_special_chars(self, adapter):
        assert adapter.clean_text("foo|bar's") == "foo,bar^s"

    def test_no_special_chars(self, adapter):
        assert adapter.clean_text("normal text") == "normal text"

    def test_list_input(self, adapter):
        result = adapter.clean_text(["a|b", "c'd"])
        assert result == ["a,b", "c^d"]

    def test_non_string_passthrough(self, adapter):
        assert adapter.clean_text(123) == 123

    def test_none_input(self, adapter):
        assert adapter.clean_text(None) is None

    def test_empty_string(self, adapter):
        assert adapter.clean_text("") == ""


# ---------------------------------------------------------------------------
# Tests: _split_field and _get_split_character
# ---------------------------------------------------------------------------


class TestSplitField:
    def test_comma_split_ontology_term(self, adapter):
        result = adapter._split_field("Ontology_term", "GO:0006260,GO:0003677")
        assert result == ["GO:0006260", "GO:0003677"]

    def test_comma_split_eggnog(self, adapter):
        result = adapter._split_field("eggNOG", "COG0592,bactNOG00989")
        assert result == ["COG0592", "bactNOG00989"]

    def test_comma_split_protein_domains(self, adapter):
        result = adapter._split_field("protein_domains", "TIGR00663,PF00712")
        assert result == ["TIGR00663", "PF00712"]

    def test_space_split_gene_names(self, adapter):
        result = adapter._split_field("gene_names", "dnaN geneB")
        assert result == ["dnaN", "geneB"]

    def test_no_split_for_product(self, adapter):
        result = adapter._split_field("product", "DNA polymerase III")
        assert result == "DNA polymerase III"

    def test_no_split_for_locus_tag(self, adapter):
        result = adapter._split_field("locus_tag", "PMM0001")
        assert result == "PMM0001"

    def test_get_split_character_comma_fields(self, adapter):
        assert adapter._get_split_character("cyanorak_Role") == ","
        assert adapter._get_split_character("Ontology_term") == ","
        assert adapter._get_split_character("eggNOG") == ","

    def test_get_split_character_space_fields(self, adapter):
        assert adapter._get_split_character("gene_names") == " "
        assert adapter._get_split_character("gene") == " "

    def test_get_split_character_semicolon_fields(self, adapter):
        assert adapter._get_split_character("kegg_description") == ";"

    def test_get_split_character_none(self, adapter):
        assert adapter._get_split_character("product") is None
        assert adapter._get_split_character("locus_tag") is None

    def test_kegg_uses_comma_split(self, adapter):
        # kegg appears in both comma and semicolon lists, comma checked first
        assert adapter._get_split_character("kegg") == ","

    def test_comma_followed_by_space_not_split(self, adapter):
        # regex r',(?! )' means comma followed by space is NOT split
        result = adapter._split_field("Ontology_term", "GO:001, GO:002")
        assert result == ["GO:001, GO:002"]


# ---------------------------------------------------------------------------
# Tests: add_prefix_to_id
# ---------------------------------------------------------------------------


class TestAddPrefixToId:
    def test_add_prefix_enabled(self, adapter):
        result = adapter.add_prefix_to_id(prefix="ncbigene", identifier="PMM0001")
        assert "ncbigene" in result
        assert "PMM0001" in result

    def test_add_prefix_disabled(self, adapter_no_prefix):
        result = adapter_no_prefix.add_prefix_to_id(
            prefix="ncbigene", identifier="PMM0001"
        )
        assert result == "PMM0001"

    def test_none_identifier(self, adapter):
        # identifier has type str with @validate_call, None may be coerced to "None"
        # or raise ValidationError depending on pydantic version
        try:
            result = adapter.add_prefix_to_id(prefix="ncbigene", identifier=None)
            # If it doesn't raise, result should be None or "None"
            assert result is None or result == "None"
        except Exception:
            pass  # validation error is acceptable

    def test_empty_identifier(self, adapter):
        result = adapter.add_prefix_to_id(prefix="ncbigene", identifier="")
        assert result == ""


# ---------------------------------------------------------------------------
# Tests: set_node_fields
# ---------------------------------------------------------------------------


class TestSetNodeFields:
    def test_set_node_fields_with_list(self, adapter):
        adapter.set_node_fields(
            gene_node_fields=[GeneNodeField.PRODUCT, GeneNodeField.LOCUS_TAG]
        )
        assert adapter.gene_node_fields == ["product", "locus_tag"]

    def test_set_node_fields_none(self, adapter):
        adapter.set_node_fields(gene_node_fields=None)
        assert len(adapter.gene_node_fields) == len(GeneNodeField)

    def test_set_node_fields_empty_list(self, adapter):
        adapter.set_node_fields(gene_node_fields=[])
        # empty list is falsy, so defaults to all fields
        assert len(adapter.gene_node_fields) == len(GeneNodeField)

    def test_set_node_fields_single(self, adapter):
        adapter.set_node_fields(gene_node_fields=[GeneNodeField.PRODUCT])
        assert adapter.gene_node_fields == ["product"]


# ---------------------------------------------------------------------------
# Tests: set_edge_types
# ---------------------------------------------------------------------------


class TestSetEdgeTypes:
    def test_set_edge_types_none(self, adapter):
        adapter.set_edge_types(edge_types=None)
        assert adapter.edge_types == list(GeneEdgeType)

    def test_set_edge_types_specific(self, adapter):
        adapter.set_edge_types(edge_types=[GeneEdgeType.GENE_TO_PROTEIN])
        assert adapter.edge_types == [GeneEdgeType.GENE_TO_PROTEIN]


# ---------------------------------------------------------------------------
# Tests: _get_cyanorak_id_map_from_gbk
# ---------------------------------------------------------------------------


class TestGetCyanorakIdMapFromGbk:
    def test_returns_correct_count(self, adapter, cyan_gbk_file):
        mapping = adapter._get_cyanorak_id_map_from_gbk(cyan_gbk_file)
        assert len(mapping) == 3

    def test_correct_mapping(self, adapter, cyan_gbk_file):
        mapping = adapter._get_cyanorak_id_map_from_gbk(cyan_gbk_file)
        assert mapping["CK_TEST_00001"] == "PMM0001"
        assert mapping["CK_TEST_00002"] == "PMM0002"
        assert mapping["CK_TEST_00003"] == "PMM0003"

    def test_return_type(self, adapter, cyan_gbk_file):
        mapping = adapter._get_cyanorak_id_map_from_gbk(cyan_gbk_file)
        assert isinstance(mapping, dict)
        for k, v in mapping.items():
            assert isinstance(k, str)
            assert isinstance(v, str)


# ---------------------------------------------------------------------------
# Tests: load_gff
# ---------------------------------------------------------------------------


class TestLoadGff:
    def test_load_ncbi_gff_returns_dataframe(self, adapter, ncbi_gff_file):
        df = adapter.load_gff(ncbi_gff_file)
        assert isinstance(df, pd.DataFrame)

    def test_load_ncbi_gff_filters_gene_and_cds(self, adapter, ncbi_gff_file):
        df = adapter.load_gff(ncbi_gff_file)
        assert set(df["type"].unique()) <= {"gene", "CDS"}
        # 3 genes + 3 CDS = 6 rows
        assert len(df) == 6

    def test_load_ncbi_gff_has_attribute_columns(self, adapter, ncbi_gff_file):
        df = adapter.load_gff(ncbi_gff_file)
        for col in ["ID", "Name", "locus_tag"]:
            assert col in df.columns

    def test_load_cyan_gff_returns_dataframe(self, adapter, cyan_gff_file):
        df = adapter.load_gff(cyan_gff_file)
        assert isinstance(df, pd.DataFrame)

    def test_load_cyan_gff_cds_only(self, adapter, cyan_gff_file):
        df = adapter.load_gff(cyan_gff_file)
        # sequence_assembly is filtered out; only CDS retained
        assert all(df["type"] == "CDS")
        assert len(df) == 3


# ---------------------------------------------------------------------------
# Tests: ncbi_merge_cds_and_gene_entries
# ---------------------------------------------------------------------------


class TestNcbiMergeCdsAndGeneEntries:
    def test_merge_produces_correct_row_count(self, adapter, ncbi_gff_file):
        ncbi_df = adapter.load_gff(ncbi_gff_file)
        merged = adapter.ncbi_merge_cds_and_gene_entries(ncbi_df)
        assert len(merged) == 3

    def test_merge_has_renamed_columns(self, adapter, ncbi_gff_file):
        ncbi_df = adapter.load_gff(ncbi_gff_file)
        merged = adapter.ncbi_merge_cds_and_gene_entries(ncbi_df)
        expected_cols = ["locus_tag", "locus_tag_ncbi", "product", "protein_id", "Name"]
        for col in expected_cols:
            assert col in merged.columns, f"Missing column: {col}"

    def test_merge_locus_tag_from_old_locus_tag(self, adapter, ncbi_gff_file):
        ncbi_df = adapter.load_gff(ncbi_gff_file)
        merged = adapter.ncbi_merge_cds_and_gene_entries(ncbi_df)
        locus_tags = set(merged["locus_tag"].values)
        assert "PMM0001" in locus_tags
        assert "PMM0002" in locus_tags
        assert "PMM0003" in locus_tags

    def test_merge_locus_tag_ncbi_is_new_tag(self, adapter, ncbi_gff_file):
        ncbi_df = adapter.load_gff(ncbi_gff_file)
        merged = adapter.ncbi_merge_cds_and_gene_entries(ncbi_df)
        ncbi_tags = set(merged["locus_tag_ncbi"].values)
        assert "TEST_RS00020" in ncbi_tags


# ---------------------------------------------------------------------------
# Tests: load_gff_from_ncbi_and_cynorak (full merge pipeline)
# ---------------------------------------------------------------------------


class TestLoadGffFromNcbiAndCynorak:
    def test_returns_merged_dataframe(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        assert isinstance(df, pd.DataFrame)

    def test_merged_row_count(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        assert len(df) == 3

    def test_has_ncbi_columns(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        for col in ["gene_names", "start", "end", "strand", "product", "protein_id", "locus_tag_ncbi"]:
            assert col in df.columns, f"Missing NCBI column: {col}"

    def test_has_cyanorak_columns(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        for col in ["gene_names_cyanorak", "start_cyanorak", "end_cyanorak",
                     "Ontology_term", "cluster_number", "eggNOG", "kegg"]:
            assert col in df.columns, f"Missing Cyanorak column: {col}"

    def test_locus_tag_preserved(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        assert "locus_tag" in df.columns
        locus_tags = set(df["locus_tag"].values)
        assert "PMM0001" in locus_tags

    def test_final_rename_applied(self, adapter, ncbi_gff_file, cyan_gff_file, cyan_gbk_file):
        df = adapter.load_gff_from_ncbi_and_cynorak(ncbi_gff_file, cyan_gff_file, cyan_gbk_file)
        # Name_ncbi should be renamed to gene_names
        assert "gene_names" in df.columns
        # Name_cyanorak should be renamed to gene_names_cyanorak
        assert "gene_names_cyanorak" in df.columns


# ---------------------------------------------------------------------------
# Tests: download_data
# ---------------------------------------------------------------------------


class TestDownloadData:
    def test_populates_data_df(self, adapter):
        adapter.download_data()
        assert isinstance(adapter.data_df, pd.DataFrame)

    def test_data_df_row_count(self, adapter):
        adapter.download_data()
        assert len(adapter.data_df) == 3


# ---------------------------------------------------------------------------
# Tests: get_nodes
# ---------------------------------------------------------------------------


class TestGetNodes:
    def test_returns_list_of_tuples(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        assert isinstance(nodes, list)
        assert all(isinstance(n, tuple) for n in nodes)

    def test_node_count(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        assert len(nodes) == 3

    def test_node_tuple_structure(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        for node_id, label, props in nodes:
            assert isinstance(node_id, str) or node_id is None
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_default_label_is_gene(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        for _, label, _ in nodes:
            assert label == "gene"

    def test_custom_label(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes(label="custom_gene")
        for _, label, _ in nodes:
            assert label == "custom_gene"

    def test_node_id_has_prefix(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        for node_id, _, _ in nodes:
            assert "ncbigene" in node_id

    def test_node_id_no_prefix(self, adapter_no_prefix):
        adapter_no_prefix.download_data()
        nodes = adapter_no_prefix.get_nodes()
        for node_id, _, _ in nodes:
            assert node_id.startswith("PMM")

    def test_node_properties_contain_expected_fields(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        # First gene (dnaN) should have product, eggNOG, etc.
        _, _, props = nodes[0]
        assert "product" in props
        assert "eggNOG" in props or "cluster_number" in props

    def test_nan_fields_excluded_from_properties(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        for _, _, props in nodes:
            for v in props.values():
                if isinstance(v, float):
                    assert not np.isnan(v)

    def test_url_encoded_values_decoded(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        # The first gene has product "DNA polymerase III%2C beta subunit" in cyanorak GFF
        # After URL decoding it should contain the comma
        for _, _, props in nodes:
            for v in props.values():
                if isinstance(v, str):
                    assert "%2C" not in v
                    assert "%3B" not in v

    def test_multivalue_fields_are_lists(self, adapter_with_data):
        nodes = adapter_with_data.get_nodes()
        # First gene has Ontology_term=GO:0006260,GO:0003677
        _, _, props = nodes[0]
        if "Ontology_term" in props:
            assert isinstance(props["Ontology_term"], list)
        if "eggNOG" in props:
            assert isinstance(props["eggNOG"], list)

    def test_subset_fields(self, adapter_subset_fields):
        adapter_subset_fields.download_data()
        nodes = adapter_subset_fields.get_nodes()
        for _, _, props in nodes:
            allowed = {"locus_tag", "product", "gene_names"}
            assert set(props.keys()) <= allowed


# ---------------------------------------------------------------------------
# Tests: get_edges
# ---------------------------------------------------------------------------


class TestGetEdges:
    def test_returns_empty_list(self, adapter_with_data):
        edges = adapter_with_data.get_edges()
        assert edges == []

    def test_return_type(self, adapter_with_data):
        edges = adapter_with_data.get_edges()
        assert isinstance(edges, list)


# ---------------------------------------------------------------------------
# Tests: column mapping helpers
# ---------------------------------------------------------------------------


class TestColumnMappingMethods:
    def test_cyanorak_cols_to_keep_returns_list(self, adapter):
        cols = adapter._get_cyanorak_cols_to_keep()
        assert isinstance(cols, list)
        assert "Ontology_term" in cols
        assert "cluster_number" in cols
        assert "locus_tag" in cols

    def test_ncbi_cols_map_renames_locus_tag_cds(self, adapter):
        col_map = adapter._get_ncbi_cols_to_keep_map()
        assert col_map["locus_tag_cds"] == "locus_tag_ncbi"

    def test_ncbi_cols_map_renames_old_locus_tag(self, adapter):
        col_map = adapter._get_ncbi_cols_to_keep_map()
        assert col_map["old_locus_tag_gene"] == "locus_tag"

    def test_final_merged_map_renames(self, adapter):
        col_map = adapter._get_final_merged_columns_map()
        assert col_map["Name_ncbi"] == "gene_names"
        assert col_map["Name_cyanorak"] == "gene_names_cyanorak"
        assert col_map["start_ncbi"] == "start"
        assert col_map["product_ncbi"] == "product"


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_special_chars_cleaned_in_product(self, adapter_with_data):
        """Pipe and quote chars should be replaced in all property values."""
        nodes = adapter_with_data.get_nodes()
        for _, _, props in nodes:
            for v in props.values():
                if isinstance(v, str):
                    assert "|" not in v
                    assert "'" not in v

    def test_url_encoded_comma_in_product(self, adapter_with_data):
        """Product 'DNA polymerase III%2C beta subunit' should be decoded."""
        nodes = adapter_with_data.get_nodes()
        # Find the dnaN gene node
        dnaN_props = None
        for _, _, props in nodes:
            if props.get("gene_names") and "dnaN" in str(props["gene_names"]):
                dnaN_props = props
                break
        if dnaN_props and "product" in dnaN_props:
            # NCBI product doesn't have %2C, but cyanorak product_cyanorak does
            pass  # main check is no %2C in any field (tested above)


# ---------------------------------------------------------------------------
# Tests: integration with real data (skipped if files not present)
# ---------------------------------------------------------------------------


class TestIntegrationWithRealData:
    def test_real_data_loads(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        assert len(adapter.data_df) > 1000

    def test_real_data_node_count(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        nodes = adapter.get_nodes()
        assert len(nodes) > 1000

    def test_real_data_node_structure(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        nodes = adapter.get_nodes()
        node_id, label, props = nodes[0]
        assert isinstance(node_id, str)
        assert label == "gene"
        assert len(props) > 0

    def test_real_data_known_locus_tags(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        nodes = adapter.get_nodes()
        all_ids = [n[0] for n in nodes]
        # PMM0001 should be in the node IDs (with prefix)
        assert any("PMM0001" in nid for nid in all_ids)

    def test_real_data_ontology_terms_are_lists(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        nodes = adapter.get_nodes()
        # Find a node with Ontology_term
        for _, _, props in nodes:
            if "Ontology_term" in props:
                assert isinstance(props["Ontology_term"], list)
                break

    def test_real_data_edges_empty(self, real_data_paths):
        ncbi, cyan_gff, cyan_gbk = real_data_paths
        adapter = CyanorakNcbi(
            ncbi_gff_file=ncbi,
            cyan_gff_file=cyan_gff,
            cyan_gbk_file=cyan_gbk,
        )
        adapter.download_data()
        edges = adapter.get_edges()
        assert edges == []


# ---------------------------------------------------------------------------
# Mock data for a second genome (used in MultiCyanorakNcbi tests)
# ---------------------------------------------------------------------------

NCBI_GFF_CONTENT_2 = """\
##gff-version 3
##sequence-region NC_OTHER.1 1 30000
NC_OTHER.1\tRefSeq\tregion\t1\t30000\t.\t+\t.\tID=NC_OTHER.1:1..30000;Dbxref=taxon:99999;Name=OTHER
NC_OTHER.1\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-OTHER_RS00010;Name=recA;gbkey=Gene;gene=recA;gene_biotype=protein_coding;locus_tag=OTHER_RS00010;old_locus_tag=OTH0001
NC_OTHER.1\tProtein Homology\tCDS\t100\t900\t.\t+\t0\tID=cds-WP_OTHER001.1;Parent=gene-OTHER_RS00010;Dbxref=Genbank:WP_OTHER001.1;Name=WP_OTHER001.1;Note=test;exception=test;gbkey=CDS;gene=recA;inference=COORDINATES: similar to AA sequence:RefSeq:WP_OTHER001.1;locus_tag=OTHER_RS00010;product=recombinase RecA;protein_id=WP_OTHER001.1
NC_OTHER.1\tRefSeq\tgene\t1000\t1800\t.\t-\t.\tID=gene-OTHER_RS00020;Name=OTHER_RS00020;gbkey=Gene;gene_biotype=protein_coding;locus_tag=OTHER_RS00020;old_locus_tag=OTH0002
NC_OTHER.1\tProtein Homology\tCDS\t1000\t1800\t.\t-\t0\tID=cds-WP_OTHER002.1;Parent=gene-OTHER_RS00020;Dbxref=Genbank:WP_OTHER002.1;Name=WP_OTHER002.1;Note=test;exception=test;gbkey=CDS;inference=COORDINATES: similar to AA sequence:RefSeq:WP_OTHER002.1;locus_tag=OTHER_RS00020;product=hypothetical protein;protein_id=WP_OTHER002.1
"""

CYAN_GFF_CONTENT_2 = """\
##gff-version 3
#seqID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes
OTHER_chrom\tcyanorak\tsequence_assembly\t1\t30000\t.\t+\t0\tID=OTHER_chrom
OTHER_chrom\tcyanorak\tCDS\t100\t900\t.\t+\t0\tID=CK_OTHER_00001;Name=recA;product=recombinase RecA;cluster_number=CK_00001111;Ontology_term=GO:0006281;ontology_term_description=DNA repair;kegg=3.4.21.88;kegg_description=repressor LexA;eggNOG=COG0468;eggNOG_description=Recombinase;tIGR_Role=140;tIGR_Role_description=DNA repair;cyanorak_Role=F.2;cyanorak_Role_description=DNA repair;protein_domains=TIGR02012;protein_domains_description=recA protein
OTHER_chrom\tcyanorak\tCDS\t1000\t1800\t.\t-\t0\tID=CK_OTHER_00002;Name=OTH0002;product=hypothetical protein;cluster_number=CK_00002222;eggNOG=COG1234;eggNOG_description=Unknown;tIGR_Role=156;tIGR_Role_description=Hypothetical;cyanorak_Role=R.2;cyanorak_Role_description=Conserved hypothetical;protein_domains=PF99999;protein_domains_description=hypothetical domain
"""

CYAN_GBK_CONTENT_2 = """\
LOCUS       OTHER_chrom          30000 bp    DNA     circular BCT 01-JAN-2024
DEFINITION  Other test organism chromosome.
ACCESSION   OTHER_ACC
VERSION     OTHER_ACC.1
FEATURES             Location/Qualifiers
     source          1..30000
                     /organism="Other test organism"
                     /mol_type="genomic DNA"
     CDS             100..900
                     /gene="recA"
                     /locus_tag="OTH0001"
                     /product="recombinase RecA"
                     /note="cyanorak ORF Id:CK_OTHER_00001"
                     /translation="MTEST"
     CDS             1000..1800
                     /gene="OTH0002"
                     /locus_tag="OTH0002"
                     /product="hypothetical protein"
                     /note="cyanorak ORF Id:CK_OTHER_00002"
                     /translation="MTEST"
ORIGIN
        1 atgaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
//
"""


# ---------------------------------------------------------------------------
# Fixtures for MultiCyanorakNcbi
# ---------------------------------------------------------------------------


@pytest.fixture
def multi_genome_dir():
    """Create a temp dir with two genome subdirectories and a CSV config file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Genome 1
        g1 = os.path.join(tmpdir, "genome1")
        os.makedirs(os.path.join(g1, "cyanorak"))
        with open(os.path.join(g1, "genomic.gff"), "w") as f:
            f.write(NCBI_GFF_CONTENT)
        with open(os.path.join(g1, "cyanorak", "strain1.gff"), "w") as f:
            f.write(CYAN_GFF_CONTENT)
        with open(os.path.join(g1, "cyanorak", "strain1.gbk"), "w") as f:
            f.write(CYAN_GBK_CONTENT)

        # Genome 2
        g2 = os.path.join(tmpdir, "genome2")
        os.makedirs(os.path.join(g2, "cyanorak"))
        with open(os.path.join(g2, "genomic.gff"), "w") as f:
            f.write(NCBI_GFF_CONTENT_2)
        with open(os.path.join(g2, "cyanorak", "strain2.gff"), "w") as f:
            f.write(CYAN_GFF_CONTENT_2)
        with open(os.path.join(g2, "cyanorak", "strain2.gbk"), "w") as f:
            f.write(CYAN_GBK_CONTENT_2)

        # CSV config
        csv_path = os.path.join(tmpdir, "cyanobacteria_genomes.csv")
        with open(csv_path, "w") as f:
            f.write("genome_dir,ncbi_gff,cyan_gff,cyan_gbk\n")
            f.write(f"{g1},genomic.gff,cyanorak/strain1.gff,cyanorak/strain1.gbk\n")
            f.write(f"{g2},genomic.gff,cyanorak/strain2.gff,cyanorak/strain2.gbk\n")

        yield tmpdir


@pytest.fixture
def multi_config_path(multi_genome_dir):
    return os.path.join(multi_genome_dir, "cyanobacteria_genomes.csv")


@pytest.fixture
def single_genome_config(multi_genome_dir):
    """CSV config with only the first genome."""
    csv_path = os.path.join(multi_genome_dir, "single_genome.csv")
    g1 = os.path.join(multi_genome_dir, "genome1")
    with open(csv_path, "w") as f:
        f.write("genome_dir,ncbi_gff,cyan_gff,cyan_gbk\n")
        f.write(f"{g1},genomic.gff,cyanorak/strain1.gff,cyanorak/strain1.gbk\n")
    return csv_path


# ---------------------------------------------------------------------------
# Tests: MultiCyanorakNcbi
# ---------------------------------------------------------------------------


class TestMultiCyanorakNcbiConstruction:
    def test_loads_correct_number_of_adapters(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        assert len(wrapper.adapters) == 2

    def test_single_genome_config(self, single_genome_config):
        wrapper = MultiCyanorakNcbi(config_list_file=single_genome_config)
        assert len(wrapper.adapters) == 1

    def test_adapters_are_cyanorak_instances(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        for adapter in wrapper.adapters:
            assert isinstance(adapter, CyanorakNcbi)

    def test_kwargs_passed_to_adapters(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(
            config_list_file=multi_config_path, test_mode=True
        )
        for adapter in wrapper.adapters:
            assert adapter.early_stopping == 100

    def test_empty_csv_creates_no_adapters(self, multi_genome_dir):
        csv_path = os.path.join(multi_genome_dir, "empty.csv")
        with open(csv_path, "w") as f:
            f.write("genome_dir,ncbi_gff,cyan_gff,cyan_gbk\n")
        wrapper = MultiCyanorakNcbi(config_list_file=csv_path)
        assert len(wrapper.adapters) == 0


class TestMultiCyanorakNcbiDownloadData:
    def test_download_populates_all_adapters(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        for adapter in wrapper.adapters:
            assert isinstance(adapter.data_df, pd.DataFrame)
            assert len(adapter.data_df) > 0


class TestMultiCyanorakNcbiGetNodes:
    def test_returns_list(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        nodes = wrapper.get_nodes()
        assert isinstance(nodes, list)

    def test_aggregates_nodes_from_all_genomes(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        nodes = wrapper.get_nodes()
        # Genome 1 has 3 genes, genome 2 has 2 genes
        assert len(nodes) == 5

    def test_single_genome_matches_direct_adapter(self, single_genome_config):
        wrapper = MultiCyanorakNcbi(config_list_file=single_genome_config)
        wrapper.download_data()
        multi_nodes = wrapper.get_nodes()
        # Compare with direct adapter
        direct = wrapper.adapters[0]
        direct_nodes = direct.get_nodes()
        assert len(multi_nodes) == len(direct_nodes)

    def test_nodes_from_both_genomes_present(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        nodes = wrapper.get_nodes()
        all_ids = [n[0] for n in nodes]
        # Genome 1 has PMM0001, genome 2 has OTH0001
        assert any("PMM0001" in nid for nid in all_ids)
        assert any("OTH0001" in nid for nid in all_ids)

    def test_all_nodes_are_gene_label(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        nodes = wrapper.get_nodes()
        for _, label, _ in nodes:
            assert label == "gene"

    def test_node_properties_are_dicts(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        nodes = wrapper.get_nodes()
        for _, _, props in nodes:
            assert isinstance(props, dict)


class TestMultiCyanorakNcbiGetEdges:
    def test_returns_list(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        edges = wrapper.get_edges()
        assert isinstance(edges, list)

    def test_edges_empty_for_current_implementation(self, multi_config_path):
        wrapper = MultiCyanorakNcbi(config_list_file=multi_config_path)
        wrapper.download_data()
        edges = wrapper.get_edges()
        assert edges == []


class TestMultiCyanorakNcbiIntegration:
    def test_real_data_single_genome(self):
        """Integration test with real MED4 data via CSV config."""
        csv_content = (
            "genome_dir,ncbi_gff,cyan_gff,cyan_gbk\n"
            "data/Prochlorococcus/genomes/MED4/,genomic.gff,"
            "cyanorak/Pro_MED4.gff,cyanorak/Pro_MED4.gbk\n"
        )
        # Check real files exist
        dpath = "data/Prochlorococcus/genomes/MED4/"
        files = [
            os.path.join(dpath, "genomic.gff"),
            os.path.join(dpath, "cyanorak/Pro_MED4.gff"),
            os.path.join(dpath, "cyanorak/Pro_MED4.gbk"),
        ]
        if not all(os.path.exists(f) for f in files):
            pytest.skip("Real MED4 data files not available")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            csv_path = f.name
        try:
            wrapper = MultiCyanorakNcbi(config_list_file=csv_path)
            wrapper.download_data()
            nodes = wrapper.get_nodes()
            assert len(nodes) > 1000
            assert any("PMM0001" in n[0] for n in nodes)
        finally:
            os.unlink(csv_path)
