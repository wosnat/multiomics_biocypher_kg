"""Unit tests for the CyanoRak gbff -> (gff + gbk) converter.

The converter bridges a CyanoRak GenBank export (ProMIT1327_* locus tags, every
annotation in CDS /note qualifiers) onto canonical NCBI locus tags via a
diamond protein-match table, so the standard ``build_gene_mapping`` CyanoRak
merge can join by locus_tag. These tests pin the note->GFF-attribute mapping
(the part most likely to silently corrupt CK clusters / roles).
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation

from multiomics_kg.download.convert_cyanorak_gbff_to_gff import (
    build_gff_attributes,
    _strip_code,
    _enc,
    _code_list,
)


def _cds(notes, gene="None", product="hypothetical protein"):
    return SeqFeature(
        FeatureLocation(0, 99, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["ProMIT1327_00001"], "gene": [gene],
                    "product": [product], "note": notes,
                    "translation": ["MAAA"]},
    )


def test_id_and_cluster_number():
    feat = _cds([
        "cyanorak ORF Id:CK_Pro_MIT1327_00001",
        "cyanorak cluster number:CK_00057378",
    ])
    a = build_gff_attributes(feat)
    assert a["ID"] == "CK_Pro_MIT1327_00001"
    assert a["cluster_number"] == "CK_00057378"  # the headline CyanoRak value


def test_multi_role_is_comma_joined():
    """A single space-separated role note must become comma-separated, matching
    the real CyanoRak GFF (e.g. ``D.1.5,R.2``) so the config delimiter splits it."""
    feat = _cds([
        "cyanorak ORF Id:CK_Pro_MIT1327_00007",
        "Cyanorak role:D.1.5 R.2",
        "Cyanorak role def.:D.1.5 R.2=Phosphorus",
        "TIGR role:156",
        "TIGR role def.:156=Hypothetical proteins / Conserved",
    ])
    a = build_gff_attributes(feat)
    assert a["cyanorak_Role"] == "D.1.5,R.2"
    assert a["cyanorak_Role_description"] == "Phosphorus"
    assert a["tIGR_Role"] == "156"
    assert a["tIGR_Role_description"] == "Hypothetical proteins / Conserved"


def test_go_terms_get_prefix_and_join():
    feat = _cds([
        "cyanorak ORF Id:CK_Pro_MIT1327_00002",
        "GO terms Molecular function:0003824",
        "GO terms Biological process:0008152",
    ])
    a = build_gff_attributes(feat)
    assert a["Ontology_term"] == "GO:0003824,GO:0008152"


def test_eggnog_and_domains_flattened():
    feat = _cds([
        "cyanorak ORF Id:CK_Pro_MIT1327_00003",
        "COG:COG0499 NOG39254",
        "cyaNOG:cyaNOG07342",
        "pFam:PF13181 PF13374",
        "Interpro Domain:IPR001234",
    ])
    a = build_gff_attributes(feat)
    assert a["eggNOG"] == "COG0499,NOG39254,cyaNOG07342"
    assert a["protein_domains"] == "PF13181,PF13374,IPR001234"


def test_name_skips_none_sentinel_uses_gene():
    none_feat = _cds(["cyanorak ORF Id:CK_Pro_MIT1327_00001"], gene="None")
    assert "Name" not in build_gff_attributes(none_feat)
    named = _cds(["cyanorak ORF Id:CK_Pro_MIT1327_00001"], gene="yrrB_1")
    assert build_gff_attributes(named)["Name"] == "yrrB_1"


def test_strip_code():
    assert _strip_code("157=Unknown function / General") == "Unknown function / General"
    assert _strip_code("no equals sign") == "no equals sign"


def test_enc_encodes_structural_chars_keeps_spaces():
    # commas/semicolons/equals percent-encoded; spaces stay literal (CyanoRak style)
    assert _enc("DNA replication, recombination") == "DNA replication%2C recombination"
    assert _enc("a=b;c") == "a%3Db%3Bc"


def test_code_list_splits_multiple_notes_and_spaces():
    notes = ["pFam:PF1 PF2", "pFam:PF3"]
    assert _code_list(notes, "pFam:") == ["PF1", "PF2", "PF3"]
