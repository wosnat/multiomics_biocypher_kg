"""End-to-end test: YAML transform pipeline produces resolved fields."""
from __future__ import annotations

import json
import sqlite3
import textwrap
from pathlib import Path

import pytest


EGGNOG_LINES = textwrap.dedent("""\
    ##
    #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
    WP_001.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\tR00299\t-\t-\t3.A.1.1.1,99.X.99\tGH13_1,XX99\t-\t-
""")


@pytest.fixture(autouse=True)
def patch_metabolism_caches(monkeypatch, tmp_path):
    """Tiny resolver/tcdb/cazy fixtures matching the eggNOG test data above."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY(source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY(name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY(source, value, mnxr_id));
    """)
    conn.execute("INSERT INTO reactions VALUES ('MNXR101234', '', '', '', 'B', NULL)")
    conn.execute("INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00299', 'MNXR101234')")
    conn.commit()
    conn.close()

    from multiomics_kg.utils import metabolite_utils as mu
    monkeypatch.setattr(mu, "DEFAULT_DB_PATH", db)

    from multiomics_kg.download.utils import annotation_transforms as at
    monkeypatch.setattr(at, "_RESOLVER_CONN", None)

    # TCDB: 3.A.1.1.1 valid, 99.X.99 invalid
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    # CAZy: GH13_1 valid (and GH13 / GH for closure), XX99 invalid
    from multiomics_kg.utils import cazy_utils as cu
    cazy_path = tmp_path / "cazy_hierarchy.json"
    cazy_path.write_text(json.dumps({"GH": {}, "GH13": {}, "GH13_1": {}}))
    monkeypatch.setattr(cu, "DEFAULT_PATH", cazy_path)
    monkeypatch.setattr(cu, "_CACHE", None)

    yield


def test_yaml_transforms_produce_resolved_fields(tmp_path, monkeypatch):
    """Run the YAML pipeline against a tiny eggNOG file and assert merged fields.

    - kegg_reactions: R00299 → MNXR101234
    - transporter_classification: 3.A.1.1.1 valid, 99.X.99 dropped
    - cazy_ids: GH13_1 valid, XX99 dropped
    - No literal 'None' strings (regression test for the framework filter fix)
    """
    import yaml
    from multiomics_kg.download import build_gene_annotations as bga

    data_dir = tmp_path / "MED4_data"
    eggnog_dir = data_dir / "eggnog"
    eggnog_dir.mkdir(parents=True)
    (eggnog_dir / "MED4.emapper.annotations").write_text(EGGNOG_LINES)
    (data_dir / "gene_mapping.csv").write_text("locus_tag,protein_id\nPMM0001,WP_001.1\n")

    project_root = Path(__file__).parent.parent
    config = yaml.safe_load((project_root / "config/gene_annotations_config.yaml").read_text())

    row = {
        "strain_name":     "MED4",
        "preferred_name":  "Prochlorococcus MED4",
        "data_dir":        str(data_dir),
        "ncbi_taxon_id":   "59919",
    }
    bga.process_strain(row, config, force=True, pfam_data=None)

    merged = json.loads((data_dir / "gene_annotations_merged.json").read_text())
    gene = merged["PMM0001"]

    assert gene["kegg_reactions"] == ["MNXR101234"]
    assert gene["transporter_classification"] == ["3.A.1.1.1"]
    assert gene["cazy_ids"] == ["GH13_1"]
    # Regression: None must not leak as the literal string "None"
    assert "None" not in gene["transporter_classification"]
    assert "None" not in gene["cazy_ids"]
    assert "None" not in gene["kegg_reactions"]


def test_apply_transform_filters_none_from_list_path(monkeypatch):
    """Regression: _apply_transform's list path must drop None-returning transforms.

    Without this fix, a None-returning transform applied to a list value via
    a passthrough/single resolver would leak `None` into the merged dict, then
    later get serialized as the literal string 'None' downstream.
    """
    from multiomics_kg.download.build_gene_annotations import AnnotationBuilder

    # Minimal stub builder; _apply_transform doesn't depend on instance state
    builder = AnnotationBuilder.__new__(AnnotationBuilder)

    # validate_cazy is registered in _TRANSFORMS and returns None for invalid IDs.
    # The patch_metabolism_caches fixture (autouse) has already loaded a CAZy hierarchy
    # with GH13_1 valid and XX99 invalid.
    result = builder._apply_transform("validate_cazy", ["GH13_1", "XX99"])

    assert result == ["GH13_1"]
    assert None not in result
    assert "None" not in result  # belt-and-suspenders against str() coercion


def test_per_strain_metabolism_report_written(tmp_path, monkeypatch):
    """After process_strain runs, step2_metabolism_report.json sits next to gene_annotations_merged.json."""
    import yaml
    from multiomics_kg.download import build_gene_annotations as bga

    data_dir = tmp_path / "MED4_data"
    eggnog_dir = data_dir / "eggnog"
    eggnog_dir.mkdir(parents=True)
    (eggnog_dir / "MED4.emapper.annotations").write_text(EGGNOG_LINES)
    (data_dir / "gene_mapping.csv").write_text("locus_tag,protein_id\nPMM0001,WP_001.1\n")

    project_root = Path(__file__).parent.parent
    config = yaml.safe_load((project_root / "config/gene_annotations_config.yaml").read_text())

    row = {
        "strain_name":     "MED4",
        "preferred_name":  "Prochlorococcus MED4",
        "data_dir":        str(data_dir),
        "ncbi_taxon_id":   "59919",
    }
    bga.process_strain(row, config, force=True, pfam_data=None)

    report_path = data_dir / "step2_metabolism_report.json"
    assert report_path.exists()
    report = json.loads(report_path.read_text())

    assert report["strain"] == "MED4"
    assert report["gene_count"] == 1

    kr = report["kegg_reactions"]
    assert kr["raw_total"] == 1
    assert kr["resolved_total"] == 1
    assert kr["resolved_unique_mnxr"] == 1
    assert kr["unresolved_unique"] == 0

    tc = report["transporter_classification"]
    assert tc["raw_total"] == 2
    assert tc["validated_total"] == 1
    assert "99.X.99" in tc["invalid_examples"]

    cz = report["cazy_ids"]
    assert cz["raw_total"] == 2
    assert cz["validated_total"] == 1
    assert "XX99" in cz["invalid_examples"]
