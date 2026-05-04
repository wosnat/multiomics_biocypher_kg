"""Tests for Phase 2 metabolomics extension of step 6 (build_kegg_metabolism_xrefs)."""
import json
import sqlite3
from pathlib import Path

import pandas as pd
import yaml


def _make_stub_resolver(tmp_path: Path) -> sqlite3.Connection:
    """Build a tiny SQLite resolver that mimics the MNX schema.

    Just enough for resolve_metabolite + mnxm_to_primary_id to work for our cases.
    """
    db_path = tmp_path / "stub.db"
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT, mass REAL, inchikey TEXT);
        CREATE TABLE compound_aliases (mnxm_id TEXT, source TEXT, value TEXT);
        CREATE TABLE compound_names (mnxm_id TEXT, name_normalized TEXT);
        INSERT INTO compounds VALUES ('MNXM41', 'D-Glucose', 'C6H12O6', 180.06, 'WQZGKKKJIJFFOK-GASJEMHNSA-N');
        INSERT INTO compounds VALUES ('MNXM58', 'Pyruvate', 'C3H3O3', 87.03, 'LCTONWCANYUPML-UHFFFAOYSA-M');
        INSERT INTO compound_aliases VALUES ('MNXM41', 'kegg.compound', 'C00031');
        INSERT INTO compound_aliases VALUES ('MNXM58', 'kegg.compound', 'C00022');
        INSERT INTO compound_names VALUES ('MNXM41', 'glucose');
        INSERT INTO compound_names VALUES ('MNXM58', 'pyruvate');
    """)
    conn.commit()
    return conn


def _make_paperconfig(
    tmp_path: Path, csv_path: Path, *, with_aliases: bool = False
) -> tuple[Path, dict]:
    """Create a minimal paperconfig with one metabolite_assays_table entry."""
    pc_path = tmp_path / "paperconfig.yaml"
    cfg = {
        "publication": {
            "papername": "Test 2026",
            "doi": "10.1/test",
            "experiments": {
                "exp1": {"organism": "Prochlorococcus MIT9303", "compartment": "whole_cell"},
            },
            "supplementary_materials": {
                "tab_a": {
                    "type": "metabolite_assays_table",
                    "filename": str(csv_path),
                    "experiment": "exp1",
                    "organism": "Prochlorococcus MIT9303",
                    "name_col": "compound",
                    "id_col": "KEGG_ID",
                    "id_type": "kegg.compound",
                    "assays": [{
                        "metric_type": "x",
                        "value_kind": "numeric",
                        "field_description": "d",
                        "sample_columns": [{"replicate_columns": ["c"]}],
                    }],
                }
            },
        }
    }
    if with_aliases:
        cfg["publication"]["supplementary_materials"]["tab_a"]["aliases_file"] = "metabolite_aliases.yaml"
        (tmp_path / "metabolite_aliases.yaml").write_text(
            'MysteryAcid: "kegg.compound:C00022"\n'
        )
    pc_path.write_text(yaml.safe_dump(cfg))
    return pc_path, cfg


def test_harvest_paper_metabolites_id_col_direct(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_metabolites

    csv = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["GlucoseRow"], "KEGG_ID": ["C00031"], "c": ["1"]}).to_csv(csv, index=False)
    pc_path, cfg = _make_paperconfig(tmp_path, csv)
    conn = _make_stub_resolver(tmp_path)

    harvest = _harvest_paper_metabolites([(pc_path, cfg)], conn)

    assert harvest["alias_to_primary"]["GlucoseRow"] == "kegg.compound:C00031"
    assert harvest["resolution_methods"]["GlucoseRow"] == "kegg_direct"
    assert "C00031" in harvest["paper_kegg_cpds"]
    assert harvest["paper_non_kegg_compounds"] == {}
    assert harvest["unresolved"] == []


def test_harvest_paper_metabolites_aliases_override_beats_resolver(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_metabolites

    csv = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["MysteryAcid"], "KEGG_ID": [""], "c": ["1"]}).to_csv(csv, index=False)
    pc_path, cfg = _make_paperconfig(tmp_path, csv, with_aliases=True)
    conn = _make_stub_resolver(tmp_path)

    harvest = _harvest_paper_metabolites([(pc_path, cfg)], conn)

    assert harvest["alias_to_primary"]["MysteryAcid"] == "kegg.compound:C00022"
    assert harvest["resolution_methods"]["MysteryAcid"] == "alias_override"


def test_harvest_paper_metabolites_name_match_via_resolver(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_metabolites

    csv = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["glucose"], "KEGG_ID": [""], "c": ["1"]}).to_csv(csv, index=False)
    pc_path, cfg = _make_paperconfig(tmp_path, csv)
    conn = _make_stub_resolver(tmp_path)

    harvest = _harvest_paper_metabolites([(pc_path, cfg)], conn)

    assert harvest["alias_to_primary"]["glucose"] == "kegg.compound:C00031"
    assert harvest["resolution_methods"]["glucose"] == "name_match"


def test_harvest_paper_metabolites_unresolved(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_metabolites

    csv = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["UnknownX"], "KEGG_ID": [""], "c": ["1"]}).to_csv(csv, index=False)
    pc_path, cfg = _make_paperconfig(tmp_path, csv)
    conn = _make_stub_resolver(tmp_path)

    harvest = _harvest_paper_metabolites([(pc_path, cfg)], conn)

    assert "UnknownX" in harvest["unresolved"]
    assert harvest["resolution_methods"]["UnknownX"] == "unresolved"
    assert "UnknownX" not in harvest["alias_to_primary"]


def test_fold_paper_metabolites_unions_evidence_sources():
    # kegg_data["compounds"] is keyed by BARE C-numbers in production (the
    # `kegg.compound:` prefix is added at adapter time when emitting Metabolite
    # nodes). The fold function must look up using bare keys. The earlier
    # version of this test used prefixed keys and silently passed despite the
    # production bug — see commit 1049a9f for the fix.
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _fold_paper_metabolites_into_kegg_data
    kegg_data = {
        "compounds": {
            "C00031": {
                "name": "Glucose",
                "evidence_sources": ["metabolism"],
            },
            "C00022": {
                "name": "Pyruvate",
                "evidence_sources": [],
            },
        },
        "additional_compounds": {},
    }
    paper_kegg_cpds = {"C00031", "C00022"}
    paper_non_kegg = {
        "chebi:16865": {"name": "GABA", "formula": "C4H9NO2", "mass": 103.06, "inchikey": None,
                        "mnxm_id": "MNXM61", "chebi_id": "16865"},
    }
    _fold_paper_metabolites_into_kegg_data(
        kegg_data,
        paper_kegg_cpds=paper_kegg_cpds,
        paper_non_kegg_compounds=paper_non_kegg,
    )
    assert sorted(kegg_data["compounds"]["C00031"]["evidence_sources"]) == ["metabolism", "metabolomics"]
    assert kegg_data["compounds"]["C00022"]["evidence_sources"] == ["metabolomics"]
    assert kegg_data["additional_compounds"]["chebi:16865"]["evidence_sources"] == ["metabolomics"]
    assert kegg_data["additional_compounds"]["chebi:16865"]["name"] == "GABA"


def test_write_metabolite_id_mapping_creates_three_tier_skeleton(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _write_metabolite_id_mapping

    harvest = {
        "alias_to_primary": {"Glucose": "kegg.compound:C00031", "GABA": "kegg.compound:C00334"},
        "resolution_methods": {"Glucose": "kegg_direct", "GABA": "alias_override"},
        "unresolved": ["MysteryX"],
        "per_paper": {"/tmp/p.yaml": {"resolved": 2, "unresolved": 1, "total": 3}},
    }
    out = tmp_path / "metabolite_id_mapping.json"
    _write_metabolite_id_mapping(harvest, out)
    data = json.loads(out.read_text())
    assert data["specific_lookup"] == {}
    assert data["multi_lookup"] == {}
    assert data["conflicts"] == {}
    assert data["name_lookup"]["Glucose"] == ["kegg.compound:C00031"]
    assert data["name_lookup"]["GABA"] == ["kegg.compound:C00334"]
    assert "compounds" in data
    assert "Glucose" in data["compounds"]["kegg.compound:C00031"]["aliases"]
    assert data["unresolved"] == ["MysteryX"]
    assert data["schema_version"] == 1
