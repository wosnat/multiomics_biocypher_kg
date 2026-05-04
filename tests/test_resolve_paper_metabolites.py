"""Tests for step 7 (resolve_paper_metabolites)."""
import json

import pandas as pd


def test_resolve_paper_metabolites_writes_resolved_csv(tmp_path):
    from multiomics_kg.download.resolve_paper_metabolites import resolve_paper_metabolites_for_entry

    csv_path = tmp_path / "metab.csv"
    pd.DataFrame({
        "compound": ["Glucose", "GABA", "MysteryX"],
        "KEGG_ID": ["C00031", "", ""],
    }).to_csv(csv_path, index=False)

    mapping = {
        "specific_lookup": {},
        "multi_lookup": {},
        "name_lookup": {
            "Glucose": ["kegg.compound:C00031"],
            "GABA": ["kegg.compound:C00334"],
        },
        "conflicts": {},
        "compounds": {},
    }

    entry = {
        "type": "metabolite_assays_table",
        "filename": str(csv_path),
        "name_col": "compound",
        "id_col": "KEGG_ID",
        "id_type": "kegg.compound",
    }

    out_csv, report = resolve_paper_metabolites_for_entry(entry, mapping)
    assert out_csv == csv_path.with_name("metab_resolved.csv")
    df = pd.read_csv(out_csv, dtype=str, keep_default_na=False)
    rows = {r["compound"]: (r["metabolite_id"], r["resolution_method"]) for _, r in df.iterrows()}
    assert rows["Glucose"] == ("kegg.compound:C00031", "kegg_direct")
    assert rows["GABA"] == ("kegg.compound:C00334", "name_match")
    assert rows["MysteryX"] == ("", "unresolved")
    assert report["resolved"] == 2
    assert report["unresolved"] == 1
    assert report["total_rows"] == 3
    assert report["method_counts"]["kegg_direct"] == 1
    assert report["method_counts"]["name_match"] == 1
    assert report["method_counts"]["unresolved"] == 1
    # Report file is also written
    report_path = csv_path.with_name("metab_resolution_report.json")
    assert report_path.exists()
    saved = json.loads(report_path.read_text())
    assert saved["resolved"] == 2


def test_resolve_paper_metabolites_ambiguous_picks_first_sorted(tmp_path):
    from multiomics_kg.download.resolve_paper_metabolites import resolve_paper_metabolites_for_entry

    csv_path = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["AmbName"], "KEGG_ID": [""]}).to_csv(csv_path, index=False)

    mapping = {
        "name_lookup": {
            "AmbName": ["kegg.compound:C00200", "kegg.compound:C00100"],
        },
    }

    entry = {
        "filename": str(csv_path),
        "name_col": "compound",
        "id_col": "KEGG_ID",
        "id_type": "kegg.compound",
    }
    _, report = resolve_paper_metabolites_for_entry(entry, mapping)
    df = pd.read_csv(csv_path.with_name("metab_resolved.csv"), dtype=str, keep_default_na=False)
    row = df.iloc[0].to_dict()
    assert row["metabolite_id"] == "kegg.compound:C00100"   # sorted first
    assert row["resolution_method"] == "ambiguous_multi_id"


def test_resolve_paper_metabolites_id_col_with_prefix(tmp_path):
    """User-provided id_col cell already containing a prefix (e.g. 'kegg.compound:C00031') should pass through."""
    from multiomics_kg.download.resolve_paper_metabolites import resolve_paper_metabolites_for_entry

    csv_path = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["X"], "KEGG_ID": ["kegg.compound:C00022"]}).to_csv(csv_path, index=False)
    mapping = {"name_lookup": {}}
    entry = {
        "filename": str(csv_path),
        "name_col": "compound",
        "id_col": "KEGG_ID",
        "id_type": "kegg.compound",
    }
    _, report = resolve_paper_metabolites_for_entry(entry, mapping)
    df = pd.read_csv(csv_path.with_name("metab_resolved.csv"), dtype=str, keep_default_na=False)
    assert df.iloc[0]["metabolite_id"] == "kegg.compound:C00022"
    assert df.iloc[0]["resolution_method"] == "kegg_direct"


def test_resolve_paper_metabolites_skip_rows(tmp_path):
    """skip_rows must be honored when reading the source CSV."""
    from multiomics_kg.download.resolve_paper_metabolites import resolve_paper_metabolites_for_entry

    csv_path = tmp_path / "metab.csv"
    csv_path.write_text("# comment\n# another comment\ncompound,KEGG_ID\nGlucose,C00031\n")

    mapping = {"name_lookup": {}}
    entry = {
        "filename": str(csv_path),
        "name_col": "compound",
        "id_col": "KEGG_ID",
        "id_type": "kegg.compound",
        "skip_rows": 2,
    }
    _, report = resolve_paper_metabolites_for_entry(entry, mapping)
    df = pd.read_csv(csv_path.with_name("metab_resolved.csv"), dtype=str, keep_default_na=False)
    assert df.iloc[0]["compound"] == "Glucose"
    assert df.iloc[0]["metabolite_id"] == "kegg.compound:C00031"
