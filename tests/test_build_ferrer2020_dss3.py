"""Ferrer 2020 builder: blank new tags / proteins that NCBI pairs with another old tag."""
import importlib.util
from pathlib import Path

import pandas as pd

BUILDER = Path(__file__).resolve().parents[1] / "scripts/build_modified_csv/build_ferrer2020_dss3_rnaseq.py"


def _builder():
    spec = importlib.util.spec_from_file_location("build_ferrer2020_dss3_rnaseq", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GM = pd.DataFrame({
    "locus_tag_ncbi": ["SPO_RS21615", "SPO_RS21620", "SPO_RS07235"],
    "protein_id": ["WP_011242129.1", "WP_011242130.1", "WP_999.1"],
    "old_locus_tags": ["SPOA0352", "SPOA0353", None],
})


def _rows(*rows):
    return pd.DataFrame(rows, columns=["locus_tag", "new_locus_tag", "protein_id"])


def test_dropped_orf_loses_the_next_genes_new_tag_and_protein():
    # SPOA0351 was dropped by RefSeq; the xlsx gives it SPOA0352's RS tag and protein.
    out = _builder().blank_contradicting_ids(
        _rows(["SPOA0351", "SPO_RS21615", "WP_011242129.1"]), GM)
    assert out.loc[0, "new_locus_tag"] == ""
    assert out.loc[0, "protein_id"] == ""


def test_shifted_row_with_known_old_tag_still_blanked():
    # SPOA0352 carrying SPOA0353's RS tag: the original one-way check.
    out = _builder().blank_contradicting_ids(
        _rows(["SPOA0352", "SPO_RS21620", "WP_011242130.1"]), GM)
    assert out.loc[0, "new_locus_tag"] == ""
    assert out.loc[0, "protein_id"] == ""


def test_consistent_rows_and_unpaired_new_tags_kept():
    out = _builder().blank_contradicting_ids(_rows(
        ["SPOA0352", "SPO_RS21615", "WP_011242129.1"],
        ["SPO1418", "SPO_RS07235", "WP_999.1"],   # RS tag with no NCBI old tag: no contradiction
    ), GM)
    assert out["new_locus_tag"].tolist() == ["SPO_RS21615", "SPO_RS07235"]
    assert out["protein_id"].tolist() == ["WP_011242129.1", "WP_999.1"]
