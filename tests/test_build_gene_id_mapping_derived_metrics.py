"""Unit tests for extract_rows_from_derived_metrics_table in build_gene_id_mapping.

A derived_metrics_table entry uses id_columns identically to a csv entry,
but without the statistical_analyses indirection. The extractor harvests
(id_val, id_type) pairs per row for GeneIdGraph seeding.
"""
from pathlib import Path

import pandas as pd

from multiomics_kg.download.build_gene_id_mapping import (
    extract_rows_from_derived_metrics_table,
)


def test_extract_rows_basic(tmp_path, monkeypatch):
    """name_col + id_columns yields one row per non-empty CSV row."""
    from multiomics_kg.download import build_gene_id_mapping as bgm

    csv = tmp_path / "s4a.csv"
    csv.write_text(
        "NCBI ID_2,Gene Name,flag\n"
        "PMN2A_RS00015,dnaN,Y\n"
        "PMN2A_RS00020,,Y\n"
    )
    # Patch PROJECT_ROOT so the extractor's relative-path resolution lands in tmp_path.
    monkeypatch.setattr(bgm, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": csv.name,
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp1",
        "name_col": "NCBI ID_2",
        "id_columns": [
            {"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"},
            {"column": "Gene Name", "id_type": "gene_name"},
        ],
        "metrics": [{
            "metric_type": "periodic_in_axenic_LD",
            "value_kind": "boolean",
            "value_col": "flag",
            "true_tokens": ["Y"],
        }],
    }
    rows = extract_rows_from_derived_metrics_table(entry, "Biller 2018", "dm_s4a")
    # Row 1: (PMN2A_RS00015, locus_tag_ncbi), (dnaN, gene_name) — gene_name is second-order
    # Row 2: (PMN2A_RS00020, locus_tag_ncbi) — gene_name blank → not emitted
    assert len(rows) == 2
    pairs_row_1 = {p for p in rows[0][0]}
    assert ("PMN2A_RS00015", "locus_tag_ncbi") in pairs_row_1
    assert ("dnaN", "gene_name") in pairs_row_1
    pairs_row_2 = {p for p in rows[1][0]}
    assert ("PMN2A_RS00020", "locus_tag_ncbi") in pairs_row_2
    assert not any(pt == "gene_name" for _, pt in pairs_row_2)


def test_extract_rows_missing_id_columns_returns_empty(tmp_path, monkeypatch):
    from multiomics_kg.download import build_gene_id_mapping as bgm

    csv = tmp_path / "s4a.csv"
    csv.write_text("NCBI ID_2,flag\nPMN2A_RS00015,Y\n")
    monkeypatch.setattr(bgm, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": csv.name,
        "name_col": "NCBI ID_2",
        # no id_columns → extractor returns no useful pairs
    }
    rows = extract_rows_from_derived_metrics_table(entry, "X", "dm_x")
    assert rows == []
