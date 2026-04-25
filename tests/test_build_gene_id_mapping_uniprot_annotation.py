"""Tests for uniprot_annotation_string id_type expansion in step-3 row extraction.

The extractor lives in gene_id_utils; this file verifies that
build_gene_id_mapping calls it at the right place for each of the four
extraction functions (id_translation, csv, gene_clusters, derived_metrics).
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from multiomics_kg.download.build_gene_id_mapping import (
    extract_rows_from_csv_table,
    extract_rows_from_derived_metrics_table,
    extract_rows_from_gene_clusters_table,
    extract_rows_from_id_translation,
)


def _write_csv(path: Path, rows: list[dict], header: list[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)


@pytest.fixture
def annot_string():
    return ("Q31DF2_PROM9 Type II secretion system protein-like protein "
            "OS=Prochlorococcus marinus (strain MIT 9312) GN=PMT9312_0032 PE=4 SV=1")


def _flatten_pairs(rows):
    """Helper: flatten (row_pairs, source_label) tuples into a single list of pairs."""
    return [pair for row_pairs, _ in rows for pair in row_pairs]


def test_id_translation_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "id_translation.csv"
    _write_csv(csv_path, [{"Locus": "PMT9312_0032", "Annot": annot_string}],
               ["Locus", "Annot"])
    entry = {
        "type": "id_translation",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "id_columns": [
            {"column": "Locus", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    # build_gene_id_mapping reads paths relative to PROJECT_ROOT; patch via cwd
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_id_translation(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs
    # The raw long string must NOT be emitted with id_type=uniprot_annotation_string
    assert not any(t == "uniprot_annotation_string" for _, t in pairs)


def test_id_translation_no_match_emits_no_extra_pairs(tmp_path):
    csv_path = tmp_path / "id_translation.csv"
    _write_csv(csv_path, [{"Locus": "PMT9312_0032", "Annot": "hypothetical protein"}],
               ["Locus", "Annot"])
    entry = {
        "type": "id_translation",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "id_columns": [
            {"column": "Locus", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_id_translation(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    # Locus pair still emitted; Annot column contributes nothing
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert all(t in {"locus_tag"} for _, t in pairs)


def test_derived_metrics_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "dm.csv"
    _write_csv(csv_path,
               [{"Protein ID": "Q31DF2", "UniProt Annotation": annot_string,
                 "value_col": "1.5"}],
               ["Protein ID", "UniProt Annotation", "value_col"])
    entry = {
        "type": "derived_metrics_table",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "name_col": "Protein ID",
        "id_columns": [
            {"column": "Protein ID", "id_type": "uniprot_accession"},
            {"column": "UniProt Annotation", "id_type": "uniprot_annotation_string"},
        ],
        "metrics": [],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_derived_metrics_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("Q31DF2", "uniprot_accession") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs


def test_csv_table_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "csv_table.csv"
    _write_csv(csv_path,
               [{"Protein ID": "Q31DF2", "UniProt Annotation": annot_string,
                 "logfc": "1.5"}],
               ["Protein ID", "UniProt Annotation", "logfc"])
    entry = {
        "type": "csv",
        "filename": str(csv_path.relative_to(tmp_path)),
        "id_columns": [
            {"column": "Protein ID", "id_type": "uniprot_accession"},
            {"column": "UniProt Annotation", "id_type": "uniprot_annotation_string"},
        ],
        "statistical_analyses": [
            {"id": "test_de", "experiment": "exp", "name_col": "Protein ID",
             "logfc_col": "logfc"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_csv_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs


def test_gene_clusters_expands_uniprot_annotation(tmp_path, annot_string):
    csv_path = tmp_path / "clusters.csv"
    _write_csv(csv_path,
               [{"gene_id": "PMT9312_0032", "Annot": annot_string, "cluster": "C1"}],
               ["gene_id", "Annot", "cluster"])
    entry = {
        "type": "gene_clusters",
        "filename": str(csv_path.relative_to(tmp_path)),
        "organism": "Prochlorococcus MIT9312",
        "gene_id_col": "gene_id",
        "cluster_col": "cluster",
        "id_columns": [
            {"column": "gene_id", "id_type": "locus_tag"},
            {"column": "Annot", "id_type": "uniprot_annotation_string"},
        ],
    }
    import multiomics_kg.download.build_gene_id_mapping as bgim
    monkey_root = bgim.PROJECT_ROOT
    bgim.PROJECT_ROOT = tmp_path
    try:
        rows = extract_rows_from_gene_clusters_table(entry, "TestPaper", "tk")
    finally:
        bgim.PROJECT_ROOT = monkey_root
    pairs = _flatten_pairs(rows)
    assert ("PMT9312_0032", "locus_tag") in pairs
    assert ("Q31DF2_PROM9", "uniprot_entry_name") in pairs
    assert ("PMT9312_0032", "gene_name") in pairs
