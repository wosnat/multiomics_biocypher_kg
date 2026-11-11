"""Duplicate table rows that resolve to one gene in one analysis are logged.

The adapter builds one edge id per (publication, analysis, gene), so BioCypher keeps a
single row; the adapter and the paperconfig validator now say so instead of staying silent.
"""
import logging
import sys
from pathlib import Path

import pandas as pd

from multiomics_kg.adapters.omics_adapter import _log_duplicate_gene_rows

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_paperconfig import _duplicate_resolved_genes_warning  # noqa: E402


def test_no_duplicates_is_silent(caplog):
    rows = {"e1": [("ncbigene:A", "A", 1.0, "significant")],
            "e2": [("ncbigene:B", "B", -1.0, "not significant")]}
    with caplog.at_level(logging.WARNING):
        assert _log_duplicate_gene_rows(rows, "an1", "t.csv") == 0
    assert "more than one row" not in caplog.text


def test_duplicates_counted_and_disagreement_reported(caplog):
    rows = {
        "e1": [("ncbigene:A", "probe1", 1.0, "significant"),
               ("ncbigene:A", "probe2", 0.2, "not significant")],   # disagree
        "e2": [("ncbigene:B", "B", 2.0, "significant"),
               ("ncbigene:B", "B", 2.0, "significant")],             # identical repeat
        "e3": [("ncbigene:C", "C", 1.0, "significant")],
    }
    with caplog.at_level(logging.WARNING):
        assert _log_duplicate_gene_rows(rows, "an1", "t.csv") == 2
    assert "an1: 2 genes have more than one row in t.csv (1 with differing" in caplog.text
    assert "probe1" in caplog.text and "probe2" in caplog.text


def _resolved(tmp_path, tags, fc=None):
    src = tmp_path / "table.csv"
    src.write_text("x\n")
    cols = {"resolved_locus_tag": tags}
    if fc is not None:
        cols["logFC"] = fc
    pd.DataFrame(cols).to_csv(tmp_path / "table_resolved.csv", index=False)
    return str(src)


def test_validator_warns_on_duplicate_resolved_tags(tmp_path):
    fn = _resolved(tmp_path, ["PMM0001", "PMM0002", "PMM0001", None, None])
    msg = _duplicate_resolved_genes_warning(fn, "supp_1")
    assert msg is not None and msg.startswith("supp_1: 1 genes") and "PMM0001" in msg


def test_validator_silent_without_duplicates_or_resolved_file(tmp_path):
    assert _duplicate_resolved_genes_warning(_resolved(tmp_path, ["A", "B", None, None]), "s") is None
    assert _duplicate_resolved_genes_warning(str(tmp_path / "missing.csv"), "s") is None


def test_validator_ignores_rows_without_a_fold_change(tmp_path):
    # A footnote row resolved to a gene (Beliaev 2014) has no logFC: the adapter skips it.
    fn = _resolved(tmp_path, ["A", "B", "A"], fc=["1.5", "-2", None])
    assert _duplicate_resolved_genes_warning(fn, "s", ["logFC"]) is None
    fn = _resolved(tmp_path, ["A", "B", "A"], fc=["1.5", "-2", "0.3*"])
    assert "1 genes" in _duplicate_resolved_genes_warning(fn, "s", ["logFC"])
