"""Unit tests for multiomics_kg/download/utils/cli.py.

Coverage
--------
- add_common_args: argument parsing smoke test
- load_config: valid YAML, missing file (sys.exit)
- load_genome_rows: all rows, strain filter, non-matching filter (sys.exit)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.utils.cli import (
    add_common_args,
    load_config,
    load_genome_rows,
)


# ─── add_common_args ──────────────────────────────────────────────────────────

class TestAddCommonArgs:
    def _make_parser(self, default_config: Path) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        add_common_args(parser, default_config)
        return parser

    def test_defaults(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        parser = self._make_parser(cfg)
        args = parser.parse_args([])
        assert args.strains is None
        assert args.force is False
        assert args.config == str(cfg)

    def test_strains_parsed_as_list(self, tmp_path):
        parser = self._make_parser(tmp_path / "config.yaml")
        args = parser.parse_args(["--strains", "MED4", "MIT9313"])
        assert args.strains == ["MED4", "MIT9313"]

    def test_force_flag(self, tmp_path):
        parser = self._make_parser(tmp_path / "config.yaml")
        args = parser.parse_args(["--force"])
        assert args.force is True

    def test_config_override(self, tmp_path):
        parser = self._make_parser(tmp_path / "default.yaml")
        custom = tmp_path / "custom.yaml"
        args = parser.parse_args(["--config", str(custom)])
        assert args.config == str(custom)


# ─── load_config ──────────────────────────────────────────────────────────────

class TestLoadConfig:
    def test_loads_valid_yaml(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        cfg.write_text("key: value\nnumber: 42\n")
        result = load_config(cfg)
        assert result == {"key": "value", "number": 42}

    def test_returns_dict(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        cfg.write_text("a: 1\nb: 2\n")
        assert isinstance(load_config(cfg), dict)

    def test_missing_file_exits(self, tmp_path):
        missing = tmp_path / "nonexistent.yaml"
        with pytest.raises(SystemExit):
            load_config(missing)

    def test_missing_file_prints_error(self, tmp_path, capsys):
        missing = tmp_path / "nonexistent.yaml"
        with pytest.raises(SystemExit):
            load_config(missing)
        captured = capsys.readouterr()
        assert "nonexistent.yaml" in captured.err


# ─── load_genome_rows ─────────────────────────────────────────────────────────

GENOME_CSV_HEADER = "strain_name,data_dir,ncbi_accession,taxid\n"
GENOME_CSV_ROWS = [
    "MED4,data/Prochlorococcus/genomes/MED4,GCF_000011465.1,59919\n",
    "MIT9313,data/Prochlorococcus/genomes/MIT9313,GCF_000011485.1,74546\n",
    "AS9601,data/Prochlorococcus/genomes/AS9601,GCF_000015645.1,385988\n",
]


def _write_csv(path: Path, rows: list[str]) -> None:
    path.write_text(GENOME_CSV_HEADER + "".join(rows))


class TestLoadGenomeRows:
    def test_returns_all_rows(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, GENOME_CSV_ROWS)
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows()
        assert len(rows) == 3
        assert rows[0]["strain_name"] == "MED4"

    def test_strain_filter(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, GENOME_CSV_ROWS)
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows(strains=["MED4", "MIT9313"])
        assert [r["strain_name"] for r in rows] == ["MED4", "MIT9313"]

    def test_single_strain_filter(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, GENOME_CSV_ROWS)
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows(strains=["AS9601"])
        assert len(rows) == 1
        assert rows[0]["strain_name"] == "AS9601"

    def test_nonmatching_strain_exits(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, GENOME_CSV_ROWS)
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            with pytest.raises(SystemExit):
                load_genome_rows(strains=["NOSUCHSTRAIN"])

    def test_nonmatching_strain_prints_error(self, tmp_path, capsys):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, GENOME_CSV_ROWS)
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            with pytest.raises(SystemExit):
                load_genome_rows(strains=["NOSUCHSTRAIN"])
        assert "NOSUCHSTRAIN" in capsys.readouterr().err

    def test_skips_comment_lines(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            GENOME_CSV_HEADER
            + "# this is a comment\n"
            + GENOME_CSV_ROWS[0]
        )
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows()
        assert len(rows) == 1

    def test_skips_rows_missing_strain_name(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            GENOME_CSV_HEADER
            + ",data/Prochlorococcus/genomes/X,GCF_X,12345\n"  # no strain_name
            + GENOME_CSV_ROWS[0]
        )
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows()
        assert len(rows) == 1
        assert rows[0]["strain_name"] == "MED4"

    def test_skips_rows_missing_data_dir(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            GENOME_CSV_HEADER
            + "NODIRSTR,,GCF_X,12345\n"  # no data_dir
            + GENOME_CSV_ROWS[0]
        )
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows()
        assert len(rows) == 1
        assert rows[0]["strain_name"] == "MED4"

    def test_rows_contain_expected_keys(self, tmp_path):
        csv_path = tmp_path / "genomes.csv"
        _write_csv(csv_path, [GENOME_CSV_ROWS[0]])
        with patch("multiomics_kg.download.utils.cli.GENOMES_CSV", csv_path):
            rows = load_genome_rows()
        assert "strain_name" in rows[0]
        assert "data_dir" in rows[0]
        assert "ncbi_accession" in rows[0]
