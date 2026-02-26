"""Unit tests for multiomics_kg/download/download_genome_data.py.

Patching strategy
-----------------
The standalone download functions (_ncbi_download_genome, _cyanorak_download_file)
use *lazy* imports (``from X import Y`` inside the function body) to avoid loading
heavy dependencies at import time.

For pypath curl calls, patch the attribute on the *source* module:
  ``pypath.share.curl.Curl``
  ``pypath.share.curl.cache_off``

The UniProt step (step3_uniprot) delegates to
``multiomics_kg.download.download_uniprot.download_uniprot`` via a lazy import.
Patch it on the *source* module:
  ``multiomics_kg.download.download_uniprot.download_uniprot``

For step-level tests, patch the standalone functions directly:
  ``multiomics_kg.download.download_genome_data._ncbi_download_genome``
  ``multiomics_kg.download.download_genome_data._cyanorak_download_file``
  ``multiomics_kg.download.download_uniprot.download_uniprot``

For ``PROJECT_ROOT`` (a module-level ``Path`` constant used by step 3 to build
cache paths) patch it as ``multiomics_kg.download.download_genome_data.PROJECT_ROOT``.
Steps 1 and 2 construct ``data_dir`` as ``str(PROJECT_ROOT / genome["data_dir"])``;
passing an *absolute* path in ``genome["data_dir"]`` makes pathlib ignore
``PROJECT_ROOT``, so no patching is needed for those steps.
"""

import json
import logging
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure the project root is importable when running pytest from any cwd
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiomics_kg.download.download_genome_data import (
    _cyanorak_download_file,
    _get_org_group,
    _ncbi_download_genome,
    _read_genomes_csv,
    step1_ncbi,
    step2_cyanorak,
    step3_uniprot,
    step4_eggnog,
)


# ── fixtures ──────────────────────────────────────────────────────────────────

GENOMES_CSV_CONTENT = """\
ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade
GCF_000011465.1,Pro_MED4,59919,MED4,cache/data/Prochlorococcus/genomes/MED4/,HLI
# comment line should be ignored
GCF_000011485.1,Pro_MIT9313,74547,MIT9313,cache/data/Prochlorococcus/genomes/MIT9313/,LLIV
GCF_000014585.1,,28108,MIT1002,cache/data/Alteromonas/genomes/MIT1002/,
GCF_000015000.1,,28108,EZ55,cache/data/Alteromonas/genomes/EZ55/,
"""


@pytest.fixture
def tmp(tmp_path):
    """Return a temp directory as Path."""
    return tmp_path


@pytest.fixture
def genomes_csv(tmp):
    p = tmp / "cyanobacteria_genomes.csv"
    p.write_text(GENOMES_CSV_CONTENT)
    return p


@pytest.fixture
def genomes(genomes_csv):
    return _read_genomes_csv(genomes_csv)


# ── patch targets ─────────────────────────────────────────────────────────────

_CURL_CLS  = "pypath.share.curl.Curl"
_CURL_OFF  = "pypath.share.curl.cache_off"
_REQUESTS_GET = "requests.get"

_NCBI_DL     = "multiomics_kg.download.download_genome_data._ncbi_download_genome"
_CYANORAK_DL = "multiomics_kg.download.download_genome_data._cyanorak_download_file"
_UNIPROT_DL  = "multiomics_kg.download.download_uniprot.download_uniprot"


# ── _read_genomes_csv ─────────────────────────────────────────────────────────

class TestReadGenomesCsv:
    def test_parses_all_data_rows(self, genomes):
        assert len(genomes) == 4

    def test_skips_comment_lines(self, genomes):
        strains = [g["strain_name"] for g in genomes]
        assert "# comment line should be ignored" not in strains

    def test_row_keys(self, genomes):
        g = genomes[0]
        for key in ("ncbi_accession", "cyanorak_organism", "ncbi_taxon_id", "strain_name", "data_dir"):
            assert key in g

    def test_empty_cyanorak_organism(self, genomes):
        alteromonas = next(g for g in genomes if g["strain_name"] == "MIT1002")
        assert alteromonas["cyanorak_organism"] == ""

    def test_values_correct(self, genomes):
        med4 = genomes[0]
        assert med4["ncbi_accession"] == "GCF_000011465.1"
        assert med4["cyanorak_organism"] == "Pro_MED4"
        assert med4["ncbi_taxon_id"] == "59919"
        assert med4["strain_name"] == "MED4"


# ── _get_org_group ────────────────────────────────────────────────────────────

class TestGetOrgGroup:
    @pytest.mark.parametrize("data_dir,expected", [
        ("cache/data/Prochlorococcus/genomes/MED4/", "Prochlorococcus"),
        ("cache/data/Prochlorococcus/genomes/MED4",  "Prochlorococcus"),
        ("cache/data/Synechococcus/genomes/CC9311/",  "Synechococcus"),
        ("cache/data/Alteromonas/genomes/MIT1002/",   "Alteromonas"),
    ])
    def test_standard_paths(self, data_dir, expected):
        assert _get_org_group(data_dir) == expected

    def test_trailing_slash_stripped(self):
        assert _get_org_group("cache/data/Alteromonas/genomes/MIT1002/") == "Alteromonas"

    def test_no_trailing_slash(self):
        assert _get_org_group("cache/data/Alteromonas/genomes/MIT1002") == "Alteromonas"


# ── _ncbi_download_genome ─────────────────────────────────────────────────────

class TestNcbiDownloadGenome:
    def test_returns_false_when_all_files_exist(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        for fname in ["genomic.gff", "protein.faa", "genomic.gbff"]:
            (data_dir / fname).write_text("content")

        result = _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=False)
        assert result is False

    def test_returns_true_and_creates_files(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_c = MagicMock()
        mock_c.result = {
            "ncbi_dataset/data/GCF/genomic.gff":  "GFF content",
            "ncbi_dataset/data/GCF/protein.faa":  "FASTA content",
            "ncbi_dataset/data/GCF/genomic.gbff": "GBFF content",
        }

        with patch(_CURL_CLS, return_value=mock_c):
            result = _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=False)

        assert result is True
        assert (data_dir / "genomic.gff").read_text() == "GFF content"
        assert (data_dir / "protein.faa").read_text() == "FASTA content"
        assert (data_dir / "genomic.gbff").read_text() == "GBFF content"

    def test_partial_cache_triggers_download(self, tmp):
        """Only GFF exists — should still download."""
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        (data_dir / "genomic.gff").write_text("content")

        mock_c = MagicMock()
        mock_c.result = {
            "genomic.gff": "new gff",
            "protein.faa": "new faa",
            "genomic.gbff": "new gbff",
        }

        with patch(_CURL_CLS, return_value=mock_c):
            result = _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=False)

        assert result is True

    def test_force_removes_existing_files_and_redownloads(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        for fname in ["genomic.gff", "protein.faa", "genomic.gbff"]:
            (data_dir / fname).write_text("old")

        mock_c = MagicMock()
        mock_c.result = {
            "genomic.gff":  "new gff",
            "protein.faa":  "new faa",
            "genomic.gbff": "new gbff",
        }

        with patch(_CURL_CLS, return_value=mock_c), \
             patch(_CURL_OFF, return_value=MagicMock()):
            result = _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=True)

        assert result is True
        assert (data_dir / "genomic.gff").read_text() == "new gff"

    def test_raises_connection_error_on_failed_download(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_c = MagicMock()
        mock_c.result = None

        with patch(_CURL_CLS, return_value=mock_c):
            with pytest.raises(ConnectionError):
                _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=False)

    def test_creates_data_dir_if_missing(self, tmp):
        data_dir = tmp / "new" / "path" / "MED4"
        assert not data_dir.exists()

        mock_c = MagicMock()
        mock_c.result = {"genomic.gff": "g", "protein.faa": "p", "genomic.gbff": "b"}

        with patch(_CURL_CLS, return_value=mock_c):
            _ncbi_download_genome("GCF_000011465.1", str(data_dir), force=False)

        assert data_dir.exists()


# ── _cyanorak_download_file ───────────────────────────────────────────────────

class TestCyanorakDownloadFile:
    """Tests for _cyanorak_download_file.

    NOTE: The Cyanorak web server (bioinformatics.psb.ugent.be) is intermittently
    unavailable — it can return connection errors or throttle requests without warning.
    This is a known upstream issue. In production, use ``--skip-cyanorak`` to reuse
    cached files when the server is down.  These unit tests mock the HTTP layer so they
    are not affected by live server availability.
    """
    def test_returns_false_when_file_exists(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        cyan_dir = data_dir / "cyanorak"
        cyan_dir.mkdir(parents=True)
        (cyan_dir / "Pro_MED4.gff").write_text("content")

        result = _cyanorak_download_file("Pro_MED4", str(data_dir), "gff", force=False)
        assert result is False

    def test_returns_true_and_creates_gff(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = "GFF file content"

        with patch(_REQUESTS_GET, return_value=mock_resp):
            result = _cyanorak_download_file("Pro_MED4", str(data_dir), "gff", force=False)

        assert result is True
        assert (data_dir / "cyanorak" / "Pro_MED4.gff").read_text() == "GFF file content"

    def test_returns_true_and_creates_gbk(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = "GBK file content"

        with patch(_REQUESTS_GET, return_value=mock_resp):
            result = _cyanorak_download_file("Pro_MED4", str(data_dir), "gbk", force=False)

        assert result is True
        assert (data_dir / "cyanorak" / "Pro_MED4.gbk").read_text() == "GBK file content"

    def test_force_removes_existing_file_and_redownloads(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        cyan_dir = data_dir / "cyanorak"
        cyan_dir.mkdir(parents=True)
        gff = cyan_dir / "Pro_MED4.gff"
        gff.write_text("old")

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = "new content"

        with patch(_REQUESTS_GET, return_value=mock_resp):
            _cyanorak_download_file("Pro_MED4", str(data_dir), "gff", force=True)

        assert gff.read_text() == "new content"

    def test_raises_connection_error_on_failed_download(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_resp = MagicMock()
        mock_resp.ok = False
        mock_resp.status_code = 404

        with patch(_REQUESTS_GET, return_value=mock_resp):
            with pytest.raises(ConnectionError):
                _cyanorak_download_file("Pro_MED4", str(data_dir), "gff", force=False)

    def test_creates_cyanorak_subdir_if_missing(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)

        mock_c = MagicMock()
        mock_c.result = "content"

        with patch(_CURL_CLS, return_value=mock_c):
            _cyanorak_download_file("Pro_MED4", str(data_dir), "gff", force=False)

        assert (data_dir / "cyanorak").is_dir()


# ── step1_ncbi ────────────────────────────────────────────────────────────────

class TestStep1Ncbi:
    def _genome(self, tmp, strain="MED4", accession="GCF_000011465.1"):
        return {"strain_name": strain, "ncbi_accession": accession,
                "data_dir": str(tmp / "genomes" / strain)}

    def test_calls_download_for_each_genome(self, tmp):
        genomes = [
            self._genome(tmp, "MED4", "GCF_1"),
            self._genome(tmp, "MIT9313", "GCF_2"),
        ]
        with patch(_NCBI_DL, return_value=True) as mock_dl:
            step1_ncbi(genomes, force=False)
        assert mock_dl.call_count == 2

    def test_passes_correct_accession_and_data_dir(self, tmp):
        genome = self._genome(tmp)
        with patch(_NCBI_DL, return_value=True) as mock_dl:
            step1_ncbi([genome], force=False)
        accession, data_dir, force = mock_dl.call_args[0]
        assert accession == "GCF_000011465.1"
        assert data_dir == str(tmp / "genomes" / "MED4")
        assert force is False

    def test_passes_force_flag(self, tmp):
        genome = self._genome(tmp)
        with patch(_NCBI_DL, return_value=True) as mock_dl:
            step1_ncbi([genome], force=True)
        _, _, force = mock_dl.call_args[0]
        assert force is True

    def test_logs_skip_when_download_returns_false(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genome = self._genome(tmp)
        with patch(_NCBI_DL, return_value=False):
            step1_ncbi([genome], force=False)
        assert "SKIP" in caplog.text

    def test_logs_ok_when_download_returns_true(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genome = self._genome(tmp)
        with patch(_NCBI_DL, return_value=True):
            step1_ncbi([genome], force=False)
        assert "OK" in caplog.text


# ── step2_cyanorak ────────────────────────────────────────────────────────────

class TestStep2Cyanorak:
    def _genome(self, tmp, strain="MED4", cyan_org="Pro_MED4"):
        return {"strain_name": strain, "cyanorak_organism": cyan_org,
                "data_dir": str(tmp / "genomes" / strain)}

    def test_skips_strain_without_cyanorak_organism(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genome = self._genome(tmp, cyan_org="")
        with patch(_CYANORAK_DL) as mock_dl:
            step2_cyanorak([genome], force=False)
        mock_dl.assert_not_called()
        assert "SKIP" in caplog.text

    def test_calls_download_for_both_extensions(self, tmp):
        genome = self._genome(tmp)
        with patch(_CYANORAK_DL, return_value=True) as mock_dl:
            step2_cyanorak([genome], force=False)
        exts = [call[0][2] for call in mock_dl.call_args_list]
        assert sorted(exts) == ["gbk", "gff"]

    def test_passes_correct_organism_and_data_dir(self, tmp):
        genome = self._genome(tmp)
        with patch(_CYANORAK_DL, return_value=True) as mock_dl:
            step2_cyanorak([genome], force=False)
        first_call = mock_dl.call_args_list[0][0]
        assert first_call[0] == "Pro_MED4"
        assert first_call[1] == str(tmp / "genomes" / "MED4")

    def test_passes_force_flag(self, tmp):
        genome = self._genome(tmp)
        with patch(_CYANORAK_DL, return_value=True) as mock_dl:
            step2_cyanorak([genome], force=True)
        forces = [call[0][3] for call in mock_dl.call_args_list]
        assert all(f is True for f in forces)

    def test_logs_skip_when_both_return_false(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genome = self._genome(tmp)
        with patch(_CYANORAK_DL, return_value=False):
            step2_cyanorak([genome], force=False)
        assert "SKIP" in caplog.text

    def test_logs_ok_when_either_returns_true(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genome = self._genome(tmp)
        with patch(_CYANORAK_DL, side_effect=[True, False]):
            step2_cyanorak([genome], force=False)
        assert "OK" in caplog.text


# ── step3_uniprot ─────────────────────────────────────────────────────────────

class TestStep3Uniprot:
    def _make_genomes(self):
        return [
            {"strain_name": "MED4",    "ncbi_taxon_id": "59919",
             "data_dir": "cache/data/Prochlorococcus/genomes/MED4/"},
            {"strain_name": "MIT1002", "ncbi_taxon_id": "28108",
             "data_dir": "cache/data/Alteromonas/genomes/MIT1002/"},
            {"strain_name": "EZ55",    "ncbi_taxon_id": "28108",  # same taxid as MIT1002
             "data_dir": "cache/data/Alteromonas/genomes/EZ55/"},
        ]

    def test_deduplicates_by_taxid(self, tmp):
        """Three genomes with two unique taxids → two UniProt downloads."""
        genomes = self._make_genomes()
        with patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_DL, return_value=True) as mock_dl:
            step3_uniprot(genomes, force=False)
        assert mock_dl.call_count == 2

    def test_skips_when_download_returns_false(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genomes = [self._make_genomes()[0]]
        with patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_DL, return_value=False):
            step3_uniprot(genomes, force=False)
        assert "SKIP" in caplog.text

    def test_logs_ok_when_download_returns_true(self, tmp, caplog):
        caplog.set_level(logging.INFO)
        genomes = [self._make_genomes()[0]]
        with patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_DL, return_value=True):
            step3_uniprot(genomes, force=False)
        assert "OK" in caplog.text

    def test_passes_correct_taxid_and_cache_dir(self, tmp):
        genomes = [self._make_genomes()[0]]  # MED4, taxid=59919
        with patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_DL, return_value=True) as mock_dl:
            step3_uniprot(genomes, force=False)
        taxid, cache_dir, force = mock_dl.call_args[0]
        assert taxid == 59919
        assert cache_dir == tmp / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        assert force is False

    def test_correct_org_group_path_for_alteromonas(self, tmp):
        genomes = [self._make_genomes()[1]]  # MIT1002
        with patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_DL, return_value=True) as mock_dl:
            step3_uniprot(genomes, force=False)
        _, cache_dir, _ = mock_dl.call_args[0]
        assert cache_dir == tmp / "cache" / "data" / "Alteromonas" / "uniprot" / "28108"


# ── download_uniprot (standalone module) ──────────────────────────────────────

_FETCH_RAW   = "multiomics_kg.download.download_uniprot.fetch_raw_uniprot"
_PREPROCESS  = "multiomics_kg.download.download_uniprot.preprocess_uniprot_data"


class TestDownloadUniprot:
    """Tests for the standalone download_uniprot() function."""

    def _import(self):
        from multiomics_kg.download.download_uniprot import download_uniprot
        return download_uniprot

    def test_returns_false_when_both_files_exist(self, tmp):
        cache_dir = tmp / "uniprot" / "59919"
        cache_dir.mkdir(parents=True)
        (cache_dir / "uniprot_raw_data.json").write_text("{}")
        (cache_dir / "uniprot_preprocess_data.json").write_text("{}")

        download_uniprot = self._import()
        result = download_uniprot(59919, cache_dir, force=False)
        assert result is False

    def test_returns_true_and_saves_json_files(self, tmp):
        cache_dir = tmp / "uniprot" / "59919"
        raw_data = {"cc_function": {"P00001": "test"}}

        download_uniprot = self._import()
        with patch(_FETCH_RAW, return_value=(raw_data, {"P00001"})), \
             patch(_PREPROCESS, return_value=(raw_data, set())), \
             patch(_CURL_OFF, return_value=MagicMock()):
            result = download_uniprot(59919, cache_dir, force=False)

        assert result is True
        assert (cache_dir / "uniprot_raw_data.json").exists()
        assert (cache_dir / "uniprot_preprocess_data.json").exists()
        saved = json.loads((cache_dir / "uniprot_raw_data.json").read_text())
        assert "cc_function" in saved

    def test_calls_fetch_and_preprocess(self, tmp):
        cache_dir = tmp / "uniprot" / "59919"
        raw_data: dict = {}

        download_uniprot = self._import()
        with patch(_FETCH_RAW, return_value=(raw_data, set())) as mock_fetch, \
             patch(_PREPROCESS, return_value=(raw_data, set())) as mock_pre, \
             patch(_CURL_OFF, return_value=MagicMock()):
            download_uniprot(59919, cache_dir, force=False)

        mock_fetch.assert_called_once()
        mock_pre.assert_called_once()

    def test_force_overwrites_existing_cache(self, tmp):
        cache_dir = tmp / "uniprot" / "59919"
        cache_dir.mkdir(parents=True)
        (cache_dir / "uniprot_raw_data.json").write_text('{"old": true}')
        (cache_dir / "uniprot_preprocess_data.json").write_text('{"old": true}')

        new_data = {"cc_function": {"P00002": "new"}}
        download_uniprot = self._import()
        with patch(_FETCH_RAW, return_value=(new_data, {"P00002"})), \
             patch(_PREPROCESS, return_value=(new_data, set())), \
             patch(_CURL_OFF, return_value=MagicMock()):
            result = download_uniprot(59919, cache_dir, force=True)

        assert result is True
        saved = json.loads((cache_dir / "uniprot_raw_data.json").read_text())
        assert "cc_function" in saved

    def test_rev_false_passed_to_fetch(self, tmp):
        cache_dir = tmp / "uniprot" / "59919"
        raw_data: dict = {}

        download_uniprot = self._import()
        with patch(_FETCH_RAW, return_value=(raw_data, set())) as mock_fetch, \
             patch(_PREPROCESS, return_value=(raw_data, set())), \
             patch(_CURL_OFF, return_value=MagicMock()):
            download_uniprot(59919, cache_dir, force=False, rev=False)

        _, kwargs = mock_fetch.call_args
        assert kwargs.get("rev") is False

    def test_creates_cache_dir_if_missing(self, tmp):
        cache_dir = tmp / "new" / "path" / "59919"
        assert not cache_dir.exists()
        raw_data: dict = {}

        download_uniprot = self._import()
        with patch(_FETCH_RAW, return_value=(raw_data, set())), \
             patch(_PREPROCESS, return_value=(raw_data, set())), \
             patch(_CURL_OFF, return_value=MagicMock()):
            download_uniprot(59919, cache_dir, force=False)

        assert cache_dir.exists()


# ── step4_eggnog ──────────────────────────────────────────────────────────────

class TestStep4Eggnog:
    def _genome(self, tmp, strain="MED4", taxid="59919", create_faa=True):
        data_dir = tmp / "cache" / "data" / "Prochlorococcus" / "genomes" / strain
        data_dir.mkdir(parents=True)
        if create_faa:
            (data_dir / "protein.faa").write_text(">WP_001\nMACKEREL\n")
        return {
            "strain_name": strain,
            "ncbi_taxon_id": taxid,
            "data_dir": str(data_dir.relative_to(tmp)),
        }, data_dir

    def test_skips_when_eggnog_data_dir_not_set(self, tmp, capsys):
        genome, _ = self._genome(tmp)
        with patch.dict(os.environ, {}, clear=True), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp):
            step4_eggnog([genome], force=False, cpu=1)
        # No subprocess should be called — just early return
        # (error is logged, not raised)

    def test_skips_when_eggnog_db_missing(self, tmp, capsys):
        genome, _ = self._genome(tmp)
        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(tmp / "nonexistent")}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp):
            step4_eggnog([genome], force=False, cpu=1)

    def test_skips_strain_without_protein_faa(self, tmp):
        genome, data_dir = self._genome(tmp, create_faa=False)
        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run") as mock_run:
            step4_eggnog([genome], force=False, cpu=1)

        mock_run.assert_not_called()

    def test_skips_when_annotation_file_exists(self, tmp):
        genome, data_dir = self._genome(tmp)
        eggnog_dir = data_dir / "eggnog"
        eggnog_dir.mkdir()
        anno_file = eggnog_dir / "MED4.emapper.annotations"
        anno_file.write_text("# comment\nPMM0001\tresult\n")

        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run") as mock_run:
            step4_eggnog([genome], force=False, cpu=1)

        mock_run.assert_not_called()

    def test_runs_emapper_when_annotation_missing(self, tmp):
        genome, data_dir = self._genome(tmp)
        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()
        anno_file = data_dir / "eggnog" / "MED4.emapper.annotations"

        def fake_run(cmd, **kwargs):
            # Create annotation file to simulate successful run
            anno_file.parent.mkdir(parents=True, exist_ok=True)
            anno_file.write_text("# comment\nPMM0001\tresult\n")
            return MagicMock(returncode=0)

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", side_effect=fake_run) as mock_run:
            step4_eggnog([genome], force=False, cpu=2)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "-i" in cmd
        assert "--cpu" in cmd
        assert "2" in cmd
        assert "--data_dir" in cmd
        assert str(eggnog_db) in cmd

    def test_force_reruns_even_if_annotation_exists(self, tmp):
        genome, data_dir = self._genome(tmp)
        eggnog_dir = data_dir / "eggnog"
        eggnog_dir.mkdir()
        anno_file = eggnog_dir / "MED4.emapper.annotations"
        anno_file.write_text("# old\nPMM0001\told\n")

        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        def fake_run(cmd, **kwargs):
            anno_file.write_text("# new\nPMM0001\tnew\n")
            return MagicMock(returncode=0)

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", side_effect=fake_run) as mock_run:
            step4_eggnog([genome], force=True, cpu=1)

        mock_run.assert_called_once()

    def test_reports_failure_on_nonzero_exit(self, tmp, capsys):
        genome, data_dir = self._genome(tmp)
        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", return_value=MagicMock(returncode=1)):
            step4_eggnog([genome], force=False, cpu=1)

        out = capsys.readouterr().out
        assert "FAILED" in out

    def test_passes_eggnog_db_to_subprocess(self, tmp):
        genome, data_dir = self._genome(tmp)
        eggnog_db = tmp / "my_eggnog_db"
        eggnog_db.mkdir()
        anno_file = data_dir / "eggnog" / "MED4.emapper.annotations"

        def fake_run(cmd, **kwargs):
            anno_file.parent.mkdir(parents=True, exist_ok=True)
            anno_file.write_text("# comment\nPMM0001\tresult\n")
            return MagicMock(returncode=0)

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("multiomics_kg.download.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", side_effect=fake_run) as mock_run:
            step4_eggnog([genome], force=False, cpu=1)

        cmd = mock_run.call_args[0][0]
        idx = cmd.index("--data_dir")
        assert cmd[idx + 1] == str(eggnog_db)
