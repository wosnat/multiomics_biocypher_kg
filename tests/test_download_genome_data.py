"""Unit tests for scripts/download_genome_data.py.

Patching strategy
-----------------
The step functions use *lazy* imports (``from X import Y`` inside the
function body) to avoid loading heavy dependencies (pypath, adapters) when
the script is invoked with ``--help``.  Because the name ``Y`` is bound
inside the function rather than at module level, we cannot patch
``scripts.download_genome_data.Y``.  Instead we patch the attribute on the
*source* module (e.g. ``multiomics_kg.adapters.cyanorak_ncbi_adapter.CyanorakNcbi``).
Python's import system caches modules in ``sys.modules``, so patching the
source attribute is picked up by subsequent ``from … import`` statements.

For ``PROJECT_ROOT`` (a module-level ``Path`` constant used by step 3 to
build cache paths) we patch it directly as
``scripts.download_genome_data.PROJECT_ROOT``.

Steps 1 and 2 construct the genome's ``data_dir`` as
``str(PROJECT_ROOT / genome["data_dir"])``.  Passing an *absolute* path in
``genome["data_dir"]`` makes pathlib ignore ``PROJECT_ROOT`` (absolute path
wins in ``Path.__truediv__``), so no ``PROJECT_ROOT`` patch is needed for
those steps.
"""

import json
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure the project root is importable when running pytest from any cwd
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.download_genome_data import (
    _get_org_group,
    _read_genomes_csv,
    step1_ncbi,
    step2_cyanorak,
    step3_uniprot,
    step4_eggnog,
)


# ── fixtures ─────────────────────────────────────────────────────────────────────

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


# ── _read_genomes_csv ────────────────────────────────────────────────────────────

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


# ── _get_org_group ───────────────────────────────────────────────────────────────

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


# ── step1_ncbi ───────────────────────────────────────────────────────────────────

# Patch targets for lazy imports inside step1_ncbi / step2_cyanorak:
_CYAN_CLS   = "multiomics_kg.adapters.cyanorak_ncbi_adapter.CyanorakNcbi"
_CURL_OFF   = "pypath.share.curl.cache_off"   # only needed for force=True tests
# Patch target for lazy imports inside step3_uniprot:
_UNIPROT_CLS = "multiomics_kg.adapters.uniprot_adapter.Uniprot"
_CURL_OFF_UP = "pypath.share.curl.cache_off"  # step3 always calls cache_off()


class TestStep1Ncbi:
    """step1_ncbi: NCBI genome download.

    genome["data_dir"] is set to an *absolute* tmp path, so pathlib's
    ``PROJECT_ROOT / absolute_path`` resolves to ``absolute_path`` (absolute
    wins).  No PROJECT_ROOT patching needed for steps 1/2.
    """

    def test_skips_when_all_files_exist(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        for fname in ["genomic.gff", "protein.faa", "genomic.gbff"]:
            (data_dir / fname).write_text("content")

        genome = {"strain_name": "MED4", "ncbi_accession": "GCF_000011465.1",
                  "data_dir": str(data_dir)}

        with patch(_CYAN_CLS) as mock_cls:
            step1_ncbi([genome], force=False)

        mock_cls.assert_not_called()

    def test_calls_download_when_files_missing(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        (data_dir / "genomic.gff").write_text("content")  # only one of three

        genome = {"strain_name": "MED4", "ncbi_accession": "GCF_000011465.1",
                  "data_dir": str(data_dir)}
        mock_adapter = MagicMock()

        with patch(_CYAN_CLS, return_value=mock_adapter):
            step1_ncbi([genome], force=False)

        mock_adapter._download_ncbi_genome.assert_called_once()

    def test_force_removes_existing_files_and_redownloads(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        files = []
        for fname in ["genomic.gff", "protein.faa", "genomic.gbff"]:
            p = data_dir / fname
            p.write_text("old content")
            files.append(p)

        genome = {"strain_name": "MED4", "ncbi_accession": "GCF_000011465.1",
                  "data_dir": str(data_dir)}
        mock_adapter = MagicMock()

        with patch(_CYAN_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF, return_value=MagicMock()):
            step1_ncbi([genome], force=True)

        for p in files:
            assert not p.exists(), f"Expected {p.name} to be removed before re-download"
        mock_adapter._download_ncbi_genome.assert_called_once()

    def test_processes_multiple_genomes(self, tmp):
        mock_adapters = []
        genomes = []
        for strain in ["MED4", "MIT9313"]:
            data_dir = tmp / "genomes" / strain
            data_dir.mkdir(parents=True)
            genomes.append({"strain_name": strain, "ncbi_accession": f"GCF_{strain}",
                            "data_dir": str(data_dir)})
            mock_adapters.append(MagicMock())

        with patch(_CYAN_CLS, side_effect=mock_adapters):
            step1_ncbi(genomes, force=False)

        for adapter in mock_adapters:
            adapter._download_ncbi_genome.assert_called_once()

    def test_adapter_instantiated_with_correct_args(self, tmp):
        data_dir = tmp / "genomes" / "MED4"
        data_dir.mkdir(parents=True)
        genome = {"strain_name": "MED4", "ncbi_accession": "GCF_000011465.1",
                  "data_dir": str(data_dir)}
        mock_adapter = MagicMock()

        with patch(_CYAN_CLS, return_value=mock_adapter) as mock_cls:
            step1_ncbi([genome], force=False)

        mock_cls.assert_called_once_with(
            ncbi_accession="GCF_000011465.1",
            data_dir=str(data_dir),
        )


# ── step2_cyanorak ───────────────────────────────────────────────────────────────

class TestStep2Cyanorak:
    def test_skips_strain_without_cyanorak_organism(self, tmp):
        data_dir = tmp / "genomes" / "MIT1002"
        data_dir.mkdir(parents=True)
        genome = {"strain_name": "MIT1002", "cyanorak_organism": "",
                  "data_dir": str(data_dir)}

        with patch(_CYAN_CLS) as mock_cls:
            step2_cyanorak([genome], force=False)

        mock_cls.assert_not_called()

    def test_skips_when_both_files_exist(self, tmp):
        cyan_org = "Pro_MED4"
        data_dir = tmp / "genomes" / "MED4"
        cyan_dir = data_dir / "cyanorak"
        cyan_dir.mkdir(parents=True)
        (cyan_dir / f"{cyan_org}.gff").write_text("content")
        (cyan_dir / f"{cyan_org}.gbk").write_text("content")

        genome = {"strain_name": "MED4", "cyanorak_organism": cyan_org,
                  "data_dir": str(data_dir)}

        with patch(_CYAN_CLS) as mock_cls:
            step2_cyanorak([genome], force=False)

        mock_cls.assert_not_called()

    def test_downloads_when_gbk_missing(self, tmp):
        cyan_org = "Pro_MED4"
        data_dir = tmp / "genomes" / "MED4"
        cyan_dir = data_dir / "cyanorak"
        cyan_dir.mkdir(parents=True)
        (cyan_dir / f"{cyan_org}.gff").write_text("content")  # GFF present, GBK absent

        genome = {"strain_name": "MED4", "cyanorak_organism": cyan_org,
                  "data_dir": str(data_dir)}
        mock_adapter = MagicMock()

        with patch(_CYAN_CLS, return_value=mock_adapter):
            step2_cyanorak([genome], force=False)

        mock_adapter._download_cyanorak_gff.assert_called_once()
        mock_adapter._download_cyanorak_gbk.assert_called_once()

    def test_force_removes_files_and_redownloads(self, tmp):
        cyan_org = "Pro_MED4"
        data_dir = tmp / "genomes" / "MED4"
        cyan_dir = data_dir / "cyanorak"
        cyan_dir.mkdir(parents=True)
        gff = cyan_dir / f"{cyan_org}.gff"
        gbk = cyan_dir / f"{cyan_org}.gbk"
        gff.write_text("old")
        gbk.write_text("old")

        genome = {"strain_name": "MED4", "cyanorak_organism": cyan_org,
                  "data_dir": str(data_dir)}
        mock_adapter = MagicMock()

        with patch(_CYAN_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF, return_value=MagicMock()):
            step2_cyanorak([genome], force=True)

        assert not gff.exists()
        assert not gbk.exists()
        mock_adapter._download_cyanorak_gff.assert_called_once()
        mock_adapter._download_cyanorak_gbk.assert_called_once()


# ── step3_uniprot ────────────────────────────────────────────────────────────────

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

    def _mock_adapter(self, data=None):
        m = MagicMock()
        m.data = data or {"cc_function": {"P00001": "test"}}
        m.uniprot_ids = {"P00001"}
        return m

    def test_deduplicates_by_taxid(self, tmp):
        """Three genomes with two unique taxids → two UniProt downloads."""
        genomes = self._make_genomes()
        call_count = [0]

        def make_adapter(**kwargs):
            call_count[0] += 1
            return self._mock_adapter()

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, side_effect=make_adapter), \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=False)

        assert call_count[0] == 2  # one per unique (org_group, taxid)

    def test_skips_when_both_json_files_exist(self, tmp):
        genomes = [self._make_genomes()[0]]  # MED4 only
        cache_dir = tmp / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        cache_dir.mkdir(parents=True)
        (cache_dir / "uniprot_raw_data.json").write_text("{}")
        (cache_dir / "uniprot_preprocess_data.json").write_text("{}")

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS) as mock_cls:
            step3_uniprot(genomes, force=False)

        mock_cls.assert_not_called()

    def test_downloads_and_saves_json_files(self, tmp):
        genomes = [self._make_genomes()[0]]
        mock_adapter = self._mock_adapter()

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=False)

        cache_dir = tmp / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        assert (cache_dir / "uniprot_raw_data.json").exists()
        assert (cache_dir / "uniprot_preprocess_data.json").exists()
        raw = json.loads((cache_dir / "uniprot_raw_data.json").read_text())
        assert "cc_function" in raw

    def test_calls_preprocess_methods(self, tmp):
        genomes = [self._make_genomes()[0]]
        mock_adapter = self._mock_adapter(data={})

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=False)

        mock_adapter._download_uniprot_data.assert_called_once()
        mock_adapter._preprocess_uniprot_data.assert_called_once()
        mock_adapter._preprocess_organisms.assert_called_once()

    def test_force_overwrites_existing_cache(self, tmp):
        genomes = [self._make_genomes()[0]]
        cache_dir = tmp / "cache" / "data" / "Prochlorococcus" / "uniprot" / "59919"
        cache_dir.mkdir(parents=True)
        (cache_dir / "uniprot_raw_data.json").write_text('{"old": true}')
        (cache_dir / "uniprot_preprocess_data.json").write_text('{"old": true}')

        mock_adapter = self._mock_adapter(data={"cc_function": {"P00002": "new"}})
        mock_adapter.uniprot_ids = {"P00002"}

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=True)

        mock_adapter._download_uniprot_data.assert_called_once()
        raw = json.loads((cache_dir / "uniprot_raw_data.json").read_text())
        assert "cc_function" in raw

    def test_correct_org_group_path_for_alteromonas(self, tmp):
        """Alteromonas taxid 28108 → cache/data/Alteromonas/uniprot/28108/."""
        genomes = [self._make_genomes()[1]]  # MIT1002
        mock_adapter = self._mock_adapter(data={})

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, return_value=mock_adapter), \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=False)

        expected = tmp / "cache" / "data" / "Alteromonas" / "uniprot" / "28108"
        assert (expected / "uniprot_raw_data.json").exists()

    def test_rev_false_passed_to_adapter(self, tmp):
        """Verify Uniprot is instantiated with rev=False for full coverage."""
        genomes = [self._make_genomes()[0]]
        mock_adapter = self._mock_adapter(data={})

        with patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch(_UNIPROT_CLS, return_value=mock_adapter) as mock_cls, \
             patch(_CURL_OFF_UP, return_value=MagicMock()):
            step3_uniprot(genomes, force=False)

        _, kwargs = mock_cls.call_args
        assert kwargs.get("rev") is False


# ── step4_eggnog ─────────────────────────────────────────────────────────────────

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
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp):
            step4_eggnog([genome], force=False, cpu=1)
        # No subprocess should be called — just early return
        # (error is logged, not raised)

    def test_skips_when_eggnog_db_missing(self, tmp, capsys):
        genome, _ = self._genome(tmp)
        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(tmp / "nonexistent")}), \
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp):
            step4_eggnog([genome], force=False, cpu=1)

    def test_skips_strain_without_protein_faa(self, tmp):
        genome, data_dir = self._genome(tmp, create_faa=False)
        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
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
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
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
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
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
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", side_effect=fake_run) as mock_run:
            step4_eggnog([genome], force=True, cpu=1)

        mock_run.assert_called_once()

    def test_reports_failure_on_nonzero_exit(self, tmp, capsys):
        genome, data_dir = self._genome(tmp)
        eggnog_db = tmp / "eggnog_db"
        eggnog_db.mkdir()

        with patch.dict(os.environ, {"EGGNOG_DATA_DIR": str(eggnog_db)}), \
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
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
             patch("scripts.download_genome_data.PROJECT_ROOT", tmp), \
             patch("subprocess.run", side_effect=fake_run) as mock_run:
            step4_eggnog([genome], force=False, cpu=1)

        cmd = mock_run.call_args[0][0]
        idx = cmd.index("--data_dir")
        assert cmd[idx + 1] == str(eggnog_db)
