"""--normalize leaves the disk-only raw IPS cache gzipped, whichever form it found."""
import gzip
import importlib.util
import json
from pathlib import Path

import pytest

RUNNER = Path(__file__).resolve().parents[1] / ".claude/skills/interproscan-run/run_interproscan.py"


@pytest.fixture(scope="module")
def runner():
    spec = importlib.util.spec_from_file_location("run_interproscan", RUNNER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _strain_dir(tmp_path: Path, raw: dict, gz: bool) -> Path:
    data_dir = tmp_path / "S1"
    out = data_dir / "interproscan"
    out.mkdir(parents=True)
    name = out / ("S1.interproscan.raw.json" + (".gz" if gz else ""))
    opener = gzip.open if gz else open
    with opener(name, "wt") as f:
        json.dump(raw, f)
    return data_dir


@pytest.mark.parametrize("gz", [False, True])
def test_normalize_leaves_raw_gzipped(runner, tmp_path, gz):
    data_dir = _strain_dir(tmp_path, {"interproscan-version": "5.78-109.0", "results": []}, gz)
    status, _ = runner.normalize_strain("S1", data_dir)
    out = data_dir / "interproscan"
    assert status == "NORMALIZED"
    assert (out / "S1.interproscan.raw.json.gz").exists()
    assert not (out / "S1.interproscan.raw.json").exists()
    with gzip.open(out / "S1.interproscan.raw.json.gz", "rt") as f:
        assert json.load(f)["results"] == []
