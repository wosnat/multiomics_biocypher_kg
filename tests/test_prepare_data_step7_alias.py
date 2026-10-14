import subprocess
import sys


def test_step7_module_is_a_noop():
    r = subprocess.run([sys.executable, "-m", "multiomics_kg.download.resolve_paper_metabolites", "--force"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "folded into step 6" in r.stdout


def test_prepare_data_default_steps_drop_7():
    text = open("scripts/prepare_data.sh").read()
    assert 'STEPS="0 9 1 2 3 4 5 6 8 10"' in text
    assert 'STEPS="9 1 2 3 4 5 6 8 10"' in text
    assert "folded into step 6" in text
