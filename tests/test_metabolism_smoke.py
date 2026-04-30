"""Integration smoke test: run prepare_data step 6 + adapter against the real cache.

Slow (~2 min). Skipped unless invoked with -m slow or "slow or kg".
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parent.parent
KEGG_DATA = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_data.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"
XREFS_OUT = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_metabolism_xrefs.json"


@pytest.mark.slow
def test_step6_smoke_run_against_real_cache():
    if not KEGG_DATA.exists():
        pytest.skip(f"{KEGG_DATA} missing — run KEGG download first")
    if not RESOLVER_DB.exists():
        pytest.skip(f"{RESOLVER_DB} missing — run prepare_data step 0 sub-step 7 first")
    # Verify kegg_data.json has reaction data (not just BRITE KO/pathway data)
    import json as _json
    kegg_data = _json.loads(KEGG_DATA.read_text())
    if "reaction_to_compounds" not in kegg_data:
        pytest.skip(
            "kegg_data.json lacks reaction_to_compounds — re-run KEGG download "
            "(uv run python -m multiomics_kg.download.download_genome_data --steps 6 --force)"
        )

    result = subprocess.run(
        [sys.executable, "-m", "multiomics_kg.download.build_kegg_metabolism_xrefs", "--force"],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    assert result.returncode == 0, f"step 6 failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    assert XREFS_OUT.exists()
    data = json.loads(XREFS_OUT.read_text())
    assert 2000 <= len(data["reactions"]) <= 6000
    assert 1000 <= len(data["compounds"]) <= 5000
    # All keys are R-numbers / C-numbers
    assert all(k.startswith("R") for k in data["reactions"])
    assert all(k.startswith("C") for k in data["compounds"])


@pytest.mark.slow
def test_metabolism_adapter_against_real_xrefs():
    if not XREFS_OUT.exists():
        pytest.skip(f"{XREFS_OUT} missing — run step 6 first")
    # Check that the xrefs file has compound data (not just reactions from partial KEGG cache)
    import json as _json
    xrefs_check = _json.loads(XREFS_OUT.read_text())
    if len(xrefs_check.get("compounds", {})) == 0:
        pytest.skip(
            "kegg_metabolism_xrefs.json has no compounds — re-run step 6 after full KEGG download"
        )
    from multiomics_kg.adapters.metabolism_adapter import MultiMetabolismAdapter

    adapter = MultiMetabolismAdapter(
        genome_config_file=str(PROJECT_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"),
    )
    nodes = list(adapter.get_nodes())
    assert any(n[1] == "reaction" for n in nodes)
    assert any(n[1] == "metabolite" for n in nodes)

    edges = list(adapter.get_edges())
    labels = {e[3] for e in edges}
    assert labels == {"gene_catalyzes_reaction", "reaction_has_metabolite", "reaction_in_kegg_pathway"}
