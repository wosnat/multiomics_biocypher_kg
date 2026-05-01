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
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"
KEGG_DATA_OUT = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_data.json"


@pytest.mark.slow
def test_step6_smoke_run_against_real_cache():
    if not RESOLVER_DB.exists():
        pytest.skip(f"{RESOLVER_DB} missing — run `bash scripts/refresh_mnx.sh` first")

    result = subprocess.run(
        [sys.executable, "-m", "multiomics_kg.download.build_kegg_metabolism_xrefs", "--force"],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    assert result.returncode == 0, f"step 6 failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    assert KEGG_DATA_OUT.exists()
    data = json.loads(KEGG_DATA_OUT.read_text())

    # New nested structure: top-level keys for each entity type.
    assert set(data.keys()) >= {
        "kos", "pathways", "subcategories", "categories", "reactions", "compounds",
    }
    assert 2000 <= len(data["reactions"]) <= 6000
    assert 1000 <= len(data["compounds"]) <= 5000
    # All keys are R-numbers / C-numbers
    assert all(k.startswith("R") for k in data["reactions"])
    assert all(k.startswith("C") for k in data["compounds"])
    # Sanity-check pathway pruning: every compound.pathways entry must be in pathways set.
    pathway_set = set(data["pathways"])
    for cpd in data["compounds"].values():
        for p in cpd.get("pathways", []):
            assert p in pathway_set, f"unevidenced pathway {p} in compound list"


@pytest.mark.slow
def test_metabolism_adapter_against_real_xrefs():
    if not KEGG_DATA_OUT.exists():
        pytest.skip(f"{KEGG_DATA_OUT} missing — run step 6 first")
    # Check that the cache has compound data (not just a partial KEGG cache)
    data = json.loads(KEGG_DATA_OUT.read_text())
    if len(data.get("compounds", {})) == 0:
        pytest.skip(
            "kegg_data.json has no compounds — re-run step 6 after full KEGG download"
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
