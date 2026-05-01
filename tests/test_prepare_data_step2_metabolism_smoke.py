"""End-to-end smoke test for Phase 1.1B against real cache files.

Assumes the heavy MNX resolver + TCDB/CAZy hierarchies have already been built
via `bash scripts/refresh_mnx.sh` (one-time, ~30 min, ~2.5 GB SQLite). The test
only re-runs the MED4 gene-annotation merge and validates its metabolism report.

If any required fixture is missing, the test is skipped.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


REQUIRED_FIXTURES = [
    "cache/data/mnx/metabolite_resolver.db",
    "cache/data/mnx/metabolite_id_mapping_report.json",
    "cache/data/tcdb/tcdb_hierarchy.json",
    "cache/data/cazy/cazy_hierarchy.json",
    "cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations",
    "cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv",
]


@pytest.mark.slow
def test_phase_1_1b_full_pipeline_med4():
    """Run step 2 for MED4 against pre-built MNX cache. Assert metabolism report sane."""
    project_root = Path(__file__).parent.parent
    for rel in REQUIRED_FIXTURES:
        if not (project_root / rel).exists():
            pytest.skip(f"Missing fixture: {rel} (run `bash scripts/refresh_mnx.sh` first)")

    report = json.loads((project_root / "cache/data/mnx/metabolite_id_mapping_report.json").read_text())
    assert report["compound_count"] >= 100_000

    # Run the gene annotation merge for MED4
    subprocess.run(
        ["uv", "run", "python", "-m",
         "multiomics_kg.download.build_gene_annotations", "--strains", "MED4", "--force"],
        cwd=project_root, check=True,
    )

    # 3. Assert metabolism report shows expected sanity values for MED4
    strain_dir = project_root / "cache/data/Prochlorococcus/genomes/MED4"
    report_path = strain_dir / "step2_metabolism_report.json"
    assert report_path.exists()
    metabolism = json.loads(report_path.read_text())

    assert metabolism["strain"] == "MED4"
    assert metabolism["gene_count"] >= 1500  # MED4 has 1838 genes
    # 526 KEGG_Reaction-annotated genes; spec lower bound >= 400 unique kept R-numbers.
    assert metabolism["kegg_reactions"]["kept_unique_r_numbers"] >= 400
    # 121 KEGG_TC-annotated genes; spec lower bound 50.
    assert metabolism["transporter_classification"]["validated_unique"] >= 50
    # CAZy coverage on MED4 is sparse (1.3%); spec lower bound 5.
    assert metabolism["cazy_ids"]["validated_unique"] >= 5
