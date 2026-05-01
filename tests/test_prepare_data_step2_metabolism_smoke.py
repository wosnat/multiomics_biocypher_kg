"""End-to-end smoke test for Phase 1.1B against real cache files.

Assumes:
- cache/data/mnx/, cache/data/tcdb/ are populated (sub-step 6 has run)
- cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations exists
- cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv exists

If any of those is missing, the test is skipped.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


REQUIRED_FIXTURES = [
    "cache/data/mnx/chem_prop.tsv",
    "cache/data/mnx/chem_xref.tsv",
    "cache/data/mnx/reac_prop.tsv",
    "cache/data/mnx/reac_xref.tsv",
    "cache/data/tcdb/families.tsv",
    "cache/data/tcdb/superfamilies.tsv",
    "cache/data/tcdb/substrates.tsv",
    "cache/data/Prochlorococcus/genomes/MED4/eggnog/MED4.emapper.annotations",
    "cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv",
]


@pytest.mark.slow
def test_phase_1_1b_full_pipeline_med4():
    """Run scripts/refresh_mnx.sh + step 2 against real cache. Assert metabolism report sane."""
    project_root = Path(__file__).parent.parent
    for rel in REQUIRED_FIXTURES:
        if not (project_root / rel).exists():
            pytest.skip(f"Missing fixture: {rel}")

    # 1. Build resolver + hierarchies (was sub-step 7; now scripts/refresh_mnx.sh)
    subprocess.run(
        ["bash", "scripts/refresh_mnx.sh", "--force"],
        cwd=project_root, check=True,
    )

    assert (project_root / "cache/data/mnx/metabolite_resolver.db").exists()
    assert (project_root / "cache/data/tcdb/tcdb_hierarchy.json").exists()
    assert (project_root / "cache/data/cazy/cazy_hierarchy.json").exists()
    report = json.loads((project_root / "cache/data/mnx/metabolite_id_mapping_report.json").read_text())
    assert report["compound_count"] >= 100_000

    # 2. Run the gene annotation merge for MED4
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
