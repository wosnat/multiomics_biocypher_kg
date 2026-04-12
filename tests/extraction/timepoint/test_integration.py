"""Integration test covering extract → merge → remap with a mocked LLM."""
import json
from pathlib import Path
from unittest.mock import MagicMock

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extract import extract_one_paper
from multiomics_kg.extraction.timepoint.merge import merge_one_paper
from multiomics_kg.extraction.timepoint.remap import remap_value


def test_full_flow_extract_merge_remap(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Fake 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: X
      organism: Prochlorococcus MED4
      treatment_type: [carbon]
      omics_type: RNASEQ
  supplementary_materials:
    tbl:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_a
          experiment: exp1
          name_col: gene
          logfc_col: log2fc
""")
    (tmp_path / "data").mkdir(exist_ok=True)
    (tmp_path / "data" / "fake.pdf").write_text("x")

    # LLM emits other:heat_acclimated
    fake_llm = MagicMock()
    fake_llm.call.return_value = ({
        "analyses": [{
            "analysis_id": "DE_a",
            "timepoint": "acclimated",
            "timepoint_hours": None,
            "growth_phase": "other:heat_acclimated",
            "self_assessment": "high",
            "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "Methods"}],
            "source_figures": [],
        }]
    }, {"input_tokens": 50, "output_tokens": 100, "model": "fake"})

    # 1. Extract
    extract_one_paper(tmp_path, llm_client=fake_llm)

    # 2. Merge (writes other:heat_acclimated verbatim)
    merge_one_paper(tmp_path, force=False)

    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    assert data["publication"]["supplementary_materials"]["tbl"][
        "statistical_analyses"][0]["growth_phase"] == "other:heat_acclimated"

    # 3. Remap to canonical (simulating enum promotion). For the test we
    #    bypass the is_valid_growth_phase check by remapping to an
    #    existing canonical value.
    remap_value([tmp_path], "other:heat_acclimated", "acclimated_steady_state")

    with open(pc) as f:
        data = yaml.load(f)
    assert data["publication"]["supplementary_materials"]["tbl"][
        "statistical_analyses"][0]["growth_phase"] == "acclimated_steady_state"

    ext = json.loads((tmp_path / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "acclimated_steady_state"
    assert "Remapped from" in ext["analyses"][0]["assessment_notes"]
