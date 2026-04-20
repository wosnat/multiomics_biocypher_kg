"""End-to-end Biller 2018 integration test for ObservationsAdapter.

Runs against the real paperconfig at
`data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`.
Asserts expected DerivedMetric + binding-edge + measurement-edge counts
match the Plan 2 budget.

Requires Plan 1's _resolved.csv files to be present (run prepare_data.sh
step 4 if they're missing).
"""
from pathlib import Path

import pytest

from multiomics_kg.adapters.observations_adapter import ObservationsAdapter


BILLER_PC = Path("data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml")
BILLER_DOI = "10.1128/mSystems.00040-18"


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present (Plan 1 author step)",
)
def test_biller_2018_dm_nodes():
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    nodes = adapter.get_nodes()
    dm_nodes = [n for n in nodes if n[1] == "derived_metric"]
    # 4 entries x 2/2/2/1 metrics = 7 DerivedMetric nodes
    assert len(dm_nodes) == 7

    metric_types = sorted({props["metric_type"] for _, _, props in dm_nodes})
    assert metric_types == [
        "darkness_survival_class",
        "periodic_in_axenic_LD",
        "periodic_in_axenic_extended_darkness",
        "periodic_in_coculture_LD",
        "periodic_in_coculture_extended_darkness",
    ]
    # Sanity: 4 distinct DMs carry value_kind=boolean from NATL2A axenic + coculture;
    # 2 boolean from MIT1002; 1 categorical from S5
    by_kind = {}
    for _, _, props in dm_nodes:
        by_kind.setdefault(props["value_kind"], 0)
        by_kind[props["value_kind"]] += 1
    assert by_kind == {"boolean": 6, "categorical": 1}


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_binding_edges():
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    by_type = {}
    for e in edges:
        by_type.setdefault(e[3], 0)
        by_type[e[3]] += 1

    # 7 DMs -> 7 pub + 7 exp + 7 org binding edges
    assert by_type.get("publication_has_derived_metric", 0) == 7
    assert by_type.get("experiment_has_derived_metric", 0) == 7
    assert by_type.get("derived_metric_belongs_to_organism", 0) == 7


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_flag_edge_count_in_band():
    """Boolean measurement edges - expected band 3,500 - 4,700."""
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    flag_edges = [e for e in edges if e[3] == "derived_metric_flags_gene"]
    assert 3500 <= len(flag_edges) <= 4700, (
        f"derived_metric_flags_gene count {len(flag_edges)} outside expected band"
    )


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_classifies_edge_count():
    """Categorical measurement edges - S5 has 258 resolved rows (verified)."""
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
    }
    edges = adapter.get_edges()
    cls_edges = [e for e in edges if e[3] == "derived_metric_classifies_gene"]
    # Tolerance: 250-270 rows depending on resolver edge cases
    assert 250 <= len(cls_edges) <= 270, (
        f"derived_metric_classifies_gene count {len(cls_edges)} outside expected band"
    )


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_no_numeric_edges():
    """Biller 2018 has no numeric metrics - DM quantifies_gene count must be 0."""
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    q_edges = [e for e in edges if e[3] == "derived_metric_quantifies_gene"]
    assert len(q_edges) == 0
