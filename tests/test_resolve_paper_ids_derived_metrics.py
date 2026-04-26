"""Unit tests for resolve_derived_metrics_entry in resolve_paper_ids.

The resolver writes <stem>_resolved.csv next to the source, adding
a resolved_locus_tag column and a resolution_method column.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from multiomics_kg.download import resolve_paper_ids as rpi
from multiomics_kg.utils.gene_id_utils import MappingData


def _seed_mapping(genome_dir: Path) -> None:
    genome_dir.mkdir(parents=True, exist_ok=True)
    gene_id_mapping = {
        "version": 2,
        "strain": "NATL2A",
        "specific_lookup": {
            "PMN2A_RS00015": "PMN2A_RS00015",
            "PMN2A_RS00020": "PMN2A_RS00020",
        },
        "multi_lookup": {},
        "conflicts": {},
        "genes": {
            "PMN2A_RS00015": {"tier1_ids": ["PMN2A_RS00015"], "tier2_ids": [], "tier3_ids": []},
            "PMN2A_RS00020": {"tier1_ids": ["PMN2A_RS00020"], "tier2_ids": [], "tier3_ids": []},
        },
        "stats": {"n_genes": 2, "n_specific": 2, "n_multi": 0, "n_conflicts": 0, "passes": 1},
    }
    (genome_dir / "gene_id_mapping.json").write_text(json.dumps(gene_id_mapping))


def test_resolve_derived_metrics_entry_writes_resolved_csv(tmp_path, monkeypatch):
    src = tmp_path / "s4a.csv"
    src.write_text(
        "NCBI ID_2,flag\n"
        "PMN2A_RS00015,Y\n"
        "PMN2A_RS00020,\n"
        "UNKNOWN_ID,Y\n"
    )

    genome_dir = tmp_path / "cache/data/Prochlorococcus/genomes/NATL2A"
    _seed_mapping(genome_dir)

    monkeypatch.setattr(
        "multiomics_kg.utils.gene_id_utils.get_genome_dir",
        lambda organism, project_root: str(genome_dir),
    )
    monkeypatch.setattr(rpi, "get_genome_dir",
                        lambda organism, project_root=None: str(genome_dir))
    monkeypatch.setattr(rpi, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": str(src.relative_to(tmp_path)),
        "organism": "Prochlorococcus NATL2A",
        "experiment": "exp1",
        "name_col": "NCBI ID_2",
        "id_columns": [{"column": "NCBI ID_2", "id_type": "locus_tag_ncbi"}],
        "metrics": [{"metric_type": "periodic_in_axenic_LD", "value_kind": "boolean",
                     "value_col": "flag", "true_tokens": ["Y"]}],
    }

    result = rpi.resolve_derived_metrics_entry("Biller 2018", "dm_s4a", entry, force=True)
    assert result is not None
    assert not result.get("skipped")

    resolved = tmp_path / "s4a_resolved.csv"
    assert resolved.exists()

    out = pd.read_csv(resolved)
    assert "resolved_locus_tag" in out.columns
    assert "resolution_method" in out.columns
    assert out.loc[out["NCBI ID_2"] == "PMN2A_RS00015", "resolved_locus_tag"].iloc[0] == "PMN2A_RS00015"
    assert pd.isna(out.loc[out["NCBI ID_2"] == "UNKNOWN_ID", "resolved_locus_tag"].iloc[0])


def test_resolve_derived_metrics_entry_resolves_non_canonical_locus_tag_column(tmp_path, monkeypatch):
    """Regression: name_col='locus_tag' must still resolve when values are
    a non-canonical assembly's locus tags (e.g. RefSeq TX50_RS* for MED4 vs
    canonical PMM*). Earlier code skipped this case as an optimization,
    producing dangling Derived_metric_quantifies_gene edges at import time.
    """
    src = tmp_path / "table_s4.csv"
    src.write_text(
        "locus_tag,avg_reads\n"
        "TX50_RS05485,12615\n"   # alt-assembly id → maps to PMM1943
        "PMM0869,10732\n"        # already canonical
        "BOGUS,1\n"              # unresolvable
    )

    genome_dir = tmp_path / "cache/data/Prochlorococcus/genomes/MED4"
    genome_dir.mkdir(parents=True, exist_ok=True)
    (genome_dir / "gene_id_mapping.json").write_text(json.dumps({
        "version": 2,
        "strain": "MED4",
        "specific_lookup": {
            "TX50_RS05485": "PMM1943",
            "PMM1943": "PMM1943",
            "PMM0869": "PMM0869",
        },
        "multi_lookup": {},
        "conflicts": {},
        "genes": {
            "PMM1943": {"tier1_ids": ["PMM1943", "TX50_RS05485"], "tier2_ids": [], "tier3_ids": []},
            "PMM0869": {"tier1_ids": ["PMM0869"], "tier2_ids": [], "tier3_ids": []},
        },
        "stats": {"n_genes": 2, "n_specific": 3, "n_multi": 0, "n_conflicts": 0, "passes": 1},
    }))

    monkeypatch.setattr(rpi, "get_genome_dir",
                        lambda organism, project_root=None: str(genome_dir))
    monkeypatch.setattr(rpi, "PROJECT_ROOT", tmp_path)

    entry = {
        "type": "derived_metrics_table",
        "filename": str(src.relative_to(tmp_path)),
        "organism": "Prochlorococcus MED4",
        "experiment": "exp1",
        "name_col": "locus_tag",  # the trap: column named locus_tag, values aren't all canonical
        "id_columns": [{"column": "locus_tag", "id_type": "locus_tag_ncbi"}],
        "metrics": [{"metric_type": "vesicle_dna_avg_read_coverage", "value_kind": "numeric",
                     "value_col": "avg_reads", "rankable": "true", "has_p_value": "false"}],
    }

    result = rpi.resolve_derived_metrics_entry("Biller 2014", "s4", entry, force=True)
    assert result is not None
    assert not result.get("skipped"), f"resolver wrongly skipped: {result.get('reason')}"

    resolved = tmp_path / "table_s4_resolved.csv"
    assert resolved.exists()
    out = pd.read_csv(resolved)
    assert out.loc[out["locus_tag"] == "TX50_RS05485", "resolved_locus_tag"].iloc[0] == "PMM1943"
    assert out.loc[out["locus_tag"] == "PMM0869", "resolved_locus_tag"].iloc[0] == "PMM0869"
    assert pd.isna(out.loc[out["locus_tag"] == "BOGUS", "resolved_locus_tag"].iloc[0])
