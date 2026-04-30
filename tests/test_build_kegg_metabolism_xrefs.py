"""Unit tests for step 6: build_kegg_metabolism_xrefs."""
from __future__ import annotations

import json

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def test_collect_gene_reachable_kegg_ids(tmp_path, monkeypatch):
    # Two synthetic strains
    s1 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MED4"
    s2 = tmp_path / "cache" / "data" / "Prochlorococcus" / "genomes" / "MIT9301"
    _write_strain_annotations(s1, {
        "PMM0001": {"kegg_reactions": ["R00200", "R00010"]},
        "PMM0002": {"kegg_reactions": ["R00010"]},
        "PMM0003": {},  # no reactions
    })
    _write_strain_annotations(s2, {
        "P9301_001": {"kegg_reactions": ["R12345"]},
    })

    # Stub: load_genome_rows returns rows pointing at the two strains
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [
        {"data_dir": str(s1)},
        {"data_dir": str(s2)},
    ])

    # Stub kegg_data.json with a tiny reaction_to_compounds map
    kegg_data = {
        "reaction_names": {
            "R00200": "pyruvate kinase reaction",
            "R00010": "phosphofructokinase",
            "R12345": "test rxn",
        },
        "compound_names": {
            "C00074": "phosphoenolpyruvate",
            "C00008": "ADP",
            "C00031": "D-glucose",
            "C00002": "ATP",
            "C99999": "test cpd",
        },
        "reaction_to_compounds": {
            "R00200": ["C00074", "C00008"],
            "R00010": ["C00031", "C00002"],
            "R12345": ["C99999"],
            "R_unused": ["C_should_not_appear"],
        },
        "reaction_to_pathways": {
            "R00200": ["ko00010"],
        },
    }

    rxns, cpds = bx._collect_gene_reachable_ids(kegg_data)

    # Only gene-reachable R-numbers
    assert rxns == {"R00200", "R00010", "R12345"}
    # C-numbers reachable from those reactions
    assert cpds == {"C00074", "C00008", "C00031", "C00002", "C99999"}
    # R_unused / C_should_not_appear are pruned
    assert "R_unused" not in rxns
    assert "C_should_not_appear" not in cpds


def test_collect_handles_missing_kegg_reactions_field(tmp_path, monkeypatch):
    s = tmp_path / "strain"
    _write_strain_annotations(s, {
        "g1": {},  # no kegg_reactions key
        "g2": {"kegg_reactions": []},  # empty list
        "g3": {"kegg_reactions": ["R00001"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(s)}])
    kegg_data = {"reaction_to_compounds": {"R00001": ["C00001"]}}
    rxns, cpds = bx._collect_gene_reachable_ids(kegg_data)
    assert rxns == {"R00001"}
    assert cpds == {"C00001"}
