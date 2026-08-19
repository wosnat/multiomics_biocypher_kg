"""Unit tests for the NCBIfam adapter (NcbifamFamily nodes + gene edges + the
InterPro bridge). Mirrors ``tests/test_interpro_adapter.py`` structure."""

from __future__ import annotations

import json
from pathlib import Path

from multiomics_kg.adapters.ncbifam_adapter import (
    NcbifamAnnotationAdapter,
    MultiNcbifamAdapter,
    _ncbifam_node_id,
)


def _write_strain(tmp_path: Path, genes: dict, calls: dict, strain: str = "TESTSTRAIN") -> Path:
    """Write a minimal merged-annotations JSON + faceted calls.json for one strain."""
    genome_dir = tmp_path / strain
    (genome_dir / "interproscan").mkdir(parents=True)
    with open(genome_dir / "gene_annotations_merged.json", "w", encoding="utf-8") as fh:
        json.dump(genes, fh)
    with open(genome_dir / "interproscan" / f"{strain}.interproscan.calls.json", "w", encoding="utf-8") as fh:
        json.dump(calls, fh)
    return genome_dir


# ── node id helper ──────────────────────────────────────────────────────────

def test_ncbifam_node_id_uses_colon_curie_form():
    # House-minted colon CURIE (KG-SYNC-002): `ncbifam` is unregistered
    # everywhere (bioregistry/identifiers.org/Biolink), but the ontology's
    # peers (tcdb/interpro/pfam/merops.*) all use colon CURIEs, so consumers
    # get one id grammar.
    assert _ncbifam_node_id("TIGR00198") == "ncbifam:TIGR00198"
    assert _ncbifam_node_id("NF002735") == "ncbifam:NF002735"


# ── per-strain gene edges ────────────────────────────────────────────────────

def test_gene_edge_carries_evalue_and_score(tmp_path):
    """Single-library rule: NCBIfam edges keep BOTH evalue and score (unlike
    interpro's Gene_has_interpro_entry, which omits score)."""
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "ncbifam_ids": ["TIGR00198"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": "IPR010987",
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = NcbifamAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    eid, src, tgt, label, props = edges[0]
    assert label == "gene_has_ncbifam_family"
    assert src == "ncbigene:LT001"
    assert tgt == "ncbifam:TIGR00198"
    assert props["start"] == 5
    assert props["end"] == 480
    assert props["evalue"] == 1e-200
    assert props["score"] == 850.5


def test_gene_edge_best_row_selected_by_min_evalue(tmp_path):
    """Two facet rows for the SAME accession on the same protein => pick the
    lower-evalue row's evidence; still exactly ONE edge for the (gene, acc)."""
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "ncbifam_ids": ["NF002735"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "NF002735", "name": "psbI", "ipr": None,
                     "start": 1, "end": 38, "evalue": 1e-10, "score": 40.0},
                    {"accession": "NF002735", "name": "psbI", "ipr": None,
                     "start": 2, "end": 39, "evalue": 1e-30, "score": 60.0},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = NcbifamAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    props = edges[0][4]
    assert props["evalue"] == 1e-30  # the lower one
    assert props["start"] == 2
    assert props["end"] == 39
    assert props["score"] == 60.0


def test_gene_edge_fail_soft_when_facet_row_missing(tmp_path):
    """Merged ncbifam_ids is authority: an accession absent from the calls.json
    facet rows (merge/calls skew) still gets an edge, fail-soft with {} props."""
    genes = {
        "LT001": {
            "protein_id": "WP_999999999.1",
            "ncbifam_ids": ["TIGR00198"],
        }
    }
    genome_dir = _write_strain(tmp_path, genes, {})
    a = NcbifamAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    assert edges[0][4] == {}


def test_gene_edge_sparse_props_omit_nulls(tmp_path):
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "ncbifam_ids": ["NF002735"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "NF002735", "name": "x", "ipr": None,
                     "start": None, "end": None, "evalue": None, "score": None},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = NcbifamAnnotationAdapter(genome_dir)
    props = list(a.get_edges())[0][4]
    assert props == {}


def test_gene_edge_skips_falsy_accession_entries(tmp_path):
    """An empty-string/None entry in ncbifam_ids (data glitch) must not emit an
    edge to a bare `ncbifam_` node id -- mirrors the `if acc:` filter already
    applied in get_all_ncbifam_ids."""
    genes = {
        "LT001": {
            "protein_id": "WP_1.1",
            "ncbifam_ids": ["", "TIGR00198", None],
        }
    }
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": None,
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = NcbifamAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    assert edges[0][2] == "ncbifam:TIGR00198"


def test_gene_with_no_ncbifam_ids_yields_no_edges(tmp_path):
    genes = {"LT001": {"protein_id": "WP_000000001.1", "ncbifam_ids": []}}
    genome_dir = _write_strain(tmp_path, genes, {})
    a = NcbifamAnnotationAdapter(genome_dir)
    assert list(a.get_edges()) == []


def test_gene_with_no_protein_id_yields_no_edges(tmp_path):
    genes = {"LT001": {"protein_id": "", "ncbifam_ids": ["TIGR00198"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    a = NcbifamAnnotationAdapter(genome_dir)
    assert list(a.get_edges()) == []


def test_test_mode_caps_at_100(tmp_path):
    genes = {
        f"LT{i:04d}": {"protein_id": f"WP_{i:09d}.1", "ncbifam_ids": ["TIGR00198"]}
        for i in range(150)
    }
    genome_dir = _write_strain(tmp_path, genes, {})
    a = NcbifamAnnotationAdapter(genome_dir, test_mode=True)
    edges = list(a.get_edges())
    assert len(edges) == 100


# ── get_all_ncbifam_ids() ───────────────────────────────────────────────────

def test_get_all_ncbifam_ids(tmp_path):
    genes = {
        "LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198", "NF002735"]},
        "LT002": {"protein_id": "WP_2.1", "ncbifam_ids": ["TIGR00198"]},
        "LT003": {"protein_id": "WP_3.1"},
    }
    genome_dir = _write_strain(tmp_path, genes, {})
    a = NcbifamAnnotationAdapter(genome_dir)
    assert a.get_all_ncbifam_ids() == {"TIGR00198", "NF002735"}


# ── MultiNcbifamAdapter: nodes ───────────────────────────────────────────────

def _cfg(tmp_path: Path, genome_dir: Path) -> Path:
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(genome_dir) + "\n")
    return cfg


def test_node_shape_from_reference(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {
        "TIGR00198": {
            "name": "catalase/peroxidase HPI",
            "family_type": "equivalog",
            "gene_symbol": "katG",
            "description": "A curated narrative.",
        }
    }
    nodes = {nid: (label, props) for nid, label, props in m.get_nodes()}
    assert "ncbifam:TIGR00198" in nodes
    label, props = nodes["ncbifam:TIGR00198"]
    assert label == "ncbifam family"
    assert props["name"] == "catalase/peroxidase HPI"
    assert props["ncbifam_id"] == "TIGR00198"
    assert props["family_type"] == "equivalog"
    assert props["gene_symbol"] == "katG"
    assert props["description"] == "A curated narrative."
    assert props["level"] == 0


def test_node_sparse_props_omitted_when_absent(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["NF002735"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {
        "NF002735": {"name": "photosystem II reaction center protein I",
                     "family_type": "equivalog"},
    }
    props = {nid: p for nid, _l, p in m.get_nodes()}["ncbifam:NF002735"]
    assert "gene_symbol" not in props
    assert "description" not in props


def test_node_clean_str_applied(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["NF002735"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {
        "NF002735": {"name": "it's a pipe|test", "family_type": "equivalog",
                     "description": "another 'quote' | pipe"},
    }
    props = {nid: p for nid, _l, p in m.get_nodes()}["ncbifam:NF002735"]
    assert "'" not in props["name"]
    assert "|" not in props["name"]
    assert "'" not in props["description"]
    assert "|" not in props["description"]


def test_node_retired_accession_fallback_from_calls(tmp_path):
    """Accession observed on a gene but absent from the reference (retired) still
    gets a node, with `name` pulled from any strain's NCBIFAM facet row and
    `family_type` omitted."""
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR99999"]}}
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR99999", "name": "retired family name",
                     "ipr": None, "start": 1, "end": 10, "evalue": 1e-5, "score": 20.0},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {}  # accession retired -- absent from the current TSV
    props = {nid: p for nid, _l, p in m.get_nodes()}["ncbifam:TIGR99999"]
    assert props["name"] == "retired family name"
    assert "family_type" not in props
    assert props["ncbifam_id"] == "TIGR99999"
    assert props["level"] == 0


def test_get_nodes_observed_only_no_ancestors(tmp_path):
    """NcbifamFamily is flat/observed-only -- no ancestor walk (unlike InterPro)."""
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {
        "TIGR00198": {"name": "katG", "family_type": "equivalog"},
        "TIGR00357": {"name": "unrelated, unobserved family", "family_type": "equivalog"},
    }
    node_ids = {nid for nid, _l, _p in m.get_nodes()}
    assert node_ids == {"ncbifam:TIGR00198"}


def test_download_data_raises_when_reference_missing(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": []}}
    genome_dir = _write_strain(tmp_path, genes, {})
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)), cache_root=tmp_path / "nope")
    try:
        m.download_data()
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as exc:
        assert "ncbifam_reference.json" in str(exc)


# ── MultiNcbifamAdapter: interpro bridge ────────────────────────────────────

def test_bridge_edges_pruned_to_injected_interpro_kept_ids(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": "IPR010987",
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(
        str(_cfg(tmp_path, genome_dir)), interpro_kept_ids={"IPR010987"},
    )
    m._reference = {"TIGR00198": {"name": "katG", "family_type": "equivalog"}}
    bridge_edges = [e for e in m.get_edges() if e[3] == "ncbifam_family_in_interpro_entry"]
    assert len(bridge_edges) == 1
    eid, src, tgt, label, props = bridge_edges[0]
    assert src == "ncbifam:TIGR00198"
    assert tgt == "interpro:IPR010987"
    assert props == {}


def test_bridge_edges_dropped_when_acc_has_no_emitted_node(tmp_path):
    """Source-side dangling guard: calls.json can carry an (acc -> ipr) facet pair
    for an accession that is NOT in the strain's merged ncbifam_ids (the
    pre-Task-18 reality -- merged seeds are empty everywhere today, but calls.json
    facet rows exist independently). get_nodes() never emits a node for such an
    acc, so the bridge must skip it too, even when the target ipr is kept."""
    genes = {"LT001": {"protein_id": "WP_1.1"}}  # no ncbifam_ids at all
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00621", "name": "some family", "ipr": "IPR011344",
                     "start": 1, "end": 100, "evalue": 1e-50, "score": 200.0},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(
        str(_cfg(tmp_path, genome_dir)), interpro_kept_ids={"IPR011344"},
    )
    m._reference = {}
    bridge_edges = [e for e in m.get_edges() if e[3] == "ncbifam_family_in_interpro_entry"]
    assert bridge_edges == []


def test_bridge_edges_none_when_interpro_kept_ids_not_injected(tmp_path):
    """interpro_kept_ids=None (default) => no bridge edges emitted at all
    (injection contract, mirrors pfam_node_ids in interpro_adapter)."""
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": "IPR010987",
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {"TIGR00198": {"name": "katG", "family_type": "equivalog"}}
    bridge_edges = [e for e in m.get_edges() if e[3] == "ncbifam_family_in_interpro_entry"]
    assert bridge_edges == []


def test_bridge_edges_dropped_when_target_not_in_kept_set(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": "IPR010987",
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(
        str(_cfg(tmp_path, genome_dir)), interpro_kept_ids={"IPR999999"},  # different entry
    )
    m._reference = {"TIGR00198": {"name": "katG", "family_type": "equivalog"}}
    bridge_edges = [e for e in m.get_edges() if e[3] == "ncbifam_family_in_interpro_entry"]
    assert bridge_edges == []


def test_get_edges_includes_gene_edges_via_delegation(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "ncbifam_ids": ["TIGR00198"]}}
    calls = {
        "WP_1.1": {
            "libraries": {
                "NCBIFAM": [
                    {"accession": "TIGR00198", "name": "katG", "ipr": None,
                     "start": 5, "end": 480, "evalue": 1e-200, "score": 850.5},
                ]
            }
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    m = MultiNcbifamAdapter(str(_cfg(tmp_path, genome_dir)))
    m._reference = {"TIGR00198": {"name": "katG", "family_type": "equivalog"}}
    gene_edges = [e for e in m.get_edges() if e[3] == "gene_has_ncbifam_family"]
    assert len(gene_edges) == 1
    assert gene_edges[0][1] == "ncbigene:LT001"
    assert gene_edges[0][2] == "ncbifam:TIGR00198"
