"""Unit tests for the InterPro adapter (faceted calls.json consumption, precomputed
per-entry evidence, ancestor walk) plus a fast real-cache sanity check on MED4."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.adapters.interpro_adapter import (
    InterproAnnotationAdapter,
    MultiInterproAnnotationAdapter,
    _ec_node_id,
    _interpro_node_id,
    kept_ids,
)

MED4_DIR = Path("cache/data/Prochlorococcus/genomes/MED4")


def _write_strain(tmp_path: Path, genes: dict, calls: dict, strain: str = "TESTSTRAIN") -> Path:
    """Write a minimal merged-annotations JSON + faceted calls.json for one strain."""
    genome_dir = tmp_path / strain
    (genome_dir / "interproscan").mkdir(parents=True)
    with open(genome_dir / "gene_annotations_merged.json", "w", encoding="utf-8") as fh:
        json.dump(genes, fh)
    with open(genome_dir / "interproscan" / f"{strain}.interproscan.calls.json", "w", encoding="utf-8") as fh:
        json.dump(calls, fh)
    return genome_dir


# ── house rule R1: emitted `libraries` / `evalue_library` are lowercase and
# declared ─────────────────────────────────────────────────────────────────
#
# The calls.json artifact stores InterProScan's NATIVE (UPPERCASE) library
# names — see the "PFAM" keys elsewhere in this file. R1 normalization
# happens here, at the adapter's emission boundary
# (`InterproAnnotationAdapter.get_edges`), not inside the parser.

def test_edge_libraries_and_evalue_library_are_lowercase_and_declared(tmp_path):
    from multiomics_kg.utils.controlled_vocab import VOCAB

    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "interpro_entries": ["IPR000001"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "interpro_entries": {
                "IPR000001": {
                    "type": "FAMILY",
                    "start": 1,
                    "end": 36,
                    "evalue": 4.1e-18,
                    "evalue_library": "PFAM",
                    "libraries": ["HAMAP", "PFAM"],
                    "match_count": 2,
                }
            },
            "libraries": {},
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = InterproAnnotationAdapter(genome_dir)
    props = list(a.get_edges())[0][4]
    assert props["libraries"] == ["hamap", "pfam"]
    assert props["evalue_library"] == "pfam"
    for value in props["libraries"]:
        VOCAB.check("Gene_has_interpro_entry", "libraries", value)
    VOCAB.check("Gene_has_interpro_entry", "libraries", props["evalue_library"])


def test_node_interpro_type_is_lowercase_and_declared(tmp_path):
    from multiomics_kg.utils.controlled_vocab import VOCAB

    genes = {"LT001": {"protein_id": "WP_000000001.1", "interpro_entries": ["IPR000001"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(genome_dir) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = {
        "IPR000001": {"parent": None, "level": 0, "type": "HOMOLOGOUS_SUPERFAMILY", "name": "fam"},
    }
    m._observed_ids = lambda: {"IPR000001"}
    nodes = {nid: props for nid, _label, props in m.get_nodes()}
    interpro_type = nodes["interpro:IPR000001"]["interpro_type"]
    assert interpro_type == "homologous_superfamily"
    VOCAB.check("InterproEntry", "interpro_type", interpro_type)


# ── pure: ancestor pruning ────────────────────────────────────────────────────

REF = {
    "IPR000004": {"parent": "IPR000002", "level": 2},
    "IPR000002": {"parent": "IPR000001", "level": 1},
    "IPR000001": {"parent": None, "level": 0},
    "IPR000685": {"parent": None, "level": 0},
}


def test_kept_ids_pulls_full_ancestor_chain():
    kept = kept_ids({"IPR000004"}, REF)
    assert kept == {"IPR000004", "IPR000002", "IPR000001"}  # + ancestors


def test_kept_ids_entry_absent_from_reference_kept_as_is():
    kept = kept_ids({"IPR999999"}, REF)
    assert kept == {"IPR999999"}  # no ancestors, but itself kept (node still emitted)


def test_kept_ids_cycle_safe():
    cyclic = {"A": {"parent": "B"}, "B": {"parent": "A"}}
    assert kept_ids({"A"}, cyclic) == {"A", "B"}  # terminates


# ── faceted calls.json → edge props (Task 12 core) ────────────────────────────

def test_edge_props_from_rollup_with_evalue(tmp_path):
    """A normal rollup entry: evalue present => evalue_library present, no score key."""
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "interpro_entries": ["IPR000001"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "interpro_entries": {
                "IPR000001": {
                    "type": "FAMILY",
                    "start": 1,
                    "end": 36,
                    "evalue": 4.1e-18,
                    "evalue_library": "PFAM",
                    "libraries": ["HAMAP", "PFAM"],
                    "match_count": 2,
                }
            },
            "libraries": {},
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = InterproAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    eid, src, tgt, label, props = edges[0]
    assert label == "gene_has_interpro_entry"
    assert src == "ncbigene:LT001"
    assert tgt == "interpro:IPR000001"
    assert props["start"] == 1
    assert props["end"] == 36
    assert props["evalue"] == 4.1e-18
    assert props["evalue_library"] == "pfam"
    assert props["libraries"] == ["hamap", "pfam"]
    assert props["match_count"] == 2
    assert "score" not in props


def test_edge_props_no_evalue_omits_evalue_library(tmp_path):
    """Profile/pattern-only entries: evalue null => neither evalue nor evalue_library key."""
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "interpro_entries": ["IPR000001"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "interpro_entries": {
                "IPR000001": {
                    "type": "HOMOLOGOUS_SUPERFAMILY",
                    "start": 1,
                    "end": 35,
                    "evalue": None,
                    "evalue_library": None,
                    "libraries": ["SUPERFAMILY"],
                    "match_count": 1,
                }
            },
            "libraries": {},
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = InterproAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    props = edges[0][4]
    assert "evalue" not in props
    assert "evalue_library" not in props
    assert "score" not in props
    assert props["match_count"] == 1
    assert props["libraries"] == ["superfamily"]


def test_edge_fail_soft_when_merged_seed_missing_from_rollup(tmp_path):
    """Merged interpro_entries is authority: an entry absent from the calls.json
    rollup (merge/calls skew) still gets an edge, fail-soft with empty evidence."""
    genes = {
        "LT001": {
            "protein_id": "WP_000000001.1",
            "interpro_entries": ["IPR000001", "IPR999999"],
        }
    }
    calls = {
        "WP_000000001.1": {
            "interpro_entries": {
                "IPR000001": {
                    "type": "FAMILY",
                    "start": 1,
                    "end": 10,
                    "evalue": 1e-10,
                    "evalue_library": "PFAM",
                    "libraries": ["PFAM"],
                    "match_count": 1,
                }
            },
            "libraries": {},
        }
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = InterproAnnotationAdapter(genome_dir)
    edges = {e[2]: e[4] for e in a.get_edges()}
    assert edges["interpro:IPR000001"]["match_count"] == 1
    assert edges["interpro:IPR999999"] == {"match_count": 0, "libraries": []}


def test_edge_fail_soft_when_protein_has_no_calls_at_all(tmp_path):
    """Protein entirely missing from calls.json (e.g. no calls.json ran) still
    emits fail-soft edges for every merged-authority entry."""
    genes = {
        "LT001": {
            "protein_id": "WP_999999999.1",
            "interpro_entries": ["IPR000001"],
        }
    }
    genome_dir = _write_strain(tmp_path, genes, {})
    a = InterproAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    assert edges[0][4] == {"match_count": 0, "libraries": []}


def test_gene_with_no_interpro_entries_yields_no_edges(tmp_path):
    genes = {"LT001": {"protein_id": "WP_000000001.1", "interpro_entries": []}}
    genome_dir = _write_strain(tmp_path, genes, {})
    a = InterproAnnotationAdapter(genome_dir)
    assert list(a.get_edges()) == []


def test_gene_with_no_protein_id_yields_no_edges(tmp_path):
    genes = {"LT001": {"protein_id": "", "interpro_entries": ["IPR000001"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    a = InterproAnnotationAdapter(genome_dir)
    assert list(a.get_edges()) == []


# ── PFAM bridge from the PFAM facet ────────────────────────────────────────────

def test_get_pfam_to_interpro_reads_pfam_facet(tmp_path):
    calls = {
        "WP_000000001.1": {
            "interpro_entries": {},
            "libraries": {
                "PFAM": [
                    {"accession": "PF02532.1", "ipr": "IPR003686", "start": 1, "end": 36},
                ],
                "NCBIFAM": [
                    {"accession": "NF002735", "ipr": None, "start": 1, "end": 38},
                ],
            },
        }
    }
    genome_dir = _write_strain(tmp_path, {}, calls)
    a = InterproAnnotationAdapter(genome_dir)
    mapping = a.get_pfam_to_interpro()
    # version suffix stripped, NCBIFAM (non-PFAM library / null ipr) excluded
    assert mapping == {"PF02532": "IPR003686"}


# ── node description (Task 6 sparse field) ─────────────────────────────────────

def test_node_carries_description_when_reference_provides_it(tmp_path):
    genes = {"LT001": {"protein_id": "WP_000000001.1", "interpro_entries": ["IPR000001"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(genome_dir) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = {
        "IPR000001": {"parent": None, "level": 0, "type": "FAMILY", "name": "fam",
                       "description": "A narrative description."},
    }
    m._observed_ids = lambda: {"IPR000001"}
    nodes = {nid: props for nid, _label, props in m.get_nodes()}
    assert nodes["interpro:IPR000001"]["description"] == "A narrative description."


def test_node_omits_description_when_reference_has_none(tmp_path):
    genes = {"LT001": {"protein_id": "WP_000000001.1", "interpro_entries": ["IPR000001"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(genome_dir) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = {
        "IPR000001": {"parent": None, "level": 0, "type": "FAMILY", "name": "fam"},
    }
    m._observed_ids = lambda: {"IPR000001"}
    nodes = {nid: props for nid, _label, props in m.get_nodes()}
    assert "description" not in nodes["interpro:IPR000001"]


# ── kept_node_accessions() (Task 14 consumes this) ─────────────────────────────

def test_kept_node_accessions_matches_observed_plus_ancestors(tmp_path):
    genes = {"LT001": {"protein_id": "WP_000000001.1", "interpro_entries": ["IPR000004"]}}
    genome_dir = _write_strain(tmp_path, genes, {})
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(genome_dir) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = REF
    m._observed_ids = lambda: {"IPR000004"}
    assert m.kept_node_accessions() == {"IPR000004", "IPR000002", "IPR000001"}


# ── fast real-cache sanity (MED4) ─────────────────────────────────────────────

@pytest.mark.skipif(not MED4_DIR.exists(), reason="MED4 cache not present")
def test_med4_edges_wellformed():
    a = InterproAnnotationAdapter(MED4_DIR)
    edges = list(a.get_edges())
    assert len(edges) > 1000
    for eid, src, tgt, label, props in edges[:500]:
        assert label == "gene_has_interpro_entry"
        assert src.startswith("ncbigene:")
        assert tgt.startswith("interpro:IPR")
        assert "score" not in props
        assert "match_count" in props
        assert "libraries" in props
        if "start" in props and "end" in props:
            assert props["start"] <= props["end"]        # valid envelope
        if "evalue" in props:
            assert props["evalue"] >= 0
            assert "evalue_library" in props
        else:
            assert "evalue_library" not in props


# ── Layer A: InterproEntry → EC/CAZy router edges (design 2026-08-10) ──────────

def test_observed_ec_node_ids_wellformed():
    a = InterproAnnotationAdapter(MED4_DIR)
    obs = a.observed_ec_node_ids()
    assert obs, "MED4 has EC-annotated genes"
    assert all(nid.startswith("ec:") for nid in obs)


def test_layer_a_related_ec_pruned_and_marked(tmp_path):
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")
    m = MultiInterproAnnotationAdapter(
        str(cfg), pfam_node_ids=set(), ec_node_ids=_all_med4_ec_node_ids(),
    )

    strain = m._strain_adapters[0]
    observed = strain.observed_ec_node_ids()
    # pick a real MED4 EC token
    some_ec = next(
        ec for g in strain._genes.values() for ec in (g.get("ec_numbers") or [])
        if _ec_node_id(ec) in observed
    )
    # FAMILY single-EC → not ambiguous; DOMAIN multi-EC (2nd EC unobserved) → ambiguous + pruned
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": [some_ec]},
        "IPR900002": {"type": "DOMAIN", "name": "dom", "ec_numbers": [some_ec, "9.9.9.9"]},
    }
    m._observed_ids = lambda: {"IPR900001", "IPR900002"}

    ec_edges = [e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"]
    by_src = {e[1]: e[4] for e in ec_edges}
    fam_id = _interpro_node_id("IPR900001")
    dom_id = _interpro_node_id("IPR900002")
    # Layer A edges carry no properties (deleted 2026-08-16, spec §3 R3): the old
    # `ambiguous`/`source_db` fields are gone; multiplicity/type are derivable
    # from the graph instead (out-degree + n.interpro_type).
    assert by_src[fam_id] == {}
    assert by_src[dom_id] == {}
    # the unobserved 9.9.9.9 was pruned (dangling-proof): only one edge per entry
    assert sum(1 for e in ec_edges if e[1] == dom_id) == 1


def _all_med4_ec_node_ids() -> set[str]:
    """Every EC node id MED4's genes reference — stands in for the real Expasy
    node set in tests that only need 'the endpoint exists'."""
    return InterproAnnotationAdapter(MED4_DIR).observed_ec_node_ids() | {
        _ec_node_id("9.9.9.9")
    }


def test_layer_a_related_ec_pruned_to_existing_ec_nodes(tmp_path):
    """An EC a gene carries but Expasy has NO node for (obsolete InterPro xref,
    e.g. 1.2.8.1) must not produce a Layer-A router edge — it would dangle."""
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")

    strain = InterproAnnotationAdapter(MED4_DIR)
    observed = strain.observed_ec_node_ids()
    some_ec = next(
        ec for g in strain._genes.values() for ec in (g.get("ec_numbers") or [])
        if _ec_node_id(ec) in observed
    )
    # The gene set carries `some_ec` AND the obsolete one, but only `some_ec`
    # has an EcNumber node.
    obsolete = "1.2.8.1"
    m = MultiInterproAnnotationAdapter(
        str(cfg), pfam_node_ids=set(), ec_node_ids={_ec_node_id(some_ec)},
    )
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": [some_ec]},
        "IPR900003": {"type": "FAMILY", "name": "obs", "ec_numbers": [obsolete]},
    }
    m._observed_ids = lambda: {"IPR900001", "IPR900003"}
    # Pretend a gene carries the obsolete EC, so `observed_ec` alone would keep it.
    for a in m._strain_adapters:
        a.observed_ec_node_ids = lambda: observed | {_ec_node_id(obsolete)}

    ec_edges = [e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"]
    targets = {e[2] for e in ec_edges}
    assert _ec_node_id(some_ec) in targets
    assert _ec_node_id(obsolete) not in targets


def test_layer_a_related_ec_suppressed_without_ec_node_ids(tmp_path):
    """No EC node set injected → no Layer-A EC edges (mirrors pfam_node_ids)."""
    cfg = tmp_path / "genomes.csv"
    cfg.write_text("data_dir\n" + str(MED4_DIR) + "\n")
    m = MultiInterproAnnotationAdapter(str(cfg), pfam_node_ids=set())
    m._reference = {
        "IPR900001": {"type": "FAMILY", "name": "fam", "ec_numbers": ["2.7.7.7"]},
    }
    m._observed_ids = lambda: {"IPR900001"}
    assert not [
        e for e in m.get_edges() if e[3] == "interpro_entry_related_to_ec_number"
    ]


def test_layer_a_edges_carry_no_properties():
    """rev 4/5: ambiguous is derivable and source_db is a constant."""
    import pathlib
    src = pathlib.Path("multiomics_kg/adapters/interpro_adapter.py").read_text()
    assert '"ambiguous"' not in src
    assert '"source_db"' not in src
