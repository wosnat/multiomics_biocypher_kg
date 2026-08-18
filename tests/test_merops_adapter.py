"""Unit tests for the MEROPS adapter (MeropsFamily nodes + scored gene edges +
is-a hierarchy) and the Phase-2 vocabulary helpers in utils/merops_diamond.
Mirrors ``tests/test_ncbifam_adapter.py`` structure."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiomics_kg.adapters.merops_adapter import (
    MeropsAnnotationAdapter,
    MultiMeropsAnnotationAdapter,
    _merops_node_id,
)
from multiomics_kg.utils.merops_diamond import (
    call_class,
    catalytic_type_word,
    classify_code,
    edge_target_code,
    family_type,
    parse_clan_txt,
    type_example_names,
)


# ── vocabulary helpers ───────────────────────────────────────────────────────

def test_classify_code_levels():
    assert classify_code("SC") == (0, "merops_clan")
    assert classify_code("S14") == (1, "merops_family")
    assert classify_code("M103") == (1, "merops_family")
    assert classify_code("S08A") == (2, "merops_subfamily")
    assert classify_code("S08.001") is None  # identifiers are not node codes
    assert classify_code("") is None


def test_catalytic_type_words():
    assert catalytic_type_word("S14") == "serine"
    assert catalytic_type_word("M01") == "metallo"
    assert catalytic_type_word("PA") == "mixed"       # mixed-type clan
    assert catalytic_type_word("U32") == "unknown"
    assert catalytic_type_word("N06") == "asparagine_lyase"
    assert catalytic_type_word("I39") is None         # inhibitors: no catalytic type


def test_family_type():
    assert family_type("S14") == "peptidase"
    assert family_type("I39") == "inhibitor"
    assert family_type("IN") == "inhibitor"           # inhibitor clan


def test_call_class_verdicts():
    assert call_class({"entry_type": "peptidase", "best_hit_kind": "holotype"}) == "peptidase"
    assert call_class({"entry_type": "peptidase", "best_hit_kind": "putative"}) == "peptidase"
    assert call_class(
        {"entry_type": "peptidase", "best_hit_kind": "nonpeptidase_homolog"}
    ) == "nonpeptidase_homolog"
    # inhibitor family wins over hit kind
    assert call_class(
        {"entry_type": "inhibitor", "best_hit_kind": "nonpeptidase_homolog"}
    ) == "inhibitor"


def test_edge_target_code_tier23_uses_called_code():
    assert edge_target_code({"level_kind": "merops_family", "code": "S14"}) == "S14"
    assert edge_target_code({"level_kind": "merops_subfamily", "code": "S08A"}) == "S08A"


def test_edge_target_code_tier1_attaches_at_subfamily_else_family():
    # identifier-level call with a subfamily claim → subfamily node
    assert edge_target_code(
        {"level_kind": "merops_id", "code": "S26.014", "subfamily": "S26A", "family": "S26"}
    ) == "S26A"
    # identifier-level call in a family without subfamilies → family node
    assert edge_target_code(
        {"level_kind": "merops_id", "code": "M75.001", "subfamily": None, "family": "M75"}
    ) == "M75"


def test_edge_target_code_malformed_returns_none():
    assert edge_target_code({"level_kind": "merops_family", "code": "S08.001"}) is None
    assert edge_target_code({"level_kind": "merops_family", "code": None}) is None


def test_parse_clan_txt_skips_unassigned_sentinels():
    text = (
        "A-\t\tFamilies not assigned to clans\t\t\tpeptidase\tYes\n"
        "AA\t\tWater nucleophile; two Asp\t\tRef\tpeptidase\tYes\n"
        "IN\t\tInhibitor clan\t\tRef\tinhibitor\tYes\n"
    )
    clans = parse_clan_txt(text)
    assert "A-" not in clans
    assert clans["AA"] == {"description": "Water nucleophile; two Asp", "family_type": "peptidase"}
    assert clans["IN"]["family_type"] == "inhibitor"


def test_type_example_names_prefers_001_and_skips_dead_homologs():
    text = "\n".join([
        ">MER1 - some late enzyme (Org) [M01.021]#M01#{u}~s X~",
        ">MER2 - aminopeptidase N (Homo sapiens) [M01.001]#M01#{u}~s X~",
        ">MER3 - dead relative (Org) [C26.965]#C26#{u}~s X~",
        ">MER4 - putative thing (Org) [C26.A01]#C26#{u}~s X~",
        ">MER5 - subtilisin Carlsberg (Bacillus licheniformis) [S08.001]#S08A#{u}~s X~",
    ])
    names = type_example_names(text)
    assert names["M01"] == "aminopeptidase N"          # .001 wins, organism stripped
    assert names["C26"] == "putative thing"            # 9xx never names; putative fallback
    assert names["S08"] == "subtilisin Carlsberg"
    assert names["S08A"] == "subtilisin Carlsberg"     # subfamily named via its token


# ── node ids ─────────────────────────────────────────────────────────────────

def test_merops_node_id_registered_prefixes():
    assert _merops_node_id("SC") == "merops.clan:SC"
    assert _merops_node_id("S14") == "merops.family:S14"
    assert _merops_node_id("S08A") == "merops.family:S08A"


# ── fixtures ─────────────────────────────────────────────────────────────────

def _candidate(**over) -> dict:
    base = {
        "code": "S14", "family": "S14", "subfamily": None, "clan": "SK",
        "catalytic_type": "S", "entry_type": "peptidase",
        "level_kind": "merops_family", "tier": 3, "confidence_score": 0.258,
        "identity": 43.1, "qcov": 85.5, "scov": 81.2, "evalue": 6.09e-41,
        "length": 188, "consensus_n": 13, "consensus_agreement": "family",
        "best_hit_id": "S14.008", "best_hit_mernum": "MER0004857",
        "best_hit_kind": "holotype", "homolog_hit_fraction": 0.0,
    }
    base.update(over)
    return base


def _write_strain(tmp_path: Path, genes: dict, calls: dict, strain: str = "TESTSTRAIN") -> Path:
    genome_dir = tmp_path / strain
    (genome_dir / "merops").mkdir(parents=True)
    with open(genome_dir / "gene_annotations_merged.json", "w", encoding="utf-8") as fh:
        json.dump(genes, fh)
    with open(genome_dir / "merops" / f"{strain}.merops.calls.json", "w", encoding="utf-8") as fh:
        json.dump(calls, fh)
    return genome_dir


def _write_reference(tmp_path: Path) -> Path:
    ref = {
        "families": {
            "S14": {"name": "endopeptidase Clp", "clan": "SK", "family_type": "peptidase"},
            "S08": {"name": "subtilisin Carlsberg", "clan": "SB", "family_type": "peptidase"},
            "U32": {"name": "collagenase U32", "clan": None, "family_type": "peptidase"},
            "I39": {"name": "alpha-2-macroglobulin", "clan": "IL", "family_type": "inhibitor"},
        },
        "subfamily_names": {"S08A": "subtilisin Carlsberg"},
        "clans": {
            "SK": {"description": "Ser/Lys or Ser/His catalytic dyad", "family_type": "peptidase"},
            "SB": {"description": "Ser/His/Asp subtilisin fold", "family_type": "peptidase"},
            "IL": {"description": "Inhibitor clan IL", "family_type": "inhibitor"},
        },
    }
    path = tmp_path / "merops_reference.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(ref, fh)
    return path


def _write_genome_config(tmp_path: Path, genome_dirs: list[Path]) -> Path:
    path = tmp_path / "genomes.csv"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("strain_name,data_dir\n")
        for d in genome_dirs:
            fh.write(f"{d.name},{d}\n")
    return path


# ── per-strain gene edges ────────────────────────────────────────────────────

def test_gene_edge_carries_call_class_and_evidence(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]}}
    calls = {"WP_1.1": {"calls": [_candidate()]}}
    genome_dir = _write_strain(tmp_path, genes, calls)
    a = MeropsAnnotationAdapter(genome_dir)
    edges = list(a.get_edges())
    assert len(edges) == 1
    eid, src, tgt, label, props = edges[0]
    assert label == "gene_has_merops_family"
    assert src == "ncbigene:LT001"
    assert tgt == "merops.family:S14"
    assert props["call_class"] == "peptidase"
    assert props["tier"] == 3
    assert props["confidence_score"] == 0.258
    assert props["identity"] == 43.1
    assert props["qcov"] == 85.5
    assert props["evalue"] == 6.09e-41
    assert props["consensus_n"] == 13
    assert props["best_hit_id"] == "S14.008"
    assert props["best_hit_kind"] == "holotype"
    # deliberately dropped from the edge
    assert "scov" not in props and "length" not in props
    assert "homolog_hit_fraction" not in props and "consensus_agreement" not in props


def test_gene_edge_nonpeptidase_homolog_class(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["C26"]}}
    calls = {"WP_1.1": {"calls": [_candidate(
        code="C26", family="C26", clan="PC", best_hit_id="C26.965",
        best_hit_kind="nonpeptidase_homolog", homolog_hit_fraction=1.0,
    )]}}
    a = MeropsAnnotationAdapter(_write_strain(tmp_path, genes, calls))
    (_, _, _, _, props), = a.get_edges()
    assert props["call_class"] == "nonpeptidase_homolog"


def test_multi_domain_protein_fans_out(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S09", "S33"]}}
    calls = {"WP_1.1": {"calls": [
        _candidate(code="S09", family="S09", clan="SC"),
        _candidate(code="S33", family="S33", clan="SC"),
    ]}}
    a = MeropsAnnotationAdapter(_write_strain(tmp_path, genes, calls))
    targets = {e[2] for e in a.get_edges()}
    assert targets == {"merops.family:S09", "merops.family:S33"}


def test_tier1_call_attaches_at_family_with_identifier_on_edge(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["M75.001"]}}
    calls = {"WP_1.1": {"calls": [_candidate(
        code="M75.001", family="M75", subfamily=None, level_kind="merops_id",
        tier=1, best_hit_id="M75.001", identity=72.9,
    )]}}
    a = MeropsAnnotationAdapter(_write_strain(tmp_path, genes, calls))
    (_, _, tgt, _, props), = a.get_edges()
    assert tgt == "merops.family:M75"
    assert props["best_hit_id"] == "M75.001"   # identifier survives on the edge
    assert props["tier"] == 1


def test_no_calls_json_yields_no_edges(tmp_path):
    genome_dir = tmp_path / "EMPTY"
    genome_dir.mkdir()
    with open(genome_dir / "gene_annotations_merged.json", "w", encoding="utf-8") as fh:
        json.dump({"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]}}, fh)
    a = MeropsAnnotationAdapter(genome_dir)
    assert list(a.get_edges()) == []


# ── Multi orchestrator: nodes + hierarchy ────────────────────────────────────

@pytest.fixture
def multi(tmp_path):
    genes = {
        "LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]},
        "LT002": {"protein_id": "WP_2.1", "merops_ids": ["S08A"]},
        "LT003": {"protein_id": "WP_3.1", "merops_ids": ["U32"]},
        "LT004": {"protein_id": "WP_4.1", "merops_ids": ["I39"]},
    }
    calls = {
        "WP_1.1": {"calls": [_candidate()]},
        "WP_2.1": {"calls": [_candidate(
            code="S08A", family="S08", subfamily="S08A", clan="SB",
            level_kind="merops_subfamily", tier=2, identity=51.0,
        )]},
        "WP_3.1": {"calls": [_candidate(code="U32", family="U32", clan=None, catalytic_type="U")]},
        "WP_4.1": {"calls": [_candidate(
            code="I39", family="I39", clan="IL", catalytic_type=None,
            entry_type="inhibitor",
        )]},
    }
    genome_dir = _write_strain(tmp_path, genes, calls)
    config = _write_genome_config(tmp_path, [genome_dir])
    m = MultiMeropsAnnotationAdapter(
        genome_config_file=str(config),
        reference_path=_write_reference(tmp_path),
    )
    m.download_data()
    return m


def test_nodes_observed_plus_ancestors_only(multi):
    nodes = {nid: props for nid, _, props in multi.get_nodes()}
    assert set(nodes) == {
        "merops.clan:SK", "merops.clan:SB", "merops.clan:IL",
        "merops.family:S14", "merops.family:S08", "merops.family:S08A",
        "merops.family:U32", "merops.family:I39",
    }
    for props in nodes.values():
        assert props["level"] in (0, 1, 2)


def test_node_properties(multi):
    nodes = {nid: props for nid, _, props in multi.get_nodes()}
    s14 = nodes["merops.family:S14"]
    assert s14 == {
        "merops_id": "S14", "level": 1, "level_kind": "merops_family",
        "family_type": "peptidase", "catalytic_type": "serine",
        "name": "endopeptidase Clp",
    }
    clan = nodes["merops.clan:SK"]
    assert clan["level"] == 0 and clan["level_kind"] == "merops_clan"
    assert clan["name"] == "SK"
    assert clan["description"] == "Ser/Lys or Ser/His catalytic dyad"
    subfam = nodes["merops.family:S08A"]
    assert subfam["level"] == 2 and subfam["name"] == "subtilisin Carlsberg"
    inhib = nodes["merops.family:I39"]
    assert inhib["family_type"] == "inhibitor"
    assert "catalytic_type" not in inhib  # sparse — absent, not empty/None
    # clan-unassigned family stays a root but keeps its level
    assert nodes["merops.family:U32"]["level"] == 1


def test_hierarchy_edges_child_to_parent(multi):
    isa = {(e[1], e[2]) for e in multi.get_edges()
           if e[3] == "merops_family_is_a_merops_family"}
    assert isa == {
        ("merops.family:S14", "merops.clan:SK"),
        ("merops.family:S08", "merops.clan:SB"),
        ("merops.family:S08A", "merops.family:S08"),
        ("merops.family:I39", "merops.clan:IL"),
        # U32 has no clan → root, no parent edge
    }


def test_all_gene_edge_targets_are_nodes(multi):
    node_ids = {nid for nid, _, _ in multi.get_nodes()}
    gene_edges = [e for e in multi.get_edges() if e[3] == "gene_has_merops_family"]
    assert len(gene_edges) == 4
    assert {e[2] for e in gene_edges} <= node_ids


def test_inhibitor_edge_call_class(multi):
    edges = {e[2]: e[4] for e in multi.get_edges() if e[3] == "gene_has_merops_family"}
    assert edges["merops.family:I39"]["call_class"] == "inhibitor"


def test_missing_reference_fails_loudly(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]}}
    calls = {"WP_1.1": {"calls": [_candidate()]}}
    genome_dir = _write_strain(tmp_path, genes, calls)
    config = _write_genome_config(tmp_path, [genome_dir])
    m = MultiMeropsAnnotationAdapter(
        genome_config_file=str(config),
        reference_path=tmp_path / "nope.json",
    )
    with pytest.raises(FileNotFoundError, match="prepare_data"):
        m.download_data()
