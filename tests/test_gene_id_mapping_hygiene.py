"""Gene-ID mapping hygiene bundle (plans/gene_id_mapping_hygiene.md, 2026-08-28).

1. GFF ``Name`` typing — a gene symbol / protein accession in ``Name`` must not
   become a Tier-1 locus tag (KT2440 tnpB x12 merge, MED4 rrf/frr merge).
2. Runaway guard — build_diagnostic_report flags genes hoarding Tier-1 ids.
3. ``gene-`` prefix heuristic — ``gene-PMM0950`` resolves through ``PMM0950``.
"""
import pytest

from multiomics_kg.download import build_gene_id_mapping as bgm
from multiomics_kg.download.gene_id_graph import (
    RUNAWAY_TIER1_MIN,
    GeneIdGraph,
)
from multiomics_kg.utils.gene_id_utils import _heuristic_candidates


# ── 1. GFF Name typing ───────────────────────────────────────────────────────

def _gff(tmp_path, rows):
    p = tmp_path / "a.gff"
    p.write_text("##gff-version 3\n" + "".join(
        f"chr\tRefSeq\t{feat}\t1\t2\t.\t+\t.\t{attrs}\n" for feat, attrs in rows))
    return p


def _rows_for(tmp_path, monkeypatch, rows):
    monkeypatch.setattr(bgm, "PROJECT_ROOT", tmp_path)
    path = _gff(tmp_path, rows)
    return bgm.extract_rows_from_annotation_gff({"filename": path.name}, "paper", "gff")


def test_gff_name_equal_to_gene_symbol_is_tier3(tmp_path, monkeypatch):
    rows = _rows_for(tmp_path, monkeypatch, [
        ("gene", "ID=gene-PP_RS16230;Name=tnpB;gene=tnpB;locus_tag=PP_RS16230;old_locus_tag=PP_3113"),
    ])
    pairs = dict(rows[0][0])
    assert pairs["tnpB"] == "gene_name"
    assert pairs["PP_RS16230"] == "locus_tag_ncbi"
    assert pairs["PP_3113"] == "old_locus_tag"


def test_gff_name_equal_to_protein_id_is_tier2(tmp_path, monkeypatch):
    rows = _rows_for(tmp_path, monkeypatch, [
        ("CDS", "ID=cds-WP_1.1;Name=WP_1.1;gene=dnaN;locus_tag=TX50_RS00020;protein_id=WP_1.1"),
    ])
    pairs = rows[0][0]
    assert ("WP_1.1", "protein_id_refseq") in pairs
    assert ("WP_1.1", "locus_tag_ncbi") not in pairs


def test_gff_name_equal_to_locus_tag_not_duplicated(tmp_path, monkeypatch):
    rows = _rows_for(tmp_path, monkeypatch, [
        ("gene", "ID=gene-PMM0001;Name=PMM0001;locus_tag=PMM0001;old_locus_tag=PMED4_00011"),
    ])
    pairs = rows[0][0]
    assert pairs.count(("PMM0001", "locus_tag_ncbi")) == 1


def test_gff_name_that_is_neither_stays_tier1(tmp_path, monkeypatch):
    rows = _rows_for(tmp_path, monkeypatch, [
        ("gene", "ID=gene-X;Name=P9301_05911;locus_tag=P9301_RS02890"),
    ])
    assert ("P9301_05911", "locus_tag_ncbi") in rows[0][0]


def test_shared_symbol_no_longer_merges_genes(tmp_path, monkeypatch):
    """Two IS copies with Name=tnpB stay two genes (the KT2440 M8001_03425 case)."""
    rows = _rows_for(tmp_path, monkeypatch, [
        ("gene", "ID=gene-A;Name=tnpB;gene=tnpB;locus_tag=A_RS1;old_locus_tag=A_1"),
        ("gene", "ID=gene-B;Name=tnpB;gene=tnpB;locus_tag=B_RS1;old_locus_tag=B_1"),
    ])
    g = GeneIdGraph()
    g.add_anchor("A_1")
    g.add_anchor("B_1")
    g.process_all_rows(rows)
    assert g.specific_lookup["A_RS1"] == "A_1"
    assert g.specific_lookup["B_RS1"] == "B_1"
    assert sorted(g.multi_lookup.get("tnpB", [])) == ["A_1", "B_1"]
    assert not g.conflicts


# ── 2. Runaway guard ─────────────────────────────────────────────────────────

def _graph_with_counts(counts):
    g = GeneIdGraph()
    for i, n in enumerate(counts):
        lt = f"G{i:04d}"
        g.add_anchor(lt)
        for k in range(n):
            g.add_id_for_gene(lt, f"{lt}_alt{k}", "old_locus_tag", "src")
    return g


def test_runaway_guard_flags_hoarder():
    g = _graph_with_counts([6] * 20 + [40])
    rep = g.build_diagnostic_report()
    assert rep["tier1_count_median"] == 6
    assert rep["tier1_count_max"] == 40
    assert [r["locus_tag"] for r in rep["runaway_genes"]] == ["G0020"]
    assert any(w.startswith("[RUNAWAY]") and "G0020 (40)" in w for w in rep["warnings"])


def test_runaway_guard_quiet_on_normal_spread():
    g = _graph_with_counts([6] * 20 + [9, 10])
    rep = g.build_diagnostic_report()
    assert rep["runaway_genes"] == []
    assert not any(w.startswith("[RUNAWAY]") for w in rep["warnings"])


def test_runaway_guard_floor_protects_small_medians():
    """3x a median of 3 is 9 — below the floor, so 11 ids is not a runaway."""
    g = _graph_with_counts([3] * 20 + [RUNAWAY_TIER1_MIN - 1])
    assert g.build_diagnostic_report()["runaway_genes"] == []


# ── 3. gene- prefix heuristic ────────────────────────────────────────────────

@pytest.mark.parametrize("raw,expected", [
    ("gene-PMM0950", "PMM0950"),
    ("cds-WP_011132156.1", "WP_011132156.1"),
    ("rna-TX50_RS01725", "TX50_RS01725"),
    ("gene-PMED4_ncRNA_Yfr11", "PMED4_ncRNA_Yfr11"),
])
def test_gff_prefix_stripped(raw, expected):
    assert expected in _heuristic_candidates(raw)


def test_no_prefix_no_candidate():
    assert "PMM0950" not in _heuristic_candidates("PMM0950")
    assert _heuristic_candidates("gene-") == []


# ── 4. Pass 2b: heuristic candidates against multi_lookup singletons ─────────

def test_heuristic_candidate_resolves_via_multi_singleton():
    from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
    md = MappingData(specific_lookup={}, multi_lookup={"AAV93747.1": ["SPO0429"]},
                     conflicts={}, locus_tags={"SPO0429"})
    lt, how = resolve_row({"NCBI reference": "AAV93747"}, "NCBI reference",
                          [{"column": "NCBI reference", "id_type": "protein_id_refseq"}], md)
    assert (lt, how) == ("SPO0429", "heuristic_multi:NCBI reference")


def test_heuristic_candidate_ambiguous_in_multi_is_not_resolved():
    from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
    md = MappingData(specific_lookup={}, multi_lookup={"AHB86005.1": ["SPO0344a", "SPO0344b"]},
                     conflicts={}, locus_tags={"SPO0344a", "SPO0344b"})
    lt, how = resolve_row({"c": "AHB86005"}, "c", [], md)
    assert lt is None and how == "ambiguous"


def test_gene_prefixed_id_resolves_through_bare_tag():
    from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
    md = MappingData(specific_lookup={"PMM0950": "PMM0950"}, multi_lookup={},
                     conflicts={"gene-PMM0950": ["PMM0950", "PMM0236"]}, locus_tags={"PMM0950"})
    lt, how = resolve_row({"GeneID": "gene-PMM0950"}, "GeneID", [], md)
    assert (lt, how) == ("PMM0950", "heuristic:GeneID")


# ── 5. Pass 3a: unique gene_name claimant beats synonym-only claimants ───────

def test_named_claimant_wins_ambiguous_symbol():
    from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
    md = MappingData(specific_lookup={}, multi_lookup={"atpB": ["A9601_16401", "A9601_16581"]},
                     conflicts={}, locus_tags={"A9601_16401", "A9601_16581"},
                     named_lookup={"atpB": ["A9601_16401"]})
    assert resolve_row({"g": "atpB"}, "g", [], md) == ("A9601_16401", "multi_named:g")


def test_two_named_claimants_stay_ambiguous():
    from multiomics_kg.utils.gene_id_utils import MappingData, resolve_row
    md = MappingData(specific_lookup={}, multi_lookup={"psbA": ["A", "B"]}, conflicts={},
                     locus_tags={"A", "B"}, named_lookup={"psbA": ["A", "B"]})
    assert resolve_row({"g": "psbA"}, "g", [], md) == (None, "ambiguous")


def test_load_mapping_v2_builds_named_lookup(tmp_path):
    import json
    from multiomics_kg.utils.gene_id_utils import load_mapping_v2
    (tmp_path / "gene_id_mapping.json").write_text(json.dumps({
        "version": 2, "specific_lookup": {}, "conflicts": {},
        "multi_lookup": {"atpB": ["G1", "G2"], "rbcL": ["G3"]},
        "genes": {
            "G1": {"tier1_ids": [], "tier2_ids": [], "tier3_ids": [{"id": "atpB", "type": "gene_name"}]},
            "G2": {"tier1_ids": [], "tier2_ids": [], "tier3_ids": [{"id": "atpB", "type": "gene_synonym"}]},
            "G3": {"tier1_ids": [], "tier2_ids": [], "tier3_ids": [{"id": "rbcL", "type": "gene_name"}]},
        }}))
    md = load_mapping_v2(tmp_path)
    assert md.named_lookup == {"atpB": ["G1"]}   # singleton rbcL needs no tie-break entry


# ── 6. A paper column cannot promote an annotation-known gene symbol to Tier 1 ─

def test_known_gene_name_typed_tier1_by_paper_is_demoted():
    g = GeneIdGraph()
    g.add_anchor("G1"); g.add_anchor("G2")
    g.add_id_for_gene("G1", "rplF", "gene_name", "annotation")       # seeding: rplF is a symbol
    g.add_id_for_gene("G2", "Q318J8", "uniprot_accession", "annotation")
    # biller-2022-shaped row: "Gene Number" typed locus_tag holds the symbol,
    # anchored through the accession onto G2
    g.process_all_rows([([("rplF", "locus_tag"), ("Q318J8", "uniprot_accession")], "paper/table")])
    assert "rplF" not in g.specific_lookup                       # never Tier 1
    assert sorted(g.multi_lookup["rplF"]) == ["G1", "G2"]        # paper claim kept as Tier 3
    assert not g.conflicts
    rep = g.build_diagnostic_report()
    assert rep["tier1_demoted_known_names"] == 1


def test_unknown_token_typed_tier1_still_promoted():
    g = GeneIdGraph()
    g.add_anchor("G1")
    g.process_all_rows([([("P9301_05911", "old_locus_tag"), ("G1", "locus_tag")], "paper/table")])
    assert g.specific_lookup["P9301_05911"] == "G1"
