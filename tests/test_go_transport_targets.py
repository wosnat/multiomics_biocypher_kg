from tests.metabolite_identity_helpers import make_chebi, term
from multiomics_kg.utils.metabolite_identity.go_transport import (
    go_substrate_strings, assemble_go_targets, rewrite_go_targets, target_distribution)


def test_go_substrate_strings_primary_and_skipped():
    chebi = make_chebi(CHEBI_29195=term("cyanate", alt=["CHEBI:11111"]))
    out, skipped = go_substrate_strings(
        {"GO:1": ["CHEBI:11111", "CHEBI:99999"], "GO:9": ["CHEBI:29195"]}, {"GO:1"}, chebi)
    assert out == {"GO:1": [("CHEBI:11111", "CHEBI:29195;cyanate")]}
    assert skipped == ["CHEBI:99999"]


def test_kg_term_without_go_plus_entry_is_skipped():
    chebi = make_chebi(CHEBI_29195=term("cyanate"))
    out, skipped = go_substrate_strings({"GO:1": ["CHEBI:29195"]}, {"GO:1", "GO:0000000"}, chebi)
    assert "GO:0000000" not in out and skipped == []


def _ts(*pairs):
    return {"GO:1": list(pairs)}


def test_assemble_drops_filtered_and_keeps_k():
    ts = _ts(("CHEBI:15378", "CHEBI:15378;H+"), ("CHEBI:29103", "CHEBI:29103;K+"),
             ("CHEBI:29195", "CHEBI:29195;cyanate"))
    s2p = {"CHEBI:15378;H+": "kegg.compound:C00080", "CHEBI:29103;K+": "kegg.compound:C00238",
           "CHEBI:29195;cyanate": "kegg.compound:C01417"}
    subs, roles, drops = assemble_go_targets(ts, s2p, {}, set(), {"kegg.compound:C00080"})
    assert set(subs["GO:1"]) == {"kegg.compound:C00238", "kegg.compound:C01417"}
    assert drops["filtered"] == 1 and not roles


def test_role_without_node_dropped():
    ts = {"GO:1": [("CHEBI:1", "CHEBI:1;a"), ("CHEBI:2", "CHEBI:2;b")]}
    rm = {"CHEBI:1;a": "chebi.role:99999", "CHEBI:2;b": "chebi.role:26672"}
    subs, roles, drops = assemble_go_targets(ts, {}, rm, {"chebi.role:26672"}, set())
    assert roles == {"GO:1": {"chebi.role:26672": ["CHEBI:2"]}}
    assert drops["role_without_node"] == 1 and not subs


def test_unresolved():
    ts = {"GO:1": [("CHEBI:1", "CHEBI:1;a")]}
    subs, roles, drops = assemble_go_targets(ts, {}, {}, set(), set())
    assert drops["unresolved"] == 1 and not subs


def test_two_axioms_one_node_one_target():
    ts = _ts(("CHEBI:29195", "CHEBI:29195;cyanate"), ("CHEBI:28024", "CHEBI:28024;cyanic acid"))
    s2p = {"CHEBI:29195;cyanate": "kegg.compound:C01417",
           "CHEBI:28024;cyanic acid": "kegg.compound:C01417"}
    subs, _, _ = assemble_go_targets(ts, s2p, {}, set(), set())
    assert subs["GO:1"] == {"kegg.compound:C01417": ["CHEBI:28024", "CHEBI:29195"]}


def test_rewrite_merges_and_post_filter():
    gs = {"GO:1": {"A": ["CHEBI:2"], "B": ["CHEBI:1"], "X": ["CHEBI:3"]}}
    out, n = rewrite_go_targets(gs, {"B": "A", "X": "Z"}, {"Z"})
    assert out == {"GO:1": {"A": ["CHEBI:1", "CHEBI:2"]}} and n == 1


def test_target_distribution():
    targets = {"GO:1": {"T1": ["c"]}, "GO:2": {"T1": ["c"], "T2": ["c"]}, "GO:3": {"T2": ["c"]}}
    genes = {"s1": [{"GO:1"}, {"GO:3"}], "s2": [{"GO:2"}, {"GO:9"}]}
    rows = target_distribution(genes, targets)
    assert rows == [
        {"target": "T1", "genes": 2, "gene_share": 0.6667, "terms": 2, "term_share": 0.6667, "strains": 2},
        {"target": "T2", "genes": 2, "gene_share": 0.6667, "terms": 2, "term_share": 0.6667, "strains": 2},
    ]
