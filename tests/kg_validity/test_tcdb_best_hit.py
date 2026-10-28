"""KG validity tests for TCDB diamond best hits after the TCDB cleanup.

Best hits are their own `Gene_has_tcdb_family` edges (`tcdb_call = 'diamond_best_hit'`,
one per best-hit system, `rank` 1 = top); the six `best_hit_*` properties are gone and
best-hit substrates are the system node's own substrate edges.

Specs: docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md (selection
rule best_hit_v1) and docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md §5.3/§9.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
GENOMES_CSV = REPO / "data" / "Prochlorococcus" / "genomes" / "cyanobacteria_genomes.csv"
TCDB_PRUNED = REPO / "cache/data/tcdb/tcdb_pruned.json"


@lru_cache(maxsize=1)
def _expected_edges_by_call() -> dict[str, int]:
    """Gene→TCDB edges the merge inputs predict, per tcdb_call (spec §9).

    Mirrors the adapter: eggNOG ids from gene_annotations_merged.json, diamond tier
    calls and best-hit systems (incl. near-ties and call-less proteins) from each
    strain's calls.json; each TC id lands on itself when kept, else on its
    seed_aliases anchor when that is kept, else is dropped; within one call kind,
    duplicates of a (gene, landed node) collapse to one edge.
    """
    if not TCDB_PRUNED.exists():
        pytest.skip(f"{TCDB_PRUNED.relative_to(REPO)} missing; run prepare_data step 6")
    pruned = json.loads(TCDB_PRUNED.read_text())
    kept = set(pruned["kept_tcdb_ids"])
    aliases = pruned.get("seed_aliases", {})

    def land(tc):
        if not tc:
            return None
        if tc in kept:
            return tc
        anchor = aliases.get(tc)
        return anchor if anchor in kept else None

    keys: set[tuple[str, str, str]] = set()
    lines = [ln for ln in GENOMES_CSV.read_text().splitlines() if not ln.lstrip().startswith("#")]
    for row in csv.DictReader(lines):
        d = REPO / (row.get("data_dir") or "").strip()
        merged_path = d / "gene_annotations_merged.json"
        if not merged_path.exists():
            continue
        calls_path = d / "tcdb" / f"{d.name}.tcdb.calls.json"
        calls = json.loads(calls_path.read_text()) if calls_path.exists() else {}
        for lt, gene in json.loads(merged_path.read_text()).items():
            for tc in gene.get("tcdb_eggnog_ids") or []:
                if (t := land(tc)):
                    keys.add((lt, t, "eggnog"))
            rec = calls.get((gene.get("protein_id") or "").strip()) or {}
            for cand in rec.get("calls") or []:
                if (t := land(cand.get("tcid"))):
                    keys.add((lt, t, "diamond_tier"))
            best = rec.get("best_hit")
            for hit in (best.get("hits") or []) if isinstance(best, dict) else []:
                if (t := land(hit.get("tcid"))):
                    keys.add((lt, t, "diamond_best_hit"))
    return dict(Counter(k[2] for k in keys))


@pytest.mark.kg
def test_edge_counts_by_call_match_merge_inputs(run_query):
    """Exact: graph edges per tcdb_call == the merge-input prediction (2026-10-23:
    eggnog 17,069 / diamond_tier 41,242 / diamond_best_hit 17,512)."""
    expected = _expected_edges_by_call()
    live = {r["k"]: r["n"] for r in run_query(
        "MATCH ()-[r:Gene_has_tcdb_family]->() RETURN r.tcdb_call AS k, count(*) AS n")}
    print(f"expected {expected}; live {live}")
    assert live == expected


@pytest.mark.kg
def test_best_hits_at_most_three_per_gene(run_query):
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->()
        WITH g, count(r) AS k WHERE k > 3
        RETURN count(g) AS n""")[0]["n"]
    assert n == 0


@pytest.mark.kg
def test_best_hits_target_systems_of_classes_1_to_4(run_query):
    """best_hit_v1 keeps TC classes 1-4 only; a best hit lands on a tc_system unless
    a retired id was re-anchored via seed_aliases (then it is not direct)."""
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->(t:TcdbFamily)
        RETURN count(r) AS n,
               sum(CASE WHEN NOT split(t.tcdb_id, '.')[0] IN ['1', '2', '3', '4'] THEN 1 ELSE 0 END) AS bad_class,
               sum(CASE WHEN t.level_kind <> 'tc_system' AND r.tcdb_gate = 'direct' THEN 1 ELSE 0 END) AS bad_direct,
               sum(CASE WHEN t.level_kind = 'tc_system' THEN 1 ELSE 0 END) AS on_system""")[0]
    assert row["n"] > 0 and row["bad_class"] == 0 and row["bad_direct"] == 0, row
    assert row["on_system"] >= 0.99 * row["n"], row


@pytest.mark.kg
def test_call_less_proteins_get_best_hit_edges(run_query):
    """The ~188 proteins with a best hit but no tier call (top HSP < 50 aa) now
    carry a diamond_best_hit edge (they had none). 698 over 117 strains
    (2026-10-28, Ford 2021's 74 Cyanorak genomes)."""
    n = run_query("""
        MATCH (g:Gene) WHERE EXISTS { (g)-[:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->() }
          AND NOT EXISTS { (g)-[:Gene_has_tcdb_family {tcdb_call: 'diamond_tier'}]->() }
        RETURN count(g) AS n""")[0]["n"]
    assert 550 <= n <= 850, f"{n} genes with best hits but no tier call (expected ~700)"


@pytest.mark.kg
@pytest.mark.parametrize("locus_tag,tcid,substrate", [
    ("PMM0370", "3.A.1.16.2", "cyanate"),   # cynA
    ("PMM0372", "3.A.1.16.2", "cyanate"),   # cynD — tier call is only 3.A.1
    ("PMM0263", "1.A.11.2.5", "ammon"),     # amt1
    ("PMM0970", "3.A.1.4.4", "urea"),       # urtA
    ("PMM0973", "3.A.1.4.4", "urea"),       # urtD
])
def test_med4_nitrogen_best_hits(run_query, locus_tag, tcid, substrate):
    """The rank-1 best hit names the system, whose OWN substrate edges carry the
    nitrogen source (no copy on the gene edge)."""
    rows = run_query("""
        MATCH (g:Gene {locus_tag: $lt})-[r:Gene_has_tcdb_family {tcdb_call: 'diamond_best_hit'}]->(t:TcdbFamily)
        WHERE r.rank = 1
        RETURN t.tcdb_id AS tc, t.name AS name,
               [(t)-[:Tcdb_family_transports_metabolite|Tcdb_family_transports_chemical_role]->(m) | m.name] AS subs""",
                     lt=locus_tag)
    assert len(rows) == 1, rows
    r = rows[0]
    assert r["tc"] == tcid
    assert substrate in r["name"].lower(), r["name"]   # a real, untruncated TCDB name
    assert any(substrate in (s or "").lower() for s in r["subs"]), r["subs"]
