"""KG validity tests for TCDB diamond best-hit details (best_hit_v1).

Spec: docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
GENOMES_CSV = REPO / "data" / "Prochlorococcus" / "genomes" / "cyanobacteria_genomes.csv"

BH_PROPS = ("best_hit_tcids", "best_hit_names", "best_hit_identity",
            "best_hit_substrate_tcids", "best_hit_substrate_ids", "best_hit_substrate_names")


def _expected_edge_count() -> tuple[int, int]:
    """(genes expected to carry the record, attach edges the adapter prune drops).

    A gene carries it when its protein has a calls.json best_hit with attach_tcid and
    that TC id is one of the gene's diamond ids; the edge survives when its target is a
    kept TC id or seed-aliased onto one.
    """
    pruned = json.loads((REPO / "cache/data/tcdb/tcdb_pruned.json").read_text())
    kept = set(pruned["kept_tcdb_ids"])
    aliases = pruned.get("seed_aliases", {})
    lines = [ln for ln in GENOMES_CSV.read_text().splitlines() if not ln.lstrip().startswith("#")]
    expected = dropped = 0
    for row in csv.DictReader(lines):
        d = REPO / (row.get("data_dir") or "").strip()
        calls_path = d / "tcdb" / f"{d.name}.tcdb.calls.json"
        merged_path = d / "gene_annotations_merged.json"
        if not calls_path.exists() or not merged_path.exists():
            continue
        calls = json.loads(calls_path.read_text())
        for gene in json.loads(merged_path.read_text()).values():
            best = calls.get((gene.get("protein_id") or "").strip(), {}).get("best_hit") or {}
            tc = best.get("attach_tcid")
            if not tc or tc not in (gene.get("tcdb_diamond_ids") or []):
                continue
            if tc in kept or aliases.get(tc) in kept:
                expected += 1
            else:
                dropped += 1
    return expected, dropped


@pytest.mark.kg
def test_best_hit_props_co_occur_and_align(run_query):
    row = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WHERE r.best_hit_tcids IS NOT NULL OR r.best_hit_names IS NOT NULL
           OR r.best_hit_identity IS NOT NULL OR r.best_hit_substrate_ids IS NOT NULL
        RETURN count(r) AS n,
          sum(CASE WHEN r.best_hit_tcids IS NULL OR r.best_hit_names IS NULL
                     OR r.best_hit_identity IS NULL
                     OR size(r.best_hit_names) <> size(r.best_hit_tcids) THEN 1 ELSE 0 END) AS bad_head,
          sum(CASE WHEN r.best_hit_substrate_ids IS NOT NULL AND (
                     r.best_hit_substrate_tcids IS NULL OR r.best_hit_substrate_names IS NULL
                     OR size(r.best_hit_substrate_tcids) <> size(r.best_hit_substrate_ids)
                     OR size(r.best_hit_substrate_names) <> size(r.best_hit_substrate_ids)
                     OR any(t IN r.best_hit_substrate_tcids WHERE NOT t IN r.best_hit_tcids))
                   THEN 1 ELSE 0 END) AS bad_subs,
          sum(CASE WHEN size(r.best_hit_tcids) > 3 THEN 1 ELSE 0 END) AS too_many
    """)[0]
    assert row["n"] > 0
    assert row["bad_head"] == 0 and row["bad_subs"] == 0 and row["too_many"] == 0, row


@pytest.mark.kg
def test_best_hit_only_on_diamond_edges(run_query):
    n = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->()
        WHERE r.best_hit_tcids IS NOT NULL AND NOT 'tcdb_diamond' IN r.sources
        RETURN count(r) AS n""")[0]["n"]
    assert n == 0


@pytest.mark.kg
def test_best_hit_on_at_most_one_edge_per_gene(run_query):
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
        WHERE r.best_hit_tcids IS NOT NULL
        WITH g, count(r) AS k WHERE k > 1
        RETURN count(g) AS n""")[0]["n"]
    assert n == 0


@pytest.mark.kg
def test_best_hit_edge_count_matches_artifacts(run_query):
    expected, dropped = _expected_edge_count()
    print(f"expected best-hit edges {expected}; attach edges dropped by the prune {dropped}")
    n = run_query("""
        MATCH ()-[r:Gene_has_tcdb_family]->() WHERE r.best_hit_tcids IS NOT NULL
        RETURN count(r) AS n""")[0]["n"]
    assert n == expected


@pytest.mark.kg
@pytest.mark.parametrize("locus_tag,tcid,substrate", [
    ("PMM0370", "3.A.1.16.2", "cyanate"),   # cynA
    ("PMM0372", "3.A.1.16.2", "cyanate"),   # cynD — tier call is only 3.A.1
    ("PMM0263", "1.A.11.2.5", "ammon"),     # amt1
    ("PMM0970", "3.A.1.4.4", "urea"),       # urtA
    ("PMM0973", "3.A.1.4.4", "urea"),       # urtD
])
def test_med4_nitrogen_best_hits(run_query, locus_tag, tcid, substrate):
    rows = run_query("""
        MATCH (g:Gene {locus_tag: $lt})-[r:Gene_has_tcdb_family]->()
        WHERE r.best_hit_tcids IS NOT NULL
        RETURN r.best_hit_tcids AS tcids, r.best_hit_names AS names,
               r.best_hit_substrate_names AS subs""", lt=locus_tag)
    assert len(rows) == 1, rows
    r = rows[0]
    assert r["tcids"][0] == tcid
    assert r["names"][0] != tcid                         # a real TCDB name, not the fallback
    assert any(substrate in s.lower() for s in r["subs"] or []), r["subs"]
