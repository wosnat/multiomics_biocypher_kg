"""Validation report for TCDB diamond best-hit details (best_hit_v1).

Spec: docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md §6.
The rule is frozen; nothing here tunes it. MED4 is training data (the rule was measured
on it); every other strain is hold-out. Agreement with the comparators is an upper
bound — they share homology with the TCDB calls.

Part A reads the committed artifacts (calls.json, .tcdb.tsv, merged annotations).
Part B queries the live graph (bolt://localhost:7687) after a rebuild.

Run:  uv run python scripts/tcdb_best_hit_validation.py [--artifacts-only] > report.md
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from multiomics_kg.utils.tcdb_diamond import parse_diamond_row, parse_tcdb_subject_id  # noqa: E402

GENOMES_CSV = REPO / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
TRAINING = {"MED4"}
NITROGEN_STRAINS = ["MED4", "MIT9313", "WH8102", "HOT1A3"]
N_KEYWORDS = ("ammon", "urea", "cyanate", "nitrate", "nitrite", "amino acid", "peptide",
              "glutamate", "glutamine", "arginine", "spermidine", "putrescine", "purine",
              "pyrimidine", "nucleoside", "taurine", "glycine betaine")
# Fixed before running (spec §6): product keywords for reciprocal-check failures.
FAIL_CATEGORIES = [
    ("chaperone", ("chaperone", "heat shock", "dnak", "clpb", "clpc", "groe")),
    ("protease", ("protease", "peptidase", "ftsh")),
    ("photosynthesis/respiration", ("photosystem", "chlorophyll", "cytochrome b6", "nad(p)h",
                                    "nadh dehydrogenase", "ndh")),
    ("transporter-like", ("transport", "permease", "abc", "channel", "porter", "efflux")),
]


def _strains():
    lines = [ln for ln in GENOMES_CSV.read_text().splitlines() if not ln.lstrip().startswith("#")]
    for row in csv.DictReader(lines):
        d = REPO / (row.get("data_dir") or "").strip()
        if (d / "tcdb" / f"{d.name}.tcdb.calls.json").exists():
            yield row["strain_name"], d


def _category(product: str) -> str:
    p = (product or "").lower()
    for name, kws in FAIL_CATEGORIES:
        if any(k in p for k in kws):
            return name
    return "other"


def part_a() -> None:
    print("## A. Artifacts (calls.json / .tcdb.tsv / merged annotations)\n")
    print("| Strain | set | proteins with call | best hit | no call | split (>1 system) | "
          "reciprocal fail |")
    print("|---|---|---|---|---|---|---|")
    fail_cats: dict[str, Counter] = {"training": Counter(), "hold-out": Counter()}
    tot = Counter()
    for strain, d in _strains():
        calls = json.loads((d / "tcdb" / f"{d.name}.tcdb.calls.json").read_text())
        merged = json.loads((d / "gene_annotations_merged.json").read_text())
        product = {g.get("protein_id"): g.get("product") or "" for g in merged.values()}
        # gate-passing proteins (same gate as select_best_hits) minus best hits = reciprocal fails
        gated = set()
        for line in open(d / "tcdb" / f"{d.name}.tcdb.tsv"):
            r = parse_diamond_row(line)
            if not r or r["qcov"] < 75 or r["scov"] < 75:
                continue
            p = parse_tcdb_subject_id(r["subject_id"], truncate_to=5)
            if p and p[1][0] in "1234" and len(p[1].split(".")) == 5:
                gated.add(r["query_id"])
        best = {q: r["best_hit"] for q, r in calls.items() if "best_hit" in r}
        with_call = sum(1 for r in calls.values() if r.get("calls"))
        no_call = sum(1 for q in best if not calls[q].get("calls"))
        split = sum(1 for b in best.values() if len(b["tcids"]) > 1)
        fails = gated - set(best)
        bucket = "training" if strain in TRAINING else "hold-out"
        for q in fails:
            fail_cats[bucket][_category(product.get(q, ""))] += 1
        print(f"| {strain} | {bucket} | {with_call} | {len(best)} | {no_call} | {split} | {len(fails)} |")
        tot.update(best=len(best), no_call=no_call, split=split, fails=len(fails))
    print(f"\nTotals: {tot['best']} best hits, {tot['no_call']} without a call, "
          f"{tot['split']} split, {tot['fails']} reciprocal-check failures.\n")
    print("Reciprocal-check failures by product keyword (list fixed before running):\n")
    for bucket, c in fail_cats.items():
        print(f"- {bucket}: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))
    print()


def _q(session, cypher, **kw):
    return session.run(cypher, **kw).data()


def _ident_bin(x: float) -> str:
    return "<40" if x < 40 else ("40-70" if x < 70 else ">=70")


def part_b() -> None:
    from neo4j import GraphDatabase
    print("## B. Live graph\n")
    # slice 1's GO target filter, as canonical node ids (written by step 6)
    filter_ids = set(json.loads((REPO / "cache/data/kegg/kegg_data.json").read_text())
                     .get("go_transport_excluded_targets", []))
    drv = GraphDatabase.driver("bolt://localhost:7687", auth=None)
    with drv.session() as s:
        rows = _q(s, """
            MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
            WHERE r.best_hit_tcids IS NOT NULL
            OPTIONAL MATCH (g)-[e:Gene_has_tcdb_family]->(et:TcdbFamily) WHERE 'eggnog' IN e.sources
            WITH g, r, t, collect(et.tcdb_id) AS egg
            OPTIONAL MATCH (g)-[:Gene_has_kegg_ko]->(:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory)
            WITH g, r, t, egg, count(CASE WHEN b.tree_code = 'ko02000' THEN 1 END) > 0 AS brite_transporter
            OPTIONAL MATCH (g)-[:Gene_has_cyanorak_role]->(cr:CyanorakRole)
            WITH g, r, t, egg, brite_transporter, collect(toLower(cr.name)) AS roles
            OPTIONAL MATCH (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function]->()
                           -[:Molecular_function_transports_metabolite|Biological_process_transports_metabolite]->(m)
            RETURN g.locus_tag AS lt, g.organism_name AS org, g.gene_name AS gene, g.product AS product,
                   r.best_hit_tcids AS tcids, r.best_hit_names AS names, r.best_hit_identity AS ident,
                   r.best_hit_substrate_ids AS sub_ids, r.best_hit_substrate_names AS sub_names,
                   t.tcdb_id AS edge_tc, egg, brite_transporter, roles,
                   collect(DISTINCT m.id) AS go_subs""")
        pf = _q(s, """
            MATCH (g:Gene)-[r:Gene_has_tcdb_family]->() WHERE r.best_hit_tcids IS NOT NULL
            WITH g, r, split(r.best_hit_tcids[0], '.') AS p
            WITH g, p[0] + '.' + p[1] + '.' + p[2] AS fam
            OPTIONAL MATCH (f:TcdbFamily {tcdb_id: fam})-[:Tcdb_family_has_pfam_domain]->(pfam:Pfam)
            WITH g, fam, collect(pfam) AS bridged
            OPTIONAL MATCH (g)-[:Gene_has_pfam]->(gp:Pfam)
            WITH g, bridged, collect(gp) AS own
            RETURN g.locus_tag AS lt,
                   CASE WHEN size(bridged) = 0 THEN 'no_bridge'
                        WHEN any(x IN own WHERE x IN bridged) THEN 'corroborated'
                        ELSE 'uncorroborated' END AS pfam""")
        node_ids = {r["id"] for r in _q(s, "MATCH (n) WHERE n:Metabolite OR n:ChemicalRole RETURN n.id AS id")}
        rollup = {r["tc"]: r["n"] for r in _q(s, """
            MATCH (t:TcdbFamily)-[x:Tcdb_family_transports_metabolite]->() RETURN t.tcdb_id AS tc, count(x) AS n""")}
        orgs = {r["org"]: r["strain"] for r in _q(s, """
            MATCH (o:OrganismTaxon) RETURN o.preferred_name AS org, o.strain_name AS strain""")}
    pfam_of = {r["lt"]: r["pfam"] for r in pf}

    def bucket(r):
        return "training" if any(k in (r["org"] or "") for k in TRAINING) else "hold-out"

    stats: dict[str, Counter] = defaultdict(Counter)
    disagreements = []
    no_node = in_filter = total_ids = 0
    for r in rows:
        b = bucket(r)
        top = r["tcids"][0]
        fam = ".".join(top.split(".")[:3])
        egg = r["egg"]
        if not egg:
            e = "no_eggnog"
        elif top in egg:
            e = "same"
        elif any(top.startswith(x + ".") for x in egg):
            e = "best_hit_descendant"
        elif any(".".join(x.split(".")[:3]) == fam for x in egg):
            e = "same_family_other_branch"
        else:
            e = "different_family"
            disagreements.append(r)
        stats[b][f"eggnog:{e}"] += 1
        stats[b][f"brite_ko02000:{r['brite_transporter']}"] += 1
        stats[b][f"pfam:{pfam_of.get(r['lt'], 'no_bridge')}"] += 1
        if r["roles"]:
            stats[b][f"cyanorak_transport_role:{any('transport' in x for x in r['roles'])}"] += 1
        subs = set(r["sub_ids"] or [])
        go = set(r["go_subs"] or [])
        stats[b]["go:" + ("no_go" if not go else ("no_best_hit_substrates" if not subs
                          else ("overlap" if subs & go else "disjoint")))] += 1
        stats[b][f"identity:{_ident_bin(r['ident'])}"] += 1
        stats[b][f"split:{len(r['tcids']) > 1}"] += 1
        stats[b][f"fallback_attach:{not top.startswith(r['edge_tc'] + '.') and r['edge_tc'] != top}"] += 1
        for sid in subs:
            total_ids += 1
            no_node += sid not in node_ids
            in_filter += sid in filter_ids
    print(f"Edges with best-hit details: {len(rows)}\n")
    keys = sorted({k for c in stats.values() for k in c})
    print("| measure | training (MED4) | hold-out |")
    print("|---|---|---|")
    for k in keys:
        print(f"| {k} | {stats['training'][k]} | {stats['hold-out'][k]} |")
    print(f"\nDistinct (edge, substrate id) pairs: {total_ids}; ids with no Metabolite/ChemicalRole "
          f"node: {no_node}; ids in slice 1's GO target filter: {in_filter}.")
    print("\n`fallback_attach:True` counts records whose top hit lies outside the edge's TC id "
          "(cross-family top, expected 0).\n")
    print("KO-name substrate overlap: deferred (no structured KO → substrate source).\n")

    print("### Nitrogen sanity table\n")
    print("| strain | gene | best hit | name | identity | best-hit substrates | tier-rollup substrates |")
    print("|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: (r["org"] or "", r["lt"])):
        strain = orgs.get(r["org"], "")
        if strain not in NITROGEN_STRAINS:
            continue
        subs = r["sub_names"] or []
        if not any(k in s.lower() for s in subs for k in N_KEYWORDS):
            continue
        print(f"| {strain} | {r['lt']} {r['gene'] or ''} | {', '.join(r['tcids'])} | "
              f"{(r['names'][0] or '')[:50]} | {r['ident']:.0f} | {'; '.join(sorted(set(subs)))[:80]} | "
              f"{rollup.get(r['edge_tc'], 0)} |")

    print("\n### eggNOG different-family disagreements — random sample of 30 (seed 0) for manual review\n")
    print("| gene | product | best hit | name | identity | eggNOG TC |")
    print("|---|---|---|---|---|---|")
    for r in random.Random(0).sample(disagreements, min(30, len(disagreements))):
        print(f"| {r['lt']} {r['gene'] or ''} | {(r['product'] or '')[:40]} | {r['tcids'][0]} | "
              f"{(r['names'][0] or '')[:45]} | {r['ident']:.0f} | {', '.join(r['egg'])} |")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artifacts-only", action="store_true")
    a = ap.parse_args()
    print("# TCDB best-hit validation (best_hit_v1)\n")
    part_a()
    if not a.artifacts_only:
        part_b()
