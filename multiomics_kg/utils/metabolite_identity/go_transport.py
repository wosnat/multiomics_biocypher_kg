"""Pure helpers assembling GO transport-substrate targets (no network, no I/O)."""
from collections import Counter, defaultdict

from multiomics_kg.utils.chebi_utils import ChebiIndex


def go_substrate_strings(axioms: dict[str, list[str]], kg_terms: set[str],
                         chebi: ChebiIndex):
    """-> ({go_id: [(go_chebi_id, "CHEBI:<primary n>;<name>")]}, skipped_chebi_ids)."""
    out: dict[str, list[tuple[str, str]]] = {}
    skipped: list[str] = []
    for go_id, cids in axioms.items():
        if go_id not in kg_terms:
            continue
        for cid in cids:
            primary = chebi.primary(cid)
            t = chebi.term(primary) if primary else None
            if not primary or not t:
                skipped.append(cid)
                continue
            out.setdefault(go_id, []).append((cid, f"{primary};{t['name']}"))
    return out, skipped


def assemble_go_targets(term_strings, string_to_primary, role_map, role_ids, excluded_primaries):
    substrates: dict = defaultdict(lambda: defaultdict(set))
    roles: dict = defaultdict(lambda: defaultdict(set))
    drops: Counter = Counter()
    for go_id, pairs in term_strings.items():
        for cid, s in pairs:
            if s in role_map:
                rid = role_map[s]
                if rid in role_ids:
                    roles[go_id][rid].add(cid)
                else:
                    drops["role_without_node"] += 1
                continue
            prim = string_to_primary.get(s)
            if prim is None:
                drops["unresolved"] += 1
            elif prim in excluded_primaries:
                drops["filtered"] += 1
            else:
                substrates[go_id][prim].add(cid)
    fin = lambda d: {g: {t: sorted(c) for t, c in m.items()} for g, m in d.items()}
    return fin(substrates), fin(roles), drops


def rewrite_go_targets(go_substrates, canonical_of, excluded_canonical):
    out: dict = {}
    dropped = 0
    for go_id, targets in go_substrates.items():
        merged: dict[str, set] = defaultdict(set)
        for t, cids in targets.items():
            c = canonical_of.get(t, t)
            if c in excluded_canonical:
                dropped += 1
                continue
            merged[c].update(cids)
        if merged:
            out[go_id] = {t: sorted(c) for t, c in merged.items()}
    return out, dropped


def target_distribution(gene_terms_by_strain, go_targets):
    total_genes = 0
    genes: Counter = Counter()
    terms: Counter = Counter()
    strains: dict[str, set] = defaultdict(set)
    for targets in go_targets.values():
        for t in targets:
            terms[t] += 1
    for strain, gene_list in gene_terms_by_strain.items():
        for gterms in gene_list:
            hit = [g for g in gterms if g in go_targets]
            if not hit:
                continue
            total_genes += 1
            for t in {t for g in hit for t in go_targets[g]}:
                genes[t] += 1
                strains[t].add(strain)
    n_terms = len(go_targets)
    rows = [{"target": t, "genes": genes[t],
             "gene_share": round(genes[t] / total_genes, 4) if total_genes else 0.0,
             "terms": terms[t],
             "term_share": round(terms[t] / n_terms, 4) if n_terms else 0.0,
             "strains": len(strains[t])} for t in terms]
    rows.sort(key=lambda r: (-r["genes"], r["target"]))
    return rows
