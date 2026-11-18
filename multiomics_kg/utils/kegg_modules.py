"""KEGG module definitions: parse a module's formula into steps and score one organism against it.

Grammar (KEGG Mapper convention):
  * a space at the top level separates **steps**; a step is covered when its expression is true;
  * ``,`` separates alternatives (OR); parentheses group, and inside them a space means AND;
  * ``+`` joins the subunits of one complex (all required); ``-K…`` marks an optional subunit;
  * ``--`` is a reaction KEGG has no KO for; as a whole step it is dropped from the step list
    (otherwise the module could never be complete), as an alternative it is never true;
  * a definition on several lines is the module's separate parts — the lines are joined as
    further steps;
  * ``M0xxxx`` names another module, true when that module is complete for the organism.

EC fallback (``evaluate_organism``'s ``ec_candidates``): an EC-only gene may stand in for a KO
that is not inside a complex, each such gene fills at most one step of a module (maximum
bipartite matching), and only over the steps the KO evidence leaves uncovered.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_TOKEN = re.compile(r"K\d{5}|M\d{5}|--|[(),+\- ]")


@dataclass(frozen=True)
class Step:
    index: int          # 1-based, over the steps left after dropping whole-step ``--``
    definition: str     # the step's formula text, verbatim


@dataclass(frozen=True)
class ModuleResult:
    steps_total: int
    steps_present: int
    steps_present_via_ec: int
    missing_step_indices: list[int] = field(default_factory=list)
    missing_step_definitions: list[str] = field(default_factory=list)
    missing_kos: list[str] = field(default_factory=list)
    ec_filled_step_indices: list[int] = field(default_factory=list)
    ec_genes: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        missing = self.steps_total - self.steps_present
        return {0: "complete", 1: "one_step_missing", 2: "two_steps_missing"}.get(missing, "partial")


# --- parsing -------------------------------------------------------------------------------

def _tokenize(text: str) -> list[str]:
    tokens = _TOKEN.findall(text)
    if "".join(tokens).replace(" ", "") != text.replace(" ", ""):
        raise ValueError(f"unparseable module definition: {text!r}")
    return tokens


def _split_top(tokens: list[str], sep: str) -> list[list[str]]:
    parts, cur, depth = [], [], 0
    for t in tokens:
        if t == "(":
            depth += 1
        elif t == ")":
            depth -= 1
        if t == sep and depth == 0:
            parts.append(cur)
            cur = []
        else:
            cur.append(t)
    parts.append(cur)
    return [p for p in parts if any(t != " " for t in p)]


def _strip(tokens: list[str]) -> list[str]:
    while tokens and tokens[0] == " ":
        tokens = tokens[1:]
    while tokens and tokens[-1] == " ":
        tokens = tokens[:-1]
    if tokens and tokens[0] == "(" and tokens[-1] == ")":
        depth = 0
        for i, t in enumerate(tokens):
            depth += (t == "(") - (t == ")")
            if depth == 0 and i < len(tokens) - 1:
                return tokens          # "(A) (B)": outer parens are not one group
        return _strip(tokens[1:-1])
    return tokens


# Tree nodes: ("or", [..]) ("and", [..]) ("cx", [(optional, node), ..]) ("ko", K) ("mod", M) ("none",)
def _parse(tokens: list[str]):
    tokens = _strip(tokens)
    if not tokens:
        raise ValueError("empty module expression")
    alts = _split_top(tokens, ",")
    if len(alts) > 1:
        return ("or", [_parse(a) for a in alts])
    terms = _split_top(tokens, " ")
    if len(terms) > 1:
        return ("and", [_parse(t) for t in terms])
    groups, buf, depth, optional = [], [], 0, False
    for t in (t for t in tokens if t != " "):
        if t == "(":
            depth += 1
        elif t == ")":
            depth -= 1
        if depth == 0 and t in "+-" and buf:
            groups.append((optional, buf))
            buf, optional = [], t == "-"
            continue
        if depth == 0 and t == "-" and not buf:
            optional = True
            continue
        buf.append(t)
    if buf:
        groups.append((optional, buf))
    if len(groups) > 1:
        return ("cx", [(o, _parse(g)) for o, g in groups])
    only = groups[0][1]
    if len(only) == 1:
        t = only[0]
        if t == "--":
            return ("none",)
        return ("ko", t) if t[0] == "K" else ("mod", t)
    return _parse(only)


def _step_trees(definition: str) -> list[tuple[Step, tuple]]:
    tokens = _tokenize(definition.replace("\n", " "))
    out = []
    for part in _split_top(tokens, " "):
        tree = _parse(part)
        if tree[0] == "none":
            continue
        out.append((Step(len(out) + 1, "".join(part).strip()), tree))
    return out


def parse_steps(definition: str) -> list[Step]:
    """The module's steps, multi-line definitions joined, whole-step ``--`` dropped."""
    return [s for s, _ in _step_trees(definition)]


def _tree(step: Step):
    return _parse(_tokenize(step.definition))


def _kos(node, in_complex=False, optional=False, out=None):
    """(KO, in_complex, optional) for every KO in the tree."""
    out = [] if out is None else out
    kind = node[0]
    if kind == "ko":
        out.append((node[1], in_complex, optional))
    elif kind in ("or", "and"):
        for c in node[1]:
            _kos(c, in_complex, optional, out)
    elif kind == "cx":
        for opt, c in node[1]:
            _kos(c, True, optional or opt, out)
    return out


def _mods(node, out=None):
    out = [] if out is None else out
    if node[0] == "mod":
        out.append(node[1])
    elif node[0] in ("or", "and"):
        for c in node[1]:
            _mods(c, out)
    elif node[0] == "cx":
        for _, c in node[1]:
            _mods(c, out)
    return out


def ko_roles(steps: list[Step]) -> dict[str, tuple[list[int], str]]:
    """KO → (sorted step indices it serves, ``required`` | ``optional``)."""
    idx: dict[str, set[int]] = {}
    role: dict[str, set[bool]] = {}
    for s in steps:
        for ko, _, opt in _kos(_tree(s)):
            idx.setdefault(ko, set()).add(s.index)
            role.setdefault(ko, set()).add(opt)
    conflict = sorted(k for k, r in role.items() if len(r) > 1)
    if conflict:
        raise ValueError(f"KO required in one step and optional in another: {conflict}")
    return {k: (sorted(idx[k]), "optional" if True in role[k] else "required") for k in idx}


# --- evaluation ----------------------------------------------------------------------------

def _true(node, kos, ec_kos, mod_ok, in_complex=False) -> bool:
    kind = node[0]
    if kind == "none":
        return False
    if kind == "mod":
        return mod_ok(node[1])
    if kind == "ko":
        return node[1] in kos or (not in_complex and node[1] in ec_kos)
    if kind == "or":
        return any(_true(c, kos, ec_kos, mod_ok, in_complex) for c in node[1])
    if kind == "and":
        return all(_true(c, kos, ec_kos, mod_ok, in_complex) for c in node[1])
    return all(_true(c, kos, ec_kos, mod_ok, True) for opt, c in node[1] if not opt)


def _match(candidates: dict[int, set[str]]) -> dict[int, str]:
    """Maximum bipartite matching step → gene (augmenting paths, deterministic order)."""
    owner: dict[str, int] = {}

    def augment(step: int, seen: set[str]) -> bool:
        for gene in sorted(candidates[step]):
            if gene in seen:
                continue
            seen.add(gene)
            if gene not in owner or augment(owner[gene], seen):
                owner[gene] = step
                return True
        return False

    for step in sorted(candidates):
        augment(step, set())
    return {step: gene for gene, step in owner.items()}


def evaluate_organism(modules: dict[str, list[Step]], kos: set[str],
                      ec_candidates: dict[str, set[str]]) -> dict[str, ModuleResult]:
    """Score every module in ``modules`` for one organism.

    ``kos`` are the organism's KOs; ``ec_candidates`` maps a KO to the locus tags of EC-only
    genes that may stand in for it. Referenced modules are evaluated first; a reference cycle
    or a reference to a module not in ``modules`` raises ``ValueError``.
    """
    trees = {m: [(s, _tree(s)) for s in steps] for m, steps in modules.items()}
    ec_kos = {k for k, genes in ec_candidates.items() if genes}
    results: dict[str, ModuleResult] = {}
    active: set[str] = set()

    def result(m: str) -> ModuleResult:
        if m in results:
            return results[m]
        if m not in trees:
            raise ValueError(f"module reference to unknown module {m}")
        if m in active:
            raise ValueError(f"module reference cycle through {m}")
        active.add(m)
        results[m] = _evaluate(trees[m])
        active.discard(m)
        return results[m]

    def ko_complete(m: str) -> bool:
        r = result(m)
        return r.status == "complete" and r.steps_present_via_ec == 0

    def any_complete(m: str) -> bool:
        return result(m).status == "complete"

    def _evaluate(step_trees) -> ModuleResult:
        ko_ok = [_true(t, kos, set(), ko_complete) for _, t in step_trees]
        # Steps a reference alone makes true (no gene consumed): EC pass without EC genes.
        ref_ok = [not a and _true(t, kos, set(), any_complete) for a, (_, t) in zip(ko_ok, step_trees)]
        need = {}
        for i, (s, t) in enumerate(step_trees):
            if ko_ok[i] or ref_ok[i] or not _true(t, kos, ec_kos, any_complete):
                continue
            genes = set()
            for ko, in_cx, _ in _kos(t):
                if not in_cx:
                    genes |= ec_candidates.get(ko, set())
            need[s.index] = genes
        matched = _match(need)
        ec_fill: dict[int, str] = dict(matched)
        for i, (s, t) in enumerate(step_trees):
            if ref_ok[i]:
                refs = [r for r in _mods(t) if any_complete(r) and not ko_complete(r)]
                ec_fill[s.index] = sorted(refs)[0]
        present = [ko_ok[i] or s.index in ec_fill for i, (s, _) in enumerate(step_trees)]
        missing = [(s, t) for (s, t), p in zip(step_trees, present) if not p]
        missing_kos = sorted({ko for _, t in missing for ko, _, opt in _kos(t) if not opt})
        filled = sorted(ec_fill)
        return ModuleResult(
            steps_total=len(step_trees),
            steps_present=sum(present),
            steps_present_via_ec=len(filled),
            missing_step_indices=[s.index for s, _ in missing],
            missing_step_definitions=[s.definition for s, _ in missing],
            missing_kos=missing_kos,
            ec_filled_step_indices=filled,
            ec_genes=[ec_fill[i] for i in filled],
        )

    for m in modules:
        result(m)
    return {m: results[m] for m in modules}
