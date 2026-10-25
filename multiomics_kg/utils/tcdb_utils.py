"""TCDB hierarchy accessor API + curated Pfam→TC map.

Reads cache/data/tcdb/tcdb_hierarchy.json (built by `multiomics_kg.download.build_kegg_metabolism_xrefs` in prepare_data step 6).
"""
from __future__ import annotations

import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

DEFAULT_PATH = Path("cache/data/tcdb/tcdb_hierarchy.json")
_CACHE: dict[str, dict] | None = None

# TCDB's curated cross-references. Both source URLs embed the data in HTML, so
# `download_tcdb_xref_map` strips non-data lines before caching.
PFAM_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/pfam.py"
PFAM_MAP_PATH = Path("cache/data/tcdb/raw/tcdb_pfam_map.tsv")
GO_MAP_URL = "https://www.tcdb.org/cgi-bin/projectv/public/go.py"
GO_MAP_PATH = Path("cache/data/tcdb/raw/tcdb_go_map.tsv")

# Bridge edges are only emitted onto tc_family (level 2) or deeper. A 5-part TCID
# whose nearest surviving ancestor is a tc_class / tc_subclass would produce an
# edge like "this domain relates to Channels and Pores" — true but useless.
MIN_BRIDGE_LEVEL = 2


def load_tcdb() -> dict[str, dict]:
    """Load the TCDB hierarchy JSON. Cached at module level."""
    global _CACHE
    if _CACHE is None:
        with open(DEFAULT_PATH, encoding="utf-8") as f:
            _CACHE = json.load(f)
    return _CACHE


def is_valid_tcdb(tc_id: str) -> bool:
    """Return True if `tc_id` exists as a key in the hierarchy."""
    if not tc_id:
        return False
    return tc_id in load_tcdb()


def tcdb_ancestors(tc_id: str) -> list[str]:
    """Return root-to-parent ancestor chain for `tc_id`. Empty for unknown / root."""
    h = load_tcdb()
    if tc_id not in h:
        return []
    chain: list[str] = []
    parent = h[tc_id].get("parent")
    while parent is not None:
        chain.append(parent)
        parent = h.get(parent, {}).get("parent")
    return list(reversed(chain))


# ============================================================================
# Transport direction (spec 2026-10-22 §5.1) — TCDB-stated only
# ============================================================================

# Whole words incl. plural / agent-noun forms; \b also splits hyphenated tokens
# ("Drug Exporter-1" matches exporter). These DECIDE a direction, so they stay narrow.
_UPTAKE_RE = re.compile(
    r"\b(?:uptake|import|imports|importer|importers|influx)\b", re.IGNORECASE)
# Secretion is export by definition ("General secretory pathway", "Type I secretion system");
# eukaryotic secretory pathways name an organelle and are caught by _ORGANELLE_RE.
_EXPORT_RE = re.compile(
    r"\b(?:export|exports|exporter|exporters|efflux|extrusion|secretion|secretions|secretory)\b",
    re.IGNORECASE)
# Broader verb forms, used only to VETO a direction from a node's description.
_UPTAKE_ANY_RE = re.compile(
    r"\b(?:uptake|import|imports|imported|importing|importer|importers|influx)\b",
    re.IGNORECASE)
_EXPORT_ANY_RE = re.compile(
    r"\b(?:export|exports|exported|exporting|exporter|exporters|efflux|effluxes"
    r"|extrusion|extrudes?|extruded|extruding|secretion|secretions|secretory|secretes?"
    r"|secreted|secreting)\b", re.IGNORECASE)
# The curator was unsure. Anywhere before the direction word in the name counts
# ("Putative heavy metal (Me 2+) exporter", "Probable pyruvate:H+ importer"); a
# hedge AFTER it does not ("uptake permease, AmtB (may need ...)" stays uptake).
_HEDGE_RE = re.compile(
    r"\b(?:may|might|could|likely|possibl[ey]|potential(?:ly)?|putative|probabl[ey]"
    r"|presum\w*|predicted|suggested|proposed|unknown|uncharacteri[sz]ed|hypothetical"
    r"|fanciful|incorrectly|regulators?)\b|\?", re.IGNORECASE)
# The node is named by resemblance only. Weak right after the direction word too
# ("Lipid A exporter homologue"), and a hedged/generic name never inherits.
_MEMBER_RE = re.compile(
    r"\b(?:homolog\w*|paralog\w*|\w+-like|like protein|DUF\d+|UPF\d+"
    r"|family (?:member|protein))\b|\?", re.IGNORECASE)
# Organelle anywhere in the name: import into a peroxisome, export from the nucleus,
# loading of vesicles, "uptake into vacuoles" are not uptake/export of a bacterial cell.
_ORGANELLE_RE = re.compile(
    r"\b(?:perox[iy]som\w*|mitochondri\w*|chloroplast\w*|plastid\w*|apicoplast\w*"
    r"|thylakoid\w*|nuclear|nucleus|lysosom\w*|vacuol\w*|tonoplast\w*|golgi"
    r"|endoplasmic|endosom\w*|vesic\w*|synap\w*|hydrogenosom\w*|glycosom\w*"
    r"|organell\w*|lamellar bod\w*)\b", re.IGNORECASE)
_ER_RE = re.compile(r"\bER\b")  # case-sensitive: the endoplasmic reticulum
# The protein takes part in a process named by the direction word, rather than moving the
# substrate itself: "K+ channel, Kir4.2, involved in bile secretion", "controls salivary fluid
# secretion", "affects secretion of proteins". Anywhere before the direction word.
_ROLE_BEFORE_RE = re.compile(
    r"\b(?:involved in|controls?|affects?|required for|contributes? to|participates? in"
    r"|regulates?)\b", re.IGNORECASE)
# "not involved in the export", "non-exporter", "unable to import".
_NEGATION_BEFORE_RE = re.compile(
    r"\b(?:not|no|non|never|unable|lacks?|without)\b(?:\W+\w+){0,3}\W*$", re.IGNORECASE)
# Clause end after the direction word (for the homologue check only).
_CLAUSE_END_RE = re.compile(r"[,;(\[]")
# A description that says the direction is two-way or contested vetoes it, wherever
# it says so. ("Unknown function" is not here: descriptions use it for accessory
# proteins; only the node's own name counts for that, see _UNKNOWN_FUNCTION_RE.)
_DESC_VETO_RE = re.compile(
    r"\bboth\s+(?:uptake|import|influx)\w*\s+and\s+(?:export|efflux|extrusion)\w*"
    r"|\bboth\s+(?:export|efflux|extrusion)\w*\s+and\s+(?:uptake|import|influx)\w*"
    r"|\bbidirectional\w*|\bboth directions\b|\bunder debate\b|\bcontroversial\b"
    r"|\bdisputed\b|\bhypothes\w*", re.IGNORECASE)
_UNKNOWN_FUNCTION_RE = re.compile(
    r"\b(?:unknown function|function unknown|uncharacteri[sz]ed)\b", re.IGNORECASE)
# The description's first sentence is TCDB's full system name (the stored name is
# truncated: "The antimicrobial peptide" for "... HBD-2 exporter, YejABEF"). Only it
# is read for the opposite direction word; later sentences discuss other proteins
# and assays.
_FIRST_SENTENCE_RE = re.compile(r"(?<=[A-Za-z0-9\)\]])\.\s+(?=[A-Z(])")
# Description phrases that use a direction word about something else.
_DESC_IGNORE_RE = re.compile(r"\bsignal (?:sequence|peptide)s? for export\b", re.IGNORECASE)
_ABC = "3.A.1"
# ABC subfamilies 3.A.1.1-99 are importers (binding-protein-dependent + ECF), so
# the number alone is enough there. 3.A.1.100+ is NOT a reliable export claim:
# it holds ABC-F ribosome-protection proteins (.120/.121), FtsEX (.140), gliding
# motility (.132), uncharacterized families, and the eukaryotic ABCA-G groups
# (.2xx, grouped by fold), some of them importers (BacA .203.4). Those need a name.
_ABC_UPTAKE_MAX = 99

DIRECTION_OVERRIDES_PATH = Path("config/tcdb_direction_overrides.yaml")
_OVERRIDE_SCOPES = {"subtree"}
_OVERRIDES_CACHE: dict[str, dict] | None = None


def load_direction_overrides(path: Path | str | None = None) -> dict[str, dict]:
    """``{tc_id: {scope, reason}}`` from ``config/tcdb_direction_overrides.yaml``.

    ``scope: subtree`` — the node and every descendant are ``unknown``. Overrides
    only ever REMOVE a direction; there is no way to assert one.
    """
    global _OVERRIDES_CACHE
    if path is None and _OVERRIDES_CACHE is not None:
        return _OVERRIDES_CACHE
    import yaml

    p = Path(path) if path is not None else DIRECTION_OVERRIDES_PATH
    raw = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    out: dict[str, dict] = {}
    for tc_id, entry in raw.items():
        tc_id = str(tc_id)
        if not isinstance(entry, dict) or entry.get("scope") not in _OVERRIDE_SCOPES \
                or not entry.get("reason"):
            raise ValueError(f"{p}: {tc_id}: need scope in {sorted(_OVERRIDE_SCOPES)} and a reason")
        out[tc_id] = {"scope": entry["scope"], "reason": str(entry["reason"])}
    if path is None:
        _OVERRIDES_CACHE = out
    return out


def _lineage(tc_id: str, hierarchy: dict) -> list[str]:
    """``tc_id`` and its ancestors, nearest first."""
    out, node = [], tc_id
    while node is not None and node in hierarchy and node not in out:
        out.append(node)
        node = hierarchy[node].get("parent")
    return out


def _weak_direction_word(name: str, match: re.Match) -> bool:
    before, after = name[:match.start()], name[match.end():]
    after_clause = _CLAUSE_END_RE.split(after, maxsplit=1)[0]
    return bool(_HEDGE_RE.search(before) or _NEGATION_BEFORE_RE.search(before)
                or _ROLE_BEFORE_RE.search(before)
                or _ORGANELLE_RE.search(name) or _ER_RE.search(name)
                or _MEMBER_RE.search(after_clause))


def _name_direction(name: str) -> tuple[str, str | None]:
    """Direction stated by one TCDB name: ``(direction, word_seen)``.

    ``direction`` is ``uptake`` / ``export`` / ``unknown``. ``word_seen`` is the
    direction the name stated even when it was rejected as weak (``"both"`` for a
    name stating both), so the ABC branch can tell "the name says nothing" from
    "the name says export, tentatively".
    """
    up, out = _UPTAKE_RE.search(name), _EXPORT_RE.search(name)
    if up and out:
        return "unknown", "both"
    if not (up or out):
        return "unknown", None
    word = "uptake" if up else "export"
    if _weak_direction_word(name, up or out):
        return "unknown", word
    return word, word


def _description_vetoes(tc_id: str, hierarchy: dict, direction: str) -> bool:
    """The node's own TCDB description contradicts or doubts ``direction``."""
    desc = _DESC_IGNORE_RE.sub("", hierarchy.get(tc_id, {}).get("description") or "")
    if not desc:
        return False
    if _DESC_VETO_RE.search(desc):
        return True
    opposite = _EXPORT_ANY_RE if direction == "uptake" else _UPTAKE_ANY_RE
    return bool(opposite.search(_FIRST_SENTENCE_RE.split(desc, maxsplit=1)[0]))


def tcdb_direction(tc_id: str, hierarchy: dict,
                   overrides: dict | None = None) -> tuple[str, str | None]:
    """Transport direction TCDB states for ``tc_id``: ``(direction, basis)``.

    Conservative: ``unknown`` whenever the evidence is weak, conflicting or
    contested; a false ``unknown`` is cheap, a wrong direction is not. A direction
    is never inherited from a family or subfamily NAME: a family is named after
    its characterised member, and its relatives (reached by homology) may go the
    other way or have no known function.

    1. ABC numbering (basis ``abc_numbering``), every level from the subfamily
       down: ``3.A.1.<n>`` with ``n <= 99`` is ``uptake`` (binding-protein-
       dependent and ECF importers) unless a name from the node up to its
       subfamily says export or both (even tentatively), the node's own name says
       its function is unknown, or its description doubts it. ``n >= 100`` has no
       numbering direction (ABC-F, FtsEX, uncharacterized and eukaryotic ABCA-G
       families live there).
    2. Own name (basis ``tcdb_name``), ``tc_system`` nodes only: the system's own
       name with exactly one direction word ({uptake, import(s), importer(s),
       influx} / {export(s), exporter(s), efflux, extrusion, secretion, secretory},
       whole word; secretion is export) decides, unless it is weak: a hedge word
       anywhere before the direction word (putative, probable, may, uncharacterized,
       ``?``, ...) or a role phrase ("involved in", "controls", "affects",
       "required for"), a negation just before it, an organelle anywhere in the name, or a resemblance word right
       after it ("exporter homologue"). Family, subfamily, class and subclass
       names never give a direction, to themselves or to their systems.
    3. The node's own TCDB description can only veto: both directions,
       bidirectional, under debate, disputed, hypothesis (anywhere), or the
       opposite direction word in its first sentence (TCDB's full system name).
    4. ``config/tcdb_direction_overrides.yaml`` ``subtree`` entries make the node
       and its descendants ``unknown``.
    """
    if overrides is None:
        overrides = load_direction_overrides()
    lineage = _lineage(tc_id, hierarchy)
    if any(n in overrides for n in lineage) \
            or (tc_id not in hierarchy and any(tc_id.startswith(k + ".") for k in overrides)):
        return "unknown", None
    own = hierarchy.get(tc_id, {})
    parts = tc_id.split(".")
    if len(parts) >= 4 and ".".join(parts[:3]) == _ABC and parts[3].isdigit() \
            and int(parts[3]) <= _ABC_UPTAKE_MAX:
        subfamily = ".".join(parts[:4])
        names = [hierarchy[n].get("name") or "" for n in lineage[:lineage.index(subfamily) + 1]] \
            if subfamily in lineage else [own.get("name") or ""]
        if any(_name_direction(n)[1] in ("export", "both") for n in names) \
                or _UNKNOWN_FUNCTION_RE.search(own.get("name") or "") \
                or _description_vetoes(tc_id, hierarchy, "uptake"):
            return "unknown", None
        return "uptake", "abc_numbering"
    if own.get("level_kind", "tc_system" if len(parts) == 5 else None) != "tc_system":
        return "unknown", None
    direction, _ = _name_direction(own.get("name") or "")
    if direction == "unknown" or _description_vetoes(tc_id, hierarchy, direction):
        return "unknown", None
    return direction, "tcdb_name"


# ============================================================================
# Curated Pfam → TC map (backs the Pfam_in_tcdb_family bridge edge)
# ============================================================================
#
# TCDB publishes a curated table of which Pfam domains are associated with which
# transporter families. It is the one piece of diamond↔eggNOG cross-checking that
# is NOT derivable from the graph, so it is modelled as a bridge edge
# (Pfam -> TcdbFamily), mirroring Pfam_in_interpro_entry.
#
# Previously downloaded + consumed inside the /tcdb-diamond runner, which made the
# runner depend on downstream annotation state. It now belongs to prepare_data
# step 6 alongside the other TCDB reference TSVs.


def _valid_pfam(token: str) -> bool:
    return token.startswith("PF") and len(token) >= 7 and token[2:7].isdigit()


def _valid_go(token: str) -> bool:
    return token.startswith("GO:") and token[3:].isdigit()


_XREF_KINDS = {
    "pfam": (PFAM_MAP_URL, PFAM_MAP_PATH, _valid_pfam),
    "go": (GO_MAP_URL, GO_MAP_PATH, _valid_go),
}


def download_tcdb_xref_map(kind: str, dest: Path | None = None, force: bool = False) -> Path:
    """Download a TCDB cross-reference map to a clean 3-column TSV.

    ``kind`` is ``"pfam"`` or ``"go"``. Writes ``<xref_id> \\t <TC_id> \\t <name>``.
    Skipped when the file exists and ``force`` is False.
    """
    if kind not in _XREF_KINDS:
        raise ValueError(f"unknown TCDB xref kind {kind!r}; expected one of {sorted(_XREF_KINDS)}")
    url, default_path, is_valid = _XREF_KINDS[kind]
    out = Path(dest) if dest is not None else default_path
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not force:
        return out

    with urllib.request.urlopen(url) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    with open(out, "w", encoding="utf-8") as fh:
        for line in text.splitlines():
            parts = line.rstrip().split("\t")
            if len(parts) < 3 or not is_valid(parts[0]) or parts[1].count(".") < 2:
                continue
            fh.write("\t".join(parts[:3]) + "\n")
    return out


def load_tcdb_xref_map(kind: str, path: Path | None = None) -> dict[str, set[str]]:
    """Load a TCDB cross-reference map into ``{xref_id: {full_tcid, ...}}``.

    TCIDs are kept at their published (5-part) depth — the roll-up to a surviving
    node happens in `build_tcdb_bridges`, which needs the original depth to record
    provenance. Returns an empty dict when the file is absent.
    """
    if kind not in _XREF_KINDS:
        raise ValueError(f"unknown TCDB xref kind {kind!r}; expected one of {sorted(_XREF_KINDS)}")
    _, default_path, is_valid = _XREF_KINDS[kind]
    p = Path(path) if path is not None else default_path
    if not p.exists():
        return {}
    result: dict[str, set[str]] = defaultdict(set)
    with open(p, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not is_valid(parts[0]) or parts[1].count(".") < 2:
                continue
            result[parts[0]].add(parts[1])
    return dict(result)


def build_tcdb_bridges(
    xref_to_tcids: dict[str, set[str]],
    hierarchy: dict[str, dict],
    kept: set[str],
    min_level: int = MIN_BRIDGE_LEVEL,
) -> dict[str, dict[str, list[str]]]:
    """Roll a published xref→TCID map onto the pruned node set.

    TCDB publishes these maps at 5-part `tc_system` depth, but the pruned
    hierarchy keeps only gene-anchored TCIDs + ancestors plus the substrate-bearing
    systems of kept families (range-only) — so some published TCIDs have no node.
    Each is walked up to its nearest **kept** node (itself when kept); ``kept``
    is the full kept set, gene-anchored and range-only alike.
    Targets shallower than ``min_level`` are dropped as uninformative.

    Returns ``{tc_id: {xref_id: [original_tcids, ...]}}`` — the original
    pre-roll-up TCIDs are retained so the edge can record exactly what TCDB
    curated, and roll-up precision is not silently lost.

    **Semantics of the result**: an entry means *"TCDB's curated reference
    proteins for this transport family carry this domain / this GO annotation"*.
    That is a statement about the composition of the TC family. It reads
    correctly outward from the TcdbFamily; it does **not** license the reverse
    inference that anything carrying the xref belongs to the family.
    """
    bridges: dict[str, dict[str, list[str]]] = {}
    for xref_id, tcids in xref_to_tcids.items():
        for tcid in tcids:
            anchor: str | None = tcid if tcid in kept else None
            if anchor is None:
                cur = hierarchy.get(tcid, {}).get("parent")
                while cur is not None:
                    if cur in kept:
                        anchor = cur
                        break
                    cur = hierarchy.get(cur, {}).get("parent")
            if anchor is None:
                continue
            if hierarchy.get(anchor, {}).get("level", 0) < min_level:
                continue
            bridges.setdefault(anchor, {}).setdefault(xref_id, []).append(tcid)
    return {
        tc: {x: sorted(set(v)) for x, v in sorted(xrefs.items())}
        for tc, xrefs in sorted(bridges.items())
    }


# ============================================================================
# Retired TCDB ids → the ids TCDB moved them to (post-3a P2)
# ============================================================================
#
# TCDB retires a family by renumbering its systems; the retired family page then
# only lists lines like "9.A.40.1.2 has been moved to 1.C.126.1.2". eggNOG and the
# diamond library still carry the old ids, so step 6 re-anchors those seeds onto
# the current ids. TCDB sometimes reused a number for two different moves; such ids
# are ambiguous and are only aliased through config/tcdb_moved_overrides.yaml.

_MOVED_RE = re.compile(
    r"\b(\d+\.[A-Za-z]\.\d+(?:\.\d+)*)\s+has been moved to\s+(\d+\.[A-Za-z]\.\d+(?:\.\d+)*)")


def _norm_tcid(tc: str) -> str:
    """TC ids carry one upper-case letter (``1.a.17.4.1`` → ``1.A.17.4.1``)."""
    return tc.strip().upper()


def parse_moved_tcids(page_texts) -> dict[str, set[str]]:
    """``{old id: {new ids}}`` from TCDB page HTML/text "<old> has been moved to <new>"
    lines. Tags are stripped first; ids are upper-cased. An old id with more than one
    target (TCDB reused the number) keeps every target, so callers can see the conflict."""
    out: dict[str, set[str]] = defaultdict(set)
    for text in page_texts:
        flat = re.sub(r"<[^>]*>", " ", text)
        for old, new in _MOVED_RE.findall(flat):
            out[_norm_tcid(old)].add(_norm_tcid(new))
    return dict(out)


def _prefix(tc: str, n_parts: int) -> str:
    return ".".join(tc.split(".")[:n_parts])


def current_tcids(families: set[str], systems: set[str]) -> set[str]:
    """Ids TCDB still lists: 3-part families (``families.tsv``), 5-part systems
    (``acc2tcid.tsv`` / ``superfamilies.tsv``) and the 4-part subfamilies of those systems."""
    out = set(families) | set(systems)
    out |= {_prefix(s, 4) for s in systems if s.count(".") >= 4}
    return out


def still_current_moves(raw_moves: dict, current: set[str] | None) -> list[str]:
    """Parsed old ids that TCDB still lists as current (never aliased)."""
    if not current:
        return []
    return sorted(old for old in raw_moves if old in current)


def resolve_moved_aliases(
    raw_moves: dict[str, set[str]],
    hierarchy: dict,
    overrides: dict[str, str] | None = None,
    current: set[str] | None = None,
) -> tuple[dict[str, str], list[str]]:
    """Turn parsed moves into seed aliases ``{old id: new id}``. Returns
    ``(aliases, ambiguous)`` where ``ambiguous`` lists the ids with conflicting targets.

    * An id with exactly one target that is in ``hierarchy`` is aliased to it.
    * A coarser id (subfamily / family prefix of moved systems) is aliased to target
      prefix T at the same depth only when EVERY system the hierarchy holds under it
      was moved, none ambiguously, and all their targets share T. Otherwise no alias.
    * ``overrides`` (curated) win over everything parsed. An override on a coarse id
      also decides its ambiguous descendants by prefix substitution (override
      ``9.B.99.1 -> 1.B.93.1`` maps ``9.B.99.1.3 -> 1.B.93.1.3``), falling back to the
      override target itself when the substituted id is not in ``hierarchy``.
    * ``current`` (ids TCDB still lists, from `current_tcids`): a parsed old id in it is
      treated as not moved — never aliased, and it blocks coarse aliases above it.
    A target missing from ``hierarchy`` gives no alias (the caller falls back to
    ancestor truncation).
    """
    aliases: dict[str, str] = {}
    current = current or set()
    raw_moves = {old: t for old, t in raw_moves.items() if old not in current}
    ambiguous = sorted(old for old, targets in raw_moves.items() if len(targets) > 1)
    single = {old: next(iter(t)) for old, t in raw_moves.items() if len(t) == 1}
    for old, new in single.items():
        if new in hierarchy:
            aliases[old] = new
    coarse: set[str] = set()
    for old in raw_moves:
        n = old.count(".") + 1
        for depth in range(3, n):
            coarse.add(_prefix(old, depth))
    for c in sorted(coarse - current):
        depth = c.count(".") + 1
        systems = sorted({t for t in hierarchy if t.startswith(c + ".")
                          and t.count(".") >= depth} | {o for o in raw_moves
                                                         if o.startswith(c + ".")})
        leaves = [s for s in systems if not any(o.startswith(s + ".") for o in systems)]
        if not leaves or any(s not in single for s in leaves):
            continue
        targets = {_prefix(single[s], depth) for s in leaves}
        if len(targets) != 1:
            continue
        target = targets.pop()
        if target.count(".") + 1 == depth and target in hierarchy:
            aliases[c] = target
    for old, new in (overrides or {}).items():
        aliases[old] = new
        for amb in ambiguous:
            if amb.startswith(old + ".") and amb not in (overrides or {}):
                sub = new + amb[len(old):]
                aliases[amb] = sub if sub in hierarchy else new
    return dict(sorted(aliases.items())), ambiguous


def load_moved_overrides(path: Path, hierarchy: dict) -> dict[str, str]:
    """``{old id: target}`` from ``config/tcdb_moved_overrides.yaml``
    (``{old_id: {target, reason}}``). Missing file → ``{}``. Raises ValueError on an
    entry without exactly a non-empty ``target`` and ``reason``, or a target that is not
    in ``hierarchy``."""
    import yaml
    path = Path(path)
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a mapping")
    out: dict[str, str] = {}
    for old, e in data.items():
        if not isinstance(e, dict) or set(e) != {"target", "reason"} \
                or not e["target"] or not e["reason"]:
            raise ValueError(f"{path}: entry {old!r} needs exactly a non-empty target "
                             f"and reason, got {e!r}")
        target = _norm_tcid(str(e["target"]))
        if target not in hierarchy:
            raise ValueError(f"{path}: {old!r} target {target!r} is not in the TCDB hierarchy")
        out[_norm_tcid(str(old))] = target
    return out
