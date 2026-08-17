"""Pure-Python parsing of the InterPro reference release files.

No filesystem, no network — the download + caching orchestrator lives in
``multiomics_kg/download/build_interpro_reference.py`` (prepare_data step 9).
This module is the unit-testable core: it turns the two InterPro
``current_release`` text files into the committed reference dict consumed by
``interpro_adapter`` at KG-build time.

Two source files (InterPro FTP ``current_release/``):

- ``entry.list`` — tab-separated, one header line, then
  ``ENTRY_AC <tab> ENTRY_TYPE <tab> ENTRY_NAME``. Covers *every* InterPro entry
  (names + types). ``ENTRY_TYPE`` is Title_case (``Family``,
  ``Homologous_superfamily``, ``Conserved_site`` …); we uppercase it to match the
  ``interpro_type`` strings InterProScan writes into calls.json (``FAMILY``,
  ``HOMOLOGOUS_SUPERFAMILY``, ``CONSERVED_SITE`` …).

- ``ParentChildTreeFile.txt`` — the *is-a* hierarchy only (not the cross-type
  "contains/found-in" relationship). Each line is ``(--)*IPRxxxxxx::Name::…`` where
  the count of leading ``--`` pairs is the depth. A line's parent is the nearest
  preceding line at depth-1. Most entries never appear here (they are parentless
  roots → level 0, parent None).

- ``interpro2go`` — entry→GO mappings, one per line:
  ``InterPro:IPRxxxxxx <name> > GO:<go name> ; GO:xxxxxxx``. Comment lines start
  with ``!``. ~30K mappings over ~14.8K entries.

- ``interpro.xml.gz`` — the full entry XML; the only source of entry→pathway
  (``<db_xref db="METACYC"|"REACTOME" …/>``), entry→**EC** (``db="EC"``) and
  entry→**CAZy** (``db="CAZY"``) cross-references. InterProScan emits GO and
  pathways but never EC or CAZy, so the reference is the *only* source of those.
  Streamed line-by-line by the caller, so this module takes an *iterable of lines*
  (:func:`parse_entry_db_xrefs`, one pass for all three).

  **There are no KEGG pathway xrefs in InterPro** — ``KEGG`` is a legacy token in
  ``interpro.dtd``'s allowed-``db`` list but carries zero entries (verified
  against the full release, 2026-08-06). Reactome xrefs are species-expanded
  (~507K, one per Reactome pathway *per organism*) and are excluded by default as
  noise for bacterial genomes; MetaCyc (~80K) is the useful bacterial set.

Why these files rather than re-running InterProScan with ``--goterms
--pathways``: the xrefs InterProScan writes into ``calls.json`` come from
``signature.entry`` and are therefore *entry-level* — every match of the same
``IPR`` accession carries an identical GO/pathway set. Looking them up from the
reference release is equivalent to, and ~27 wallclock-hours cheaper than, a
42-strain re-scan.

Combined result: ``{IPRxxxxxx: {"name", "type", "parent", "level"[, "go_terms"]
[, "pathways"][, "ec_numbers"][, "cazy_ids"][, "description"]}}`` — ``name``/
``type`` from entry.list (authoritative, untruncated), ``parent``/``level``
from the tree (``None``/``0`` when the entry is not in the tree).
``go_terms``, ``pathways``, ``ec_numbers``, ``cazy_ids`` and ``description``
are all **sparse**: the key is absent rather than an empty list/string when
the entry has none, which keeps the committed JSON small (only ~27% of
entries carry GO, ~39% of *observed* entries carry EC, far fewer CAZy).
``ec_numbers``/``cazy_ids`` are raw InterPro tokens; normalization is the
consumer's job. ``description`` is the first paragraph of the entry's
``<abstract>`` (plain text, tags stripped, capped to 400 chars) — see
:func:`clean_abstract` / :func:`parse_entry_descriptions`, also parsed from
``interpro.xml`` in the same streaming pass as the db_xrefs.

See ``docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md``.
"""

from __future__ import annotations

import html
import re
from collections.abc import Iterable

# ENTRY_TYPE strings as they appear in entry.list → the UPPER form InterProScan
# emits in calls.json. We normalize via ``.upper()`` which already maps every
# known value correctly (``Homologous_superfamily`` → ``HOMOLOGOUS_SUPERFAMILY``);
# this set is the allow-list used to flag unexpected new types.
KNOWN_INTERPRO_TYPES = {
    "FAMILY",
    "DOMAIN",
    "HOMOLOGOUS_SUPERFAMILY",
    "REPEAT",
    "CONSERVED_SITE",
    "ACTIVE_SITE",
    "BINDING_SITE",
    "PTM",
}

_TREE_LINE_RE = re.compile(r"^(?P<dashes>(?:--)*)(?P<acc>IPR\d{6})::")

# `InterPro:IPR000003 Retinoid X receptor > GO:DNA binding ; GO:0003677`
_IPR2GO_RE = re.compile(r"^InterPro:(?P<acc>IPR\d{6})\s.*?;\s*(?P<go>GO:\d{7})\s*$")

_XML_ENTRY_OPEN_RE = re.compile(r'<interpro\s+id="(IPR\d{6})"')
_XML_ENTRY_CLOSE = "</interpro>"
_XML_DB_XREF_RE = re.compile(r'<db_xref[^>]*\bdb="([A-Z]+)"[^>]*\bdbkey="([^"]+)"')
_XML_ABSTRACT_OPEN = "<abstract"
_XML_ABSTRACT_CLOSE = "</abstract>"

# clean_abstract() building blocks. Real InterPro abstracts pretty-print each
# citation as `[` <newline> `<cite idref="PUB..."/>` <newline> `]` (one or more
# <cite> tags, comma-separated) — after tag stripping the brackets are empty
# and must be swept up too. A handful of entries additionally carry a literal
# `[cite:...]` text marker (not a tag) left over from InterPro's authoring
# tooling; that is swept separately since it never matches a real XML tag.
_TAG_RE = re.compile(r"<[^>]+>")
_CITE_MARKER_RE = re.compile(r"\[cite:[^\]]*\]", re.IGNORECASE)
_EMPTY_BRACKET_RE = re.compile(r"\[\s*(?:,\s*)*\]")
_WHITESPACE_RE = re.compile(r"\s+")
_SPACE_BEFORE_PUNCT_RE = re.compile(r"\s+([.,;:)])")

# InterPro's uppercase `db` token → the `databaseName` casing InterProScan writes
# into calls.json, so both artifacts speak one vocabulary (`MetaCyc:PWY-1042`).
PATHWAY_DB_NAMES = {"METACYC": "MetaCyc", "REACTOME": "Reactome"}

# Reactome xrefs are species-expanded (one per Reactome pathway per organism,
# ~507K rows) and are near-pure noise for marine bacteria — MetaCyc only.
DEFAULT_PATHWAY_DBS = ("METACYC",)


def normalize_type(raw_type: str | None) -> str | None:
    """Uppercase an entry.list ENTRY_TYPE to the calls.json ``interpro_type`` form."""
    if not raw_type:
        return None
    return raw_type.strip().upper().replace(" ", "_")


def parse_entry_list(text: str) -> dict[str, dict[str, str]]:
    """Parse ``entry.list`` → ``{IPRxxxxxx: {"type", "name"}}``.

    Skips the header line (``ENTRY_AC``) and any malformed rows. ``type`` is the
    uppercased/normalized entry type; ``name`` is the entry description.
    """
    out: dict[str, dict[str, str]] = {}
    for line in text.splitlines():
        line = line.rstrip("\n")
        if not line or line.startswith("ENTRY_AC"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        acc = parts[0].strip()
        if not acc.startswith("IPR"):
            continue
        entry_type = normalize_type(parts[1])
        name = parts[2].strip()
        out[acc] = {"type": entry_type or "", "name": name}
    return out


def parse_parent_child_tree(text: str) -> dict[str, dict]:
    """Parse ``ParentChildTreeFile.txt`` → ``{IPRxxxxxx: {"parent", "level"}}``.

    Depth = number of leading ``--`` pairs. A line's parent is the most recent
    line seen at ``depth - 1``. Only entries that appear in the tree are returned;
    the caller defaults every other entry to ``level 0`` / ``parent None``.
    """
    out: dict[str, dict] = {}
    # depth → the accession most recently opened at that depth (the running path)
    path: dict[int, str] = {}
    for line in text.splitlines():
        m = _TREE_LINE_RE.match(line)
        if not m:
            continue
        depth = len(m.group("dashes")) // 2
        acc = m.group("acc")
        parent = path.get(depth - 1) if depth > 0 else None
        out[acc] = {"parent": parent, "level": depth}
        path[depth] = acc
        # invalidate any deeper path entries — we've moved to a new subtree
        for d in [d for d in path if d > depth]:
            del path[d]
    return out


def parse_interpro2go(text: str) -> dict[str, list[str]]:
    """Parse ``interpro2go`` → ``{IPRxxxxxx: [GO:0003677, …]}`` (sorted, deduped).

    Comment lines (``!``) and anything not matching the
    ``InterPro:<acc> … ; GO:<id>`` shape are skipped. One entry maps to many GO
    terms, each on its own line.
    """
    acc_to_go: dict[str, set[str]] = {}
    for line in text.splitlines():
        m = _IPR2GO_RE.match(line.rstrip("\n"))
        if not m:
            continue
        acc_to_go.setdefault(m.group("acc"), set()).add(m.group("go"))
    return {acc: sorted(gos) for acc, gos in acc_to_go.items()}


def parse_entry_db_xrefs(
    lines: Iterable[str],
    include_dbs: Iterable[str],
) -> dict[str, dict[str, list[str]]]:
    """Stream ``interpro.xml`` lines → ``{IPRxxxxxx: {DB: [raw dbkey, …]}}``.

    Generic entry-level ``<db_xref>`` extractor. Takes an *iterable of lines* (not
    a string) so the caller can stream the 42 MB gzipped XML without decompressing
    it into memory. Only ``db`` tokens in *include_dbs* (upper-cased) are kept;
    raw ``dbkey`` values are returned **verbatim** (no label prefix), grouped by
    the upper-case ``db`` token. This is the single XML pass that feeds pathways
    (via :func:`parse_pathway_xrefs`), EC numbers (``db="EC"``), and CAZy families
    (``db="CAZY"``) — the two latter as raw ids (``1.1.1.1`` / ``GH13``).

    Cross-references are attributed to the most recently opened ``<interpro
    id="…">`` element and dropped outside any entry, so header/footer blocks
    cannot leak in. EC/CAZy are entry-level metadata; member-signature db_xrefs
    (PFAM, PROSITE, …) inside ``<member_list>`` are simply not in *include_dbs*.
    """
    wanted = {db.upper() for db in include_dbs}
    acc_to_dbs: dict[str, dict[str, set[str]]] = {}
    current: str | None = None
    for line in lines:
        m = _XML_ENTRY_OPEN_RE.search(line)
        if m:
            current = m.group(1)
        if current is not None:
            for db, key in _XML_DB_XREF_RE.findall(line):
                if db in wanted:
                    acc_to_dbs.setdefault(current, {}).setdefault(db, set()).add(key)
        if _XML_ENTRY_CLOSE in line:
            current = None
    return {
        acc: {db: sorted(keys) for db, keys in dbs.items()}
        for acc, dbs in acc_to_dbs.items()
    }


def parse_pathway_xrefs(
    lines: Iterable[str],
    include_dbs: Iterable[str] = DEFAULT_PATHWAY_DBS,
) -> dict[str, list[str]]:
    """Stream ``interpro.xml`` lines → ``{IPRxxxxxx: ["MetaCyc:PWY-1042", …]}``.

    Thin wrapper over :func:`parse_entry_db_xrefs` that renders each kept xref as
    ``{DatabaseName}:{key}`` via :data:`PATHWAY_DB_NAMES`, matching the ``DB:id``
    form InterProScan writes into calls.json.
    """
    raw = parse_entry_db_xrefs(lines, include_dbs)
    out: dict[str, list[str]] = {}
    for acc, dbs in raw.items():
        pws: set[str] = set()
        for db, keys in dbs.items():
            label = PATHWAY_DB_NAMES.get(db, db)
            pws.update(f"{label}:{key}" for key in keys)
        if pws:
            out[acc] = sorted(pws)
    return out


def clean_abstract(html_text: str, cap: int = 400) -> str:
    """First paragraph of an InterPro entry ``<abstract>``, plain text, capped.

    Takes the text up to the first ``</p>`` (or the whole string if there is
    none), strips every XML/HTML tag (``<p>``, ``<i>``, ``<sub>``,
    ``<cite idref="..."/>`` …), unescapes HTML entities, sweeps up the empty
    ``[ ]``/``[ , ]`` citation brackets left behind once ``<cite>`` tags are
    removed plus the rare literal ``[cite:...]`` text marker, collapses
    whitespace, tidies stray whitespace before punctuation, and truncates to
    *cap* characters.
    """
    text = html_text.split("</p>", 1)[0]
    text = _TAG_RE.sub("", text)
    text = html.unescape(text)
    text = _CITE_MARKER_RE.sub("", text)
    text = _EMPTY_BRACKET_RE.sub("", text)
    text = _WHITESPACE_RE.sub(" ", text).strip()
    text = _SPACE_BEFORE_PUNCT_RE.sub(r"\1", text)
    return text[:cap]


def parse_entry_descriptions(lines: Iterable[str]) -> dict[str, str]:
    """Stream ``interpro.xml`` lines → ``{IPRxxxxxx: description}``.

    Captures the first ``<abstract>`` block per entry (attributed the same
    way :func:`parse_entry_db_xrefs` attributes ``<db_xref>`` elements — to
    the most recently opened ``<interpro id="…">``), stops accumulating as
    soon as the first ``</p>`` (or ``</abstract>`` for a paragraph-less
    abstract) is seen, and runs :func:`clean_abstract` over the accumulated
    raw markup. Entries with no ``<abstract>`` or an empty one are omitted —
    sparse, matching ``go_terms``/``pathways``/``ec_numbers``/``cazy_ids``.
    """
    out: dict[str, str] = {}
    current: str | None = None
    in_abstract = False
    buf: list[str] = []
    for line in lines:
        m = _XML_ENTRY_OPEN_RE.search(line)
        if m:
            current = m.group(1)
            in_abstract = False
            buf = []
        if current is not None and current not in out:
            if not in_abstract and _XML_ABSTRACT_OPEN in line:
                in_abstract = True
            if in_abstract:
                buf.append(line)
                if "</p>" in line or _XML_ABSTRACT_CLOSE in line:
                    desc = clean_abstract("".join(buf))
                    if desc:
                        out[current] = desc
                    in_abstract = False
        if _XML_ENTRY_CLOSE in line:
            current = None
    return out


def build_reference(
    entry_list_text: str,
    tree_text: str,
    go_map: dict[str, list[str]] | None = None,
    pathway_map: dict[str, list[str]] | None = None,
    ec_map: dict[str, list[str]] | None = None,
    cazy_map: dict[str, list[str]] | None = None,
    description_map: dict[str, str] | None = None,
) -> dict[str, dict]:
    """Combine the release files into the committed reference dict.

    ``{IPRxxxxxx: {"name", "type", "parent", "level"[, "go_terms"][, "pathways"]
    [, "ec_numbers"][, "cazy_ids"][, "description"]}}`` — name/type from
    entry.list; parent/level from the tree (``None``/``0`` when not in the
    tree); ``go_terms`` from *go_map* (``interpro2go``); ``pathways`` /
    ``ec_numbers`` / ``cazy_ids`` / ``description`` from the ``interpro.xml``
    stream (*pathway_map* / *ec_map* / *cazy_map* / *description_map*).

    ``ec_numbers`` / ``cazy_ids`` are stored **raw** (``1.1.1.1`` / ``GH13``) —
    ``normalize_ec`` and bare-3-level EC normalization are applied by the
    consumer (Phase-2 gene-annotation enrichment), not here, so the reference
    stays a faithful copy of InterPro. ``description`` is the entry's first
    abstract paragraph, plain text, capped to 400 chars (see
    :func:`clean_abstract` / :func:`parse_entry_descriptions`).

    All five xref/description fields are **sparse** — the key is omitted
    entirely for entries with none, rather than carrying an empty
    list/string. Consumers must use ``meta.get("ec_numbers", [])`` /
    ``meta.get("description", "")``.

    Entries present only in the tree (should not happen, but guarded) get an
    empty name/type.
    """
    entries = parse_entry_list(entry_list_text)
    tree = parse_parent_child_tree(tree_text)
    go_map = go_map or {}
    pathway_map = pathway_map or {}
    ec_map = ec_map or {}
    cazy_map = cazy_map or {}
    description_map = description_map or {}

    ref: dict[str, dict] = {}
    for acc, meta in entries.items():
        t = tree.get(acc)
        ref[acc] = {
            "name": meta["name"],
            "type": meta["type"],
            "parent": t["parent"] if t else None,
            "level": t["level"] if t else 0,
        }
    # Tree-only accessions (defensive): keep them with empty name/type so the
    # hierarchy never references a missing node.
    for acc, t in tree.items():
        if acc not in ref:
            ref[acc] = {"name": "", "type": "", "parent": t["parent"], "level": t["level"]}

    # Sparse xref fields, applied only to entries that exist in the reference —
    # a GO/pathway mapping for a retired accession is dropped, never resurrected
    # as a nameless node.
    for acc, gos in go_map.items():
        if acc in ref and gos:
            ref[acc]["go_terms"] = sorted(set(gos))
    for acc, pws in pathway_map.items():
        if acc in ref and pws:
            ref[acc]["pathways"] = sorted(set(pws))
    for acc, ecs in ec_map.items():
        if acc in ref and ecs:
            ref[acc]["ec_numbers"] = sorted(set(ecs))
    for acc, czs in cazy_map.items():
        if acc in ref and czs:
            ref[acc]["cazy_ids"] = sorted(set(czs))
    for acc, desc in description_map.items():
        if acc in ref and desc:
            ref[acc]["description"] = desc
    return ref
