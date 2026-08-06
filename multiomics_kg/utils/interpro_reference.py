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

Combined result: ``{IPRxxxxxx: {"name", "type", "parent", "level"}}`` — ``name``
and ``type`` from entry.list (authoritative, untruncated), ``parent``/``level``
from the tree (``None``/``0`` when the entry is not in the tree).

See ``docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md``.
"""

from __future__ import annotations

import re

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


def build_reference(entry_list_text: str, tree_text: str) -> dict[str, dict]:
    """Combine the two files into the committed reference dict.

    ``{IPRxxxxxx: {"name", "type", "parent", "level"}}`` — name/type from
    entry.list; parent/level from the tree (``None``/``0`` when not in the tree).
    Entries present only in the tree (should not happen, but guarded) get an
    empty name/type.
    """
    entries = parse_entry_list(entry_list_text)
    tree = parse_parent_child_tree(tree_text)

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
    return ref
