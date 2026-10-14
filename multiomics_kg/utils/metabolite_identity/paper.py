"""Paper-metabolite resolution (spec §3.6, §3.6a, §6.5, §6.6). Stages D (harvest) and J."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd
import yaml

from multiomics_kg.utils.metabolite_identity.names import is_stereo_qualified, name_key
from multiomics_kg.utils.paperconfig_utils import iter_metabolite_assays_tables

STEREO_ASSIGNMENTS = ("as_reported", "dominant_form_assumed")
_NO_ID = {"0", "-", "na", "n/a", "none"}
_ID_METHOD = {"kegg.compound": "kegg_direct", "chebi": "chebi_direct", "mnx": "mnx_direct"}


@dataclass(frozen=True)
class Override:
    id: str
    stereo_assignment: str | None
    reason: str
    raw_name: str


@dataclass
class PaperRow:
    paper: str
    table: str
    row_index: int
    name: str
    direct_id: str | None
    id_method: str | None

    @property
    def name_key(self) -> str:
        return name_key(self.name)


@dataclass
class PaperTable:
    paper: str
    table: str
    csv_path: Path | None
    skip_rows: int
    rows: list[PaperRow]
    name_col: str = "compound"

    @property
    def key(self) -> str:
        return f"{self.paper}|{self.table}"


def load_paper_overrides(pc_path: Path, cfg: dict) -> dict[str, Override]:
    files = sorted({e.get("aliases_file") for _, e in iter_metabolite_assays_tables(cfg)
                    if e.get("aliases_file")})
    out: dict[str, Override] = {}
    for f in files:
        path = pc_path.parent / f
        if not path.exists():
            raise FileNotFoundError(f"{pc_path}: aliases_file {f!r} not found at {path}")
        doc = yaml.safe_load(path.read_text())
        if doc is None:
            doc = {}
        if not isinstance(doc, dict):
            raise ValueError(f"{path}: alias file must be a mapping of name -> id, "
                             f"got {type(doc).__name__}")
        for raw, val in doc.items():
            if not isinstance(raw, str):
                # YAML 1.1 reads bare NO / No / Yes / On / Off keys as booleans (and digits as ints).
                raise ValueError(f"{path}: alias key {raw!r} is not a string; quote the key in the "
                                 f"yaml (e.g. \"NO\": ...)")
            if isinstance(val, dict):
                if not val.get("id"):
                    raise ValueError(f"{path}: alias {raw!r} has no 'id'")
                ov = Override(str(val["id"]), val.get("stereo_assignment"), str(val.get("reason", "")), str(raw))
            else:
                ov = Override(str(val), None, "", str(raw))
            if ov.stereo_assignment not in (None, *STEREO_ASSIGNMENTS):
                raise ValueError(f"{path}: {raw!r} stereo_assignment {ov.stereo_assignment!r} "
                                 f"not in {STEREO_ASSIGNMENTS}")
            k = name_key(str(raw))
            if k in out and out[k].id != ov.id:
                raise ValueError(f"{path}: {raw!r} and {out[k].raw_name!r} share name_key {k!r} "
                                 f"but point to {ov.id} vs {out[k].id}")
            out[k] = ov
    return out


def harvest_paper_tables(paperconfigs, project_root: Path):
    tables, overrides = [], {}
    for pc_path, cfg in paperconfigs:
        try:
            paper = str(Path(pc_path).resolve().relative_to(project_root.resolve()))
        except ValueError:
            paper = str(pc_path)
        overrides[paper] = load_paper_overrides(Path(pc_path), cfg)
        for entry_key, entry in iter_metabolite_assays_tables(cfg):
            csv = entry.get("filename")
            if not csv:
                continue
            csv_full = Path(csv) if Path(csv).is_absolute() else project_root / csv
            skip = int(entry.get("skip_rows", 0) or 0)
            df = pd.read_csv(csv_full, dtype=str, keep_default_na=False, skiprows=skip or None)
            name_col = entry.get("name_col") or "compound"
            id_col, id_type = entry.get("id_col") or "", entry.get("id_type") or ""
            for col in (name_col, id_col):
                if col and col not in df.columns:
                    raise ValueError(f"{pc_path} [{entry_key}]: column {col!r} not in {csv_full} "
                                     f"(columns: {list(df.columns)})")
            rows = []
            for i, row in df.iterrows():
                name = (row.get(name_col) or "").strip()
                if not name:
                    continue
                direct, method = None, None
                if id_col and id_type:
                    raw = (row.get(id_col) or "").strip()
                    if raw and raw.lower() not in _NO_ID:
                        direct = raw if ":" in raw else f"{id_type}:{raw}"
                        method = _ID_METHOD.get(id_type, "id_direct")
                rows.append(PaperRow(paper, entry_key, int(i), name, direct, method))
            tables.append(PaperTable(paper, entry_key, csv_full, skip, rows, name_col))
    return tables, overrides


@dataclass
class Decision:
    metabolite_id: str | None
    resolution_method: str
    stereo_assignment: str | None


@dataclass
class StereoContext:
    families: dict[str, frozenset[str]]
    unspecified: set[str]
    is_anomer: Callable[[str, str], bool]
    catalyst_counts: dict[str, int]

    def variants(self, node: str) -> frozenset[str]:
        fam = self.families.get(node, frozenset({node}))
        return frozenset(m for m in fam if m == node or not self.is_anomer(m, node))


def collapse_anomers(candidates, ctx: StereoContext) -> set[str]:
    """Anomer guard (spec §3.6): drop a candidate that is only an anomer of another candidate
    that is the anomer-unspecified (stereo-unspecified) node, e.g. beta-GlcNAc next to GlcNAc."""
    cands = set(candidates)
    return {m for m in cands
            if not any(u != m and u in ctx.unspecified and ctx.is_anomer(m, u) for u in cands)}


def dominant_form(candidates, ctx: StereoContext):
    """Dominant-stereo-form rule (spec §3.6): the family member with the most catalysts.

    The pool is every candidate's stereo family (unspecified node included), minus members
    that are only anomers of a candidate (anomer guard). Refuses (winner ``None``) when the
    winner has 0 catalysts, ties the runner-up, or has < 2x the runner-up's count; the info
    dict is returned either way so drift / ambiguity reports can show the counts.
    """
    candidates = collapse_anomers(candidates, ctx)
    pool = set()
    for c in candidates:
        pool |= ctx.variants(c)
    pool = {m for m in pool if m in candidates or not any(ctx.is_anomer(m, c) for c in candidates)}
    if not pool:
        return None, {"winner": None, "counts": {"winner": 0, "runner_up": 0}, "runner_up_id": None}
    ranked = sorted(pool, key=lambda m: (-ctx.catalyst_counts.get(m, 0), m))
    win = ranked[0]
    run = ranked[1] if len(ranked) > 1 else None
    wc, rc = ctx.catalyst_counts.get(win, 0), ctx.catalyst_counts.get(run, 0) if run else 0
    info = {"winner": win, "counts": {"winner": wc, "runner_up": rc}, "runner_up_id": run}
    ok = wc > 0 and wc >= 2 * rc and (run is None or wc != rc)
    return (win if ok else None), info


@dataclass
class ResolutionResult:
    decisions: dict = field(default_factory=dict)
    lock: dict = field(default_factory=dict)
    drift: list = field(default_factory=list)
    ambiguous: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)
    collisions: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    stale_overrides: list = field(default_factory=list)
    stale_lock: list = field(default_factory=list)
    lock_remapped: list = field(default_factory=list)
    name_lookup: dict = field(default_factory=dict)


def resolve_paper_rows(tables, *, overrides, name_candidates, id_index, ctx, lock, today):
    """Stage J: resolve every paper row in the spec §3.6 order.

    Per row: (1) the paper's alias override, (2) the direct id (through ``id_index``),
    (3) the frozen lock, (4) the name match / dominant-form rule. Nothing here raises:
    ambiguous names, unresolvable override targets and legacy overrides on stereo targets are
    collected in ``ambiguous`` / ``errors`` for stage P to fail on. Stale lock entries stay in
    the returned lock (a human deletes them); new freezes are added with ``frozen_on=today``.
    """
    res = ResolutionResult()
    lock = {"names": {k: dict(v) for k, v in sorted(lock.get("names", {}).items())},
            "direct_ids": {k: dict(v) for k, v in sorted(lock.get("direct_ids", {}).items())}}
    # Lock keys whose entry is stale: kept in the lock for a human to delete (spec §6.6),
    # never applied, never overwritten.
    dead: dict[str, set[str]] = {"names": set(), "direct_ids": set()}
    # Re-map lock keys / targets through the current aliases (spec §3.6a).
    for section in ("names", "direct_ids"):
        for key in list(lock[section]):
            entry = lock[section][key]
            tgt = id_index.get(entry["id"])
            if tgt is None:
                res.stale_lock.append({"section": section, "key": key, "id": entry["id"],
                                       "why": "target gone"})
                dead[section].add(key)
                continue
            if tgt != entry["id"]:
                res.lock_remapped.append({"section": section, "key": key, "what": "target",
                                          "old": entry["id"], "new": tgt})
                entry["id"] = tgt
            new_key = id_index.get(key, key) if section == "direct_ids" else key
            if new_key != key:
                if new_key in lock[section]:
                    # Two old keys merged into one canonical id: keep the surviving one,
                    # leave the other where it is for a human.
                    res.stale_lock.append({"section": section, "key": key, "id": entry["id"],
                                           "why": f"key merged into {new_key}, already locked"})
                    dead[section].add(key)
                    continue
                res.lock_remapped.append({"section": section, "key": key, "what": "key",
                                          "old": key, "new": new_key})
                del lock[section][key]
                lock[section][new_key] = entry
    used_overrides: set[tuple[str, str]] = set()

    def live(section, key):
        return key in lock[section] and key not in dead[section]

    def has_variants(node):
        return len(ctx.variants(node)) > 1

    def error(msg):
        if msg not in res.errors:
            res.errors.append(msg)

    def rule(section, key, candidates):
        win, info = dominant_form(candidates, ctx)
        if live(section, key):
            entry = lock[section][key]
            # win is None when the margin fell below 2x / tie / no catalysis (dominant_form refuses).
            if win != entry["id"]:
                drift = {"section": section, "key": key, "frozen": entry["id"],
                         "recomputed": win, "counts": info["counts"]}
                if drift not in res.drift:
                    res.drift.append(drift)
            return entry["id"]
        if win is None:
            return None
        if key not in dead[section]:
            lock[section][key] = {"id": win, "frozen_on": today, "counts": info["counts"],
                                  "runner_up_id": info["runner_up_id"]}
        return win

    for t in tables:
        for row in t.rows:
            k = row.name_key
            key3 = (row.paper, row.table, row.row_index)
            ov = overrides.get(row.paper, {}).get(k)
            if ov is not None:
                used_overrides.add((row.paper, k))
                tgt = id_index.get(ov.id)
                if tgt is None:
                    error(f"{row.paper}: alias {ov.raw_name!r} → {ov.id} resolves to no compound")
                    continue
                sa = ov.stereo_assignment
                if sa is None:
                    if has_variants(tgt):
                        error(f"{row.paper}: alias {ov.raw_name!r} → {tgt} has stereo variants; "
                              f"add an explicit stereo_assignment")
                        continue
                    sa = "as_reported"
                res.decisions[key3] = Decision(tgt, "alias_override", sa)
                continue
            if row.direct_id:
                tgt = id_index.get(row.direct_id)
                if tgt is None:
                    res.unresolved.append({"paper": row.paper, "table": row.table, "name": row.name,
                                           "direct_id": row.direct_id})
                    continue
                if tgt in ctx.unspecified and has_variants(tgt):
                    chosen = rule("direct_ids", tgt, {tgt})
                    if chosen is None:
                        res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                              "candidates": sorted(ctx.variants(tgt))})
                        continue
                    method = "dominant_stereo_form" if chosen != tgt else row.id_method
                    res.decisions[key3] = Decision(chosen, method, "dominant_form_assumed")
                else:
                    res.decisions[key3] = Decision(tgt, row.id_method, "as_reported")
                continue
            cands = collapse_anomers({id_index.get(c, c) for c in name_candidates.get(k, set())}, ctx)
            qualified = is_stereo_qualified(row.name)
            if live("names", k):
                chosen = rule("names", k, cands or {lock["names"][k]["id"]})
                # Same method rule as the freezing path below, so build 1 and build 2 agree.
                method = ("name_match" if len(cands) == 1 and chosen in cands
                          else "dominant_stereo_form")
                res.decisions[key3] = Decision(chosen, method, "dominant_form_assumed")
                continue
            if not cands:
                res.unresolved.append({"paper": row.paper, "table": row.table, "name": row.name})
                continue
            fams = {ctx.families.get(c, frozenset({c})) for c in cands}
            if len(cands) == 1:
                (c,) = cands
                if not has_variants(c) or qualified:
                    res.decisions[key3] = Decision(c, "name_match", "as_reported")
                elif c in ctx.unspecified:
                    chosen = rule("names", k, {c})
                    if chosen is None:
                        res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                              "candidates": sorted(ctx.variants(c))})
                        continue
                    method = "dominant_stereo_form" if chosen != c else "name_match"
                    res.decisions[key3] = Decision(chosen, method, "dominant_form_assumed")
                else:
                    res.decisions[key3] = Decision(c, "name_match", "dominant_form_assumed")
                continue
            if len(fams) == 1 and not qualified:
                chosen = rule("names", k, cands)
                if chosen is not None:
                    res.decisions[key3] = Decision(chosen, "dominant_stereo_form", "dominant_form_assumed")
                    continue
            res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                  "candidates": sorted(cands)})

    # Same-table collisions (spec §3.6): distinct name_keys on one node need both overridden.
    for t in tables:
        by_node: dict[str, set[str]] = defaultdict(set)
        for row in t.rows:
            d = res.decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                by_node[d.metabolite_id].add(row.name_key)
        for node, keys in sorted(by_node.items()):
            if len(keys) < 2:
                continue
            ovs = overrides.get(t.paper, {})
            if all(k in ovs and id_index.get(ovs[k].id) == node for k in keys):
                continue
            res.collisions.append({"paper": t.paper, "table": t.table, "node": node, "names": sorted(keys)})

    seen_names = {row.name_key for t in tables for row in t.rows}
    seen_ids = {id_index.get(row.direct_id) for t in tables for row in t.rows if row.direct_id}
    for k in sorted(set(lock["names"]) - seen_names - dead["names"]):
        res.stale_lock.append({"section": "names", "key": k, "id": lock["names"][k]["id"], "why": "no row"})
    for k in sorted(set(lock["direct_ids"]) - seen_ids - dead["direct_ids"]):
        res.stale_lock.append({"section": "direct_ids", "key": k, "id": lock["direct_ids"][k]["id"], "why": "no row"})
    for paper, ovs in sorted(overrides.items()):
        for k, ov in sorted(ovs.items()):
            if (paper, k) not in used_overrides:
                res.stale_overrides.append({"paper": paper, "name_key": k, "id": ov.id})
    overridden = {k for ovs in overrides.values() for k in ovs}
    for k, cands in sorted(name_candidates.items()):
        if k in overridden:
            continue
        canon = collapse_anomers({id_index.get(c, c) for c in cands}, ctx)
        if len(canon) == 1:
            res.name_lookup[k] = next(iter(canon))
    res.lock = lock
    return res
