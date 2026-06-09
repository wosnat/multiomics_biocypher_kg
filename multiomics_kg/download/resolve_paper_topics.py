#!/usr/bin/env python3
"""Resolve discussed-topic mentions to canonical Gene / KEGG-pathway nodes.

Stage 2 of the discuss-edges pipeline (extract -> resolve -> adapter), wired as
``prepare_data.sh`` step 8. For each paper with a ``publication_topics/topics.json``
(Stage 1, written by the manual /extract-discussed-topics LLM step), resolves:

  - GENE mentions to strain-specific locus tags via the per-strain
    ``gene_id_mapping`` (identifiers first via Tier-1 ``specific_lookup``,
    gene_name as a Tier-3 singleton fallback). Strain-aware: a mention
    attributed to a specific strain resolves there; an ``all`` / ``unspecified``
    mention resolves against every paper strain and emits one record per strain
    that resolves (bounded to the paper's strains, never all 23).
  - KEGG PATHWAY mentions to ``ko*`` pathway ids via a GLOBAL lookup built from
    ``kegg_data.json['pathways']`` — the exact set the kegg adapter turns into
    ``KeggTerm`` pathway nodes, so anything that resolves is guaranteed to be a
    real node (no dangling edges by construction).

Writes ``publication_topics/topics_resolved.json`` + ``resolution_report.txt`` in the
paper's subfolder. The Stage 3 adapter reads only the resolved ids and emits edges.

Usage:
  uv run python -m multiomics_kg.download.resolve_paper_topics [--force] [--papers NAME ...]
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path

from multiomics_kg.download.utils.paths import PROJECT_ROOT
from multiomics_kg.extraction.topics.extraction_utils import (
    find_all_papers,
    has_topics,
    load_topics,
    migrate_legacy,
    report_path,
    resolved_path,
    topics_path,
)
from multiomics_kg.utils.gene_id_utils import (
    MappingData,
    get_genome_dir,
    load_mapping_v2,
    resolve_row,
)

logger = logging.getLogger(__name__)

DEFAULT_KEGG_DATA = PROJECT_ROOT / "cache" / "data" / "kegg" / "kegg_data.json"


# ─── KEGG pathway lookup (global) ───────────────────────────────────────────────


def normalize_kegg_pathway_id(raw: str | None) -> str | None:
    """Normalize an LLM-supplied KEGG pathway id to the ``koNNNNN`` key form.

    Accepts ``ko00910`` / ``map00910`` / ``path:map00910`` / ``00910`` — the
    prefix is irrelevant since pathway maps share the numeric id across ko/map.
    Returns None if no 4-5 digit pathway number is present.
    """
    if not raw:
        return None
    m = re.search(r"(\d{4,5})", str(raw))
    if not m:
        return None
    return "ko" + m.group(1).zfill(5)


def build_pathway_lookup(kegg_data: dict) -> tuple[dict[str, str], dict[str, str]]:
    """Build (id_lookup, name_lookup) from kegg_data.json['pathways'].

    id_lookup: normalized koNNNNN -> koNNNNN (the canonical key, == node id base)
    name_lookup: lowercased pathway name -> koNNNNN
    Both keyed only on pathways that will become KeggTerm nodes.
    """
    pathways = kegg_data.get("pathways", {})
    id_lookup: dict[str, str] = {}
    name_lookup: dict[str, str] = {}
    for pw_id, info in pathways.items():
        norm = normalize_kegg_pathway_id(pw_id)
        if norm:
            id_lookup[norm] = pw_id
        name = (info or {}).get("name", "")
        if name:
            name_lookup[name.strip().lower()] = pw_id
    return id_lookup, name_lookup


def resolve_pathway(
    mention: dict,
    id_lookup: dict[str, str],
    name_lookup: dict[str, str],
) -> tuple[str | None, str | None, str]:
    """Resolve one pathway mention. Returns (ko_key | None, name | None, method)."""
    # id-first
    norm = normalize_kegg_pathway_id(mention.get("kegg_pathway_id"))
    if norm and norm in id_lookup:
        ko = id_lookup[norm]
        return ko, _pathway_name(ko, name_lookup), "id"
    # normalized-name fallback
    surface = (mention.get("surface_form") or "").strip().lower()
    if surface and surface in name_lookup:
        ko = name_lookup[surface]
        return ko, _pathway_name(ko, name_lookup), "name"
    return None, None, "unresolved"


def _pathway_name(ko: str, name_lookup: dict[str, str]) -> str | None:
    for name, k in name_lookup.items():
        if k == ko:
            return name
    return None


# ─── Gene mention resolution (strain-aware) ─────────────────────────────────────


def _load_strain_mapping(organism: str, cache: dict[str, MappingData | None]) -> MappingData | None:
    """Load (and cache) the MappingData for an organism, or None if unavailable."""
    if organism in cache:
        return cache[organism]
    md: MappingData | None = None
    genome_dir = get_genome_dir(organism, str(PROJECT_ROOT))
    if genome_dir:
        md = load_mapping_v2(genome_dir)
        if md is None:
            # Legacy fallback (mirrors resolve_paper_ids).
            from multiomics_kg.utils.gene_id_utils import build_id_lookup
            lookup, locus_tags, _ = build_id_lookup(genome_dir)
            if lookup is not None:
                md = MappingData(
                    specific_lookup=lookup, multi_lookup={}, conflicts={},
                    locus_tags=locus_tags or set(), version=0,
                )
    cache[organism] = md
    return md


def _resolve_value(value: str, mapping_data: MappingData) -> tuple[str | None, str]:
    """Resolve a single identifier/name value to a locus_tag via resolve_row."""
    return resolve_row({"v": value}, "v", [], mapping_data)


def _resolve_in_strain(mention: dict, mapping_data: MappingData) -> list[tuple[str, str, str]]:
    """Resolve a gene mention against one strain's mapping.

    Resolves EACH identifier independently and unions the distinct locus tags —
    so a gene family with N member ids (e.g. CCRG-2 → 12 PMT* tags) yields N
    targets, while synonyms of the same gene collapse to one. gene_name (or
    surface form) is a fallback only when NO identifier resolved.

    Returns a list of (locus_tag, method, source_value), distinct by locus_tag.
    """
    identifiers = [str(i).strip() for i in (mention.get("identifiers") or []) if str(i).strip()]
    resolved: dict[str, tuple[str, str]] = {}  # locus_tag -> (method, source)

    for ident in identifiers:
        lt, method = _resolve_value(ident, mapping_data)
        if lt and lt not in resolved:
            resolved[lt] = (method, ident)

    if not resolved:
        name_val = (mention.get("gene_name") or "").strip() or (mention.get("surface_form") or "").strip()
        if name_val:
            lt, method = _resolve_value(name_val, mapping_data)
            if lt:
                resolved[lt] = (method, name_val)

    return [(lt, m, src) for lt, (m, src) in resolved.items()]


def _match_strain(mention_strain: str, paper_strains: list[str]) -> str | None:
    """Match a mention's strain to one of the paper's strains (exact, then ci)."""
    if not mention_strain:
        return None
    for s in paper_strains:
        if s == mention_strain:
            return s
    low = mention_strain.lower()
    for s in paper_strains:
        if s.lower() == low:
            return s
    return None


def resolve_gene_mention(
    mention: dict,
    paper_strains: list[str],
    mapping_cache: dict[str, MappingData | None],
) -> list[dict]:
    """Resolve a gene mention into one or more strain-specific records.

    - specific strain -> one record (resolved or unresolved)
    - 'all' / 'unspecified' / unknown -> resolve in each paper strain; emit one
      record per strain that RESOLVES. If none resolve anywhere, emit a single
      unresolved record (for audit) with the original strain label.
    """
    raw_strain = (mention.get("strain") or "").strip()
    matched = _match_strain(raw_strain, paper_strains)

    def _record(strain: str, locus_tag: str | None, method: str, source: str = "") -> dict:
        return {
            "surface_form": mention.get("surface_form", ""),
            "gene_name": mention.get("gene_name", ""),
            "identifiers": mention.get("identifiers", []),
            "prominence": mention.get("prominence", "peripheral"),
            "evidence": mention.get("evidence", ""),   # carried onto the edge (Stage 3)
            "mention_strain": raw_strain,
            "strain": strain,
            "locus_tag": locus_tag,
            "resolved_from": source,                   # audit only (not on the edge)
            "resolution_method": method,               # audit only (not on the edge)
        }

    # Specific, recognized strain → resolve there only.
    if matched:
        md = _load_strain_mapping(matched, mapping_cache)
        if md is None:
            return [_record(matched, None, "no_mapping")]
        hits = _resolve_in_strain(mention, md)
        if hits:
            return [_record(matched, lt, m, src) for lt, m, src in hits]
        return [_record(matched, None, "unresolved")]

    # Generic / unknown → try every paper strain.
    records: list[dict] = []
    for strain in paper_strains:
        md = _load_strain_mapping(strain, mapping_cache)
        if md is None:
            continue
        for lt, m, src in _resolve_in_strain(mention, md):
            records.append(_record(strain, lt, m, src))
    if records:
        return records
    # Nothing resolved anywhere — keep one unresolved record for audit.
    return [_record(raw_strain or "all", None, "unresolved")]


# ─── Truncation quantification (best-effort visibility, not a fix) ──────────────

# A "truncated" id is a bare locus-tag suffix with no strain-letter prefix —
# e.g. "0246", "0994", "0246*". These cannot resolve without the strain prefix.
_TRUNCATED_ID_RE = re.compile(r"^\d{1,6}[*+]?$")


def is_truncated_id(value: str) -> bool:
    return bool(_TRUNCATED_ID_RE.match(str(value).strip()))


def _unresolved_reason(rec: dict) -> str:
    """Classify why an unresolved gene record failed (for the report)."""
    ids = [str(i).strip() for i in (rec.get("identifiers") or []) if str(i).strip()]
    has_name = bool((rec.get("gene_name") or "").strip())
    if not ids and not has_name:
        return "descriptive_only"   # surface form is a function/enzyme phrase
    if ids and all(is_truncated_id(i) for i in ids) and not has_name:
        return "truncated_id_only"  # recoverable in principle via strain prefix
    return "lookup_miss"            # had a name / full id but no mapping match


# ─── Per-paper driver ───────────────────────────────────────────────────────────


def resolve_paper_topics(
    paper_dir: Path,
    paper_strains: list[str],
    id_lookup: dict[str, str],
    name_lookup: dict[str, str],
    mapping_cache: dict[str, MappingData | None],
) -> dict | None:
    """Resolve one paper's topics. Returns a stats dict, or None if no topics."""
    migrate_legacy(paper_dir)
    data = load_topics(paper_dir)
    if not data:
        return None

    raw_genes = data.get("genes", [])

    gene_records: list[dict] = []
    gene_methods: Counter = Counter()
    for m in raw_genes:
        for rec in resolve_gene_mention(m, paper_strains, mapping_cache):
            gene_records.append(rec)
            gene_methods[rec["resolution_method"]] += 1

    pathway_records: list[dict] = []
    pathway_methods: Counter = Counter()
    for m in data.get("pathways", []):
        ko, name, method = resolve_pathway(m, id_lookup, name_lookup)
        pathway_records.append({
            "surface_form": m.get("surface_form", ""),
            "kegg_pathway_id_raw": m.get("kegg_pathway_id", ""),
            "prominence": m.get("prominence", "peripheral"),
            "evidence": m.get("evidence", ""),    # carried onto the edge (Stage 3)
            "pathway_id": ko,           # raw koNNNNN key; Stage 3 -> kegg.pathway:ko*
            "pathway_name": name,
            "resolution_method": method,          # audit only (not on the edge)
        })
        pathway_methods[method] += 1

    # Truncation quantification over the RAW captured identifiers.
    all_ids = [str(i).strip() for g in raw_genes for i in (g.get("identifiers") or []) if str(i).strip()]
    truncated_ids = [i for i in all_ids if is_truncated_id(i)]
    unresolved_reasons = Counter(
        _unresolved_reason(r) for r in gene_records if not r["locus_tag"]
    )

    meta = dict(data.get("metadata", {}))
    n_genes_resolved = sum(1 for r in gene_records if r["locus_tag"])
    n_paths_resolved = sum(1 for r in pathway_records if r["pathway_id"])
    meta["resolution"] = {
        "gene_records": len(gene_records),
        "genes_resolved": n_genes_resolved,
        "gene_methods": dict(gene_methods),
        "pathway_records": len(pathway_records),
        "pathways_resolved": n_paths_resolved,
        "pathway_methods": dict(pathway_methods),
        "paper_strains": paper_strains,
        "captured_identifiers": len(all_ids),
        "truncated_identifiers": len(truncated_ids),
        "truncated_examples": sorted(set(truncated_ids))[:10],
        "unresolved_reasons": dict(unresolved_reasons),
    }

    out = {"metadata": meta, "genes": gene_records, "pathways": pathway_records}
    resolved_path(paper_dir).parent.mkdir(parents=True, exist_ok=True)
    resolved_path(paper_dir).write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    _write_resolution_report(paper_dir, meta, gene_records, pathway_records)

    return {
        "paper": meta.get("paper", paper_dir.name),
        "out_path": resolved_path(paper_dir),
        "gene_records": len(gene_records),
        "genes_resolved": n_genes_resolved,
        "pathway_records": len(pathway_records),
        "pathways_resolved": n_paths_resolved,
    }


def _write_resolution_report(
    paper_dir: Path, meta: dict, gene_records: list[dict], pathway_records: list[dict],
) -> None:
    """Write a per-paper diagnostic report (mirrors the Stage-4 resolve report)."""
    r = meta["resolution"]
    g_tot, g_res = r["gene_records"], r["genes_resolved"]
    p_tot, p_res = r["pathway_records"], r["pathways_resolved"]
    g_pct = 100 * g_res / max(g_tot, 1)
    p_pct = 100 * p_res / max(p_tot, 1)

    lines = [
        f"Topics resolution report: {meta.get('paper', paper_dir.name)}",
        f"Strains: {r['paper_strains']}",
        "",
        f"Genes:    {g_res}/{g_tot} resolved ({g_pct:.1f}%)",
        f"  methods: " + ", ".join(f"{k}={v}" for k, v in sorted(r["gene_methods"].items())),
        f"Pathways: {p_res}/{p_tot} resolved ({p_pct:.1f}%)",
        f"  methods: " + ", ".join(f"{k}={v}" for k, v in sorted(r["pathway_methods"].items())),
        "",
        f"Captured identifiers: {r['captured_identifiers']}  "
        f"(truncated/bare-number: {r['truncated_identifiers']}"
        + (f", e.g. {r['truncated_examples']}" if r["truncated_examples"] else "")
        + ")",
        f"Unresolved gene reasons: {r['unresolved_reasons']}",
        "    descriptive_only = function/enzyme phrase, no symbol or id (paper style)",
        "    truncated_id_only = bare locus number, recoverable via strain prefix",
        "    lookup_miss       = had a name/full id but no mapping match",
        "",
    ]

    unresolved = [g for g in gene_records if not g["locus_tag"]]
    if unresolved:
        lines.append(f"Unresolved gene mentions ({len(unresolved)}):")
        for g in unresolved:
            label = g.get("gene_name") or g.get("surface_form") or "?"
            ids = g.get("identifiers") or []
            lines.append(
                f"  - {label}  [{_unresolved_reason(g)}]  strain={g.get('strain')}"
                + (f"  ids={ids}" if ids else "")
            )
        lines.append("")

    unresolved_pw = [p for p in pathway_records if not p["pathway_id"]]
    if unresolved_pw:
        lines.append(f"Unresolved pathway mentions ({len(unresolved_pw)}):")
        for p in unresolved_pw:
            lines.append(
                f"  - {p.get('surface_form')}  (raw id={p.get('kegg_pathway_id_raw') or '—'})"
            )
        lines.append("")

    report_path(paper_dir).write_text("\n".join(lines), encoding="utf-8")


def load_kegg_data(path: Path = DEFAULT_KEGG_DATA) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"kegg_data.json not found at {path}. Run prepare_data step 6 first."
        )
    return json.loads(path.read_text())


# ─── Main ───────────────────────────────────────────────────────────────────────


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true",
                        help="Re-resolve even if publication_topics_resolved.json is newer")
    parser.add_argument("--papers", nargs="+", metavar="NAME",
                        help="Only papers whose name contains one of these (case-insensitive)")
    args = parser.parse_args()

    kegg_data = load_kegg_data()
    id_lookup, name_lookup = build_pathway_lookup(kegg_data)
    logger.info("KEGG pathway lookup: %d ids, %d names", len(id_lookup), len(name_lookup))

    mapping_cache: dict[str, MappingData | None] = {}
    results: list[dict] = []

    for paper_dir, pub, _pdf, strains in find_all_papers():
        paper = pub.get("papername", paper_dir.name)
        if args.papers and not any(p.lower() in paper.lower() for p in args.papers):
            continue
        if not has_topics(paper_dir):
            continue

        migrate_legacy(paper_dir)
        rp, tp = resolved_path(paper_dir), topics_path(paper_dir)
        if (not args.force and rp.exists() and tp.exists()
                and rp.stat().st_mtime >= tp.stat().st_mtime):
            logger.info("%s: up to date, skipping", paper)
            continue

        r = resolve_paper_topics(paper_dir, strains, id_lookup, name_lookup, mapping_cache)
        if r:
            results.append(r)
            logger.info(
                "%s: genes %d/%d resolved, pathways %d/%d resolved",
                paper, r["genes_resolved"], r["gene_records"],
                r["pathways_resolved"], r["pathway_records"],
            )

    if not results:
        print("No topics resolved (all up to date or nothing to do).")
        return

    g_res = sum(r["genes_resolved"] for r in results)
    g_tot = sum(r["gene_records"] for r in results)
    p_res = sum(r["pathways_resolved"] for r in results)
    p_tot = sum(r["pathway_records"] for r in results)
    print(f"\nResolved {len(results)} papers: "
          f"genes {g_res}/{g_tot}, pathways {p_res}/{p_tot}")


if __name__ == "__main__":
    main()
