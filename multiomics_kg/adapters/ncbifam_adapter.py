"""NCBIfam ontology adapter — InterProScan Phase-2 KG integration (Task 13 of
the InterPro multi-ontology redesign).

Yields:
- NcbifamFamily nodes (observed-only, flat — no hierarchy: the source TSV
  carries no parent/child column, so an ``is_a`` closure is impossible;
  hierarchical context comes via the ``Ncbifam_family_in_interpro_entry``
  bridge instead), with name/family_type/gene_symbol/description from the
  committed reference cache (``cache/data/ncbifam/ncbifam_reference.json``,
  built by prepare_data step 9). Retired accessions absent from the reference
  still get a minimal node (name recovered from any strain's calls.json
  NCBIFAM facet row; ``family_type`` omitted) so their gene edges never
  dangle — same defensive fallback as ``interpro_adapter``.
- Gene_has_ncbifam_family edges (one per gene×accession) carrying evidence
  read DIRECTLY from the strain's faceted Phase-1 calls.json
  ``libraries.NCBIFAM`` rows: start/end envelope, evalue, AND score (both
  kept — single library, homogeneous HMMER bit scale; unlike InterPro's
  Gene_has_interpro_entry, which omits score because entries can be reached
  through multiple heterogeneous-scale member DBs). When an accession has
  more than one facet row on the same protein, the row with the lowest
  evalue (None sorts last) wins. Merged ``ncbifam_ids`` (a per-gene list of
  accessions) is the authority for which edges exist; when an authority
  accession is missing from the calls.json facet rows (skew between merge
  and calls artifacts), the edge is still emitted fail-soft with ``{}``.
- Ncbifam_family_in_interpro_entry bridge edges (NcbifamFamily → InterproEntry)
  connecting the new NCBIfam layer to the existing InterPro layer, derived
  from the calls.json NCBIFAM facet rows' ``ipr`` field. Dangling-proof:
  emitted only when the target InterPro entry is a member of the injected
  ``interpro_kept_ids`` set (Task 12's ``kept_node_accessions()`` — pfam bridge
  precedent); ``None`` (not injected) → emit no bridge edges at all.

Two-class shape mirrors ``interpro_adapter`` (per-strain edges + multi
orchestrator that owns nodes + the ontology-bridge edges). See
``docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md``.
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.adapters.interpro_adapter import _interpro_node_id

logger = logging.getLogger(__name__)


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _ncbifam_node_id(acc: str) -> str:
    # `ncbifam` is NOT a registered bioregistry prefix (normalize_prefix -> None
    # in this project's venv), so normalize_curie would fall back to the
    # underscore form anyway -- spelled out explicitly here (psortb/signalp
    # precedent) rather than relying on that fallback silently.
    return f"ncbifam_{acc}"


def _best_facet_row(rows: list[dict]) -> dict:
    """Pick the row with the lowest evalue (None sorts last)."""
    return min(rows, key=lambda r: (r.get("evalue") is None, r.get("evalue")))


class NcbifamAnnotationAdapter:
    """Per-strain adapter: Gene_has_ncbifam_family edges with precomputed evidence.

    Reads the merged JSON (locus_tag → {protein_id, ncbifam_ids}) for the gene
    set + protein_id join (authority for which edges exist), and the strain's
    faceted Phase-1 calls.json (WP_ → {libraries: {NCBIFAM: [...]}}) for the
    per-accession evidence — no aggregation happens at merge time for this
    field, so this adapter reads the raw facet rows directly (interpro_adapter
    precedent).
    """

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)
        self._calls = self._load_calls()

    def _load_calls(self) -> dict[str, dict]:
        strain = self.genome_dir.name
        path = self.genome_dir / "interproscan" / f"{strain}.interproscan.calls.json"
        if not path.exists():
            return {}
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return {str(k).strip(): v for k, v in data.items() if isinstance(v, dict)}

    def get_all_ncbifam_ids(self) -> set[str]:
        """Distinct NCBIfam accessions observed on this strain's genes."""
        ids: set[str] = set()
        for gene in self._genes.values():
            for acc in gene.get("ncbifam_ids") or []:
                if acc:
                    ids.add(acc)
        return ids

    def get_ncbifam_to_interpro(self) -> dict[str, str]:
        """acc → IPR mapping from this strain's NCBIFAM facet rows (for the bridge)."""
        out: dict[str, str] = {}
        for call in self._calls.values():
            for row in (call.get("libraries") or {}).get("NCBIFAM") or []:
                acc = row.get("accession")
                ipr = row.get("ipr")
                if acc and ipr:
                    out[acc] = ipr
        return out

    def get_facet_names(self) -> dict[str, str]:
        """acc → name observed in this strain's NCBIFAM facet rows — fallback
        naming for accessions retired out of the current reference TSV."""
        out: dict[str, str] = {}
        for call in self._calls.values():
            for row in (call.get("libraries") or {}).get("NCBIFAM") or []:
                acc = row.get("accession")
                name = row.get("name")
                if acc and name and acc not in out:
                    out[acc] = name
        return out

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            protein_id = (gene.get("protein_id") or "").strip()
            accs = gene.get("ncbifam_ids") or []      # merged = authority
            if not protein_id or not accs:
                continue
            call = self._calls.get(protein_id) or {}
            rows_by_acc: dict[str, list[dict]] = {}
            for row in (call.get("libraries") or {}).get("NCBIFAM") or []:
                acc = row.get("accession")
                if acc:
                    rows_by_acc.setdefault(acc, []).append(row)
            for acc in accs:
                rows = rows_by_acc.get(acc) or []
                props: dict = {}
                if rows:
                    best = _best_facet_row(rows)
                    for k in ("start", "end", "evalue", "score"):
                        if best.get(k) is not None:
                            props[k] = best[k]
                yield (f"{locus_tag}-has_ncbifam-{acc}", _gene_node_id(locus_tag),
                       _ncbifam_node_id(acc), "gene_has_ncbifam_family", props)
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(
            f"NcbifamAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→NCBIfam edges"
        )


class MultiNcbifamAdapter:
    """Multi-strain orchestrator: owns NcbifamFamily nodes + the InterPro bridge.

    Node set is observed-only (union of every strain's ``get_all_ncbifam_ids()``)
    — flat, no ancestor walk (unlike InterPro/TCDB/BRITE pruning: the source TSV
    has no hierarchy column). ``interpro_kept_ids`` is the InterPro adapter's
    kept node-id set (``MultiInterproAnnotationAdapter.kept_node_accessions()``)
    — injected by ``create_knowledge_graph`` so the bridge never references a
    non-existent InterproEntry node; ``None`` → emit no bridge edges.
    """

    def __init__(
        self,
        genome_config_file: str,
        cache_root: str | Path = "cache/data",
        interpro_kept_ids: set[str] | None = None,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self.cache_root = Path(cache_root)
        # None = InterPro kept-id set not provided → emit NO bridge edges (can't
        # guarantee the InterproEntry endpoints exist). A provided set (even
        # empty) prunes to it.
        self.interpro_kept_ids = interpro_kept_ids
        self._reference: dict[str, dict] = {}
        self._strain_adapters: list[NcbifamAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                NcbifamAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiNcbifamAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """Load the committed NCBIfam reference cache (built by prepare_data step 9)."""
        path = self.cache_root / "ncbifam" / "ncbifam_reference.json"
        if not path.exists():
            raise FileNotFoundError(
                f"NCBIfam reference cache missing: {path}. "
                f"Run `bash scripts/prepare_data.sh --steps 9 --force` first."
            )
        with open(path, encoding="utf-8") as fh:
            self._reference = json.load(fh)
        logger.info(f"MultiNcbifamAdapter: loaded {len(self._reference)} reference entries")

    def _observed_ids(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_ncbifam_ids()
        return ids

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._reference:
            self.download_data()
        observed = self._observed_ids()
        facet_names: dict[str, str] | None = None
        for acc in sorted(observed):
            ref = self._reference.get(acc)
            if ref is None:
                # Retired accession: absent from the current reference TSV.
                # Fall back to any strain's calls.json facet row name so the
                # gene edge referencing it never dangles.
                if facet_names is None:
                    facet_names = {}
                    for adapter in self._strain_adapters:
                        for a, name in adapter.get_facet_names().items():
                            facet_names.setdefault(a, name)
                    facet_names = facet_names or {}
                name = facet_names.get(acc, "")
                logger.warning(f"NCBIfam accession {acc} not in reference (retired?); emitting minimal node")
                props = {"name": _clean_str(name), "ncbifam_id": acc, "level": 0}
            else:
                props = {
                    "name": _clean_str(ref.get("name")),
                    "ncbifam_id": acc,
                    "family_type": _clean_str(ref.get("family_type")),
                    "level": 0,
                }
                gene_symbol = _clean_str(ref.get("gene_symbol"))
                if gene_symbol:
                    props["gene_symbol"] = gene_symbol
                description = _clean_str(ref.get("description"))
                if description:
                    props["description"] = description
            yield _ncbifam_node_id(acc), "ncbifam family", props
        logger.info(
            f"MultiNcbifamAdapter.get_nodes: {len(observed)} NcbifamFamily nodes (observed-only, flat)"
        )

    def get_edges(self):
        if not self._reference:
            self.download_data()

        # 1. NcbifamFamily → InterproEntry bridge edges (dangling-proof)
        acc_to_ipr: dict[str, str] = {}
        for adapter in self._strain_adapters:
            acc_to_ipr.update(adapter.get_ncbifam_to_interpro())
        bridge = 0
        for acc, ipr in sorted(acc_to_ipr.items()):
            # dangling-proof: no set provided → no bridges; set provided → require membership
            if self.interpro_kept_ids is None or ipr not in self.interpro_kept_ids:
                continue
            yield (
                f"{acc}-in_interpro-{ipr}",
                _ncbifam_node_id(acc),
                _interpro_node_id(ipr),
                "ncbifam_family_in_interpro_entry",
                {},
            )
            bridge += 1

        # 2. Gene → NcbifamFamily edges via per-strain delegation
        gene = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene += 1
        logger.info(
            f"MultiNcbifamAdapter.get_edges: {bridge} interpro-bridge, {gene} gene edges"
        )
