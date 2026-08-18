"""InterPro ontology adapter — InterProScan Phase-2 KG integration.

Yields:
- InterproEntry nodes (observed IPR entries + their is-a ancestors), with
  name/type/level/description from the committed reference cache
  (``cache/data/interpro/interpro_reference.json``, built by prepare_data step 9).
  ``description`` is sparse (Task 6): present only when the reference entry
  carries a non-empty narrative.
- Interpro_entry_is_a_interpro_entry parent edges (child → parent), pruned to kept
  entries (mirrors TCDB/BRITE subhierarchy pruning — never dangles).
- Gene_has_interpro_entry edges (one per gene×entry) carrying evidence read
  DIRECTLY from the strain's faceted Phase-1 calls.json per-protein
  ``interpro_entries`` rollup (precomputed at normalize time — no aggregation
  here): start/end envelope, evalue (+ ``evalue_library`` naming which member DB
  it came from), member-DB libraries, match_count. NO ``score``. Scored ontology
  edge (the Changes_expression_of / psortb / signalp precedent). NO e-value
  filtering — each member DB pre-applies its own curated threshold (see the
  design spec §1). Merged ``interpro_entries`` (a per-gene list of accessions) is
  the authority for which edges exist; when an authority entry is missing from
  the calls.json rollup (skew between merge and calls artifacts), the edge is
  still emitted fail-soft with ``{match_count: 0, libraries: []}``.
- Pfam_in_interpro_entry bridge edges (Pfam → InterproEntry) connecting the
  existing eggNOG Pfam layer to the InterPro layer, derived from the calls.json
  PFAM signature→entry mapping. Dangling-proof: emitted only when the source Pfam
  node exists (injected ``pfam_node_ids``, BRITE-``known_ko_ids`` style) and the
  target entry is kept.
- **Layer A** (design 2026-08-10): Interpro_entry_related_to_ec_number /
  _related_to_cazy_family cross-references from the reference's entry-level EC/CAZy
  xrefs. A deliberately WEAK, recall-biased ROUTER (weak ``related_to`` verb) that
  homes the multi-EC / DOMAIN ECs and fold CAZy that Layer B refuses to stamp on
  genes — read backward it is low-precision; never use it to assign gene function.
  These edges carry NO properties (deleted 2026-08-16, spec §3 R3): the old
  ``ambiguous`` flag is APPROXIMATELY derivable from the graph — a consumer writes
  ``WITH n, count(r) AS k WHERE k > 1 OR n.interpro_type <> 'FAMILY'`` — and
  ``source_db`` was a hardcoded literal (``"interpro.xml"``), which the edge type
  already says. The derivation is not exact, because edges are pruned to
  gene-referenced ECs *after* the flag would have been computed: of 5,021 entries
  carrying a Layer-A EC edge it differs from the old flag on 105 — **42** read as
  specific though the entry lists competing alternatives no gene in this corpus
  carries (IPR001461 "Aspartic peptidase A1" lists 20 distinct ECs, 1 survives),
  **32** drop a class+member nuance (``2.6.1.-`` plus ``2.6.1.50``), and **31**
  are cases the derivation gets *right* where the old flag was a false positive
  (raw tokens collapsing to one EC after obsolete-number remapping, e.g.
  ``3.3.1.1`` → ``3.13.2.1``). Treat it as a floor on multiplicity, not a
  faithful replay. Dangling-proof by pruning to the EC/CAZy nodes the gene edges
  already created (self-computed from the merged JSONs — no injection).

Two-class shape mirrors cazy_adapter (per-strain edges + multi orchestrator that
owns nodes/hierarchy). See
``docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md``.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.adapters.cazy_adapter import _cazy_node_id, _most_specific_id
from multiomics_kg.download.utils.annotation_transforms import _TRANSFORMS

logger = logging.getLogger(__name__)

_BARE_EC_RE = re.compile(r"^\d+\.\d+\.\d+$")  # 3.4.21 (no 4th field) — needs `.-`


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _interpro_node_id(interpro_id: str) -> str:
    return normalize_curie(f"interpro:{interpro_id}") or f"interpro_{interpro_id}"


def _pfam_node_id(pfam_accession: str) -> str:
    return normalize_curie(f"pfam:{pfam_accession}") or f"pfam_{pfam_accession}"


def _ec_node_id(ec_number: str) -> str:
    """EC node id — must match functional_annotation_adapter / ec_adapter."""
    return normalize_curie(f"eccode:{ec_number}") or f"eccode_{ec_number}"


def _normalize_ec_token(raw: str) -> list[str]:
    """Normalise a raw InterPro EC token the SAME way Layer B does, so a Layer-A
    edge can be pruned against the EC nodes the gene edges created (bare 3-level →
    ``.-``; ``normalize_ec`` remaps obsolete numbers). May return 0..n tokens."""
    raw = (raw or "").strip()
    if _BARE_EC_RE.match(raw):
        raw = raw + ".-"
    fn = _TRANSFORMS.get("normalize_ec")
    out = fn(raw) if fn else raw
    toks = out if isinstance(out, list) else [out]
    return [t for t in (str(x).strip() for x in toks) if t]


def kept_ids(observed: set[str], reference: dict[str, dict]) -> set[str]:
    """Observed entries + every is-a ancestor (walk the reference parent chain).

    Cycle-safe. Ensures hierarchy edges never dangle (TCDB/BRITE pruning pattern). Pure.
    """
    kept: set[str] = set()
    for acc in observed:
        cur: str | None = acc
        seen: set[str] = set()
        while cur and cur not in seen:
            seen.add(cur)
            kept.add(cur)
            ref = reference.get(cur)
            cur = ref.get("parent") if ref else None
    return kept


class InterproAnnotationAdapter:
    """Per-strain adapter: Gene_has_interpro_entry edges with precomputed evidence.

    Reads the merged JSON (locus_tag → {protein_id, interpro_entries}) for the gene
    set + protein_id join (authority for which edges exist), and the strain's
    faceted Phase-1 calls.json (WP_ → {interpro_entries: {IPR: {...}}, ...}) for
    the per-entry evidence rollup already computed at normalize time — no
    aggregation happens here.
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

    def get_all_interpro_ids(self) -> set[str]:
        """Distinct IPR entries observed on this strain's genes."""
        ids: set[str] = set()
        for gene in self._genes.values():
            for acc in gene.get("interpro_entries") or []:
                if acc:
                    ids.add(acc)
        return ids

    def observed_ec_node_ids(self) -> set[str]:
        """EC node ids that this strain's gene→EC edges reference (dangling-proof
        target set for the Layer-A InterproEntry→EcNumber router)."""
        out: set[str] = set()
        for gene in self._genes.values():
            for ec in gene.get("ec_numbers") or []:
                if ec:
                    out.add(_ec_node_id(ec))
        return out

    def observed_cazy_node_ids(self) -> set[str]:
        """CAZy node ids that this strain's gene→CAZy edges reference."""
        out: set[str] = set()
        for gene in self._genes.values():
            for tok in gene.get("cazy_ids") or []:
                spec = _most_specific_id(tok)
                if spec:
                    out.add(_cazy_node_id(spec))
        return out

    def get_pfam_to_interpro(self) -> dict[str, str]:
        """PF* → IPR* mapping from this strain's PFAM facet rows (for the bridge)."""
        out: dict[str, str] = {}
        for call in self._calls.values():
            for row in (call.get("libraries") or {}).get("PFAM") or []:
                sig = row.get("accession")
                ipr = row.get("ipr")
                if sig and ipr:
                    out[sig.split(".")[0]] = ipr
        return out

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            protein_id = (gene.get("protein_id") or "").strip()
            entries = gene.get("interpro_entries") or []      # merged = authority
            if not protein_id or not entries:
                continue
            rollups = (self._calls.get(protein_id) or {}).get("interpro_entries") or {}
            for acc in entries:
                ent = rollups.get(acc)
                props: dict = {"match_count": 0, "libraries": []}
                if ent:
                    props = {"match_count": ent.get("match_count") or 0,
                             "libraries": [_clean_str(x) for x in ent.get("libraries") or []]}
                    for k in ("start", "end"):
                        if ent.get(k) is not None:
                            props[k] = ent[k]
                    if ent.get("evalue") is not None:
                        props["evalue"] = ent["evalue"]
                        props["evalue_library"] = _clean_str(ent.get("evalue_library"))
                yield (f"{locus_tag}-has_interpro-{acc}", _gene_node_id(locus_tag),
                       _interpro_node_id(acc), "gene_has_interpro_entry", props)
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(
            f"InterproAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→InterPro edges"
        )


class MultiInterproAnnotationAdapter:
    """Multi-strain orchestrator: owns InterproEntry nodes + hierarchy + Pfam bridge.

    Prunes the InterPro node set to entries observed on any strain PLUS their is-a
    ancestors (so hierarchy edges never dangle). ``pfam_node_ids`` is the KG's
    global set of raw Pfam accessions (PF*) with a node — injected by
    ``create_knowledge_graph`` (BRITE-``known_ko_ids`` precedent) so the Pfam bridge
    never references a non-existent Pfam node. ``ec_node_ids`` is the same
    guarantee for the Layer-A EC router: the EC node ids
    ``MultiEcAnnotationAdapter`` emits, needed because a gene's ``ec_numbers`` can
    name an obsolete/invalid EC that Expasy has no node for.
    """

    def __init__(
        self,
        genome_config_file: str,
        cache_root: str | Path = "cache/data",
        pfam_node_ids: set[str] | None = None,
        ec_node_ids: set[str] | None = None,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self.cache_root = Path(cache_root)
        # None = Pfam node set not provided → emit NO bridge edges (can't guarantee
        # the Pfam endpoints exist). A provided set (even empty) prunes to it.
        self.pfam_node_ids = (
            None if pfam_node_ids is None
            else {p.split(".")[0] for p in pfam_node_ids}
        )
        # Same contract for Layer-A EC edges: None → emit none.
        self.ec_node_ids = ec_node_ids
        self._reference: dict[str, dict] = {}
        self._strain_adapters: list[InterproAnnotationAdapter] = []
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
                InterproAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiInterproAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """Load the committed InterPro reference cache (built by prepare_data step 9)."""
        path = self.cache_root / "interpro" / "interpro_reference.json"
        if not path.exists():
            raise FileNotFoundError(
                f"InterPro reference cache missing: {path}. "
                f"Run `bash scripts/prepare_data.sh --steps 9 --force` first."
            )
        with open(path, encoding="utf-8") as fh:
            self._reference = json.load(fh)
        logger.info(f"MultiInterproAnnotationAdapter: loaded {len(self._reference)} reference entries")

    def _observed_ids(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_interpro_ids()
        return ids

    def _kept_ids(self, observed: set[str]) -> set[str]:
        return kept_ids(observed, self._reference)

    def kept_node_accessions(self) -> set[str]:
        """Full pruned InterProEntry node-id set (observed + is-a ancestors).

        Public accessor for downstream consumers (Task 14) that need the exact
        set of InterPro accessions this adapter will emit nodes for, without
        duplicating the observed-ids + ancestor-walk logic.
        """
        if not self._reference:
            self.download_data()
        return self._kept_ids(self._observed_ids())

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._reference:
            self.download_data()
        observed = self._observed_ids()
        kept = self._kept_ids(observed)
        for acc in sorted(kept):
            ref = self._reference.get(acc)
            if ref is None:
                # Defensive: observed entry absent from reference (retired ID?).
                # Emit a minimal node so its gene edge never dangles.
                logger.warning(f"InterPro entry {acc} not in reference; emitting minimal node")
                props = {"name": "", "interpro_id": acc, "interpro_type": "", "level": 0}
            else:
                props = {
                    "name": _clean_str(ref.get("name")),
                    "interpro_id": acc,
                    "interpro_type": _clean_str(ref.get("type")),
                    "level": int(ref.get("level") or 0),
                }
                description = _clean_str(ref.get("description"))
                if description:
                    props["description"] = description
            yield _interpro_node_id(acc), "interpro entry", props
        logger.info(
            f"MultiInterproAnnotationAdapter.get_nodes: {len(kept)} InterproEntry nodes "
            f"({len(observed)} observed + {len(kept) - len(observed)} ancestors)"
        )

    def get_edges(self):
        if not self._reference:
            self.download_data()
        observed = self._observed_ids()
        kept = self._kept_ids(observed)

        # 1. is-a hierarchy edges (child → parent), both endpoints kept
        hier = 0
        for acc in sorted(kept):
            ref = self._reference.get(acc)
            parent = ref.get("parent") if ref else None
            if parent and parent in kept:
                yield (
                    f"{acc}-is_a-{parent}",
                    _interpro_node_id(acc),
                    _interpro_node_id(parent),
                    "interpro_entry_is_a_interpro_entry",
                    {},
                )
                hier += 1

        # 2. Pfam → InterproEntry bridge edges (dangling-proof)
        pf_to_ipr: dict[str, str] = {}
        for adapter in self._strain_adapters:
            pf_to_ipr.update(adapter.get_pfam_to_interpro())
        bridge = 0
        for pf, ipr in sorted(pf_to_ipr.items()):
            if ipr not in kept:
                continue
            # dangling-proof: no set provided → no bridges; set provided → require membership
            if self.pfam_node_ids is None or pf not in self.pfam_node_ids:
                continue
            yield (
                f"{pf}-in_interpro-{ipr}",
                _pfam_node_id(pf),
                _interpro_node_id(ipr),
                "pfam_in_interpro_entry",
                {},
            )
            bridge += 1

        # 3. Layer A — InterproEntry → EC / CAZy cross-references (design 2026-08-10).
        #    A deliberately WEAK, recall-biased ROUTER: read outward from a gene's
        #    known family it corroborates; read backward (carries EC → is that enzyme)
        #    it is low-precision. NEVER use it to assign a gene its function — that is
        #    what Gene_catalyzes_ec_number (Layer B) is for. These edges carry NO
        #    properties — a consumer APPROXIMATELY derives the old `ambiguous` flag
        #    as `count(r) > 1 OR n.interpro_type <> 'FAMILY'` (a floor, not a replay:
        #    the prune below runs after the flag would have been computed — see the
        #    module docstring for the 42/32/31 breakdown). Pruned to EC/CAZy nodes
        #    the gene edges already created. For EC that set is further intersected
        #    with the injected `ec_node_ids`, because a gene's ec_numbers can name an
        #    obsolete/invalid EC (InterPro xref) that Expasy has no node for —
        #    MultiEcAnnotationAdapter drops those gene edges, so Layer A must drop
        #    them too or it reintroduces the dangling target.
        observed_ec: set[str] = set()
        observed_cazy: set[str] = set()
        for adapter in self._strain_adapters:
            observed_ec |= adapter.observed_ec_node_ids()
            observed_cazy |= adapter.observed_cazy_node_ids()
        observed_ec = set() if self.ec_node_ids is None else observed_ec & self.ec_node_ids
        la_ec = la_cazy = 0
        for acc in sorted(kept):
            ref = self._reference.get(acc)
            if not ref:
                continue
            ecs = ref.get("ec_numbers") or []
            if ecs:
                for raw in ecs:
                    for norm in _normalize_ec_token(raw):
                        nid = _ec_node_id(norm)
                        if nid not in observed_ec:
                            continue
                        yield (
                            f"{acc}-related_ec-{norm}",
                            _interpro_node_id(acc),
                            nid,
                            "interpro_entry_related_to_ec_number",
                            {},
                        )
                        la_ec += 1
            czs = ref.get("cazy_ids") or []
            if czs:
                for cz in czs:
                    spec = _most_specific_id(cz)
                    if not spec:
                        continue
                    nid = _cazy_node_id(spec)
                    if nid not in observed_cazy:
                        continue
                    yield (
                        f"{acc}-related_cazy-{cz}",
                        _interpro_node_id(acc),
                        nid,
                        "interpro_entry_related_to_cazy_family",
                        {},
                    )
                    la_cazy += 1

        # 4. Gene → InterproEntry edges via per-strain delegation
        gene = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene += 1
        logger.info(
            f"MultiInterproAnnotationAdapter.get_edges: {hier} hierarchy, {bridge} pfam-bridge, "
            f"{la_ec} related-EC, {la_cazy} related-CAZy, {gene} gene edges"
        )
