"""MEROPS peptidase ontology adapter — observed-only (CAZy pattern).

Yields:
- MeropsFamily nodes (clan + family + subfamily, only codes reachable from
  calls observed in any strain's `<strain>.merops.calls.json`). Names come
  from the committed `cache/data/merops/merops_reference.json` (prepare_data
  step 9): family.txt type-example names with scan-library holotype fallback,
  clan.txt descriptions.
- Merops_family_is_a_merops_family parent edges (subfamily → family,
  family → clan; families MEROPS leaves clan-unassigned are roots).
- Gene_has_merops_family scored edges — one per calls.json candidate, attached
  at the called depth (tier-1 identifier calls attach at subfamily/family;
  the identifier survives as edge `best_hit_id`). Edge `call_class` is the
  read-first verdict: peptidase | inhibitor | nonpeptidase_homolog.

Two-class shape mirrors cazy_adapter; the calls.json-for-edge-evidence read
mirrors tcdb_adapter / interpro_adapter. Single evidence source
(merops-diamond). Design:
docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.controlled_vocab import VOCAB
from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.merops_diamond import (
    call_class,
    catalytic_type_word,
    classify_code,
    edge_target_code,
    family_type,
)

logger = logging.getLogger(__name__)

DEFAULT_REFERENCE_PATH = Path("cache/data/merops/merops_reference.json")


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _merops_node_id(code: str) -> str:
    """merops.clan:SC for clans, merops.family:S14 / merops.family:S08A else.

    Both prefixes are registered in bioregistry (verified), so the colon CURIE
    is the canonical form; the underscore fallback should never fire.
    """
    classified = classify_code(code)
    prefix = "merops.clan" if classified and classified[0] == 0 else "merops.family"
    return normalize_curie(f"{prefix}:{code}") or f"merops_{code}"


def _pfam_node_id(acc: str) -> str:
    return normalize_curie(f"pfam:{acc}") or f"pfam_{acc}"


# Edge-evidence fields copied verbatim from each candidate (call_class is derived).
_EDGE_EVIDENCE_FIELDS = (
    "tier", "confidence_score", "identity", "qcov", "evalue", "consensus_n",
)


class MeropsAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_merops_family edges.

    Reads the merged JSON (locus_tag → {protein_id, merops_ids}) for the gene
    set + protein_id join, and the strain's Phase-1 calls.json (WP_ → calls)
    for the per-candidate evidence folded onto each edge.
    """

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)
        self._calls = self._load_calls()

    def _load_calls(self) -> dict[str, dict]:
        strain = self.genome_dir.name
        path = self.genome_dir / "merops" / f"{strain}.merops.calls.json"
        if not path.exists():
            return {}
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return {str(k).strip(): v for k, v in data.items() if isinstance(v, dict)}

    def get_all_observed_codes(self) -> set[str]:
        """Distinct edge-target codes (family/subfamily) observed on this strain."""
        codes: set[str] = set()
        for rec in self._calls.values():
            for cand in rec.get("calls", []):
                target = edge_target_code(cand)
                if target:
                    codes.add(target)
        return codes

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            protein_id = (gene.get("protein_id") or "").strip()
            if not protein_id or not gene.get("merops_ids"):
                continue
            rec = self._calls.get(protein_id)
            if not rec:
                continue
            for cand in rec.get("calls", []):
                target = edge_target_code(cand)
                if target is None:
                    logger.debug(
                        f"MEROPS malformed candidate {cand.get('code')!r} on {locus_tag} "
                        f"(strain {self.genome_dir.name}), skipping"
                    )
                    continue
                props: dict = {
                    "call_class": VOCAB.check(
                        "Gene_has_merops_family", "call_class", call_class(cand)),
                    "best_hit_id": _clean_str(cand.get("best_hit_id")),
                    "best_hit_kind": VOCAB.check(
                        "Gene_has_merops_family", "best_hit_kind",
                        cand.get("best_hit_kind")),
                }
                for field in _EDGE_EVIDENCE_FIELDS:
                    if cand.get(field) is not None:
                        props[field] = cand[field]
                yield (
                    f"{locus_tag}-has_merops-{target}",
                    _gene_node_id(locus_tag),
                    _merops_node_id(target),
                    "gene_has_merops_family",
                    props,
                )
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(
            f"MeropsAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→MEROPS edges"
        )


class MultiMeropsAnnotationAdapter:
    """Multi-strain orchestrator: owns MeropsFamily nodes + parent edges."""

    def __init__(
        self,
        genome_config_file: str,
        reference_path: str | Path = DEFAULT_REFERENCE_PATH,
        pfam_node_ids: set[str] | None = None,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self.reference_path = Path(reference_path)
        self.pfam_node_ids = pfam_node_ids
        self._reference: dict | None = None
        self._strain_adapters: list[MeropsAnnotationAdapter] = []
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
                MeropsAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiMeropsAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """Load the committed reference — fails loudly when missing."""
        if not self.reference_path.exists():
            raise FileNotFoundError(
                f"MEROPS reference cache missing: {self.reference_path} — "
                f"run `bash scripts/prepare_data.sh --steps 10` "
                f"(uv run python -m multiomics_kg.download.build_merops_reference)."
            )
        with open(self.reference_path, encoding="utf-8") as fh:
            self._reference = json.load(fh)

    def _ref(self) -> dict:
        if self._reference is None:
            self.download_data()
        return self._reference  # type: ignore[return-value]

    def _all_observed_codes(self) -> set[str]:
        codes: set[str] = set()
        for adapter in self._strain_adapters:
            codes |= adapter.get_all_observed_codes()
        return codes

    def _parent_of(self, code: str) -> str | None:
        """Parent code: subfamily → its family; family → its clan (reference);
        clan → None. Clan lookups tolerate families absent from the reference."""
        classified = classify_code(code)
        if classified is None:
            return None
        level, _ = classified
        if level == 2:
            return code[:-1]  # S08A → S08
        if level == 1:
            return (self._ref()["families"].get(code) or {}).get("clan")
        return None

    def _expand_with_ancestors(self, observed: set[str]) -> set[str]:
        out: set[str] = set()
        for code in observed:
            cur: str | None = code
            while cur is not None and cur not in out:
                out.add(cur)
                cur = self._parent_of(cur)
        return out

    def _node_props(self, code: str, level: int, level_kind: str) -> dict:
        ref = self._ref()
        props = {
            "merops_id": code,
            "level": level,
            "level_kind": level_kind,
            "family_type": family_type(code),
        }
        # sparse (tcdb-superfamily pattern): inhibitors have no catalytic type —
        # omit the key so the Neo4j property is absent, not an empty string
        ct = catalytic_type_word(code)
        if ct is not None:
            props["catalytic_type"] = ct
        if level == 0:
            clan_meta = ref["clans"].get(code) or {}
            props["name"] = _clean_str(code)
            if clan_meta.get("description"):
                props["description"] = _clean_str(clan_meta["description"])
            if clan_meta.get("family_type"):
                props["family_type"] = clan_meta["family_type"]
        elif level == 1:
            fam_meta = ref["families"].get(code) or {}
            props["name"] = _clean_str(fam_meta.get("name") or code)
            if fam_meta.get("family_type"):
                props["family_type"] = fam_meta["family_type"]
            for key, val in (ref.get("cleavage", {}).get(code) or {}).items():
                if key == "cleavage_summary":
                    props[key] = _clean_str(val)
                elif key == "cleavage_p1_residues":
                    props[key] = [_clean_str(v) for v in val]
                else:
                    props[key] = val
        else:
            props["name"] = _clean_str(ref.get("subfamily_names", {}).get(code) or code)
        return props

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        all_codes = self._expand_with_ancestors(self._all_observed_codes())
        count = 0
        for code in sorted(all_codes):
            classified = classify_code(code)
            if classified is None:
                continue
            level, level_kind = classified
            yield _merops_node_id(code), "merops family", self._node_props(code, level, level_kind)
            count += 1
        logger.info(f"MultiMeropsAnnotationAdapter.get_nodes: {count} MeropsFamily nodes")

    def get_edges(self):
        all_codes = self._expand_with_ancestors(self._all_observed_codes())

        # 1. Parent edges: subfamily → family, family → clan. Every parent is a
        # node by construction (_expand_with_ancestors walks to the root), so
        # no filtering is needed — cazy precedent.
        parent_count = 0
        for code in sorted(all_codes):
            parent = self._parent_of(code)
            if parent is None:
                continue
            yield (
                f"{code}-parent-{parent}",
                _merops_node_id(code),
                _merops_node_id(parent),
                "merops_family_is_a_merops_family",
                {},
            )
            parent_count += 1

        # 2. Family→Pfam bridge (MEROPS interpro.txt, family-level only).
        # Dangling-proof: emitted only for injected, existing Pfam node ids
        # (TCDB-bridge precedent) — pfam_node_ids=None -> no bridge edges.
        bridge_count = 0
        if self.pfam_node_ids is not None:
            for fam, pfams in sorted(self._ref().get("pfam_bridge", {}).items()):
                if fam not in all_codes:
                    continue
                for pf, n in sorted(pfams.items()):
                    pf_id = _pfam_node_id(pf)
                    if pf_id not in self.pfam_node_ids:
                        continue
                    yield (
                        f"{fam}-has_pfam-{pf}",
                        _merops_node_id(fam),
                        pf_id,
                        "merops_family_has_pfam_domain",
                        {"member_id_count": n},
                    )
                    bridge_count += 1

        # 3. Gene→MEROPS edges via per-strain delegation. Every edge target is
        # in the node set by construction (same edge_target_code both places).
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(
            f"MultiMeropsAnnotationAdapter.get_edges: {parent_count} parent, "
            f"{bridge_count} bridge, {gene_count} gene edges"
        )
