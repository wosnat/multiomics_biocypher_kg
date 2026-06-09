"""Adapter for publication "discusses" edges (Stage 3 of the discuss-edges pipeline).

Reads each paper's resolved topic file (`publication_topics/topics_resolved.json`,
written by `resolve_paper_topics`) and emits literature-index edges:

  - ``Publication_discusses_gene``          (Publication → Gene)
  - ``Publication_discusses_kegg_pathway``  (Publication → KeggTerm)

Pure edge adapter — the Publication / Gene / KeggTerm nodes already exist.
Edge source = the paper's Publication node id (`doi:{doi}`, with the DOI taken from
the paperconfig or the PDF-extraction cache, mirroring the omics + cluster adapters).
Targets: gene → ``ncbigene:{locus_tag}``; pathway → ``kegg.pathway:ko*``.
Each edge carries ``prominence`` (central|peripheral) + the extraction ``evidence``
quote. Only resolved mentions are emitted; resolution info stays in the json.
"""
import json
import logging
from pathlib import Path

from multiomics_kg.adapters.cluster_adapter import _load_pdf_cache
from multiomics_kg.extraction.topics.extraction_utils import migrate_legacy, resolved_path
from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.paperconfig_utils import (
    get_paper_name,
    get_publication,
    load_all_paperconfigs,
    load_paperconfig,
)

logger = logging.getLogger(__name__)


def _clean_str(value: str) -> str:
    """Strip BioCypher-breaking characters from a string property."""
    return str(value).replace("'", "^").replace("|", "")


class PublicationTopicsAdapter:
    """Emits discuss-edges for one paperconfig's resolved topics file."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.paper_dir = Path(config_file).parent
        self.doi = self._extract_doi()

    def _extract_doi(self) -> str:
        """DOI from paperconfig, falling back to the PDF-extraction cache.

        Same precedence the omics/cluster adapters use so the source id matches
        the existing Publication node.
        """
        pub = get_publication(self.config)
        doi = pub.get("doi", "")
        if doi:
            return doi
        pdf_path = pub.get("papermainpdf", "")
        if pdf_path:
            cached = _load_pdf_cache().get(pdf_path, {})
            doi = cached.get("publication", {}).get("doi", "")
        return doi or ""

    def _load_resolved(self) -> dict:
        migrate_legacy(self.paper_dir)
        p = resolved_path(self.paper_dir)
        if not p.exists():
            return {}
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            logger.warning("Failed to parse resolved topics: %s", p)
            return {}

    def get_edges(self):
        """Yield (edge_id, source_id, target_id, label, properties) tuples."""
        if not self.doi:
            logger.warning("No DOI for %s — skipping discuss-edges", self.paper_name)
            return
        data = self._load_resolved()
        if not data:
            return
        pub_id = normalize_curie(f"doi:{self.doi}") or f"doi:{self.doi}"

        # Gene edges — dedup by target, prominence central-wins.
        seen_genes: dict[str, dict] = {}
        for g in data.get("genes", []):
            lt = g.get("locus_tag")
            if not lt:
                continue
            target = normalize_curie(f"ncbigene:{lt}") or f"ncbigene:{lt}"
            prom = g.get("prominence", "peripheral")
            if target in seen_genes:
                if prom == "central":
                    seen_genes[target]["prominence"] = "central"
                continue
            seen_genes[target] = {
                "prominence": _clean_str(prom),
                "evidence": _clean_str(g.get("evidence", "")),
            }
        for target, props in seen_genes.items():
            edge_id = f"{self.doi}::discusses_gene::{target}"
            yield (edge_id, pub_id, target, "publication_discusses_gene", props)

        # Pathway edges — dedup by target, prominence central-wins.
        seen_paths: dict[str, dict] = {}
        for p in data.get("pathways", []):
            ko = p.get("pathway_id")
            if not ko:
                continue
            target = normalize_curie(f"kegg.pathway:{ko}") or f"kegg.pathway:{ko}"
            prom = p.get("prominence", "peripheral")
            if target in seen_paths:
                if prom == "central":
                    seen_paths[target]["prominence"] = "central"
                continue
            seen_paths[target] = {
                "prominence": _clean_str(prom),
                "evidence": _clean_str(p.get("evidence", "")),
            }
        for target, props in seen_paths.items():
            edge_id = f"{self.doi}::discusses_pathway::{target}"
            yield (edge_id, pub_id, target, "publication_discusses_kegg_pathway", props)


class MultiPublicationTopicsAdapter:
    """Reads paperconfig list(s) and delegates to PublicationTopicsAdapter instances."""

    def __init__(self, config_list_file, test_mode: bool = False, **kwargs):
        self.adapters: list[PublicationTopicsAdapter] = []
        list_files = config_list_file if isinstance(config_list_file, list) else [config_list_file]
        for pc_path, _config in load_all_paperconfigs([Path(lf) for lf in list_files]):
            self.adapters.append(
                PublicationTopicsAdapter(config_file=str(pc_path), test_mode=test_mode)
            )

    def download_data(self, cache: bool = True) -> None:
        """No external data — resolved files are produced by prepare_data."""
        return None

    def get_nodes(self):
        """No nodes — Publication / Gene / KeggTerm already exist."""
        return []

    def get_edges(self):
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges())
        return edges
