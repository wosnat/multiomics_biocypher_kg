"""Gene ID Mapping Graph — iterative convergence.

Builds a two-tier lookup from locus_tag anchors and alternative IDs collected
from genome annotations (NCBI, Cyanorak, UniProt) and paper supplementary
tables (id_translation, annotation_gff, csv id_columns).

Tier classification
-------------------
Tier 1 (gene-unique): locus_tag, locus_tag_ncbi, locus_tag_cyanorak,
    old_locus_tag, alternative_locus_tag, jgi_id, probeset_id,
    uniprot_entry_name (after stripping _ORGANISM suffix),
    cds_fna_id (lcl|<accession>_cds_<protein_id>_<n> from cds_from_genomic.fna)
    → goes into specific_lookup (1:1)
    → conflicts are data errors (flagged in conflicts dict)

Tier 2 (protein-level): protein_id_refseq, protein_id, uniprot_accession
    → goes into multi_lookup (1:many; paralogs may share the same protein)
    → NOT a conflict when multiple genes share one Tier 2 ID

Tier 3 (generic): gene_name, gene_synonym, gene_oln, em_preferred_name, *
    → goes into multi_lookup (1:many by design)
    → never triggers convergence; used for resolution only as last resort

Algorithm
---------
1. Seed: for each locus_tag, add its annotation IDs directly to the graph.
2. Iterate: process all paper-source rows repeatedly until no new ID is added
   to specific_lookup or multi_lookup (typically 2–3 passes).
3. Build: call to_json_structure() to get the v2 gene_id_mapping.json payload.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

# ─── Tier classification ──────────────────────────────────────────────────────

TIER1_TYPES: frozenset[str] = frozenset({
    "locus_tag",
    "locus_tag_ncbi",
    "locus_tag_cyanorak",
    "old_locus_tag",
    "alternative_locus_tag",
    "jgi_id",
    "probeset_id",
    "uniprot_entry_name",
    "cds_fna_id",  # lcl|<accession>_cds_<protein_id>_<n> from cds_from_genomic.fna
})

TIER2_TYPES: frozenset[str] = frozenset({
    "protein_id_refseq",
    "protein_id",
    "uniprot_accession",
})

# Everything else → Tier 3 (generic, many-to-many by design)


def get_id_tier(id_type: str) -> int:
    """Return 1, 2, or 3 for the given id_type."""
    if id_type in TIER1_TYPES:
        return 1
    if id_type in TIER2_TYPES:
        return 2
    return 3


def normalize_id(id_val: str, id_type: str) -> list[str]:
    """Return candidate normalized forms of id_val (raw form first).

    For uniprot_entry_name: also return the form with _ORGANISM stripped.
    Example: "DNAA_PROM0" → ["DNAA_PROM0", "DNAA"]
    """
    id_val = str(id_val).strip()
    if not id_val or id_val.lower() in ("nan", ""):
        return []
    candidates: list[str] = [id_val]
    if id_type == "uniprot_entry_name":
        idx = id_val.rfind("_")
        if idx > 0:
            stripped = id_val[:idx]
            if stripped and stripped not in candidates:
                candidates.append(stripped)
    return candidates


# ─── Main graph class ─────────────────────────────────────────────────────────


class GeneIdGraph:
    """Iterative gene ID equivalence graph for one organism/strain.

    All operations are O(n) or better; no Union-Find needed at these dataset
    sizes (~2K genes, ~20K alt-IDs).
    """

    def __init__(self) -> None:
        # Primary outputs
        self.specific_lookup: dict[str, str] = {}       # Tier 1 id → locus_tag (1:1)
        self.multi_lookup: dict[str, list[str]] = {}    # Tier 2+3 id → [locus_tag, ...]
        self.conflicts: dict[str, list[str]] = {}       # Tier 1 id → [lt1, lt2] (data errors)

        # Per-gene ID collections (for gene_id_mapping.json "genes" section)
        self._genes: dict[str, dict[str, list]] = {}

        # Processing statistics
        self._stats: dict[str, Any] = {
            "passes": 0,
            "new_tier1_total": 0,
            "new_tier2_total": 0,
            "unresolved_rows_per_source": defaultdict(int),
        }

    # ── Seeding ───────────────────────────────────────────────────────────────

    def add_anchor(self, locus_tag: str) -> None:
        """Register a canonical locus_tag (called once per gene during seeding)."""
        self.specific_lookup[locus_tag] = locus_tag
        if locus_tag not in self._genes:
            self._genes[locus_tag] = {"tier1_ids": [], "tier2_ids": [], "tier3_ids": []}

    def add_id_for_gene(
        self, locus_tag: str, id_val: str, id_type: str, source: str
    ) -> None:
        """Add a known alt-ID for a known gene (seeding phase).

        Handles normalization (e.g., uniprot_entry_name suffix stripping).
        Silently skips values that are empty, NaN, or identical to the anchor.
        """
        for candidate in normalize_id(id_val, id_type):
            tier = get_id_tier(id_type)
            self._add_mapping(locus_tag, candidate, id_type, tier, source)

    # ── Paper-source processing ───────────────────────────────────────────────

    def process_row(
        self,
        row_ids: list[tuple[str, str]],
        source_name: str,
    ) -> bool:
        """Process one row from a paper source.

        Args:
            row_ids: list of (id_val, id_type) pairs from this row.
            source_name: label used in diagnostic output (paper/table name).

        Returns:
            True if any new entry was added to specific_lookup or multi_lookup.
        """
        anchor = self._find_anchor(row_ids)
        if anchor is None:
            self._stats["unresolved_rows_per_source"][source_name] += 1
            return False

        changed = False
        for id_val, id_type in row_ids:
            tier = get_id_tier(id_type)
            for candidate in normalize_id(id_val, id_type):
                if self._add_mapping(anchor, candidate, id_type, tier, source_name):
                    changed = True
            # Whitespace-split compound values (e.g. "dnaA PMM0001" in gene_name column)
            if " " in id_val and tier <= 2:
                for token in id_val.split():
                    token = token.strip()
                    if token and self._add_mapping(anchor, token, id_type, tier, source_name):
                        changed = True
        return changed

    def process_all_rows(
        self,
        all_rows: list[tuple[list[tuple[str, str]], str]],
    ) -> int:
        """Iterate over all rows until convergence.

        Args:
            all_rows: list of (row_ids, source_name) tuples from all sources.

        Returns:
            Number of passes taken to converge.
        """
        passes = 0
        while True:
            passes += 1
            # Reset unresolved counts each pass so we only report the final state
            self._stats["unresolved_rows_per_source"] = defaultdict(int)
            any_change = False
            for row_ids, source_name in all_rows:
                if self.process_row(row_ids, source_name):
                    any_change = True
            if not any_change:
                break
        self._stats["passes"] = passes
        return passes

    # ── Private helpers ───────────────────────────────────────────────────────

    def _find_anchor(self, row_ids: list[tuple[str, str]]) -> str | None:
        """Find a locus_tag anchor from this row's IDs.

        Priority:
        1. Tier 1 normalized match in specific_lookup
        2. Whitespace-split tokens from any field → specific_lookup
        3. Tier 2 singleton match in multi_lookup
        """
        # Phase 1: Tier 1 exact + normalized match
        for id_val, id_type in row_ids:
            if get_id_tier(id_type) == 1:
                for candidate in normalize_id(id_val, id_type):
                    if candidate in self.specific_lookup:
                        return self.specific_lookup[candidate]

        # Phase 2: whitespace-split tokens (handles "gene_name locus_tag" compound)
        for id_val, _ in row_ids:
            if " " in id_val:
                for token in id_val.split():
                    token = token.strip()
                    if token and token in self.specific_lookup:
                        return self.specific_lookup[token]

        # Phase 3: Tier 2 singleton
        for id_val, id_type in row_ids:
            if get_id_tier(id_type) == 2:
                matches = self.multi_lookup.get(id_val)
                if matches and len(matches) == 1:
                    return matches[0]

        return None

    def _add_mapping(
        self,
        anchor: str,
        id_val: str,
        id_type: str,
        tier: int,
        source: str,
    ) -> bool:
        """Add id_val → anchor mapping. Returns True if a new entry was created.

        Handles conflicts for Tier 1, expected multi for Tier 2, and stores
        per-gene ID records.  Tier 3 additions never return True (they don't
        trigger convergence iteration because they can't serve as anchors).
        """
        if not id_val:
            return False
        # Skip the anchor's self-mapping — it is added by add_anchor()
        if id_val == anchor and tier == 1 and id_type == "locus_tag":
            return False

        if tier == 1:
            return self._add_tier1(anchor, id_val, id_type, source)
        elif tier == 2:
            return self._add_tier2(anchor, id_val, id_type, source)
        else:
            self._add_tier3(anchor, id_val, id_type)
            return False  # Tier 3 never triggers convergence

    def _add_tier1(self, anchor: str, id_val: str, id_type: str, source: str) -> bool:
        existing = self.specific_lookup.get(id_val)
        if existing is not None:
            if existing == anchor:
                return False  # already there
            # Tier 1 conflict: same ID maps to two different locus_tags
            if id_val not in self.conflicts:
                self.conflicts[id_val] = [existing]
            if anchor not in self.conflicts[id_val]:
                self.conflicts[id_val].append(anchor)
            return False  # Don't overwrite; keep original mapping

        self.specific_lookup[id_val] = anchor
        gene = self._genes.setdefault(anchor, {"tier1_ids": [], "tier2_ids": [], "tier3_ids": []})
        existing_ids = {e["id"] for e in gene["tier1_ids"]}
        if id_val not in existing_ids:
            gene["tier1_ids"].append({"id": id_val, "type": id_type, "sources": [source]})
        self._stats["new_tier1_total"] += 1
        return True

    def _add_tier2(self, anchor: str, id_val: str, id_type: str, source: str) -> bool:
        if id_val not in self.multi_lookup:
            self.multi_lookup[id_val] = [anchor]
            gene = self._genes.setdefault(anchor, {"tier1_ids": [], "tier2_ids": [], "tier3_ids": []})
            existing_ids = {e["id"] for e in gene["tier2_ids"]}
            if id_val not in existing_ids:
                gene["tier2_ids"].append({"id": id_val, "type": id_type, "sources": [source]})
            self._stats["new_tier2_total"] += 1
            return True
        elif anchor not in self.multi_lookup[id_val]:
            self.multi_lookup[id_val].append(anchor)
            gene = self._genes.setdefault(anchor, {"tier1_ids": [], "tier2_ids": [], "tier3_ids": []})
            existing_ids = {e["id"] for e in gene["tier2_ids"]}
            if id_val not in existing_ids:
                gene["tier2_ids"].append({"id": id_val, "type": id_type, "sources": [source]})
            return True  # New anchor added → triggers convergence
        return False

    def _add_tier3(self, anchor: str, id_val: str, id_type: str) -> None:
        if id_val not in self.multi_lookup:
            self.multi_lookup[id_val] = [anchor]
        elif anchor not in self.multi_lookup[id_val]:
            self.multi_lookup[id_val].append(anchor)
        gene = self._genes.setdefault(anchor, {"tier1_ids": [], "tier2_ids": [], "tier3_ids": []})
        existing_ids = {e["id"] for e in gene["tier3_ids"]}
        if id_val not in existing_ids:
            gene["tier3_ids"].append({"id": id_val, "type": id_type})

    # ── Output ────────────────────────────────────────────────────────────────

    def to_json_structure(self, organism: str, strain: str) -> dict:
        """Build the gene_id_mapping.json v2 payload."""
        # Exclude locus_tag self-mappings from specific_lookup in output
        specific_out = {k: v for k, v in self.specific_lookup.items() if k != v}
        return {
            "version": 2,
            "organism": organism,
            "strain": strain,
            "stats": {
                "n_genes": len(self._genes),
                "n_specific": len(specific_out),
                "n_multi": len(self.multi_lookup),
                "n_conflicts": len(self.conflicts),
                "passes": self._stats["passes"],
            },
            "genes": self._genes,
            "specific_lookup": specific_out,
            "multi_lookup": self.multi_lookup,
            "conflicts": self.conflicts,
        }

    def build_diagnostic_report(self) -> dict:
        """Build a per-ID-type diagnostic report for reclassification guidance.

        Returns a dict with per-id_type statistics:
          - how many IDs of this type ended up in specific_lookup (unique)
          - how many in multi_lookup (multi-match)
          - how many in conflicts (data errors)
        """
        type_stats: dict[str, dict[str, int]] = defaultdict(
            lambda: {"unique": 0, "multi": 0, "conflict": 0}
        )

        for gene_entry in self._genes.values():
            for id_rec in gene_entry.get("tier1_ids", []):
                t = id_rec["type"]
                id_val = id_rec["id"]
                if id_val in self.conflicts:
                    type_stats[t]["conflict"] += 1
                else:
                    type_stats[t]["unique"] += 1

            for id_rec in gene_entry.get("tier2_ids", []):
                t = id_rec["type"]
                id_val = id_rec["id"]
                n = len(self.multi_lookup.get(id_val, []))
                if n > 1:
                    type_stats[t]["multi"] += 1
                else:
                    type_stats[t]["unique"] += 1

            for id_rec in gene_entry.get("tier3_ids", []):
                t = id_rec["type"]
                id_val = id_rec["id"]
                n = len(self.multi_lookup.get(id_val, []))
                if n > 1:
                    type_stats[t]["multi"] += 1
                else:
                    type_stats[t]["unique"] += 1

        # Add reclassification warnings
        warnings = []
        for id_type, stats in type_stats.items():
            total = stats["unique"] + stats["multi"] + stats["conflict"]
            if total == 0:
                continue
            tier = get_id_tier(id_type)
            multi_pct = 100 * stats["multi"] / total
            if tier == 1 and stats["conflict"] > 0:
                warnings.append(
                    f"[CONFLICT] Tier 1 id_type '{id_type}' has {stats['conflict']} "
                    f"data conflicts — check annotation quality."
                )
            if tier == 1 and multi_pct > 10:
                warnings.append(
                    f"[RECLASSIFY?] id_type '{id_type}' declared Tier 1 but "
                    f"{multi_pct:.0f}% of IDs match multiple genes — consider Tier 2."
                )
            if tier == 2 and stats["unique"] > 0:
                pass  # Expected: most Tier 2 singletons are fine

        return {
            "per_id_type": dict(type_stats),
            "warnings": warnings,
            "unresolved_rows_per_source": dict(self._stats["unresolved_rows_per_source"]),
        }
