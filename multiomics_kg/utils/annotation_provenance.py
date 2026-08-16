"""Edge-level provenance + confidence for gene→ontology annotation edges.

Shared by the GO / EC / Pfam / CAZy adapters. Turns the per-token maps that
``build_gene_annotations`` writes into ``gene_annotations_merged.json`` —
``<field>_source`` (``{token: [source, …]}``) and ``<field>_evidence``
(``{token: evidence}``, sparse) — into the three edge properties every
gene→ontology edge now carries:

- ``sources`` (str[]) — who asserted the annotation
  (``ncbi|cyanorak|uniprot|eggnog|interproscan``).
- ``evidence`` (str) — inference strength: ``curated`` > ``signature`` (direct
  Pfam HMM hit) > ``family_inferred`` > ``domain_inferred``. Tokens no source
  labelled inferred default to ``curated`` (they came from a curated source).
- ``evidence_score`` (float in [0,1]) — advisory, never a filter; a ready sort
  key. Multiply by 3 and round to recover the signal count.

See ``docs/superpowers/specs/2026-08-10-interpro-two-layer-integration-design.md``
§5.3. The label stays a coarse relation; confidence lives here.
"""

from __future__ import annotations

_CURATED_SOURCES = {"ncbi", "cyanorak", "uniprot", "eggnog"}

_SIGNAL_COUNT = 3   # module constant, mirrored in controlled_vocabularies.yaml


def annotation_edge_props(gene: dict, field: str, token: str) -> dict:
    """Build ``{sources, evidence, evidence_score}`` for one (gene, term) edge.

    *field* is the merged-annotation field the term lives in (``go_terms``,
    ``ec_numbers``, ``pfam_ids``, ``cazy_ids``); *token* is the term id.
    """
    src_map = gene.get(f"{field}_source") or {}
    ev_map = gene.get(f"{field}_evidence") or {}
    sources = list(src_map.get(token) or [])
    evidence = ev_map.get(token) or "curated"

    # evidence_score (advisory): three independent +1 signals.
    score = 0
    # (1) corroborated by >=2 *independent* sources. Pfam-via-eggNOG and
    #     Pfam-via-InterPro are the SAME signal (InterPro integrates Pfam), so
    #     collapse that dependent pair before counting (spec §9.2).
    effective = set(sources)
    if field == "pfam_ids" and {"eggnog", "interproscan"} <= effective:
        effective.discard("interproscan")
    if len(effective) >= 2:
        score += 1
    # (2) at least one high-trust assertion (curated source or a direct signature).
    if evidence in ("curated", "signature"):
        score += 1
    # (3) not a bare domain-level inference.
    if evidence != "domain_inferred":
        score += 1

    props: dict = {
        "evidence": evidence,
        # R4: normalized to [0,1] so a Pfam score and a TCDB score are not
        # arithmetically comparable-but-wrong. round(x * signal_count)
        # recovers the raw count; NEVER truncate (0.333 * 3 = 0.999).
        "evidence_score": round(score / _SIGNAL_COUNT, 3),
    }
    if sources:
        props["sources"] = sources
    return props
