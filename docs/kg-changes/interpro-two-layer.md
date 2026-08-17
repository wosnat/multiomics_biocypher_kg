# InterPro two-layer integration — what changed

> **SUPERSEDED (2026-08-17)** by the multi-ontology redesign — see
> [`docs/kg-changes/interpro-multi-ontology.md`](interpro-multi-ontology.md). This document describes the
> PREVIOUS integration (predates the NCBIfam split, the `evalue_library` edge-property change, and the 9-bucket
> `annotation_quality` resolution). The Layer A/B provenance mechanics described below are still accurate and
> unchanged by the redesign — only the surrounding InterPro edge shape and bucket accounting moved. Kept for
> history.

**Design:** `docs/superpowers/specs/2026-08-10-interpro-two-layer-integration-design.md`
**Commits:** Phase 1 `88b20f24`, Phase 2 `b30356cd`, Phase 3 `6c962f50`
**Status:** KG-side complete. MCP/explorer surfacing is a separate follow-up (§ MCP below).

InterPro was promoted from a bookkeeping source (a bare `interpro_entries` list)
to a **contributing functional annotator**, via two complementary layers plus
edge-level provenance. Nothing curated is discarded — cross-references that are
unsafe to stamp on genes are routed to a marked ontology→ontology router instead.

## Layer B — richer gene annotation (gene→ontology)

InterPro entry-level xrefs now propagate into the existing gene ontology fields,
**noise-gated and type-aware** (an entry only donates a term where it is a
specific claim):

| Field | Filter | Aggregate net-new (42 strains) |
|---|---|---|
| `pfam_ids` | direct PFAM signature hits (no inference) | **+13,977** (133,249 corroborating) |
| `go_terms` | FAMILY + DOMAIN (**fold excluded**) | +45,226 |
| `ec_numbers` | **FAMILY + single-EC** only; bare 3-level `3.4.21`→`3.4.21.-` | +9,538 |
| `cazy_ids` | FAMILY + DOMAIN (fold excluded) | +642 |

Multi-EC FAMILY entries and DOMAIN-level ECs are **not** stamped on genes (each
would be ~1/N precise) — they go to Layer A. `alternate_functional_descriptions`
gains `[interpro] <entry name>` entries (FAMILY/DOMAIN only).

## Edge-level provenance + confidence (gene→ontology)

`Gene_involved_in_biological_process` / `_located_in_cellular_component` /
`_enables_molecular_function`, `Gene_catalyzes_ec_number`, `Gene_has_pfam`,
`Gene_has_cazy_family` now carry three properties (`annotation_provenance.py`):

- **`sources`** (str[]) — `ncbi|cyanorak|uniprot|eggnog|interpro`.
- **`evidence`** (str) — strength: `curated` > `signature` (direct Pfam HMM) >
  `family_inferred` > `domain_inferred`. The single field an LLM/ORA should read
  to tell a curated fact from a domain guess. The label stays coarse on purpose.
- **`evidence_score`** (int 0–3, advisory, never a filter) — +1 for ≥2
  *independent* sources (eggNOG-Pfam and InterPro-Pfam count as one — §circularity),
  +1 curated/signature, +1 not a bare domain inference.

Backing this: the gene-annotation merge now records a **per-token provenance map**
(`<field>_source: {token: [source,…]}`) for union fields (`go_terms`,
`ec_numbers`, `pfam_ids`, `cazy_ids`) — `_resolve_union` honours `track_source`
(previously a silent no-op), and `<field>_evidence` maps carry the InterPro
inference strength. `_has_source_label` / `_compute_contributing_sources` read the
map shape.

## Layer A — ontology→ontology router (the "keep everything" half)

Two new edges home the xrefs Layer B refuses:

- `Interpro_entry_related_to_ec_number` (InterproEntry → EcNumber) — **6,961**
- `Interpro_entry_related_to_cazy_family` (InterproEntry → CazyFamily) — **122**

Deliberately weak `related_to` verb + `ambiguous` bool (`true` = multi-term entry
or non-FAMILY type) + `source_db`. **A recall-biased ROUTER, not an annotation** —
read outward from a gene's known family it corroborates; read backward (carries EC
→ therefore is that enzyme) it is low-precision. **Never assign gene function from
these** — that is `Gene_catalyzes_ec_number`. Same contract as
`Publication_discusses_*` and the TCDB bridges.

Dangling-proof by **pruning to the EC/CAZy nodes the gene edges already created**.
Consequence: a DOMAIN/multi-EC entry's EC only gets a router edge if that EC exists
somewhere in the graph. DOMAIN-ECs no gene carries (~9K of the reference's 16K) are
**not** materialised — the simplicity/completeness tradeoff vs. injecting orphan EC
nodes. Revisit with node injection if ontology-level ORA needs the full set.

### EC pruning correction (2026-08-12)

The original pruning target — the EC set self-computed from the merged JSONs, with
no cross-adapter injection — assumed every EC a gene carries has an `EcNumber`
node. It does not. `EcNumber` nodes are the **Expasy** hierarchy (7,337 ids), while
InterPro's entry-level EC xrefs include obsolete numbers (`1.2.8.1` — deleted
before Expasy's transfer map covers it) and invalid ones (`2.8.3.183` on IPR017821;
Expasy's `2.8.3.*` stops at `.28`). `normalize_ec` remaps neither, so Layer B
propagated both onto genes, and `neo4j-admin import` skipped 9 edges: 5
`Gene_catalyzes_ec_number` + 4 `Interpro_entry_related_to_ec_number`.

Fix (both halves needed — Layer A prunes to what the gene edges create, so the gene
edges have to be clean first):

- `MultiEcAnnotationAdapter.all_ec_node_ids()` exposes the Expasy node universe;
  `get_edges` prunes gene→EC edges to it and logs each dropped EC number.
- `create_knowledge_graph` injects that set as
  `MultiInterproAnnotationAdapter(ec_node_ids=…)`, and Layer A intersects
  `observed_ec` with it. `None` → emit no Layer-A EC edges, the same contract as
  `pfam_node_ids`.

No count changes (`Gene_catalyzes_ec_number` 69,026; Layer A 6,854 distinct) —
import was already dropping these. Separately, Layer A emits ~103 duplicate edge
ids (two raw tokens on one entry normalising to the same EC) that dedup at import;
harmless and pre-existing.

## Reference cache (step 9)

`interpro_reference.json` now also carries sparse `ec_numbers` (10,849 entries) +
`cazy_ids` (227) parsed from `interpro.xml.gz` db_xrefs (`db="EC"` / `db="CAZY"`),
via a generic `parse_entry_db_xrefs` one-pass extractor. Stored raw;
normalization is the consumer's. The InterPro reference is now lazily loaded
inside gene-annotation build (step 2), mirroring the Pfam precedent.

## Not changed / deferred

- **Circularity:** eggNOG-Pfam and InterPro-Pfam are the *same* signal (InterPro
  integrates Pfam); `evidence_score` collapses that pair — don't count them as two
  independent sources downstream.
- **`informative_annotation_types` / `annotation_quality`:** InterPro-inferred
  GO/EC/CAZy do count (they enrich genes); the `evidence` edge property preserves
  the curated-vs-inferred distinction. Confirm bucket semantics on next post-import
  review.
- **Deferred (own specs):** NCBIfam `hmm_PGAP.tsv` GO/EC + TigrRole heterotroph
  roles; residue-level sites + PANTHER subfamily GO (raw JSON is perishable);
  MetaCyc pathway layer. Layer A GO was dropped (redundant — is-a ancestors add 5
  entries).

## MCP / explorer follow-up (out of scope here — spec §7)

- Optional `source_filter` / `evidence_filter` (curated-only vs include inferred)
  on ontology-edge tools, reading the new edge `sources` / `evidence`.
- Surface `sources` + `evidence` in `gene_ontology_terms` / `gene_overview`.
- ORA over InterPro must stratify by `(interpro_type, level)` — `interpro_type`
  primary. Expose the 2-hop `gene → entry → EC|CAZy` (`related_to`, `ambiguous`)
  as an explicit opt-in "candidate/router" mode, never default annotation.
