# InterPro multi-ontology redesign — design spec

**Date:** 2026-08-17
**Status:** approved design, pending implementation plan
**Supersedes (on completion):** `docs/kg-changes/interproscan-extension.md`, `docs/kg-changes/interpro-two-layer.md` (and their underlying specs `2026-07-22-interproscan-domains-design.md` phase-2 portions, `2026-08-10-interpro-two-layer-integration-design.md`)

## 1. Motivation

Three drivers:

1. **The Phase-1 calls.json is IPR-centric and patchy.** Member-DB signatures
   (NCBIfam, HAMAP, PANTHER, Gene3D, …) are buried as rows in a flat
   `matches[]` list; per-protein `go_terms`/`pathways` aggregates carry no
   attribution. Using InterProScan as a *multi-ontology source* — each member
   DB feeding its respective ontology where one exists in the KG, or standing
   ready as an optional new ontology — is not natively supported by the
   artifact shape.
2. **UniProt is removing most of our strains from its database.** UniProt was
   the main `curated`-grade source for GO, product names, and functional
   descriptions. InterProScan runs locally against our own FASTAs — it cannot
   be cleaned out from under us — and becomes the annotation backbone for
   affected strains. HAMAP (UniProt's own microbial ruleset, run locally by
   the scan) and NCBIfam (curated product names + gene symbols) can recover
   naming that UniProt deletes.
3. **The existing two-layer integration grew bolted-on.** Layer A/B semantics
   are sound but were retrofitted onto an artifact that fights them. Full
   restart: nothing InterPro-touched is treated as settled; prior pieces are
   reused only where they fit the new uniform structure. The current
   `InterproEntry` content in the KG is not consumed by the explorer, so
   scrapping and rebuilding is safe.

## 2. Decision log

| Decision | Choice |
|---|---|
| Scrap scope | **Full restart** — calls.json format, merge sources, Layer A/B propagation, `InterproEntry` ontology all redesigned; provenance infrastructure (track_source, evidence maps, `annotation_provenance.py`) is kept as infrastructure |
| Ontology scope | **Core + NCBIfam**: GO → existing GO terms; PFAM → existing `Pfam` nodes; IPR\* → rebuilt `InterproEntry`; NCBIfam → **new `NcbifamFamily` ontology**. All other member DBs stay as evidence facets in calls.json, elevatable later by one merge-config line |
| EC / CAZy / MetaCyc xrefs | **Re-include EC + CAZy** into gene fields with the proven noise gates, rebuilt cleanly at merge time; MetaCyc stays unused; **Reactome dropped everywhere** (noise for marine bacteria; survives only in gitignored raw.json) |
| Architecture | **Approach 1 — faceted pure-signal calls.json + merge-time xref enrichment.** calls.json records only what the tool said, faceted by member DB; reference-derived enrichment happens at merge time from central reference files |
| Central-data principle | Ontology details (names, types, hierarchies, xrefs) live in dedicated **central reference files**; per-strain artifacts carry only per-protein scan results |
| Naming recovery | **Yes** — HAMAP descriptions + NCBIfam product names feed `alternate_functional_descriptions`; NCBIfam `gene_symbol` is a lowest-priority `gene_name` fallback |
| GO gate | **FAMILY + DOMAIN** entries donate GO (folds/superfamilies excluded), labeled `family_inferred` / `domain_inferred` — tighten at query time, not build time (205 MED4 proteins get their only InterPro GO from DOMAIN entries; e.g. hypothetical + MFS domain → transmembrane transporter activity) |
| EC gate | **Single-EC FAMILY only** (multi-EC families are wrong (N−1)/N of the time; domain-ECs don't generalize). Refused xrefs park on ontology-level router edges |
| Direct-attribution rule | Gene-level Pfam/NCBIfam edges come **only from direct HMM hits** (`libraries.PFAM` / `libraries.NCBIFAM`); entry-mediated sibling signatures are never stamped on genes — the entry↔signature association exists only as ontology-level bridges. IPR attachment is direct by construction (exact matched entry only; ancestors are nodes + is-a edges, never gene edges). GO/EC/CAZy are inherently entry-mediated — handled by gates + `evidence` labels |
| Framing | add-a-tool / integrate-a-tool template family: Phase-1 redo under `/interproscan-run` (re-parse of cached `raw.json`, **no re-scan**), Phase 2 through the gene-annotation-merge front door |
| prepare_data step numbers | Deferred to the implementation plan |

## 3. Artifact layer

### 3.1 Per-strain `<strain>.interproscan.calls.json` (new format)

Rebuilt from the cached, gitignored `<strain>.interproscan.raw.json` (present
for all 42 strains) via a `--normalize` re-parse. One record per protein
(WP_ key), faceted by member DB:

```json
"WP_002805124.1": {
  "md5": "…", "match_count": 5,
  "libraries": {
    "PFAM":    [{"accession": "PF02532",  "name": "Photosystem II reaction centre I protein",
                 "ipr": "IPR003686", "start": 1, "end": 36, "evalue": 4.1e-18, "score": 76.3}],
    "NCBIFAM": [{"accession": "NF002735", "name": "photosystem II reaction center protein I",
                 "ipr": null, "start": 1, "end": 38, "evalue": 3.3e-23, "score": 92.7}],
    "HAMAP": ["…"], "PANTHER": ["…"], "SUPERFAMILY": ["…"]
  },
  "interpro_entries": {
    "IPR003686": {"type": "FAMILY", "libraries": ["HAMAP", "PANTHER", "PFAM"],
                  "match_count": 3, "start": 1, "end": 36, "evalue": 4.1e-18, "score": 76.3}
  },
  "go_terms": {"GO:0015979": ["IPR003686"], "GO:0009523": ["IPR003686"]}
}
```

Rules:

- **`libraries` is the primary structure.** All 17 member DBs appear
  *sparsely* (key present only when that DB matched). Signatures that don't
  integrate into an IPR entry keep `ipr: null` — nothing lost. Any member DB
  becomes a merge source later by declaring its facet.
- **`interpro_entries`** is the per-protein rollup the KG adapter consumes:
  aggregated envelope (min start / max end), best evalue (min, nullable), best
  score (max, nullable), contributing libraries, match count. Entry names /
  hierarchy are NOT duplicated here (central reference is authoritative);
  `type` is kept inline as a jq-inspection convenience only.
- **`go_terms` carries entry attribution** — which IPR entries donated each
  term. This is scan-emitted tool output (`--goterms`), so it belongs in the
  per-strain file.
- **Dropped:** the `pathways` block (entry-level, derivable, ~162K Reactome
  refs of per-strain bloat) and the per-strain `entry_xrefs.json` sidecar
  (replaced by the central reference; files deleted).
- `skill_summary.json` keeps its QC role (input_proteins, calls_made,
  parse_failures, per-library distribution, sentinel_rate, wallclock,
  image digest), updated to the new schema.

### 3.2 Central reference files

| File | Status | Contents |
|---|---|---|
| `cache/data/interpro/interpro_reference.json` | exists (prepare_data step 9) — kept as-is | IPR → name, type, parent, level, sparse `go_terms` / `pathways` (MetaCyc-only) / `ec_numbers` / `cazy_ids` |
| `cache/data/ncbifam/ncbifam_reference.json` | **new** — `build_ncbifam_reference.py` | NCBIfam/TIGR accession → product name, `family_type` (equivalog / subfamily / domain / …), `gene_symbol`, sparse EC/GO xrefs. Source: NCBI's published `hmm_PGAP.tsv` (exact columns confirmed at implementation). NCBIfam is flat — no hierarchy; `family_type` is the specificity signal |
| `cache/data/pfam/pfam_reference.json` | exists — unchanged | Pfam names, clans |
| `cache/data/go_terms/go_namespace_cache.json` | exists — unchanged | GO names, namespaces, DAG |

## 4. Merge / prepare_data integration

The `interproscan` source in `config/gene_annotations_config.yaml` grows from
one bookkeeping field to a multi-field source reading the new facets:

| Merged field | How | Evidence label |
|---|---|---|
| `interpro_entries` | passthrough — keys of `interpro_entries` | — (bookkeeping + routing) |
| `ncbifam_ids` | **new field** — accessions from `libraries.NCBIFAM` | `signature` |
| `pfam_ids` | `interproscan` joins the union as a 4th source (`libraries.PFAM` accessions) | `signature` |
| `go_terms` | `interproscan` joins the union — keys of the `go_terms` facet, gated FAMILY + DOMAIN via reference entry types (attribution data is in the artifact) | `family_inferred` / `domain_inferred` |
| `ec_numbers` | merge-time enrichment from `interpro_reference.json`: gene's IPR entries → entry EC xrefs, gated **FAMILY + single-EC** | `family_inferred` |
| `cazy_ids` | same enrichment path, gated FAMILY + DOMAIN (fold excluded) | `family_inferred` / `domain_inferred` |
| `alternate_functional_descriptions` | **naming recovery**: HAMAP signature descriptions (`[hamap] …`) from `libraries.HAMAP` + NCBIfam product names (`[ncbifam] …`) from the reference | — |
| `gene_name` | NCBIfam `gene_symbol` as **lowest-priority fallback** (after ncbi / cyanorak / uniprot / eggnog — never overrides an existing name) | — |

Mechanics carried over unchanged (infrastructure, not the scrapped part):
per-token `<field>_source` / `<field>_evidence` provenance maps,
`track_source`, `annotation_provenance.py` edge-property derivation,
`contributing_sources`, the `interproscan` DataSource node.

prepare_data flow (step numbers assigned in the implementation plan):

1. Phase-1 re-normalize — `interproscan-run --normalize` re-parses cached
   raw.json → new-format calls.json for all 42 strains (signalp `--normalize`
   precedent; no re-scan; machines without raw.json use the committed
   calls.json as-is).
2. NCBIfam reference build — new module `build_ncbifam_reference.py`
   downloads `hmm_PGAP.tsv` → `ncbifam_reference.json` (global central file,
   InterPro-step-9 pattern; `--force` / `--refetch-raw` flags).
3. Merge — `build_gene_annotations.py` (step 2) consumes the fields above.
4. InterPro reference — step 9 stays as-is.

## 5. KG schema

### 5.1 Nodes

| Node | Status | Notes |
|---|---|---|
| `InterproEntry` | rebuilt, same shape | Observed entries + is-a ancestors (pruned). ID `interpro:IPR*`; name/type/level from central reference |
| `NcbifamFamily` | **new** | Observed-only (~4,957 nodes: 2,204 TIGR\*, 2,753 NF\*), flat. ID `ncbifam:TIGR01234` / `ncbifam:NF002735` (underscore fallback if `ncbifam` is not a bioregistry prefix — check at implementation). Properties: `name` (product name), `ncbifam_id`, `family_type`, `gene_symbol` (sparse), `level` (always 0) |
| `Pfam`, GO terms, `EcNumber`, `CazyFamily` | existing, untouched | Gain InterProScan as an additional evidence source on their gene edges |

### 5.2 Gene → ontology edges

| Edge | Status | Evidence properties |
|---|---|---|
| `Gene_has_interpro_entry` | rebuilt, same evidence shape | `start`/`end`, `evalue` (nullable), `score` (nullable), `libraries` (str[]), `match_count`. No e-value cutoff (member DBs pre-apply curated thresholds) |
| `Gene_has_ncbifam_family` | **new** (~68K matches → deduped per (gene, accession)) | `start`/`end`, `evalue`, `score`. Single-source — no `sources[]` |
| `Gene_has_pfam` | existing | `sources` gains `interpro` on direct scan hits; `evidence: signature`. eggNOG-Pfam + InterPro-Pfam still count as ONE source in `evidence_score` (circularity guard) |
| GO ×3, `Gene_catalyzes_ec_number`, `Gene_has_cazy_family` | existing | `sources` / `evidence` / `evidence_score` via provenance maps. EC edges stay pruned to `all_ec_node_ids()` (Expasy universe — dangling-proof) |

### 5.3 Ontology → ontology edges

| Edge | Status |
|---|---|
| `Interpro_entry_is_a_interpro_entry` | rebuilt from reference (is-a, child → parent) |
| `Pfam_in_interpro_entry` | kept — Pfam ↔ InterPro overlap link, never touches genes |
| `Ncbifam_family_in_interpro_entry` | **new**, same pattern — from the `ipr` field on NCBIFAM facet entries; dangling-proof via injected node sets |
| `Interpro_entry_related_to_ec_number` / `_related_to_cazy_family` | kept — router semantics (`ambiguous` flag, `source_db`), homes gate-refused xrefs, pruned to EC/CAZy nodes the gene edges create |

### 5.4 Post-import + routing

- `NcbifamFamily`: computed `gene_count` (direct), `organism_count`; scalar
  indexes (`ncbifam_family_id_idx`, level) + `ncbifamFamilyFullText` on
  `name`, `gene_symbol`.
- `Gene.annotation_types` += `interpro` (as before) **and `ncbifam`**;
  `Gene.interpro_entry_count` kept; **`Gene.ncbifam_family_count` added**.
- `informative_annotation_types` / `annotation_quality` buckets: deliberately
  untouched (redundant with pfam/go/ec; revisit separately).
- NCBIfam's own EC/GO xrefs: stored in the reference but **not propagated to
  genes this round** — future enrichment source in the same pattern.

## 6. Validation & testing

Unit tests (fast suite):

- New pure parser `multiomics_kg/utils/interproscan.py` (raw → faceted
  calls.json): IPR-integrated + non-integrated matches, GO attribution,
  sparse facets, dropped-pathways rule.
- `build_ncbifam_reference.py` parser: fixture `hmm_PGAP.tsv` rows →
  reference entries.
- Merge tests per new/changed field with provenance-map assertions.
- Rewritten `interpro_adapter` tests + new NCBIfam adapter tests: node/edge
  shape, dangling-proofing, string sanitization.

Pipeline validation (integrate-a-tool gate):

1. Re-normalize all 42 strains → `git diff --stat` on calls.json + jq
   spot-checks (spot-check table updated in interproscan-run SKILL.md).
2. `prepare_data --steps 2 --force` (+ new reference step) → per-strain
   `gene_annotations_merged.json` delta review.
3. Full Docker rebuild → `import.report` clean (no skipped relationships).
4. `/omics-edge-snapshot` before/after — expression edges byte-identical.
5. `pytest -m kg` + updated KG-validity assertions (NcbifamFamily counts,
   edge properties, provenance fields) + snapshot regeneration.

Expected count movements (judge the rebuild against predictions):
`NcbifamFamily` ~4,957 nodes + ~60–68K `Gene_has_ncbifam_family` edges new;
`Gene_has_interpro_entry` ~unchanged (~397K); `Gene_has_pfam` source mix
shifts (more `interpro`-tagged); GO/EC/CAZy edges approximately reproduce the
old Layer-B gains; per-strain `entry_xrefs.json` files deleted.

## 7. Scrap / migration list

- Old-format calls.json → overwritten by re-normalize.
- Per-strain `entry_xrefs.json` → deleted (central reference replaces).
- `interpro_adapter.py` internals rewritten for the new format; NCBIfam
  adapter added (new module or same file — implementation choice).
- `interproscan-run` SKILL.md updated: new Output Schema, `--normalize`
  workflow, refreshed spot checks.
- `docs/kg-changes/interproscan-extension.md` + `interpro-two-layer.md` get
  superseded-by headers pointing at the new kg-changes doc (not deleted).
- CLAUDE.md: adapters list, Neo4j labels, key graph facts, data locations.

## 8. Out of scope

- Elevating any further member DB (PANTHER, HAMAP-as-ontology, CDD, folds…) —
  the faceted artifact makes each a one-line merge change later.
- NCBIfam EC/GO xref propagation to genes.
- `informative_annotation_types` / `annotation_quality` bucket changes.
- MCP / explorer surfacing (source/evidence filters, router-mode traversals).
- Any re-scan of InterProScan itself (raw.json is reused as-is).
- MetaCyc/Reactome pathway modeling.
