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

**Research framing.** The driving question is the *Prochlorococcus–Alteromonas
interaction mechanism*, with four candidate mechanisms: cross-feeding,
ammonification (N remineralization), ROS detoxification, and exoenzymatic
proteolytic activity. All four live disproportionately on Alteromonas genes —
the least-annotated organisms in the KG and the ones UniProt is dropping.
This redesign upgrades exactly that substrate: precise family identity for
ROS-detox enzymes (NCBIfam catalase/peroxidase families, per-strain
presence/absence), EC → Reaction → ammonia paths for ammonification,
transporter/enzyme function hints on hypothetical genes for cross-feeding
(DOMAIN-level GO), and one leg of the secreted-protease intersection query
(family identity × SignalP × PSORTb × coculture expression) for
exoproteolysis. MEROPS and antiSMASH are the named follow-ups this framing
points at (§8).

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
| Quality buckets | `ncbifam` joins as the **9th source bucket**; **no `interpro` bucket** (conduit routing evidence into go/ec/pfam/cazy, not an independent evidence kind — bucket-scale circularity); `has_any_edge` gains both new edge types (bug-fix: interpro/ncbifam-only genes read `catch_all_only`, not `no_evidence`). Measured: 185 / 706 gene movements, §4.2 |
| Framing | add-a-tool / integrate-a-tool template family: Phase-1 redo under `/interproscan-run` (re-parse of cached `raw.json`, **no re-scan**), Phase 2 through the gene-annotation-merge front door |
| prepare_data step numbers | Deferred to the implementation plan |

## 3. Artifact layer

### 3.1 Per-strain `<strain>.interproscan.calls.json` (new format)

**Production (Phase-1, outside prepare_data):** rebuilt from the cached,
gitignored `<strain>.interproscan.raw.json` (present for all 42 strains) via
`interproscan-run --normalize` — a re-parse, no re-scan (signalp
`--normalize` precedent). A skill invocation like the scan itself; the
resulting calls.json is committed, and prepare_data starts from it. Machines
without raw.json use the committed calls.json as-is; only normalization
changes require raw. One record per protein (WP_ key), faceted by member DB:

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
                  "match_count": 3, "start": 1, "end": 36,
                  "evalue": 4.1e-18, "evalue_library": "PFAM"}
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
  aggregated envelope (min start / max end), best evalue (min, nullable) with
  `evalue_library` attribution, contributing libraries, match count. Entry
  names / hierarchy are NOT duplicated here (central reference is
  authoritative); `type` is kept inline as a jq-inspection convenience only.
- **Cross-library numeric rule: count, don't combine.** Scores are a scale
  zoo (HMMER bits for Pfam/NCBIfam, HAMAP/PROSITE normalized profile units,
  none at all for PANTHER/SUPERFAMILY/patterns) — the rollup carries **no
  `score`** (the old max-across-libraries mixed incomparable scales;
  per-match scores stay in `libraries`, interpretable within their own
  library). E-values are comparable in *kind* (lower = better) but computed
  by different algorithms against different DB sizes, so min-evalue is kept
  as a labeled "strongest single-library support" heuristic — evidence-only,
  never a filter or arithmetic input. Cross-library confidence is expressed
  by counting: `size(libraries)` = independent corroboration (distinct model
  families agreeing); `match_count` = total matches grouped into this entry
  (includes within-library repetition, e.g. repeat domains). The
  protein-level `match_count` counts all matches including non-integrated
  (`ipr: null`) ones. There is no "IPR selection" — every match carries its
  signature's curated IPR assignment; the rollup groups matches by IPR.
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
| `cache/data/interpro/interpro_reference.json` | exists (prepare_data step 9) — **extended**: `description` added | IPR → name, type, parent, level, **`description` (new: plain-texted first paragraph of the curated `<abstract>` from `interpro.xml.gz`, which step 9 already streams; capped ~500 chars; sparse)**, sparse `go_terms` / `pathways` (MetaCyc-only) / `ec_numbers` / `cazy_ids`. Size policy: cap ~400 chars (functional gist); expected +12–18 MB committed. **This does NOT reopen the pruned-artifact decision** — but if the built file exceeds ~25 MB, the documented fallback is pruning the *descriptions only* to observed entries + ancestors: legal because descriptions are adapter-only metadata (the merge never reads them) and the observed set derives from the committed calls.json corpus, which exists BEFORE prepare_data — no ordering violation. Fallback cost: couples the reference to the strain corpus (new strain → re-run reference step, or fail-soft to missing descriptions), which is why it is the fallback, not the default. **MetaCyc parks here and stops** — not in calls.json, not merged, not in the KG (the KG's pathway layer is KEGG-native and InterPro has zero KEGG xrefs); dormant option for a future elevation |
| `cache/data/ncbifam/ncbifam_reference.json` | **new** — `build_ncbifam_reference.py` | accession → product name, `family_type`, `gene_symbol`, `description` (curated `comment` column — 64% of observed families carry one), `gene_synonyms`, `pmids`, sparse EC/GO xrefs. Source: NCBI's `hmm_PGAP.tsv` (24 columns, verified 2026-08-17). Flat — no hierarchy, **data-verified**: no parent/child column exists in the TSV; `family_type` is a per-node specificity *label*, not a link (even `exception` doesn't record which family it overrides), so `is_a` edges are impossible — hierarchical context comes via the `Ncbifam_family_in_interpro_entry` bridge instead. **family_type vocab, observed distribution (4,947 of our 4,957 accessions in current TSV):** `equivalog` 3,496 (71% — the conserved-function precision type; validates the elevation), `subfamily`/`superfamily` 934, `domain`/`repeat`/`signature`/`_domain` variants 329, `exception` 74 (more-specific override), `hypoth_equivalog`(+`_domain`) 126 (unknown function → §5.5), `PfamEq`/`PfamAutoEq` **0 observed** (InterProScan's NCBIFAM excludes Pfam-wrapped models — no Pfam duplication). The 10 retired accessions absent from the current TSV get name-fallback from the calls.json facet `name`, `family_type` null |
| `cache/data/pfam/pfam_reference.json` | exists — unchanged | Pfam names, clans |
| `cache/data/go_terms/go_namespace_cache.json` | exists — unchanged | GO names, namespaces, DAG |

**Deliberately NO committed pruned artifact** (do not apply the KEGG/TCDB
pattern here). Pruning exists but stays in-memory at adapter time (InterPro:
observed + is-a ancestors; NCBIfam: observed-only, flat). Every reason
KEGG/TCDB earned a committed pruned cache is absent: full references are
already small committed files (vs huge network-fetched raw), pruning is a
millisecond dict filter (vs needing the ~4 GB MNX resolver), one adapter
consumes each (vs three), and the pruned set is a trivial projection of
(reference × merged seeds) — a committed copy would only add a third
skew-able artifact. Structural blocker besides: references must exist BEFORE
the merge (enrichment), while pruning needs the merge's OUTPUT — a committed
pruned file would force an extra post-merge step recomputable at build for
free. Revisit only if several adapters ever need the same pruned tree.
Nothing downloads at merge time — step 2 lazily *loads* committed
references; only the reference-build steps download (raw inputs gitignored).

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
| `alternate_functional_descriptions` | **naming recovery**: HAMAP signature descriptions (`[hamap] …`) from `libraries.HAMAP` + NCBIfam product names (`[ncbifam] …`) from the reference. **Dedup rule:** skip a token that case-insensitively matches the existing `product` or an already-present description — recovery is for disagreement and absence, not echo (NCBI products are themselves PGAP-derived from NCBIfam, see §4.1) | — |
| `gene_name` | NCBIfam `gene_symbol` as **lowest-priority fallback** (after ncbi / cyanorak / uniprot / eggnog — never overrides an existing name) | — |

### 4.1 Relationship to NCBI/PGAP annotation

RefSeq genomes are annotated by PGAP, whose evidence library is the same
`hmm_PGAP` = NCBIfam collection the scan runs — so NCBI `product` strings are
often *conclusions derived from* NCBIfam hits (the GFF `inference=…HMM:NF*`
attribute records the winner sparsely; 54 CDS in MED4; deliberately not
parsed — the scan supersedes it). Running NCBIfam ourselves recovers the
evidence behind those conclusions: all matches (not just PGAP's best), with
e-values/coordinates, on one current library version across all 42 strains
(NCBI annotations are frozen at each assembly's annotation date and mutually
inconsistent in vintage). Consequences: strong product↔NCBIfam-name agreement
is expected (corroboration); disagreements flag annotation-vintage drift;
recovery value concentrates on `hypothetical protein` genes and
UniProt-sourced fields.

### 4.2 Clean-integration invariants

- Shared gene→ontology edges (GO ×3, EC, CAZy, Pfam): one edge per
  (gene, term); InterPro joins via the existing union + `sources[]` /
  `evidence` provenance — no parallel InterPro-specific edge types.
  eggNOG-Pfam + InterPro-Pfam continue to count as ONE independent source in
  `evidence_score` (circularity guard).
- `Gene_has_ncbifam_family` is single-source (verified: eggNOG emits no TIGR
  tokens in its PFAMs column) — no `sources[]`, psortb/signalp precedent.
- `annotation_quality` / `informative_annotation_types` bucket resolution
  (measured against the live KG, 2026-08-17; 15,002 genes at
  no_evidence/catch_all_only):
  - **`ncbifam` becomes the 9th source bucket** (predicate: informative
    non-hypoth `Gene_has_ncbifam_family` edge, §5.5 filter). Independent
    curated library → counting pfam+ncbifam as two sources is honest
    corroboration. Residual lift: 185 low-quality genes climb, plus
    legitimate single→multi promotions. Follows the documented bucket
    maintenance procedure (post-import ×2, CLAUDE.md, bucket-count test).
  - **NO `interpro` bucket — decided, not deferred.** InterPro is a
    *conduit* routing evidence into the go/ec/pfam/cazy buckets, not an
    independent evidence kind: ~85% of genes carry entries, and an entry
    derived from a gene's sole Pfam hit would count the same HMM hit twice,
    mass-promoting single→multi spuriously (bucket-scale circularity).
  - **`has_any_edge` gains both new edge types** (bug-fix): a gene whose
    only annotation is an InterPro entry or NCBIfam family currently reads
    `no_evidence`; it must read `catch_all_only`. The 706 interpro-only
    low-quality genes land there — honest for what are mostly fold-only
    entries with no functional xref.
  - InterPro-sourced GO/EC/CAZy tokens continue to flow into the existing
    go/ec/cazy buckets (the gene genuinely has that evidence).
- The `interproscan` DataSource node's `info_types` is updated to the new
  contribution list (interpro_entries, ncbifam, pfam, go_terms, ec, cazy,
  naming recovery).
- **Merged-JSON authority rule (calls.json is an evidence sidecar).**
  Adapters emit edges ONLY for entries present in the merged seed fields
  (`interpro_entries` / `ncbifam_ids`); calls.json only decorates those edges
  with per-match evidence. Skew (calls.json re-normalized without re-running
  the step-2 merge) therefore fails soft — missing decorations, never
  phantom edges — and is loudly flagged by a dedicated consistency unit test
  (§6), not silently ignored.

Mechanics carried over unchanged (infrastructure, not the scrapped part):
per-token `<field>_source` / `<field>_evidence` provenance maps,
`track_source`, `annotation_provenance.py` edge-property derivation,
`contributing_sources`, the `interproscan` DataSource node.

**Implementation principle: declarative-first.** Express every new field
through `gene_annotations_config.yaml` (sources, union/passthrough rules,
`transform:` functions, `track_source`) as far as the machinery reaches —
`ncbifam_ids`, the pfam/go union additions, and simple facet extraction are
plain config. Bespoke Python is reserved for what the per-value transform
model genuinely cannot express (reference-lookup gating for EC/CAZy/GO
types, the naming-recovery dedup-against-product rule) — and those follow
the `enrich_pfam_fields` post-merge-enrichment precedent, each as a small
named function, with the field *wiring* still declared in the YAML. No
ad-hoc merge logic buried in `build_gene_annotations.py`.

**Ordering principle — enrichment concentrates in step 2.** All gene-level
enrichment (unions, gates, naming recovery) happens in the step-2 merge, so
`gene_annotations_merged.json` is the single seed for everything downstream:
ontology pruning (InterproEntry = observed + is-a ancestors; NcbifamFamily =
observed-only; TCDB/KEGG step-6 precedent) and Gene routing all read the
merged fields, never re-deriving from calls.json. Adapters return to
calls.json only for per-match **edge evidence** and to central references
only for **node metadata**. Both references are **committed** (raw downloads
under `cache/data/{interpro,ncbifam}/raw/` gitignored) — a fresh checkout
builds the KG with no downloads; reference steps rerun only on upstream
releases.

**Two separate pipelines.** Calls.json regeneration is Phase-1 (a
`/interproscan-run --normalize` skill invocation, §3.1) and is NOT a
prepare_data step — prepare_data starts from the committed calls.json, the
same way it starts from committed eggNOG/psortb/signalp artifacts.

prepare_data flow, in dependency order (exact step numbers assigned in the
implementation plan, but the ORDER is normative — **reference builds run
before the merge**; the current "step 9, independent of 0–8" numbering is
wrong as a dependency graph and only works because the reference is
committed; a full `--force --refetch-raw` rebuild in numeric order would
merge against a stale reference):

1. Central reference builds — `build_interpro_reference.py` (current step 9,
   **renumbered to before the merge**) + new `build_ncbifam_reference.py`
   (`hmm_PGAP.tsv` → `ncbifam_reference.json`; `--force` / `--refetch-raw`
   flags, same contract).
2. Merge — `build_gene_annotations.py` (step 2) consumes calls.json facets +
   both references; all gene-level enrichment lands here.
3. Downstream pruning/adapters read `gene_annotations_merged.json` (+
   calls.json for edge evidence, references for node metadata).

## 5. KG schema

### 5.1 Nodes

| Node | Status | Notes |
|---|---|---|
| `InterproEntry` | rebuilt, same shape + `description` | Observed entries + is-a ancestors (pruned). ID `interpro:IPR*`; name/type/level/**description** (sparse — truncated curated abstract) from central reference |
| `NcbifamFamily` | **new** | Observed-only (~4,957 nodes: 2,204 TIGR\*, 2,753 NF\*), flat. ID `ncbifam:TIGR01234` / `ncbifam:NF002735` (underscore fallback if `ncbifam` is not a bioregistry prefix — check at implementation). Properties: `name` (product name), `ncbifam_id`, `family_type`, `gene_symbol` (sparse), `description` (sparse — curated `comment` from the reference, 64% of observed), `level` (always 0). TIGRFAM accessions (`TIGR*`) are absorbed unmodified — enables the future `Ncbifam_family_has_tigr_role` bridge (§8) without migration |
| `Pfam`, GO terms, `EcNumber`, `CazyFamily` | existing, untouched | Gain InterProScan as an additional evidence source on their gene edges |

### 5.2 Gene → ontology edges

| Edge | Status | Evidence properties |
|---|---|---|
| `Gene_has_interpro_entry` | rebuilt | `start`/`end`, `evalue` (nullable) + `evalue_library`, `libraries` (str[]), `match_count`. **No `score` property** (the old one mixed incomparable scales — §3.1 count-don't-combine rule; `size(libraries)` is the corroboration signal). No e-value cutoff (member DBs pre-apply curated thresholds) |
| `Gene_has_ncbifam_family` | **new** (~68K matches → deduped per (gene, accession)) | `start`/`end`, `evalue`, `score` — both kept: single library, homogeneous HMMER bit scale. Single-source — no `sources[]` |
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

- `InterproEntry` keeps its full rollup set (unchanged from current):
  `gene_count` (DIRECT — the correct ORA count), `organism_count`,
  `member_count` (direct child entries), `is_promiscuous` (gene_count ≥
  1000; ~22 ubiquitous entries). Indexes unchanged
  (`interpro_entry_{id,type,level}_idx`); `interproEntryFullText` extended to
  `name`, `description` (OrthologGroup precedent).
- `NcbifamFamily`: computed `gene_count` (direct), `organism_count` —
  flat-ontology pattern (psortb/signalp precedent). **No `member_count`**
  (no hierarchy) and **no `is_promiscuous`** (mean ~14 genes/family; curated
  families lack the rollup-inflation problem the flag exists for). Scalar
  indexes `ncbifam_family_id_idx`, `ncbifam_family_type_idx` (family_type is
  the stratification key, the `interpro_entry_type_idx` analog),
  `ncbifam_family_level_idx` (convention) + `ncbifamFamilyFullText` on
  `name`, `gene_symbol`, `description` (OrthologGroup precedent).
- `Gene.annotation_types` += `interpro` (as before) **and `ncbifam`**;
  `Gene.interpro_entry_count` kept; **`Gene.ncbifam_family_count` added**.
- `informative_annotation_types` / `annotation_quality`: `ncbifam` joins as
  the 9th source bucket; NO `interpro` bucket (conduit, not evidence kind);
  `has_any_edge` gains both new edge types — full resolution + measurements
  in §4.2.
- NCBIfam's own EC/GO xrefs: stored in the reference but **not propagated to
  genes this round** — future enrichment source in the same pattern.

### 5.5 Term informativeness (`is_uninformative` marking)

Both new ontologies plug into the existing post-import F1.1 mechanism
(`is_uninformative = 'true'` sentinel, absent otherwise; vocabulary in
`config/uninformative_terms.yaml`, stamped by post-import Cypher — GO-roots /
COG-S / role / KEGG-KO precedent):

- **`InterproEntry`** — name-pattern rule (mirrors the KEGG KO rule): flag
  entries named `Protein of unknown function*` / `Domain of unknown function
  DUF*` / `Uncharacterised protein family*`.
- **`NcbifamFamily`** — typed rule (cleaner than patterns, an advantage of
  NCBIfam's curation): flag `family_type` ∈ {`hypoth_equivalog`,
  `hypoth_equivalog_domain`} (126 observed nodes, ~2.5%) + the same
  name-pattern fallback (`hypothetical` / `uncharacterized` / `DUF`).

This defines the "InterPro term-informativeness filter" whose absence was the
recorded reason `informative_annotation_types` folding was deferred. The
folding is now RESOLVED (§4.2): the NCBIfam flags gate the new 9th bucket's
predicate directly; the InterPro flags feed explorer/MCP display (no
`interpro` bucket exists to gate — decided, not deferred). Recorded deviation: F1.1's
guiding principle leaves Pfam DUFs unflagged ("conveys family membership");
flagging InterPro/NCBIfam DUF-equivalents deviates in letter, not spirit —
the flag is a sentinel, never a filter on edge emission, and DUF-family
membership stays fully queryable.

### 5.6 Gene-node delta (complete enumeration)

`go_terms` / `ec_numbers` / `pfam_ids` / `cazy_ids` are merged-JSON fields
that materialize as EDGES — the InterPro contributions there never touch the
Gene node. What actually changes on `Gene`:

| Property | Change |
|---|---|
| `annotation_types` | += `ncbifam` token (`interpro` already present) |
| `ncbifam_family_count` | **new** routing count |
| `interpro_entry_count` | kept |
| `gene_name` | may newly be *filled* (never overwritten) by the NCBIfam `gene_symbol` fallback; `gene_summary` shifts with it where that happens |
| `alternate_functional_descriptions` | gains `[hamap]` / `[ncbifam]` entries (deduped, §4) → flows into `geneFullText` |
| `annotation_state` / `annotation_quality` / `informative_annotation_types` | no mechanism change; membership shifts indirectly (§4.2, §6) |
| `contributing_sources` | unchanged (`interproscan` label already existed) |
| everything else (`product`, `function_description`, `protein_family`, identifiers…) | untouched — NCBIfam `gene_synonyms` / `pmids` stay in the reference, deliberately not propagated |


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
- **Calls↔merge consistency test** (`test_interproscan_consistency.py`, fast
  suite): for every strain with both files committed, join genes to proteins
  and assert set-equality between the merged seed fields and the calls.json
  facets — `interpro_entries` == calls entry keys, `ncbifam_ids` == NCBIFAM
  facet accessions (both are ungated passthroughs, so equality is exact;
  gated fields like GO are excluded by design). Catches stale-merge skew in
  either direction; skips strains missing either file.

Pipeline validation (integrate-a-tool gate):

1. Re-normalize all 42 strains → `git diff --stat` on calls.json + jq
   spot-checks (spot-check table updated in interproscan-run SKILL.md).
2. `prepare_data --steps 2 --force` (+ new reference step) → per-strain
   `gene_annotations_merged.json` delta review.
3. Full Docker rebuild → `import.report` clean (no skipped relationships).
4. `/omics-edge-snapshot` before/after — expression edges byte-identical.
5. `pytest -m kg` + updated KG-validity assertions (NcbifamFamily counts,
   edge properties, provenance fields) + snapshot regeneration.

**Spot checks** (three layers; all pre-verified against committed data
2026-08-17; the calls.json layer goes into the interproscan-run SKILL.md
spot-check table, the KG layer into `pytest -m kg` assertions):

*Layer 1 — calls.json (`jq`):*

| Check | Expectation |
|---|---|
| MED4 `WP_002805124.1` (PsbI) | `libraries.PFAM` has `PF02532` with `ipr: IPR003686`; `libraries.NCBIFAM` has `NF002735` with `ipr: null`; `go_terms["GO:0015979"]` attributed to `IPR003686` |
| EZ55 `WP_156086936.1` / MIT1002 `WP_014977393.1` | `libraries.NCBIFAM` includes `TIGR00198` (catalase-peroxidase katG) |
| MED4, whole file | **no** protein carries `TIGR00198` (Black Queen) |
| MED4 `WP_011131653.1` (MsrB) | `libraries.NCBIFAM` includes `TIGR00357` |

*Layer 2 — `gene_annotations_merged.json`:*

| Check | Expectation |
|---|---|
| MED4 ATP synthase c (`WP_002805169.1`) | `ec_numbers` contains `7.1.2.2` with `interpro` in its source map (single-EC FAMILY **accept**) |
| MED4 argininosuccinate lyase (`WP_011131650.1`) | the 5 fumarate-lyase-family ECs from `IPR000362` are **not** interpro-attributed in `ec_numbers` (multi-EC **refuse**; the correct 4.3.2.1 may still arrive from other sources) |
| MsrB | `ncbifam_ids` contains `TIGR00357` |
| naming dedup | at least one strain-level count of skipped `[ncbifam]` tokens > 0 (echo of `product`), asserted non-zero |

*Layer 3 — KG (Cypher, `pytest -m kg`):*

| Check | Expectation |
|---|---|
| katG | `(g)-[:Gene_has_ncbifam_family]->(:NcbifamFamily {ncbifam_id:'TIGR00198'})` exists for EZ55 + MIT1002 genes; **zero** *Prochlorococcus* genes have it (extends the existing katG Black-Queen test to the NCBIfam layer) |
| KT2440 `WP_010953880.1` (hypothetical + MFS domain) | GO edge to `GO:0022857` with `evidence: 'domain_inferred'`, `'interpro'` ∈ `sources` |
| `TIGR00198` node | `family_type: 'equivalog'`, `description` non-null |
| informativeness | ≥1 `NcbifamFamily` with `family_type: 'hypoth_equivalog'` carries `is_uninformative: 'true'`; DUF-named `InterproEntry` flagged |

Expected count movements (judge the rebuild against predictions):
`NcbifamFamily` ~4,957 nodes + ~60–68K `Gene_has_ncbifam_family` edges new;
`Gene_has_interpro_entry` ~unchanged (~397K); `Gene_has_pfam` source mix
shifts (more `interpro`-tagged); GO/EC/CAZy edges approximately reproduce the
old Layer-B gains; `annotation_state`/`annotation_quality` moves per the §4.2
measurements — ~185 low-quality genes climb via the new `ncbifam` bucket,
~706 interpro-only genes move `no_evidence` → `catch_all_only` via the
`has_any_edge` fix, single→multi promotions only where pfam + ncbifam
genuinely co-occur (no mass promotion — the rejected `interpro` bucket was
the mass-promotion risk); per-strain `entry_xrefs.json` files deleted.

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
- **Named follow-ups for the interaction-mechanism research question** (future
  `/add-a-tool` candidates, not this design): **MEROPS** peptidase
  classification (exoproteolytic-activity shortlists — the one real ontology
  gap; BRITE ko01002 + InterPro FAMILY entries approximate it meanwhile) and
  **antiSMASH** BGC detection (siderophore/vitamin cross-feeding currency;
  per-genome region calls).
- **`Ncbifam_family_has_tigr_role` bridge** (follow-up): TIGR\* families →
  existing `TigrRole` nodes via JCVI's frozen `TIGRFAMS_ROLE_LINK` archive.
  Would give non-Cyanorak strains (all Alteromonas) role-category access via
  the 2-hop `gene → NcbifamFamily → TigrRole`, plus a Cyanorak-role QC
  cross-check. Partial by nature (`NF*` families have no roles; taxonomy
  frozen ~2014). Node-level unification of NcbifamFamily and TigrRole is
  rejected — identity ontology vs category taxonomy.
- NCBIfam EC/GO xref propagation to genes.
- MCP / explorer surfacing (source/evidence filters, router-mode traversals).
- Any re-scan of InterProScan itself (raw.json is reused as-is).
- MetaCyc/Reactome pathway modeling.
