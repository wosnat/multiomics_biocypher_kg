# Changelog

All notable changes to the multi-omics knowledge graph are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions use the KG release scheme `X.Y.Z[-(alpha|beta|rc).N]` and are tagged
`kg-X.Y.Z…`.

**Process (accumulate-then-cut):** log notable changes under **[Unreleased]** as
they land. At release time, `/release-kg` *cuts* `[Unreleased]` into a dated
version section, stamps the same version onto `Schema_info.version`, and renders
the GitHub Release notes from that section. The changelog is the source of truth;
the GitHub Release is a rendering of one section. See `plans/alpha_release.md` §2.3.

**Authoring conventions for preflight subsections.** Each version section MAY
include two preflight-facing subsections that are extracted at build time and
stamped onto `Schema_info.release_highlights` / `Schema_info.release_breaking`.
MCP/explorer clients surface them when a user connects, so they answer two
preflight questions:

- `### Highlights` — "what can I now ask that I couldn't before?" One bullet
  per user-facing capability or data layer added in this release. Cap at ~5
  bullets. Plain prose, no internal jargon.
- `### Breaking` — "did anything change meaning under me?" One bullet per
  silently-changed semantics (redefinitions, repurposed fields, default flips).
  Errors/removals are also fine here. **Omit the subsection entirely when there
  are no breaking items** — do not leave an empty `### Breaking` heading
  (renders as a blank bullet on the client).

Both are extracted verbatim (markdown). The rest of the version section
(`### Added` / `### Changed` / `### Fixed`) is unchanged in role.

**Authoring convention for data changes.** The graph's *content* is part of the
release, not just its code, so each version section MAY include a `### Data`
subsection — placed after `### Breaking` and before `### Added` — logging every
publication or dataset added, re-wired, superseded, or corrected. One bullet per
paper/dataset, naming: the paper key as it appears under
`data/Prochlorococcus/papers_and_supp/` (+ DOI for new papers), the organism(s)
and omics type, what the experiments actually compare, and the scale (number of
experiments / analyses / metrics, or edge counts if measured). A correction to an
existing paper says what the KG asserted before and what it asserts now, because
that changes query results. Strain additions belong here too — in `### Data`
regardless of whether code landed with them, with the code side logged separately
under `### Added` / `### Changed`. Unlike `### Highlights` / `### Breaking`, `### Data` is **not**
extracted onto `Schema_info` — omit the heading entirely when a release adds no
data. `/paperconfig`, `/add-a-strain` and `/integrate-a-tool` each end by writing
here, and `/release-kg` preflight warns when paperconfigs changed since the last
tag with nothing logged.

## [Unreleased]

### Highlights

- **TCDB node names (T6).** 916 of the 1,515 `TcdbFamily` nodes rendered as
  their bare TC id — 60% of the transporter ontology was reachable only by
  exact id. Names are now scraped from tcdb.org into a committed
  `cache/data/tcdb/tcdb_names.json` (browse.php + one `result.php` page per
  kept family; 351 families, 0 failures) and folded in by step 6. Unnamed
  drops **916 → 487**, all genuinely unnamed upstream (482 subfamilies TCDB
  never named, 2 specificities, 3 retired families) — those keep the bare id
  by design (`t.name = t.tcdb_id` is the fallback marker). `tc_specificity`
  nodes gain a sparse 400-char `description` (284 nodes; `name` is the
  citation-stripped first sentence, 150-char capped); `tcdbFamilyFullText`
  now covers `name, tcdb_id, superfamily, description`, so "sugar porter" /
  "ammonia" searches reach transporter nodes for the first time. Two class
  names were wrong vs upstream and are fixed (`Channels/Pores`, `Accessory
  Factors Involved in Transport`); the hardcoded `_TC_CLASS_NAMES` table is
  deleted. New standalone scraper
  `multiomics_kg/download/scrape_tcdb_names.py` (resumable, 2.5 s cadence,
  scoped to kept families); step 6 warns with a top-up command when a strain
  onboarding introduces families the scrape hasn't covered. See
  `docs/kg-changes/tcdb-two-source-upgrade.md` ("Node names") and
  `docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md`.

- **InterPro domains.** Every gene now carries direct InterProScan domain/family
  calls: ask "which InterPro entries does gene X have?", "which genes carry the
  RuBisCO large-subunit family?", or filter genes by protein domain — a
  method-independent cross-check on the existing eggNOG/Pfam annotations, with
  domain coordinates and member-DB evidence on each edge, and a curated
  description on every entry node (searchable full-text).
- **NCBIfam families — a new ontology.** Genes carry direct hits from NCBI's
  curated prokaryotic family HMMs (4,957 `NcbifamFamily` nodes, 67K edges,
  47K genes): function-precise identity ("catalase-peroxidase katG", not just
  "peroxidase domain"), comparable across all 42 strains on one library
  version. The katG Black-Queen check works at family level now: present in
  every Alteromonas, absent from every Prochlorococcus. NCBIfam product names
  and HAMAP descriptions also backfill gene naming where UniProt (which is
  removing most of our strains) goes dark.
- **InterPro now enriches gene function, with provenance.** InterPro domain/
  family calls feed each gene's GO, EC, CAZy and Pfam annotations (thousands of
  genes gain a first functional term), and every GO/EC/Pfam/CAZy edge now says
  *who* asserted it (`sources`) and *how strong* the claim is (`evidence`:
  curated vs domain-inferred). Ask "which genes have this EC — from curated
  sources only?" or tell a domain-inferred guess from a reviewed annotation.
- **Transporter substrates now say whether they can be trusted.** Every gene with
  a transporter call carries a `resolved` / `family_inferred` verdict, so you can
  tell "this gene moves nitrate" from "this gene is *some* ABC transporter, and
  ABC transporters collectively move 554 things". Substrate counts per metabolite
  and per gene were corrected to match, and every transported metabolite now
  reports how many distinct transporter systems move it.
- **What *Alteromonas* releases into the medium.** A new publication (Lu 2026)
  adds the EZ55 exudate proteome — which proteins are detected in the >50 kDa
  cell-free supernatant versus the matching whole-cell lysate, for the ancestral
  strain and two lineages evolved under elevated pCO₂. Ask which EZ55 proteins
  are exuded, and whether that set shifted with pCO₂ adaptation.
- **Proteases and their inhibitors are now classified.** Every strain's proteome
  was scanned against MEROPS, the authoritative peptidase classification: ask
  "which subtilisin-family (S08) peptidases does this Alteromonas carry?", count
  an organism's proteases honestly, or join with the signal-peptide and
  localization layers for secreted exoproteases. Each call says up front whether
  it is a real peptidase, a protease *inhibitor*, or a catalytically dead
  look-alike (`call_class`), so dead homologs never inflate protease counts.
- **Ask the graph what a filter accepts, instead of guessing.** Every value a
  property or edge can take — which "evidence" strengths are possible, which
  sources contributed an annotation, which categories a transporter's substrate
  breadth falls into — is now published as data in the graph itself, so a tool
  built against the KG can look the allowed values up rather than hard-coding a
  list that silently goes stale when the graph changes.

### Breaking

- **`Gene_has_interpro_entry` no longer carries `score`.** Member-DB scores are
  incomparable scales (HMMER bits vs HAMAP/PROSITE profile units), so the old
  max-across-libraries value was meaningless whenever two libraries were
  involved. The edge now carries `evalue` (best/min, nullable) + new
  `evalue_library` (which member DB produced it); cross-library confidence is
  `size(libraries)` / `match_count`, by counting not arithmetic.
  `Gene_has_ncbifam_family` (single library, homogeneous HMMER scale) keeps
  both `evalue` and `score`.
- **`annotation_state` semantics tightened.** The quality source buckets grew
  8 → 9 (`ncbifam` added; deliberately NO `interpro` bucket — it is a conduit
  routing evidence into the go/ec/pfam/cazy buckets, and counting it would
  double-count the same HMM hit). `has_any_edge` now includes
  `Gene_has_interpro_entry` + `Gene_has_ncbifam_family`, so 420 genes whose
  only annotation was an InterPro/NCBIfam edge moved `no_evidence` →
  `catch_all_only` (they were misread before); ~184 more climbed into
  informative states via the new bucket. Whole-KG: `no_evidence` 12,481 →
  12,061, `catch_all_only` 5,752 → 5,988.
- **`<strain>.interproscan.calls.json` format changed completely** (faceted:
  per-member-DB `libraries` facets + `interpro_entries` rollups + attributed
  `go_terms`; no `matches[]` list, no pathways) and the per-strain
  `entry_xrefs.json` sidecars were **deleted** (the central
  `interpro_reference.json` replaces them). Any external consumer of the old
  artifact shape must migrate.
- **`metabolite_count` is now the catalysis arm only** on `Gene`, `Metabolite`
  (`gene_count`) and `OrganismTaxon`. It previously unioned catalysis with
  transport, which mixed a p90-of-11 signal with a p90-of-554 one; 23,137 genes
  had transport evidence only, so their number was entirely the inflated arm. The
  transport arm moved to `transported_metabolite_count` (`Gene`, `OrganismTaxon`)
  and `transporter_gene_count` (`Metabolite`). Readers wanting the old union must
  now add the two.
- **`Metabolite.transporter_count` changed definition** from "distinct
  `tc_specificity` nodes" to "distinct transporter systems at maximal depth". It
  was reading 0 for 83% of transported metabolites, so any `transporter_count > 0`
  filter was silently excluding most of them.
- **Substrate queries must use `Tcdb_family_transports_metabolite.substrate_depth
  = 'most_specific'`**, not `level_kind = 'tc_specificity'` — the latter now matches
  only 466 of 11,263 substrate edges.
- **Controlled-vocabulary alignment pass** (`docs/kg-changes/vocabulary-contract.md`).
  `TcdbFamily.is_promiscuous` and `InterproEntry.is_promiscuous` are **deleted** — both
  only ever restated a threshold over a count the node already publishes
  (`level >= 2 AND metabolite_count >= 50`, `gene_count >= 1000` respectively); query
  the predicate directly instead. `evidence_score` is now a **float in `[0,1]`** (was
  int 0–3) on every gene→ontology annotation edge, with paired `signal_count` /
  `signals` so `round(score × signal_count)` recovers the old integer.
  `Gene_has_tcdb_family.tcdb_evidence_score` is renamed `evidence_score` (also now
  float `[0,1]`, was int 0–5) and `Gene.tcdb_best_evidence_score` is renamed
  `tcdb_evidence_score_max`; its three component booleans
  (`agrees_across_sources`/`pfam_corroborated`/`go_corroborated`) are renamed
  `source_agreement`/`pfam_support`/`go_support` and are now meaningful-pair strings,
  not native `bool`. `Tcdb_family_transports_metabolite.substrate_depth` values rename
  `deepest`/`ancestor` → `most_specific`/`inherited`. `sources` values rename
  `interpro` → `interproscan` and `diamond` → `tcdb_diamond` (every value is now the
  `id` of a `DataSource` node). The Layer-A router edges
  (`Interpro_entry_related_to_ec_number` / `_cazy_family`) lose their `ambiguous`
  (native `bool`, uniformly `false` in every shipped build — BioCypher never
  round-trips adapter-emitted `bool`) and `source_db` (a hardcoded literal) properties
  entirely; they now carry no properties. See the doc for full derivation recipes and
  a migration checklist.
- **The Lu EZ55 exudate study moved from its bioRxiv preprint to the published
  AEM paper, and every node id it owns changed with it.** alpha.5 and alpha.6
  shipped it as `Lu 2025` / `10.1101/2025.05.28.656624`; it is now `Lu 2026` /
  `10.1128/aem.00798-26`. Because `Publication`, `Experiment`
  (`{doi}_{experiment_key}`) and `DerivedMetric` ids are all DOI-derived, the
  old ids no longer resolve — anything pinned to them returns empty rather than
  stale. The two **numeric** metric types the preprint reported
  (`exoproteome_detection_replicates`, `whole_cell_detection_replicates`:
  detection counts 0–3 across the three LTPE strains) are **removed from
  `KNOWN_METRIC_TYPES`** and replaced by 6 **boolean** per-strain detection
  metrics, so the same evidence is now queried as presence/absence per named
  lineage instead of a count. A `metric_type` filter on either old name now
  matches nothing.

### Data

- **Lu 2026** (new publication, `10.1128/aem.00798-26`) — *Alteromonas macleodii*
  EZ55 exudate proteomics. Two `compartment` experiments over the same cultures:
  the >50 kDa cell-free supernatant (EV-enriched exudate, `omics_type:
  EXOPROTEOMICS`, `compartment: exoproteome`) and the matching crude cell lysate
  (`PROTEOMICS`, `whole_cell`). Table S1 is a detection table, not a fold-change
  table, so it is wired as two `derived_metrics_table` entries carrying **6
  boolean DerivedMetrics** over 602 proteins — presence/absence in each
  compartment for each of the three LTPE strains (LTPE26 ancestor, LTPE397 and
  LTPE403 evolved 500 generations at 400 / 800 ppm pCO₂; distinct evolved
  lineages, not replicates). **Supersedes the bioRxiv preprint of the same
  study**, which shipped in alpha.5 and alpha.6 as `Lu 2025`
  (`10.1101/2025.05.28.656624`) and is now de-registered — see `### Breaking`.
  The registered-paper count is therefore unchanged at 36. Note the current
  `Lu 2025` directory is **a different paper** (ISME `10.1093/ismejo/wrae259`,
  the same group's pCO₂ coevolution survey): `48ca7103` reused the path after
  the preprint was renamed away. It is staged, not registered.
- **Weissberg 2025** — 4 new RNA-seq contrasts wired onto the existing
  experiments: HOT1A3 coculture-with-MED4 vs axenic at day 18 and day 31,
  MED4 long-term starvation (days 60 + 89 vs log-exponential day 7), and MED4
  coculture-with-HOT1A3 vs axenic at day 18. All resolve via `locus_tag_ncbi`.
- **Per-strain tool artifacts regenerated across all 42 strains** — InterProScan
  re-run with `--goterms --pathways` enabled by default, and the TCDB
  `calls.json` regenerated without derived fields. Artifact refreshes only; the
  schema-side consequences are in `### Added` (InterPro two-layer integration,
  TCDB two-source upgrade).
- **biller 2016 — corrected contrast semantics on the two MIT1002 experiments**
  (filed against this repo by the downstream `multiomics_analysis` consumer,
  ticket 3). Both are *within-coculture* time contrasts
  (24 h and 48 h vs 12 h after addition), but `control` read
  "Co-culture with *Prochlorococcus* NATL2A", dropping the reference timepoint —
  a reader querying the KG alone misread the denominator as a generic coculture
  state or as t0. The 12 h reference is restored, and `experimental_context`
  now records it (indexed in `experimentFullText`, so it is searchable rather
  than recoverable only from the ingestion config). `coculture` also moves from
  `treatment_type` to `background_factors` and the treatment becomes
  `growth_phase`: the study has no axenic *Alteromonas* arm, so a
  `treatment_type = coculture` filter was surfacing two experiments with no
  coculture-vs-axenic handle at all. **`Tests_coculture_with` and
  `coculture_partner` are unaffected** — both are gated on `treatment_organism`.
  A sweep of all 36 registered paperconfigs / 174 experiments found no other
  collapsed time-contrast reference. `validate_paperconfig.py` now accepts a partner
  organism in either `treatment_type` or `background_factors` and warns only
  when it is accounted for in neither.

### Added

- **Annotation-state distribution baseline tool.**
  `tests/kg_validity/capture_annotation_state.py` (`--save` / `--compare`,
  omics-edge-snapshot pattern) captures the Gene `annotation_state` /
  `annotation_quality` / `annotation_types` / `informative_annotation_types`
  distributions (global + per-organism) into the committed
  `annotation_state_baseline.json`, so bucket-movement claims across rebuilds
  (e.g. the 2026-08-17 has_any_edge fix's −420 no_evidence) are reproducible
  from artifacts. First baseline captured from the 2026-08-18 build: 124,751
  genes — no_evidence 12,061 (9.7%) / catch_all_only 5,988 (4.8%) /
  informative_single 12,642 (10.1%) / informative_multi 94,060 (75.4%).
- **`ControlledVocabulary` nodes — the value sets a property or edge can take,
  published as data** (design
  `docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md`, consumer doc
  `docs/kg-changes/vocabulary-contract.md`). One node per (label-or-edge-type,
  property) pair — `applies_to`, `property`, `value_type`, `closed`, `values`,
  `description`, plus sparse `min_value`/`max_value`/`signal_count`/`signals` for
  numeric/score vocabularies — sourced from `config/controlled_vocabularies.yaml`
  (THE source of truth) via `multiomics_kg/utils/controlled_vocab.py` and emitted
  by the node-only `controlled_vocabulary_adapter.py` (`DataSource`-adapter
  pattern). A four-gate test suite (adapter unit checks, `--test`-build CSV
  scan, `slow`-build CSV scan, live-graph) checks that no `closed` vocabulary
  ever emits an undeclared value; the loader's `VOCAB.check()` is wired into
  only one adapter (`tcdb_adapter.py`) today, so this is a detection net across
  test runs, not a build-time guard on every emitter. The reverse
  direction — everything declared was actually observed — is a separate,
  opt-in coverage check (`exhaustive: true`).
  `Schema_info.controlled_vocabularies_hash` (sha256 over the emitted set) lets a
  consumer detect drift between releases instead of discovering it through a
  wrong answer. Five house rules now govern every vocabulary the KG ships: **R1**
  lowercase `snake_case` for KG-minted values, external database terms preserved
  verbatim; **R1b** namespace a value only when it collides across labels; **R2**
  every `sources` value joins a `DataSource` node via `id = 'data_source:' + value`; **R3** don't materialize a
  threshold over an already-stored count; **R4** one score name per concept, on
  one `[0,1]` float scale. **R5 — no native `bool`.** A two-state fact is now
  required to be a categorical string naming both states meaningfully (the
  existing `OrthologGroup.has_cross_genus_members: cross_genus | single_genus`
  is the precedent), because BioCypher does not round-trip adapter-emitted
  `bool` — the one place it shipped was silently `false` on every edge in every
  build ever deployed. `value_type` in the contract admits `string`,
  `string_array`, `float`, `int`, `bool_string`, never `bool`.
- **MEROPS peptidase ontology** (merops-diamond Phase 2; design
  `docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md`, release
  notes `docs/kg-changes/merops-extension.md`). 155 `MeropsFamily` nodes
  (clan → family → subfamily, observed-only; ids `merops.clan:SC` /
  `merops.family:S14`; catalytic type as full words, inhibitor families typed
  `family_type='inhibitor'`) + ~4.2K scored `Gene_has_merops_family` edges
  (`call_class` peptidase|inhibitor|nonpeptidase_homolog, tcdb-parity `tier`,
  `confidence_score`, alignment stats, `best_hit_id`/`best_hit_kind`) + 108
  `Merops_family_is_a_merops_family` edges. New `merops_diamond` merge source
  (`merops_ids`, 9th DataSource node), committed
  `cache/data/merops/merops_reference.json` (prepare_data step 9 sub-builder),
  post-import rollups (`gene_count`/`peptidase_gene_count`/`organism_count`/
  `member_count`), Gene routing `merops_family_count` + `merops_classes`,
  tier-gated `annotation_types`/`informative_annotation_types` += `merops`
  (no annotation_quality bucket). Indexes + `meropsFamilyFullText`.
  **2026-08-18 follow-up** (design
  `docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`):
  new `Merops_family_has_pfam_domain` bridge edge (MeropsFamily → Pfam,
  family-level only, 156 live edges, `member_id_count` property) built
  from MEROPS's own curated `interpro.txt` family→Pfam map, giving the
  single-source `Gene_has_merops_family` edge an independent corroboration
  signal via new post-import property `pfam_support`
  (`corroborated`|`uncorroborated`, same shape as the TCDB analog); three new
  sparse family-level `MeropsFamily` cleavage-specificity properties from
  `Substrate_search.txt` (`cleavage_p1_residues`, `cleavage_summary`,
  `known_cleavage_count` — e.g. "cleaves after Lys (36%) / Arg (34%) - 39567
  known cleavages"), `cleavage_summary` folded into `meropsFamilyFullText`;
  the corresponding MEROPS→GO bridge was evaluated and rejected on
  measurement (all-kingdom rollup too noisy, see `plans/backlog.md`); 5 new
  controlled-vocabulary entries complete the merops set of 10 in
  `config/controlled_vocabularies.yaml`; MEROPS reference build moved from
  prepare_data step 9 to a new step 10 (own log), and prepare_data gained a
  `--rebuild` flag that reruns every derived step (`9 1 2 3 4 5 6 7 8 10`)
  with `--force` in one call.

- **InterPro two-layer integration** (design
  `docs/superpowers/specs/2026-08-10-interpro-two-layer-integration-design.md`).
  - **Layer B** — InterPro entry xrefs propagate into `go_terms` (+45K),
    `ec_numbers` (+9.5K), `cazy_ids` (+642), `pfam_ids` (+14K net-new), noise-gated
    (GO/CAZy: FAMILY+DOMAIN, fold excluded; EC: FAMILY+single-EC only; Pfam:
    direct HMM hits). `alternate_functional_descriptions` gains `[interpro]` names.
  - **Edge provenance** — `sources` (str[] `ncbi|cyanorak|uniprot|eggnog|interproscan`;
    `interpro` renamed `interproscan` — see `### Breaking`),
    `evidence` (`curated`>`signature`>`family_inferred`>`domain_inferred`),
    `evidence_score` (now a float `[0,1]`, advisory — see `### Breaking`) on all six gene→ontology edge types
    (`Gene_involved_in_biological_process`/`_located_in_cellular_component`/
    `_enables_molecular_function`, `Gene_catalyzes_ec_number`, `Gene_has_pfam`,
    `Gene_has_cazy_family`). Backed by per-token `<field>_source`/`<field>_evidence`
    maps in `gene_annotations_merged.json`.
  - **Layer A** — `Interpro_entry_related_to_ec_number` (~6,961) +
    `Interpro_entry_related_to_cazy_family` (~122): a recall-biased **router**
    (weak `related_to` verb; carries no properties — the `ambiguous` bool and
    `source_db` it originally shipped with were both deleted, see `### Breaking`)
    homing the multi-EC/DOMAIN ECs and fold CAZy Layer B refuses to stamp on genes.
    **Not an annotation** — never assign gene function from it.
  - Reference cache (step 9) gains sparse `ec_numbers`/`cazy_ids` per entry.
  - See `docs/kg-changes/interpro-two-layer.md`.
- **InterProScan InterPro-entry ontology** (`/integrate-a-tool` Phase 2, then
  rebuilt by the multi-ontology redesign — see below). `InterproEntry` nodes
  (12,999; `interpro:IPR*`) with `interpro_type`, is-a `level`, and sparse
  curated `description` (400-char abstract; in `interproEntryFullText`);
  `Gene_has_interpro_entry` edges (397,342; 102,895 genes, ~85%) carrying
  domain envelope + best e-value with `evalue_library` attribution + member-DB
  `libraries` + `match_count` (no `score` — see Breaking);
  `Interpro_entry_is_a_interpro_entry` hierarchy; `Pfam_in_interpro_entry`
  bridge linking the eggNOG Pfam layer to InterPro (a link, not a merge).
  Post-import `gene_count`/`organism_count`/`member_count`/`is_promiscuous`;
  scalar + full-text indexes. `Gene.annotation_types` gains `'interpro'` +
  `Gene.interpro_entry_count`. **ORA over InterPro must stratify by
  `(interpro_type, level)`** (type primary). `prepare_data` **step 9**
  (central references, ordered before the merge) + committed
  `interpro_reference.json`; `interproscan` `DataSource` node. No e-value
  cutoff (evidence-only). See `docs/kg-changes/interpro-multi-ontology.md`.
- **NCBIfam family ontology** (multi-ontology redesign,
  `docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`).
  4,957 observed-only `NcbifamFamily` nodes (2,204 TIGR\* + 2,753 NF\*; flat —
  no hierarchy exists in the data; `family_type` is the specificity label,
  71% `equivalog`) with `name`/`gene_symbol`/`description` from the new
  committed `ncbifam_reference.json` (`hmm_PGAP.tsv`, prepare_data step 9).
  67,459 `Gene_has_ncbifam_family` edges (direct HMM hits only; `evalue` +
  `score`, single homogeneous scale) + 2,630 `Ncbifam_family_in_interpro_entry`
  bridges (double-sided dangling guard). Node IDs `ncbifam_TIGR*`/`ncbifam_NF*`
  (underscore — not a bioregistry prefix). Post-import
  `gene_count`/`organism_count`, `Gene.ncbifam_family_count`,
  `annotation_types` + informative buckets gain `'ncbifam'`,
  `is_uninformative` on 195 unknown-function families (126 via the typed
  `hypoth_equivalog` rule — a new third rule kind in
  `config/uninformative_terms.yaml` — + 69 via name patterns), 3 scalar
  indexes + `ncbifamFamilyFullText`.
- **Naming recovery from local scans.** `[ncbifam]`/`[hamap]` tokens join
  `alternate_functional_descriptions` (skipped when they merely echo the
  product — NCBI products are themselves PGAP-derived from NCBIfam), and
  NCBIfam `gene_symbol` fills `gene_name` as a lowest-priority fallback
  (never overwrites; `gene_name_source: ncbifam`). Compensates for UniProt
  removing most project strains from its database.
- **Calls↔merge consistency guard** (`tests/test_interproscan_consistency.py`):
  per-strain exact set-equality between the merged seeds
  (`interpro_entries`, `ncbifam_ids`) and the calls.json facets — a
  re-normalize without re-running the step-2 merge now fails loudly instead
  of silently skewing. Plus a static 9-bucket gate
  (`tests/test_annotation_quality_buckets.py`) pinning the bucket list and
  `.cypher`/`.sh` agreement (the bucket-count test CLAUDE.md referenced but
  which never existed).

### Changed

- **Faceted InterProScan artifacts.** All 42 strains' calls.json re-normalized
  (from cached raw output, no re-scan) to the multi-ontology format: sparse
  per-member-DB `libraries` facets (all 17 DBs, `ipr: null` preserved for
  non-integrated hits), per-protein `interpro_entries` rollups (min-evalue +
  `evalue_library`; count-don't-combine — no cross-library score), `go_terms`
  with donating-entry attribution, no pathways (Reactome/MetaCyc dropped from
  per-strain files; MetaCyc parks in the central reference). New
  `interproscan-run --normalize` re-parse mode (NORMALIZED/NO_RAW/FAILED,
  batch-safe).
- **GO/EC/CAZy gates re-implemented donor-attributed** in the step-2 merge:
  a GO transfers iff ≥1 donating entry is FAMILY/DOMAIN (evidence
  `family_inferred`/`domain_inferred`); EC stays single-EC-FAMILY-only; the
  gates now read the artifact's attribution instead of re-deriving.
- **InterPro reference descriptions ship observed-only** (full corpus would be
  27.1 MB > the 25 MB gate): descriptions cover the 12,999 observed entries +
  ancestors. Corpus-coupling consequence: onboarding a strain ⇒ re-run
  prepare_data step 9, else new entries fail-soft to missing descriptions.
- **prepare_data step 9 = central references** (InterPro + NCBIfam builders)
  and the default step order is now dependency-ordered (`0 9 1 2 …`) — the
  references must exist before the step-2 merge that consumes them.

### Fixed

- **InterPro redesign deferred cleanups (all 7 from
  `plans/interpro_redesign_backlog.md`).** NCBIfam `is_uninformative` DUF
  name-pattern aligned with `config/uninformative_terms.yaml`
  (`.*DUF\d.*` → `.*\bDUF\d.*` in both post-import scripts; verified a no-op
  on the live graph — same 9 nodes match). `DataSource` `interproscan`
  `info_types` now lists the Layer-B enrichment contributions
  (go_terms/ec_numbers/cazy_ids/pfam_ids/alternate_functional_descriptions/
  gene_name) via a new `info_types_extra` YAML override (spec §4.2).
  `interproscan-run --normalize` `sentinel_rate` denominator now matches scan
  mode (FASTA header count, falling back to call count). New multi-xref
  fan-out test pins the parser's shallow-copy contract. Stale mid-branch
  comments trimmed; `enrich_interpro_fields` docstring overclaim reworded;
  `acc_to_ipr` last-write-wins documented in code.
- **Duplicated KEGG reaction cross-references in `kegg_data.json`.** KEGG's
  `/link/pathway/reaction` endpoint serves **both** prefix forms for every link
  (`path:map00220` *and* `path:rn00220` — 19,775 of each, an exact pairing), and
  `_parse_reaction_to_pathways` normalizes both to the same `ko` id, so a plain
  append stored every pathway twice: 2,105 of 2,375 reactions carried 6,408
  duplicate entries. `_parse_reaction_to_compounds` had the same missing dedup
  against ~165 literally duplicated upstream rows (a compound on both sides of a
  reaction, e.g. H+ `C00080`), affecting 45 reactions. Both parsers now dedup
  order-preservingly. **No graph impact** — `Reaction_in_kegg_pathway` (6,408)
  and `Reaction_has_metabolite` (10,149) were already one-edge-per-pair, since
  edge ids collapse at import; this was cache bloat (6,408 spurious lines per
  rebuild) that would have misled any consumer counting occurrences rather than
  distinct values.

- **Dangling gene→EC and Layer-A EC edges.** `EcNumber` nodes are the Expasy
  hierarchy (7,337 ids), but InterPro's entry-level EC xrefs include obsolete
  (`1.2.8.1`) and invalid (`2.8.3.183`) numbers that `normalize_ec` cannot remap,
  so Layer B propagated them onto genes and `neo4j-admin import` skipped 9 edges
  (5 `Gene_catalyzes_ec_number` + 4 `Interpro_entry_related_to_ec_number`).
  `MultiEcAnnotationAdapter` now prunes gene→EC edges to `all_ec_node_ids()` and
  logs what it drops; that set is injected into `MultiInterproAnnotationAdapter`
  as `ec_node_ids` so Layer A prunes against it too (`None` → no Layer-A EC
  edges, mirroring `pfam_node_ids`). No count changes — import was already
  dropping these. `output/import.report` is now empty.

- **`Protein.catalytic_activities` was being split apart by the array delimiter.**
  `build_protein_annotations` sanitised `|` and `'` only on the scalar
  (`_resolve_passthrough`) path, never on the list (`_resolve_passthrough_list`)
  one, and `uniprot_adapter` was the only adapter not applying `_clean_str` at the
  yield point. 239 pipe-bearing UniProt values therefore reached the CSV where
  BioCypher's `|` array delimiter re-split them, producing 385 fragment elements
  across `Protein.catalytic_activities` (e.g. `"Release of an N-terminal amino
  acid, Xaa-"` + `"-Yaa-, in which Xaa is preferably Leu…"` as two entries), and
  1,568 Protein nodes retained raw `'`. Sanitisation now runs on every token in
  both paths, with `_clean_value` in the adapter as defense in depth;
  `protein_annotations.json` regenerated for all 39 taxids (3,002 lines, pure
  character substitution). Gene nodes were never affected.

- **Fabricated UniProt release on ~89K edges.** `uniprot_adapter` hardcoded
  `data_version = "2024_03"` and stamped it on every `Gene_encodes_protein`
  (44,646) and `Protein_belongs_to_organism` (44,515) edge, while the underlying
  `uniprot_raw_data.json` was actually fetched in 2026 — nothing in the download
  path captures a release at all. The `version` property is dropped rather than
  guessed; `source` / `licence` are unchanged. (It never reached Protein *nodes*:
  `version`/`source`/`licence` are not in the `protein` schema whitelist.)

- **TCDB substrate rollup had no depth marker, inflating three scalars.** Step 6
  materialises every descendant's substrates onto every ancestor (deliberately —
  it is what keeps an ancestor substrate-reachable after its leaves are pruned),
  but nothing distinguished "this node is the transporter system" from "this node
  is an ancestor of one". Consequences, all fixed together and folded into the
  existing unreleased TCDB upgrade contract
  ([`docs/kg-changes/tcdb-two-source-upgrade.md`](docs/kg-changes/tcdb-two-source-upgrade.md) §7):
  - `Tcdb_family_transports_metabolite` gains **`substrate_depth`** (`'most_specific'` |
    `'inherited'` — see `### Breaking` for the 2026-08-18 rename from `'deepest'`/
    `'ancestor'`). A (node, substrate) fact — a node can be most-specific
    for one substrate and inherited for another. Categorical string, not bool.
  - `Metabolite.transporter_count` read **0 for 1,218 of 1,462 (83%)** transported
    metabolites: it filtered `level_kind = 'tc_specificity'`, but only 466 of
    11,263 substrate edges now sit there after the ancestor-only prune. Now counts
    distinct `substrate_depth = 'most_specific'` sources — non-zero for all 1,462.
  - The transport arm now counts each gene's **deepest TC attachments only**.
    6,950 genes are annotated at both an ancestor and its own descendant (e.g.
    `3.A.1` *and* `3.A.1.14`) and were inheriting the superfamily's whole rollup
    anyway; p90 drops 554 → 97 while keeping 26,813 of 26,894 genes (99.7%).
  - New `Gene.transport_substrate_resolution` (`'resolved'` 28,405 / p90 33 ·
    `'family_inferred'` 1,671 / p90 554), sparse — absent when the gene has no
    TCDB edge. Answers the friction logged by the Alteromonas coculture analysis
    ("carry a confident-vs-inferred flag on every substrate tag"). Deliberately
    **not** tier-gated: 11,871 `resolved` genes are tier-3-only narrow `2.A.x`
    carriers, exactly what a tier gate would wrongly discard.
  - **Bridge-map fetch failures no longer degrade silently.** `tcdb_pruned.json` is
    committed, so a failed tcdb.org fetch used to be caught, logged at WARNING, and
    written out with empty `pfam_bridge`/`go_bridge` — durably dropping ~8.4K bridge
    edges. Step 6 now raises when the download fails with no cached TSV, and when a
    rebuilt bridge comes out empty beside a populated committed one. Verified the
    happy path is untouched: a full `--force` re-run reproduces `tcdb_pruned.json`
    and `kegg_data.json` byte-identically.
  - **`--test` no longer produces a broken TCDB layer.** `MultiTcdbAnnotationAdapter`
    was the only ontology adapter capping its *node* output at 100 (cazy, interpro,
    psortb and signalp all cap the per-gene edge loop and emit the full ontology),
    so gene/parent/substrate/bridge edges pointed at ~1,400 unemitted nodes and
    `skip_bad_relationships: true` dropped them silently. The ontology is bounded
    reference data; only the 53K gene edges needed capping.
  - **Removed the adapter's dead `_TC_CLASS_NAMES` copy.** `build_tcdb_hierarchy`
    (step 6) owns the authoritative table and names every class, so the fallback
    never fired — and could not have helped the one unnamed class, TC `6`, because
    the adapter's copy had no entry for it. Verified no-op: node-name hashes across
    all 1,515 `TcdbFamily` nodes are unchanged.
  - **`consensus_collapse` now enforces its 5-part assumption** instead of assuming
    it. `parse_tcdb_subject_id` accepts 3-5 part TCIDs and list slicing does not pad,
    so a group of 4-part hits would have matched at depth 5, reported `"5_part"`,
    inflated the agreement weight to 1.0 and mislabelled a subfamily call as
    `tc_specificity`. Latent only — TCDB ships 5-part headers throughout
    (40,520/40,520 candidates across 42 strains).

- **Adapter INFO logging was being discarded.** Twelve adapters use stdlib
  `logging.getLogger(__name__)`; with no handler configured, the root logger's
  `lastResort` fallback emits WARNING and above only, so every INFO diagnostic
  they produced was dropped from the build log — node/edge tallies, TCDB
  seed-alias remap counts, and the substrate depth breakdown. Only warnings
  survived, which is why the safety nets looked healthy while the accounting was
  mute. `create_knowledge_graph.configure_logging()` attaches a handler to the
  `multiomics_kg` package logger with `propagate = False`. Deliberately not
  `logging.basicConfig()`: BioCypher's logger sets `propagate = True` *and*
  attaches its own StreamHandler, so a root handler would print every BioCypher
  record twice (verified: 1x each after the fix).

- **Duplicate `Protein_belongs_to_organism` edges.** The adapter deduped both
  gene- and organism-edges on `(locus_tag, ncbi_acc)`, but the organism edge
  depends only on the assembly — so the 80 proteins mapping to several locus tags
  within one assembly emitted ~131 duplicate organism edges per build, absorbed
  silently by BioCypher's `Deduplicator` (with a `Duplicate edge type
  Protein_belongs_to_organism found` warning). Now deduped on the correct key.
  No graph change: the deduplicator was already collapsing them.

## [0.1.0-alpha.6] - 2026-06-13

### Highlights
- **Preflight release summary in your MCP/explorer client.** Each KG release
  can now ship two short markdown bullets — what's new and what changed
  meaning — that your client surfaces the moment it connects. Lets you
  answer "what can I now ask that I couldn't before?" and "did anything
  silently change under me?" before you start querying, instead of
  discovering a redefined field through a wrong answer.
- **Two more Prochlorococcus strains: MIT1327 (LLIV) and MIT1314 (HLII).**
  The last two strains of the Soussan 2025 N/P-starvation panel are now in the
  graph, so cross-strain queries and ortholog comparisons span the full panel.
  MIT1327 includes curated CyanoRak ortholog clusters and roles.

### Added
- **Prochlorococcus MIT1327 + MIT1314** — the two Soussan 2025 N/P-starvation
  panel strains previously deferred for lack of a public assembly, completing
  the 15-strain panel (47 OrganismTaxon nodes / 40 genome strains total).
  Genome-only (no DE paper). MIT1327 (`GCF_001632125.1`, taxid 1801626, clade
  LLIV) carries full CyanoRak annotation (CK ortholog clusters + CyanoRak/TIGR
  roles) on 2308 genes, bridged from a CyanoRak-team gbff via diamond
  protein-sequence matching because the strain has no public CyanoRak server
  export — see the new reusable converter
  `multiomics_kg/download/convert_cyanorak_gbff_to_gff.py`. MIT1314
  (`GCF_034093315.1`, taxid 3096220, clade HLII) is NCBI/eggNOG-only. Both
  carry eggNOG, PSORTb, SignalP, and tcdb-diamond annotation.
- New optional string properties on the `Schema_info` node:
  `release_highlights` and `release_breaking`. Sourced from `### Highlights`
  and `### Breaking` subsections inside each version's CHANGELOG entry,
  extracted at build time by `extract_preflight_subsection()` in
  `release_kg.py`, and stamped by `post-import.sh` (Group 4). Absent or
  empty subsections → real `null` property (not empty-string), so legacy
  KGs and no-subsection releases are indistinguishable on the wire.
  Existing `kg_release_info` callers pick both up automatically (they
  already return `Schema_info { .* }`).
- CHANGELOG preamble: authoring conventions for the two preflight
  subsections (~5 bullets each, omit `### Breaking` entirely when nothing
  to flag rather than leaving an empty heading).
- New string property `deployment_role` on the `Schema_info` node: the KG
  **self-declares** its deployment identity (`local-dev` / `staging` /
  `production`) so the explorer-side `kg_release_info` preflight can echo it
  instead of guessing from port heuristics. Sourced from `KG_DEPLOYMENT_ROLE`
  (forwarded to the `post-process` container by all three compose files),
  stamped by `post-import.sh` (Group 4); defaults to `local-dev` on a plain
  `docker compose up`. `/release-kg` sets `staging` for its staging stack and
  `production` for the Track A (`--target local`) deploy. The explorer-side
  consumer is separate, out-of-repo work.

### Fixed
- `/release-kg --target local`: two **re-deploy** bugs surfaced cutting
  kg-0.1.0-alpha.6 — the first release to redeploy over a *live* alpha stack
  (alpha.5 was the bootstrap cut, so neither path had run against existing
  containers before).
  - **Blue/green build collided with the live `alpha-deploy`.** Both the
    transient `kg-alpha-build` verify project and the live `kg-alpha` project
    sourced `docker-compose.alpha.yml`, which hardcoded `container_name:
    alpha-deploy`. Container names are daemon-global, so the build project
    couldn't create its verify-deploy while the live blue container was up
    (`Conflict. The container name "/alpha-deploy" is already in use`). Fixed by
    parameterizing all four alpha `container_name`s with `${ALPHA_CONTAINER_SUFFIX:-}`;
    the build project sets `-build` (→ `alpha-deploy-build` etc.) while the live
    flip leaves it unset and owns the canonical `alpha-deploy`.
  - **Stale `staging-deploy` locked the Neo4j volume on staging re-runs.** A
    prior release leaves `staging-deploy` `Up`; its volume lock made the next
    run's `neo4j-admin import` fail with "The database is in use." Phase 5 now
    tears the staging stack down (volume preserved) before building.
- `/release-kg --target local`: two first-cut deploy bugs surfaced while cutting
  kg-0.1.0-alpha.5 to the lab box.
  - **Live flip re-stamped `Schema_info` to `0.0.0-dev`.** `_alpha_flip_live_deploy`
    brought the live `alpha-deploy` up with `docker compose -p kg-alpha up -d deploy`;
    `deploy`'s `depends_on` chain re-ran build→import→post-process in the live
    project, and because that env intentionally omits `KG_RELEASE_VERSION`, the
    re-import re-stamped the release graph with the dev fallback version. Fixed by
    adding `--no-deps` so the live deploy just *mounts* the already-built/verified
    `kg-alpha-<color>` volume. (Also added to the rollback path.)
  - **Explorer provisioning aborted the deploy on a fresh volume.**
    `_alpha_provision_explorer_user` ran `CREATE USER explorer IF NOT EXISTS …` then
    `ALTER USER explorer SET PASSWORD …` to the same value; on a fresh blue/green
    volume the `ALTER`-to-identical-password is rejected by Neo4j ("Old and new
    password cannot be the same"), crashing the deploy after the graph was already
    live. Replaced with a single idempotent `CREATE OR REPLACE USER explorer …`.
- `/release-kg` SKILL.md: un-stubbed the `--target local` documentation (it was
  still described as raising `NotImplementedError`); documented the two flip/provision
  gotchas above.
- `scripts/alpha_firewall.sh`: new operator helper to apply the `DOCKER-USER`
  allowlist restricting the alpha ports (`17474`/`17687`) to a confirmed lab subnet
  (plan §2.6). Parameterized on the CIDR, idempotent (`-C` check before `-I`), and
  refuses to run against the campus-wide `/16`.

## [0.1.0-alpha.5] - 2026-06-09

### Added
- **Publication "discusses" edges** — a narrative literature index linking each
  publication to the genes and KEGG pathways it discusses in prose (regulators,
  model genes, pathways named in text), distinct from the supplementary DE-table
  expression data. Recall-biased *router* ("which papers discuss gene/pathway X?"),
  best-effort, not exhaustive. New relationship types:
  - `Publication_discusses_gene` (Publication → Gene)
  - `Publication_discusses_kegg_pathway` (Publication → KeggTerm, pathway-level)

  Both carry `prominence` (`central` | `peripheral`) and the extraction `evidence`
  quote. ~1,099 gene + ~140 pathway edges across 40 publications.

  Three-stage pipeline (`plans/publication_discusses_edges.md`):
  - **Extract** — `multiomics_kg/extraction/topics/` + `/extract-discussed-topics`
    skill: full-PDF LLM extraction with 15-page chunking and reference-page skipping
    (fits per-request token caps; large PDFs no longer fail). Writes
    `publication_topics/topics.json`. Strain-aware (attributes each gene to one of the
    paper's strains); captures verbatim locus tags; self-reports an
    `uncaptured_identifiers` triage signal.
  - **Resolve** — `prepare_data.sh` **step 8**
    (`multiomics_kg/download/resolve_paper_topics.py`): genes resolved
    per-strain via `gene_id_mapping` (identifiers-first; gene families fan out to one
    edge per member; `all`/`unspecified` mentions resolve in each paper strain);
    pathways via a global `kegg_data.json` lookup (dangling-proof — only resolves to
    KeggTerm nodes that exist). Writes `topics_resolved.json` + a `resolution_report.txt`
    (per-paper stats, method breakdown, truncated-id count, `unresolved_reasons` tally).
  - **Adapter** — `multiomics_kg/adapters/publication_topics_adapter.py`: pure edge
    adapter; source = the paper's Publication node id (DOI from paperconfig or
    PDF-extraction cache); targets `ncbigene:{locus_tag}` / `kegg.pathway:ko*`. Wired
    into `create_knowledge_graph.py` after the omics adapter.
- Post-import rollups: `Publication.discussed_gene_count` / `discussed_pathway_count`,
  `Gene.discussed_in_publication_count` (in both `post-import.cypher` and `.sh`).
- Two `config/schema_config.yaml` edge types for the above.
- Tests: 44 unit tests (resolution, chunking/merge, adapter) + 8 KG-validity tests
  (`tests/kg_validity/test_discuss_edges.py`: endpoint correctness/no-dangling,
  prominence enum, evidence presence, post-import rollups).
- Shared `multiomics_kg/extraction/pdf_utils.py` helpers: `upload_pdf` (lifted from
  cluster extraction), `count_pages`, `find_references_page`, `page_chunks`,
  `write_page_range_pdf`.

### Changed

### Fixed

## [0.1.0-alpha.4] - 2026-06-08

### Added
- Track A infrastructure scaffolding (no behavior yet — these wire
  into `/release-kg --target local` in a follow-up):
  - `docker-compose.alpha.yml` — alpha-stack compose override
    matching the plan §2.4 sketch. Renames containers `alpha-build` /
    `alpha-import` / `alpha-post-process` / `alpha-deploy`, forwards
    `KG_*` env vars to post-process, drops the Biochatter UI, and
    hands the data volume off to a `${ALPHA_DATA_VOLUME}`-selected
    external name (default `kg-alpha-blue`) for the blue/green flip.
  - `.env.alpha.example` — committed template for the operator-only
    secrets (`ALPHA_BIND_IP`, `NEO4J_AUTH`, `ALPHA_EXPLORER_PASSWORD`).
  - `scripts/alpha_up.sh` / `scripts/alpha_down.sh` — operator
    wrappers around the `-p kg-alpha -f … -f …` invocation. Read the
    active color from `.alpha_active_color` (written by the release
    flow) and refuse to run if `.env.alpha` or the marker is missing.
  - `.gitignore` adds `.env.alpha` and `.alpha_active_color`.
- `/release-kg --target local` now implements the Track A lab-box
  deploy (was a stub). Replaces `deploy_local_stub` with the real
  blue/green flip orchestration in `release_kg.py`:
  - Reads `.env.alpha` for `ALPHA_BIND_IP`, `NEO4J_AUTH`,
    `ALPHA_EXPLORER_PASSWORD`; refuses to run on unfilled
    `REPLACE_WITH_…` / `<...>` placeholders.
  - Determines active/inactive colors from `.alpha_active_color`
    (first cut bootstraps into `blue`).
  - Builds the alpha-inactive color via a transient
    `kg-alpha-build` compose project on temp localhost ports;
    verifies Schema_info via `docker exec`. Pytest L2 is skipped
    here because Phase 5 already ran it against the same tagged
    KG on the staging stack (same tag + env = same KG).
  - Flips: stops the live alpha-deploy (if any), brings it up on
    the new color bound to `ALPHA_BIND_IP`. Best-effort rollback
    to the previous color on flip failure.
  - Provisions the shared `explorer` read login idempotently via
    `CREATE USER explorer IF NOT EXISTS … + ALTER … CHANGE NOT
    REQUIRED`.
  - Warns (does not block) if the `DOCKER-USER` chain has no rule
    for `:17474`/`:17687` — the firewall allowlist needs sudo on
    the lab box.
  - Updates `.alpha_active_color`, prints an operator summary
    (Bolt URI, browser URI, distribute-credentials reminder).
  - 18 new unit tests cover the deterministic pieces — color
    rotation, marker round-trip, env-alpha parsing, validation
    (missing keys, unfilled placeholders, malformed NEO4J_AUTH).

### Changed
- Decision (2026-06-06): the alpha runs on **Track A** — the lab box
  at `132.75.249.47` (leadership choice). Aura (was Track B) archived;
  `/release-kg` drops the `--target aura` backend. Plan
  (`plans/alpha_release.md`) and tester guide
  (`docs/kg_mcp_guide.md` §2) reframed accordingly. Remaining
  agnostic-vs-local design split stays intact; `--target local`
  is now implemented (see Added above) and this release is the
  first cut deployed to the lab box via it.

### Fixed

## [0.1.0-alpha.3] - 2026-06-06

### Added
- `/release-kg` Phase 5 now gates the release on the KG validity suite
  passing against the staging stack (`uv run pytest tests/kg_validity/
  --neo4j-url bolt://localhost:27687 -q`, ~73 s, 1012 assertions).
  Catches structural / semantic regressions before Phase 7 publishes a
  GitHub Release. `--skip-kg-tests` flag bypasses for emergencies. On
  failure, the staging stack is left running so the operator can
  inspect.
- `/release-kg` Phase 7 now embeds a "What changed since kg-X.Y.Z"
  diff block in the GitHub Release notes. Compares the prior published
  release's `metadata.json` to the current build's metadata: headline
  Schema_info count deltas (papers / experiments / genes / organisms /
  expression edges) + per-publication expression-edge changes (new /
  changed / removed). The per-publication detail catches the
  net-zero-but-paper-A-lost-paper-B-gained regression class that
  Schema_info totals would hide. `metadata.json` now also carries the
  full `per_publication_edges: {doi: count}` map for the *next*
  release's diff. Soft-fails on first-ever release / older prior
  releases without `metadata.json` — release still publishes.

### Changed
- `/release-kg`'s default `--mcp-min` value now reads from
  `[tool.release-kg].mcp_min_version` in `pyproject.toml` (hard
  fallback `0.1.0` if the file/section is missing). Previously the
  default was a Python constant inside the skill script. Repo-wide
  bumps of the cross-repo explorer-MCP contract now live in
  declarative config; per-release override via `--mcp-min` is
  unchanged.
- Bumped `[tool.release-kg].mcp_min_version` from `0.1.0` to
  `0.1.0a1` (PEP 440 alpha). Matches the current
  `multiomics_explorer` version, so the explorer's
  `kg_release_info()` compat check passes against builds off this
  branch (string equality on `mcp_min_version`).

### Fixed

## [0.1.0-alpha.2] - 2026-06-02

### Added

### Changed

### Fixed
- `docker-compose.staging.yml` now forwards `KG_RELEASE_VERSION` and
  `KG_GIT_*` env vars from the compose process into the `post-process`
  container, so post-import.sh Group 4 stamps `Schema_info.version` with
  the tagged version instead of silently falling back to `0.0.0-dev`.
  The 0.1.0-alpha.1 cut built and deployed successfully but stamped the
  wrong version because of this gap, which is why the tag exists but no
  GitHub Release was published for it; this release is the first cut
  with a correctly-stamped staging verify.

## [0.1.0-alpha.1] - 2026-06-02

### Added
- `Schema_info` release metadata, stamped at post-import (every build): `version`,
  `built_at`, `git_sha`, `git_sha_short`, `git_branch`, `git_dirty`,
  `mcp_min_version`, `release_notes_url`, plus computed counts (`paper_count`,
  `experiment_count`, `gene_count`, `organism_count`, `expression_edge_count`).
  Dev builds stamp `0.0.0-dev`; releases stamp the tagged version. Added in both
  `scripts/post-import.sh` (Group 4) and `scripts/post-import.cypher`.

### Changed

### Fixed
