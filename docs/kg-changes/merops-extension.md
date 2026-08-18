# MEROPS peptidase ontology (MeropsFamily nodes)

**Date:** 2026-08-17
**Spec:** [`docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md`](../superpowers/specs/2026-08-17-merops-kg-integration-design.md)
**Phase 1:** [`/merops-diamond`](../../.claude/skills/merops-diamond/SKILL.md) ([design](../superpowers/specs/2026-08-17-merops-diamond-design.md))
**Track:** 3A ontology · hierarchical (clan → family → subfamily) · scored · functional (tier-gated)

## What's changing

The KG gains the authoritative protease classification. `/merops-diamond`'s
per-strain `<strain>.merops.calls.json` (diamond blastp vs. the 5,009-sequence
MEROPS scan library, release 12.5) now flows through the gene-annotation merge
as the `merops_diamond` source (field `merops_ids: str[]`, join `protein_id`),
and a new `merops_adapter` materializes an observed-only `MeropsFamily`
ontology with scored `Gene_has_merops_family` edges. Single evidence source —
eggNOG does not emit MEROPS ids and InterProScan has no MEROPS member DB.

The central design concern is **making three easy-to-miss distinctions
un-missable to an LLM reader**, at each read-point:

1. **On the edge — `call_class`** (`peptidase` | `inhibitor` |
   `nonpeptidase_homolog`): MEROPS files catalytically dead relatives inside
   peptidase families under `.9xx` identifiers; a gene whose best hits are
   those is fold evidence, NOT protease evidence. Threshold-free (derived from
   `best_hit_kind` + the I-family flag), `expression_status` pattern.
2. **On the Gene — `merops_classes: str[]`**: distinct `call_class` values, so
   `gene_details` never shows a bare count that reads as "N proteases" when
   the calls are dead homologs or inhibitors.
3. **On the family node — `peptidase_gene_count` vs `gene_count`**: the
   metabolite-count-split lesson. C26/C44-style families are dominated by dead
   relatives; the gap between the two counts shows it before anyone does ORA
   or per-organism protease counting with the wrong number.

## New node label: `MeropsFamily`

- **155 nodes** (observed-only + ancestors, CAZy pruning): 41 `merops_clan`
  (level 0) · 97 `merops_family` (level 1) · 17 `merops_subfamily` (level 2).
- IDs: `merops.clan:SC` / `merops.family:S14` / `merops.family:S08A` — both
  prefixes registered in bioregistry → colon CURIEs.
- **No identifier-level (S01.001) nodes** — only 6 tier-1 calls exist in the
  whole corpus; they attach at subfamily/family and the identifier survives as
  edge `best_hit_id` (for tier-1 calls identical to the called code).
- Properties (adapter): `name` (family/subfamily: MEROPS type-example name —
  family.txt where curated, else derived from the scan library's `.001`/lowest
  holotype, organism suffix stripped; 100% of observed families named; clans:
  code), `description` (sparse, clans only — clan.txt fold/mechanism text),
  `merops_id`, `level`, `level_kind`, `family_type`
  (`peptidase`|`inhibitor`), `catalytic_type` (**full words**: serine,
  cysteine, metallo, aspartic, threonine, glutamic, asparagine_lyase, mixed,
  unknown; null on inhibitors — one-letter MEROPS codes are insider jargon).
- Computed (post-import): `gene_count` + `organism_count` (subtree `*0..`,
  CAZy pattern), `peptidase_gene_count` (subtree, `call_class='peptidase'`
  edges only), `member_count` (direct children).
- Node metadata source: **committed `cache/data/merops/merops_reference.json`**
  (594 families / 139 subfamily names / 123 clans), built by prepare_data
  **step 10** (moved from step 9 on 2026-08-18 — see the dated section below)
  (`build_merops_reference.py`) from `family.txt` + `clan.txt` +
  `merops_scan.lib` — reads `$MEROPS_DATA_DIR/DB/` when present (Phase-1
  install), else self-downloads (~2 MB). KG build itself needs zero network.

## New edge types

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_merops_family` | Gene → MeropsFamily | ~4.2K (≤4,254 candidates; one per candidate, multi-domain proteins fan out) | `call_class`, `tier` (1–3), `confidence_score` (0–1), `identity`, `qcov`, `evalue`, `consensus_n`, `best_hit_id`, `best_hit_kind` (`holotype`\|`putative`\|`nonpeptidase_homolog`) |
| `Merops_family_is_a_merops_family` | MeropsFamily → MeropsFamily (child → parent) | 108 | none |

- Edges attach at the **called depth** (tier-truncated: tier 1 = identifier
  claim → attaches at subfamily/family; tier 2 = subfamily; tier 3 = family).
  Families MEROPS leaves clan-unassigned are hierarchy roots.
- Expected `call_class` split (42 strains): peptidase 3,425 ·
  nonpeptidase_homolog 777 · inhibitor 52. Tier shape: ~92% tier 3 (remote
  homology, conservative-by-design — not a defect), ~8% tier 2, 6 tier-1.
- Deliberately NOT on the edge (stay in committed calls.json): `scov`,
  `length`, `consensus_agreement`, `homolog_hit_fraction` (redundant once
  `call_class` exists).

## Cross-node rollups

- **Functional, TIER-GATED (tcdb precedent):** `'merops'` folds into
  `Gene.annotation_types` and `Gene.informative_annotation_types` **only for
  genes with a tier ≤ 2 edge** (~330 genes) — tier-3 calls stay findable
  through the edges without inflating annotation coverage.
- **NO `annotation_quality` bucket** (bucket count stays 9): single evidence
  source, tier-3-dominated — the Phase-1 spec's explicit call; BRITE precedent
  for informative-without-bucket. Revisit if a second MEROPS evidence source
  lands (`interpro.txt` Pfam bridge / Foldseek).
- `Gene_has_merops_family` **does** count toward `has_any_edge`
  (`catch_all_only` vs `no_evidence`) — interpro/tcdb-tier-3 precedent.
- Gene routing signals (post-import): `merops_family_count` (int, ALL edges,
  ungated), `merops_classes` (str[], sorted distinct `call_class` values).
- `contributing_sources` += `'merops_diamond'` (build-time, via the merge);
  9th `DataSource` node (`data_source:merops_diamond`, "MEROPS (diamond)").
- `OrganismTaxon` / `Publication`: no per-ontology rollup (unchanged).

## Indexes

Scalar: `merops_family_level_idx`, `merops_family_level_kind_idx`,
`merops_family_id_idx` (on `merops_id`). Full-text: `meropsFamilyFullText` on
`name`, `merops_id`, `description`.

## Validation results

- Unit: `tests/test_merops_adapter.py` (vocab helpers, node/edge shape,
  tier-1 attachment, call_class, reference-missing failure) — 27 tests; full
  fast suite 2,318 passed.
- MED4 real-data smoke: 55 nodes / 33 is-a / 59 gene edges (= 59 candidates);
  ClpP (PMM1313) → `merops.family:S14` `{call_class: peptidase, tier: 3,
  best_hit_id: S14.008}`.
- `/omics-edge-snapshot`: `pre_merops_integration` baseline captured;
  compare after the Docker rebuild (expression edges must be unchanged).
- `pytest -m kg`: `tests/kg_validity/test_merops.py` added (node counts +
  level distribution, CURIE form, S14/C26 spot checks, hierarchy sanity,
  call_class/tier/score ranges, tier-gate consistency, rollup invariants,
  MED4-vs-MIT1002 lifestyle gradient) — runs against the rebuilt graph.
- `snapshot_data.json`: regenerate after the rebuild.

## What this enables

- "Every subtilisin-family gene across Alteromonas, ranked by confidence":
  `MATCH (g:Gene)-[r:Gene_has_merops_family]->(m:MeropsFamily {merops_id:'S08'})
  WHERE r.call_class = 'peptidase' RETURN g, r ORDER BY r.confidence_score DESC`
- Honest per-organism protease counting via `peptidase_gene_count` /
  `call_class`, joined with the SignalP + PSORTb layers for the secreted
  exoprotease story (heterotroph DOM degradation).
- Protease-inhibitor biology via `family_type = 'inhibitor'`.

## What does NOT change

- Phase-1 `/merops-diamond` runner + artifacts untouched.
- No existing node/edge/property modified — additions-only, EXCEPT:
  `annotation_state`/`annotation_quality` can shift `no_evidence` →
  `catch_all_only` for genes whose only functional edge is a merops call
  (deliberate, has_any_edge extension).
- **MCP / explorer integration is out of scope** — registering `merops` in the
  explorer's `ONTOLOGY_CONFIG` + ontology enum is an explorer-repo follow-up.
  Until then, `run_cypher` reaches everything.
- Deferred (recorded): MEROPS `GO_annotation.txt` GO bridges — measured and
  rejected 2026-08-18, see the dated section below and `plans/backlog.md`.
  The `interpro.txt` Pfam bridge and `Substrate_search.txt`
  cleavage-specificity node properties this line originally deferred landed
  2026-08-18 — see below.

---

## 2026-08-18 — Pfam bridge + pfam_support + cleavage specificity

**Spec:** [`docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`](../superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md)
**Backlog items closed:** "MEROPS cross-ontology bridges", "MEROPS cleavage
specificity as node properties" (`plans/backlog.md`)

This follow-up gives the single-source MEROPS layer (a) an independent
corroboration signal via MEROPS's own curated family→Pfam map, and (b)
functional depth on the family nodes via curated cleavage-site specificity.

### New edge: `Merops_family_has_pfam_domain`

- `MeropsFamily → Pfam`, **family-level attachment only** (level-1 nodes, 95
  of the 97 observed families) — subfamily attachment adds nothing at a
  median of 1 Pfam per family.
- Source: MEROPS's own `interpro.txt` (182 MB, one row per curated UniProt
  member of a MEROPS subfamily/identifier with its Pfam domains). The
  committed reference's `pfam_bridge` block covers the **full MEROPS
  catalog** — 321 families / 440 distinct (family, Pfam) pairs — not
  pre-pruned to our observed families; the adapter prunes at edge-build time
  to observed `MeropsFamily` nodes + the injected `pfam_node_ids`
  (`MultiPfamAnnotationAdapter.all_pfam_ids()`, TCDB-bridge dangling-proof
  pattern; `pfam_node_ids=None` → no bridge edges). Against our 97 observed
  families, `interpro.txt` collapses to **~183 distinct (family, Pfam)
  pairs** (62 families have exactly 1 Pfam, 18 have 2; the S09/C26 tail is
  shared-fold biology — α/β-hydrolase, GATase — not noise). **Measured live
  edge count: 156** — the 27-pair gap is observed-family pairs whose Pfam
  has no node in the graph (pruned by the injected `pfam_node_ids`, as
  designed).
- Property: `member_id_count` (int) — how many distinct curated MEROPS
  identifiers back the (family, Pfam) pair. **Deliberate deviation** from
  TCDB's `curated_tcids` list: MEROPS backing lists run to hundreds of
  identifiers (every characterized S01 member), so the count carries the
  support signal without a bloated array.
- **Verb: `has`** (user-confirmed) — same claim shape as
  `Tcdb_family_has_pfam_domain`: *"MEROPS's curated member proteins of this
  family carry this domain"*, a factual composition statement. `related_to`
  stays reserved for weak recall-biased routers (InterPro Layer A).
- **Read-direction warning (documented, not encoded):** forward (gene's
  family known → do its domains agree?) is the sound direction and the only
  one `pfam_support` uses. Backward (has PF00117 → is a C26 peptidase) is
  unsafe in exactly the shared-fold families where `call_class` already
  flags dead homologs (`nonpeptidase_homolog`) — the same biology, two
  coherent guards.
- Spot check: `merops.family:S14` → `PF00574` (ClpP).

### `pfam_support` on `Gene_has_merops_family`

- Post-import property, R5 string pair: `corroborated` | `uncorroborated`.
  Same name and shape as the TCDB analog
  (`Gene_has_tcdb_family.pfam_support`) — no native bool, per the vocabulary
  contract.
- Definition (sound direction only): `corroborated` ⇔ the gene has a
  `Gene_has_pfam` edge to a Pfam that this edge's family — the target node,
  or its level-1 ancestor when the edge attaches at a subfamily — bridges to
  via `Merops_family_has_pfam_domain`.
- **Advisory, never a filter.** No composite score: with only two signals
  (`tier`, `pfam_support`) a synthesized score adds nothing over reading
  both. The no-`annotation_quality`-bucket decision for MEROPS is unchanged
  — corroboration is not a second independent evidence source.
- **Measured `pfam_support` split by call_class (2026-08-18 rebuild):**
  `peptidase` 3,285 corroborated / 142 uncorroborated (95.9%) · `inhibitor`
  52 / 0 (100%) · `nonpeptidase_homolog` 768 / 10 (**98.7%** corroborated).
  The design predicted corroboration would be *depressed* on dead-homolog
  edges — it is not, and the reason is instructive: **the Pfam bridge
  corroborates the FOLD assignment, which catalytically dead relatives share
  with the active members of their family.** So `pfam_support` and
  `call_class` are orthogonal, complementary signals, not redundant ones:
  `pfam_support = corroborated` says "the family/fold placement is right";
  `call_class = nonpeptidase_homolog` says "but it is probably not an active
  peptidase." Read both — a corroborated dead homolog is still not a
  protease.

### Cleavage specificity — three sparse `MeropsFamily` properties

On family-level (level-1) nodes; **sparse** — absent when MEROPS has no
cleavage data for the family (≤88 of 97 observed families covered; the
committed reference's `cleavage` block itself spans 256 families across the
full MEROPS catalog, not pre-pruned), never an empty value.

| Property | Type | Content |
|---|---|---|
| `cleavage_p1_residues` | str[] | Top P1 residues (residue N-terminal of the cut — Schechter–Berger nomenclature, the dominant specificity determinant) with ≥10% share of curated cleavages, max 3, frequency order. Closed to the 20 standard amino-acid three-letter codes (the raw file's modified-residue codes, e.g. `TyI`, are filtered out so the vocabulary can be closed). |
| `cleavage_summary` | str | Readable sentence, e.g. `"cleaves after Lys (36%) / Arg (34%) / Glu (11%) - 39567 known cleavages (14% physiological)"` for S01 — no `\|` or `'` (the clean-string constraint applies to computed strings too). Joins `meropsFamilyFullText` (index drop+recreate). |
| `known_cleavage_count` | int | Number of curated cleavage records backing the profile ("known cleavages" is MEROPS's own phrasing; the `CLE*` record ids literally mean cleavage). S01: 39,567. |

Source: MEROPS's own `Substrate_search.txt` (29 MB, ~100K curated cleavage
records — enzyme identifier, substrate, position, P4→P4' residues,
physiological/non-physiological/synthetic flag), aggregated over **all**
record types (synthetic assays are legitimate specificity evidence; the
physiological share is *reported* in the summary, not used as a filter).

**Honesty caveats carried by the design** (read before treating this as a
per-gene claim):

1. Support varies four orders of magnitude across families (39,567 events
   for S01 down to a handful) — hence `known_cleavage_count` is exposed
   rather than hidden.
2. Family-level aggregation blends members — S01 blends trypsin's Lys/Arg
   specificity with chymotrypsin's aromatic-residue specificity. The
   property reads "what this family is known to cut," **never** "what your
   gene cuts."
3. Cleavage records are not independent evidence — one degradomics paper can
   contribute many rows to the same family's profile.

### Vocabulary registration

`config/controlled_vocabularies.yaml` now carries **10** MEROPS entries
total. **5 were already registered** by a parallel change before this work
landed — `Gene_has_merops_family.call_class`, `Gene_has_merops_family.best_hit_kind`,
`MeropsFamily.level_kind`, `MeropsFamily.family_type`,
`MeropsFamily.catalytic_type`. **This work adds the remaining 5**, completing
the set:

| Entry | Shape |
|---|---|
| `Gene_has_merops_family.tier` | int, 1–3 (NOT sparse — unlike tcdb, every merops edge carries tier) |
| `Gene_has_merops_family.pfam_support` | string, closed: `corroborated`, `uncorroborated` |
| `MeropsFamily.cleavage_p1_residues` | string_array, closed: the 20 standard amino-acid three-letter codes, sparse |
| `Merops_family_has_pfam_domain.member_id_count` | int, min 1 |
| `Gene.merops_classes` | string_array, closed (3 call_class values), default `[]` |

### Measured GO rejection

A parallel MEROPS→GO bridge (`GO_annotation.txt`, ~1 MB, per-identifier GO
terms from member proteins across **all kingdoms**) was evaluated and
**rejected** (user decision, 2026-08-18) — not implemented:

- 81/97 observed families have GO data, but the all-kingdom member rollup
  gives a **median of 19, max of 389** GO terms per family, including
  biologically-wrong eukaryote-only terms.
- The completeness win is small: only **338** of 4,197 merops-called genes
  lack GO entirely (and only 1 is weakly annotated,
  `annotation_state <= catch_all_only`) — versus the **1,311** genes that
  justified TCDB's GO bridge.
- Recorded in `plans/backlog.md` as a rejected-on-measurement item, revisit
  only with a concrete use case + a filtering design (e.g. ≥ N supporting
  identifiers).

### prepare_data: step 9 → step 10, `--rebuild`

- **MEROPS reference moved from step 9 to a new step 10.** Step 9's identity
  is "central references consumed by the step-2 merge" (the reason it runs
  before step 2 in the default order); the MEROPS reference is consumed only
  at KG-build time (`merops_adapter`), never by the merge, so keeping it in
  step 9 muddied that dependency story — and now adds a 182 MB download to a
  step people run for merge freshness. Default `STEPS` order:
  `"0 9 1 2 3 4 5 6 7 8 10"`. `--steps 10` reruns MEROPS reference alone;
  own log `logs/prepare_data_step10.log`.
- **New `--rebuild` flag**: reruns every derived step with `--force` —
  equivalent to `--steps 9 1 2 3 4 5 6 7 8 10 --force` (dependency order: 9
  before 2; step 0 downloads deliberately excluded). Mutually exclusive with
  `--steps`; composes with `--strains` / `--skip-cyanorak` / `--refetch-raw`
  as usual.
