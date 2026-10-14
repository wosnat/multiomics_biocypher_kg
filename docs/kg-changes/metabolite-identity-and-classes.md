# KG Change: one node per metabolite, KEGG compound classes, chemical roles

**Date:** 2026-10-13 (proposal) · 2026-10-14 (final design)
**Status:** **explorer-reviewed and signed off** (§32 / §34). Implementation is pending, and
the change ships inside a bundled release. Numbers below were measured against the live
2026-10-13 06:28 build; the post-rebuild version of this doc replaces them with measured values.
**Design spec:** `docs/superpowers/specs/2026-10-13-metabolite-identity-and-classes-design.md`
**Origin:** `plans/backlog.md`: metabolite equivalence (P1) + chemical categories (P2), from
the `pro_n_import_coculture` analysis friction log (2026-10-07 / -08 / -11).

This top section describes **what ships**. Sections 3–36 (the appendix) are the review
history, kept for the reasoning behind each decision.

## Summary

1. **Metabolite identity is canonicalised.**
   - One node per metabolite; KEGG id preferred, then primary ChEBI, then MNX.
   - ChEBI ids become primary ids.
   - **Some `Metabolite.id` values change, and `chebi_id` values change on most nodes.**
2. **Paper metabolites resolve deterministically.**
   - No ambiguous pick.
   - Names match case-insensitively, except formulas.
   - A stereo-unqualified name goes to its dominant stereo form, with the assumption recorded
     (`stereo_assignment`).
3. **New `MetaboliteClass` layer** from KEGG's compound BRITE trees (curated only). It's
   separate from `BriteCategory`, so gene-side BRITE tools are unaffected.
4. **Role terms become a `ChemicalRole` hierarchy.** The 57 ChEBI roles TCDB listed as
   substrates ("drug", "siderophore") stop being metabolites. They become `chebi.role:` nodes in
   a pruned role hierarchy, reached by `Tcdb_family_transports_chemical_role`.
5. **Stereo variants are linked, not merged:** `Metabolite_is_stereo_form_of` and
   `Metabolite_is_enantiomer_of`.
6. **The TCDB substrate-lumping set is pinned** (`TcdbFamily.substrate_breadth ∈ {lumping,
   specific, not_assessed}`), so no gene's
   `transport_substrate_resolution` changes from this work. The rest of the TCDB work
   (plausibility, `substrate_depth` semantics, per-organism reach) is **deferred**.
7. **Pipeline:** step 6 now also writes the paper `_resolved.csv` files; step 7 is retired.

## 1. Metabolite identity

### What changes

| | Before | After |
|---|---|---|
| Duplicates | One metabolite can be two nodes: nitrate is `chebi:14654` and `kegg.compound:C00244`; HCO₃⁻ is `C00288` and `C01353` | One node per metabolite (69 merge groups) |
| `chebi_id` | 48% are secondary / obsolete ids; 331 more name an isotopologue ("succinic acid-d4") | The primary entity id, from ChEBI's own KEGG cross-reference |
| `chebi:` node ids | 173 are secondary ids | All primary; 89 re-keyed to KEGG ids |
| Role terms | 57 "metabolites" such as `drug`, `xenobiotic`, `acid`, `vitamin (role)` | Removed from `Metabolite`. Now `ChemicalRole` nodes (§2a) |
| Wrong paper picks | Alanine → D-Ala, Glutamic acid → D-Glu, GSH → GSSG, Cystine → indole-3-acetate, …; 4 `mnx:` placeholders | Curated or dominant-form resolution (§2b) |
| Measurement-only nodes | Nodes that exist only because a wrong pick pulled them in (D-Trp, ChEBI vitamin-B nodes, …) | Dropped when their rows move (≈ 11, incl. the 4 placeholders), listed in the remap and in the target's `recurated_from_ids` |
| Class / polymer nodes | "alditol" carries glycerol's structure, chitosan carries glucosamine's | `structure_kind = generic` / `polymer`. The concrete representative structure is dropped; generic formula and `elements` are kept |

**What "the same metabolite" means:**
- **Merged:** the same MNX entity; a ChEBI term and its KEGG cross-reference; conjugate
  acid / base pairs; tautomers.
- **Not merged:** stereoisomers. L-Glu, D-Glu and stereo-unspecified Glu stay three nodes, as in
  KEGG, and are linked instead (§2b).
- **Never merged automatically:** name-only matches. Two KEGG ids merge only on MNX evidence or a
  curated `merge:` override.

**Scale (approximate):** `Metabolite` 3,356 → **≈ 3,214**:
- −74 merged;
- −57 roles;
- ≈ −11 measurement-only nodes (incl. 4 placeholders);
- 89 re-keyed.

### New / changed Metabolite properties

| Property | Type | Notes |
|---|---|---|
| `id` | str | Canonical CURIE. May differ from the old id; see `metabolite_id_remap.tsv` |
| `kegg_compound_id`, `chebi_id`, `mnxm_id`, `hmdb_id` | str | **Scalars, kept, bare**, holding the canonical entity's value (E2). `chebi_id` tie-break: the charge-0 member, then the lowest number. On `chebi:`-keyed nodes `chebi_id` always equals the id's local part |
| `kegg_compound_ids`, `chebi_ids` | str[] | **new**, same bare format as the scalar, scalar value first |
| `mnxm_ids` | str[] | **new, sparse**: present when the node spans more than one MNX entity |
| `all_identifiers` | str[] | **new, dense, CURIE-only**: one canonical form per identifier, incl. every merged-away or re-keyed id. Forms (prefix lower case, local part normalised): `kegg.compound:C01353` (upper case) · `chebi:14654` (primary and secondary, digits, no leading zeros) · `mnx:MNXM732398` · `hmdb:HMDB0002878` (7-digit form only; 5-digit ids padded at build) · `inchikey:<27-char key>` (upper case). Each entry is on exactly one node across `Metabolite` ∪ `ChemicalRole`, and InChIKey is unique after the B2 scrub (the build fails on a collision). **Unlike `Gene.all_identifiers`**, which is bare. Also in `metaboliteFullText` |
| `recurated_from_ids` | str[] | **new, sparse, no identity meaning**: ids of dropped placeholder or measurement-only nodes whose paper rows now point here (e.g. `kegg.compound:C00525` D-Trp on L-Trp). Never in `all_identifiers` |
| `identity_basis` | str[] | **new, sparse, controlled**: on merged **and re-keyed** nodes (a re-keyed node carries `[chebi_kegg_xref]`); ⊆ `{mnx_same_entity, chebi_kegg_xref, chebi_conjugate, chebi_tautomer, manual}` |
| `structure_kind` | str | **new, dense, controlled**: `defined` \| `generic` \| `polymer` |
| `formula`, `mass`, `inchikey`, `smiles`, `elements` | — | Concrete representative values are dropped from class / polymer nodes that carried them. Generic formulas and `elements` are kept; ≈ 49 formula-less class nodes gain both from ChEBI. `mass` is never present on a generic formula; `elements` ⊆ the periodic table |
| `name` | str | KEGG name on KEGG-keyed nodes, ChEBI name on `chebi:` nodes |
| `kegg_compound_categories` | str[] | **new, sparse, post-import**: br08001 top-level names (`Peptides`, `Carbohydrates`, …). Never inferred |
| `evidence_sources` | str[] | `metabolomics` now only on the nodes paper rows finally resolve to; ≈ 7 nodes lose it |

### Migration notes for the explorer

- **Resolve user ids via CURIE normalisation:**
  1. normalise the input to its CURIE;
  2. try `Metabolite.id`;
  3. then try `IN all_identifiers`;
  4. on a miss, try `recurated_from_ids` (warn "removed; its paper measurements now live on <id>");
  5. on a `chebi:` miss, try **`IN ChemicalRole.all_identifiers`** (warn "this is a chemical
     role"). Old role-metabolite ids are mostly ChEBI **secondary** ids (`chebi:10004` →
     `chebi.role:33229`), so a prefix swap would miss them.
- **Counts shift slightly.** Merges deduplicate, roles leave and measurement-only nodes drop,
  which affects:
  - `Metabolite.organism_count` / `transporter_count`;
  - `Gene` / `OrganismTaxon` `transported_metabolite_count`;
  - `Organism_has_metabolite`;
  - `Reaction_has_metabolite`;
  - `Tcdb_family_transports_metabolite` (−392 role edges, moved to the role edge type, plus
    deduplication on merged targets);
  - `TcdbFamily.metabolite_count`, which loses its roles (e.g. `3.A.1` −33);
  - `measured_*` → 0 on the ≈ 7 nodes that lose the `metabolomics` tag.
- **Merged nodes:** list rollups (`pathway_*`, `evidence_sources`, `organism_names`) are unioned;
  counts are recomputed, never summed.
- **`name` values change** from MNX's lowercase alias to the KEGG / ChEBI name. Goldens and
  search-score order move.
- **Unchanged:** `transport_substrate_resolution` (pinned lumping set) and expression.
- **Assay values are untouched.** Assay edges point at canonical / re-resolved nodes; every move
  is in `metabolite_assay_repoints.tsv`.
- **`list_metabolites`:**
  - the no-formula count drops by ~49;
  - `metabolite_elements` filters only gain matches;
  - `metaboliteFullText` covers name + formula + identifiers.
- **The KEGG pathway set may gain pathway-level `KeggTerm` nodes with no gene**, from re-keyed
  compounds' pathway links. Gene-anchored tools are unaffected.
- **Answers that change:** transport ↔ catalysis joins (e.g. "does MED4 both take up and reduce
  nitrate?") now work on one node, and measured amino acids sit on their L-forms.

## 2. MetaboliteClass: KEGG compound BRITE

### Schema

**Node `MetaboliteClass`.**
- Id: `kegg.brite:<tree_code>.A<i>.B<j>…`. Same positional scheme as `BriteCategory`; the `br08xxx` tree codes cannot collide.
- Properties:

| Property | Type | Notes |
|---|---|---|
| `name` | str | |
| `tree` | str | snake_case tree name |
| `tree_code` | str | e.g. `br08001` |
| `level` | int | 0 = broadest |
| `level_kind` | str | `kegg_compound_class` \| `kegg_compound_subclass` \| `kegg_compound_group` \| `kegg_compound_subgroup` |
| `metabolite_count` | int, post-import | **subtree** count |
| `direct_metabolite_count` | int, post-import | |
| `organism_count` | int, post-import | |

**Edges:**
- `Metabolite_in_metabolite_class` (Metabolite → MetaboliteClass). Carries `sources: ['kegg']` (joins `data_source:kegg`) and `evidence: 'curated'`. A compound listed under several classes has several edges.
- `Metabolite_class_is_a_metabolite_class` (child → parent).

**Indexes:**
- `metabolite_class_tree_idx`, `metabolite_class_level_idx`, `metabolite_class_level_kind_idx`
- full-text `metaboliteClassFullText` on `name`

### Trees and size

| Tree | Covers | Class nodes |
|---|---|---|
| br08001 Compounds with biological roles: Peptides (incl. amino acids, amines), Carbohydrates, Nucleic acids, Vitamins and cofactors, Lipids, Organic acids, Hormones and transmitters, Antibiotics, Steroids | 328 metabolites | 121 |
| br08002 Lipids | 517 | 163 |
| br08003 Phytochemical compounds | 143 | 76 |
| br08021 Glycosides · br08007 Pesticides · br08009 Natural toxins · br08006 Endocrine disruptors · br08005 Bioactive peptides | 95 | 71 |
| **Union** | **867 metabolites, about 27% of the ≈ 3,214 nodes** | **≈ 431** |

≈ 1,130 `Metabolite_in_metabolite_class` and ≈ 420 `Metabolite_class_is_a_metabolite_class` edges.
- Memberships on an alias KEGG id land on its canonical node.
- `tree` is a closed vocabulary of the 8 tree names.
- Ids share the `kegg.brite:` prefix with gene BRITE (see §36, question F).

### Things to know before using it

- **Coverage is partial by design.** KEGG lists representative compounds, not every member of a category. Fumarate, ethanol, carnitine, ribulose-5-P, copper and selenate are unclassified. A metabolite without a class means "KEGG doesn't classify it", **not** "not in this category". Enrichment needs a background of *classified* metabolites.
- **The categories are KEGG's, not textbook chemistry.**
  - Amines, putrescine and dopamine sit under **Peptides** (as do amino acids).
  - myo-inositol is under **Carbohydrates**.
  - Retinol is under **Vitamins**.
  - Tool docs should say so.
- **There are no inferred classes.** A ChEBI-derived mapping was measured against BRITE and only 80% of compounds shared at least one category, with heavy over-assignment, so it was left out.
- **It's separate from gene BRITE.** `ontology='brite'`, `brite_tree` filter values and `genes_by_ontology` are unchanged and won't see these trees.

## 2a. ChemicalRole: role-typed TCDB substrates

- **Meaning:** a TCDB entry like "drug" says **"this family transports something TCDB calls a
  drug"**: an extra substrate given by role, not a role of the family or of its metabolites. The
  edge type carries that meaning.
- **Node `ChemicalRole`:** a pruned ChEBI role hierarchy, ≈ 111 nodes.
  - The 57 TCDB roles plus ancestors.
  - Level 0 = ChEBI's three branches (`biological role`, `chemical role`, `application`); no
    `role` root.
  - Id `chebi.role:<primary ChEBI number>`, a house-minted prefix that never collides with
    `chebi:` metabolites.
  - Properties `name`, `level` (shortest distance, no best-effort marker), `level_kind ∈
    {chebi_role_branch, chebi_role}`.
  - Post-import subtree counts, both counted **distinct** over the DAG, each with a `direct_*`
    version:
    - `transporter_count` counts sources whose role edge is `substrate_depth = 'most_specific'`;
    - `tcdb_family_count` counts sources at any depth.
  - **`all_identifiers`** (dense): its own `chebi.role:` id plus the `chebi:` CURIEs of its primary
    and secondary ChEBI ids, which include every removed role-metabolite id.
- **Edges:**
  - `Chemical_role_is_a_chemical_role` (child → parent);
  - `Tcdb_family_transports_chemical_role` (TcdbFamily → ChemicalRole, ≈ 392). It mirrors the
    metabolite edge: same rollup, same `substrate_depth`.
- **Not propagated:** no Gene / OrganismTaxon / Metabolite count reaches a role.
- **Search and indexes:** `chemicalRoleFullText` on `name`; `chemical_role_id_idx` /
  `_level_idx` / `_level_kind_idx`. In `kg_schema` from day one; not in `list_filter_values`;
  the siderophore route goes through `run_cypher` for now.

## 2b. Also new

- **`Metabolite_is_stereo_form_of`** (stereo-specified → unspecified, directed, ≈ 60) and
  **`Metabolite_is_enantiomer_of`** (symmetric, stored once lower → higher id, query undirected,
  ≈ 25).
  - Both carry `sources: ['chebi']`, `evidence: 'curated'`.
  - Both are emitted only between nodes that survive the measurement-only drop. ≈ 60 / ≈ 25 are
    counts **before** that drop (L-/D-Trp is lost); the final counts come with the post-rebuild
    doc.
  - No anomer pair is an enantiomer edge. No rollup crosses either edge.
- **`stereo_assignment ∈ {as_reported, dominant_form_assumed}`** on `Assay_quantifies_metabolite`
  and `Assay_flags_metabolite`. It's dense on every assay edge and records **what the paper said**.
  - `as_reported`: the paper named the form, or the metabolite has no stereo variant other than
    anomers.
  - `dominant_form_assumed`: an unqualified name (or a stereo-unspecified direct id) landed on one
    member of a stereo family, **whatever chose it**: the rule, the lock, an override, or a single
    exact match. That includes the 3 amino-acid placeholders, every case-fold-converged row, and
    the case where the unspecified node itself wins.
  - Anomer-only families (GlcNAc, MurNAc) keep their target `as_reported`.
- **The dominant-form rule.** A stereo-unqualified paper name, a stereo-unspecified direct id, or
  an exact match onto a stereo-unspecified node goes to the stereo form with the most catalysing
  reactions.
  - It needs a ≥ 2× margin.
  - It's frozen on first computation, and drift is reported.
  - The anomer guard means it never moves onto an α/β anomer.
  - A name **states a form** when it carries a descriptor token anywhere in the name, bounded by
    start / `-` / space / `(` / `,`, with every token scanned. Descriptors:
    - `L-` / `D-` / `DL-`;
    - `(R)-`, `(S)-`, `(R,S)-`, `(R,R)-`, `(2S)-`, `(2R,3S)-`;
    - `(+)-`, `(−)-` (U+2212) and ASCII `(-)-`;
    - `meso-`;
    - `alpha-` / `beta-` / `α-` / `β-` followed by `D-` / `L-`.

    Bare `R-` / `S-` don't count (`S-Adenosyl…` is a sulfur locant). E / Z isomers are separate
    KEGG ids.
  - **Order per row:** paper override → direct id → lock → match → rule. An override beats the
    row's direct id.
  - **Ambiguous names with no override fail the build; unmatched names are unresolved** (no edge,
    reported). Today 1 unmatched name, 6 rows.
  - Today **7 paper names move**, all onto L-forms: Asn, Ala, Glu, Gln, Trp, Met, Phe. 2
    stereo-unspecified direct ids qualify; 0 of them move.
- **Case-insensitive paper-name matching**, except names that parse as formulas (`Co` ≠ `CO`).
- **`TcdbFamily.substrate_breadth ∈ {lumping, specific, not_assessed}`** (dense): the pinned
  lumping set, today's 13 (`2.A.1`, `2.A.1.1`, `2.A.1.2`, `2.A.1.3`, `2.A.3`, `2.A.6`, `2.A.6.2`,
  `2.A.7`, `2.A.66`, `2.A.66.1`, `3.A.1`, `3.A.1.1`, `3.A.3`).
  - `not_assessed` applies at level < 2.
  - `transport_substrate_resolution` reads `substrate_breadth = 'lumping'`.
  - **The threshold query `level >= 2 AND metabolite_count >= 50` no longer reproduces the set**
    (`3.A.3`: 47 molecules, still `lumping`), so filter on `substrate_breadth`.
  - Frozen ids never flip silently; a newly kept family is classified by the threshold and
    appended (`lumping_set_added`).
  - A documented R3 exception, sunset with the deferred TCDB spec.
- **`metaboliteFullText` = `name` + `formula` + `all_identifiers`.** The live index is `name`
  only, so `formula` is new too. A KG test checks that bare-id queries (`C01353`, `MNXM732398`,
  `HMDB0002878`) still hit. If the analyser doesn't produce bare tokens, a search-only property of
  bare local parts is added to the index. New `metabolite_structure_kind_idx` and
  `metabolite_class_id_idx`.
- **New DataSource nodes** `data_source:kegg` and `data_source:chebi` (ChEBI CC BY 4.0
  attribution).
- **Controlled vocabularies:**
  - `Metabolite.identity_basis` / `.structure_kind` / `.kegg_compound_categories`;
  - `MetaboliteClass.level_kind` / `.tree`;
  - `ChemicalRole.level_kind`;
  - `TcdbFamily.substrate_breadth`;
  - `Tcdb_family_transports_chemical_role.substrate_depth`;
  - `sources` / `evidence` on the three new curated edge types;
  - `stereo_assignment` on both assay edges.

  The hash changes.
- **Two release TSVs** (with the post-rebuild doc):
  - `metabolite_id_remap.tsv`: old id, new id, basis ∈ merge basis | `rekey_chebi_kegg_xref` |
    `removed_role` (empty new id; a `role_id` column gives the `chebi.role:` id) |
    `removed_placeholder` / `remeasured_elsewhere` (new id = the node the rows now point to) |
    `chebi_primary_id` | fallback bases for ids none of those explain: the old id is now an
    identifier of another node (`rekey_chebi_kegg_xref` / `chebi_primary_id`; approximate —
    conjugate / zwitterion relocations land here too), the old entry's MNX id now sits on another
    node (`isotopologue_rejected` | `class_term_rejected` | `mnx_same_entity`), the old id's TCDB
    substrate string now resolves to another node (`isotopologue_rejected` | `class_term_rejected` |
    `chebi_conjugate` | `chebi_primary_id`, read from the old ChEBI term), or
    `no_longer_reachable` (empty new id; 0 rows in this build). Written by `scripts/metabolite_identity_tsvs.py`
    from the identity report (whose `baseline` is the previous `kegg_data.json`).
  - `metabolite_assay_repoints.tsv`: paper, table, paper name, old id, new id,
    `stereo_assignment`, reason. The reason is one of `alias_override` |
    `dominant_stereo_form` | `case_fold_convergence` | `direct_id_canonicalised` |
    `placeholder_recurated` | `remeasured_elsewhere` | `merged` (only when stage E merged, or
    stage D re-keyed, the old id into the new node) | `name_re_resolved` (the row's name now
    resolves elsewhere; the old node may still exist); the first applicable wins.
- **The post-rebuild doc also reports:**
  - the per-name stereo table: the 36 B5 names plus the case-fold convergences, with old id,
    new id, `stereo_assignment`, and `catalyst_gene_count` / `transporter_gene_count` on the new
    node;
  - the expected metabolism deltas, so `/omics-edge-snapshot` drops aren't read as regressions;
  - the `structure_kind` counts (and how many carry transport edges);
  - the final stereo / enantiomer edge counts;
  - the number of 0-gene pathway `KeggTerm` nodes added.

---

# Appendix: review history (2026-10-13 / 14)

The sections below are the original proposal's open questions and the full KG ↔ explorer
exchange, in order. Where they disagree with the top of this doc, **the top is authoritative.**

## 3. Questions for the explorer

1. **Alias resolution.** Is resolving through `all_identifiers` enough? Or do you want an explicit `retired_ids` list, separate from current cross-references?
2. **Tool surface for classes.** Is `kegg_compound_categories` on `list_metabolites` rows (plus a filter) enough for v1? Or should a `metabolites_by_class` / `MetaboliteClass` ontology mode in `search_ontology` / `ontology_term_details` come at the same time?
3. **Scope.** Should all 8 trees ship, or only br08001 + br08002? Narrowing loses 107 metabolites, mostly xenobiotic transport substrates.
4. **`structure_kind = generic`.** Should `list_metabolites` hide generic nodes ("amino acid", "Peptide") by default, or only label them?
5. **Breaking-change timing.** Does the id remap need a release boundary or an `mcp_min_version` bump?

## 4. Explorer response (2026-10-13)

**Verdict: go.** This fixes backlog D1/D2 at the source. There are three things to settle before
implementation (E1–E3) and some smaller asks. Code references are to `multiomics_explorer` main
at `81f2d3e`.

### Blocking: settle before implementation

- **E1. Retired *prefixed* ids are not resolved by the explorer today.**
  - `_is_metabolite_alias` (`api/functions.py`) passes any `prefix:` input through verbatim. The only exception is upper-case `CHEBI:NNN`.
  - So after this change, `chebi:14654` and `kegg.compound:C01353` would land in `not_found` with no hint.
  - Explorer fix: if a prefixed id misses on `Metabolite.id`, look it up in `all_identifiers`. Results are reported in `resolved_aliases`, as for bare ids.
  - For this the KG needs to **specify the exact string form of each `all_identifiers` entry**, and to **include every retired node id in its full prefixed form**. That means the 92 re-keyed ids, the ~79 merged-away ids and the secondary `chebi:` ids.
  - Preferred form: the node-id form (`kegg.compound:C01353`, `chebi:14654`, `mnx:MNXM…`) plus the bare/xref forms the resolver already accepts (`C01353`, `CHEBI:14654`, `HMDB…`, `MNXM…`). Then one membership test covers every input shape.
  - No exact-match index is needed: ~3.2k nodes is a trivial scan. `metaboliteFullText` is the wrong tool for this, because Lucene tokenises on `:`.
- **E2. Fate of the scalar xref properties.**
  - The resolver (`build_resolve_metabolite_aliases`) and the `list_metabolites` xref filters use scalar equality: `m.kegg_compound_id IN $…`, `m.chebi_id`, `m.hmdb_id`, `m.mnxm_id`.
  - Please confirm that `kegg_compound_id`, `hmdb_id` and `mnxm_id` stay as scalars holding the canonical entity's value.
  - Please state the element format of the new arrays. `chebi_ids` should match `chebi_id`: raw numeric, no `CHEBI:`.
  - The explorer will then switch the filters to array membership.
  - Note: `list_metabolites` already has params named `kegg_compound_ids` / `chebi_ids` / `mnxm_ids`. Same name as the new properties, same meaning, so that's fine.
- **E3. Assay edges that collapse onto one node.**
  - "Assay edges re-point to canonical nodes" leaves one case open: an assay that measured two ids which now merge.
  - Re-pointing could give the assay **two `Assay_quantifies_metabolite` edges to the same metabolite**. That shows up as duplicate rows in `metabolites_by_quantifies_assay`, `metabolites_by_flags_assay` and `assays_by_metabolite`, and inflates `total_metabolite_count` and `measured_assay_count`.
  - Please state the rule (keep both edges with a `source_metabolite_id`, or refuse to merge measured pairs) and how many assay × metabolite pairs it affects.

### Answers to §3

1. **Alias resolution:** `all_identifiers` is enough, with the conditions in E1. A separate `retired_ids` isn't needed: the explorer reports the input → canonical mapping in `resolved_aliases` whatever the cause.
   - One ask: please ship an **old-id → new-id mapping file** (TSV) with the post-rebuild doc. The explorer has hard-coded `chebi:` / `kegg.compound:` ids in ~12 files (goldens, examples, tool docs, `docs/backlog.md`), and the absorb step needs the mapping to update them.
2. **Class tool surface, v1:** yes to `kegg_compound_categories` on `list_metabolites` rows, plus a filter and a `by_kegg_compound_category` envelope.
   - A `MetaboliteClass` mode in `search_ontology` / `ontology_term_details` is a later change. `ONTOLOGY_CONFIG` assumes gene-anchored terms (`gene_rel`, gene counts, trust axes), so metabolite-anchored terms need their own design.
   - Until then, the other 7 trees can be reached through `run_cypher` only, and the docs will say so.
3. **Scope:** ship all 8 trees. The KG-side cost is small, and the 107 metabolites you would lose are exactly the TCDB xenobiotic substrates that transport questions reach.
4. **`structure_kind = generic`:** label it, don't hide it.
   - Hiding by default would silently drop TCDB substrate rows ("amino acid", "peptide" are common transport substrates). That contradicts the "never default-filter real data" convention.
   - The explorer will add `structure_kind` per row, a `by_structure_kind` envelope breakdown and an opt-in filter.
   - Related: generic nodes losing `elements` means they drop out of `metabolite_elements=['N']` on `metabolites_by_gene` and out of `list_metabolites(elements=…)`. That is correct, but it's a visible behaviour change, and the explorer will document it.
   - Please give counts for `generic` / `polymer` nodes, and how many of them have transport edges.
5. **Breaking-change timing:** bump `mcp_min_version`, and ship the change as an official release that stamps `breaking_changes`. (The alpha.7 cut stamped none; that is explorer backlog 1.2.) Proposed order:
   1. The explorer ships the E1 resolver change first. It reads `coalesce(m.all_identifiers, [])`, so it is a no-op against the current KG.
   2. The KG rebuild sets `mcp_min_version` to that explorer version.
   3. Explorer absorb.

### Smaller asks / notes

- **Controlled vocabularies:** please add `ControlledVocabulary` nodes for `Metabolite.structure_kind`, `Metabolite.identity_basis`, `Metabolite.kegg_compound_categories`, `MetaboliteClass.level_kind` and `TcdbFamily.substrate_roles`, following the vocabulary contract. The hash will change, and the explorer will re-pin during absorb.
- **Stereo split leaves part of D1 open.**
  - The D1 friction also included generic "Cysteine" C00736 (TCDB) vs "L-Cysteine" C00097 (reactions). That was 20 of 49 zero-enzyme cells.
  - Keeping stereoisomers apart is correct, but the transport ↔ catalysis join still fails for these.
  - Not blocking. Could the review report (or a later change) add a non-merging link such as `stereo_parent_id` or a `Metabolite_is_stereo_form_of` edge? Then the explorer could at least warn "0 catalysts; L-form C00097 has N".
- **Numbers to reconcile:** §1 says 92 `chebi:` nodes are re-keyed, but the migration notes say "141 of 635 become KEGG ids". If the gap is the merged-away nodes (92 + 49?), please say so.
- **Rollups after a merge:** please confirm that `pathway_ids` / `pathway_names` / `evidence_sources` / `organism_names` are unioned, and that `catalyst_gene_count` / `transporter_gene_count` / `measured_*` are **recomputed** on the merged node, not summed.
- **`name` change (MNX lowercase → KEGG/ChEBI name):** fine. It moves `list_metabolites` goldens and search-score order. The explorer will absorb this.
- **`substrate_roles`:** no explorer surface planned for now. It shows up in `ontology_term_details(verbose=True)` properties without any change, which is enough.

## 5. KG reply (2026-10-14)

All three blocking items are settled in the spec (§3.4, §3.11). Answering the review surfaced
one more bug, which is now in scope. Two additions came out of the review.

### E1: retired prefixed ids

Agreed as proposed. `all_identifiers` is dense on every `Metabolite`. It holds every
identifier of every merged member, **including every retired node id**, in every form the
resolver accepts:

| Identifier | Forms in `all_identifiers` |
|---|---|
| KEGG compound | `kegg.compound:C01353`, `C01353` |
| ChEBI (primary and secondary) | `chebi:14654`, `CHEBI:14654` |
| MNX | `mnx:MNXM732398`, `MNXM732398` |
| HMDB | `HMDB0002878` |
| InChIKey | full key |

- A bare ChEBI number is **not** included, because it collides with other namespaces. Bare
  numbers resolve through `chebi_ids`.
- A KG test asserts that every form appears on exactly one node.
- We will ship `metabolite_id_remap.tsv` (old id → new id + basis) with the post-rebuild
  version of this doc.

### E2: scalar cross-reference properties

Confirmed. `kegg_compound_id`, `chebi_id`, `mnxm_id` and `hmdb_id` **stay scalars** holding the
canonical entity's value, in today's formats. The new arrays use the same element format as
their scalar:
- `kegg_compound_ids` are bare `C…`;
- `chebi_ids` are raw numeric with no `CHEBI:`, `chebi_id` first;
- `mnxm_ids` are `MNXM…`.

### E3: assay edges onto one node

- **Merges cause 0 collisions.** None of the 69 merge groups contains two measured
  metabolites.
- **Guard for future data:** step 6 refuses a merge that would put two distinct rows of one
  assay table on one node. Step 7 fails loudly when two distinct names of one table resolve to
  one id without an explicit override.
- **But the check found an existing collision from name resolution:** reduced
  **"Glutathione" (GSH, `C00051`) resolves to GSSG (`C00127`)** in Biller 2022 and Kujawinski
  2023.
  - So 6 (assay, metabolite, condition, time point) slots in the current graph carry two edges
    with different values.
  - It's the only distinct-name collision across all paper tables. Every other multi-edge
    pair is a legitimate per-condition or per-time-point edge.
  - It's fixed by the new item below.

### New in scope: ambiguous paper-metabolite names (bug B5)

36 of 149 paper metabolites were resolved as `ambiguous_multi_id`, meaning the code silently
took one of several candidates. Many picks are wrong:

| Paper name | Current id | Problem |
|---|---|---|
| Alanine | `C00133` | D-Alanine |
| Glutamic acid | `C00217` | D-Glutamate |
| Tryptophan | `C00525` | D-Tryptophan |
| Glutamine | `C00819` | D-Glutamine |
| Cystine | `C00954` | indole-3-acetate (shared with "indole 3-acetic acid") |
| Glutathione | `C00127` | GSSG |
| Arginine, cysteine, serine | `mnx:` placeholder nodes | not linked to reactions |

Fix:
- Ambiguous matches are never applied. All 36 are curated as `alias_override` entries, each
  with a reason.
- Stereo-unqualified amino acids go to KEGG's stereo-unspecified entry, because LC-MS doesn't
  resolve stereo, and the new stereo link below connects them to the L-form.
- A new uncurated ambiguous name fails step 7.

**What this means for the explorer:**
- Several measured metabolites change node.
- `assays_by_metabolite` / `metabolites_by_*_assay` goldens that touch amino acids or
  glutathione will move.
- The measured ↔ catalysis join starts working for amino acids.

### New in scope: `Metabolite_is_stereo_form_of` (your D1 ask)

- A curated, **non-merging** edge, stereo-specified → unspecified/broader form (L-Cysteine
  `C00097` → Cysteine `C00736`; D-Cysteine `C00793` → Cysteine `C00736`).
- Carries `sources: ['chebi']` and `evidence: 'curated'`.
- Rule: a direct ChEBI `is_a` between the two terms **and** the same InChIKey connectivity block.
  Both are required.
- **60 edges**, of which **38 connect a transport-only node to a catalysis node**: your "0
  catalysts; L-form has N" case.
- No rollup crosses it. Counts stay per node, and the hop is explicit.

### Smaller items

- **92 vs 141 reconciled.** 141 `chebi:` nodes have a ChEBI→KEGG cross-reference: 49 merge
  into an existing KEGG node, 92 are re-keyed. The full group computation refines this to
  **74 nodes merged away + 89 re-keyed + 57 removed → ≈ 3,225 nodes**. The guard rejected 4
  ChEBI-only KEGG↔KEGG links, including water `C00001` ↔ hydroxide `C01328`, which must stay
  apart.
- **Rollups after a merge:** confirmed. `pathway_ids` / `pathway_names` / `evidence_sources`
  / `organism_names` are unioned. `catalyst_gene_count` / `transporter_gene_count` /
  `transporter_count` / `measured_*` are recomputed from the merged node's edges, never
  summed.
- **`generic` / `polymer` counts.** A draft rule gives ≈ 736 generic-or-polymer nodes, 382 of
  them with transport edges. The final rule (separating `polymer` from `generic`) is settled
  in the plan, and measured numbers come with the post-rebuild doc.
- **Controlled vocabularies.** All five you listed, plus `sources` / `evidence` on
  `Metabolite_in_metabolite_class` and `Metabolite_is_stereo_form_of`.
- **Release: held.** This change ships inside a larger release with other items, so it won't
  be tagged on its own. Your proposed order still holds for that release:
  1. Explorer E1 resolver first (a no-op on the current KG).
  2. The KG release bumps `mcp_min_version` to that explorer version and stamps
     `breaking_changes`.
  3. Explorer absorb.

## 6. Explorer follow-up (2026-10-14)

E1, E2 and E3 are closed. The cross-reference formats, the `exactly one node` test, the remap
TSV, the merge guard and the rollup rule all work for the explorer. B5 and the stereo edge are
welcome. There is one real concern (F1) and a few small clarifications.

### F1. Losing `elements` on ≈ 736 generic/polymer nodes would hurt N-import questions

- 736 nodes is ~23% of all metabolites, and 382 of them have transport edges. Generic transport
  substrates such as "amino acid", "peptide" and "polyamine" are exactly what an N-import
  question reaches.
- With no `elements`, they drop out of `metabolites_by_gene(metabolite_elements=['N'])` and
  `list_metabolites(elements=['N'])`. That is the filter the `pro_n_import_coculture` analysis
  depends on.
- **Ask:** treat `elements` differently from the other dropped properties. Keep `elements` on
  generic nodes when composition is known for the class, taking it from ChEBI's generic formula
  with the `R` removed (e.g. amino acid `C2H4NO2R` → `{C, H, N, O}`).
- Drop `formula`, `mass`, `inchikey` and `smiles` as planned. They really are wrong for a class.
- Where no class formula exists (most polymers), leave `elements` absent.
- If you'd rather not, say so. The explorer would then document the gap and add a warning when
  an element filter drops transport rows on generic nodes. But the KG fix is the better one.

### F2. B5 may break the direct measured ↔ catalysis join it is meant to fix

- Measured stereo-unqualified amino acids move to KEGG's *stereo-unspecified* entry (e.g.
  Alanine `C01401`). In KEGG, reactions mostly sit on the L-forms, and the current wrong D-picks
  sometimes have them (D-Ala `C00133`: `alr`, `ddl`).
- So on the node itself, measured Ala/Glu/Gln/Trp may go from *some* catalysts to *none*. The
  join then only works through the `Metabolite_is_stereo_form_of` hop.
- The decision is right ("LC-MS doesn't resolve stereo"). But "the measured ↔ catalysis join
  starts working" is true only through that hop. Please reword it in the post-rebuild doc.
- **Ask:** in the post-rebuild doc, for each of the 36 re-curated names, list the new id with
  `catalyst_gene_count` / `transporter_gene_count` on the new node, and the same counts summed
  over its stereo forms. The explorer will use the list to size the work and to build test
  cases.
- **Explorer side (not a KG ask):** the explorer will make the stereo hop visible in the same
  release. A sparse per-row hint on `genes_by_metabolite` / `assays_by_metabolite` /
  `list_metabolites` along the lines of "0 catalysts; stereo form C00041 has N". Possibly also
  an opt-in `include_stereo_forms`. Shipping this alongside the KG change means amino-acid
  questions never get worse in between.

### Small clarifications

- **Bare secondary ChEBI numbers.** `chebi_ids` holds primary ids only, so a bare *secondary*
  number won't resolve through it. No KG change needed: the explorer will also test
  `'chebi:' + key IN all_identifiers`. Please just confirm that secondary ids are present in
  the `chebi:` form (§5 table suggests yes).
- **`mnx:` placeholder nodes** for arginine / cysteine / serine (B5 table). Are they deleted
  once their assay edges are re-pointed, or kept as orphans? If deleted, please add them to
  `metabolite_id_remap.tsv`.
- **Remap TSV scope.** Please add a second section, or a sibling file, for the B5 *edge*
  re-points: paper, table, paper name, old id, new id. They aren't node-id changes, so the
  node remap won't show them. The explorer needs them to classify the `assays_by_metabolite` /
  `metabolites_by_*_assay` golden diffs during absorb.

### Release

Understood that it's held and bundled. The explorer will land the E1 resolver change (now
including the bare-secondary-ChEBI path) ahead of the bundle. It is a no-op against the
current KG. F2's stereo hint waits for the rebuild, because it needs the edge to exist.

## 7. KG reply to the follow-up (2026-10-14)

### F1: `elements` on generic nodes (accepted; the measured fix goes further)

Measuring the 736 nodes showed that the original rule ("drop structure on generic nodes") was
too broad. The actual bug is narrower:

| Group | Nodes | With transport edges | Today | After |
|---|---|---|---|---|
| Generic **formula** already (`C2H4NO2R`, `*`, `(…)n`) | 422 | 109 | `formula` + `elements` | **unchanged**: a generic formula describes the class correctly |
| Class term with a **concrete representative** structure (alditol = glycerol, chitin = GlcNAc) | 20 | 17 | concrete formula / mass / InChIKey / SMILES + `elements` | concrete values **dropped**; `elements` **kept** (class composition = representative's); formula refilled from ChEBI's generalized formula where present (9) |
| No formula at all | 316 | 275 | nothing | **49 gain** a generic formula + `elements` from ChEBI's `generalized_empirical_formula` (44 of them transport substrates) |

So **no metabolite loses `elements`**, and N-filter recall goes up. The spec now has this as an
acceptance criterion and a KG test: "the set of metabolites with `elements` only grows".

### F2: stereo-unspecified landing

Agreed, and reworded in the spec. On the node itself, a measured stereo-unqualified amino acid
may have fewer direct catalysts than today's wrong D-pick. The measured ↔ catalysis join works
**through `Metabolite_is_stereo_form_of`**. The post-rebuild doc will carry the per-name table
you asked for, one row per re-curated name:
- the old id and the new id;
- `catalyst_gene_count` / `transporter_gene_count` on the new node;
- the same counts summed over its stereo forms.

Shipping your stereo hint in the same release sounds right.

### Small clarifications

- **Secondary ChEBI ids** are present in `all_identifiers` in both `chebi:NNN` and `CHEBI:NNN`
  form. Your `'chebi:' + key` path works.
- **`mnx:` placeholder nodes.** There are 4, not 3. The fourth is "Cytosolic Ribose 5
  Phosphate", a compartment-prefixed name. Their only edges are assay and
  `Organism_has_metabolite`, so after re-curation they are **deleted** and listed in
  `metabolite_id_remap.tsv` with basis `removed_placeholder`.
- **Remap scope.** There will be two files:
  - `metabolite_id_remap.tsv`: node-id changes. Columns: old id, new id, basis ∈ merge basis |
    `rekey_chebi_kegg_xref` | `removed_role` | `removed_placeholder`.
  - `metabolite_assay_repoints.tsv`: B5 edge re-points. Columns: paper, table, paper name,
    old id, new id, override reason.

## 8. Explorer close-out (2026-10-14)

All items are resolved, and the explorer has no further blocking asks. The F1 rework is better
than our ask: since `elements` only grows, `metabolites_by_gene(metabolite_elements=['N'])`
recall rises instead of falling. Before implementation, please fix a few small consistency
points in this doc and the spec:

- **§1 is now stale on generic nodes.** It needs updating in three places:
  - the "Generic nodes" row of the *What changes* table;
  - the `formula, mass, inchikey, smiles, elements` row of the property table ("Now absent on `generic`/`polymer` nodes");
  - the migration note on `with_chebi` / `no_formula` / `metabolite_elements` ("will move").

  All three still describe the dropped rule. They should match §7: only the 20 representative-structure nodes lose concrete values, and `elements` never shrinks. The `list_metabolites` no-formula count *drops* by ~49.
- **`mass` on generic-formula nodes.** For the 422 nodes with an `R` / `*` / `(…)n` formula, please confirm that `mass` is absent. A mass computed from a formula containing `R` is meaningless, and `list_metabolites` filters on mass.
- **`elements` parsing.** Please confirm that `R`, `*` and polymer `n` never appear as `elements` entries. A KG test like "`elements` ⊆ the periodic table" would cover it.
- **Remap TSV rows with no target.** For `removed_role` rows, please leave `new id` **empty** rather than inventing a target, and the explorer will treat them as "removed". For `removed_placeholder` rows, please say whether `new id` is the re-curated target.

The explorer plan is unchanged:
1. The E1 resolver change, with the `'chebi:' + key` path, lands ahead of the bundled release.
2. After the rebuild, the absorb step uses the two TSVs and the per-name stereo table.
3. The stereo hint ships in the same release as the KG change.

## 9. KG reply to the close-out (2026-10-14)

- **§1 consistency.** All three rows now match §7: the "Generic nodes" row, the property-table
  row, and the `with_chebi` / `no_formula` / `metabolite_elements` migration note. The
  no-formula count drops by ~49.
- **`mass` on generic-formula nodes.** Confirmed absent: 0 of the 420 such nodes carry a mass
  in the live graph. It is now an explicit rule, covering ChEBI-filled formulas too, with a KG
  test.
- **`elements` parsing.** Confirmed: the live graph holds 42 distinct symbols, all real
  elements, with no `R`, `*` or `n`. Added as an explicit rule: an allow-list after `chemparse`,
  with `R` treated as a placeholder before tokenising so a generic `Rn` is never read as radon.
  KG test: `elements` ⊆ periodic table.
- **Remap rows with no target.**
  - `removed_role` rows have an **empty** `new id`.
  - `removed_placeholder` rows carry the re-curated target in `new id`: the node that the
    placeholder's assay rows now point to. The same re-point also appears, per paper row, in
    `metabolite_assay_repoints.tsv`.

The spec is updated to match. Next step on the KG side: spec sign-off, then the implementation
plan.

## 10. KG question: prefixed-only `all_identifiers`? (2026-10-14, open, not decided)

This reopens the string form agreed in E1 / §5. Nothing is changed in the spec yet; we would
like the explorer's view first.

**Current agreement (§5).** `all_identifiers` holds both the prefixed and the bare/xref form
(`kegg.compound:C01353` + `C01353`, `chebi:14654` + `CHEBI:14654`, `mnx:MNXM…` + `MNXM…`), so
one membership test matches any raw input.

**Alternative under consideration.** `all_identifiers` holds **one canonical CURIE per
identifier**, in the same convention as node ids:
- `kegg.compound:C01353`
- `chebi:14654` (primary and secondary)
- `mnx:MNXM…`
- `hmdb:HMDB…`
- `inchikey:…`

The explorer normalises user input to that form before the membership test, using the patterns
it already has (`_METABOLITE_ALIAS_PATTERNS`):
- bare `C\d{5}` → `kegg.compound:`;
- `MNXM…` → `mnx:`;
- `CHEBI:NNN` / bare digits → `chebi:`;
- `HMDB…` → `hmdb:`;
- a 27-character InChIKey → `inchikey:`.

**Why we're considering it:**
- One form per identifier. The KG doesn't encode input-normalisation choices into the data.
- The list isn't doubled.
- The upper-case `CHEBI:` duplicate and the bare-number exclusion rule go away.
- The uniqueness test ("each entry on exactly one node") becomes a plain duplicate check.
- It matches the house CURIE convention used for node ids (`normalize_curie`).

**What it costs the explorer:**
- A normalisation step before lookup, instead of matching raw input.
- The E1 resolver change in flight would have to adopt it before the bundled release.

**Not changing either way:**
- The scalar xref properties (`kegg_compound_id`, `chebi_id`, `mnxm_id`, `hmdb_id`) and the new
  arrays (`kegg_compound_ids`, `chebi_ids`, `mnxm_ids`) stay **bare**, as agreed in E2.
  Prefixing them would be a separate vocabulary-contract change with filter migrations and no
  new capability.
- `mnx:` stays house-minted. The registered bioregistry prefix is `metanetx.chemical`, but
  switching would rename node ids.

**Questions for the explorer:**
1. Is normalise-then-match acceptable in the resolver, or is matching raw input against the
   data a deliberate design choice on your side (e.g. shared with `resolve_gene`)?
2. Does any other consumer read `all_identifiers` raw, such as goldens, `run_cypher` examples or
   docs, in a way that relies on bare forms being present?
3. Do you prefer one canonical form, or the current both-forms agreement?

## 11. Explorer answer on prefixed-only `all_identifiers` (2026-10-14)

§9 is accepted as is, and all of §1–§9 is closed.

**§10: go with one canonical CURIE per identifier.** Answers to your three questions:

1. **Normalise-then-match is fine.**
   - Matching raw input was a convenience, not a design choice. The resolver already sorts every input by shape (`_METABOLITE_ALIAS_PATTERNS`), so mapping the shape to a CURIE costs one line per pattern and keeps the rules in one place.
   - Rework cost is zero: the E1 resolver change isn't written yet, so it will adopt this form from the start.
   - It doesn't follow `resolve_gene`, but that's fine. `resolve_gene` matches raw, case-insensitively (`toLower`), against `Gene.all_identifiers`, because locus tags and RefSeq ids are unprefixed and namespace-unambiguous.
   - Metabolites are the reverse case. A bare `14654` or a bare InChIKey only means something with a namespace, so prefixed entries are the right fit.
2. **No consumer reads `Metabolite.all_identifiers` today.** The property doesn't exist yet. The only `all_identifiers` hits in the explorer are on `Gene` (`resolve_gene`, `gene_details` goldens, unit tests), and they are unaffected.
   - **Side effect:** the same property name will follow two conventions. `Gene.all_identifiers` holds bare ids; `Metabolite.all_identifiers` holds CURIEs.
   - That's acceptable, but please say so in the property's schema description so `kg_schema` / `run_cypher` users see it. The explorer will add a line to `docs://guide/conventions`.
3. **Preference: one canonical form.**

**Conditions:**
- **Normalise the local part in the KG too, not just the prefix.** The explorer will upper-case KEGG / HMDB / MNX / InChIKey local parts and lower-case the prefix before matching. Please store them in exactly that form.
  - **HMDB:** store the 7-digit form only (`hmdb:HMDB0002878`). The explorer will pad legacy 5-digit input (`HMDB02878`) before matching.
  - **ChEBI:** no leading zeros.
- **InChIKey becomes a new input shape.** The explorer doesn't accept InChIKeys today and will add the pattern `^[A-Z]{14}-[A-Z]{10}-[A-Z]$` → `inchikey:`. Please state how the KG handles one InChIKey shared by several nodes (e.g. merged tautomers). Is that impossible after the merge, or covered by the "exactly one node" test? If it can happen, the explorer will expand to all matches and warn, as it already does for xref collisions.
- **Full-text search.** `all_identifiers` is indexed in `metaboliteFullText`. Please add a KG test that `list_metabolites(search_text='C01353')`-style Lucene queries (a bare local part) still hit with prefixed-only entries. With the standard analyser, `kegg.compound:C01353` should yield a `c01353` token, but this needs checking rather than assuming.
- **No change to E2:** the scalars and the `*_ids` arrays stay bare, as you wrote.

## 12. KG decision: prefixed-only `all_identifiers` (2026-10-14, decided)

Decided: **one canonical CURIE per identifier.** The spec (§3.4) is updated. Your conditions:

- **Local-part normalisation in the KG.** Accepted exactly as you specified:
  - lower-case prefix;
  - KEGG / MNX / HMDB / InChIKey local parts upper case;
  - HMDB in the 7-digit form only (`hmdb:HMDB0002878`; legacy 5-digit ids are padded at build
    time);
  - ChEBI as digits with no leading zeros;
  - InChIKey stored as `inchikey:<27-character key>`.
- **Shared InChIKey.** It can't happen after the change, and the build fails if it does.
  - After merging, 5 InChIKeys would still sit on two nodes. All are the representative-
    structure bug (B2): a class or polymer carrying its monomer's key (alditol / glycerol,
    (1,4-α-D-glucosyl)n / glucose, chitosan / glucosamine, PHB / 3-hydroxybutanoate, chitin /
    GlcNAc).
  - The B2 fix removes the key from the generic/polymer side, so InChIKey is unique by
    construction.
  - The "each entry on exactly one node" test covers it, and any leftover collision fails the
    build. You don't need the expand-and-warn path for InChIKey; keeping it as a defensive
    fallback is harmless.
- **Full-text with bare local parts.** Added as a KG test: `C01353`, `MNXM732398` and
  `HMDB0002878` queries must hit through `metaboliteFullText`, verified against the live
  analyser. If the standard analyser doesn't produce the bare token, the bare local parts go
  into a separate search-only property, and `all_identifiers` stays CURIE-only.
- **Two conventions under one property name.** Noted. The `Metabolite.all_identifiers` schema
  description will say "CURIEs, unlike `Gene.all_identifiers`".
- **E2 unchanged:** the scalars and the `*_ids` arrays stay bare.

## 13. Explorer sign-off (2026-10-14)

§12 is accepted as written, and the explorer has no open questions. **The explorer signs off on
this spec.**

One note on the full-text fallback in §12: a separate search-only property is fine, as long as
it is in `metaboliteFullText`. `list_metabolites(search_text=…)` uses only that index, so the
explorer needs no change either way.

Explorer work items:

1. **Before the bundled release:** the E1 resolver.
   - Turn each input shape into its CURIE: KEGG, ChEBI (`CHEBI:` / bare digits, leading zeros stripped), MNX, HMDB (5→7-digit pad), InChIKey (new).
   - Check `Metabolite.id` first, then `'<curie>' IN coalesce(m.all_identifiers, [])`.
   - Report the mapping in `resolved_aliases`.
   - Keep the expand-and-warn path as a defensive fallback.
   - It is a no-op against the current KG.
2. **Same release as the KG change:** `structure_kind` (per row, `by_structure_kind`, an opt-in filter) and `kegg_compound_categories` (per row, a filter, `by_kegg_compound_category`) on `list_metabolites`, plus the stereo-form hint (F2).
3. **Absorb:**
   - update goldens and examples from `metabolite_id_remap.tsv` and `metabolite_assay_repoints.tsv`;
   - re-pin the vocabulary hash;
   - switch the `list_metabolites` xref filters to `*_ids` array membership;
   - document the CURIE vs bare `all_identifiers` asymmetry in `docs://guide/conventions`.

## 13a. KG change to F2: dominant stereo form for unqualified paper names (2026-10-14, decided)

This **reverses** the B5 default agreed in §5 / §7 (stereo-unqualified amino acids → KEGG's
stereo-unspecified entry). It simplifies your F2 work.

**Why.** The catalysts sit on the specified forms. The unspecified nodes have almost none:

| Amino acid | L-form catalyst genes | Unspecified node | D-form |
|---|---|---|---|
| Alanine | 520 | 18 | 103 |
| Glutamate | 1,905 | 69 | 101 |
| Methionine | 509 | 0 | — |
| Asparagine | 126 | 0 | 0 |
| Serine | 436 | 0 | 0 |
| Tryptophan | 207 | — | 0 |

Their large transporter counts come only from TCDB's generic substrate listings. The current
measured data is inconsistent too: Lys, Pro, Leu, Ile, Thr and Val already sit on L-forms;
Ala, Glu, Gln and Trp sit on D-forms; Asn, Met and Phe sit on unspecified nodes.

**Rule.** When a paper name **does not specify stereo**, it resolves to the member of its stereo
family with the **most gene-reachable catalysing reactions** in the build. That is the form our
organisms' metabolism uses: L for amino acids, usually D for sugars.
- The family is the same strict one as the stereo link: a ChEBI `is_a` / enantiomer link plus
  the same InChIKey connectivity block.
- A tie, or a family with no catalysis, needs a curated override; nothing is guessed.
- **Never applied** when the paper names a form ("D-Alanine" stays D), or when an override
  exists (e.g. a chiral-method paper).

**Measured over all 149 paper metabolites:**

| Outcome | Count |
|---|---|
| Already on the most-catalysed form | 24 |
| No stereo variant in the KG | 118 |
| **Move** | **7**, all onto the L-form |

The 7 moves: Asparagine → L-Asn, D-Ala → L-Ala, D-Glu → L-Glu, D-Gln → L-Gln, D-Trp → L-Trp,
Methionine → L-Met, Phenylalanine → L-Phe. Arginine / cysteine / serine (the `mnx:`
placeholders) go to L-Arg / L-Cys / L-Ser by override.

**New edge property.** `Assay_quantifies_metabolite` and `Assay_flags_metabolite` gain a dense
**`stereo_assignment ∈ {as_reported, dominant_form_assumed}`**, registered as a controlled
vocabulary. A consumer can tell "the paper said L-Ala" from "the KG assigned L-Ala because the
paper said Alanine".

**What it means for the explorer:**
- **The measured ↔ catalysis join works on the node itself** for amino acids, with no stereo
  hop. The F2 hint ("0 catalysts; stereo form has N") is still useful, but only for TCDB's
  generic substrates (e.g. Cysteine `C00736` with 0 catalysts vs L-Cysteine `C00097` with 652).
  It is no longer needed for measured amino acids.
- The per-name table in the post-rebuild doc carries `stereo_assignment` per row. The summed
  over stereo forms column is dropped, because it is no longer needed.
- Surfacing `stereo_assignment` on `assays_by_metabolite` / `metabolites_by_*_assay` rows is
  your call. We'd suggest it, because it is the only place this assumption is recorded.

**Questions:**
1. Does this change your F2 plan (stereo hint for TCDB generic substrates only)?
2. Do you want a rollup on the `Metabolite` node as well (e.g. a sparse
   `measured_stereo_assignments: str[]` of the distinct values on its assay edges), or is the
   edge enough?

## 14. Explorer answer on the dominant-stereo-form change (2026-10-14)

*(Housekeeping: there are now two sections numbered 13. Please renumber the KG one to 13a, or
renumber this one, whichever you prefer.)*

**Agreed.** Assigning the dominant form is the better default. Unqualified LC-MS "Alanine" in
these cultures is overwhelmingly L. The current mix of L, D and unspecified is worse than
either consistent rule. `stereo_assignment` keeps the assumption visible on every row.
Three conditions follow.

### Conditions

- **C1. Stable across builds.** "Most gene-reachable catalysing reactions *in the build*" depends
  on the data. Adding organisms or KO annotations could flip a sugar's winner between releases,
  and a measured value would then silently change node.
  - Please **freeze** each computed assignment into the curated override table (paper name →
    id, reason `dominant_form`) the first time it is computed.
  - On later builds, compare the computed winner with the frozen one. If they differ, report it
    for review rather than re-pointing silently.
  - Ties and no-catalysis families already go to override; this extends the same rule to
    *changes*.
- **C2. Require a margin, not just "no tie".** All 7 measured moves are clear (e.g. L-Ala 520 vs
  D-Ala 103). But the rule will meet sugars and hydroxy acids where the forms are close.
  - Please require a minimum dominance before auto-assigning, e.g. the winner has ≥ 2× the
    runner-up. Anything under the threshold goes to override.
  - D-forms are real biology here: peptidoglycan D-Ala and D-Glu. Close counts are where a wrong
    guess would matter.
- **C3. `stereo_assignment` should record what the paper *said*, not how the KG got there.** The
  3 `mnx:` placeholders (Arginine / Cysteine / Serine) go to L-forms by override, but the paper
  still didn't name a form.
  - Please confirm that those edges carry `dominant_form_assumed`, not `as_reported`.
  - Likewise, any future override for an unqualified name should be `dominant_form_assumed`.
    `as_reported` is only for a paper that named the form, or for an achiral metabolite.

### Answers

1. **F2 plan:** yes, it narrows.
   - The stereo hint is now only for the transport-side generic substrates (e.g. Cysteine
     `C00736` vs L-Cysteine `C00097`).
   - It lives on `genes_by_metabolite`, `metabolites_by_gene` and `list_metabolites`. The assay
     tools don't get it.
   - It still needs `Metabolite_is_stereo_form_of`, so it still ships in the same release as the
     KG change.
2. **Node rollup: the edge is enough.** Please don't add `measured_stereo_assignments`. The
   explorer will surface the property from the edge:
   - **`stereo_assignment` as a per-row field** in compact mode on `metabolites_by_quantifies_assay`,
     `metabolites_by_flags_assay` and `assays_by_metabolite`. It's short, and it qualifies the
     identity of the row, so it belongs in compact rows rather than only in verbose output.
   - **A `by_stereo_assignment` breakdown** in the envelope of all three tools.
   - **A warning** when the result contains any `dominant_form_assumed` rows: "N rows assume
     the dominant stereo form; the paper did not specify".

   The node rollup would only matter for a "which measured metabolites are stereo-assumed"
   question at the `list_metabolites` level. The explorer can answer that from the assay tools
   without it. We'll revisit if that question comes up.

## 15. KG reply on C1–C3 (2026-10-14)

Renumbered: the KG dominant-stereo section is now **13a**. All three conditions are adopted as
written, in spec §3.10:

- **C1 (stable across builds).** The first assignment is frozen into the curated override table
  (`reason: dominant_form`, `frozen_on`, the winner / runner-up counts). Later builds recompute
  and compare. A changed winner, or a margin that drops below threshold, is reported as
  `dominant_form_drift` in the identity report and never re-points silently.
- **C2 (margin).** Auto-assignment requires the winner to have ≥ 2× the runner-up's
  catalysing-reaction count. Anything closer, a tie, or no catalysis goes to override. The 7
  current moves all clear it (the lowest is L-Ala 520 vs D-Ala 103, ≈ 5×).
- **C3 (what the paper said).** Confirmed: `dominant_form_assumed` covers *every* unqualified
  name resolved to one member of a stereo family, whether by the rule or by an override,
  including the arginine / cysteine / serine `mnx:` placeholders. `as_reported` is only for a
  name that states the form, or a metabolite with no stereo variant in the KG. Overrides carry an
  explicit `stereo_assignment` field, and step 7 fails when one is missing.

The node rollup (`measured_stereo_assignments`) is dropped, as you asked. Your plans for the
per-row field, the envelope breakdown and the warning are noted.

## 16. KG note: `substrate_roles` contract + one side effect (2026-10-14, FYI)

The full contract is in spec §3.5. The parts that touch the explorer:

- **`TcdbFamily.substrate_roles: str[]`.**
  - Sparse (absent, never `[]`).
  - Values are ChEBI role names, verbatim (e.g. `drug`, `siderophore`, `xenobiotic`,
    `vitamin (role)`).
  - Rolled up exactly like substrates, so class-level nodes carry long lists (33 on `3.A.1`).
  - Registered as an **open** controlled vocabulary.
  - Today: 392 (family, role) pairs on 116 families, from 57 role terms that stop being
    `Metabolite` nodes.
- **`tcdbFamilyFullText` gains `substrate_roles`** (drop + recreate), so `search_ontology` over
  TCDB finds "siderophore" / "drug" families.
- **Side effect on a stored gene property.** P-type ATPase `3.A.3` drops from
  `metabolite_count` 50 to 49, one of its 50 being a role. It therefore leaves the
  substrate-lumping set (`level ≥ 2 AND metabolite_count ≥ 50`).
  - 66 genes have `3.A.3` as their deepest TCDB attachment, and some of them move
    `transport_substrate_resolution` from `family_inferred` to `resolved`.
  - Goldens touching P-type ATPase genes (Mg²⁺ / Cu⁺ / Ca²⁺ / K⁺ pumps) may move.
  - Docs that list the 7 lumping families should drop `3.A.3`; the list becomes 6 families + 6
    subfamilies.
  - The exact gene count comes with the post-rebuild doc.

## 17. Explorer reply on §15 / §16 (2026-10-14)

**§15:** all three conditions are adopted as written, so the stereo work is closed on our side.

**§16 `substrate_roles`:** the contract is fine. Having it sparse, verbatim, an open vocabulary and
in `tcdbFamilyFullText` works for us. `search_ontology(ontology='tcdb', search_text='siderophore')`
picks it up with no explorer change.

### One objection: don't let role removal change the lumping set (`3.A.3`)

- `3.A.3` leaving the lumping set is a **threshold artefact**, not new information. The family
  lists the same substrates as before. One of its 50 is now stored as a role instead of a
  metabolite, and that alone moves it across `metabolite_count ≥ 50`.
- The result reads the wrong way. Genes whose deepest TCDB attachment is the P-type ATPase
  superfamily would become `resolved`, which tells a user "drill into the substrates". But their
  breadth is still ~49 cations and other substrates: that is reachability, not capability. This is
  exactly the over-claim that `family_inferred` exists to prevent.
- **Ask:** measure substrate breadth for the lumping rule as
  `metabolite_count + size(coalesce(substrate_roles, []))`.
  - Roles are still TCDB substrate listings, so counting them is semantically consistent.
  - It keeps the lumping set the same as today: `3.A.3` stays, and no gene's
    `transport_substrate_resolution` changes from this change.
  - If you'd rather keep `metabolite_count` alone, please pin the lumping set explicitly
    (a curated list plus a report when the computed set differs). Then a membership change is a
    decision, not a side effect. Either way, please add a KG test that this change doesn't alter
    the lumping set.
- This also means the docs need no "6 + 6" edit. In fact the explorer docs never list the
  lumping families: `tools.py`, `metabolites_by_gene.yaml` and `tcdb.yaml` describe the rule, not
  the members.

### Absorb notes (no KG action)

- Removing roles shrinks `3.A.1`'s substrate list (33 roles). So the empirical figures in
  `docs://analysis/metabolites` will move: 13 MED4 `family_inferred` genes, 554 transport rows
  each, 312 inherited / 242 `most_specific`. The explorer will re-measure them during absorb.

## 18. KG FYI: spec restructured as a data-flow spec (2026-10-14)

The spec was rewritten around inputs → prepare-data → files → schema → adapters → post-import.
That surfaced a few gaps, which are now filled. Three of them are visible to the explorer; none
is blocking:

1. **Paper-supplied ids are now canonicalised.** Step 7 used to copy `kegg_direct` /
   `chebi_direct` ids from paper tables verbatim, bypassing identity. They now go through a new
   `id_aliases` map in `metabolite_id_mapping.json`.
   - An assay row whose paper id gets merged away lands on the canonical node instead of
     dangling.
   - These re-points appear in `metabolite_assay_repoints.tsv` with reason
     `direct_id_canonicalised`.
   - Their `stereo_assignment` is `as_reported`.
2. **`metaboliteFullText` also gains `formula`.** The live index covers `name` only; CLAUDE.md
   wrongly lists `formula`. Recreated, it covers `name`, `formula` and `all_identifiers`.
   - Formula-like `list_metabolites(search_text=…)` queries can now hit, and score order can
     shift.
   - Search goldens may move.
3. **New entries in listings:**
   - DataSource nodes `data_source:kegg` and `data_source:chebi`, which back the `sources` values
     on the class-membership and stereo edges (R2);
   - a closed vocabulary entry `MetaboliteClass.tree` (the 8 tree names).
   - Both show up in data-source / vocabulary listings and are covered by the hash re-pin you
     already plan.

Internal only, no explorer impact:
- **C1 storage.** Frozen dominant-form picks are kept in a committed lock file
  (`cache/data/metabolomics/dominant_stereo_forms.json`) instead of the curated paper alias
  files, so the build never writes into curated inputs. Behaviour is as agreed: frozen on first
  computation, drift reported, never re-pointed silently.
- **Alias format.** The per-paper `metabolite_aliases.yaml` format is extended to
  `{id, stereo_assignment, reason}`.
- **Step 7.** Its own `ambiguous_multi_id` branch is removed too, with the same outcome as the
  step-6 fix.
- **Placeholder.** The fourth `mnx:` placeholder ("Cytosolic Ribose 5 Phosphate") goes to
  D-ribose 5-phosphate. The node estimate is now ≈ 3,221.

## 19. KG reply on the `3.A.3` objection (2026-10-14, decided)

Agreed that the flip is an artefact. Your first fix (counting roles in the breadth) **isn't
enough on its own**, though. The merges also deduplicate `3.A.3`'s list (two of its substrates
collapse into one node), so it reaches 48 even with the role counted. It is the only family in
today's 13-node lumping set that moves.

**Decision: your fallback, a pinned set** (spec §3.4a):
- `tcdb_pruned.json` carries a committed `lumping_tcdb_ids`, frozen as today's 13: `2.A.1`,
  `2.A.1.1`, `2.A.1.2`, `2.A.1.3`, `2.A.3`, `2.A.6`, `2.A.6.2`, `2.A.7`, `2.A.66`, `2.A.66.1`,
  `3.A.1`, `3.A.1.1`, `3.A.3`.
- Each build recomputes the threshold set and reports `lumping_set_drift`; it never applies a
  change silently. Changing the pin is a human edit.
- New dense **`TcdbFamily.substrate_breadth ∈ {lumping, specific}`**, a closed controlled
  vocabulary.
- `transport_substrate_resolution` now reads `substrate_breadth = 'lumping'` instead of the
  inline `level ≥ 2 AND metabolite_count ≥ 50`.
- Recorded as an R3 exception. A pinned, curated set isn't a restatement of
  `metabolite_count`; after this change they can disagree (`3.A.3`: 47 molecules, still
  `lumping`).
- **No gene's `transport_substrate_resolution` changes from this spec.** KG tests: the
  `lumping` nodes equal the pin, and the resolution distribution equals the pre-change baseline.
- No docs edit about "6 + 6": noted that your docs describe the rule, not the members. If you
  want, they can now point at `substrate_breadth` instead of restating the threshold.
- Your absorb note (the `3.A.1` figures in `docs://analysis/metabolites` move because of role
  removal) still applies.

## 20. Explorer reply on §18 / §19 (2026-10-14)

**§19: accepted.** A pinned `lumping_tcdb_ids`, `substrate_breadth`, drift reporting and the two
KG tests are better than our first fix.
- The explorer will point its docs at `substrate_breadth` instead of restating the threshold:
  - `transport_substrate_resolution` field descriptions in `tools.py`;
  - `metabolites_by_gene` / `genes_by_metabolite` YAML;
  - `inputs/ontologies/tcdb.yaml`.
- `substrate_breadth` surfaces in `ontology_term_details(verbose=True)` properties without
  change. Adding it to TCDB rows elsewhere is a later, separate decision.

**§18: items 2 and 3 are fine.** Notes for absorb:
- **Full-text `formula` (item 2).** You're right that the live index is `name` only. The
  explorer's own `list_metabolites` docs say so explicitly ("`search_text` is a Lucene search
  over the metabolite name only — NOT the formula"), and that line becomes wrong once the index
  is recreated. The explorer will rewrite it during absorb: the index covers name, formula and
  identifiers. `elements` is still the right filter for composition.
  - Also, an old explorer design spec repeats the `name + formula` claim. It's a historical
    document and won't be edited.
- **Item 3:** the new DataSource nodes and `MetaboliteClass.tree` are covered by the planned
  hash re-pin.

### One question on item 1: paper-supplied *stereo-unspecified* ids

Canonicalising `kegg_direct` / `chebi_direct` ids is right, and so is `as_reported` when the
paper's id names a specific form. But there is a gap:
- Suppose a paper table supplies a **stereo-unspecified** id directly, e.g. Alanine `C01401`.
- That id doesn't assert a form; it records that the paper didn't specify one.
- Under §18 it stays on `C01401` (18 catalysts) as `as_reported`. The same compound given only
  by *name* goes to L-Ala (520) as `dominant_form_assumed`.
- So which node a measurement lands on would depend on whether the paper's table had an id
  column, not on the biology. That is the inconsistency 13a was meant to remove.

**Ask:** when a paper-direct id is the stereo-unspecified member of a stereo family (it is the
target of `Metabolite_is_stereo_form_of` edges), treat it like an unqualified name:
- apply the dominant-form rule, with the C1 freeze and the C2 margin;
- set `stereo_assignment = 'dominant_form_assumed'`;
- list it in `metabolite_assay_repoints.tsv`.

Please report how many paper-direct ids this affects. It may well be 0 today, in which case a
KG test that would catch the first one is enough.

## 21. KG reply on §20: paper-supplied stereo-unspecified ids (2026-10-14, decided)

**Adopted** (spec §3.6 / §5.3 / §6.5):
- A paper-direct id that is the stereo-unspecified member of a family goes through the same
  dominant-form rule as an unqualified name, with the C1 freeze and the C2 margin.
- When it moves: `dominant_form_assumed`, listed in `metabolite_assay_repoints.tsv`.
- Step 6 pre-computes the decision (`direct_id_resolution` in `metabolite_id_mapping.json`), so
  step 7 stays a lookup.

**Measured: 2 of 64 paper-direct ids** are stereo-unspecified family targets:
- N-Acetyl-D-glucosamine `C00140`
- N-Acetylmuramate `C02713`

**0 move.** Both already carry the catalysts (65 vs 0 for β-GlcNAc; 98 vs 0 for aldehydo-MurNAc).

**New guard that this surfaced: anomers.** Both cases differ from their variants only at the
anomeric carbon (α/β). Anomers interconvert in solution and LC-MS doesn't separate them, so the
rule must **never** move a measurement onto an anomer. Glucose shows the risk: D-glucose has 284
catalysts, α-D-glucose 183, β-D-glucose 192; today only the 2× margin happens to refuse a move.
- The rule now excludes anomer-only differences as destinations. Anomer stereo edges are still
  emitted.
- So GlcNAc / MurNAc stay `as_reported`: the paper's id named the D-form, and the anomer isn't a
  stereo choice.
- A KG test catches the first stereo-unspecified direct id that would move.

**Also decided on the KG side, FYI:** ChEBI *roles* on `Metabolite` (a `has_role` facet: cofactor,
siderophore, osmolyte, …) are **out of this spec**.
- Measured: after removing organism-observation roles ("*E. coli* metabolite") and pharmacology,
  only a few dozen compounds carry ecologically useful roles.
- A clean version needs a curated allow-list.
- It's filed in `plans/backlog.md`. It would use the same vocabulary as
  `TcdbFamily.substrate_roles`.

## 22. KG change: role-typed substrates become edges, not `substrate_roles` (2026-10-14, decided)

This **replaces the §16 `TcdbFamily.substrate_roles` contract**. Spec §3.4.

**Why.** A TCDB role entry ("drug", "siderophore") means **"this family transports something
TCDB calls a drug"**: an extra substrate given by role instead of by molecule. It is not a role
of the family, and not a role of the family's metabolite substrates. A property named
`substrate_roles` reads naturally as the latter. It also had two structural limits:
- it couldn't carry the per-pair `substrate_depth` the metabolite edges have;
- it would collide with a future molecule-level `has_role` that uses the same vocabulary but
  makes a different claim.

**New model:**
- **Node `ChemicalRole`** (57 today).
  - Id: the primary ChEBI CURIE, e.g. `chebi:33229` "vitamin (role)".
  - Properties: `name` (verbatim), and `role_branches: str[]` ⊆ {`biological role`,
    `chemical role`, `application`}, a closed vocabulary. 10 roles sit in more than one branch.
  - Post-import: `transporter_count` (most-specific sources, the `Metabolite.transporter_count`
    rule) and `tcdb_family_count`.
- **Edge `Tcdb_family_transports_chemical_role`** (TcdbFamily → ChemicalRole, 392 today on 116
  families).
  - It mirrors `Tcdb_family_transports_metabolite`: the same rollup and the same
    `substrate_depth ∈ {most_specific, inherited}`.
  - The meaning is in the edge type.
- **Full-text `chemicalRoleFullText`** on `name`, plus `chemical_role_id_idx`.
  `tcdbFamilyFullText` is **unchanged** (the §16 plan to add `substrate_roles` to it is dropped).
- **Not propagated.** No Gene / OrganismTaxon / Metabolite count reaches a `ChemicalRole`, and
  the lumping pin is unaffected (KG test).
- **`ChemicalRole` and `Metabolite` ids never coincide**, although both use `chebi:` (KG test).
- **Future molecule roles** (the backlogged `Metabolite.chebi_roles`) become
  `Metabolite_has_chemical_role` → the **same** `ChemicalRole` nodes, through a different edge
  type. "Has the role" and "transports something with the role" stay distinct but joinable.

**What changes for the explorer vs §16/§17:**
- **"Which TCDB families transport siderophores / drugs" changes route.**
  - Was: `search_ontology(ontology='tcdb', search_text='siderophore')` through
    `tcdbFamilyFullText`.
  - Now: a 1-hop query from the `ChemicalRole` node
    (`(:ChemicalRole {name:'siderophore'})<-[:Tcdb_family_transports_chemical_role]-(t)`),
    found through `chemicalRoleFullText`. Your call whether to surface it as a tool mode now or
    leave it to `run_cypher`.
- **Two removed metabolites get a new home, not a hole.** `metabolite_id_remap.tsv` rows for the
  57 roles keep basis `removed_role`, and now give the `ChemicalRole` id in a new `role_id`
  column; `new id` stays empty.
- **New label, edge type and vocabularies** (`ChemicalRole.role_branches`,
  `Tcdb_family_transports_chemical_role.substrate_depth`), covered by the planned hash re-pin.
  `TcdbFamily.substrate_roles` is **not** registered.

**Question:** does the explorer want `ChemicalRole` in `kg_schema` / `list_filter_values` from
day one, or only reachable through `run_cypher` until a tool needs it?

## 23. KG additions: role hierarchy + enantiomer edge (2026-10-14, decided)

Two refinements to §22 and the stereo model (spec §3.4 / §3.5).

**1. `ChemicalRole` is a pruned ChEBI role hierarchy, not 57 flat nodes.**
- **Contents:** the 57 TCDB roles plus their `is_a` ancestors, ≈ 111 nodes, linked by
  `Chemical_role_is_a_chemical_role` (child → parent).
- **Levels:** the level-0 nodes are ChEBI's three role branches, `biological role`,
  `chemical role` and `application`. The single `role` root isn't emitted. `level` is the
  SHORTEST `is_a` distance to a branch (a DAG, so best-effort like GO): it is not monotone
  along every `is_a` edge — a child can sit at the same or a lower level than one of its
  parents. `level_kind ∈ {chebi_role_branch, chebi_role}`.
- **`role_branches` (§22) is dropped:** a role's branch is its level-0 ancestor, and 10 roles
  have several.
- **Counts:** `transporter_count` / `tcdb_family_count` are **subtree** counts
  (`Chemical_role_is_a_chemical_role*0..`, the house convention), with `direct_*` versions.
- **Why:** "families transporting siderophores" now also catches **phytosiderophore**, and
  "insecticide" catches **neonicotinoid insecticide**; both are among today's 57.

**2. New `Metabolite_is_enantiomer_of`** (≈ 25 edges, `sources: ['chebi']`,
`evidence: 'curated'`).
- **Symmetric, stored once** (lower → higher canonical id). **Query it undirected:**
  `-[:Metabolite_is_enantiomer_of]-`.
- **Needed because 15 of the 25 enantiomer pairs have no stereo-unspecified node between them**
  (L-/D-Val, -Pro, -Asp, -Thr, -Orn, -Trp, (S)-/(R)-malate). The dominant-form rule already
  treats them as one stereo family, but the graph didn't connect them.
- **Separate from `is_stereo_form_of`**, so the two meanings stay distinct: "a more specific form
  of" (directed) vs "the mirror image of" (symmetric).

**What changes for the explorer:**
- **The F2 stereo hint:** a "0 catalysts; other form has N" lookup should follow
  **both** edge types:
  - `is_stereo_form_of` in either direction (specific ↔ unspecified);
  - `is_enantiomer_of` undirected.

  Otherwise the D-Trp / L-Trp-type pairs are missed.
- **Vocabularies:** new `ChemicalRole.level_kind`, plus `Metabolite_is_enantiomer_of.sources` /
  `.evidence`. `ChemicalRole.role_branches` is not registered. All covered by the hash re-pin.
- The §22 question (`ChemicalRole` in `kg_schema` / `list_filter_values` from day one?) still
  stands.

## 24. Explorer reply on §21–§23 (2026-10-14)

**§21: accepted.** The anomer guard is a good catch. Glucose shows that the 2× margin alone was
protecting us by luck. The 0-move result and the drift test close the item.

**§22 / §23: the `ChemicalRole` model is better than `substrate_roles`. Accepted.**
- Putting the meaning in the edge type, giving it the same `substrate_depth` as the metabolite
  edges, and leaving room for a future `Metabolite_has_chemical_role` on the same nodes: all
  three are right.
- The pruned `is_a` hierarchy also works: level 0 is the three branches, `level_kind`,
  best-effort DAG levels like GO, subtree counts with `direct_*` versions.
- One counting detail: with a DAG, subtree counts must count **distinct** transporters and
  families over `*0..`. A role reachable through two parents must not be counted twice. Please
  include a KG test, the same one GO has.

### One ask: give `ChemicalRole` its own id prefix (`chebi.role:`)

- Today both `ChemicalRole` and `Metabolite` use `chebi:`. The "ids never coincide" test keeps
  them apart in the data, but not in **routing**.
- The explorer dispatches on self-prefixed ids. `ontology_term_details` routes `go:`, `tcdb:`,
  `merops.family:`, `kegg.pathway:`, `tigr.role:` and others to a label by prefix, and the
  metabolite resolver treats `chebi:` as a metabolite namespace.
  - A role id passed to a metabolite tool would be resolved as a metabolite and land in
    `not_found`.
  - A future role mode in the ontology tools couldn't tell a role id from a metabolite id
    without probing both labels.
- **Proposal:** `chebi.role:33229`. This matches the house pattern of `<source>.<kind>:`
  (`kegg.pathway:`, `merops.family:`, `tigr.role:`).
  - The ChEBI number still sits in the local part.
  - The `role_id` column in `metabolite_id_remap.tsv` then makes the change of kind explicit:
    `chebi:33229` (removed Metabolite) → `chebi.role:33229`.
  - Deciding now is cheap: nothing reads these ids yet.
- **Explorer side either way:** when a `chebi:` metabolite id misses, the resolver will check
  `ChemicalRole` and warn: "`chebi:33229` is a chemical role (vitamin (role)), not a
  metabolite; it was removed in <release>". This covers users holding old ids.

### Answer to the §22 question

- **`kg_schema`: yes, from day one.** It introspects the live graph, so the new label, edges and
  indexes appear automatically. The explorer's `config/schema_baseline.yaml` gets updated during
  absorb.
- **`list_filter_values`: no.** Its `filter_type`s are the values some explorer tool *filters
  on*, and no tool filters on `ChemicalRole` yet. Role names are an open set of ~111, so
  `chemicalRoleFullText` is the right way to find them, not a value listing. `ChemicalRole.level_kind`
  will join `list_filter_values` when a tool exposes it.
- **Tool route for "which families transport siderophores":** `run_cypher` for now.
  - The explorer will add a few-shot example to `kg/queries.py` and a short note in
    `docs://analysis/metabolites`. The example: full-text → `Chemical_role_is_a_chemical_role*0..`
    → `Tcdb_family_transports_chemical_role` → TcdbFamily → genes via `genes_by_ontology(ontology='tcdb', term_ids=…)`.
  - A first-class tool mode, likely a role arm on `genes_by_metabolite` or a `ChemicalRole`
    entry in the ontology tools, waits until real questions ask for it. That is another reason
    to settle the prefix now.

### §23.2 enantiomer edge

Agreed. The F2 stereo hint will follow `Metabolite_is_stereo_form_of` in both directions and
`Metabolite_is_enantiomer_of` undirected, and will
also take the two-hop sibling path through an unspecified node (specific → unspecified → other
specific, e.g. D-Ala → Ala → L-Ala).

Please confirm that anomer pairs are **not** `is_enantiomer_of` edges. They are diastereomers,
and they stay on `is_stereo_form_of` per §21.

## 25. KG FYI: paper-metabolite names match case-insensitively, except formulas (2026-10-14)

Spec §3.6. No explorer action beyond the absorb you already plan.

**Rule:**
- Paper-name keys (alias overrides, `name_lookup`, the dominant-form lock, the review lists) are
  case-folded and whitespace-normalised.
- **Exception:** a name that parses entirely as element symbols with counts or a charge (`CO2`,
  `NH4+`, `Co`, `NO`) keeps its exact case (`Co` cobalt ≠ `CO` carbon monoxide). For such names
  the MNX resolver's lowercased name step is skipped.
- Acronyms (`NAD`, `ADP`, `GABA`) are folded normally.

**Why it matters for you:** case-sensitivity was splitting one compound across papers. For
example, "Glutamic acid" and "glutamic acid" landed on D-Glu and L-Glu; Arginine and Serine landed
on an L-form in one paper and an `mnx:` placeholder in another. **Eight compounds converge** on one
node each, and every moved row is in `metabolite_assay_repoints.tsv`. None of today's paper names
is a true formula, so the exception is a guard for future papers.

## 26. KG reply on §24 (2026-10-14, decided)

All three adopted (spec §3.4 / §3.5 / §3.6):

- **`chebi.role:` prefix.** `ChemicalRole` ids are `chebi.role:<primary ChEBI number>`, e.g.
  `chebi.role:33229` "vitamin (role)". The prefix is house-minted and follows the
  `<source>.<kind>:` pattern.
  - KG test: every `ChemicalRole.id` uses it, with a primary ChEBI local part.
  - The remap TSV's `role_id` column shows the change of kind: `chebi:33229` → `chebi.role:33229`.
  - Your resolver fallback ("this is a chemical role, not a metabolite") is a good complement.
- **DAG-distinct subtree counts.** Confirmed: `transporter_count` / `tcdb_family_count` count
  **distinct** families over `Chemical_role_is_a_chemical_role*0..`, so a family reached through
  two parents counts once. There is a KG test in the GO shape, plus subtree ≥ direct.
- **Anomers are not enantiomer edges.** Confirmed and tested. `Metabolite_is_enantiomer_of` comes
  only from ChEBI `is_enantiomer_of`, and anomers (diastereomers) stay on `is_stereo_form_of`.
  Measured: 0 of the 25 pairs are anomers.
  - **Side finding that tightens the anomer guard:** the one pair with a `beta-` in its name,
    beta-Citronellol ↔ (−)-Citronellol, is a genuine enantiomer pair; the "β" is part of a
    trivial name.
  - So anomer detection can't rely on a name prefix alone. It uses ChEBI's anomer annotation, or
    an `alpha-`/`beta-` prefix **plus** a ring-sugar / glycoside ChEBI class.
  - beta-Citronellol is a negative test case.

Your §22 answers are noted: `kg_schema` from day one, not `list_filter_values`, `run_cypher` plus a
few-shot example for the siderophore route, and the two-hop sibling path in the stereo hint.

## 27. Explorer sign-off, final (2026-10-14)

§25 and §26 are accepted, and the explorer has no open questions. **This supersedes the §13
sign-off.** The citronellol negative case is a good addition.

One check for §25: the 8 converged compounds include unqualified names such as "glutamic acid".
They should end up on the dominant form, with `stereo_assignment = 'dominant_form_assumed'` on
every row the case-fold moved, not only on rows the stereo rule moved. Assuming the step-6
decision runs after case-folding, this already holds. If it does, please don't change anything.

### Consolidated explorer work list (replaces the §13 list)

1. **Before the bundled release (a no-op on the current KG): the metabolite resolver.**
   - Normalise each input shape to its CURIE: KEGG; ChEBI (`CHEBI:` / bare digits, leading zeros stripped); MNX; HMDB (5→7-digit pad); InChIKey (new).
   - Match on `Metabolite.id` first, then `IN coalesce(m.all_identifiers, [])`.
   - Report the mapping in `resolved_aliases`, and keep the expand-and-warn fallback.
   - When a `chebi:` id misses, check `ChemicalRole` and warn that it's a role, not a metabolite.
2. **Same release as the KG change:**
   - `list_metabolites`: `structure_kind` and `kegg_compound_categories`, each with a per-row field, a filter and an envelope breakdown.
   - The three assay tools: `stereo_assignment` per row (compact), a `by_stereo_assignment` breakdown, and a warning when any row is `dominant_form_assumed`.
   - The stereo hint on `genes_by_metabolite` / `metabolites_by_gene` / `list_metabolites`. It follows `is_stereo_form_of` in both directions, `is_enantiomer_of` undirected, and the two-hop sibling path.
3. **Absorb:**
   - Goldens and examples, from `metabolite_id_remap.tsv` and `metabolite_assay_repoints.tsv`.
   - Re-pin the vocabulary hash and update `schema_baseline.yaml` (new labels `MetaboliteClass` / `ChemicalRole` and their edges).
   - Switch the `list_metabolites` xref filters to `*_ids` array membership.
   - Docs:
     - the CURIE vs bare `all_identifiers` asymmetry (`docs://guide/conventions`);
     - `metaboliteFullText` now covers name + formula + ids (`list_metabolites.yaml`);
     - `transport_substrate_resolution` text points at `substrate_breadth`;
     - KEGG-not-textbook class caveats.
   - Re-measure the `3.A.1` figures in `docs://analysis/metabolites`.
   - Add the `ChemicalRole` siderophore route as a few-shot example in `kg/queries.py`.

## 28. KG confirmation (2026-10-14)

Your §25 check holds as specified, with no change:
- `name_key` case-folding applies wherever a paper name is keyed, from the step-6 harvest
  onwards, and the dominant-form decision (step 6.9) runs after it.
- `stereo_assignment` records what the paper *said* (§3.6), not how the row got there. So every
  row with an unqualified name such as "glutamic acid" lands on the dominant form with
  `dominant_form_assumed`, whether case-folding, the stereo rule, or both moved it.
- The repoint TSV lists each such row.

Your consolidated work list (§27) is noted. On the KG side the next step is spec sign-off, then
the implementation plan.

## 29. KG-side spec review: explorer-facing changes (2026-10-14, for feedback)

We reviewed the spec for correctness, completeness and consistency after all the decisions above.
Most fixes are internal (pipeline ordering, file shapes, tests). The ones below touch the explorer.
**Please flag anything that conflicts with your §27 work list.**

1. **The dominant-form rule also covers exact single matches onto a stereo-unspecified node.**
   - An unqualified "Glutamate" that matches only `C00302` (69 catalysts) now goes to L-Glu
     `C00025` (1,905), with `dominant_form_assumed`. Until now only ambiguous names and
     paper-supplied ids were covered.
   - It's the same reasoning as §20/§21: an unspecified node doesn't assert a form.
   - **Measured: 0 such rows today**, so it's a guard for future papers. Any hit would appear in
     `metabolite_assay_repoints.tsv`.
2. **"Stereo-qualified name" is now defined**, and it drives `stereo_assignment` (C3).
   - A paper name *states a form* if it contains a stereo descriptor: `L-` / `D-` / `DL-`,
     `(R)-` / `(S)-` / `(+)-` / `(−)-`, `meso-`, or `alpha-` / `beta-` followed by `D-` / `L-`.
   - Otherwise it's unqualified. `beta-Citronellol` counts as unqualified; the β there is part of
     the trivial name.
   - `as_reported` now reads: the paper named the form, or the metabolite has no stereo variant
     **other than anomers**. This is the wording that makes GlcNAc / MurNAc consistent.
3. **The lumping pin isn't fully frozen for new families (refines §19).**
   - Families that were kept when the pin was frozen never change silently; disagreements are
     reported as drift, as agreed.
   - A TC family that becomes kept *later*, e.g. with a new strain, is classified by the threshold
     and appended to the pin (`lumping_set_added`).
   - Otherwise a new 300-substrate family would default to `specific`, and its genes would read
     `resolved`.
   - For you: `substrate_breadth` can gain `lumping` members between releases, always reported.
     Existing members never flip silently.
4. **Removed `mnx:` placeholder ids are *not* in `all_identifiers`** (refines E1).
   - "Arginine hydrochloride" isn't the same entity as L-Arg, so `mnx:MNXM1565340` doesn't go onto
     L-Arg.
   - These 4 ids resolve only through `metabolite_id_remap.tsv` (basis `removed_placeholder`, new
     id = the re-curated target). The resolver will return `not_found` for them.
   - Removed role ids work as you planned: they resolve through your `ChemicalRole` fallback.
   - Goal 12 now covers merged-away and re-keyed ids only.
5. **`chebi_id` tie-break.** 375 KEGG nodes have several ChEBI entity cross-references, typically
   the acid and its anion. `chebi_id` takes the **charge-0** member, then the lowest number;
   `chebi_ids` lists them all. This affects `chebi_id` values and goldens.
6. **B2 scrub also covers two KEGG polymer nodes**: chitosan `C00734` and PHB `C11552`, which
   carried their monomer's InChIKey. They lose `inchikey` / `smiles` / concrete formula and keep
   `elements`, so they stay in element filters.
7. **The KEGG pathway set may grow.**
   - The 89 compounds re-keyed from ChEBI to KEGG bring their compound → pathway links, and that
     can add pathway-level `KeggTerm` nodes with **no gene attached** (transport-reachable, the
     Option B pruning).
   - This could surface in `ontology_landscape`, `pathway_enrichment` term universes and
     `kg_schema` counts as 0-gene pathways. The plan measures it.
   - Gene annotation state must not move (a regression check).
   - Tell us if 0-gene pathways need special handling on your side; today the KEGG layer already
     prunes to KO ∪ reaction ∪ transport-reachable pathways.
8. **`Reaction_has_metabolite` drops slightly** where merged compounds deduplicate on one reaction
   (e.g. HCO₃⁻ + carbonic acid). This affects reaction-level goldens.
9. **New `data_source:chebi` node**, alongside `data_source:kegg`. It carries ChEBI's CC BY 4.0
   attribution in `description`.
10. **New scalar indexes:** `metabolite_structure_kind_idx` (for your `structure_kind` filter) and
    `chemical_role_level_idx` / `_level_kind_idx`. These are additions only.

**Internal only (FYI):**
- re-keying moved earlier in step 6 so re-keyed compounds get KEGG enrichment;
- the assay-collision guard runs on the per-table harvest and is re-checked on the final map;
- per-paper alias overrides are kept per paper (two papers may map one name differently);
- a same-key collision in one alias yaml is an error (`CO2` / `Co2` are different keys and fine);
- no `level_is_best_effort` marker on `ChemicalRole`;
- `validate_paperconfig.py` and the `/paperconfig` skill learn the alias format.

## 30. Explorer reply on §29 (2026-10-14)

§28 confirms the §25 check. Nothing in §29 conflicts with the §27 work list. Items 1, 3, 6, 8, 9
and 10 are accepted as written. Comments on the rest:

**Item 2: the stereo-descriptor definition needs to match inside a name, not only at the start.**
- Many stereo-qualified metabolite names carry the descriptor mid-name:
  - `N-Acetyl-D-glucosamine`
  - `S-Adenosyl-L-methionine`
  - `UDP-D-glucose`
  - `Glycyl-L-leucine`
- If the rule only checks the start of the name, these read as *unqualified*. They would be
  marked `dominant_form_assumed`, and could be moved.
- Please define the descriptor as a token **anywhere** in the name, bounded by start / `-` /
  space / `(` / `,`. Please also cover:
  - **multi-centre locants:** `(2S)-`, `(2R,3S)-`, `(R,R)-`;
  - **bare `R-` / `S-`;**
  - **Greek letters:** `α-D-` / `β-L-`, in addition to `alpha-` / `beta-`;
  - **both minus signs:** U+2212 `(−)` and ASCII `(-)`.
- Please test a negative case too: the `S-` in `S-Adenosyl` is a sulfur locant, not a stereo
  descriptor. The `L-` later in that name is what qualifies it, so the result is right only if
  the matcher doesn't stop at the first hit.
- Suggested KG test cases:
  - qualified: `N-Acetyl-D-glucosamine`, `S-Adenosyl-L-methionine`, `(2R,3S)-isocitrate`
  - unqualified: `beta-Citronellol`, `alpha-Ketoglutarate`, `S-Adenosylmethionine`

**Item 4: removed `mnx:` placeholder ids.**
- Agreed that they mustn't go into `all_identifiers`. "Arginine hydrochloride" isn't L-Arg.
- But a user holding `mnx:MNXM1565340` from an earlier query would get a plain `not_found`,
  even though the measurements it carried still exist on L-Arg.
- **Small ask, not blocking:** put a sparse `recurated_from_ids: str[]` on the target node. It
  stays separate from `all_identifiers` and has no identity meaning.
- The explorer would then check it after an `all_identifiers` miss and warn instead of
  resolving: "`mnx:MNXM1565340` was a placeholder whose measurements were re-curated to
  `kegg.compound:C00062` (L-Arginine) in <release>".
- This mirrors the `ChemicalRole` fallback. If you'd rather not, `not_found` is acceptable for
  4 ids, and the explorer won't hard-code them.

**Item 5: `chebi_id` tie-break.**
- Charge-0 first, then the lowest number, is fine for KEGG nodes.
- Please confirm that on **`chebi:`-keyed nodes** `chebi_id` always equals the id's local part,
  i.e. the node id uses the same tie-break.
- The explorer's `chebi_id` field description promises this ("extracted from the ID itself").
  If the two can differ, we'll reword it instead; please tell us which.

**Item 7: 0-gene KEGG pathways need no special handling on the explorer side.**
- `ontology_landscape`, `pathway_enrichment` and `cluster_enrichment` are gene-anchored (the
  gene → term `MATCH`), so 0-gene terms never enter the term-size distributions or the
  TERM2GENE universes.
- `search_ontology` browse lists them last (`gene_count DESC`), and `min_gene_count` hides them.
- On metabolites (`pathway_ids`, `top_metabolite_pathways`) they are wanted.
- Only `kg_schema` / `kg_release_info` counts and `search_ontology` goldens move, and absorb
  covers that. Please include the measured count in the post-rebuild doc.

## 31. KG reply on §30 (2026-10-14, decided)

**Item 2: adopted, with one exception.** The spec (§3.6) now matches a descriptor **anywhere** in
the name:
- tokens bounded by start / `-` / space / `(` / `,`;
- every token is scanned (no stop at the first hit);
- multi-centre locants: `(2S)-`, `(2R,3S)-`, `(R,R)-`;
- Greek `α-D-` / `β-L-`, in addition to `alpha-` / `beta-`;
- both U+2212 `(−)` and ASCII `(-)`.

Your test cases are in, plus `UDP-D-glucose`, `Glycyl-L-leucine` and `(S)-malate` as qualified,
and `Alanine` / `Glutamate` as unqualified.

**The exception is bare `R-` / `S-`. They are not descriptors.** Counting them would contradict
your own negative case:
- a bare `S-` would make `S-Adenosylmethionine` qualified;
- in metabolite names a bare `S-` is overwhelmingly a sulfur locant (S-adenosyl, S-methyl,
  S-lactoylglutathione);
- real configurations use the parenthesised forms, which are covered.

If you know of real names that use bare `R-`/`S-` as a configuration, send them and we'll
special-case them.

**Item 4: adopted.** Sparse `Metabolite.recurated_from_ids: str[]` on the target of each removed
`mnx:` placeholder (e.g. `mnx:MNXM1565340` on L-Arg).
- It has no identity meaning and is never in `all_identifiers`.
- KG test: each placeholder id appears in exactly one node's `recurated_from_ids`.
- Your warn-instead-of-resolve fallback fits.

**Item 5: confirmed.** On `chebi:`-keyed nodes, `chebi_id` **always** equals the id's local part:
- the canonical-id rule and the `chebi_id` tie-break are the same (charge-0, then the lowest
  number);
- there's a KG test for it;
- your field description ("extracted from the ID itself") stays correct.

**Item 7: noted.** The post-rebuild doc will report the measured count of 0-gene pathway-level
`KeggTerm` nodes added.

## 32. Explorer close (2026-10-14)

§31 is accepted in full.
- **Bare `R-` / `S-`:** you're right, and our ask contradicted our own negative case. Excluding
  them is correct; we know of no counter-examples.
- **§27 work list, item 1:** add one step. After an `all_identifiers` miss, check
  `recurated_from_ids` and warn instead of resolving: "placeholder re-curated to <id>". This
  sits next to the `ChemicalRole` fallback.

The explorer has no open items on this change. Over to the KG for spec sign-off and the plan.
We'll pick it up again with the post-rebuild doc (the two TSVs, the per-name stereo table, the
generic/polymer counts and the 0-gene pathway count).

## 33. KG FYI: step 7 folded into step 6; measurement-only nodes dropped when their rows move (2026-10-14)

Two simplifications, found while tracing the paper-metabolomics flow. Spec §5.2 (6.4 / 6.9a)
and §5.3.

1. **Step 6 now writes the paper `_resolved.csv` files itself; step 7 is retired.**
   - Before, step 6 decided and step 7 wrote. With ids now moving, re-running step 6 alone would
     have left stale tables, and those give dangling assay edges that `neo4j-admin` drops
     silently. One step that decides and writes can't drift.
   - For you: nothing changes in the `_resolved.csv` columns or the edges, apart from the
     `stereo_assignment` column already agreed.
2. **Paper metabolites are added only after their final resolution.**
   - The harvest now collects candidates. Only the chosen targets become paper metabolites and
     get the `metabolomics` tag.
   - A candidate that isn't chosen and has no catalysis or transport reason to exist is never
     written. About **11 nodes** go this way today: D-Trp `C00525`, `chebi:173245` Tyrosine,
     the ChEBI vitamin-B2 / B5 nodes, SAM `chebi:142094`, cysteic acid, chitotriose, and the 4
     `mnx:` placeholders.
   - The plan re-measures this. Overrides or the rule move some, and some stay if they're still a
     paper's target.
   - Otherwise these nodes would have been left with no edges.

**What changes for the explorer:**
- **`recurated_from_ids` covers more than placeholders.** It now covers every dropped node whose
  measurements moved, e.g. `kegg.compound:C00525` D-Trp on L-Trp. Your §32 fallback ("re-curated
  to <id>") covers these with no change. They are **not** in `all_identifiers`; D-Trp isn't L-Trp.
- **New remap basis value** `remeasured_elsewhere`, with new id = the target.
- **About 7 nodes lose the `metabolomics` tag** but keep their other evidence: D-Ala, D-Glu,
  D-Gln, Asn `C16438`, Met `C01733`, Phe `C02057`, and others. Their `measured_*` rollups go to 0.
- **Node estimate ≈ 3,214.**

## 34. Explorer last look (2026-10-14)

**§33: accepted.** Having one step that decides and writes removes the stale-table / dangling-edge
risk, and not creating orphan nodes is right. The explorer needs no change for
`remeasured_elsewhere` / `recurated_from_ids`. Two small points:

- **Warning wording (explorer side, no action).** `recurated_from_ids` now holds real compounds
  (D-Trp `C00525`), not only placeholders. So the explorer's warning will say "**removed**: its
  only edges were paper measurements, now assigned to <id>". It won't say "re-curated to", which
  could read as identity.
- **Order: stereo and enantiomer edges vs dropped nodes.** D-Trp is one of the §23 enantiomer
  pairs (L-/D-Trp). Please confirm that `Metabolite_is_stereo_form_of` /
  `Metabolite_is_enantiomer_of` are emitted **after** the measurement-only drop, so no edge
  points at a dropped node. Please also report the final 60 / 25 edge counts in the post-rebuild
  doc.

### Ask before sign-off: refresh the top of this doc

The header, Summary and §1 tables still describe the 2026-10-13 proposal. A reader of the
post-rebuild doc, or of the explorer absorb, shouldn't have to replay §4–§33 to know what
shipped. Please update the following, or collapse §4–§33 into an appendix:

- **Status:** "proposed, for explorer review" → "explorer-reviewed, signed off (§32/§34)".
- **Summary item 3, the §1 role-terms row, and the line under the property table:** these still
  say `TcdbFamily.substrate_roles`. That was replaced by `ChemicalRole` +
  `Tcdb_family_transports_chemical_role` + `Chemical_role_is_a_chemical_role`, with the
  `chebi.role:` prefix (§22/§23/§26).
- **The §1 property table is missing:**
  - `recurated_from_ids`;
  - that `all_identifiers` is CURIE-only (§12);
  - the `chebi_id` tie-break (§29.5).
- **Add a short "Also new" list:**
  - `Metabolite_is_stereo_form_of`;
  - `Metabolite_is_enantiomer_of`;
  - `stereo_assignment` on both assay edges;
  - `TcdbFamily.substrate_breadth` and the lumping pin;
  - the dominant-form rule (including the anomer guard and the stereo-descriptor definition);
  - `metaboliteFullText` = name + formula + identifiers;
  - the two TSVs and their basis values.

With that, the explorer has nothing further. Thanks; this was a thorough round.

## 35. KG reply on §34 (2026-10-14)

- **Edge order: confirmed and made explicit.** Stereo families are built in memory (spec 6.8)
  for the dominant-form rule. The `Metabolite_is_stereo_form_of` / `Metabolite_is_enantiomer_of`
  lists are emitted **after** the measurement-only drop, and only between surviving nodes (new
  stage 6.9b). So L-Trp ↔ D-Trp loses its edge once D-Trp `C00525` is dropped, and no edge points
  at a dropped node. The final counts (≈ 60 / ≈ 25 before the drop) go into the post-rebuild doc.
- **Warning wording:** noted ("removed; its paper measurements now live on <id>"). The top of the
  doc uses it.
- **Top of this doc refreshed.** The header, Summary, §1 and §2 now describe what ships:
  - `ChemicalRole` replaces `substrate_roles`;
  - `recurated_from_ids`, CURIE-only `all_identifiers` and the `chebi_id` tie-break are in the
    property table;
  - new §2a covers `ChemicalRole`, and §2b is the "Also new" list you asked for;
  - everything from §3 on is now an appendix of review history, with "the top is authoritative".

## 36. KG final review (2026-10-14): three items for explorer sign-off, plus top-of-doc fixes

Two independent reviews of the spec and this doc found real issues. The KG-internal fixes (step-6
stage ordering, re-keying on the paper path, buffered writes, ChEBI cache coverage) are in the
spec. **Three items touch you**, and the top of this doc now reflects all of them:

**A. Old role ids need a lookup list, not a prefix swap.** The removed role metabolites mostly
carry TCDB's **secondary** ChEBI ids: `chebi:10004` "vitamin (role)" becomes
`chebi.role:33229`. Your §24 fallback (`chebi:` miss → try `chebi.role:` with the same number)
would miss them.
- **Decision:** `ChemicalRole` gets a dense `all_identifiers` (its own `chebi.role:` id + the
  `chebi:` CURIEs of its primary and secondary ChEBI ids). These are the same ChEBI term, so
  identity holds.
- **Your fallback:** test `IN ChemicalRole.all_identifiers` instead of a prefix swap.
- Entries are unique across `Metabolite` ∪ `ChemicalRole` (KG test).
- The remap TSV example is corrected to `chebi:10004 → chebi.role:33229`.

**B. `substrate_breadth` gains a third value, `not_assessed`.** As a dense two-state property it
would have labelled class-level nodes (`3.A`, hundreds of substrates) "specific", although the
lumping rule never assessed levels < 2. Values are now `lumping` | `specific` (level ≥ 2, not
pinned) | `not_assessed` (level < 2). It's a closed vocabulary with value descriptions.

**F. Please confirm: `MetaboliteClass` shares the `kegg.brite:` prefix with gene BRITE.** This is
the same routing concern you raised for `chebi:` in §24. Ids can't collide (`br08xxx` vs
`ko0xxxx` tree codes). Our read of your code is that BRITE is routed by ontology name, not by id
prefix, so we kept it. If anything routes on `kegg.brite:`, tell us and we'll switch to a
distinct prefix.

**Also changed in the top part** (detail you'd otherwise have to reconstruct from the appendix):
- the full `stereo_assignment` semantics, including single exact matches and the unspecified node
  winning;
- the `all_identifiers` string-form rules;
- the full stereo-descriptor list;
- the per-row resolution order (an override beats the row's direct id);
- ambiguous names fail the build, while unmatched names are unresolved;
- "7 paper **names** move", not rows;
- `substrate_breadth` detail: the 13 ids, and that the threshold no longer reproduces the set;
- the `ChemicalRole` count semantics;
- the count-shift list, including `Tcdb_family_transports_metabolite`, `TcdbFamily.metabolite_count`,
  `measured_*`, and the `name` value change;
- the closed repoint-TSV reason values;
- the per-name table columns;
- the `metaboliteFullText` fallback (and that `formula` is new to the index);
- the stereo edge counts are pre-drop;
- `identity_basis` now also on re-keyed nodes;
- `MetaboliteClass` edge counts.

Section range fixed (3–36).

## 37. Explorer final sign-off (2026-10-14)

The refreshed top of the doc (Summary, §1, §2, §2a, §2b) reads correctly as the authoritative
"what ships" description, and it matches everything agreed in the appendix. §35 is accepted.
On §36:

- **A: accepted, and thanks for catching it.** Our §24 prefix-swap fallback would have missed the
  secondary-id roles (`chebi:10004`). The explorer will test `IN ChemicalRole.all_identifiers`
  as written in §1 *Migration notes*, step 5. Uniqueness across `Metabolite` ∪ `ChemicalRole`
  is what makes the two-label lookup safe.
- **B: accepted.** `not_assessed` is better than labelling a level-1 class "specific". The
  explorer docs will tell users to filter on `substrate_breadth = 'lumping'`, and that the
  threshold no longer reproduces the set.
- **F: confirmed, keep `kegg.brite:`. Nothing in the explorer routes on it.**
  - `ontology_term_details` finds a term by matching its id against the registry labels
    (`BriteCategory` among them), not by prefix. A `MetaboliteClass` id passed there comes back
    `not_found`, which is correct until a `MetaboliteClass` mode exists.
  - `brite_tree` / `ontology='brite'` read `BriteCategory.tree`, which is label-scoped.
  - The `chebi.role:` case was different: the metabolite resolver *does* dispatch on `chebi:`.
  - If a future `MetaboliteClass` term mode is added to the registry, the label match still
    separates the two, and the `br08xxx` / `ko0xxxx` tree codes can't collide.

**Explorer sign-off: final.** The §27 work list stands, with these amendments:
- the resolver uses `recurated_from_ids` with "removed; its paper measurements now live on <id>";
- the resolver uses `ChemicalRole.all_identifiers`, not a prefix swap;
- the docs reference `substrate_breadth`, including `not_assessed`.

Next contact: the post-rebuild doc.
