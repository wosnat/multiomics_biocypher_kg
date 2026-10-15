# Metabolite identity + KEGG compound classes — design

**Date:** 2026-10-13 (restructured as a data-flow spec 2026-10-14)
**Status:** explorer **final sign-off 2026-10-14** (review doc §37, which closes §36 items A / B /
F). Decisions taken in review (labels refer to review-doc sections: **E1–E3** = identifier
resolution / scalar cross-references / assay collisions; **F1/F2** = `elements` on class nodes /
stereo landing; **C1–C3** = freeze / margin / `stereo_assignment` semantics):
- E1–E3, F1/F2;
- prefixed-only `all_identifiers`;
- the dominant-stereo-form rule with C1–C3, the anomer guard, and paper-supplied ids;
- case-insensitive name keys (formula exception);
- the role-substrate model: a `ChemicalRole` hierarchy under the `chebi.role:` prefix with
  `Tcdb_family_transports_chemical_role`, which replaced the earlier `substrate_roles` property;
- `Metabolite_is_enantiomer_of`;
- the pinned lumping set, with `substrate_breadth ∈ {lumping, specific, not_assessed}`;
- `recurated_from_ids`; bare `R-` / `S-` excluded as stereo descriptors;
- step 7 folded into step 6; measurement-only nodes dropped when their rows move; stereo edges
  emitted after that drop.

Revised after two KG-side reviews on 2026-10-14 (correctness, gaps, consistency, code
claims). **User sign-off 2026-10-14.** Next: implementation plan (writing-plans).
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture (2026-10-13)":
"Metabolite equivalence across namespaces and stereo/charge variants" (P1) and "Chemical
categories for `Metabolite`" (P2).
**Explorer review doc:** `docs/kg-changes/metabolite-identity-and-classes.md`
**Deferred (separate spec, after this lands):** TCDB substrate plausibility, `substrate_depth`
semantics, per-organism transport reach, TCDB name hygiene, lumping-threshold re-tune.

All numbers were measured on 2026-10-13/14 against the live KG (built 2026-10-13 06:28), the
committed step-6 caches, MNX (`~/tools/mnx`), ChEBI release 256 and the KEGG REST compound BRITE
files. The plan re-measures them.

**How to read this spec.**
- §1–§2: why, and what "done" means.
- §3: the **rules**, i.e. what the data must say, independent of plumbing.
- §4–§9: the **data flow**: inputs → prepare-data steps → cached / resolved files → schema →
  adapters → post-import.

---

## 1. Problem

The KG mints `Metabolite` nodes from three arms, and MNX is the only bridge between the KEGG and
ChEBI namespaces:
- catalysis (KEGG reactions);
- transport (TCDB substrates, ChEBI-keyed);
- metabolomics (paper assays).

Today: 2,717 `kegg.compound:` + 635 `chebi:` + 4 `mnx:` = 3,356 nodes.

### 1.1 Bugs found

| # | Bug | Size | Cause |
|---|---|---|---|
| B1 | ChEBI ids are **secondary/obsolete** | 1,627 of 3,042 ChEBI ids (48%): 1,454 `chebi_id` values on KEGG nodes + 173 `chebi:` node ids | the lowest-sorted MNX `compound_aliases` row is kept, and it includes secondary ids; TCDB's `substrates.tsv` carries old ids verbatim |
| B1b | `chebi_id` names an **isotopologue / odd variant** | 331 KEGG nodes disagree with ChEBI's own KEGG xref (Succinate → "succinic acid-d4", Urea → "UREA C 14", Calcium → "calcium-43(2+)") | same lowest-alias pick |
| B2 | Class/polymer nodes carry a **representative monomer structure** | 20 class nodes. 5 InChIKeys sit on two nodes for this reason (alditol / glycerol, (1,4-α-D-glucosyl)n / glucose, chitosan / glucosamine, PHB / 3-hydroxybutanoate, chitin / GlcNAc) | MNX stores a representative structure on generic entities, and we copy it |
| B3 | **ChEBI role terms are `Metabolite` nodes** | 57 ("drug", "herbicide", "xenobiotic", "acid", "base", "vitamin (role)", …), 392 transport edges | TCDB lists roles as substrates |
| B4 | **Duplicate nodes** for one metabolite | 69 merge groups (§1.2) | each arm resolves separately; MNX keeps charge forms and generic ChEBI terms apart |
| B5 | **Ambiguous paper-metabolite names are silently assigned** | 36 of 149 paper metabolites. Alanine → D-Ala `C00133`, Glutamic acid → D-Glu `C00217`, Tryptophan → D-Trp `C00525`, Glutamine → D-Gln `C00819`, Cystine → indole-3-acetate `C00954`, Glutathione (GSH) → GSSG `C00127`; Arginine / cysteine / serine → `mnx:` placeholders | step 6 `_harvest_paper_metabolites` maps MNX `xref:ambiguous` to a primary; step 7 `_resolve_row` takes `sorted(hits)[0]` |

B5 has one visible symptom in the live graph. GSH and GSSG rows in Biller 2022 and Kujawinski
2023 both land on `C00127`, so 6 (assay, metabolite, condition, time point) slots carry two edges
with different values. These are the only distinct-name collisions across all paper tables.

Invariants that already hold (no action needed):
- `organism_names` / `pathway_*` sizes match their counts;
- every `evidence_sources` tag is backed by an edge;
- no metabolite is isolated;
- `mass` is absent on all ≈ 420 nodes with a generic formula (420 by the live-graph regex `*` /
  `R` / `)n`; 422 by the stricter step-6 rule that also catches a trailing `n`; the plan reports
  the final figure);
- `elements` holds only real element symbols (42 distinct).

### 1.2 Duplicate inventory, by evidence

| Basis | Pairs | Source |
|---|---|---|
| Same MNX entity on two KEGG ids (HCO₃⁻ / carbonic acid, F6P / β-D-F6P) | 18 | MNX |
| `chebi:` node whose ChEBI term cross-references an existing KEGG node | 49 | ChEBI |
| `chebi:` node whose ChEBI term cross-references a KEGG id *not yet* a node (re-key) | 92 nodes | ChEBI |
| Conjugate acid / base (nitrate `chebi:14654` ↔ nitric acid / `C00244`) | 9 new | ChEBI |
| Tautomers | 3 | ChEBI |
| Exact-name match not covered above | 16 | heuristic: **review only** |
| InChIKey first-block match | 430 | heuristic: **rejected** (ignores stereo: glucose ~ galactose) |

**Group-level result:** union-find over the curated links, with role terms removed and the §3.1
guard applied.
- 69 merge groups, 74 nodes merged away.
- 89 `chebi:` nodes re-keyed to KEGG ids. Of the 92 cross-referencing a not-yet-node KEGG id,
  3 also merge with another member of their group, so they count as merged, not re-keyed.
- 57 role nodes and 4 `mnx:` placeholders removed.
- ≈ 11 measurement-only nodes dropped when their rows move (step 6 stage K), including the 4
  placeholders.
- **3,356 → ≈ 3,214 nodes** (re-measured by the plan).
- The guard rejects 4 ChEBI-only KEGG↔KEGG links, e.g. water `C00001` ↔ hydroxide `C01328`.

### 1.3 No compound classes

The proxies the analysis tried for compound classes were about 50% wrong: KEGG pathway, and the
TCDB family that carries the compound. Two real options were measured:
- **KEGG compound BRITE** (8 trees): curated, 2–4 levels, user-facing (br08001: Peptides /
  Carbohydrates / Nucleic acids / Vitamins and cofactors / Lipids / Organic acids / Hormones and
  transmitters / Antibiotics / Steroids). Covers **867 of 2,717 KEGG metabolites** (26% of all
  nodes). **Chosen.**
- **ChEBI-inferred categories**: against BRITE on the 586 compounds both classify, 80% share at
  least one category, 27% have the exact same set, and there are 434 extra assignments.
  **Rejected.**

---

## 2. Goals and acceptance criteria

1. **One node per metabolite.** No two `Metabolite` nodes share a KEGG id, a primary ChEBI
   entity id, an MNX id or any `all_identifiers` entry.
2. **KEGG id is the identity** wherever any curated source gives one; else primary ChEBI; else
   MNX.
3. **Every ChEBI id is primary** and names the entity, never an isotopologue or a class.
4. **Every merge is attributable** (`identity_basis`). Heuristic matches are reported, never
   applied.
5. **No role terms as `Metabolite` nodes.** A TCDB substrate given only by role is modelled as
   a `ChemicalRole` node reached by `Tcdb_family_transports_chemical_role`. Its meaning ("this
   family transports something TCDB calls a drug") is carried by the edge type, so it can't be
   misread as a role *of* the family or *of* its metabolites.
6. **No class node carries a concrete representative structure**, and **no metabolite loses
   `elements`**: the set of nodes with `elements` only grows.
7. **Curated KEGG compound classes** as a `MetaboliteClass` layer, with `sources` / `evidence`
   on every membership edge. No inferred classes.
8. **All three arms land on canonical nodes.** No dangling edge (`import.report` clean).
9. **No ambiguous paper-name match is applied automatically.** No two distinct names of one
   table land on one node unless an override says so (0 same-slot duplicate assay edges).
10. **Stereo variants are linked, not merged**: `Metabolite_is_stereo_form_of` (specific →
    unspecified, directed) and `Metabolite_is_enantiomer_of` (symmetric, stored once).
11. **Stereo-unqualified paper names go to the dominant form** only with a ≥ 2× margin. The
    assignment is frozen, drift is reported, and the assumption is recorded on the edge
    (`stereo_assignment`).
12. **Every merged-away or re-keyed id still resolves** through exactly one node's
    `all_identifiers`.
    - Removed role ids resolve to a `ChemicalRole` (`chebi.role:`), through the explorer's
      fallback.
    - Nodes dropped because their paper rows moved elsewhere (the 4 `mnx:` placeholders, and
      measurement-only nodes such as D-Trp, stage K) aren't the same entity as their new target. So
      their ids are **not** added to the target's `all_identifiers`. They appear in
      `metabolite_id_remap.tsv` and in the target's `recurated_from_ids`.

---

## 3. Rules (what the data must say)

### 3.1 Identity

**What counts as the same metabolite:**

| Relationship | Merge? | Rationale |
|---|---|---|
| Same MNX entity | **yes**, unless the two KEGG ids are stereo forms of one structure (stereo guard below) | MNX's own reconciliation |
| ChEBI term ↔ its KEGG cross-reference | **yes** | curated by ChEBI |
| Conjugate acid / base (`RO:0018033/34`) | **yes** | one physiological species; KEGG is charge-agnostic |
| Tautomer (`RO:0018036`) | **yes** | same compound |
| Stereoisomers / unspecified vs specified (Glutamate `C00302` vs L-Glu `C00025`) | **no**, linked instead (§3.5) | biologically distinct; KEGG keeps them apart |
| Isotopologue | never a node id | not a metabolite in this KG |
| Generic class vs KEGG generic ("peptide" vs `C00012`) | only via a curated link above | name-only matches go to review |

**Guard.** Two different KEGG ids merge **only** when MNX says they are one entity, or when a
curated `merge:` entry in `config/metabolite_identity_overrides.yaml` says so (recorded as
`identity_basis: [manual]`). A ChEBI-only link joining two KEGG ids is reported, not applied. Likewise, a merge that would put two distinct
rows of one assay table onto one node is refused and reported (today: 0 such groups).

**Stereo guard on MNX links** (explorer ask 2026-10-15). MNX sometimes calls a stereo-unspecified
compound and its specific form, or a generic sugar and its anomer, one entity. An MNX link between
two groups with different KEGG ids is refused (reported in `rejected_links` and
`stereo_refused_links`) when some KEGG pair across them is *stereo-distinct*: a non-class ChEBI
term of one shares an InChIKey connectivity block with the other's terms but carries a stereo layer
none of them has, and the two term names are not a cis/trans or E/Z pair (§3.6: geometric isomers
are separate compounds; KEGG lists generic and trans caffeate as two ids of one compound). It uses
ChEBI's per-KEGG-id terms, not the node's MNX InChIKey, which every member of an MNX group shares.
`manual` merges are exempt. Refused pairs are linked by §3.5 instead. A class term MNX joins to one
member, or two MNX-conflated protein-bound species, can't be seen structurally and are curated
`keep_separate:` entries.

**MNX fields on a refused pair.** Both sides still carry one `mnxm_id`, and with it the same
`inchikey`, `smiles`, `hmdb_id`, `chebi_id` (MNX's ChEBI xref) and MNX ChEBI aliases. They stay on
the side whose own ChEBI structure equals MNX's InChIKey (protonation ignored); with no single
match, on the side with more gene-reachable reactions. The other side drops them, takes `inchikey`
from its own ChEBI term and keeps `formula` / `mass` (report: `mnx_split`).

**Canonical id per group** (margin added 2026-10-15, explorer ask; same 2× as the dominant-form
rule, §6.5):
1. If the group has KEGG ids:
   1. a `canonical:` pin in `metabolite_identity_overrides.yaml`, if present;
   2. else the KEGG id with the most gene-reachable reaction references, **only if** it has
      ≥ 2× the runner-up's (1 vs 0 passes);
   3. else, only when that measure is 0 vs 0, the same ≥ 2× test on all KEGG reaction references;
   4. else step 6 fails with `canonical_margin_unresolved`, naming the members, their KEGG names,
      both counts and a ready-to-paste `canonical:` entry. A silent tiebreak would rename a hub
      metabolite unannounced; a close call is the case a human should decide.
2. Else `chebi:<primary>`, preferring the charge-0 member, then the lowest number.
3. Else `mnx:<MNXM>`.

**Heuristic matches** (16 exact-name pairs) go to the identity report. A human decision is
recorded in `config/metabolite_identity_overrides.yaml` (`merge:` / `keep_separate:` with a
`reason`). An override-applied merge records `identity_basis: [manual]`.

### 3.2 Identifiers and cross-references

| Property | Rule |
|---|---|
| `id` | canonical CURIE (§3.1) |
| `kegg_compound_id`, `chebi_id`, `mnxm_id`, `hmdb_id` | **scalars, kept, bare**, as today (`C00244`, `17632`, `MNXM732399`, `HMDB0002878`), holding the canonical entity's value. `chebi_id` comes from ChEBI's own KEGG xref first (an entity, not a class term such as "nitrates"), else MNX's non-secondary, non-isotopologue alias. **Tie-break** when several entity terms qualify (375 KEGG nodes have > 1, typically the acid and its anion): the charge-0 member, then the lowest number. **On `chebi:`-keyed nodes `chebi_id` always equals the id's local part.** The canonical-id rule (§3.1) uses the same tie-break, so the two can't differ (KG test) |
| `kegg_compound_ids`, `chebi_ids`, `mnxm_ids` | **new arrays, bare** (same format as the scalar, scalar value first). `kegg_compound_ids` / `chebi_ids` are present exactly when the matching scalar is (absent on a node with no KEGG / ChEBI id); `mnxm_ids` is sparse (present only when the group spans > 1 MNX entity) |
| `all_identifiers` | **new, dense, one canonical CURIE per identifier**, covering every member's identifiers and **every merged-away or re-keyed node id** (not removed role / placeholder ids, §2 goal 12). Forms: `kegg.compound:C01353` (upper case) · `chebi:14654` (primary *and* secondary, no leading zeros) · `mnx:MNXM732398` · `hmdb:HMDB0002878` (7-digit form only) · `inchikey:<27 chars>`. Each entry is on exactly one node. **This convention differs from `Gene.all_identifiers`** (bare ids); the schema description says so |
| `recurated_from_ids` | **new, sparse** (explorer review §30): the CURIEs of dropped nodes whose paper rows were re-curated onto this node: the `mnx:` placeholders (e.g. `mnx:MNXM1565340` on L-Arg), and measurement-only nodes dropped by stage K (e.g. `kegg.compound:C00525` D-Trp on L-Trp). It **carries no identity meaning**, is kept separate from `all_identifiers`, and exists only so a consumer can say "this old id's measurements now live here" |
| `identity_basis` | new, sparse: present on merged **and on re-keyed** nodes (a `chebi:` node re-keyed to its KEGG id carries `[chebi_kegg_xref]` even with nothing merged into it). ⊆ `{mnx_same_entity, chebi_kegg_xref, chebi_conjugate, chebi_tautomer, manual}` |
| `name` | KEGG name on KEGG-keyed nodes; the ChEBI primary name on `chebi:` nodes |

**Merged-node lists.** Step-6 list fields (`pathways`, `evidence_sources`) are unioned. Every
post-import rollup is recomputed on the merged node, never summed.

**A ChEBI term several separate KEGG nodes cross-reference** (ChEBI xrefs `C00725` lipoate to
(R)-lipoic acid as well as `C16241` (R)-lipoate; the class `C00193` "Aromatic aldehyde" to
benzaldehyde) is published on at most one node (report: `shared_chebi_terms`, term → the nodes it
is withheld from). Its owner is the single node for which it is the only entity term among its
KEGG ids' xrefs, with no class term among them: (R)-lipoic acid → `C16241`, benzaldehyde →
`C00261`, α-D-galactose → `C00984`. With no single owner (hydroquinone `C00530` / `C15603`) it is
withheld from all of them. A ChEBI term that is a written KEGG node's own xref is never published
on a sibling through MNX aliases. (Until 2026-10-15 every shared term was withheld from all nodes.)

### 3.3 Structure

| Case | Rule |
|---|---|
| `structure_kind` (new, dense) | `defined` \| `generic` \| `polymer`. `generic`: the formula carries `R` / `*`, or the ChEBI term is a class. `polymer`: the formula / name marks a repeat unit. The exact split rule is settled in the plan |
| Generic formula already present (422 nodes) | kept, with `elements` derived from it |
| Class node with **no** formula (316) | filled from ChEBI `generalized_empirical_formula` when present (≈ 49 gain one, 44 of them transport substrates), with `elements` derived from it |
| Class node with a **concrete representative** structure (B2: 20) | `inchikey` / `smiles` / concrete formula / `mass` **dropped**; formula refilled from ChEBI's generalized formula where present; **`elements` kept** (class composition = representative's) |
| `mass` | absent whenever the formula is generic |
| `elements` | ⊆ the periodic table: an allow-list applied after `chemparse`, with `R` treated as a placeholder before tokenising (a generic `Rn` is never radon) |

**B2 detection.** A node is scrubbed (its concrete `inchikey` / `smiles` / formula / `mass` are
dropped, and `elements` is kept) when **either** condition holds:
- its `structure_kind` is `generic` or `polymer`;
- its InChIKey equals that of a `defined` node. The generic side is scrubbed, never the defined
  one.

The second condition catches the KEGG polymer nodes that carry a monomer key:
- chitosan `C00734` / glucosamine;
- PHB `C11552` / 3-hydroxybutanoate;
- along with the ChEBI class nodes (alditol, (1,4-α-D-glucosyl)n, chitin).

This makes InChIKey unique by construction. Any leftover collision fails the build.

### 3.4 Role terms: role-typed substrates

**Semantics.** TCDB's substrate list for a family sometimes names a substrate only by its ChEBI
**role** ("drug", "siderophore", "xenobiotic"), not by a molecule. Such an entry says **"this
family transports something TCDB calls a drug"**: an additional substrate, specified by role.
It does **not** say:
- that the family has that role; or
- that the family's `Tcdb_family_transports_metabolite` targets have that role.

The model carries this meaning in its structure rather than in documentation alone.

**Detection.** A TCDB substrate is a role when its **primary** ChEBI term has `CHEBI:50906`
("role") among its `is_a` ancestors (after secondary → primary mapping). A role never becomes a
`Metabolite` or a `Tcdb_family_transports_metabolite` edge.

**Model:**
- **Node `ChemicalRole`: a pruned ChEBI role hierarchy** (decided 2026-10-14).
  - Contents: the 57 roles TCDB uses, plus their `is_a` ancestors up to ChEBI's three role
    branches. That is ≈ 111 nodes; the plan re-measures. The single `role` root (`CHEBI:50906`)
    is **not** emitted (BRITE no-tree-root precedent).
  - The three branches are the level-0 nodes: `biological role`, `chemical role`, `application`.
    So no separate branch property is needed; a role's branch is its level-0 ancestor. A role can
    sit under several branches (10 of the 57 do), because the hierarchy is a DAG.
  - Id: **`chebi.role:<primary ChEBI number>`**, e.g. `chebi.role:33229` "vitamin (role)" (not the
    secondary `10004` TCDB used). A distinct, house-minted prefix, decided 2026-10-14 (explorer
    review §24). It follows the `<source>.<kind>:` pattern (`tigr.role:`, `merops.family:`,
    `kegg.pathway:`), so tools that route by prefix never confuse a role with a `chebi:`
    metabolite. The ChEBI number stays in the local part.
  - Properties: `name` (ChEBI name verbatim, R1), `level` (int, the **shortest** `is_a` distance
    to a branch, 0 = branch), `level_kind ∈ {chebi_role_branch, chebi_role}`.
  - Because the hierarchy is a DAG, `level` is approximate. **No `level_is_best_effort` marker is
    added**: the presence-marker exception is closed to its two sanctioned properties.
- **Edge `Chemical_role_is_a_chemical_role`** (child → parent, ChEBI `is_a`), within the pruned
  set.
- **Why a hierarchy.** Subtree queries catch sub-roles a flat lookup misses: "siderophore"
  catches **phytosiderophore**, "insecticide" catches **neonicotinoid insecticide**; both are
  among today's 57. It is also the backbone for the backlogged molecule roles ("all cofactors"
  including coenzyme / prosthetic group).
- **Edge `Tcdb_family_transports_chemical_role`** (TcdbFamily → ChemicalRole), 392 today on 116
  families. It mirrors `Tcdb_family_transports_metabolite`:
  - **Same rollup** (the full-hierarchy `subtree_substrates` rule, filtered to kept TC ids).
  - **Same `substrate_depth`** per pair (`most_specific` | `inherited`, computed by the same
    `_compute_substrate_depth`, refactored to take the `{tc_id: ids}` map as an argument; today
    it is hard-wired to `self._subtree_substrates`). By depth today: class 91 · subclass 98 · family 120 ·
    subfamily 79 · specificity 4. Up to 33 roles per node.
  - The deferred TCDB spec revises both edge types together.
- **Kept separate from future molecule roles.**
  - The backlogged `Metabolite.chebi_roles` becomes `Metabolite_has_chemical_role` →
    the *same* `ChemicalRole` nodes, through a *different* edge type.
  - "Has a role" and "transports something with a role" stay distinct claims that still join
    on one node, e.g. "families that transport something with a role this molecule has".
- **Not propagated.** No Gene, OrganismTaxon or Metabolite rollup counts role edges:
  - `transported_metabolite_count`, `Organism_has_metabolite` and `TcdbFamily.metabolite_count`
    stay molecule-only;
  - the lumping pin (§3.4a) is unaffected.

  A KG test enforces this.
- **Id prefix.** `ChemicalRole` ids use `chebi.role:` and never `chebi:`, so they can't collide
  with or be routed as `Metabolite` ids. A KG test asserts every `ChemicalRole.id` starts with
  `chebi.role:` and that its local part is a primary ChEBI id. The plan checks that
  `normalize_curie` leaves `chebi.role:` intact.
- **`ChemicalRole.all_identifiers`** (decided 2026-10-14). The role `Metabolite` nodes being
  removed mostly carry TCDB's **secondary** ChEBI ids, e.g. `chebi:10004` "vitamin (role)",
  whose primary is `CHEBI:33229`. A prefix swap (`chebi:10004` → `chebi.role:10004`) would miss.
  So each `ChemicalRole` carries a dense `all_identifiers: str[]`:
  - its own `chebi.role:<n>`;
  - the `chebi:` CURIEs of its primary **and** secondary ChEBI ids, which includes every removed
    role-metabolite id (same ChEBI term, so identity holds).

  Consumers resolve an old `chebi:` role id by testing `IN ChemicalRole.all_identifiers`. Each
  entry is on exactly one node across `Metabolite` ∪ `ChemicalRole` (KG test).

### 3.4a Substrate-lumping set: pinned (decided 2026-10-14, explorer review §17)

**Problem.** Identity work alone would change which TCDB families count as substrate-lumping
(`level ≥ 2 AND metabolite_count ≥ 50`, the predicate behind the stored
`Gene.transport_substrate_resolution`). P-type ATPase `3.A.3` lists the same substrates as
before, but:
- one of its 50 is a role (§3.4);
- two of its substrates merge into one node.

It drops to 47 molecules (48 if the role is counted) and would leave the set. Genes whose
deepest attachment is `3.A.3` would then read `resolved` despite ~48 listed substrates. That's
a threshold artefact, not new information. Counting roles in the breadth isn't enough; the merge
alone still flips it.

**Rule.** The lumping set is **pinned**.
- Step 6 writes a committed **`lumping_pin`** block into `tcdb_pruned.json` (§6.3). It is
  **seeded from the literal list** of today's 13, not recomputed (recomputing on canonical data
  would give 12, since `3.A.3` drops to 47): `2.A.1`, `2.A.1.1`, `2.A.1.2`, `2.A.1.3`, `2.A.3`, `2.A.6`, `2.A.6.2`, `2.A.7`,
  `2.A.66`, `2.A.66.1`, `3.A.1`, `3.A.1.1`, `3.A.3`. The plan re-confirms the list against the
  live graph before freezing.
- On every build, step 6 also recomputes the threshold set from the canonical data, and treats
  the two kinds of TC id differently:
  - **Ids already in `kept_tcdb_ids` when the pin was frozen never change classification
    silently.** A threshold disagreement is reported as `lumping_set_drift` and isn't applied.
    Changing the pin for them is a human edit, the same pattern as the dominant-stereo lock
    (§6.6).
  - **A TC id that becomes kept later** (e.g. a new strain brings a new family) is classified by
    the threshold, appended to the pin, and reported as `lumping_set_added`. Otherwise a newly
    kept 300-substrate family would default to `specific`, and its genes would read `resolved`:
    the very over-claim the pin prevents.
  - The pin records the kept-id set it was frozen against (`frozen_kept_ids`), so step 6 can tell
    the two cases apart.
- `tcdb_adapter` sets a dense **`TcdbFamily.substrate_breadth ∈ {lumping, specific,
  not_assessed}`** (decided 2026-10-14):
  - `lumping`: in the pin;
  - `specific`: level ≥ 2 and not in the pin;
  - `not_assessed`: level < 2 (tc_class / tc_subclass). The old predicate never applied there,
    because substrate counts scale mechanically with level; labelling `3.A` "specific" would be
    false.
- Post-import `transport_substrate_resolution` reads `t.substrate_breadth = 'lumping'` instead of
  the inline `level ≥ 2 AND metabolite_count ≥ 50`. The rest of its logic is unchanged.

**R3: a documented exception, stated plainly.** R3 forbids materialising a threshold the node's
own count already expresses; that is why `is_promiscuous` was deleted.
- `substrate_breadth` is a **compatibility pin**: it keeps today's lumping classification stable
  while identity work changes the underlying counts.
- Auto-classifying *newly kept* ids by the threshold **is** a materialised threshold. It's kept so
  new families are judged the same way the pinned ones were.
- Consumers must filter on `substrate_breadth`, **not** re-derive it: the threshold query no
  longer reproduces the set (`3.A.3`: 47 molecules, still `lumping`).
- `vocabulary-contract.md` and the post-import comment record the exception.
- **Sunset:** the deferred TCDB spec re-tunes or replaces the pin, and must retire this
  exception.

**Result:** no gene's `transport_substrate_resolution` changes because of this spec. A KG test
asserts that the `lumping` nodes equal `lumping_tcdb_ids ∩ kept_tcdb_ids`, so a frozen id that
later leaves the kept set doesn't fail it, and the post-rebuild check compares gene
resolution counts with the pre-change baseline.

### 3.5 Stereo families and the stereo link

**Stereo family.** Two canonical nodes are linked when **both** hold:
- their primary ChEBI terms are joined by a direct `is_a` or `is_enantiomer_of`;
- their InChIKeys share the first (connectivity) block.

A **stereo family** is a **connected component** of these pairwise links. The family is
computed **after** the B2 scrub, so a class can't join one. Families are never built between
members of one merge group.

**`Metabolite_is_stereo_form_of`** runs from the stereo-specified node to the unspecified or
broader node (L-Cys `C00097` → Cys `C00736`; D-Cys `C00793` → Cys `C00736`).
- Carries `sources: ['chebi']`, `evidence: 'curated'`.
- ≈ 60 edges; 38 connect a transport-only node to a catalysis node (before the measurement-only
  drop; the final count is reported).
- Emitted only between nodes that survive the measurement-only drop (stage L).
- No rollup crosses it.

**`Metabolite_is_enantiomer_of`** (decided 2026-10-14) links two members of a stereo family that
are ChEBI `is_enantiomer_of` (L-Val ↔ D-Val, (S)- ↔ (R)-malate).
- **Mirror images only.** The edge comes solely from ChEBI's `is_enantiomer_of`. α/β **anomers are
  diastereomers**, never enantiomers, so they stay on `is_stereo_form_of` (§3.6 anomer guard).
  Measured: 0 of the 25 pairs are anomers. The one pair whose name carries a `beta-`,
  beta-Citronellol ↔ (−)-Citronellol, is a genuine enantiomer pair: "β" is part of the trivial
  name, not an anomeric prefix.
- **Symmetric, stored once:** from the lower to the higher canonical id; consumers query it
  undirected (`-[:Metabolite_is_enantiomer_of]-`). A reverse copy would only double counts.
- Carries `sources: ['chebi']`, `evidence: 'curated'`. ≈ 25 edges before the measurement-only drop
  (L-/D-Trp loses its edge once D-Trp is dropped); the final count is reported.
- **Why it's needed.** 15 of the 25 enantiomer pairs have no stereo-unspecified node between
  them (L-/D-Val, -Pro, -Asp, -Thr, -Orn, -Trp, (S)-/(R)-malate, …). Without this edge they would
  be in the same stereo family for the dominant-form rule (e.g. the D-Trp → L-Trp move) yet
  unconnected in the graph.
- A separate edge type from `is_stereo_form_of` keeps the two meanings distinct: "a more
  specific form of" vs "the mirror image of".
- No rollup crosses it.

### 3.6 Paper-metabolite name resolution

**Order per paper row** (the single order used everywhere; decided 2026-10-14):
1. **The paper's own curated override** for the row's name (`metabolite_aliases.yaml`, by
   `name_key`; carried per paper as `paper_overrides`, §6.5). It wins over everything,
   including the row's direct id, so a curator can always fix a row.
2. **Direct id column.** The id is mapped through `id_aliases` to its canonical id, so a
   paper-supplied KEGG id that was merged away still lands correctly.
   - A specified form, or a node with no stereo variant (anomers aside): take it.
   - The stereo-unspecified member of a family: continue to step 3, then step 4b.
3. **The frozen dominant-form lock** (§6.6), keyed by `name_key` for names and by canonical id
   for direct ids.
4. **Match**, by `name_key`, against the global name lookup.
   - a. One candidate, which is a specified form or has no stereo variant (anomers aside): take
     it.
   - b. One candidate that is the **stereo-unspecified** member of a family, or several
     candidates that all belong to one stereo family, **and** the name is stereo-unqualified
     (see below): apply the **dominant-form rule**. The rule's candidate set is **the whole
     family, the unspecified node included**. If the unspecified node itself wins (it has the most
     catalysis), the row stays on it. On first computation the result freezes into the lock.
   - c. Anything else with several candidates is **ambiguous**.
5. Outcome when nothing resolves:
   - **Ambiguous with no override → step 6 fails**, listing the names, so a human curates them.
   - **Unmatched** (no candidate at all) **→ unresolved, not fatal.** No `metabolite_id`, no
     edge, listed in the identity report (`unresolved_paper_names`). Today:
     `s-(1-2-dicarboxyethyl)glutathione pos` (6 Kujawinski 2023 rows). Unresolved rows carry no
     `stereo_assignment`.

**Stereo-qualified names** (revised with explorer review §30). A paper name **states a form**
when it contains a stereo descriptor **anywhere** in the name, as a token bounded by the start of
the name, `-`, a space, `(` or `,`. The matcher scans **every** token; it doesn't stop at the
first one. Descriptors:
- `L-`, `D-`, `DL-`;
- parenthesised configurations, including multi-centre locants: `(R)-`, `(S)-`, `(R,S)-`,
  `(R,R)-`, `(2S)-`, `(2R,3S)-`;
- optical rotation: `(+)-`, `(−)-`, `(-)-` (both U+2212 and the ASCII hyphen-minus);
- `meso-`;
- an anomeric or Greek prefix **followed by `D-` / `L-`**: `alpha-D-`, `beta-L-`, `α-D-`, `β-L-`.

**Geometric isomers (E / Z, cis / trans) and `D(-)` / `L(+)` forms are not descriptors.**
Geometric isomers are separate KEGG ids and are never merged or linked as stereo forms. The
combined `D(-)-` / `L(+)-` notation is caught by its `D-` / `L-` part.

**Bare `R-` / `S-` are not descriptors.** In metabolite names a bare `S-` is overwhelmingly a
sulfur locant (S-adenosyl, S-methyl, S-lactoylglutathione), and counting it would mark
`S-Adenosylmethionine` as qualified. Real configurations use the parenthesised forms above.

Test cases:
- qualified: `N-Acetyl-D-glucosamine`, `S-Adenosyl-L-methionine` (qualified by the later `L-`),
  `UDP-D-glucose`, `Glycyl-L-leucine`, `(2R,3S)-isocitrate`, `(S)-malate`;
- unqualified: `beta-Citronellol`, `alpha-Ketoglutarate`, `S-Adenosylmethionine`, `Alanine`,
  `Glutamate`.

Otherwise the name is **stereo-unqualified**. Qualified names are never moved by the
dominant-form rule (they resolve to the form they name). The same test drives
`stereo_assignment`.

**Per-paper overrides and key collisions.**
- Overrides are per paper, because two papers may mean different things by one name. The global
  name lookup only holds what no paper overrides.
- Inside one paper's alias yaml, two entries that produce the **same `name_key`** but point to
  different ids are an error.
- Entries that produce **different keys** are fine. `CO2` vs `Co2` both parse as formulas, keep
  their case, and are two keys: carbon dioxide and a cobalt term.
- **Allowed shared targets.** Two distinct `name_key`s of one table may land on one id only when
  **both** names have entries in that paper's alias yaml pointing to that id (the curator said so
  explicitly). Any other same-table collision fails step 6.
- Each paper resolves on its own: the current harvest's global dedup
  (`if name in alias_to_primary: continue`, which lets the first paper's resolution win for every
  later paper) is removed.

**An ambiguous match is never applied.**

**Name keys are case-insensitive, except for formulas** (decided 2026-10-14). Every place a paper
name is used as a key goes through **one** function, `metabolite_identity.name_key(name)`:
- the per-paper alias lookup;
- `name_lookup`;
- the dominant-form lock file;
- the identity-report review lists.

How it works:
- **Default:** case-fold, then collapse whitespace.
- **Exception, case-significant names:** a name that parses entirely as element symbols, with
  optional counts and an optional charge (`CO2`, `NH4+`, `Fe3+`, `Co`, `NO`), keeps its exact
  case, because there case is meaning (`Co` cobalt vs `CO` carbon monoxide, `NO` vs `No`, `HF` vs
  `Hf`).
  - The check uses the same periodic-table allow-list as `elements` (§3.3).
  - Acronyms such as `NAD`, `ADP`, `FAD`, `GABA` and `DMSP` don't parse (they contain a letter
    that isn't an element symbol), so they are case-folded normally.
  - A false positive such as `AMP` (Am + P) is harmless: it just isn't folded.
- **MNX resolver:** `resolve_metabolite` already lowercases every name in its name-match step
  (`_normalize_name`). For case-significant names that step is **skipped**, leaving only the
  exact alias / formula match.

Why: case-sensitivity currently splits one compound across papers. Eight compounds land on two
nodes depending on the capitalisation in the source table:
- glutamic acid (D-Glu / L-Glu);
- arginine and serine (L-form / `mnx:` placeholder);
- methionine, tryptophan, phenylalanine, glutamine;
- S-adenosyl-L-methionine (`chebi:142094` / `C00019`).

With this rule plus the B5 and dominant-form fixes, each converges on one node. Another 25
case-variant groups already agree. None of today's 211 paper names is a true formula, so the
exception is a guard. Moved rows appear in `metabolite_assay_repoints.tsv`.

**Dominant-stereo-form rule.** It applies when a paper name **doesn't specify stereo** and its
match (step 4b) is either the stereo-unspecified member of a stereo family, or several candidates
that all belong to one family. Example: "Glutamate" matching only `C00302` (69 catalysts) goes to
L-Glu `C00025` (1,905). This is the same reasoning as for paper-supplied ids: an unspecified
node doesn't assert a form.
- The name goes to the member with the most gene-reachable catalysing reactions in the build:
  L for amino acids, usually D for sugars.
- **Margin:** the winner needs ≥ 2× the runner-up's count. A smaller margin, a tie, or a family
  with no catalysis requires an override.
- **Frozen** on first computation into a committed lock file (§6). Later builds recompute and
  report `dominant_form_drift` (changed winner, or margin now < 2×); they never re-point silently.
- **Never applied** when the paper names a form ("D-Alanine") or an override exists.
- **Applies to paper-supplied ids as well** (explorer review §20). A direct id that is the
  *stereo-unspecified* member of a family (the target of `Metabolite_is_stereo_form_of`) asserts
  no form.
  - It goes through the same rule, margin and freeze.
  - When it moves, it gets `dominant_form_assumed` and is listed in
    `metabolite_assay_repoints.tsv`.
  - So a measurement lands on the same node whether the paper gave a name or an id.
- **Anomer guard.** The rule never picks a member that differs from the target only in anomeric
  configuration (α/β: β-D-GlcNAc for GlcNAc, α-D-glucose for D-glucose).
  - Anomers interconvert in solution and LC-MS doesn't resolve them, so the anomer-unspecified
    form is the correct target.
  - Anomer links still become stereo edges (§3.5); they're just never a rule destination.
  - Detection (settled in the plan): ChEBI's own anomer annotation where present; otherwise the
    names differ only by an `alpha-` / `beta-` prefix **and** the compound is a ring sugar or
    glycoside (by ChEBI class). A prefix alone is not enough: beta-Citronellol uses "β" as part of
    a trivial name. Tests: GlcNAc, MurNAc, glucose, and beta-Citronellol as a negative case.
  - Without the guard, D-glucose (284 catalysts) vs β-D-glucose (192) is refused only by the
    margin, by luck.
- Measured over 149 paper metabolites: 24 already on the dominant form, 118 have no stereo
  variant, **7 move**, all to L:
  - Asparagine `C16438` → `C00152`
  - D-Ala `C00133` → `C00041`
  - D-Glu `C00217` → `C00025`
  - D-Gln `C00819` → `C00064`
  - D-Trp `C00525` → `C00078`
  - Methionine `C01733` → `C00073`
  - Phenylalanine `C02057` → `C00079`

  The smallest margin is ≈ 5×.

  Over the 64 paper-direct ids: 2 are stereo-unspecified family targets, N-Acetyl-D-glucosamine
  `C00140` and N-Acetylmuramate `C02713`.
  - Both differ from their variants only anomerically, and both already carry the catalysts
    (65 vs 0; 98 vs 0), so **0 move**.
  - They stay `as_reported`: the paper's id named the D-form, and under the guard the anomer is
    not a stereo choice.

**Curated overrides added by this spec** (in the papers' `metabolite_aliases.yaml`):
- GSH → `C00051`; GSSG / "glutathione oxidized" → `C00127`.
- Cystine is separated from indole-3-acetic acid.
- The 4 `mnx:` placeholders: Arginine / Cysteine / Serine hydrochloride or compartment names →
  L-Arg / L-Cys / L-Ser, and "Cytosolic Ribose 5 Phosphate" → D-ribose 5-phosphate.
- The remaining ambiguous names → the form the paper names.
- **Same-entity check on every drop** (decided 2026-10-14). Before a measurement-only node is
  dropped (stage K), the plan checks whether it is actually the *same entity* as its target, e.g.
  SAM `chebi:142094` vs `C00019`. If so, the right outcome is an identity **merge** (its id goes
  into `all_identifiers`, basis via ChEBI xref or `manual`), not a drop with `recurated_from_ids`.
- Placeholder collision check (measured): no table that carries a placeholder also lists its
  re-curated L-form, so the mapping creates **0** same-table collisions. The compartment wording
  in "Chlorplastic Serine" / "Cytosolic Ribose 5 Phosphate" is dropped from the identity. The plan
  confirms the assay's own `compartment` field carries it.

**`stereo_assignment ∈ {as_reported, dominant_form_assumed}`** on assay edges records **what the
paper said**:
- `as_reported`: the paper named the form (a stereo-qualified name, or a direct id of a
  specified form), or the metabolite has no stereo variant **other than anomers**;
- `dominant_form_assumed`: an unqualified name (or a stereo-unspecified direct id) resolved to
  one member of a stereo family with more than one member (anomers aside). This applies
  **whatever chose it**: the rule, the lock, an override, **or a single exact match** (step 4a:
  "Lysine" matching only L-Lys is still an assumption). It also applies when the unspecified node
  itself wins. It includes the three amino-acid placeholders and every case-fold-converged row.
- **Anomer-only families** (GlcNAc `C00140`, MurNAc `C02713`): when the anomer guard leaves only
  the target, the row **keeps the target, `as_reported`**. Neither the "no catalysis" nor the
  "tie" clause fires.
- Unresolved rows carry no value; the property is dense on every assay **edge**, and unresolved
  rows produce no edge.

### 3.6a Curated files and pruning

Three curated or frozen files steer step 6:
- the per-paper `metabolite_aliases.yaml`;
- `config/metabolite_identity_overrides.yaml`;
- the dominant-form lock.

Pruning decides which compounds are written. **Every pruning decision takes these files into
account, and none of them can be silently ignored or silently go stale.**

| File | Effect on what is kept | Applied at | Stale entry |
|---|---|---|---|
| `metabolite_aliases.yaml` (per paper) | Every override target is a **paper candidate** at harvest (stage D), so it's included and enriched even if no reaction or transport reaches it. It's **kept** whenever a current row of that paper resolves through the override, and then it's a final target with the `metabolomics` tag. An override always wins over the global lookup and the dominant-form rule for its paper | stage D (candidate), J (resolution), K (kept as a target) | An entry whose name matches no row of its paper's tables is reported as `stale_alias_override`. Its target isn't kept on its account. A target id that can't be resolved at all (no KEGG / ChEBI / MNX entry) **fails** step 6 |
| `config/metabolite_identity_overrides.yaml` | **`merge:`** joins two nodes before any pruning (stage E): the KEGG-guard exception, `identity_basis: [manual]`, then the normal winner rule. **`keep_separate:`** blocks a merge from **every** basis, MNX included. **`canonical:`** (`{id, group, reason}`) names a merge group's canonical id over the margin rule (§3.1); an `id` outside its own `group` fails loading, two pins for one group fail step 6, and a pin whose group no longer forms is reported stale. Neither kind **adds** a node: identity overrides change identity, not reachability | stage E, before stage K and the final write | An entry naming an id that isn't in the pre-prune compound set is reported as `stale_identity_override` and skipped. A `merge:` + `keep_separate:` pair on the same two ids fails step 6 |
| `dominant_stereo_forms.json` (lock) | A frozen target is a **paper candidate** at stage D, like an alias target, so it's enriched and kept while its name still has rows | stages D, J | An entry whose name / direct id no longer occurs in any paper table is reported as `stale_lock_entry`, not deleted; deleting it is a human edit |

**Precedence when they meet:**
1. the paper's alias override;
2. the identity overrides, which decide which node an id belongs to (applied at stage E, before
   any row is resolved);
3. the row's direct id;
4. the lock;
5. the dominant-form rule.

A `keep_separate:` pair can't be undone by a later stage, and the stage-K prune never merges.

**Lock keys across builds.** Lock keys (canonical direct ids) and lock targets come from a
previous build. At stage J they are mapped through the current `id_aliases`. A target that was
merged away is remapped and reported (`lock_target_remapped`); a target that no longer exists is
reported as `stale_lock_entry`.

### 3.7 Compound classes

**New label `MetaboliteClass`**, not `BriteCategory`.
- `BriteCategory` is a gene ontology end to end: KO members, `gene_count` rollups, and the
  explorer's `ontology='brite'` / `brite_tree` / `genes_by_ontology`.
- All 8 KEGG compound trees ship, pruned to our metabolites plus ancestors:

| Tree | Our metabolites | Class nodes | Membership edges |
|---|---|---|---|
| br08001 Compounds with biological roles | 328 | 121 | 356 |
| br08002 Lipids | 517 | 163 | 536 |
| br08003 Phytochemical compounds | 143 | 76 | 143 |
| br08021 Glycosides | 38 | 9 | 38 |
| br08007 Pesticides | 33 | 31 | 33 |
| br08009 Natural toxins | 14 | 18 | 14 |
| br08006 Endocrine disrupting compounds | 8 | 9 | 8 |
| br08005 Bioactive peptides | 2 | 4 | 2 |
| **Union** | **867** | **≈ 431** | **≈ 1,130** |

- A compound listed under several classes gets one edge per class.
- Memberships on an alias KEGG id land on its canonical node.
- **Id prefix:** `kegg.brite:<br08xxx>…` shares the prefix with gene BRITE (`kegg.brite:<ko0xxxx>…`).
  Ids can't collide. **Confirmed by the explorer (review doc §37):** nothing routes on the
  `kegg.brite:` prefix; terms are matched by label.
- **`level_kind` naming:** `kegg_compound_class` / `_subclass` / `_group` / `_subgroup`, not the
  gene-BRITE `brite_class` / …. Deliberately distinct, so a `level_kind` filter never mixes the
  two layers (R1b).
- **Gene side: no change** to `annotation_types`, buckets or `has_any_edge`.

---

## 4. Data flow overview

**Inputs to step 6** (`build_kegg_metabolism_xrefs.py` + new `utils/metabolite_identity.py`):

| Input | Status |
|---|---|
| KEGG REST raw | existing |
| KEGG compound BRITE `br08xxx` JSON (8 trees) | **new** |
| TCDB raw TSVs | existing |
| MNX resolver | existing |
| `chebi.obo` | **new** |
| paperconfigs + per-paper `metabolite_aliases.yaml` | existing (alias format extended) |
| `config/metabolite_identity_overrides.yaml` | **new** |

**Outputs → consumers:**

| Produced by | File | Status | Read by | Graph content |
|---|---|---|---|---|
| step 6 | `kegg/kegg_data.json` | changed | `metabolism_adapter` (changed) | `Metabolite`, `Reaction`, reaction / pathway edges, stereo-form + enantiomer edges |
| step 6 | `kegg/brite_br08xxx.json` | **new** | `metabolite_class_adapter` (**new**) | `MetaboliteClass` + membership / is-a edges |
| step 6 | `tcdb/tcdb_pruned.json` | changed | `tcdb_adapter` (changed) | `TcdbFamily` (+ `substrate_breadth`), substrate edges, `ChemicalRole` + role edges |
| step 6 | `chebi/chebi_metabolite_xrefs.json` | **new** | step 6 only (input cache) | — |
| step 6 | `metabolomics/metabolite_id_mapping.json` | changed | nobody at build time (the record of step 6's paper decisions; diagnostics) | — |
| step 6 | `metabolomics/dominant_stereo_forms.json` | **new** (lock file) | step 6 on later builds | — |
| step 6 | `metabolomics/metabolite_identity_report.json` | **new** | nobody at build time (diagnostics, release TSVs) | — |
| **step 6** (was step 7) | `<stem>_resolved.csv` | changed | `metabolite_assay_adapter` (changed) | assay edges (+ `stereo_assignment`) |

All adapters → BioCypher CSVs → `neo4j-admin import` → **`post-import.sh` / `.cypher` (changed:
new rollups and indexes, §9)**.

- Every identity, stereo, role and class decision is made in **step 6**.
- Step 6 also writes the paper `_resolved.csv` files, from the same in-memory decisions, so the
  graph's metabolite set and the paper tables can never disagree. Step 7 is folded into step 6
  (§5.3).
- Adapters only copy fields.
- Post-import only computes rollups and indexes.

---

## 5. Prepare-data steps

### 5.1 Inputs

| Input | Status | Location | Fetch / outage behaviour |
|---|---|---|---|
| ChEBI ontology `chebi.obo.gz` (~47 MB, release 256) | **new** | `cache/data/chebi/raw/` (gitignored) | step 6, fetched when the derived cache is missing or with `--refetch-raw`. On failure with `chebi_metabolite_xrefs.json` present: warn and reuse; missing ⇒ fatal (TCDB precedent) |
| KEGG compound BRITE `br:br08001…br08021` JSON (8 trees) | **new** | `cache/data/kegg/raw/` (gitignored) | step 6, same fetch rule. Fail-loud parse (`_require_parsed`); on outage, reuse the committed `brite_br08xxx.json` |
| `config/metabolite_identity_overrides.yaml` | **new**, curated, starts empty | repo | read by step 6 |
| Per-paper `metabolite_aliases.yaml` (`aliases_file` in paperconfig) | **format extended** (§6.8) | paper directories | read by step 6 |
| KEGG REST raw, TCDB raw TSVs, MNX resolver, paperconfigs | unchanged | — | — |

### 5.2 Step 6 (`build_kegg_metabolism_xrefs.py`): ordered pipeline

New logic lives in **`multiomics_kg/utils/metabolite_identity.py`** (pure functions) and
**`multiomics_kg/utils/chebi_utils.py`** (obo parse → derived cache). Step 6 calls them in this
order:

Stages are lettered (A–P) so they can't be confused with the §6.x file sections. Each stage
uses only what earlier stages produced, and **nothing is written to disk until stage P**.

| Stage | What it does | Existing / new | Output in memory |
|---|---|---|---|
| **A** | **Load ChEBI** → build `chebi_metabolite_xrefs.json` in memory (§6.4 defines coverage). It's **rebuilt whenever the raw obo is present**; the committed file is reused only when the raw download fails | new | ChEBI lookup |
| **B** | Existing TCDB hierarchy, gene-reachable sets, rollup (`_compute_subtree_substrate_strings`), prune | existing | substrate strings per kept TC id |
| **C** | **Partition role strings out** of the substrate strings (§3.4), before MNX resolution; map each to its `chebi.role:` id; build the pruned role hierarchy (`chemical_roles`, incl. each role's `all_identifiers`) and `subtree_roles`. **Depth isn't computed here**; the adapter does it (§8.2) | new, in `_resolve_substrate_strings` | role blocks |
| **D** | **Resolve substrate strings and harvest paper candidates.** Secondary ChEBI → primary first. **ChEBI → KEGG re-key for every candidate source**: the TCDB path (`_primary_for` in `_resolve_substrate_strings`) **and** the paper path (`metabolite_utils.mnxm_to_primary_id` / the harvest), alias-override targets and lock targets. When MNX gives no KEGG id, ChEBI's own KEGG xref does. Every KEGG id therefore exists before the existing steps run, **in their current order, unchanged**: transport-reachable pathway extension, KEGG pruning + MNX enrichment (`build_pruned_kegg_data` / `_bulk_enrich_compounds`), `evidence_sources` tagging, the folds. The paper harvest (`_harvest_paper_metabolites`) **collects candidates, not decisions**: per paper row, the direct id, the per-paper override target, the lock target and all MNX / name candidates, all enriched like any compound. Its global cross-paper name dedup is removed (§3.6) | changed | enriched compound entries, plus per-row paper candidates |
| **E** | **Canonicalise**: union-find over the merge links, the KEGG guard, the MNX stereo guard, the margin winner rule (unresolved → step fails), and `metabolite_identity_overrides.yaml` (`merge:` joins, `keep_separate:` blocks every basis; stale entries reported; §3.6a). After this, **no KEGG id is introduced** that stage D didn't enrich; any that would be fails the step | new | groups + alias → canonical map (`id_aliases`) |
| **F** | **Rewrite every id to canonical**: `reactions[*].compounds` (deduplicated), compound `pathways` / `evidence_sources` (union), `subtree_substrates`, paper candidates, alias-override and lock keys / targets. Lock remaps are reported (§3.6a). Done **here**, so stages I–K count reactions and transport on canonical nodes | new | canonical graph data |
| **G** | **Structure**: `structure_kind`, the B2 scrub (both conditions), generalized-formula fill, no `mass` on a generic formula, `elements` kept (§3.3) | new | scrubbed compound entries |
| **H** | **Cross-reference fields**: `chebi_id` (ChEBI KEGG xref, tie-break), `chebi_ids`, `kegg_compound_ids`, `mnxm_ids`, `identity_basis` (merged **and** re-keyed nodes), `name`, and `all_identifiers`, built **after** the scrub, so a scrubbed class node never advertises its monomer's InChIKey (§3.2) | new; replaces the lowest-alias pick | final identity fields |
| **I** | **Stereo families**: connected components of the pairwise links, on scrubbed canonical nodes (§3.5). In memory only | new | stereo families |
| **J** | **Resolve every paper row** in the §3.6 order: paper override → direct id → lock → match → dominant-form rule (margin, anomer guard). Compute new lock entries and `dominant_form_drift`; set `stereo_assignment`; collect ambiguous and unmatched names (§3.6) | new | final per-row decisions |
| **K** | **Measured set = the final targets** (§3.6a). The `metabolomics` tag goes only on final targets. A candidate that isn't a final target **and** has no catalysis or transport reason to exist is dropped (≈ 11 today), after the same-entity check (§3.6). Each drop goes to the remap (`remeasured_elsewhere` / `removed_placeholder`) and to the target's `recurated_from_ids`, never to `all_identifiers` | new | final compound set |
| **L** | **Stereo links on surviving nodes only**: `stereo_form_of` / `enantiomer_of` lists between compounds that survive K (e.g. L-Trp ↔ D-Trp is discarded once D-Trp is dropped). Final counts go to the report | new | final stereo lists |
| **M** | **Compound classes**: map the 8 BRITE trees' members to canonical ids; prune to written metabolites + ancestors (§3.7) | new | class trees |
| **N** | TCDB Pfam / GO bridges | existing, unchanged | |
| **O** | **Lumping pin**: load `lumping_pin` (on the first run, **seed it from the literal 13**, §3.4a); recompute the threshold on canonical `subtree_substrates`; report `lumping_set_drift` for frozen ids; classify and append newly kept ids (`lumping_set_added`) | new | pin |
| **P** | **Validate, then write everything at once.** Checks that fail the step: an ambiguous paper name with no override; a disallowed same-table collision; an `all_identifiers` entry or InChIKey on two nodes (across `Metabolite` ∪ `ChemicalRole`); a dangling stereo / enantiomer / class / role target; a canonical node that needs a ChEBI entry the cache lacks; a resolved row without `stereo_assignment` (unresolved rows excepted); a merge group the canonical margin rule can't decide (`canonical_margin_unresolved`); a written node holding two stereo-distinct KEGG ids outside a `manual` merge. **Only if every check passes**, write: `kegg_data.json`, `tcdb_pruned.json`, `brite_br08xxx.json`, `chebi_metabolite_xrefs.json`, `metabolite_id_mapping.json`, the lock, the identity report, and every paper `<stem>_resolved.csv` (§5.3). A failed run writes nothing, so the committed lock and pin never hold half a build | changed | files |

Determinism: the same inputs give byte-identical outputs (sorted keys and lists), apart from
`frozen_on` dates on newly frozen entries.

### 5.3 Paper `_resolved.csv` files: now written by step 6 (step 7 folded in)

**Why.** Today step 6 decides the paper metabolites and step 7 writes the resolved tables. With
ids now moving (merges, re-keys, the dominant form), re-running step 6 alone would leave stale
`_resolved.csv` files that point at old ids, giving dangling assay edges that `neo4j-admin`
drops silently. One step that both decides and writes can't drift.

**What step 6 writes**, at stage P, from the in-memory decisions of stages J / K, one
`<stem>_resolved.csv` per `metabolite_assays_table`:
- `metabolite_id`: canonical, and always a compound that step 6 also wrote to `kegg_data.json`;
- `resolution_method`: `name_match`, `alias_override`, `dominant_stereo_form`, `kegg_direct`, …;
  the `ambiguous_multi_id` method is retired;
- **new** `stereo_assignment` (§3.6). It's empty on unresolved rows, which produce no edge.

Unresolved rows stay in the file with an empty `metabolite_id` and `resolution_method =
unresolved`, as today.

The row-level rules are unchanged from §3.6:
- `name_key` lookups (paper overrides first);
- direct ids canonicalised;
- no ambiguous pick.

**Step 6 fails loudly** (stage P, before writing anything) when:
- a name is ambiguous with no override;
- two distinct **`name_key`s** of one table resolve to one id, unless both have alias entries
  pointing to it (§3.6). Case variants of one name are one key, not a collision;
- a **resolved** row has no `stereo_assignment`.

**Unmatched names don't fail the step**: they're unresolved and reported (§3.6).

**Step 7** (`resolve_paper_metabolites.py`) is **retired**. Its row logic moves into step 6
(`metabolite_identity.write_resolved_tables`). `prepare_data.sh --steps 7` stays accepted as an
alias that prints "folded into step 6" and runs nothing, so existing commands don't break. The
`--papers` convenience goes away: re-resolving a paper means re-running step 6 (`--force`, about a
minute, no network). `metabolite_id_mapping.json` stays as a written record of step 6's paper
decisions, but no step reads it any more.

### 5.4 Other steps

Steps 0–5, 8, 9 and 10 are unchanged. Step 2's `gene_annotations_merged.json` is untouched:
metabolite identity is not a gene-annotation field.

**Re-run after merge:** `prepare_data.sh --steps 6 --force`, then the Docker rebuild.

---

## 6. Cached and resolved files

All files under `cache/data/` below are **committed**, unless marked gitignored.

### 6.1 `cache/data/kegg/kegg_data.json`: changed

- `compounds` (keyed by bare `C…`) and `additional_compounds` (keyed by full CURIE) hold
  **canonical entries only**. Merged-away keys disappear.
- Removed role entries and `mnx:` placeholders disappear.
- New or changed fields per entry:

```json
"C00244": {
  "name": "Nitrate",
  "kegg_compound_ids": ["C00244"],
  "chebi_id": "48107", "chebi_ids": ["48107", "17632"],
  "mnxm_id": "MNXM732399", "mnxm_ids": ["MNXM732399", "MNXM732398"],
  "hmdb_id": "HMDB0002878",
  "all_identifiers": ["chebi:14654", "chebi:17632", "chebi:25545", "chebi:48107",
                      "hmdb:HMDB0002878", "inchikey:GRYLNZFGIOXLOG-UHFFFAOYSA-N",
                      "kegg.compound:C00244", "mnx:MNXM732398", "mnx:MNXM732399"],
  "identity_basis": ["chebi_conjugate"],
  "structure_kind": "defined",
  "formula": "HNO3", "mass": 63.012, "inchikey": "GRYLNZFGIOXLOG-UHFFFAOYSA-N", "smiles": "…",
  "evidence_sources": ["metabolism", "metabolomics", "transport"],
  "pathways": ["ko00910", "…"]
}
```

The values are illustrative; the plan fills the real ones. Sparse keys such as `stereo_form_of`
are omitted in the example because they'd be empty here.

- `recurated_from_ids` (sparse) on the target of each dropped placeholder or measurement-only node.
- **Sparse keys are omitted when empty:** `mnxm_ids`, `identity_basis`, `recurated_from_ids`, `stereo_form_of`,
  `enantiomer_of`, and
  `formula` / `mass` / `inchikey` / `smiles` where §3.3 drops them.
- `stereo_form_of` values are **full node ids**, because a target may be a `chebi:` node.
- `enantiomer_of` (full node ids) is written **only on the member whose canonical id sorts
  lower**, so each pair appears once in the file and once in the graph.
- `reactions[*].compounds`: canonical KEGG ids, deduplicated.

### 6.2 `cache/data/kegg/brite_br08xxx.json` (8 files): new

A **new, already-pruned** shape. Unlike the existing `brite_ko0xxxx.json`, which are raw KEGG
`{name, children}` trees pruned by `brite_adapter` at build time, these are pruned in step 6,
because membership must be mapped to canonical metabolite ids. Per tree:

```json
{
  "tree_code": "br08001", "tree": "compounds_with_biological_roles",
  "nodes": [{"id": "kegg.brite:br08001.A1", "parent": null, "name": "Peptides", "level": 0}, …],
  "members": [{"metabolite_id": "kegg.compound:C00041", "class_id": "kegg.brite:br08001.A1.B1"}, …]
}
```

- `metabolite_id` values are canonical.
- Raw downloads go to `cache/data/kegg/raw/` (gitignored). On a KEGG outage step 6 reuses the
  committed pruned files, re-mapping their `metabolite_id`s through the current `id_aliases`.

### 6.3 `cache/data/tcdb/tcdb_pruned.json`: changed

- `subtree_substrates`: canonical ids, deduplicated, **role primaries removed**.
- **New block** `chemical_roles: {"chebi.role:33229": {"name": "vitamin (role)", "parents":
  ["chebi.role:…"], "level": 2, "level_kind": "chebi_role", "all_identifiers":
  ["chebi.role:33229", "chebi:33229", "chebi:10004", …]}, …}`: the pruned role hierarchy, i.e.
  the 57 TCDB roles plus `is_a` ancestors up to the three branches (no `role` root). Every
  parent is itself an entry.
- **New block** `subtree_roles: {tc_id: ["chebi.role:…", …]}`, the role counterpart of
  `subtree_substrates`: sorted, kept TC ids only, a key only where non-empty.
  Depth per pair is computed by the adapter from these two blocks and the kept hierarchy, the
  same way as for substrates.
- **New block** `lumping_pin: {"frozen_on": "…", "frozen_kept_ids": [tc_id, …],
  "lumping_tcdb_ids": [tc_id, …]}` (§3.4a):
  - seeded from the literal 13;
  - frozen ids change only by a human edit;
  - ids that become kept after `frozen_on` are classified by the threshold and appended
    automatically (and reported). `frozen_kept_ids` is not extended, so the two cases stay
    distinguishable.
- `kept_tcdb_ids`, `seed_aliases`, `pfam_bridge`, `go_bridge`: unchanged.

### 6.4 `cache/data/chebi/chebi_metabolite_xrefs.json`: new (step-6 input cache)

`{"release": "256", "terms": {"CHEBI:17632": {name, charge, alt_ids, kegg, conjugates,
tautomers, enantiomers, is_a, is_role, is_class, generalized_formula}, …}}`

- **Coverage** (decided 2026-10-14), the union of:
  - every non-obsolete ChEBI term with a KEGG cross-reference (needed for `chebi_id` and
    re-keying in the KEGG → ChEBI direction);
  - every ChEBI id seen in a TCDB substrate, a paper candidate, an MNX candidate, or a
    curated / lock target;
  - their `is_a` ancestors and their conjugate / tautomer / enantiomer partners (with the
    partners' own KEGG xrefs).
- **Freshness:**
  - rebuilt by stage A whenever the raw obo is present;
  - reused as committed only when the raw download fails;
  - in that case, ids missing from the cache are reported, and step 6 fails if a canonical
    node needs one.

  No silent corpus-coupling (the InterPro `description` lesson).
- **Read by step 6 only**, never by an adapter. Written at stage P with everything else.
- Raw `chebi.obo.gz` goes to `cache/data/chebi/raw/` (gitignored).

### 6.5 `cache/data/metabolomics/metabolite_id_mapping.json`: changed

A **written record of step 6's paper decisions** (no step reads it any more; see §5.3).
`schema_version` goes to **2**. Shape:
- `id_aliases: {any id CURIE: canonical id}`: the identity map from stage E.
- `name_lookup: {name_key: canonical id}`: the global matches, single-valued (ambiguous names
  absent).
- `paper_overrides: {paper: {name_key: {id, stereo_assignment, reason}}}`: each paper's alias
  yaml, normalised.
- **`decisions: {paper: {table: {name_key or direct id: {metabolite_id, resolution_method,
  stereo_assignment}}}}`**: the per-paper, per-table outcome. This is where "two papers map one
  name differently" lives, so there is no global `stereo_assignment` / `resolution_method` map.
- `unresolved: {paper: {table: [name_key, …]}}`; `ambiguous` is always empty in a successful
  write, because ambiguity fails the step.
- `compounds`, `conflicts`, `per_paper_stats`: kept, on canonical ids.
- The legacy `specific_lookup` / `multi_lookup` keys are **dropped**: `multi_lookup` was the
  ambiguity store, and `specific_lookup` is replaced by `name_lookup`.

Methods: `name_match` | `alias_override` | `dominant_stereo_form` | `kegg_direct` / `chebi_direct`
/ `mnx_direct` | `unresolved` (`ambiguous_multi_id` is retired).

### 6.6 `cache/data/metabolomics/dominant_stereo_forms.json`: new (lock file)

```json
{"names":      {"alanine": {"id": "kegg.compound:C00041", "frozen_on": "2026-10-xx",
                            "counts": {"winner": 520, "runner_up": 103},
                            "runner_up_id": "kegg.compound:C00133"}},
 "direct_ids": {"kegg.compound:Cxxxxx": {"id": "kegg.compound:Cyyyyy", "frozen_on": "…",
                                         "counts": {…}, "runner_up_id": "…"}}}
```

- Keys are `name_key`s for names, and canonical ids for paper-supplied stereo-unspecified ids
  (explorer C1 covers both). The `direct_ids` example is hypothetical: today 0 direct ids move.
- Keys and targets are re-mapped through the current `id_aliases` on every build (§3.6a).
- New entries are computed at stage J and written at stage P, only on a successful run.
- Read on every later build.
- Edited only by hand:
  - delete an entry to re-freeze;
  - or replace it with an alias override for that paper's name. For a direct id, the override
    goes on that row's name, since overrides win over direct ids (§3.6).
- It plays the role of the "curated override table" of explorer C1, but as a lock file, so the
  build never writes into curated paper files.

### 6.7 `cache/data/metabolomics/metabolite_identity_report.json`: new (diagnostics only)

Merge groups by basis, alias → canonical map, rejected ChEBI-only KEGG↔KEGG links, refused
assay-collision merges, `stale_alias_override` / `stale_identity_override` / `stale_lock_entry`, `roles_removed`, removed placeholders, heuristic name candidates,
`ambiguous_paper_names` (with candidates), `dominant_form_drift`, counts.
- **Never read at build time.**
- The release TSVs (§12) are generated from it.

### 6.8 Per-paper `metabolite_aliases.yaml`: format extended

```yaml
Glutathione:            {id: "kegg.compound:C00051", stereo_assignment: as_reported, reason: "GSH, reduced"}
Arginine hydrochloride: {id: "kegg.compound:C00062", stereo_assignment: dominant_form_assumed, reason: "salt form; paper gives no stereo"}
```

- The old `name: id` string form stays valid **only** when the target has no stereo variant; it
  is then `as_reported`.
- Otherwise step 6 fails and asks for an explicit `stereo_assignment`.
- This spec adds the curated entries listed in §3.6.

### 6.9 `<stem>_resolved.csv` (paper directories): changed

- Columns: `metabolite_id` (canonical), `resolution_method` (new method values, §6.5), **new**
  `stereo_assignment`.
- Rows that change:
  - the 36 re-resolved names, including the GSH split;
  - the case-fold convergences (§3.6);
  - canonicalised paper-supplied ids (merged-away or re-keyed);
  - unqualified exact matches moved by the dominant-form rule (step 4b).

  Every change is listed in `metabolite_assay_repoints.tsv` with its reason.

### 6.10 `config/metabolite_identity_overrides.yaml`: new (curated)

`merge:` / `keep_separate:` lists of `{a, b, reason}`. It starts empty, and the 16 heuristic
candidates are listed in the report for review.

---

## 7. Schema changes (`config/schema_config.yaml`)

### 7.1 Nodes

| Node | Change |
|---|---|
| `metabolite` | **add:** `kegg_compound_ids: str[]`, `chebi_ids: str[]`, `mnxm_ids: str[]`, `all_identifiers: str[]` (described as "CURIEs, unlike `Gene.all_identifiers`"), `recurated_from_ids: str[]` (sparse, no identity meaning), `identity_basis: str[]`, `structure_kind: str`, `kegg_functional_categories: str[]` (post-import; named `kegg_compound_categories` until 2026-10-15), `kegg_brite_classes: str[]` (post-import). Existing properties keep their names and types; their values change per §3.2 / §3.3 |
| `tcdb family` | **add:** `substrate_breadth: str` (dense, `lumping` \| `specific`, §3.4a) |
| **`chemical role`** (new) | `preferred_id: chebi.role` (house-minted); properties `name`, `level: int`, `level_kind`; post-import `transporter_count` / `direct_transporter_count`, `tcdb_family_count` / `direct_tcdb_family_count` (§9.1) |
| **`metabolite class`** (new) | `preferred_id: kegg.brite`; properties `name`, `tree`, `tree_code`, `level: int`, `level_kind`; post-import `metabolite_count`, `direct_metabolite_count`, `organism_count` |
| `data source` | no shape change; new `kegg` and `chebi` instances (§8.5) |

### 7.2 Edges

| Edge | Change |
|---|---|
| **`metabolite is stereo form of metabolite`** (new) | Metabolite → Metabolite, `label_as_edge: metabolite_is_stereo_form_of`, `sources: str[]`, `evidence: str` |
| **`metabolite is enantiomer of metabolite`** (new) | Metabolite → Metabolite (lower → higher id; symmetric), `label_as_edge: metabolite_is_enantiomer_of`, `sources: str[]`, `evidence: str` |
| **`metabolite in metabolite class`** (new) | Metabolite → MetaboliteClass, `label_as_edge: metabolite_in_metabolite_class`, `sources: str[]`, `evidence: str` |
| **`metabolite class is a metabolite class`** (new) | child → parent, `label_as_edge: metabolite_class_is_a_metabolite_class`, no properties |
| `assay quantifies metabolite`, `assay flags metabolite` | **add** `stereo_assignment: str` |
| **`chemical role is a chemical role`** (new) | child → parent, `label_as_edge: chemical_role_is_a_chemical_role`, no properties |
| **`tcdb family transports chemical role`** (new) | TcdbFamily → ChemicalRole, `label_as_edge: tcdb_family_transports_chemical_role`, `substrate_depth: str` |

### 7.3 Controlled vocabularies (`config/controlled_vocabularies.yaml`)

New entries:

| Entry | Values |
|---|---|
| `Metabolite.identity_basis` | closed, sparse |
| `Metabolite.structure_kind` | closed: `defined`, `generic`, `polymer` |
| `Metabolite.kegg_functional_categories` (renamed from `kegg_compound_categories`, 2026-10-15) | closed: the br08001 level-0 names |
| `Metabolite.kegg_brite_classes` (2026-10-15) | closed: the 28 `"<Tree>: <class>"` level-0 classes of the 7 other compound trees; `value_descriptions` name the tree code |
| `MetaboliteClass.level_kind` | closed: `kegg_compound_class`, `_subclass`, `_group`, `_subgroup` |
| `MetaboliteClass.tree` | closed: the 8 tree names |
| `ChemicalRole.level_kind` | closed: `chebi_role_branch`, `chebi_role` |
| `Tcdb_family_transports_chemical_role.substrate_depth` | closed: `most_specific`, `inherited` (same values as the metabolite edge) |
| `TcdbFamily.substrate_breadth` | closed: `lumping`, `specific`, `not_assessed` (pinned set, §3.4a) |
| `Metabolite_in_metabolite_class.sources` / `.evidence` | `[kegg]` / `[curated]` |
| `Metabolite_is_stereo_form_of.sources` / `.evidence` | `[chebi]` / `[curated]` |
| `Metabolite_is_enantiomer_of.sources` / `.evidence` | `[chebi]` / `[curated]` |
| `Assay_quantifies_metabolite.stereo_assignment`, `Assay_flags_metabolite.stereo_assignment` | closed: `as_reported`, `dominant_form_assumed` |

**Trust-style vocabularies need `value_descriptions`** (the 2026-08-29 rule: a closed
vocabulary describes all its values or none). That covers:
- `identity_basis`, `structure_kind`, `substrate_breadth`, `stereo_assignment`;
- the `sources` / `evidence` entries of the three new edge types;
- the `sources` values `kegg` / `chebi`.

The unit test that pins `sources` text to the DataSource registry is extended to `kegg` / `chebi`.
The remap TSV's `basis` and repoint `reason` values are release-file columns, not graph
properties. They are defined in §12, not registered as vocabularies. The vocabulary hash
changes.

---

## 8. Adapter changes

### 8.1 `metabolism_adapter.py`: changed

- **`get_nodes`**: copies the new compound fields (§6.1) onto `Metabolite`.
  - Sparse fields are dropped when absent (`_drop_nulls`).
  - `_clean_str` is applied to every string and list element.
  - `_parse_elements` gains the placeholder-`R` handling and the periodic-table allow-list (§3.3).
- **New `_metabolite_stereo_edges`**:
  - iterates `compounds` and `additional_compounds`, and for each `stereo_form_of` target emits
    `metabolite_is_stereo_form_of`, id `m2s:{src}:{tgt}`, `{sources: ['chebi'], evidence:
    'curated'}`;
  - skips and logs any target not among its emitted node ids (a guard; step 6 already
    guarantees this).
- **New `_metabolite_enantiomer_edges`**: the same pattern over `enantiomer_of`, emitting
  `metabolite_is_enantiomer_of`, id `m2e:{lower}:{higher}`. It asserts source < target, so a
  reverse duplicate can't be emitted.
- `_reaction_metabolite_edges`, `_metabolite_pathway_edges`: no code change; they read canonical
  ids.
- The `get_nodes` test-mode cap stays as is.

### 8.2 `tcdb_adapter.py`: changed

- **New:** yields:
  - `chemical role` nodes from `chemical_roles` (`name` via `_clean_str`, `level`, `level_kind`,
    `all_identifiers`);
  - `chemical_role_is_a_chemical_role` edges from each entry's `parents`;
  - `tcdb_family_transports_chemical_role` edges from `subtree_roles`, id
    `t2r:{tc_id}:{role_id}`, with `substrate_depth`.

  Dangling-proof by construction: both blocks come from one file.
- **Refactor:** `MultiTcdbAnnotationAdapter._compute_substrate_depth` takes the `{tc_id: ids}` map
  as an argument; today it is hard-wired to `self._subtree_substrates`. It is then called once for
  substrates and once for roles.
- `TcdbFamily` nodes gain `substrate_breadth`:
  - `lumping` when the id is in `lumping_pin.lumping_tcdb_ids`;
  - `not_assessed` when `level < 2`;
  - otherwise `specific`.
- The substrate-edge code is unchanged; it simply no longer receives role primaries.

### 8.3 `metabolite_assay_adapter.py`: changed

- Reads the `stereo_assignment` column of `<stem>_resolved.csv` and copies it onto
  `Assay_quantifies_metabolite` / `Assay_flags_metabolite`.
- A resolved row with an empty value is an error (step 6 guarantees one).

### 8.4 `metabolite_class_adapter.py`: new

- Reads the 8 `brite_br08xxx.json` files.
- Yields `metabolite class` nodes (`name` cleaned of `[Fig]` / HTML, `tree`, `tree_code`,
  `level`, `level_kind`).
- Yields `metabolite_class_is_a_metabolite_class` edges, and `metabolite_in_metabolite_class`
  edges with `{sources: ['kegg'], evidence: 'curated'}`.
- **Dangling-proofed** by an injected `metabolite_node_ids` set taken from the metabolism adapter
  (BRITE `known_ko_ids` precedent).
- Follows the house adapter interface; it doesn't need a Multi* wrapper, since the data is
  global.

### 8.5 `data_source_adapter.py` + config: changed

- Add a `kegg` DataSource (`data_source:kegg`, scope global, provenance reference_db,
  `info_types: [compound_classes]`).
- Add `data_source:chebi` the same way, so the stereo and enantiomer edges' `sources` join
  (R2). ChEBI is CC BY 4.0, so its `description` carries the attribution, and `Schema_info`
  records the release.
- Whether these join `gene_annotations_config.yaml` `logical_sources` or a small non-gene
  registry beside it is decided in the plan. The `sources` text-pinning test must keep passing.

### 8.6 `create_knowledge_graph.py`: changed

- Instantiate `MetaboliteClassAdapter` **after** `MultiMetabolismAdapter`, injecting the emitted
  metabolite ids.
- Pass the ChEBI release into the `Schema_info` provenance.

---

## 9. Post-import changes (`scripts/post-import.sh` and `scripts/post-import.cypher`, kept identical)

### 9.1 New computations

| Computation | Rule |
|---|---|
| `MetaboliteClass.metabolite_count` | distinct metabolites in the **subtree** (`Metabolite_class_is_a_metabolite_class*0..`) |
| `MetaboliteClass.direct_metabolite_count` | distinct metabolites with a membership edge to this node |
| `MetaboliteClass.organism_count` | distinct organisms with `Organism_has_metabolite` to a subtree metabolite |
| `Metabolite.kegg_functional_categories` | sorted distinct br08001 level-0 names reachable from the metabolite; **sparse** (renamed from `kegg_compound_categories`, 2026-10-15) |
| `Metabolite.kegg_brite_classes` | sorted distinct `"<Tree title>: <class>"` over the level-0 classes of every tree except br08001 reachable from the metabolite (tree title = `tree` with `_` → space, first letter upper case; class = KEGG's name without its 2-letter LIPID MAPS code); **sparse** (explorer ask 6, 2026-10-15) |
| `ChemicalRole.transporter_count` / `direct_transporter_count` | distinct TcdbFamily sources with a `most_specific` role edge, the `Metabolite.transporter_count` rule. Subtree version over `Chemical_role_is_a_chemical_role*0..`, counting **distinct** families (the role hierarchy is a DAG; house convention for hierarchical labels) / direct version on this node only |
| `ChemicalRole.tcdb_family_count` / `direct_tcdb_family_count` | the same pair, at any depth |

### 9.2 Indexes

**New:**
- `metabolite_class_tree_idx`, `metabolite_class_level_idx`, `metabolite_class_level_kind_idx`;
- full-text `metaboliteClassFullText` on `name`;
- `chemical_role_id_idx`, `chemical_role_level_idx`, `chemical_role_level_kind_idx` (house
  convention for ontology labels);
- `metabolite_structure_kind_idx` (the explorer filters on it);
- `metabolite_class_id_idx`, for parity with `chemical_role_id_idx`.

**Drop and recreate:**
- `metaboliteFullText` → on `name`, `formula`, `all_identifiers`. Note: the live index covers
  **`name` only** (`post-import.sh`), although CLAUDE.md claims "name, formula". So `formula` is
  an addition too. CLAUDE.md is corrected, and the CHANGELOG says so. A KG test checks that a bare
  `C01353` / `MNXM732398` / `HMDB0002878` query still hits; if the analyser doesn't produce bare
  tokens, the plan adds a search-only property of bare local parts to the index instead.
- (`tcdbFamilyFullText` unchanged.) **New** full-text `chemicalRoleFullText` on `ChemicalRole.name`
  (its scalar indexes are listed above), so `search_ontology`-style lookups for "siderophore" /
  "drug" find the role, and through it the families.

Existing scalar indexes (`metabolite_kegg_id_idx`, `metabolite_chebi_idx`, …) are unchanged.

### 9.3 Unchanged code, changed values

These recompute automatically on the canonical graph:
- Metabolite rollups (`organism_count` / `_names`, `catalyst_gene_count`,
  `transporter_gene_count`, `transporter_count`, `pathway_*`, `measured_*`);
- `Organism_has_metabolite`;
- `Gene.transported_metabolite_count` and `Gene.catalyzed_metabolite_count`;
- `OrganismTaxon.*_metabolite_count`;
- `TcdbFamily.metabolite_count`;
- `Gene.transport_substrate_resolution`: **one code change.** The inline `coalesce(t.level, 0) >=
  2 AND coalesce(t.metabolite_count, 0) >= 50` breadth test becomes
  `t.substrate_breadth = 'lumping'` (§3.4a), and the R3 comment is updated to cite the
  exception. Values are unchanged vs. the pre-change graph;
- MetaboliteAssay / Experiment / Publication metabolite rollups.

No post-import statement depends on metabolite id prefixes.

---

## 10. Expected impact (re-measured by the plan)

| Area | Change |
|---|---|
| `Metabolite` nodes | 3,356 → ≈ 3,214 (−74 merged, −57 roles, ≈ −11 measurement-only incl. the 4 placeholders; 89 re-keyed) |
| `evidence_sources` | the `metabolomics` tag moves to the final targets; ≈ 7 nodes keep their other evidence but lose the tag (D-Ala, D-Glu, D-Gln, Asn `C16438`, Met `C01733`, Phe `C02057`, …) |
| Node ids | ≥ 173 secondary `chebi:` ids remapped; 89 re-keyed + 74 merged away; ≈ 68 removed (57 roles + ≈ 11 measurement-only incl. the 4 placeholders) |
| `chebi_id` values | change on ≥ 1,785 nodes |
| Structure | 20 class nodes lose a representative structure but keep `elements`; ≈ 49 gain a generic formula and `elements` |
| Paper metabolites | 36 re-resolved (7 by the rule). The GSH / GSSG duplicate edges split |
| New edges | ≈ 60 stereo-form · ≈ 25 enantiomer · ≈ 1,130 class memberships · ≈ 420 class is-a |
| New nodes | ≈ 431 `MetaboliteClass` |
| Moved edges | 392 role `Tcdb_family_transports_metabolite` edges become `Tcdb_family_transports_chemical_role` edges (same pairs, same depth rule); merged metabolite targets deduplicate |
| New role layer | ≈ 111 `ChemicalRole` nodes (the 57 TCDB roles, replacing the 57 role `Metabolite` nodes, plus ancestors up to the 3 branches) and their `Chemical_role_is_a_chemical_role` edges |
| Counts | metabolite / organism / gene metabolite counts shift slightly. The lumping set and `transport_substrate_resolution` are **unchanged** (pinned) |
| `Reaction_has_metabolite` | small drop where merged compounds deduplicate on a reaction (e.g. HCO₃⁻ + carbonic acid) |
| KEGG pathways | the 89 re-keyed compounds bring their compound → pathway links. That can **add pathway-level `KeggTerm` nodes** (transport-reachable, Option B) with no gene attached. Measured by the plan; gene annotation state must not move |

---

## 11. Testing

**Unit tests** (`tests/test_metabolite_identity.py`, `tests/test_chebi_utils.py`, fixtures, no
network).

*Identity:*
- nitrate merges (`chebi_conjugate`);
- HCO₃⁻ merges by MNX, with the winner chosen by reaction count;
- a ChEBI-only KEGG↔KEGG link is reported, not applied; water / hydroxide stay two nodes;
- an assay-collision merge is refused;
- glucose and galactose stay apart;
- secondary → primary mapping; an isotopologue is never `chebi_id`.

*Identifiers:*
- `all_identifiers` holds CURIE forms only, each unique across nodes, and every merged-away or
  re-keyed id is
  present.

*Structure:*
- the B2 scrub, by both conditions (the class nodes and the KEGG polymers chitosan / PHB); the
  generalized-formula fill; no `mass` on a generic formula;
- `elements` ⊆ periodic table, including the `Rn` case.

*Roles:*
- no Metabolite node or metabolite edge for a role;
- a `ChemicalRole` node is created under its primary id, with its `is_a` ancestors up to a branch,
  and no `role` root;
- the "siderophore" subtree contains "phytosiderophore";
- a multi-branch role has several level-0 ancestors;
- a role edge appears on the family and its kept ancestors, with `substrate_depth` matching the
  metabolite rule on the same hierarchy;
- secondary-id roles are detected.

*Stereo:*
- the stereo family needs both conditions; the B2 scrub runs before families are built;
- an enantiomer pair with no unspecified node (L-/D-Val) gets exactly one `is_enantiomer_of` edge,
  lower → higher id.

*Paper names:*
- an ambiguous name with no override fails step 6; an **unmatched** name is unresolved (empty id,
  no edge, reported) and doesn't fail;
- a curated alias override beats the row's direct id;
- the unspecified node can win the dominant-form rule (the row stays on it);
- an anomer-only family keeps the target `as_reported`;
- "Lysine" with a single L-Lys match gets `dominant_form_assumed`;
- two names sharing a target are allowed only when both have alias entries;
- lock keys and targets are re-mapped through `id_aliases`;
- a failed run writes no file (lock and pin untouched);
- an old secondary role id (`chebi:10004`) resolves via `ChemicalRole.all_identifiers`;
- a paper-path ChEBI candidate (e.g. Tyrosine) is re-keyed and enriched at stage D;
- the alias yaml has priority;
- "Alanine" → L-Ala and "D-Alanine" stays D;
- a < 2× margin or a tie is refused;
- the first assignment is frozen, and drift is reported, not re-pointed;
- `stereo_assignment` is set for rule picks and override picks;
- a direct id gets canonicalised; a stereo-unspecified direct id goes through the rule; the anomer
  guard keeps GlcNAc / MurNAc / glucose on the anomer-unspecified form;
- step 6 fails on an unresolved ambiguity, a distinct-name collision, or a missing
  `stereo_assignment`;
- `name_key`: "Glutamic acid" and "glutamic acid" give one key; `Co` / `CO`, `NO` / `No` and `CO2`
  / `Co2` stay distinct (all allowed in one alias yaml); `NAD` / `nad` fold; the resolver skips the
  lowercased name step for `CO`;
- a same-`name_key`, different-id pair in one alias yaml is an error;
- `paper_overrides` beat the global lookup; two papers can map one name differently;
- an alias-override target that no reaction or transport reaches is still written and tagged
  while a row uses it; a stale alias entry is reported and its target is not kept on its account;
  an unresolvable target id fails step 6;
- `keep_separate:` blocks an MNX-basis merge; a stale identity-override entry is reported and
  skipped; contradictory `merge:` + `keep_separate:` fails step 6; no identity override adds a
  node;
- a frozen lock target is kept while its name still has rows; a stale lock entry is reported, not
  deleted;
- stereo-qualified detection: `L-Alanine`, `(S)-malate`, `beta-D-glucose` are qualified;
  `Alanine`, `Glutamate`, `beta-Citronellol` are not;
- an exact single match to an unspecified node (`Glutamate` → `C00302`) goes through the rule to
  `C00025`;
- `chebi_id` tie-break picks the charge-0 member;
- each dropped placeholder or measurement-only id appears in exactly one node's
  `recurated_from_ids` and in no `all_identifiers`;
- a paper candidate that isn't chosen and has no other reason to exist isn't written; the
  `metabolomics` tag sits only on final targets;
- `_resolved.csv` is written by step 6, and every `metabolite_id` in it is a written compound.

*Lumping pin:*
- the first run freezes today's set;
- a later threshold change on a frozen id is reported as drift, not applied;
- a newly kept id is auto-classified and appended (`lumping_set_added`);
- the adapter maps the pin to `substrate_breadth`.

*Classes:*
- the parser fails loudly on empty input;
- ancestors-only pruning; multi-class compounds; an alias membership lands on the canonical node.

*Pipeline:*
- byte-identical outputs on re-run.

**Adapter tests**:
- the new metabolism-adapter fields and stereo edges (with the guard);
- tcdb `ChemicalRole` nodes and role edges (with depth);
- assay `stereo_assignment` pass-through;
- the class adapter's dangling-proofing.

**KG validity** (`tests/kg_validity/`):
- no duplicate KEGG / primary-ChEBI / MNX / `all_identifiers` entry across `Metabolite`;
- 0 secondary ChEBI ids; 0 role-term metabolites; every `ChemicalRole.id` uses `chebi.role:`;
- no isolated `Metabolite`; every `metabolomics`-tagged node has an assay edge, and every node with
  an assay edge carries the tag;
- on every `chebi:`-keyed `Metabolite`, `chebi_id` equals the id's local part;
- the `ChemicalRole` subtree counts count **distinct** families over the `is_a` DAG (a family
  reachable through two parents counts once; same test shape as GO), and each subtree count is ≥
  its `direct_*` count;
- no `Metabolite_is_enantiomer_of` edge joins an anomer pair;
- no Gene / OrganismTaxon / Metabolite rollup reaches a `ChemicalRole`; no `ChemicalRole` carries
  `level_is_best_effort`;
- every `sources` value on the new edges joins a `DataSource` node (`kegg`, `chebi`);
- every re-keyed node carries the KEGG enrichment (structure where defined, pathway links where
  KEGG has them);
- no class node with a concrete `inchikey`; `elements` ⊆ periodic table; no `mass` on a generic
  formula; **every pre-change node with `elements` maps (through `metabolite_id_remap.tsv`) to a
  node with `elements`**, since node ids change and a plain superset can't hold;
- 0 same-slot duplicate assay edges; every assay edge has `stereo_assignment`;
- stereo-form and enantiomer edges carry `sources` / `evidence` and join distinct nodes; no
  enantiomer pair is stored in both directions;
- `MetaboliteClass` subtree count ≥ direct count;
- the full-text queries hit (bare id, "siderophore");
- the `lumping` `TcdbFamily` set equals `lumping_tcdb_ids ∩ kept_tcdb_ids`; every level-< 2
  family is `not_assessed`;
- `structure_kind` is dense; `identity_basis` appears only on merged or re-keyed nodes;
  `kegg_functional_categories` values ⊆ the br08001 level-0 names; `kegg_brite_classes` equals its derivation on every metabolite;
- no InChIKey on two `Metabolite` nodes (a KG assertion, beyond the step-6 check); every
  `ChemicalRole.all_identifiers` entry is unique across `Metabolite` ∪ `ChemicalRole`;
- 11 `DataSource` nodes (was 9). Any test or doc pinning 9 is updated;
- the `transport_substrate_resolution` distribution equals the pre-change baseline.

**Regression:**
- `import.report` clean;
- `/omics-edge-snapshot`: expression edges unchanged. The metabolism block **will show expected
  drops**, by design: `Metabolite`, `Reaction_has_metabolite`, `Organism_has_metabolite`,
  `Tcdb_family_transports_metabolite`. The script flags any negative delta as a regression, so
  record the expected deltas in the post-rebuild doc and compare against them instead of reading
  the flags as failures;
- assay edge counts change only as documented;
- the snapshot is regenerated;
- `capture_annotation_state.py --compare` shows 0 gene movement;
- the `Reaction_has_metabolite` and pathway-level `KeggTerm` counts move only as documented (§10).

---

## 12. Release

**Held.** This change ships inside a larger, bundled release. The sequence agreed with the
explorer:
1. The explorer E1 resolver ships first (a no-op on the current KG).
2. The KG release sets `mcp_min_version` and `breaking_changes`.
3. The explorer absorbs the new KG.

Two files, generated from the identity report, ship with the post-rebuild kg-changes doc:
- **`metabolite_id_remap.tsv`**: old id, new id, basis ∈ merge basis | `rekey_chebi_kegg_xref` |
  `removed_role` (empty new id; a `role_id` column gives the `ChemicalRole` it became, e.g.
  `chebi:10004` → `chebi.role:33229`, where TCDB's id was secondary) |
  `removed_placeholder` / `remeasured_elsewhere` (new id = the node the rows now point to).
- **`metabolite_assay_repoints.tsv`**: paper, table, paper name, old id, new id,
  `stereo_assignment`, reason. Closed set of reasons:
  - `alias_override`: a curated override (incl. GSH / GSSG, cystine);
  - `dominant_stereo_form`: the rule or the lock picked the form;
  - `case_fold_convergence`: case-insensitive keys joined the name to another paper's node;
  - `direct_id_canonicalised`: a paper-supplied id was merged away or re-keyed;
  - `placeholder_recurated`: was on an `mnx:` placeholder;
  - `remeasured_elsewhere`: was on a measurement-only node dropped at stage K;
  - `merged`: the old node merged into the new one.

  When several reasons apply, the first in this order wins.

The post-rebuild doc also lists:
- for each re-resolved name (the 36 B5 names **plus** the case-fold convergences): the old id,
  the new id, `stereo_assignment`, and `catalyst_gene_count` / `transporter_gene_count` on the new
  node. The summed-over-stereo-forms column was dropped (review doc §13a);
- the measured number of pathway-level `KeggTerm` nodes added with no gene attached (§10);
- the measured `structure_kind` counts (`defined` / `generic` / `polymer`), and how many
  `generic` / `polymer` nodes carry transport edges;
- the final `Metabolite_is_stereo_form_of` / `Metabolite_is_enantiomer_of` counts.

---

## 13. Out of scope

- TCDB plausibility, `substrate_depth`, per-organism transport reach, TCDB names, and the
  lumping-threshold re-tune (deferred spec).
- ChEBI-inferred classes, or a ChEBI class ontology in the graph.
- ChEBI roles **on `Metabolite`** (`Metabolite_has_chemical_role`): backlogged; needs a curated
  allow-list.
- A culture-level measurement node; the `spent_medium` compartment.
- Explorer / MCP changes (explorer repo).

## 14. Documentation

**`CLAUDE.md`:**
- the Metabolite bullet;
- a new MetaboliteClass bullet;
- a new ChemicalRole bullet (semantics: "transports something TCDB calls X", not a role of the
  family or its metabolites);
- the TcdbFamily bullet (`substrate_breadth` and the pinned lumping set, which
  replaces the "query the threshold" text);
- the label / relationship lists;
- the step-6 description (now also writes the paper `_resolved.csv` files) and the step-7 entry
  (retired, folded into step 6); `scripts/prepare_data.sh` (`--steps 7` alias);
- Data Locations: `cache/data/chebi/`, `brite_br08xxx.json`, the lock file, the identity report.

**Other docs:**
- `docs/kg-changes/metabolite-identity-and-classes.md`: finalised with measured numbers and the
  TSVs.
- `docs/kg-changes/tcdb-two-source-upgrade.md` §4 and `docs/kg-changes/vocabulary-contract.md`:
  the pinned lumping set and the R3 exception.
- `CHANGELOG.md`: `### Breaking` (id remaps, `chebi_id` values, removed role and placeholder
  nodes, paper re-resolutions), `### Added` (MetaboliteClass, identity properties, stereo
  edges, `stereo_assignment`, `ChemicalRole` + `Tcdb_family_transports_chemical_role`,
  `substrate_breadth` as the pinned lumping set,
  with no value change), `### Highlights` ("ask for every
  carbohydrate or peptide an organism makes, takes up or was measured with").
- `plans/backlog.md`: close both items; file the deferred TCDB spec.

**Project skills** (`.claude/skills/`):

*Use during implementation and validation:*
- `/omics-edge-snapshot`: capture a baseline **before** the rebuild and compare after (see the
  regression note in §11 for the expected metabolism deltas);
- `/cypher-queries`: post-rebuild spot checks (nitrate one node, L-Ala measured, siderophore
  route, stereo hop);
- `/release-kg`: later, for the bundled release. Its per-paper expression-regression check is
  unaffected; the release cut adds `breaking_changes` and the `mcp_min_version` bump.

*Update in the same change:*
- **`paperconfig/SKILL.md`**:
  - the "ID resolution priority" section, with the new order (direct id → paper overrides → lock
    → match → dominant-form rule) and no ambiguous pick;
  - the extended `metabolite_aliases.yaml` form `{id, stereo_assignment, reason}`;
  - the stereo-qualified name rule and case-insensitive keys;
  - "step 7 just looks up …" → step 6 writes `_resolved.csv`;
  - the perf loop `--steps 6 7 --force` → `--steps 6 --force`.
- **`add-a-strain/SKILL.md`**:
  - the step-output table: step 6 now also writes `_resolved.csv`, and step 7 is retired (the
    `--steps … 7` lines stay valid through the alias, but the text should say so);
  - after onboarding, read the identity report: `lumping_set_added`, `stale_*` entries, and new
    ChEBI-rekeyed compounds.
- **`cypher-queries/SKILL.md`**:
  - metabolite templates: lookup via `all_identifiers` (CURIE), the `ChemicalRole` siderophore
    route, `MetaboliteClass` membership, the stereo-form / enantiomer hop;
  - fix the stale `Tcdb_family_transports_metabolite` hint (use `substrate_depth`, not
    `level_kind = 'tc_specificity'`).
- **`omics-edge-snapshot`** (script + SKILL.md): add `MetaboliteClass` and `ChemicalRole` node
  counts and the new edge types (`Metabolite_is_stereo_form_of`, `Metabolite_is_enantiomer_of`,
  `Metabolite_in_metabolite_class`, `Tcdb_family_transports_chemical_role`) to the metabolism
  block.
- No change: `refresh-mnx` (MNX is still read the same way; ChEBI is fetched by step 6 itself),
  `integrate-a-tool`, the per-strain tool skills.

**Step-7 retirement, other consumers** (update in the same change):
- `tests/test_resolve_paper_metabolites.py`: port its cases to the step-6 writer tests, then
  delete;
- `multiomics_kg/adapters/metabolite_assay_adapter.py`: the docstring says `_resolved.csv` is
  "written by step 7";
- `docs/kg-changes/metabolomics-extension.md`: the step-7 description;
- `scripts/prepare_data.sh`: the header comments (`--steps 6 7`), the default `STEPS` and
  `--rebuild` lists (drop 7, keep the no-op alias);
- `CLAUDE.md`: also correct the `metaboliteFullText` field list (live: `name` only), and the
  DataSource count (9 → 11).

**Other consumers of changed formats** (update in the same change):
- `scripts/validate_paperconfig.py`: validate the extended `metabolite_aliases.yaml` form
  (`{id, stereo_assignment, reason}`), the same-`name_key` collision rule, and the
  `stereo_assignment` requirement for stereo-variant targets.
- `.claude/skills/paperconfig/SKILL.md`: document the alias form and the stereo rules for new
  metabolomics papers.
- `.gitignore`: `cache/data/chebi/raw/`.
