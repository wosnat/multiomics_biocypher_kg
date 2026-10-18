# MNX merges across stereo / anomer / class boundaries (explorer ask 2026-10-15)

Source: `multiomics_explorer/docs/kg-specs/2026-10-15-metabolite-stereo-merge-asks.md`.

## Root cause

`canonical.canonicalise` lets `mnx_same_entity` (and `manual`) links bypass the KEGG↔KEGG
guard. That is by design ("two KEGG ids merge only on MNX evidence or a curated merge"), but MNX
sometimes treats a stereo-unspecified compound, or a generic sugar and its β-anomer, as one
entity. All 21 multi-KEGG nodes come from `mnx_same_entity`. 11 of them break §1.

The winner rule is **not** the cause. `_choose` already picks the KEGG id with more references
from gene-reachable reactions (ask 3a). C05345 β-F6P really has 10 kept reactions against 9 for
C00085, and C05378 has 5 against 4. Once the stereo merges stop, that rule no longer has to pick
between these ids.

## Fix

1. **Stereo guard on MNX links** (`canonical.py`). Refuse an `mnx_same_entity` link between two
   groups with different KEGG ids when some KEGG pair across them is *stereo-distinct*. Two KEGG
   ids are stereo-distinct when:
   - one has a non-class ChEBI term whose InChIKey shares the other's connectivity block, and
   - that InChIKey's stereo layer appears on none of the other's terms in that block, and
   - the two term names are not a geometric (cis/trans, E/Z) pair. §3.6 handles those
     separately; this exemption keeps caffeate merged.

   Refused links go to `rejected_links` with basis `mnx_same_entity`.

   Prototype result: it rejects exactly the 9 stereo/anomer pairs and passes all 10 pairs the
   explorer calls fine. `manual` links stay exempt.

   Once split, stage I's existing ChEBI `is_a` walk emits `Metabolite_is_stereo_form_of` for
   C00197→C00597, C01102→C05702, C03227→C02794, C04272→C04039, C06010→C00900, C16241→C00725 and
   C05345→C00085 / C05378→C00354 (anomer edges are **still emitted**: §21 says so, and the
   explorer ask misreads it). C17345/C04644 get no edge, because ChEBI has no `is_a` between them.

2. **`keep_separate` overrides** for the two pairs a structural test can't see:
   - C00193 Aromatic aldehyde / C00261 Benzaldehyde (class → member);
   - C03688 Apo-[ACP] / C16240 [Lipoyl-carrier protein]-L-lysine (MNX conflation of
     protein-bound species).

3. **MNX-derived fields on a refused pair.** Both nodes carry the same `mnxm_id`, and so the same
   `inchikey` / `smiles` / `formula` / `mass` / `hmdb_id` / `mnx_chebi_aliases`. Stage H would raise
   `IdentityCollisionError` on the shared `mnx:` / `inchikey:` / `hmdb:` identifiers. Policy
   (**needs sign-off**):
   - the MNX fields stay on the node whose ChEBI InChIKey equals MNX's InChIKey. Protonation is
     ignored, so β-F6P keeps MNXM1107898 and its HMDB id;
   - if no node matches (C03688 / C16240), they stay on the node `_choose` would have made
     canonical;
   - the other node drops `mnxm_id`, `smiles`, `hmdb_id` and `mnx_chebi_aliases`. It takes
     `inchikey` from its primary ChEBI term when that term has one, and keeps `formula` / `mass`
     (identical across stereoisomers).

4. **Tests.**
   - Unit tests for the guard (stereo pair refused, anomer pair refused, geometric pair and
     conjugate pair allowed, `manual` exempt) and for the field split.
   - Stage-P `validate_build` check: no written node's `kegg_compound_ids` holds a stereo-distinct
     pair, so the build fails before writing.
   - `tests/kg_validity/test_metabolite_identity.py`: the 11 explorer pairs live on different nodes.

5. **Rebuild + docs.**
   - `prepare_data.sh --steps 6 --force`, regenerate the release TSVs
     (`scripts/metabolite_identity_tsvs.py`), Docker rebuild, `/omics-edge-snapshot` compare,
     `pytest -m kg`.
   - Reply in `docs/kg-changes/metabolite-identity-and-classes.md` §37 and in the explorer ask doc,
     plus a CHANGELOG `### Fixed` entry.

## Not in this change

- **Ask 3b** (prefer chloride over HCl, and the same for HBr / HF). A "prefer the charged ChEBI
  term" rule would also move `Sulfite` C00094 onto hydrogensulfite C11481 (ChEBI xrefs C00094 to
  sulfurous acid), so it needs its own rule. Proposed: a small curated `canonical:` pin list in
  `metabolite_identity_overrides.yaml`. Deferred unless wanted now.
- **Ask 5** (HMDB0000122 on β-D-glucose). MNX itself puts HMDB0000122 on MNXM1364060
  (β-D-glucose); HMDB's "D-Glucose" record carries the β-pyranose structure. CHEBI:17634
  D-glucose has no KEGG xref and no InChIKey, so no node gets it. The fix needs a curated xref
  move (an `xref_move:` override) rather than an algorithm. The anomer-xref count the explorer
  asked for is reported with this change; the move itself is a follow-up.

## Decisions (explorer reply 2026-10-15)

- MNX fields go to the node whose ChEBI InChIKey matches MNX's. The generic node keeps its own
  ChEBI ids, and a ChEBI term that is a written KEGG node's own xref is never published on a
  sibling through MNX aliases.
- Anomer stereo edges are re-created (§21).
- **3b is in this rebuild:** `canonical:` pins for `C00698` chloride, `C01324` bromide and
  `C00742` fluoride. The old ids stay in `all_identifiers` and in `metabolite_id_remap.tsv`.
- **Ask 5:** report the anomer-xref count in this rebuild. The `HMDB0000122` / `CHEBI:17634` →
  `C00031` move also goes in if the override stays small.
- **Ask 6 is in this rebuild:** `kegg_brite_classes` + its ControlledVocabulary +
  `kegg_compound_categories` → `kegg_functional_categories` (### Breaking).
- **Proposed (awaiting explorer):** the canonical choice needs a ≥ 2× reaction margin. Otherwise
  a `canonical:` pin is required, else stage P fails with `canonical_margin_unresolved`.
  Measured: 0 groups affected after the guard.
- **Agreed (explorer reply 2):**
  - the margin rule fails the build;
  - the message lists members, names, both counts and a ready-to-paste `canonical:` block, also
    written to the report;
  - a fallback measure (all KEGG reaction refs, same ≥ 2×) applies only when the gene-reachable
    counts are 0 vs 0. There are 0 such groups today.

  Design closed. Next: implement.

## Status (2026-10-15): DONE

Implemented, step 6 rerun against the release baseline (`5b13d00e^`) for the TSVs, Docker rebuilt,
`pytest -m kg` run. Write-up: `docs/kg-changes/metabolite-identity-and-classes.md` §38.

Deviations:
- the `kegg_brite_classes` vocabulary carries no `value_descriptions` (its values contain the `": "`
  separator; the loader now rejects that);
- a failed margin check reports only in the error and log, not in the report.

Addition: the shared-ChEBI-term owner rule.
Ask 5 follow-up: DONE 2026-10-15 — `xref:` override kind + 6 reviewed entries (change doc §39).

## Ask 8 (2026-10-18): DONE

The stereo-form rule is rebuilt on ChEBI chemistry: an `is_a` path of ≤ 3 steps with unchanged
formula and charge, read from each node's own identity, plus `stereo_link:` / `no_stereo_link:`
overrides. Narrowed D (a `chebi:` node's structure is its own ChEBI term's) adds
`removed_stereo_xrefs` and the remap basis `stereo_form_removed`. Stereo edges 76 → 99.
Write-up: change doc §40.
