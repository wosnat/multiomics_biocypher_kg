# Backlog

Deferred and cross-cutting work that has no home in a single plan file.

**Convention.** One bullet per item: *what*, *why it was deferred*, and a pointer
to the spec section, plan file, or code comment that created the deferral. Add
the bullet in the **same commit** that creates the deferral — that is the whole
discipline, and the reason this file lives beside the code rather than in an
issue tracker (the GitHub issues opened in Jan–Feb 2026 were all completed and
none were ever closed).

Remove a bullet when the work lands; the CHANGELOG records that it did.
Items that are one focused change get their own `plans/<topic>.md` when picked
up — this file is the index, not the plan.

---

## KG semantics

- [ ] **MEROPS GO bridge — rejected on measurement (2026-08-18), revisit only
      with a concrete use case.** All-kingdom member rollup yields median 19 /
      max 389 GO terms per family incl. eukaryote-only terms; completeness win
      is 338 genes (~8%) vs the 1,311 that justified TCDB's GO bridge. Needs a
      filtering design (e.g. >= N supporting identifiers) before it can land.
      → `docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`

- [ ] **Recalibrate `is_multi_substrate` / `is_multi_gene` thresholds.** The TCDB
      threshold (`level >= 2 AND metabolite_count >= 50`) was calibrated against
      the pre-pruning node set of ~12.9K nodes; the graph now has 1,515. Held
      back from the vocabulary-contract change so the rename's validate diff
      stays readable.
      → same spec §9.3, and the `is_promiscuous` note in `CLAUDE.md`

- [ ] **Convert the grandfathered `"true"` / `"false"` properties to meaningful
      pairs.** R5 of the vocabulary contract forbids native `bool` and deprecates
      stringified booleans, but seven released properties predate it:
      `Experiment.is_time_course`, `Experiment.reports_fold_change`,
      `DerivedMetric.rankable`, `MetaboliteAssay.rankable`,
      `DerivedMetric.has_p_value`, the DM edge `significant`, and the DM flag
      edge `value`. All are MCP-read, so this is **breaking** and needs its own
      `post-import-validate` baseline. Follow the
      `has_cross_genus_members: cross_genus | single_genus` precedent.
      → `docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md` §3 R5, §10.5

- [ ] **Orphan proteins.** ~46% of UniProt proteins were reported with neither
      `Protein_belongs_to_organism` nor `Gene_encodes_protein`. The two
      kg-validity tests documented as failing on it have now passed clean on
      **two consecutive 2026-08-17 rebuilds** (InterPro multi-ontology + MEROPS),
      so what remains is verification, not a live failure: confirm the orphan
      fraction itself is gone (vs. the tests under-asserting), then close this
      and the CLAUDE.md Known Issues entry together.
      → `plans/orphan_proteins.md`

- [ ] **InterPro MetaCyc pathway xrefs.** Populated in
      `cache/data/interpro/interpro_reference.json` (`pathways`, 5,091 entries)
      but not in the graph. InterPro ships no KEGG xrefs, so this would be a new
      pathway vocabulary rather than an extension of the KO-derived layer.
      Reactome is excluded by default (species-expanded, noisy for marine
      bacteria).
      → `docs/kg-changes/interpro-multi-ontology.md` (supersedes
      `interproscan-extension.md`)

> InterPro multi-ontology redesign follow-ups (NCBIfam MCP registration,
> naming-recovery extensions, etc.) have their own plan file:
> `plans/interpro_redesign_backlog.md`.

## Explorer / MCP coordination

- [ ] **Relationship-property index on `evidence`.** Explicitly not requested for
      this release — current edge-property filters touch
      `Tcdb_family_transports_metabolite` (11,263) and `Gene_has_tcdb_family`
      (53,763), where it does not matter. It starts to matter if the deferred W2
      workstream lands `source_filter` / `evidence_filter` over
      `Gene_involved_in_biological_process` (539,873) and `Gene_has_pfam`
      (177,453). The graph currently has 86 indexes and zero relationship-property
      indexes.
      → explorer `docs/kg-specs/2026-08-16-interpro-tcdb-asks.md` KG-IPT-008

- [ ] **MCP surfacing of InterPro two-layer provenance.** Source / evidence
      filters on the gene→ontology tools, and a 2-hop router mode over the
      Layer-A `Interpro_entry_related_to_*` edges. Explorer-side work; unblocked
      now that §7.2 of the vocabulary spec establishes the GO provenance shape is
      final this release rather than pending.
      → `docs/kg-changes/interpro-two-layer.md` §7

- [ ] **Register the MEROPS ontology in the explorer.** `MeropsFamily` /
      `Gene_has_merops_family` are live in the KG but invisible to
      `ontology_landscape` / `search_ontology` / `genes_by_ontology` until the
      explorer's `ONTOLOGY_CONFIG` + ontology enum gain a `merops` entry
      (currently only `run_cypher` reaches them). Registration should surface
      `call_class` — the guard that keeps dead homologs and inhibitors out of
      protease counts — and `peptidase_gene_count` as the default count.
      → `docs/kg-changes/merops-extension.md` "What does NOT change"
