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

- [ ] **Orphan proteins — verified still real (2026-08-18); NOT closeable.**
      Direct query on the freshly rebuilt graph: **25,441 / 67,024 proteins
      (38%)** have neither `Protein_belongs_to_organism` nor
      `Gene_encodes_protein` (down from the ~46% originally reported). The two
      kg-validity tests pass only because their thresholds were loosened to
      `< 50%` (docstrings still say "currently failing") — the "passed clean on
      two consecutive rebuilds" observation was the loosened assertion, not a
      fixed data gap. The original investigation stands: determine whether the
      WP_-join gap is pre-existing or a regression from the Feb 2026 refactor
      (`fe5c2bb`), fix or accept, then re-tighten the test thresholds.
      → `plans/orphan_proteins.md`

- [ ] **TIGRFAM→TigrRole bridge — rejected on measurement (2026-08-18), revisit
      only with a concrete cross-genus use case AND a coverage-correction
      design.** The archived JCVI role links (NCBI FTP
      `hmm/TIGRFAMs/release_15.0/{TIGRFAMS_ROLE_LINK,TIGR_ROLE_NAMES}`, frozen
      2018) would let `NcbifamFamily` TIGR* nodes bridge to the existing
      `TigrRole` nodes — role_id spaces are identical, 1,579 of our 2,204 TIGR
      families carry an informative role, 91% concordant with Cyanorak's curated
      assignments where both exist. Measured gain: 19,596 genes across the 21
      role-less organisms would get a role edge, **but 19,583 of them are
      already `informative_multi`** (zero dark-gene rescue), only 1,141
      `gene_category='Unknown'` genes recategorize (~3%), and coverage is
      hard-capped at ~20–26% per heterotroph genome (TIGRFAM-hit-bounded; the
      2,753 NF* families post-date the frozen role system and can never join)
      vs 89–94% curated coverage on Pro/Syn — so the one distinctive win,
      a shared cross-genus role axis, would be ~4× coverage-biased against the
      heterotroph side and misleading without corrected backgrounds. COG /
      KEGG / BRITE / GO already provide uniform cross-genus category layers.

- [ ] **NCBIfam→GO bridge — rejected on measurement (2026-08-18), revisit only
      with a GO-corroboration design that discounts same-scan sources.**
      `hmm_PGAP.tsv`'s curated per-family `go_terms` are already parsed into
      `cache/data/ncbifam/ncbifam_reference.json` (11,480 families, 8,279 NF* +
      3,201 TIGR*) but unused; 3,387 of our 4,957 observed families carry GO
      (~8.3K would-be bridge edges). Measured over all 47,324 genes with an
      NCBIfam edge: 120,222 candidate (gene, GO) pairs, **93% already on the
      gene** — and 77% of that corroboration is circular (the existing edge
      already lists `interproscan`, the same scan lineage that would deliver
      the bridge; eggNOG-Pfam/InterPro-Pfam one-source precedent). Informative
      additive content: 7,118 pairs (2,617 refinements + 4,501 novel; 1,429
      coarser no-ops) across 5,618 genes; dark-gene rescue (no GO at all →
      gains a specific level ≥ 4 term) is **344** — at the MEROPS-GO rejection
      line (338) and far under TCDB's justifying 1,311. Independent (non-scan)
      corroboration exists (25,978 pairs, ~24.4K single-source terms could gain
      a second voice) but has no consumer until a GO evidence-score design
      exists. The `ncbi` GFF GO source is the same curation lineage frozen at
      assembly-annotation date, which is why the marginal gain concentrates on
      old assemblies.

- [ ] **TigrRole hierarchy normalization.** The 114 `TigrRole` nodes are flat
      (`level = 0` everywhere) with the JCVI mainrole/subrole two-level scheme
      embedded in compound names ("Energy metabolism / Electron transport") —
      the only hierarchical ontology in the KG with no `is_a` edges, contra the
      unified-level convention. Split into mainrole (level 0) / subrole
      (level 1) nodes + `Tigr_role_is_a_tigr_role` edges, and flag the junk
      role nodes ("Not Found", "Unclassified", "Hypothetical proteins",
      "Unknown function", "Disrupted reading frame /") as uninformative.
      Hygiene-only — fold into the next `functional_annotation_adapter` touch
      that already forces a rebuild; not worth one on its own.
      → `multiomics_kg/adapters/functional_annotation_adapter.py`
      (`_tigr_role_node_id`, `MultiCogRoleAnnotationAdapter`)

- [ ] **InterPro MetaCyc pathway xrefs — measured 2026-08-18, benefit is thin;
      revisit only with a concrete use case.** Populated in
      `cache/data/interpro/interpro_reference.json` (`pathways`, 5,091 entries)
      but not in the graph. InterPro ships no KEGG xrefs, so this would be a new
      pathway vocabulary rather than an extension of the KO-derived layer.
      Measured over all 42 strains (124,751 genes): 23,744 genes would gain a
      MetaCyc pathway, but 76% already have a KEGG pathway (and the delivery is
      entry-level family inference, weaker than the per-gene KO layer). The
      no-KEGG win (5,836 ungated / 3,792 FAMILY-gated) is driven by fold-level
      superfamilies with heavy fan-out (p90 = 55 pathways/gene among gainers);
      the truly dark-gene rescue (no KEGG *and* no GO) is **554 ungated / 271
      FAMILY-gated** — MEROPS-GO-rejection territory (338) vs the 1,311 that
      justified TCDB's GO bridge. Practical blocker: `interpro.xml` xrefs are
      bare dbkeys (no name attribute), and MetaCyc pathway names are
      registration/license-gated, so nodes would ship as unreadable `PWY-XXXX`
      ids. Reactome remains excluded by default (species-expanded, noisy for
      marine bacteria).
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

- [ ] **File an upstream Bioregistry new-prefix request for `ncbifam`.**
      KG-SYNC-002 (2026-08-19) minted `ncbifam:` as a house colon-CURIE prefix
      for `NcbifamFamily` node ids — `ncbifam` is registered nowhere today
      (verified live against bioregistry, identifiers.org, and the Biolink
      prefix map; only `tigrfam` exists and its `^TIGR\d+$` pattern cannot hold
      NF accessions). NCBIfam is a real, active NCBI resource (InterProScan's
      member-DB name; successor of TIGRFAMs), so a registration request at
      https://github.com/biopragmatics/bioregistry (new-prefix issue template)
      would make the graph retroactively registry-correct. Suggested entry:
      pattern `^(TIGR|NF)\d+$`, homepage
      https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/.
      → explorer `docs/kg-specs/2026-08-19-presync-kg-asks.md` KG-SYNC-002 / §6
