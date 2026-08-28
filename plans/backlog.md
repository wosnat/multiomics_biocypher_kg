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

## GEO processed-supplements drop

Deferrals created by the 2026-08-19 execution pass of
`plans/geo_paperconfig_updates.md` (the "easy parts" shipped; these did not).
Section references below are into that plan file, which holds the full designs.

- [ ] **Johnson 2026b — diel periodicity metrics (§2.1B).** `media-3`
      `all_cosinor_fit` (1,872 rows) → `derived_metrics_table` with 2–3 numeric
      metrics (`diel_acrophase_rad` new, `diel_amplitude`, `peak_time_h`).
      Design fully written; deferred as a scope choice (only §2.1A DE shipped).
      → `plans/geo_paperconfig_updates.md` §2.1B

- [ ] **Johnson 2026b — MED4 iModulons as `gene_clusters` (§2.1C).** 32 modules,
      1,011 membership rows already long-form in `media-8`; per-module
      descriptions transcribed from the published `iModulons` sheet (no LLM).
      The one genuinely new mechanic of the GEO drop (multi-membership +
      signed `gene_weight` as `score_col`); check `Gene_in_gene_cluster` 1:1
      assumptions in MCP/explorer before wiring. PCC7942's 78-module model and
      the `media-5` module↔module homology ride behind it.
      → `plans/geo_paperconfig_updates.md` §2.1C–D

- [ ] **iModulon activity layer (D5).** The A matrix (32 modules × 248 samples)
      and DIMA tables have no schema slot — a measurement whose subject is a
      gene set. Needs its own spec; preferred shape is a `MetaboliteAssay`-style
      assay node. Settle together with the TF-binding-site layer, where a
      `RegulatoryModule` promotion becomes worth revisiting.
      → `plans/geo_paperconfig_updates.md` D5

- [ ] **Hackl 2023 — antisense DE arm (D8, open).** Three antisense DESeq2
      tables report DE of the antisense transcript at a sense gene's locus;
      recommended encoding is separate Experiments whose `name` /
      `table_scope_detail` say "antisense-strand transcription at this locus".
      Sense arm shipped without it.
      → `plans/geo_paperconfig_updates.md` §2.5 + D8

- [ ] **Hackl 2023 — MIT1327 island remap.** Hackl used a different MIT1327
      assembly (2.58 Mb / 23 contigs vs our `GCF_001632125.1` / 31); its 11
      islands need a contig+coordinate remap before transfer. The other 11
      strains' islands shipped.
      → `plans/geo_paperconfig_updates.md` §2.4 + Tier 3

- [ ] **munoz 2022 (GSE154594).** Environmental station-ALOHA samples on a
      custom Agilent array (`GPL28884`); needs probe→gene mapping AND a
      reference-proteome-match-style organism (`Marinobacter (MarRef v6)`
      precedent), plus a 129 MB RAW.tar download. A different kind of
      integration — needs its own plan, not a strain add.
      → `plans/geo_paperconfig_updates.md` Tier 3;
      `docs/kg-changes/reference-proteome-match-organisms.md`

- [ ] **TSS / operon / UTR / ncRNA entity layer.** No node types today. Would
      serve Voigt 2014 (TSS/operons/UTRs/ncRNAs — shipped only as reduced
      per-gene scalars, §2.6), Steglich 2010 (operon half-lives, MOESM6/7) and
      Doron 2016 (45 asRNAs, moesm38) together — three papers in hand is what
      makes it a spec rather than a one-off.
      → `plans/geo_paperconfig_updates.md` Tier 3 + §2.6

- [ ] **TF binding sites / regulator→target layer.** Johnson `media-4` (75 RpaB
      sites with coordinates, motif, q-value) + `TF_regs` on every iModulon.
      The KG has no regulator→target edge. Related to D5 — an iModulon's
      regulator and a TFBS are two views of the same missing edge; settle the
      two specs together.
      → `plans/geo_paperconfig_updates.md` Tier 3 + D5

- [ ] **B1 regression check.** Guard in `build_gene_id_mapping` (or its report
      tests) that no gene accumulates more than ~3× the median tier-1 id count —
      proposed with the B1 fix (`a84db12b`) but not in the commit.
      → `plans/geo_paperconfig_updates.md` Blocker B1

- [ ] **Runaway-mapping residual: KT2440 `M8001_03425` (40 tier-1 ids).**
      Re-check done 2026-08-28: all three `gene_id_mapping.json` files were in
      fact rebuilt by the B1 fix commit (`a84db12b`), not before it. W3-18-1
      (median 6 / max 7) and PCC7002 (6 / 8) are clean; KT2440 (median 9) has
      one gene at 40 tier-1 ids — >4× the median, the shape the B1 guard above
      would catch. Inspect which paper's `id_columns` feeds it before the next
      KT2440 paper touch.
      → `plans/geo_paperconfig_updates.md` Blocker B1 consequences

### Post-review additions (2026-08-19 subagent review of all wired paperconfigs)

- [ ] **Voigt 2014 — Table S7 conserved-TSS ortholog comparison.** Per
      MED4↔MIT9313 ortholog pair: TSS distance in each strain + delta.
      Reducible to a promoter-conservation DerivedMetric, but a two-strain
      fact sits awkwardly in DerivedMetric's single-organism shape — needs a
      small design decision before wiring. Table S8 (prochlorosin/procA TSS
      scalars) also skipped: mostly new-annotation genes absent from the gene
      layer, near-zero resolution expected.
- [ ] **Hackl 2023 — tycheposon element catalog (`pro-623-elements.tsv`).**
      937 elements with coordinates; contigs of 7 deployed strains appear
      (AS9601, MIT0604, MIT1314, MIT9303, MIT9312, MIT9313, RSP50 + excluded
      MIT1327). Could ride the same coordinate-containment machinery as the
      islands (cluster_type e.g. `tycheposon_element`) — arguably the paper's
      core contribution. mmc3 (HMM profiles) / mmc4 (element-type summary) /
      mmc5 (RT-qPCR primers) reviewed and correctly out of scope.
- [ ] **Johnson 2026b — `media-7` per-iModulon enrichment tables** (GO / COG /
      KEGG-Module / KEGG-Pathway per module). Not independently wired-worthy —
      should feed the §2.1C cluster `functional_description`s when 2.1C lands.
- [ ] **he 2022 — recover ~6 MED4 genes via a strip-`gene-` heuristic.** GEO
      GeneID values are `gene-<tag>`; 6 protein-coding rows
      (gene-PMM0220/0236/0950/1858/2002/2065) fail only because no resolution
      pass strips the `gene-` prefix before the multi-singleton lookup
      (`_heuristic_candidates` in `gene_id_utils.py`). ~0.3% gain; rerun
      step 4 for the paper after.
- [ ] **MED4 mapping quirk — `RNA_41`/`gene-RNA_41` absorbed into PMM0521
      (frr).** Likely a frr/rrf 5S-rRNA name collision during the GCA/GCF
      harvest; produces one spurious duplicate resolution in he 2022 (benign:
      not_significant edge). Fix in the mapping builder, not per-paper.

## Explorer / MCP coordination

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
