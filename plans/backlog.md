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

- [x] **TIGRFAM→TigrRole — SHIPPED 2026-08-29 as (a) an ontology bridge +
      (b) equivalog-gated gene edges.** The 2026-08-18 measurement below stands
      and is kept for the record; the gene-level rejection was **superseded**
      once the transfer was gated on `family_type = 'equivalog'` (JCVI's "same
      function in every member" type, the only one meant for unconditional role
      transfer): multi-role genes 375 → 15, conflicting mainroles 308 → 14,
      Cyanorak contradictions 7% → 5.4% (sampled — facet choices on multi-role
      proteins, not errors, so both views ship as two edges), and
      `annotation_quality` movement **0**. The coverage-correction requirement
      is met by documentation rather than by a correction factor: ORA over
      `TigrRole` must use a per-organism background of genes with ≥ 1
      `Gene_has_tigr_role` edge (~700–860 genes/heterotroph genome vs
      ~1,800–1,976 on MED4). Shipped: two-level `TigrRole` +
      `Tigr_role_is_a_tigr_role`, `Ncbifam_family_has_tigr_role` (~1.9K, router
      only for the non-equivalog families), inferred `Gene_has_tigr_role` on all
      43 organisms, `gene_category` fill-only priority 4, `[tigr_role_inferred]`
      description lines.
      → `docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md`,
      `docs/kg-changes/tigr-role-bridge.md`

      *Original 2026-08-18 rejection, for the record.* The archived JCVI role links (NCBI FTP
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

- [ ] **Multi-source `Gene_has_ncbifam_family` — measured 2026-08-28, deferred.**
      The edge is `interproscan`-only; two other sources name NCBIfam/TIGRFAM
      families and were checked because the TIGR-role merge discards them.
      Cyanorak `protein_domains` `TIGR*` tokens are **93% already in** the
      gene's InterProScan hits (1,244 extra ids over 903 genes, Pro/Syn only).
      PGAP `inference` HMM ids are **84% `NF*` `domain` / `PfamEq` models that
      InterPro's NCBIfam member DB excludes by design** — 0 of 1,860 `PfamEq`
      and 108 of 14,234 `domain` families are ever observed across 43 genomes,
      and the ids interleave numerically with observed ones, so this is a
      deliberate member-DB scope difference, **not a version gap**. Adding both
      would buy corroboration only (~480 TIGR + ~205 `NF*` PGAP confirmations)
      and no new families of consequence. If ever done, gate on `for_naming` /
      family-grade `family_type`s rather than ingesting every token.
      → `docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md` §2
      ("Dropped sources"), `docs/kg-changes/tigr-role-bridge.md`

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

- [ ] **MIT9313 — `MIT9313_genbank.tsv` `PMTid` column is stale for ~13 genes.**
      The resource table (an older RefSeq build of the same assembly) gives
      e.g. `PMT_0166 → AKG35_RS12635` and `PMT_0956/0959/0962 → AKG35_RS05225`,
      where the current GFF has `PMT_0166` on `AKG35_RS13725` and RS05225 as
      `PMT0983`; typed `old_locus_tag` (Tier 1) these surface as the Tier-1
      conflicts left on MIT9313 after 2026-08-29. The other 49 conflicts (the
      `PMT_nnnn`/`PMTnnnn` "two locus-tag families") were one gene each with a
      re-called start codon and are gone — fixed in `_position_fallback_merge`
      (stop-codon/reading-frame criterion + per-contig offsets, see
      `docs/methods_position_fallback_merge.md`). Options: drop `PMTid` from the
      TSV's `id_columns` (the GFF already carries every `PMT_` tag) or demote it
      to a lower tier. Low value; the conflict mechanism already refuses them.
      → `data/Prochlorococcus/papers_and_supp/MIT9313_resources/paperconfig.yaml`,
      `cache/data/Prochlorococcus/genomes/MIT9313/gene_id_mapping.json` `conflicts`

- [x] **TigrRole hierarchy normalization — DONE 2026-08-29**, folded into the
      TIGRFAM-role bridge above exactly as this bullet anticipated (one
      `functional_annotation_adapter` touch that already forced a rebuild).
      Subroles are `level 1` / `level_kind = 'tigr_subrole'`, new mainrole nodes
      are `level 0` / `tigr_mainrole`, joined by `Tigr_role_is_a_tigr_role`; the
      5 junk subroles keep their flags and the 3 junk mainroles gained them.
      → `docs/kg-changes/tigr-role-bridge.md`

      *Original bullet.* The 114 `TigrRole` nodes are flat
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
## Analysis-driven asks — pro_n_import_coculture (2026-10-13)

Filed from `multiomics_analysis/analyses/2026-10-06-pro_n_import_coculture/gaps_and_friction.md`
and `docs/methodology-review-2026-10.md` (round-2 review). Only the KG-side items;
the explorer tool bugs (non-unique paging sort, `gene_neighbors` window OOM, silent
`min_gene_set_size` defaults, `metabolite_elements` drop count, flat `list_metabolites`
reach) are filed in the explorer repo. Sizes below were measured on 2026-10-13 against
the committed caches, not taken from the log. Ordered by proposed priority.

**P1 — small, mechanical, unblocks the analysis's own workarounds**

- [x] **Fall back to Cyanorak coordinates for genes NCBI does not call — SHIPPED
      2026-10-13 (+ sequence from the Cyanorak translation), see
      `docs/kg-changes/cyanorak-coordinate-fallback.md`; awaiting rebuild.**
      3,837 genes across the 21 Cyanorak strains (MED4 96 / 1,973; MIT9313 517 /
      2,906; MIT9303 504; WH8109 362; CC9311 339) have no `contig`/`start`/`end`
      and no `protein_id` — every one is Cyanorak-only (`locus_tag_ncbi` null).
      `gene_mapping.csv` already carries `start_cyanorak`/`end_cyanorak`/
      `strand_cyanorak` for them, so this is a merge-time fallback in step 2 (or
      the adapter), not new data. Caveat: Cyanorak stores drafts as one
      concatenated record, so PAC1/SB-style assemblies need the per-contig offset
      remap `_contig_offsets` already uses; closed genomes are a straight copy.
      Add a `coordinate_source ∈ {ncbi, cyanorak}` property (R5) so consumers can
      tell. Unblocks `gene_neighbors` and the `gene_org_contig_start_idx` seek for
      these genes. → friction log 2026-10-07 "91 MED4 genes have no coordinates";
      `multiomics_kg/download/build_gene_mapping.py`, `cyanorak_ncbi_adapter.py`

- [x] **Populate `coculture_partner` when coculture is the background, not the
      treatment — SHIPPED 2026-10-13** as a dense list + `partner_role` on the edge
      (`docs/kg-changes/coculture-partner-list.md`; explorer dialog
      `multiomics_explorer/docs/kg-specs/2026-10-13-coculture-partner-list-explorer-asks.md`).** The adapter sets it only from `treatment_organism`
      (`omics_adapter.py` ~l.390), so the 8 Weissberg 2025 N-starvation
      experiments run *in* coculture (`background_factors: [coculture, …]`,
      MED4 + HOT1A3, RNA-seq + proteomics) and García-Fernández 2022's glucose
      proteomics have `coculture_partner` null and are invisible to a partner
      filter — the analysis first concluded "no Prochlorococcus proteomics in
      coculture", wrongly. Proposal: optional paperconfig `background_organism`
      (+ `background_taxid`) on an experiment; adapter populates
      `coculture_partner` from either field and emits `Tests_coculture_with`
      with `partner_role ∈ {treatment, background}` (new R5 vocab entry). The
      Publication `organisms` rollup then covers both. → friction log 2026-10-06
      "`coculture_partner` filter misses…"; `docs/kg-changes/experiment-stats.md`

- [x] **Metabolite equivalence across namespaces and stereo/charge variants.** CLOSED 2026-10-14 — canonical identity + stereo/enantiomer edges; see `docs/superpowers/plans/2026-10-14-metabolite-identity-and-classes.md` and `docs/kg-changes/metabolite-identity-and-classes.md`. Original note:
      The transport arm is ChEBI-keyed and the metabolism arm KEGG-keyed, and
      they do not always meet on one node, so any transport ↔ catalysis join
      silently misses. Measured over the 3,356 metabolite entries in
      `kegg_data.json`: 191 InChIKey-skeleton groups (first 14 chars; 465 nodes —
      157 all-KEGG stereo/charge variants such as Cysteine C00736 vs L-Cysteine
      C00097, Glutamate C00302 vs L-/D-, 34 mixed KEGG + chebi/mnx) plus 36
      exact-name KEGG↔ChEBI pairs the skeleton misses. Nitrate is the hard case:
      `kegg.compound:C00244` (ChEBI 25545, MNXM732399) and `chebi:14654`
      (MNXM732398) share neither InChIKey nor MNX id — a ChEBI conjugate pair —
      so it needs ChEBI's conjugate_acid/base relations, not string matching.
      Proposal: a `Metabolite_equivalent_to_metabolite` edge (symmetric, emitted
      once per pair) carrying `basis ∈ {inchikey_skeleton, chebi_conjugate,
      name}` and a `Metabolite.equivalence_group` id, computed in step 6 where
      the MNX resolver is open; consumers join on the group, never on the raw
      id. Prerequisite for the TCDB plausibility item below. → friction log
      2026-10-07 "Metabolite IDs split across namespaces", "The two nitrate
      nodes…", 2026-10-11 "Generic vs stereo-specific…";
      `multiomics_kg/download/build_kegg_metabolism_xrefs.py`

- [ ] **Per-organism most-specific transport reach on `Metabolite`** (explorer
      ask, 2026-10-13 — KG half of the "flat `list_metabolites` reach" item).
      Add post-import `transport_most_specific_organism_names: str[]` (+
      `_count`, size invariant) on `Metabolite`: organisms with ≥ 1 gene →
      `TcdbFamily` → metabolite path whose `Tcdb_family_transports_metabolite`
      edge has `substrate_depth = 'most_specific'`. Today `organism_names`
      unions the inherited superfamily (ABC `3.A.1`) substrates, so transport
      reach is flat — 1,294–1,304 of 1,462 transported metabolites for every
      organism; most-specific-only spreads from 639 (MED4) to 1,201
      (Pseudomonas KT2440). Design note: decide whether the gene hop also
      restricts to `Gene_has_tcdb_family.attachment_depth = 'most_specific'`
      (the projection `transported_metabolite_count` uses) and document which;
      the two filters are on different edges and compose. Sparse/`[]` rule to
      match `organism_names`. Explorer use: `list_metabolites` `top_organisms`
      and the organism filter for the transport arm. → `scripts/post-import.sh`
      + `.cypher` (Metabolite rollup block), `controlled_vocabularies.yaml` if
      registered; CLAUDE.md Metabolite bullet

**P2 — needs a design pass**

- [ ] **Per-arm sampling day + growth state, and a shared-model identifier for
      linked contrasts.** A paperconfig analysis carries one `timepoint`, so a
      contrast whose arms were sampled on different days gets one label:
      Weissberg 2025 `DE_med4_rnaseq_coculture_vs_axenic_day_18` pairs coculture
      d18 with axenic d14, and the KG labels that axenic arm `nutrient_limited`
      while the researcher reports it in full decline. The four MED4 RNA
      contrasts are also one fitted model over shared samples (coculture d18 −
      axenic d14 ≡ direct d18 − direct d11 to 1e-15), which no property exposes.
      Proposal: optional per-analysis `treatment_timepoint` / `control_timepoint`
      (+ `_hours`, + per-arm `growth_phase`) materialised as edge properties
      alongside the existing `time_point`, and an Experiment-level
      `shared_model_id` (free string, same value across experiments fitted
      together) so "independent contrast" is a queryable claim. → friction log
      2026-10-11 "Anchor RNA contrasts are algebraically linked", 2026-10-11
      researcher update under "The Plan phase missed three…"; `omics_adapter.py`

- [x] **Chemical categories for `Metabolite` (ChEBI is_a / has_role).** CLOSED 2026-10-14 — delivered as `MetaboliteClass` (KEGG BRITE compound trees) + `ChemicalRole` for TCDB role substrates; see `docs/superpowers/plans/2026-10-14-metabolite-identity-and-classes.md`. ChEBI `has_role` on `Metabolite` stays deferred (item below). Original note: The KG
      has no compound-class layer; KEGG pathway membership and the carrying TCDB
      family were tried as proxies and ~half the sampled assignments were wrong
      (thyroxine → "amino acids", GDP-mannose → "xenobiotic"). 3,042 of 3,356
      metabolites (91%) carry a `chebi_id`, so a GO-pattern pruned closure over
      `chebi.obo` (`is_a` + `has_role`, seeds = our metabolites, ancestors only)
      would give `ChemicalClass` nodes + `Metabolite_is_a_chemical_class` /
      `Metabolite_has_chemical_role` + `Chemical_class_is_a_chemical_class`,
      answering "Lys-Arg is a dipeptide", "ferrichrome is a siderophore". Depends
      on the equivalence item so a class attached to the KEGG node reaches its
      ChEBI twin. Explicitly flagged "revisit the analysis when it lands".
      → friction log 2026-10-08 "KG has no chemical categories" (OPEN);
      `docs/kg-changes/metabolomics-extension.md` (LipidMaps note)

- [ ] **Deferred TCDB spec (follow-up to metabolite identity + classes).** Bundle: (a) per-gene substrate plausibility (`catalysis_support`, now feasible on canonical identity) — see item below; (b) re-examine `substrate_depth` semantics after the identity merge; (c) per-organism most-specific transport reach (item above); (d) TCDB name hygiene (bare-id subfamily names, role-typed entries); (e) re-tune the lumping threshold and **retire the R3 exception** (`TcdbFamily.substrate_breadth` is a stored value restating a pinned threshold; sunsets with this spec). → `docs/kg-changes/vocabulary-contract.md`, `docs/kg-changes/tcdb-two-source-upgrade.md` §4
      (f) **Transport-reachable pathways (2026-10-15).** Step 6 adds the KEGG pathways of substrates rolled up onto ANY kept TC node (gene-annotated or ancestor-only), so re-keyed eukaryote lipids (PI3P `C04549`, PIP2 `C04637`, transported only by `tcdb:9.A` / `tcdb:9`, 0 direct genes) pull in `ko04145` Phagosome / `ko04510` Focal adhesion (+ `ko01521`, `ko05212`, `ko05332` via other re-keys; pathway KeggTerms 447 → 452, zero-gene 61 → 66). Consider requiring a gene-attached TC node (option A); `test_metabolite_pathways_only_kg_evidenced` was widened to the current rule. → `docs/kg-changes/metabolite-identity-and-classes.md`

- [ ] **TCDB substrate plausibility per gene.** Even at
      `substrate_depth = most_specific` a gene inherits every substrate its TCDB
      subfamily lists: MED4 sulfate permeases `sul1`/`sul3` list amino acids and
      nitrate, `ktrA`/`nhaS` list ammonium, the cyn system (3.A.1.16.1) reads as
      a nitrate/nitrite porter in a strain with no nitrate assimilation. Not a
      bug (NitT/TauT is genuinely bispecific; `transport_substrate_resolution`
      already marks the lumping case) — an annotation-granularity limit. Two
      candidate signals, both post-import and both cheap once equivalence
      exists: (a) `catalysis_support ∈ {corroborated, uncorroborated}` on the
      (gene, substrate) pair = the same organism has a `Gene_catalyzes_reaction`
      touching an equivalent metabolite; (b) a `substrate_breadth` int on
      `Gene_has_tcdb_family` (= `metabolite_count` of the attachment, so the
      explorer can rank "dedicated" vs "broad listing" without keyword rules).
      Second occurrence after round 1's ABC coarseness. → friction log
      2026-10-06 "TCDB most-specific substrate for MED4's cyn system",
      2026-10-08 "TCDB substrate lists describe the subfamily";
      `docs/kg-changes/tcdb-two-source-upgrade.md`

- [ ] **ChEBI roles on `Metabolite` (curated allow-list) — measured 2026-10-14, deferred from
      the metabolite identity + classes spec.** ChEBI `has_role` is a facet ("what it does / is
      used for"), not a class, and the step-6 ChEBI cache will already hold it. Raw: 1,562 of
      2,985 ChEBI-linked metabolites carry ≥ 1 direct role (508 distinct), but 2,269 assignments
      are organism-observation roles ("mouse / human / *E. coli* / yeast metabolite") and much of
      the rest is pharmacology. A rule-based filter (biological + chemical role, minus the
      "metabolite" subtree) still leaves 578 metabolites (19%) / 328 roles, dominated by drug
      roles that ChEBI files under *biological* role (antibacterial 81, anti-HIV 41, analgesic 38;
      CO₂ → "analgesic", ammonia → "NMR chemical shift reference compound"). A useful version
      needs a **curated allow-list of ~15–20 roles** (cofactor 74 rolled up, coenzyme 12,
      prosthetic group, antioxidant 52, nutrient / micronutrient 33 / 26, toxin 15, signalling
      molecule 10, siderophore 8, osmolyte 8, xenobiotic, pigment). Shape:
      `Metabolite_has_chemical_role` edges → the existing `ChemicalRole` nodes (introduced by
      the identity spec for TCDB role-typed substrates), direct assertions only,
      `sources: ['chebi']`, `evidence: 'curated'`. This is a different edge type from
      `Tcdb_family_transports_chemical_role`: "has the role" vs "transports something with the
      role", joinable on one node.
      → spec `docs/superpowers/specs/2026-10-13-metabolite-identity-and-classes-design.md` §13

**P3 — large or blocked on researcher input**

- [ ] **Culture-level measurement node type.** Cell counts, N biomass, DON and
      fluorescence have no target entity (every measurement edge points at a
      Gene or a Metabolite). Needs a node off Experiment (organism, condition,
      culture/bottle id, time, growth phase, value, unit,
      `value_kind ∈ {measured, derived}`, extraction method + uncertainty since
      most papers publish growth only in figures) plus a paperconfig entry type,
      and a sample-level link from omics samples to their culture measurements.
      Medium NH₄⁺/NO₃⁻ could ride `metabolite_assays_table` today **except** that
      `compartment: spent_medium` was pruned from the vocabulary on 2026-08-29
      (DOC-004) and must be re-added in both `controlled_vocabularies.yaml` and
      `vocab/non_de_evidence.py`. Blocked on the researcher: measured vs derived,
      the NH₄⁺ rerun reconciliation, fractional days vs omics day labels, and
      whether the Weissberg CC1A3 data is published or loaded as unpublished.
      → friction log 2026-10-12 "The anchor study's culture measurements are not
      in the KG"; `docs/kg-changes/metabolomics-extension.md` (advises a separate
      type, not an extension of `MetaboliteAssay`)

- [ ] **TCDB name/depth hygiene (minor).** `tcdb:3.A.1.4.4` is named
      "The high-affinity (" — the only name in `tcdb_names.json` ending in "("
      (citation-strip bug in `scrape_tcdb_names.py`'s first-sentence cut; fix the
      regex, re-scrape that family with `--families`). PepT (`dpp`) depth looks
      inverted (generic "peptide"/"Oligopeptide" only `inherited`, heme/EDTA/
      bradykinin `most_specific`) — verify against `subtree_substrates` before
      touching `_compute_substrate_depth`. Not bugs, document only: bare-id names
      on `4.C.1.1`/`1.A.1.x`/`2.A.88.x` are the designed fallback; `urtE`'s
      score-0 `1.B.42` is a tier-3 diamond call; inherited rows reporting the
      gene's own attachment as "family" is the explorer's projection, not the
      edge. → friction log 2026-10-07 "TCDB annotation oddities"

**Reframed — no KG change**

- KO coverage of transporters: `amt1` (PMM0263) carries K03320 and TC 1.A.11;
  what it lacks is a KEGG *pathway*, because KEGG files transporters under
  BRITE (`ko02000`), not pathway maps. Use `Kegg_term_in_brite_category` for the
  transporter ↔ enzyme link; a doc line in `docs://analysis/metabolites` is the
  fix. No KEGG reactions on di/tripeptides is upstream KEGG.
  → friction log 2026-10-08 open items 6–7

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
