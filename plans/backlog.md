# Backlog

Deferred and cross-cutting work that has no home in a single plan file.

**Convention.** One bullet per item: *what*, *why it was deferred*, and a pointer
to the spec section, plan file, or code comment that created the deferral. Add
the bullet in the **same commit** that creates the deferral — that is the whole
discipline, and the reason this file lives beside the code rather than in an
issue tracker (the GitHub issues opened in Jan–Feb 2026 were all completed and
none were ever closed).

Remove a bullet when the work lands; the CHANGELOG records that it did. (Last sweep
2026-11-11: closed and superseded items removed, `plans/interpro_redesign_backlog.md` merged in.)
Items that are one focused change get their own `plans/<topic>.md` when picked
up — this file is the index, not the plan.

---

## KG semantics

- [ ] **Eukaryotic KEGG pathways kept by the category rule (2026-11-01).** The KEGG refresh
      brought in ko04145 Phagocytosis, ko04149 Peroxisome - yeast and ko04156 Integrated stress
      response signaling; KEGG files them under Cellular Processes / Environmental Information
      Processing, which `config/kegg_pathway_exclusions.yaml` keeps (as it already keeps Endocytosis,
      Peroxisome, Autophagy - animal, Lysosome biogenesis). Decide whether to exclude subcategory
      "Transport and catabolism" eukaryotic maps by id, then rebuild step 6.

- [ ] **TCDB FASTA lags the TCDB tables (2026-11-01).** `/public/tcdb` was still the 2025-08-04
      build, so systems TCDB added since then are in the hierarchy (range-only) but no gene can be
      called onto them. Re-run `/tcdb-diamond --refresh-tcdb --force` + steps 2 6 once
      `last-modified` moves.

- [ ] **Step 6 `--force` is not reproducible on a checkout without its raw caches (found
      2026-11-11).** CLAUDE.md says `--force` rebuilds "from already-cached raw inputs (~1 min, no
      network)", but on a checkout whose gitignored raws (`cache/data/chebi/raw/chebi.obo`, the
      TCDB TSVs, `cache/data/kegg/raw/br_br08*.json`) are missing or stale, it silently downloads
      the current upstream release: a run here moved ChEBI 256 → 257, rewrote ~4.5K lines of
      `tcdb_hierarchy.json` (e.g. family 1.A.123 gone upstream) and shifted KEGG pathway lists —
      a data refresh nobody asked for (reverted). Options: fail (or require `--refetch-raw`) when a
      raw input is missing, or pin each raw to the release recorded in its committed output.
      → `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`_load_chebi`, TCDB step 1, BRITE)

- [ ] **Re-scan all strains with InterProScan 110 (2026-10-31).** The calls.json
      artifacts come from InterProScan 5.78-109.0, but the reference was refreshed
      to InterPro 110.0. The adapter now drops the 38 accessions 110.0 retired
      (1,665 gene edges, 35 auto-named PANTHER families), and accessions 110.0
      added (+854) are never called. Re-scanning (new image + data under
      `~/tools/InterProScan/`, `/interproscan-run --force` per strain, then
      `prepare_data.sh --steps 9 2 --force`) restores both. Heavy: many hours
      over 118 strains — run one strain at a time (OOM).

- [ ] **Metabolite identity review lists (2026-10-31).** Two lists in
      `cache/data/metabolomics/metabolite_identity_report.json` are written every step-6
      run but nobody has worked through them. `alias_stereo_review` (3): curated
      `metabolite_aliases.yaml` `stereo_assignment` values that differ from the name rule —
      biller 2022 "Adenosyl Methionine", Kujawinski 2023 "2,3-dihydroxypropane-1-sulfonate"
      and "n-acetyl muramic acid" (curated `dominant_form_assumed`, rule `as_reported`);
      keep or drop each override. `anomer_xref_review` (18): HMDB / ChEBI ids MNX files
      under one anomer of a generic sugar; curate with `xref:` overrides in
      `config/metabolite_identity_overrides.yaml` where wrong. Then `--steps 6 --force`.

- [ ] **43 Cyanorak-only genes still have no position (2026-10-28).** The neighbour-based
      coordinate mapping leaves genes that sit at a layout break (where Cyanorak's genome is
      rotated, reverse-complemented or re-ordered against NCBI's) unplaced: CB0101 17, CB0205 11 and
      15 elsewhere. They have `sequence` but no
      `contig`/`start`/`end`/`strand`, so `gene_neighbors` skips them. Could be placed by
      translating Cyanorak's protein against the NCBI genome. → `docs/kg-changes/cyanorak-coordinate-fallback.md`

- [ ] **Zero-count fold-change magnitudes (2026-10-31, documented, not changed).** Capovilla
      2023 (to ±21.5; only VST counts published), Anjur-Dietrich 2026 (two genes at +24 / +26;
      no count table) and Huang 2021 (unreplicated) keep extreme |log2FC| from zero-count genes,
      which dominate `rank_by_effect`. Fix needs raw count tables from the authors (then
      recompute with a one-read pseudocount, as for Beliaev 2014 / Aharonovich 2016).
      → `plans/signed_linear_fold_change.md`

- [ ] **WH8020 left out of the Ford 2021 strain batch (2026-10-28).** Cyanorak annotates a
      2.29 Mb assembly NCBI does not carry, so the NCBI ↔ Cyanorak gene-merge QC failed; its
      Ford 2021 rows are dropped. Revisit if NCBI gains that assembly or a Cyanorak-only genome
      path is added. (UH18301 has only SRA reads, no assembly: not actionable.)

- [ ] **MEROPS GO bridge — rejected on measurement (2026-08-18), revisit only
      with a concrete use case.** All-kingdom member rollup yields median 19 /
      max 389 GO terms per family incl. eukaryote-only terms; completeness win
      is 338 genes (~8%) vs the 1,311 that justified TCDB's GO bridge. Needs a
      filtering design (e.g. >= N supporting identifiers) before it can land.
      → `docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`

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
      old assemblies. The same measurement covers the alternative of folding NCBIfam
      GO into gene `go_terms` in the step-2 merge (InterPro Layer B pattern; today
      `enrich_interpro_fields` uses `ncbifam_reference.json` for naming only). NCBIfam
      **EC** propagation was never measured; do that first if this is revisited.

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

## InterPro / NCBIfam follow-ups

Merged in from `plans/interpro_redesign_backlog.md` (2026-11-11; that file is retired).
Origin: the 2026-08-17 multi-ontology redesign
(`docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`,
`docs/kg-changes/interpro-multi-ontology.md`); all triaged OK-to-defer. The NCBIfam→GO,
MetaCyc and multi-source NCBIfam items sit under KG semantics above.

- [ ] **antiSMASH BGC detection.** `/add-a-tool` candidate: siderophore / vitamin
      cross-feeding currency; per-genome region calls (a new calls.json shape, per the
      add-a-tool table).
- [ ] **Further member-DB elevation** (PANTHER, HAMAP-as-ontology, CDD, folds …). Each is one
      merge-config line + a small adapter under the faceted design. Only on demonstrated need.
- Not tasks, watch notes: the `acc_to_ipr` last-write-wins note in
  `ncbifam_adapter.get_edges` must be re-checked after the InterProScan 110 re-scan above; if
  several adapters ever need the same pruned InterPro tree, materialize a committed pruned
  artifact then (spec §3.2).

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


**P2 — needs a design pass**

- [ ] **Paperconfig items deferred from the control-timing review (2026-10-18).** Each needs
      something not available in the review: (a) **Barreto 2022 S12–S14 (800 ppm partner tables)**
      are probably the partner × pCO2 interaction deviations (author-script contrasts 7–9), not the
      deviation-at-800-ppm the paper describes (logFC anticorrelates with the 400 ppm tables,
      r = -0.61 / -0.83 / -0.97); settle with the authors or a GLM re-run from `ALL_EZ55_2.Rdata`
      (needs R). The experiments carry a caveat in `design_notes`. (b) **Aharonovich 2016**
      temperature / light (config 24 C, 10 uE vs Methods maintenance 20 C, 27 uE) — the
      Supplementary information is not in the paper directory. (c) **Anjur 2026** (was Anjur 2025)
      planktonic-vs-biofilm needs a new `treatment_type` value (e.g. `lifestyle`), a vocabulary
      change; it stays `growth_phase`. (d) **Lindell 2007** phage is the generic `Phage` (taxid
      10239); naming P-SSP7 needs a `treatment_organisms.csv` row and changes the OrganismTaxon
      count pinned by `tests/kg_validity/test_organism.py`. → `plans/control_timing_curation.md`

- [ ] **Two table rows → one gene in one analysis: BioCypher silently keeps one (found 2026-10-31,
      fold-change audit live check).** The omics adapter yields both edges with the same id
      (`{doi}_{analysis_id}_{gene}`) and BioCypher dedups them, so the graph holds one of two
      possibly different measurements, chosen by row order — and with it the edge's
      `log2_fold_change` and `significant`. 85 edges over 10 papers on the 2026-10-31 build: Coe 2024
      49 (39 of the dropped copies were significant; whether the pair disagrees is unchecked), Fang 2019 9, Biller 2018 6, Domínguez 2017 5, Ma 2022 5,
      Ferrer 2020 4, Thompson 2016 3, Al-Hosani 2015 2, He 2022 1, Kratzl 2024 1. Pre-existing (same
      counts in the 2026-10-28 snapshot). Detection DONE 2026-11-11: the adapter logs each
      collision (analysis, gene, the rows' identifiers, log2FC, significance, how many disagree) and
      `validate_paperconfig.py` warns on duplicate `resolved_locus_tag`s per table. Still open: decide
      per paper — a real duplicate row, two paralogs one ID resolves to both, or two probes per gene
      (needs an explicit rule, e.g. keep the larger |log2FC| or mark the edge). Related step-4 bug
      (review, 2026-11-11): `resolve_paper_ids` resolves a footnote row to a gene (Beliaev 2014
      `supp_table_s6`/`s8` row 1342 "1 RPKM is defined as …" → `Sputw3181_2756` via `multi: Locus
      tag`); harmless today (no logFC, the adapter skips it) but resolution should skip such rows.
      → `multiomics_kg/adapters/omics_adapter.py`,
      `<table>_resolved.csv` (`locus_tag`, `resolution_method`)

- [ ] **TCDB substrate plausibility per gene** (the one live part of the old "Deferred TCDB
      spec"; its `substrate_depth` / `substrate_breadth` / most-specific-reach / name-hygiene
      items were removed or fixed by 3a and the 2026-10-24/25 TCDB review). Since 3a a gene
      reaches only the substrates of the `tc_system` it hits directly, but a system's TCDB
      substrate list can still be broader than the gene's biology (MED4 cyn system 3.A.1.16.1
      reads as a nitrate/nitrite porter in a strain with no nitrate assimilation; NitT/TauT is
      genuinely bispecific). Candidate signal, post-import, feasible on canonical identity:
      `catalysis_support ∈ {corroborated, uncorroborated}` on the (gene, substrate) pair = the
      same organism has a `Gene_catalyzes_reaction` touching that metabolite. Probably belongs
      in the 3b design. → friction log 2026-10-06 "TCDB most-specific substrate for MED4's cyn
      system", `docs/kg-changes/tcdb-cleanup.md`

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

- [ ] **Metabolite identity follow-ups (deferred from the 2026-10-14/15 implementation reviews).**
      None blocked the merge (final review triaged all as "can wait"); grouped by what they touch.
      → plan `docs/superpowers/plans/2026-10-14-metabolite-identity-and-classes.md`;
      `docs/kg-changes/metabolite-identity-and-classes.md`
      - *Curation to confirm (researcher call).* cyanocobalamin → `C02823` overrides the Kujawinski
        export's `C08230` (hydroxocobalamin); chitotriose → `chebi:143141` (triacetyl form);
        "nicotinamide adenine dinucleotide" → NAD⁺ `C00003` while the export's "NAD pos" rows keep
        NADH `C00004`; Biller "Vitamin B3" → nicotinate `C00253` (Ingalls standard). All in the
        papers' `metabolite_aliases.yaml` with reasons; re-run step 6 after any change.
      - *Semantics.*
        - `names._GEOMETRIC` also strips ring cis/trans (relative configuration, i.e. diastereomers),
          so e.g. `C04314` cis-1,2-dihydronaphthalene-1,2-diol → `C06205` lost its stereo link;
          restrict the geometric rule to double-bond descriptors.
        - ~20 abstract KEGG compounds with no formula and no ChEBI class term (`C00028` Acceptor,
          Quinone, …) are `structure_kind = defined` by default; consider a curated list or a
          KEGG-name rule.
        - `member_basis` records the first link that touched a member, not the link that attached
          it to the canonical node (manual-first ordering can misattribute); a manual merge
          blocked transitively by `keep_separate` is silent (no report entry).
        - Remap fallback basis labels are approximate (conjugate / zwitterion relocations read
          `chebi_primary_id` / `rekey_chebi_kegg_xref`); `repoint_reason`'s `identity` / `rekeys`
          parameters are unused.
        - Committed ChEBI cache lacks `CHEBI:176566/176629/176584/10809`, so a cache-only rerun
          (raw obo unavailable) could relabel 3–4 remap rows.
        - `identifier_tokens` includes bare ChEBI numbers, so digit-only full-text queries match
          them — document in the explorer contract.
      - *Robustness.* Classes parser drops a compound leaf at a tree root silently; a substrate
        equal to the role root `CHEBI:50906` would crash `_level`.
      - *Test gaps.* `download_chebi_obo` / `download_compound_brite` / `metabolite_identity_tsvs.py`
        edge cases; stale `keep_separate`, 3+-member union, `member_basis`, tie-breaks;
        first-block-mismatch stereo link; boolean assay-edge `ValueError`; validator non-str key /
        missing `id` / cross-file `name_key` conflicts; ChemicalRole multi-parent is_a; the
        `tests/kg_validity/snapshot_data.json` re-sample dropped `data_source:merops_diamond`.
      - *Docs.* paperconfig SKILL: a single stereo-unspecified match with no ≥ 2× winner is also
        ambiguous; the direct-id path wording; prepare_data.sh step-6 `run_step` label understates
        its outputs.

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

**Reframed — no KG change**

- KO coverage of transporters: `amt1` (PMM0263) carries K03320 and TC 1.A.11;
  what it lacks is a KEGG *pathway*, because KEGG files transporters under
  BRITE (`ko02000`), not pathway maps. Use `Kegg_term_in_brite_category` for the
  transporter ↔ enzyme link; a doc line in `docs://analysis/metabolites` is the
  fix. No KEGG reactions on di/tripeptides is upstream KEGG.
  → friction log 2026-10-08 open items 6–7

## Transport programme (design `docs/superpowers/specs/2026-10-18-go-transport-substrates-design.md`)

- [ ] **Project 3: TCDB nodes/edges redesign** — split 2026-10-22 into three; 3a (TCDB cleanup) done 2026-10-23 → `docs/kg-changes/tcdb-cleanup.md`:
  - [ ] **3b gene transport layer** (needs its own design): entry via TCDB, KEGG BRITE, transport Pfam/InterPro domains, Cyanorak Q roles, GO; per-gene subunit role, purpose, direction, confidence; no per-gene compound edges. Answer keys: the two `multiomics_analysis` transport analyses (2026-10-06 N import, 2026-07-06 Alteromonas carbon) and the Ford 2021 SBP table (`data/Synechococcus/papers_and_supp/Ford 2021/`). POC + draft paper section: claude.ai artifact PKYdYsmr4YefKhCBmemAqA.
    - [ ] Reciprocal check recall cost (from project 2): ~2,990 gate-passing proteins fail it corpus-wide, ~1,300 hold-out ones with transporter-like products (paralogous ABC components in heterotrophs). Decide in the 3b design (e.g. family-level reciprocity) — measured, not tuned.
  - [ ] **Systems project** (after 3b): group genes into transport systems (genomic neighbourhood + subunit role + shared TCDB lineage, as both analyses did) and assign substrates per system from the binding protein; per-gene `Gene` → `Metabolite` transport edges belong here.
- [ ] **Catalysis-only pathway creation — needs scoping, not agreed** (raised 2026-10-23). Today a pathway enters the KG via catalysis (gene KO/reaction), transport (any TCDB substrate under kept nodes) or papers (measured compounds); the GO route already never creates pathways. 66 of 452 pathways have no gene (animal signalling via Ca²⁺/Na⁺/K⁺/Cl⁻/cAMP/DAG, antibiotic and plant secondary-metabolite biosynthesis via efflux substrates). Proposal: only catalysis creates pathways; transported/measured compounds still get `Metabolite_in_pathway` to existing ones. Open: is a measured compound's pathway useful context even with no gene? Overlaps the category filter. The 2026-10-15 Phagosome / Focal adhesion pull-in (Deferred TCDB spec item (f)) is resolved by 3a itself: dropping the class-level roll-up removes 65 off-target transport-only Metabolite nodes and exactly 3 gene-less pathways (Phagosome, Focal adhesion, chondroitin sulfate GAG biosynthesis); the other 63 gene-less pathways remained for this item. After the KEGG category filter (2026-10-24) dropped 32 of them, **31 of 267 pathways have no gene** (measured on the 3a graph): 16 Metabolism (antibiotic / plant / fungal secondary metabolites, steroid hormones, heparan sulfate, GPI anchors, bisphenol), 9 Environmental Information Processing → Signal transduction (MAPK, ErbB, Ras, mTOR, Wnt, Hedgehog, TGF-beta, VEGF, sphingolipid), 6 Cellular Processes (yeast cell cycle, cellular senescence, actin cytoskeleton; Transport and catabolism: autophagy other / animal, endocytosis).
- [ ] **TCDB substrate overrides — needs curation and scoping** (decision 2026-10-25: we do not silently fix TCDB's mistakes in code). The review found ~60 upstream errors on gene-carrying systems (findings: session scratchpad `tcdb_review/findings_substrate_biology.json` + `findings_resolver_identity.json`; re-derive if lost): wrong ChEBI ids (PstABCS 3.A.1.7.6 → triphosphate, 100 genes; DAPI given as nifedipine 3.A.1.106.5 / 2.A.7.1.11; TPP = tetraphenylphosphonium read as thiamine diphosphate 2.A.7.1.10/.11; AmtB methylammonium as tetramethylammonium; PhnCDE 3.A.1.9.1 → iron(3+); UDP-sugar for a free sugar on 3.A.1.2.22/.24/.25/.27, SemiSWEET 2.A.123.2.2; ThiXYZ 3.A.1.17.13 → riboflavin; PcaM 3.A.1.4.7; CHEBI:35274 "ammonium ion derivative" for ammonium; safranin O = phenosafranine; diaminopimelate → 7,8-diaminononanoate on CysXYZ 3.A.1.3.10; arsenopyrite on ArsAB 3.A.4.1.1), non-transported compounds (PSII O₂/plastoquinone 3.E.2.2.3; Mla tetracycline from mutant phenotype 3.A.1.27.3; fructose 6-phosphate on MurJ 2.A.66.4.x), and missing substrates (F-ATPase 3.A.2.1.10 none; tetracycline on TetA 2.A.1.2.4; Fe-pyochelin on FptA 1.B.14.1.8). Scope: a `config/tcdb_substrate_overrides.yaml` (replace / remove / add per system, reason + reference), reported in the identity report, and reported upstream to TCDB. Also ChEBI-side: "corrinoid" → C06021 "Co(I) corrinoid protein" via ChEBI's own KEGG xref (3.A.1.14.9).
- [ ] **Project 4: transport-analysis recipe** (`docs://analysis/transport`, explorer-owned).
- [ ] **Direction gap**: no project answers question 3's "in which direction"; GO target start/end location is too sparse (5 MED4 import annotations, one gene).
- [ ] **GO taxon-constraint flag**: 1,912 of 1.4M gene to GO annotations violate a GO taxon constraint; flag on the gene to GO edge.
- [ ] **GO MF to Rhea / EC / KEGG reaction xref bridge** (unmeasured).
- [ ] **Obsolete gene GO annotations** (pre-existing, found in Task 0): ~330,594 gene `go_terms` in the merged annotations reference ids not live in GO 2026-07-26 (139K obsolete with no replacement, 145K with a single `replaced_by`, 45K `consider`-only); alt-id annotations are remapped (23,462), the rest dropped at build, before and after the refresh. Source annotations (eggNOG, UniProt, Cyanorak) carry old ids; consider a `replaced_by` remap or refreshing the sources.

## Paper batch 2026-10-28 — tables that need a new concept

Deferred from the 2026-10-28 paper batch (Manck 2020, Tetu 2009, Barber 2026, Ferrer 2020,
Durham 2015, Rubin 2015, Oleza 2012, Ostrowski 2010, Oleza 2015b, Koch 2019, Koch 2020,
Doré 2020/2023, Schroer 2023). Rule applied: wire what fits existing concepts; a table that
would need a new KG concept waits here instead of being forced.

- [ ] **Cross-strain protein families (multi-organism clusters).** Manck 2020 Data Set S3:
      1,804 TonB-dependent transporters over 31 *Alteromonas* genomes in 98 sequence-similarity
      network clusters (user: TBDTs are important). `gene_clusters` / `ClusteringAnalysis` is
      single-organism, so one cluster spanning strains has no home; options: a multi-organism
      ClusteringAnalysis, or attaching to OrthologGroup-like family nodes. In-KG strains: ATCC27126,
      AD45, BS11, AltDE1 resolve today; AltDE (`MADE_00089`, our assembly has only `MADE_RS`) and
      HOT1A3 (IMG `Ga0125132_`) need diamond `id_translation`s. Would make EC673 + the other 24
      genomes worth onboarding. Same concept covers Barber 2026's TRAP sequence-similarity-network
      sheet (Supplementary Data 3). → `data/Prochlorococcus/papers_and_supp/Manck 2020/paperconfig.yaml`
- [ ] **Curated transporter → substrate tables (Schroer 2023).** Table S1: 128 *R. pomeroyi*
      DSS-3 organic-carbon influx transporters, each a multi-locus unit (`SPO0050; SPO0051; ...`)
      with a mutant-confirmed substrate. Needs a transporter-unit concept and a substrate link to
      `Metabolite` (not a free-text category). Also the natural answer key for transport 3b
      validation. (Ford 2021's SBP classes shipped as categorical metrics — this is the
      substrate-linked, experimentally confirmed version.) → `plans/transport_validation_publications.md`
- [ ] **Unnamed natural-community partner + seawater-origin treatment (Oleza 2012).** DSS-3
      whole-cell proteome, 20 h in natural seawater from 4 Marseille sites, Table S2 panels A-J
      (PatternLab AC-test, p < 0.01 and |FC| ≥ 2). A-D = natural community present vs 0.22 um
      filtered (needs a treatment organism for "natural seawater microbial community"; H-J need it
      as background); E-J = beach / ferry basin / marina vs cape seawater (anthropogenic nutrients
      + copper; no `treatment_type` fits, `chemical` is a stretch). Parsed CSV is ready
      (`oleza2012_table_s2_wide.csv`, builder `build_oleza2012_table_s2.py`; 10 contrasts, log2).
      Once both concepts exist: 4 experiments, PROTEOMICS, `control_timing: matched`. Tables S3
      (30 one-shot conditions) and S4 (fractions incl. exoproteome) also unintegrated.
- [ ] **Ortholog-group-level genomics (Doré 2020, Doré 2023).** Core/accessory genes, genomic
      islands, phyletic patterns and eCAG modules keyed by Cyanorak CLOG cluster x genome. The KG
      has no attachment at OrthologGroup level for a paper's classification.
- [ ] **Formula-level metabolite features (Koch 2019, Koch 2020).** FT-ICR-MS exometabolomes
      give molecular formulas (thousands of masses per strain/phase) with no compound identity;
      `metabolite_assays_table` needs a resolvable metabolite. Needs a formula-feature concept.
- [ ] **Non-gene genomic features (Rubin 2015).** sd04 ncRNA essentiality, sd05 upstream
      regulatory-region essentiality, sd06 essential intergenic regions: no Gene nodes to attach to.
- [ ] **Bacteroidetes ortholog level (Ferrer 2020 *Polaribacter* MED152).** Full RNA-seq table in
      coculture with *T. pseudonana*; MED152 would join only at Bacteria-level COGs (no family
      level in `ORGANISM_GROUP_LEVELS`). Revisit with a second flavobacterial paper/strain.
      (*Stenotrophomonas* SKA14: skipped, off-theme.)
- [ ] **Knocked-out gene link for `mutant` experiments.** Tetu 2009 (phoB SYNW0947, phoR SYNW0948)
      and Ostrowski 2010 (ptrA SYNW1019) name the deleted gene only in experiment text; an
      Experiment → Gene "knocks out" edge would make mutants queryable by gene.

Smaller ID follow-ups from the same batch (not concepts, just mapping work):

- [ ] Rubin 2015: ~57 PCC7942 chromosomal genes have no bridge from the old GCF_000012525 tags to
      the GCF_014698905 (`H6G84_`) assembly — needs a protein-sequence `id_translation`.
- [ ] Ferrer 2020: 48 DSS-3 SPO tags absent from the current annotation; Oleza 2015b: 14 WH7803
      YP_ accessions skipped (paralog ambiguity).
- [ ] IPG legacy-accession bridge gaps: WH8102 NP_ maps only to the GCA `CAE` accession
      (`find_assembly_match` vs the GCF returns None) and `ZP_` is not matched by
      `looks_like_legacy_protein_acc` — Oleza 2015b used diamond instead. **The two are one fix**
      (measured 2026-11-11): IPG returns retired `ZP_` accessions with only an INSDC row on the
      paired GenBank assembly (`ZP_01467679.1` → `EAU71572.1` on `GCA_000153805.1`, which
      `genomic_gca.gff` anchors to `BL107_12230`, the gene diamond chose), so adding `ZP_` to the
      regex alone does nothing. Accepting the same-number `GCA_` as the target assembly is not
      small: over the cached IPG bodies it would newly match up to ~3.7K accessions on 34
      strains (WH7803 1,161, DSS-3 861, BS11 137; upper bound — the cache is not keyed by
      strain), re-resolving paper rows through each strain's GenBank annotation. Needs a
      measured step 3 + 4 run (`/omics-edge-snapshot` before/after). No paperconfig types a
      `ZP_` column `protein_id_refseq` today.
      → `multiomics_kg/utils/ncbi_protein_xref.py` (`LEGACY_PROTEIN_ACC_RE`, `find_assembly_match`)
- [ ] Oleza 2015b and Oleza 2017: PatternLab TFold raw p-values are mapped to
      `adjusted_p_value_col` (open policy item).

## Explorer / MCP coordination

- [ ] **Explorer: treat `structure_kind = 'generic'` substrate targets as a class of compounds.** TCDB lists ChEBI structural classes as substrates (alpha-amino acid → C00045 on 76 systems, peptide, lipid, fatty acid, nucleoside …); they stay generic `Metabolite` nodes by design (classes, not placeholders), so a transport answer that lands on one should say "a class of compounds", not a single molecule. From the TCDB integration review (2026-10-25). Raised with the explorer as E1 of `multiomics_explorer/docs/kg-specs/2026-10-29-transport-small-followups-explorer-asks.md`.

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
