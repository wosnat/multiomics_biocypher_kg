# Gene groups in the KG — concept scoping (pre-design)

Status: **scoping capture ahead of a real design — this is the intended next
major step in the KG's evolution, expected soon.** No schema decisions are
made here; a full design spec (brainstorming → approaches → spec → plan) will
follow. Meanwhile the GEO-branch work (`plans/geo_paperconfig_updates.md`)
proceeds unaffected — its D4/D5 decisions (encode islands and iModulons as
`gene_clusters`) stand, and are explicitly reversible (graph rebuilds from
source).

Trigger: scoping antiSMASH via `/add-a-tool` surfaced that BGCs are the first
*region*-shaped tool output, which opened the broader question: the KG is
accumulating "group of genes" data from many directions, and its only current
answer is GeneCluster.

---

## 1. The concept, split correctly

"Groups of genes" is three distinct concepts. Conflating them is the main
design risk; the split is the main organizing device.

| Concept | Definition | Membership is | Examples |
|---|---|---|---|
| **Region** | contiguous interval on a contig | positional (coords) | genomic island, prophage, BGC region, integron, CRISPR array, plasmid/replicon |
| **System** | multi-gene functional machine with component roles | functional (role-typed) | ABC transporter system, secretion system (T6SS), defense system, CRISPR-Cas, metal-resistance operon cassette |
| **Set** | functional grouping, no spatial claim | analytical | regulon, iModulon, KEGG module, co-expression cluster (= existing `GeneCluster`) |

Some data sits across categories: a BGC is *detected* as a region but *means*
a biosynthetic system; an operon is a region that often IS a system.

## 2. Instance vs type — the two-layer anatomy

Every group source has the same anatomy, and the field already thinks this way:

- **Instance**: one detected group in one strain (MED4's island at 428–487 kb;
  HOT1A3's plasmid pAM1A3; antiSMASH region 2 of ATCC27126). Carries coords
  (when spatial), tool/paper provenance, score.
- **Type / identity**: the cross-strain identity instances point at (island
  "ISL1" recurring across all HL strains; fGI-2 recurring across
  *A. mediterranea* genomes *and* on a HOT1A3 plasmid; BGC product class
  "lanthipeptide"; defense system type "CBASS").

Evidence the type layer is real and published, not our invention:

- Doré 2020's island network modules (`Pro_GI033` = MED4 ISL1, `Pro_GI004` =
  ISL4, …) reconcile automated detection with the canonical Coleman 2006 /
  Dufresne 2008 island names — a published instance→type mapping.
- fGI-2 (Fadeev 2016) exists on *A. mediterranea* chromosomes AND on a HOT1A3
  plasmid — so type identity must NOT be hard-wired to chromosomal
  coordinates, and must survive crossing replicon kinds.

KG precedent for the pattern: `ClusteringAnalysis`/`GeneCluster` (instance
groups + per-analysis provenance), `Experiment` (instance) vs treatment_type
vocab (type), and query-time ortholog-overlap for cross-strain comparison.

### 2b. Synteny — the fabric under the type layer

Synteny (conserved gene order/neighborhood across genomes) is special: it is
both a **layer** (synteny blocks = paired regions between two genomes) and
the **mechanism** the type layer is computed from — "MED4's ISL1 corresponds
to MIT9312's ISL1" is exactly an ortholog-anchored synteny statement, and
Doré 2020's island network (gene-content similarity between islands) is a
synteny-flavored computation.

Two important observations:

- **The KG already holds the primitives.** Gene coordinates + the
  `gene_org_contig_start_idx` genomic-window index (backing `gene_neighbors`)
  + `OrthologGroup` membership together make *query-time* synteny possible
  today: neighbors of gene → ortholog groups → those groups' members'
  neighbors in another strain. No external tool required for the
  ortholog-anchored form.
- **The design choice is whether to materialize.** Options range from
  (a) leave it query-time (status quo, works for spot questions, too slow /
  too implicit for "which strains share this neighborhood" at scale), through
  (b) materialize pairwise synteny-block instances (heavy: O(strains²)
  blocks), to (c) materialize only *anchored correspondences between group
  instances* — i.e., use synteny to compute the instance→type edges for
  islands/BGCs/operons and store those, not the whole synteny fabric.
  (c) is the shape the driver papers already use.

External tools if ever needed (progressiveMauve, MCScanX, Sibelia, clinker
for BGC-level synteny plots) — but the ortholog-anchored in-house derivation
is likely sufficient and stays consistent with existing KG semantics.

## 3. Inventory — tools

| Layer | Concept | Tool(s) | Input | Fit to this KG | Cost |
|---|---|---|---|---|---|
| BGCs | region/system | **antiSMASH** (scoped, see §5) | genomic.gbff | High — interaction mechanisms, metabolomics bridge | ~10 GB Docker image, overnight batch |
| Operons / TUs | region | Operon-mapper; or distance+strand heuristic refined by our own expression data | GFF (+RNA-seq) | Foundational — the grain everything else decomposes into | Low compute; prediction quality is the risk |
| Genomic islands | region | IslandViewer/IslandPath-DIMOB — **but published data (Hackl, Doré) likely supersedes running a tool** | gbff | High — Pro island biology is core | Low |
| Prophages / IS / integrons | region | geNomad, ISEScan, IntegronFinder | genome.fna | Medium-high — viral experiments in KG | Low-medium |
| Defense systems | system | DefenseFinder or PADLOC | protein.faa | High — pairs with viral experiments + island cargo | **Cheap** (HMM scan, minutes/strain) |
| Secretion systems / pili / flagella | system | MacSyFinder2 + TXSScan | protein.faa | Medium-high — candidate Alteromonas–Pro interaction machinery | Cheap |
| Transport systems (complexes) | system | none — derivable from existing TCDB edges + operon adjacency; KEGG modules cover many ABC systems | existing KG | Medium | Zero external tooling |
| KEGG modules | set | none — KEGG definition file over KOs we already have | existing annotations | Medium-high (module completeness per organism) | **Cheapest of all** |
| Regulons | set | FIMO + RegPrecise PSSMs, or curated transfer; or data-driven from our own DE | genome + motifs | Very high biological fit (NtcA ↔ N-limitation experiments) | **Highest risk** (noisy) — own scoping later |
| Replicon typing (chromosome vs plasmid) | region (container) | none — derivable from assembly metadata (contig → replicon kind) | existing cache | High leverage per unit cost — "is this gene plasmid-borne?" | Near zero |
| Synteny | cross-genome fabric (see §2b) | none needed for ortholog-anchored form (in-house from coords + OrthologGroups); progressiveMauve/MCScanX/Sibelia/clinker if alignment-grade needed | existing KG | The mechanism behind cross-strain group identity; rarely a user-facing layer itself | Zero external tooling for form (c) |

## 4. Inventory — papers (driver datasets + ground truth)

**Pro/Syn islands & mobile elements**

- **Hackl 2023** (tycheposons) — 5,598 predicted islands, 623 genomes; 11 of
  our strains verified coordinate-compatible (~139 islands, ~6,000 membership
  edges). *Already staged on branch `data/geo-processed-supplements`*, planned
  as `gene_clusters` per that plan's D4. Also the mobile-element layer
  (tycheposons themselves).
- **Doré et al. 2020**, Front. Microbiol. 11:567431 (PDF in hand at repo root)
  — 81 genomes (28 Pro, 53 Syn/Cyanobium); per-strain islands (Supp Table S5),
  cross-strain island network modules `Pro_GI###`/`Syn_GI###` (Supp Table S6)
  mapped to canonical ISL/SVR names; CLOG pangenome (core 911 genes); clade-
  specific niche-adaptation clusters (cynA-B-D cyanate, nfeD-floT flotillin).
  **Second, independent island source over overlapping strains → cross-source
  corroboration; and the published instance→type mapping.** Supp tables not
  yet fetched/verified.
- **Coleman et al. 2006** (Science) — original MED4/MIT9312 island definition
  (ISL1–5); naming ground truth.
- **Dufresne et al. 2008** (Genome Biol) — tetranucleotide-based islands
  across Pro+Syn; naming ground truth.
- **Kettler et al. 2007** (PLoS Genet) — gained-gene island definition Doré's
  method extends.

**Alteromonas fGIs, plasmids, metal-resistance cargo**

- **López-Pérez, Ramon-Marco & Rodriguez-Valera 2017**, BMC Genomics 18:36
  (doi 10.1186/s12864-016-3461-0) — conjugative elements + plasmid networks
  across *Alteromonas*; plasmid cargo cassettes = the same blocks that
  populate chromosomal fGIs. The lab that defined the fGI concept.
- **Fadeev et al. 2016**, Front. Microbiol. 7:248 — **HOT1A3** (deployed KG
  strain, coculture partner) plasmid pAM1A3 (148 kb) carries fGI-2 (otherwise
  known only from *A. mediterranea* = AltDE/AltDE1, also deployed), with
  mer/cop/czc×2 metal-resistance operons + NiFe hydrogenase cluster.
- **Cusick et al. 2020**, AEM 86:e01831-19 — *A. macleodii* megaplasmids
  (~200 kb) with copper/mercury/zinc-cadmium resistance clusters, cop operons,
  conjugation machinery. (Strains CUKW/KCC02 — not deployed; concept-relevant.)

Ecology hook: *Prochlorococcus* is notoriously copper-sensitive; metal-
resistance cargo on *Alteromonas* mobile elements is interaction-relevant
biology, not incidental.

**Sets (already on the GEO branch)**

- **Johnson 2026b** — 32 MED4 iModulons (ICA regulatory modules), planned as
  `gene_clusters` per the GEO plan's D2/D5; regulator layer (TF_regs, RpaB
  binding sites) deferred to a Tier-3 spec there.
- **Steglich 2010** — 12 RNA-decay clusters; **Voigt 2014** — TSS
  architecture (a real TSS/operon/UTR entity layer is a named Tier-3 spec on
  the GEO plan, serving Voigt + Steglich + Doron together).

## 5. antiSMASH Phase-1 feasibility (scoped 2026-08-18, not started)

- Predicts: per-genome BGC regions (product class, member genes,
  KnownClusterBlast/MIBiG similarity). First region-shaped calls.json (list,
  not WP_-keyed dict) — the add-a-tool skill's key-choice table anticipates it.
- Input: `genomic.gbff`, cached for all 42 strains. Using gbff (not fna) means
  antiSMASH reuses NCBI gene calls → regions come out with **our locus_tags**,
  zero ID-mapping work (`--genefinding-tool none`).
- Install: Flavor A, `antismash/standalone` Docker (~9 GB download, all DBs
  bundled), version line 8.0.x. Disk verified: 65 GB free; ~65–70 GB more
  reclaimable via docker prune (avoid `docker volume prune` — Neo4j volume).
- Recommended module set: KnownClusterBlast ON (the high-value MIBiG
  comparison), generic ClusterBlast/SubClusterBlast OFF (slow, low value).
- Wallclock: ~5–30 min/strain → overnight 42-strain batch.
- Expected yield: sparse and concentrated — Pro HL strains ~1–3 regions
  (mostly the ubiquitous terpene cluster); MIT9313/MIT9303 lanthipeptides
  (prochlorosins — LanM machinery confirmed present in merged annotations);
  Alteromonas + heterotrophs richest (~5–12 each; ATCC27126 shows 41
  "siderophore" + 22 "polyketide" gene-level mentions with no cluster
  structure connecting them). Total order ~100–300 regions.
- Keep-all policy in calls.json (raw evidence, re-thresholdable) — consistent
  with the repo's "evidence, never a filter" convention.
- Phase 2 would be instance nodes (GeneCluster-like), NOT an ontology — the
  integrate-a-tool ontology/gene-property fork doesn't fit cleanly. Flagged,
  not designed.

## 6. The consumer angle — what groups mean to a KG user

*(The KG feeds LLM agents and humans doing omics interpretation. Groups only
earn their place if they make answers better, not just the graph bigger.)*

### What a user actually does with groups

1. **Context on a gene** (the most common read): "gene X sits in genomic
   island ISL2 / in a siderophore BGC / on plasmid pAM1A3 / in 3 iModulons."
   One hop, changes interpretation immediately — island/plasmid genes are
   flexible-genome, HGT-prone, often stress-responsive; BGC membership names
   a product; module membership names co-regulation.
2. **Group-level analysis**: "is this island/BGC/module coordinately DE in
   this experiment?" — the group is the statistical unit, like pathway
   enrichment but over spatial/regulatory groupings. Groups + existing
   per-gene expression edges = free co-regulation analyses.
3. **Cross-strain comparison**: "which strains carry fGI-2 / a lanthipeptide
   BGC / this defense system?" — needs the type layer; instances alone can't
   answer it.
4. **Narrative units**: papers discuss islands, BGCs and operons — not gene
   lists. Groups are the mid-scale units LLM agents need to tell the story a
   paper tells. (`Publication_discusses_*` edges could eventually target
   group types.)

### How not to be overwhelmed by multiple/conflicting groupings

Principles proposed (to be tested in the eventual design):

- **A grouping is a claim by a source, never truth.** Every group instance
  carries provenance (tool+version or paper), exactly as expression edges
  belong to Experiments. Two sources disagreeing on island boundaries is
  *data*, not a bug to resolve at ingest.
- **Corroborate, don't merge.** Keep each source's instances separate; add
  computed correspondence (overlap/containment) between them, and
  corroboration status on membership (`both_sources`/`single_source` — the
  existing TCDB two-source vocabulary). Mirrors "count, don't combine"
  (InterPro) and the hierarchical-agreement lesson (TCDB): near-miss overlap
  is still corroboration.
- **Overlap across kinds is not conflict.** An operon inside a BGC inside an
  island on a plasmid is *expected nesting*, not disagreement. Real conflict
  exists only within a kind (two island predictions for the same strain).
  The region/system/set taxonomy is the routing device that prevents the
  "15 parallel groupings" feeling: the user picks a lens by their question —
  *Where did it come from?* → island/mobile/replicon. *What machine is it
  part of?* → system. *Who is it co-regulated with?* → module/regulon.
  *What does it make?* → BGC.
- **One at-a-glance summary per gene, details on demand.** Routing signals in
  the existing style (`Gene.annotation_types` precedent): e.g. a
  `group_kinds: str[]` ("island, bgc, imodulon") + counts, so a gene page
  says "member of: 1 island (2 sources agree), 1 BGC (antiSMASH), 3
  iModulons" without dumping every instance. The default view is the summary;
  conflicts and per-source detail are one hop deeper, opt-in.
- **Don't rank sources globally; scope defaults per kind.** For islands,
  published+multi-source beats single-tool; for BGCs there is only antiSMASH;
  for modules there is only the paper. A per-kind default source (with
  everything else reachable) beats a universal precedence rule.

### Open consumer-side questions (for the design, later)

- Does the MCP layer get a `gene_groups(gene)` view and a
  `group_members(group)` view, and does group-level enrichment become a tool?
- Do groups fold into `annotation_types`/routing signals, and under what
  gating (the tier-gate/informative-gate precedents)?
- How do `Publication_discusses_*` edges interact with group types (papers
  discuss ISL4 by name)?

## 7. Relationship to current work (nothing blocked)

- **GEO branch ships as planned.** Hackl islands and iModulons land as
  `gene_clusters` (its D4/D5). Migration to a first-class group model later is
  cheap — source of truth is paperconfig + derived CSVs; the graph rebuilds.
- **antiSMASH Phase 1 can proceed independently** whenever desired — Phase-1
  artifacts (calls.json) are schema-free by design; only its Phase 2 touches
  this concept.
- **The GEO plan's Tier-3 specs** (TSS/operon/UTR/ncRNA entities; TF binding
  sites / regulator→target; iModulon activity) are all facets of this same
  concept and should be settled together with — or explicitly against — an
  eventual group design, not independently.
- **This doc's successor is a real design spec, expected soon** (approaches,
  schema, trade-offs) via the brainstorming → spec → plan path. This is
  positioned as the next major step in the KG's evolution, not a someday item.

## 8. Questions the design must answer (agenda for the spec conversation)

1. **Schema shape**: one generic instance label (`GeneGroup`/`GenomicRegion`
   with `group_kind`) vs per-kind labels (`GenomicIsland`, `BgcRegion`,
   `Replicon`, …) vs continued GeneCluster reuse — and where the
   region/system/set distinction lives (label vs property).
2. **The type layer**: node type(s) for cross-strain group identities (island
   families, BGC product classes, system types), and how instance→type edges
   are computed (published mappings like Doré S6 vs derived ortholog-anchored
   synteny, §2b option c) and provenanced.
3. **Membership semantics**: positional (region contains gene) vs role-typed
   (system component) vs weighted (module) — one edge type with kind-specific
   properties, or per-kind edge types.
4. **Multi-source handling**: corroborate-don't-merge mechanics (overlap
   edges? corroboration status on membership? per-kind default source) —
   §6's principles turned into schema.
5. **Consumer surface**: gene-page summary (routing signals, `group_kinds`),
   MCP views (`gene_groups`, `group_members`, group-level enrichment),
   full-text/index treatment, and gating into `annotation_types`.
6. **Migration path**: which existing encodings move (Hackl islands,
   iModulons as gene_clusters; D4/D5 anticipated this) and which stay
   (expression clusters are genuinely GeneClusters).
7. **Sequencing**: which layer lands first as the pilot (antiSMASH Phase 2?
   Hackl+Doré islands? replicon typing as the cheapest?) and what each
   subsequent layer reuses.
8. **Relationship to deferred specs**: the GEO plan's Tier-3 items
   (TSS/operon/UTR/ncRNA entities, TF-binding-sites/regulator→target,
   iModulon activity) — settled inside this design or explicitly out.
