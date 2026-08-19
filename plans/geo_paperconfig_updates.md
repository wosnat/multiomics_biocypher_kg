# Scope: paperconfig updates from the GEO processed-file drop

Status: **2026-08-19 execution pass COMPLETE (wiring side)** — everything the
pass scoped is wired, validated, and committed on `data/geo-processed-supplements`;
awaiting the Docker rebuild in the other clone. Shipped: 2.2 (Steglich, 94.5%),
2.4 (Hackl islands, 11 strains, 6,003 rows, 100%), 2.6 (Voigt, 99.8%),
2.3 (Doron, D7=(a), three hosts 97.9/97.1/93.9%), 2.5 sense arm (MIT0604,
97.9% via prefix id_translation), 2.1A (Johnson 2026b DE, 97.8%), the
`johnson 2026 a/` → `Johnson 2026b/` rename + triage-doc F8 correction, and
the **WH8109 onboarding** (registered + prepare_data + all six tools; loop-back
steps 1–7 pending tool completion). Measured deltas vs plan: Steglich cluster
membership is ~1.0K (not ~1.9K — the 12 clusters cover 1,102 rows); Voigt
resolution beat the estimate (old-locus-tag column resolves directly); WH8109
first measurement 2,589/2,757 = 93.9% (unresolved = JCVI features dropped by
RefSeq PGAP). New vocab registrations: 8 KNOWN_METRIC_TYPES, cluster_types
`decay_pattern`/`genomic_island`/`response_pattern`, `chemical` promoted into
the closed `Experiment.treatment_type` vocabulary (first emitter: Hackl MMC).
Test-mode build gate green (4/4) after two catches (chemical vocab; pandas
bool round-trip on Voigt's has_primary_tss). Discuss-topics extraction run
for all six new papers. Remaining for the pass: WH8109 prepare_data loop-back,
Docker rebuild + `pytest -m kg` + `/omics-edge-snapshot` + `/check-gene-ids`
(per-item acceptance numbers below), snapshot regeneration.

Previous status: Tier 1 done pending build (1.2 Huang wired; 1.1 he 2022 wired
and **unblocked** — Blocker B1 was fixed in `a84db12b` and the rebuilt MED4
mapping verified clean). Tier 2 planned 2026-08-18: six items, four promoted
out of Tier 3 on evidence — see Findings F1–F4.

**2026-08-19 execution pass (scoped with user):** today ships the easy items —
**2.2, 2.3 (D7 settled = (a)), 2.4, 2.5 sense arm only, 2.6, and 2.1A (DE
only)** — in the plan's suggested order (2.2 → 2.4 → 2.6 → 2.3 → 2.5 → 2.1A),
**plus new-strain onboarding where a wired paper needs it**: Synechococcus
**WH8109** via `/add-a-strain` (kicked off first so its background tool runs
overlap the paper work), which promotes the **Doron WH8109 arm** out of Tier 3 —
2.3 becomes three hosts. Everything else is **filed as backlog items** in
`plans/backlog.md` (section "GEO processed-supplements drop"): 2.1B diel
metrics, 2.1C iModulons (+ the PCC7942 model and `media-5` module homology),
the Hackl antisense arm (D8), the remaining Tier-3 items (munoz 2022 — not a
plain strain add, it needs a reference-proteome-match organism + probe mapping;
Hackl MIT1327 island remap; iModulon activity layer; TSS/operon/UTR/ncRNA
entities; TF-binding-site layer), the B1 regression check, and the
runaway-strain mapping re-check (W3-18-1, PCC7002, KT2440). The
`johnson 2026 a/` → `Johnson 2026b/` rename still happens today because 2.1A
needs the directory.
Branch: `data/geo-processed-supplements`
Input: commit `6acdd54f` (84 files, 7.4 MB) + the seven paper dirs from `dd594f23`
Triage doc: `docs/geo_prochlorococcus_candidates.md`

## What this covers

Turning the newly staged GEO/journal files into `paperconfig.yaml` entries. Ten
directories hold new material. They are **not** one kind of work — three tiers
plus one blocked group, gated by a single policy question (see Decisions).

---

## Tier 1 — trivial: published DE, native IDs, deployed strain

### 1.1 he 2022 — replace the DE source with the full tables

| | now | after |
|---|---|---|
| source | paper supp S1/S2 | `GSE195946_Group_*_DE_anno.txt.gz` |
| MED4 rows | 31 | 2,041 |
| NATL1A rows | 82 | 2,239 |
| `table_scope` | `significant_only` | `all_detected_genes` |
| `prefiltered` | `true` | `false` |

- IDs: `GeneID` = `gene-PMM0087` / `gene-NATL1_08781` → **tier-1 hit, verified**
  (`gene-PMM0087` → `PMM0087`; `gene-NATL1_00011` → `NATL1_00011`). No bridge work.
- Columns: `logFC`, `PValue`, **`FDR`** (the current entries use `p-value` as the
  adjusted column — the new files give a real FDR).
- **Replace, do not add**: two `csv` entries pointing at the same experiment would
  emit two `Changes_expression_of` edges per gene. The S1/S2 entries come out.
- Note: the paper's S1/S2 ("highly DE": p<0.05 AND abs(logFC)>1) are a *stricter*
  subset than GEO's own `_DE_significant_anno` (54 / 139 rows). Neither is the
  full table; both are superseded.
- Experiments/organisms unchanged. No new strain, no new schema.

**Acceptance:** `/omics-edge-snapshot` shows He 2022 rising 113 to ~4,280 edges with
no other paper changing; unresolved rate < 2%; `pytest -m kg` green;
`adjusted_p_value` non-null on all new edges.

### 1.2 Huang 2020 — new paperconfig, WH7803 phage infection

- `GSE150732_WH7803.Differential_Expression.xls`: 2,528 rows, 5 timepoints
  (15 min, 1, 3, 5, 7 h vs CT30m), paired `logFC-*` / `Pvalue-*` columns.
- IDs: `Gene Locus` is mostly native `SYNWH7803_RS#####` → **tier-1 verified**
  (`SYNWH7803_RS03435` → `SynWH7803_0680`); a minority carry a gene name
  instead (`rbcL`, `ppnK`), which also resolves (`rbcL` → `SynWH7803_0678`).
- Shape: 1 Publication, 1 Experiment (`treatment_type: [viral]`,
  `treatment_organism` = phage S-SBP1), 5 statistical_analyses (one per timepoint).
- `Pvalue-*` is a raw p-value, not adjusted → set `pvalue_col`, leave
  `adjusted_p_value_col` unset (edges get null `adjusted_p_value`, as with
  Thompson 2016).

**Acceptance:** new Publication + Experiment + ~12.6K edges; match rate > 95%;
registered in `data/Synechococcus/papers_and_supp/paperconfig_files.txt`;
`pytest -m kg` green.

---

## Tier 2 — medium: new paperconfig, several entry types

**Six items.** Four were in Tier 3 in the previous revision. The demotion of
those blockers is evidence-based — see Findings F1–F4 under Decisions. Every
item below is on a **deployed strain with native, already-resolving IDs**; none
needs `/add-a-strain` and none needs a schema change.

| # | item | strains | shape | entry types | est. edges |
|---|---|---|---|---|---|
| 2.1 | Johnson 2026b — iModulons + diel | MED4 (+PCC7942) | DE + periodicity + modules | `csv` x2, `derived_metrics_table`, `gene_clusters` | ~7.5K DE + 1.0K membership |
| 2.2 | Steglich 2010 — RNA half-lives | MED4 | per-gene scalars + decay clusters | `derived_metrics_table`, `gene_clusters` | ~5.8K metric + ~1.9K membership |
| 2.3 | Doron 2016 — cyanophage infection | WH7803, WH8102 | time-course microarray DE | `csv` x2 | ~35.1K |
| 2.4 | Hackl 2023 — genomic islands | 11 strains | genomic intervals → gene sets | `gene_clusters` x11 | ~6.0K membership |
| 2.5 | Hackl 2023 — MIT0604 MMC/UV DE | MIT0604 | unfiltered DESeq2 | `csv` + `id_translation` | ~6.0K (sense) |
| 2.6 | Voigt 2014 — TSS architecture | MED4, MIT9313 | per-gene TSS scalars | `derived_metrics_table` x2 | <=14.9K + 3.7K |

**Suggested order** — 2.2 → 2.4 → 2.6 → 2.3 → 2.5 → 2.1. Cheapest and most
self-contained first; 2.1 last because it is the only one with an unsettled
design question (D5) and the largest surface.

---

### 2.1 Johnson 2026b — GSE314951 iModulon paper

> **2026-08-19 split:** only **§A (DE, csv x2)** ships in today's pass, plus the
> directory rename it needs. §B (diel periodicity metrics) and §C (iModulons)
> are deferred to `plans/backlog.md` — same design as written here, nothing
> re-decided.

**Identity, first.** The directory `johnson 2026 a/` holds bioRxiv
`2026.04.15.718746` — *"Extreme genome reduction selectively retains modular
regulatory architecture in Prochlorococcus MED4"* (Johnson, Sadler, …, Bohutskyi;
PNNL). DOI `10.1101/2026.04.15.718746`. This is **not** the paper already sitting
in `Johnson 2026/` (`2026.01.21.700212`, the Lattice-Microbe whole-cell model,
which merely *cites* GSE314951). The triage doc attributes GSE314951 to the
modelling preprint; that attribution is wrong and should be corrected there too.

Rename `johnson 2026 a/` → `Johnson 2026b/` and strip the `media-N(3).xlsx`
browser suffixes in the same commit.

All IDs are `Refseq Locus` = `TX50_RS#####`, **verified 1,863 / 1,920 tier-1
(97.0%)** in MED4's `specific_lookup`, +14 via multi-singleton → 43 unresolved
(2.2%). MED4's canonical locus tag is `PMM*`, so `TX50_RS*` resolves as an
alternative — no bridge work.

#### A. Differential expression — `csv` x2

`media-2` sheets `Dark vs Light 2-4h (DESeq2)` and `Light vs Dark 14-16h (DESeq2)`,
1,920 rows each, each holding **two timepoints as paired columns** → 4 contrasts.

- Shape: 1 Publication, **2 Experiments** (one per perturbation direction), 4
  `statistical_analyses` (2 per csv entry, distinguished by `timepoint` /
  `timepoint_hours`).
- Experiment 1 — dark shift during the light period: `treatment_type: [darkness]`,
  `background_factors: [axenic, diel]`, timepoints 2 h and 4 h.
- Experiment 2 — light exposure during the dark period: `treatment_type: [light]`,
  `background_factors: [axenic, diel]`, timepoints 14 h and 16 h.
- Columns are `log2FC \n2H Dark vs Light` / `p-adjusted \n2H Dark vs Light` —
  note the **embedded newline** in every header. Export each sheet to a flat CSV
  with clean headers rather than quoting newlines in the paperconfig.
- `table_scope: all_detected_genes`, `prefiltered: false` — all 1,920 rows are
  present in every sheet.
- 21 °C, 12:12 L:D, Pro99, `omics_type: RNASEQ`, `test_type: DESeq2`.

#### B. Diel periodicity — `derived_metrics_table`

`media-3` `all_cosinor_fit` (1,872 rows) is the primary table: `acrophase`,
`amplitude`, `p_val`, `p_adj` per gene. Mirrors Waldbauer 2012 / Biller 2018.

Proposed metrics (all `numeric`; parent Experiment = a diel time-course
Experiment over the 36 circadian samples):

| metric_type | column | unit | rankable | registry |
|---|---|---|---|---|
| `diel_acrophase_rad` | `acrophase` | radians | false | **new** |
| `diel_amplitude` | `amplitude` | log2 | true | already registered (numeric) |
| `peak_time_h` | `cosinor_peak` (from `cosinor_padj_1e-2`) | hours | false | already registered (numeric) |

`p_adj` rides along as the metric's `adjusted_p_value` with
`has_p_value: true`, `p_value_threshold: 0.01` (the paper's own cutoff).

Do **not** also wire `fourier_padj_1e-2` / `intersect_` / `union_` as separate
metrics — they are filtered views of the same fit plus a re-analysis of Zinser
2009, and would triple-count the same gene. `fourier_padj_1e-2` carries
`peak_fourier`, which *is* a distinct estimator; if wanted it goes in as one
extra numeric metric (`peak_time_fourier_h`, since `peak_time_h` is taken by
cosinor), not as a second copy of the periodicity call.

#### C. iModulons — `gene_clusters` (the one new mechanic)

**The reshape problem in the old §2.1a does not exist.** `media-8` ships the 32
MED4 iModulons as **per-module member sheets already in long form** — one sheet
per module, one row per member gene, `gene_weight` in column B plus the full gene
annotation block. The ICA thresholding is already applied. Measured:
**32 modules, 1,011 membership rows, 737 distinct genes**, of which 726 are
tier-1 and exactly **1** is unresolved.

Build step: concatenate the 32 member sheets into one long CSV
`med4_imodulon_membership.csv` with columns `refseq_locus`, `imodulon`,
`gene_weight`. That is the entire derivation — a ~15-line script, committed
alongside the source.

```yaml
med4_imodulons:
  type: gene_clusters
  name: "MED4 iModulons (ICA, 248 samples)"
  filename: "data/.../Johnson 2026b/med4_imodulon_membership.csv"
  organism: "Prochlorococcus MED4"
  gene_id_col: "refseq_locus"
  cluster_col: "imodulon"
  cluster_type: "regulatory_module"
  cluster_method: "ICA (iModulonDB pipeline), 248 RNA-seq samples"
  score_col: "gene_weight"
  omics_type: RNASEQ
```

- **Multi-membership is safe.** `cluster_adapter` builds each edge id as
  `{cluster_id}__{gene_locus}`, and `cluster_id` already contains the module
  name, so a gene in three modules yields three distinct edges. Verified by
  reading the adapter, not assumed. The `Gene_in_gene_cluster` 1:1 worry in the
  old §2.1a is unfounded — but consumers that *assume* 1:1 should still be
  checked before the build (grep MCP/explorer for `Gene_in_gene_cluster`).
- **Per-module descriptions come from the paper, not an LLM.** The `iModulons`
  sheet carries `Class`, `Description`, `explained_variance`, `n_genes`,
  `TF_regs` (e.g. `+ lexA (PMM1262)`), `Prot_regs`, and KEGG/COG/GO enrichment
  strings per module. Hand-build `cluster_extractions/med4_imodulons.json` in the
  documented `{"metadata": …, "clusters": {key: {…}}}` shape from that sheet:
  `name` ← `name`, `functional_description` ← `Class` + `Description` +
  enrichment descriptions, `expression_dynamics` ← `explained_variance` +
  `TF_regs`. This is published metadata transcribed, so it is *more* trustworthy
  than the LLM path, and `load_extraction` does not care how the JSON was made.
- **`gene_weight` is signed.** ICA weights carry direction (a gene can be
  negatively weighted in a module). `score_col` lands the raw value on the edge;
  do not take an absolute value, and say so in the ClusteringAnalysis
  `experimental_context` so a consumer ranking by score knows to use `abs()`.

#### D. Deferred within 2.1

- **PCC7942 iModulons** (`media-6`: 78 modules, 2,706 genes, 520 samples) are a
  second, complete ICA model on a strain we deploy. Same mechanic, second
  `gene_clusters` entry, `organism: "Synechococcus PCC7942"`. Wire it **after**
  MED4 lands and only if MED4's membership edges look right — it doubles the
  new-mechanic surface for no new mechanic.
- **`media-5` BBH + iModulon homology** — MED4↔PCC7942 bidirectional best hits
  (1,018 pairs) and 19 module↔module correspondences with a Pearson r. The BBH
  layer is redundant with existing eggNOG/Cyanorak OrthologGroup bridging. The
  module↔module map has no slot (a GeneCluster→GeneCluster edge). Folded into D5.
- **`media-4`** — 75 RpaB binding-site predictions (TF, operon, coordinates,
  motif, q-value). A TF-binding-site layer. Out of scope; see Tier 3.
- **DIMA sheets** (26 rows x 2 contrasts, per-module activity difference +
  p-adj) and the **A matrix** — see **D5**.

**Acceptance (2.1):** A → 1 Publication, 2 Experiments, 4 analyses, ~7.5K
`Changes_expression_of` edges at >95% resolution. B → 1 Experiment, 2–3
DerivedMetric nodes, ~5.6K `Derived_metric_quantifies_gene` edges. C → 1
ClusteringAnalysis, 32 GeneCluster nodes, 1,011 `Gene_in_gene_cluster` edges,
every cluster carrying a non-null `functional_description`. `pytest -m kg` green;
`/omics-edge-snapshot` shows no other paper moving.

---

### 2.2 Steglich 2010 — MED4 RNA half-lives  *(was Tier 3)*

**F1 — the PowerPoint blocker is wrong.** `13059_2010_2347_MOESM5_ESM.XLS`
sheet `Tabelle1` is the whole-transcriptome half-life table: **2,043 rows**, all
with a non-null half-life. The legend's file numbering is offset by four (files
1–4 are the PPTs), so Additional file 1's table ships as MOESM5. Nothing needs
extracting from a PPT.

Columns: `old annotation` (`PMM0001`), `new annotation` (`PMED4_00001`),
`gene name and function`, `expression at 0 min in log2`, `expression` (bool),
`cluster`, `half-life time [min]` + SE, `decay rate [min]` + SE + error bounds.

**IDs.** `PMED4_*` resolves 1,928/2,042; `PMM*` resolves 1,691/2,043; taking
`new annotation` as `name_col` with `old annotation` as a fallback `id_column`
gives **1,930 / 2,043 = 94.5%**, 1,930 distinct locus tags (no collapsing). The
113 leftovers are `PMED4_pseudo_*` pseudogenes, `PMM_or0124`-style ncRNA/operon
features, and `-` placeholders — all legitimately absent from the gene layer.

> The `-` placeholder cells here are exactly the class of junk that caused
> Blocker B1. They are now rejected by the fix in `a84db12b`; do not add a
> Tier-1 `id_column` that would re-admit them.

Entries:

- `derived_metrics_table` (parent: a rifampicin RNA-decay Experiment):
  - `rna_half_life_min` — numeric, unit `min`, `rankable: true` — **new**
  - `rna_decay_rate_per_min` — numeric, unit `1/min`, `rankable: true` — **new**
  - `expression_at_t0_log2` — numeric, `rankable: true` — **new**
  - the two SE columns ride in `field_description`, not as separate metrics
    (they are precision, not measurement).
- `gene_clusters` for the 12 decay clusters (`cluster` column, values 1–12):
  `cluster_type: "decay_pattern"`, `cluster_method: "half-life clustering"`.
  Per-cluster descriptions need an extraction pass against the paper PDF — or
  ship without descriptions initially.

**Treatment-vocabulary question.** The measurement is a half-life, not a
treatment response; rifampicin is the *method*, not a perturbation whose effect
we report. `treatment_type` has no `chemical` value (`chemical` is a
`background_factors` value). Recommend `treatment_type: []` with
`background_factors: [axenic, chemical]` and the rifampicin detail in
`experimental_context`. Confirm against `validate_paperconfig.py` before writing.

**MIT9313 arm.** The paper covers MED4 *and* MIT9313, but only the MED4 table
ships as XLS. MOESM6 (operon half-lives, 220 rows) is also MED4. Any MIT9313
half-lives are inside the four PPTs — leave them out and note it in
`table_scope_detail`.

**Acceptance:** 1 Publication, 1–2 Experiments, 3 DerivedMetric nodes, ~5.8K
`Derived_metric_quantifies_gene` edges, 12 GeneCluster nodes, ~1.9K membership
edges; resolution ≥ 94%; `pytest -m kg` green.

---

### 2.3 Doron 2016 — cyanophage Syn9 infection, WH7803 + WH8102  *(was Tier 3)*

**F2 — there is no probe-mapping problem.** GSE63690 probe IDs are
`<locus_tag>_at`: `SynWH7803_1834_at`, `SYNW2135_at`, `Syncc8109_1439_at`. A
three-character suffix strip yields a native locus tag. Two of the three hosts
are **already deployed**, and the third is being onboarded in the 2026-08-19
pass (WH8109 via `/add-a-strain` — see the execution-pass note at the top):

| file | host | in KG? | host rows | resolution |
|---|---|---|---|---|
| `GSE63690_Processed_data_WH7803.txt.gz` | WH7803 | yes | 2,580 | 2,510 / 2,580 = **97.3%** |
| `GSE63690_Processed_data_WH8102.txt.gz` | WH8102 | yes | 2,585 | 2,504 / 2,585 = **96.9%** |
| `GSE63690_Processed_data_WH8109_Syn9.txt.gz` | WH8109 | **onboarding 2026-08-19** | 2,758 | measure after onboarding |

Each file also carries ~232 `BSV9_gp*` phage probes — **filter them out**; the
phage has no Gene nodes and they would become dangling edges.

**Shape.** Value columns are per-timepoint infected/control log-ratios:
WH7803 `S0.C0, S30.C30, S60.C60, S90.C90, S120.C120, S180.C180, S240.C240`
(minutes); WH8102 `SYN9.0.C0 … SYN9.4.C4` (hours); WH8109 per its own column
set (read at wiring time). So: 1 Experiment per host,
`treatment_type: [viral]`, `treatment_organism` = *Synechococcus* phage Syn9,
7 `statistical_analyses` each (one per timepoint), `timepoint_hours` from the
column label. The WH8109 entry lands **after** its `/add-a-strain` onboarding
completes far enough for gene-ID resolution (mapping built); if the onboarding
stalls on tool runs, ship WH7803 + WH8102 and add WH8109 as soon as its
mapping exists.

**The p-value caveat is the one thing to get right.** limma reports a single
omnibus `F` / `P.Value` / `adj.P.Val` per gene for the whole series — *"this
gene changed at some point"* — not a per-timepoint p. Wiring the same
`adj.P.Val` onto all 7 timepoints would make `expression_status` claim
timepoint-specific significance the paper never asserted. See **D7**.

Registration goes in `data/Synechococcus/papers_and_supp/paperconfig_files.txt`
(per D3).

Note the GSE74921 RNA-seq arm stays out: it is **phage-only** (every row
`NC_008296`), as `geo_source.txt` records. The journal supplements
(`moesm37/38/39`) are likewise phage gene tables (syn9 genes, 45 asRNAs, P-TIM40
expression clusters) — none of them host DE.

**Acceptance:** 1 Publication, 3 Experiments (2 if WH8109 onboarding hasn't
produced a mapping yet), 21 analyses (14 without WH8109), ~35.1K expression
edges from the two deployed hosts (2,510 x 7 + 2,504 x 7) plus up to ~19K from
WH8109 (2,758 x 7 at a comparable resolution rate), 0 dangling `BSV9_*` edges,
`Tests_coculture_with` → Syn9 on every experiment; resolution ≥ 96% per host.

---

### 2.4 Hackl 2023 — genomic islands as gene sets  *(new; D4)*

`mmc2.xlsx` = **5,598 predicted genomic islands across 623 genomes**, as
`(genome_id, contig_id, start, end)` intervals. Most of those genomes are SAGs
we will never hold — but **12 are our strains**, and the coordinate frames were
verified against the Hackl assembly metadata in `mmc1`:

| our strain | Hackl assembly size | our max gene end | Δ | islands | genes in islands | % of genome |
|---|---|---|---|---|---|---|
| MED4 | 1,657,990 | 1,657,817 | 173 | 10 | 525 | 26.6% |
| AS9601 | 1,669,886 | 1,669,718 | 168 | 12 | 492 | 25.2% |
| MIT9301 | 1,641,879 | 1,641,711 | 168 | 11 | 501 | 25.9% |
| MIT9312 | 1,709,204 | 1,709,069 | 135 | 12 | 516 | 26.1% |
| MIT9313 | 2,410,873 | 2,410,637 | 236 | 12 | 551 | 18.7% |
| MIT9303 | 2,682,675 | 2,682,470 | 205 | 19 | 666 | 21.4% |
| NATL1A | 1,864,731 | 1,864,543 | 188 | 14 | 750 | 33.7% |
| NATL2A | 1,842,899 | 1,842,712 | 187 | 14 | 700 | 31.6% |
| SS120 | 1,751,080 | 1,750,904 | 176 | 8 | 501 | 25.5% |
| MIT1314 | 1,704,447 | 1,704,332 | 115 | 14 | 402 | 21.3% |
| RSP50 (Hackl `RS50`) | 1,656,133 | 1,656,033 | 100 | 13 | — | — |
| **MIT1327** | 2,579,618 | **715,496** | **1,864,122** | 11 | — | **EXCLUDE** |

Every Δ is the trailing distance from the last gene to the contig end — i.e. the
same assembly. MIT1327 is a **different assembly** (23 Hackl contigs vs our 31,
and our annotation reaches only 715 kb on any single contig), so its islands
cannot be transferred without a coordinate remap. **Wire 11 strains, defer
MIT1327 to Tier 3.**

Independent biological validation on MED4: island genes are **50% hypothetical**
(54% in the largest island) against **24% genome-wide**. That is the textbook
flexible-genome signature and confirms the intervals landed on the right
sequence, not merely on a sequence of the right length.

**Modelling (D4, decided): `gene_clusters`, one analysis per strain.**

```yaml
med4_genomic_islands:
  type: gene_clusters
  name: "MED4 predicted genomic islands (Hackl 2023)"
  filename: "data/.../Hackl 2023/islands/MED4_island_membership.csv"
  organism: "Prochlorococcus MED4"
  gene_id_col: "locus_tag"
  cluster_col: "island_id"
  cluster_type: "genomic_island"
  cluster_method: "tycheposon-pipeline genomic island prediction"
```

- Membership CSV is **derived, not published**: a committed script joins `mmc2`
  intervals against each strain's `gene_annotations_merged.json` on
  `start >= island.start AND end <= island.end`. Containment (not overlap) is the
  right predicate — a gene straddling a boundary is not an island gene. Emit the
  script's per-strain counts into a report so the table above stays checkable.
- `island_id` = `{genome_id}_{contig}_{start}_{end}` so the cluster node id
  carries the coordinates, which GeneCluster has no dedicated fields for. Put the
  human-readable coordinates in the cluster `name` too.
- Descriptions: no per-island narrative exists and none should be invented. Ship
  without `cluster_extractions/`; `functional_description` stays null.
- 11 ClusteringAnalysis nodes, ~139 GeneCluster nodes, ~6,000 membership edges
  (5,604 measured across the 10 strains in the table above; RSP50 not yet
  counted — its annotation cache was not read in this pass).

**Not a core/flexible label** — see D6.

**Acceptance:** 11 ClusteringAnalysis nodes each with a
`Clusteringanalysis_belongs_to_organism` edge, ~139 GeneCluster nodes with
`member_count > 0`, ~6,000 `Gene_in_gene_cluster` edges, zero dangling; MED4's
island-gene hypothetical fraction reproducible from the graph (≈50% vs 24%
genome-wide) as a spot check.

---

### 2.5 Hackl 2023 — MIT0604 mitomycin C / UV DE  *(was Tier 3)*

**F3 — MIT0604 is already in the KG.** It is registered in
`cyanobacteria_genomes.csv` (`GCF_000757845.1`, taxid 1501268, Cyanorak
`Pro_MIT0604`) with **2,137 genes** merged and a built `gene_id_mapping.json`.
The Tier-3 note calling for `/add-a-strain` predates that onboarding.

**F4 — the ID bridge is one prefix.** The authors' `MIT0604_00001` is our
Cyanorak id `CK_Pro_MIT0604_00001` minus the `CK_Pro_` prefix. Prefixing and
looking up gives **2,003 / 2,047 = 97.9%** tier-1 hits. (A coordinate join via
`pro-deseq2.tsv`'s `start`/`end` is an independent cross-check and reaches
88.3% — lower, so use the prefix as primary and the coordinate join only to
audit the 44 misses.)

The bridge is a generated `id_translation` CSV: for every `CK_Pro_MIT0604_*` key
in `specific_lookup`, emit `(MIT0604_*, locus_tag)`. Deterministic, no diamond
run, ~10 lines — but it *is* a generated artifact, so commit the generator next
to the CSV.

**Tables.** Six DESeq2 tables, all unfiltered (~2,047 rows):
`{MMC_vs_control, UV_vs_control, MMC_vs_UV} x {sense, antisense}`.

- **Sense** (3 tables): 1 Publication, 3 Experiments (MMC and UV vs control;
  MMC-vs-UV is a treatment-vs-treatment contrast, so set `control_condition` to
  the UV arm and say so in `experimental_context`), `table_scope:
  all_detected_genes`, `prefiltered: false`, `test_type: DESeq2`, columns
  `log2FoldChange` / `padj`. Check the `treatment_type` vocab for both arms —
  neither mitomycin C nor UV has an obvious canonical value today.
- **Antisense** (3 tables): **open question D8.**
- The URL-encoded `%2C` in `geneFunction` (see Working notes) appears in
  `pro-deseq2.tsv`, not in the GEO tables we wire; still decode if any text field
  reaches a node property.

**Acceptance (sense arm):** 1 Publication, 3 Experiments, 3 analyses, ~6.0K
expression edges at ≥97% resolution; `id_translation` generator committed;
`pytest -m kg` green.

---

### 2.6 Voigt 2014 — per-gene TSS architecture  *(was Tier 3; D9)*

TSS / operons / UTRs / ncRNAs have no schema slot and that has not changed. But
Tables S4 (MED4, 4,190 rows) and S5 (MIT9313, 8,746 rows) key every TSS to a
native `oldLocusTag` (`PMM*` / `P9313_*`) and **reduce cleanly to per-gene
scalars**:

| strain | TSS rows | distinct genes | resolution |
|---|---|---|---|
| MED4 | 4,190 | 1,647 | 1,592 / 1,647 = 96.7% |
| MIT9313 | 8,746 | 2,144 | 2,124 / 2,144 = 99.1% |

TSS types present: `gTSS` (primary/genic), `aTSS` (antisense), `iTSS` (internal),
`nTSS` (orphan), `aArray` (array-only antisense).

Proposed `derived_metrics_table` per strain, from an aggregation script that
collapses the multi-row-per-gene tables:

| metric_type | value_kind | derivation |
|---|---|---|
| `has_primary_tss` | boolean | any `gTSS` row for the gene |
| `antisense_tss_count` | numeric | count of `aTSS` + `aArray` rows |
| `internal_tss_count` | numeric | count of `iTSS` rows |
| `minus10_element_score` | numeric | max `-10 element score` across the gene's gTSS rows |
| `tss_distance_to_cds` | numeric | `TSS distance` of the gene's primary gTSS |

All five are **new** `KNOWN_METRIC_TYPES` entries. Parent Experiment is a
transcriptome-architecture mapping run (`omics_type: RNASEQ`,
`test_type: "TSS mapping (454 + Solexa + tiling array)"`), one per strain.

**What this deliberately does not capture:** transcript boundaries, UTR lengths,
operon structure, and the ncRNA catalogue — the paper's actual contribution. The
reduction is a *routing* signal ("does this gene have mapped antisense
transcription?"), not a representation of the paper. Say so in every
`field_description` so nobody reads `antisense_tss_count` as a complete antisense
annotation. A real TSS layer stays in Tier 3.

**Acceptance:** 1 Publication, 2 Experiments, 10 DerivedMetric nodes, up to
~14.9K `Derived_metric_quantifies_gene` (4 numeric metrics x 3,716 resolved
genes; sparser in practice — `tss_distance_to_cds` and `minus10_element_score`
exist only for genes with a gTSS, `internal_tss_count` only for iTSS genes) +
~3.7K `Derived_metric_flags_gene` edges; aggregation script committed with its
row-count report.

---

## Tier 3 — genuinely blocked or needs its own plan

Four of the previous six Tier-3 items moved to Tier 2 above. The WH8109 arm
moved into the 2026-08-19 pass (onboarding via `/add-a-strain`, see 2.3).
What remains — **each filed as a backlog bullet in `plans/backlog.md`**:

| item | blocker | size |
|---|---|---|
| **munoz 2022 (GSE154594)** | Nothing staged but two PDFs. Environmental samples at station ALOHA on a custom Agilent array (`GPL28884`), 129 MB RAW.tar not downloaded. Needs probe→gene mapping *and* a reference-proteome-match-style organism (the `Marinobacter (MarRef v6)` precedent). A different kind of integration. | own plan |
| **Hackl 2023 — MIT1327 islands** | Hackl used a different MIT1327 assembly (2.58 Mb / 23 contigs) than our `GCF_001632125.1`. Needs a contig+coordinate remap before its 11 islands can transfer. Everything else in 2.4 is unaffected. | small, but its own verification |
| **iModulon activity layer** | See **D5**. The A matrix (32 modules x 248 samples) and DIMA tables (per-module activity difference + p-adj per contrast) have no representation. | own spec |
| **TSS / operon / UTR / ncRNA entities** | No node types. Would serve Voigt 2014 (TSS, operons, UTRs, ncRNAs), Steglich 2010 (operon half-lives, MOESM6/MOESM7) and Doron 2016 (45 asRNAs, moesm38) **together** — three papers already in hand, which is what makes it worth a spec rather than a one-off. | own spec |
| **TF binding sites / regulator→target** | Johnson `media-4`: 75 RpaB sites with operon, coordinates, motif, q-value. Plus `TF_regs` on every iModulon. A regulator→target layer the KG has nothing of today. Related to D5 — an iModulon's regulator and a TFBS are two views of the same missing edge. | own spec |

---

## Tier 4 — counts-only: NOT WIRED (decided)

Six datasets are **raw count matrices with no DE**: Thompson 2016 (GSE79359,
40 per-sample tables), Biller 2018 (GSE93197), biller 2016 (GSE73511),
coe 2024 (GSE264347), tetu 2019 (GSE118155), and GSE335434 (111 samples,
under `johnson 2026 a/`).

Every expression edge in the KG today comes from a **published** DE table. None
of these can become edges without us running DESeq2/edgeR ourselves — a new
capability, a new provenance story (`test_type` would name our own analysis),
and a reproducibility burden. Wang 2014's GSE49517 RPKM is the same question in
a different shape (per-sample RPKM, not a contrast).

**Decision D1: no.** We do not compute our own DE in this pass. The invariant
holds — every `Changes_expression_of` edge traces to a peer-reviewed analysis.
These files stay committed as reference/provenance: they document the detection
background of each experiment and let us check `table_scope` claims against the
actual gene universe. Revisit per-dataset if a specific question needs it; the
`geo_source.txt` in each directory records what is there.

**No paperconfig changes for:** Thompson 2016, Biller 2018, biller 2016,
coe 2024, tetu 2019, Wang 2014, and `johnson 2026 a/GSE335434/`.

---

## Findings that changed the plan (2026-08-18)

Each was verified against the files and the deployed caches, not inferred.

- **F1 — Steglich 2010 is not a PowerPoint problem.** The whole-transcriptome
  half-life table ships as `MOESM5_ESM.XLS` (2,043 rows, native `PMM*`/`PMED4_*`).
  The legend's file numbering is offset by four. → Tier 3 to 2.2.
- **F2 — Doron 2016 needs no probe mapping.** Probe IDs are `<locus_tag>_at`;
  WH7803 and WH8102 are deployed and resolve at 97.3% / 96.9%. Only the WH8109
  arm is blocked. → Tier 3 to 2.3 (WH8109 stays in Tier 3).
- **F3 — MIT0604 is already onboarded** (`GCF_000757845.1`, 2,137 genes,
  Cyanorak `Pro_MIT0604`). The `/add-a-strain` blocker is stale. → Tier 3 to 2.5.
- **F4 — the Hackl ID bridge is a prefix.** `MIT0604_00001` = our
  `CK_Pro_MIT0604_00001` minus `CK_Pro_`; 97.9% tier-1. No diamond, no GFF.
- **F5 — Hackl's genomic islands cover 11 deployed strains**, coordinate frames
  verified by assembly size (Δ ≤ 236 bp) and biologically validated on MED4 (50%
  hypothetical in islands vs 24% genome-wide). MIT1327 fails the frame check.
- **F6 — Johnson's iModulon M matrix never needs melting.** `media-8` already
  ships thresholded per-module member sheets (32 modules, 1,011 rows).
- **F7 — Blocker B1 is fixed and deployed** (`a84db12b`). MED4's max tier-1 id
  count is 16 against a median of 8 (was 273); `--` and `-` are no longer in
  `specific_lookup`; 42 conflicts remain (was 43). he 2022 is unblocked.
- **F8 — `johnson 2026 a` and `Johnson 2026` are two different papers.** The
  triage doc's attribution of GSE314951 to `2026.01.21.700212` is wrong.

---

## Decisions taken

- **D1 — compute DE from count matrices? NO.** See Tier 4. Six datasets stay
  reference-only.
- **D2 — Johnson scope: A + B + C.** All three pieces, including iModulons as
  `gene_clusters`. C is the only genuinely new mechanic in this pass; see 2.1C.
- **D3 — Synechococcus registration: yes, in the Synechococcus registry.**
  There are **two** paperconfig lists and `create_knowledge_graph.py` loads both —
  `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` and
  `data/Synechococcus/papers_and_supp/paperconfig_files.txt` (the latter already
  registers Tal 2009, Beliaev 2014, Kratzl 2024, Ma 2022, Oleza 2015, Oleza 2017,
  kaur 2018, Bernstein 2017). Huang 2020 and Doron 2016 go in the Synechococcus
  list. `multiomics_kg/utils/paperconfig_utils.py` holds the canonical pair in
  `PAPERCONFIG_LISTS`.
  Note: there is only **one** genome registry —
  `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` holds every strain,
  Synechococcus included.
  Side note: CLAUDE.md's "Reads three registry files" section names only the
  Prochlorococcus paperconfig list — worth correcting when we next touch that doc.
- **D4 — genomic islands are `gene_clusters`, one analysis per strain.**
  `cluster_type: genomic_island`, membership by coordinate containment,
  coordinates encoded in the cluster id and name. Zero schema change; island
  identity and coordinates are preserved, which a per-gene boolean would lose.
  Rejected: a first-class `GenomicIsland` node type — more faithful (an island is
  a locus, not a clustering result) but it is a schema change, a new adapter, new
  post-import rollups, and the `gene_clusters` encoding loses nothing we can
  currently query on.
- **D6 — core/flexible stays published-only.** `pangenome_membership` is added
  only where a paper publishes such a column (Wang 2014 today, as a categorical
  DerivedMetric). We do **not** derive a core-vs-flexible label from Hackl's
  islands, and we do not compute one from OrthologGroup coverage. Hackl's islands
  are *predicted variable regions within one genome*; pangenome core/flexible is a
  *cross-genome presence* statement. They correlate strongly and are not the same
  claim, and synthesizing one from the other would put an unpublished inference
  behind a name that already means something published. **No core-gene work in
  this pass.**
- **D9 — Voigt 2014 is wired as per-gene derived metrics.** Five metrics per
  strain (2.6), explicitly labelled as a routing signal rather than a
  representation of the paper's TSS/UTR/operon contribution, which stays Tier 3.

### D5 — iModulon activity: recommendation, not yet decided

The membership half (2.1C) is settled. The open half is **module activity**: the
A matrix (32 modules x 248 samples) and the DIMA tables (per-module activity
difference + p-adj + n significant genes, per contrast). Nothing in the schema
represents a measurement whose *subject is a gene set*.

**Recommendation: ship membership now as `GeneCluster`; give activity its own
spec later; do not create a `RegulatoryModule` node type in this pass.**

Reasoning:

- Membership semantics are **identical** to what `GeneCluster` already models —
  a gene set with per-member weights and analysis provenance. A parallel node
  type would duplicate ClusteringAnalysis/GeneCluster wholesale, including
  post-import rollups, indexes and full-text. `cluster_type: regulatory_module`
  already says "this is an ICA regulatory module".
- The KG's precedent is that **new evidence shapes** get new node types
  (MetaboliteAssay, DerivedMetric) while **new groupings of genes** reuse
  GeneCluster. Activity is a new evidence shape; membership is not.
- Migration cost if we later promote iModulons to their own type is low — the
  source of truth is the paperconfig plus a derived CSV, and the graph is rebuilt
  from source. What breaks is consumer knowledge ("iModulons are GeneClusters"),
  not data.

When the activity spec is written, the three candidate shapes are:

1. **Assay-node mirror of `MetaboliteAssay`** — a `ModuleActivityAssay` node per
   (Experiment x contrast), with `Assay_changes_activity_of → GeneCluster`
   carrying `activity_difference`, `adjusted_p_value`, `n_significant_genes`.
   Follows the closest existing precedent and keeps per-contrast provenance.
   *Preferred if activity is wanted.*
2. **`Experiment_changes_activity_of → GeneCluster`** — a direct mirror of
   `Changes_expression_of` with a cluster target. Simplest, but the KG has never
   had a measurement edge landing on anything but Gene or Metabolite.
3. **Properties on the GeneCluster node** — cheapest, but hard-codes one paper's
   contrasts into node properties. Rejected.

A first-class `RegulatoryModule` type becomes worth revisiting if and when the
regulator layer lands (TF_regs, `media-4` RpaB sites, module↔module homology in
`media-5`) — at that point a module is the hub of three edge types and no longer
looks like a clustering result. That is the Tier-3 "TF binding sites" spec, and
D5 should be settled together with it rather than separately.

### D7 — Doron omnibus FDR: **DECIDED (a)** (2026-08-19, with user)

limma gives one `adj.P.Val` per gene for the whole time series, not per
timepoint. Two readings:

- **(a) DECIDED** — set `adjusted_p_value_col` on **one** analysis per host
  (the terminal timepoint), `table_scope_detail` naming the omnibus semantics;
  the other six timepoints get a null `adjusted_p_value` (the Thompson 2016 /
  Huang 2020 precedent). Honest, but makes six timepoints look untested.
- **(b) rejected** — repeat the omnibus FDR on all 7 timepoints. Defensible as
  "this gene is significant in this series", but inflates
  `Experiment.significant_up_count` roughly sevenfold and makes per-timepoint
  `expression_status` claim something the paper never tested.

### D8 — Hackl antisense tables: open

The three antisense tables report differential expression of the *antisense
transcript at a sense gene's locus*. Emitting `Changes_expression_of` from an
antisense experiment to the sense Gene node conflates two different transcripts
at one address.

- **(a) recommended** — encode them as **separate Experiments** whose `name` and
  `table_scope_detail` state "antisense-strand transcription at this locus". The
  edge stays honest ("this experiment measured a change at this locus") and the
  data stays filterable.
- **(b)** drop the antisense arm. Loses half of what GEO adds over the GitHub
  `pro-deseq2.tsv`.
- **(c)** wait for the ncRNA/antisense entity layer (Tier 3). Correct but blocks
  the sense arm on a spec that does not exist.

The sense arm can ship regardless of how this settles.

---

## Explicitly NOT in scope

- No schema changes, no new node/edge types.
- One new strain only: **WH8109** (added to the 2026-08-19 pass by user
  request; MIT0604 turned out to be already onboarded, see F3). munoz 2022's
  organism need is NOT a plain strain add and stays out.
- No post-import Cypher changes.
- No re-running `prepare_data.sh` steps 0-2 (no new genomes).
- No DE computed by us (D1) — Tier 4 stays reference-only.
- No core/flexible derivation (D6).
- No iModulon activity, TF-binding-site, or TSS-entity layers (D5 / Tier 3).
- Steps 3/4 (`build_gene_id_mapping`, `resolve_paper_ids`) **will** need re-running
  for every paper touched.
- New `KNOWN_METRIC_TYPES` registrations **are** in scope — 2.2 adds 3, 2.6 adds
  5, 2.1B adds 1. That is a vocab edit, not a schema change.

## Working notes

- Files are committed **gzipped**. `pd.read_csv` reads `.gz` transparently, but
  `resolve_paper_ids.py` writes `<stem>_resolved.csv` — a `.txt.gz` source yields
  an awkward stem. **Decompress the specific files we wire up**; leave the rest
  gzipped.
- `GSE118155` count CSVs use **CR-only** line endings.
- `GSE195946` and the Hackl tables carry **URL-encoded** text (`%2C`); decode
  before it reaches a node property (and note `clean_text()` silently rewrites
  the pipe and apostrophe characters).
- Johnson `media-2` column headers contain **embedded newlines**; export sheets
  to flat CSVs rather than quoting them in YAML.
- Doron `moesm39` sheet names are Hebrew (`גיליון1`); read by index, not name.
- Several Excel reads need `PYTHONIOENCODING=utf-8` on this Windows shell —
  `media-3`'s `imodulon_fit` sheet and Doron's `.xls` files contain non-cp1252
  characters that crash a bare `print`.
- Three derived CSVs will be generated and committed in this pass, each with its
  generator script: Johnson iModulon membership (2.1C), Hackl island membership
  (2.4), Hackl MIT0604 id_translation (2.5). Plus two aggregation outputs for
  Voigt (2.6).

---

## Blocker B1 — runaway ID merge in the MED4 gene mapping — **FIXED**

Fixed in `a84db12b` ("fix(gene-id): reject placeholder and free-text cells as
gene-unique IDs") and verified against the rebuilt mapping on 2026-08-18:

| | before | after |
|---|---|---|
| MED4 max tier-1 ids on one gene | 273 (`PMM0236`) | 16 (`PMM0767`) |
| MED4 median tier-1 ids | 8 | 8 |
| `--` / `-` in `specific_lookup` | present | absent |
| MED4 conflicts | 43 | 42 |

The two root causes were (1) placeholder cells (`--`, `-`, `NA`) typed as Tier-1
identifiers — Wang 2014's MED4 table used `--` for 99 RNA features, merging them
all into `PMM0236` — and (2) free text split into per-word "identifiers"
(Beliaev 2014 footnote rows). Both classes are now rejected at ingest.

**Consequences for this plan:**

- **he 2022 (1.1) is unblocked** and can be built.
- The other three runaway strains reported in the original scan (W3-18-1,
  PCC7002, KT2440) should be re-checked against the same rule before their papers
  are next touched — the fix is global but their mappings need a step-3 rebuild
  to pick it up.
- Steglich 2010 (2.2) and Wang 2014 both feed `-`-style placeholder cells. Do not
  add a Tier-1 `id_column` that would re-admit them.
- A regression check that no gene accumulates more than ~3x the median tier-1 id
  count is still worth adding; it was proposed as part of B1 and is not in the
  commit.
