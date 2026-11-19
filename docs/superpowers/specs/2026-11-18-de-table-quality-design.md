# DE table quality: fix KG bugs, flag paper limitations

**Date:** 2026-11-18
**Status:** **Decided** (2026-11-19). Frozen: changes after this date are logged as dated amendments.
**Inputs:**
- Explorer ask `multiomics_explorer/docs/superpowers/specs/2026-11-18-kg-side-de-value-provenance-ask.md` (F1–F7, plus our response in its §5)
- spock analysis `2026-11-16-hli_light_stress` `gaps_and_friction.md` (assessment items 3 and 4, the "Table sanity and coverage check" library candidate) and its decisions A.8, A.13
- Explorer backlog 2d.5 (a / b), 2d.6, E1, E6, E7
- KG `plans/backlog.md` (zero-count magnitudes, Barreto, row collisions) and `plans/signed_linear_fold_change.md` (2026-10-31 fold-change audit)

## 1. Goal

The end goal is a **table sanity and coverage check** in the MCP. An analysis should learn, before it reads any target values, which DE experiments it can trust and for what. The KG precomputes everything this needs; the explorer only reads it and presents it.

Problems fall into two kinds, and they are handled differently:

| | KG bug | Paper limitation |
|---|---|---|
| What it is | The KG stores something the source does not say: an import, resolution, curation or modelling error | The source itself is limited or ambiguous; the KG stores it faithfully |
| What we do | **Fix** it, then log the fix in CHANGELOG `### Data` | **Flag** it on the Experiment, so the MCP can warn |
| After the work | Gone from the graph | Visible to every consumer, with its reason |

**Decided:** stored values are not changed by the flagging work. Fixes for KG bugs are separate, one per bug.

**Decided (2026-11-18):** the line is **B**. A *bug* is only something the KG got wrong relative to what the paper published. Anything the paper published as-is is *flagged*, never recomputed, even when the KG could correct it from the paper's own columns. Consequences:
- Aharonovich 2016's three *Prochlorococcus* tables (§3.3) are a limitation: flag them.
- Shared-probe twins (§3.4) are a platform limitation. **Removed from this spec 2026-11-19** (§10); backlog item "microarray probe-level features".
- **Decided (2026-11-18):** the recomputations already in place (Aharonovich HOT1A3, Beliaev 2014, Fuszard 2012 NATL2A) stay for now, and those experiments get a flag saying the KG recomputed the fold change from the paper's own columns.

## 2. Non-goals

- Recording how every value was made (the ask's F2–F6 provenance fields). A possible later spec. F1 `fold_change_type` is the one exception and is in this spec (§7.5).
- Deciding for the user whether to drop an experiment. The KG states the problem; the analysis decides.
- Explorer-side work (tools, warnings, docs). This spec only defines what the KG publishes.

## 3. Known cases, provisionally sorted

Sources: the backlogs and audit listed above, plus checks made on 2026-11-18 against the source files and the live graph.

### 3.1 KG bugs (fix)

| Case | Evidence | Status |
|---|---|---|
| Two table rows resolve to one gene; BioCypher keeps one of the two edges by row order | 85 edges over 10 papers (Coe 2024: 49, of which 39 of the dropped copies were significant). Detection added 2026-11-11 | Open: decide per paper |
| Aharonovich 2016 experiment metadata: config says 24 °C, 10 µE; the Methods give 20 °C, 27 µE for maintenance | KG backlog P2 (b); Supplementary information missing | Open: needs the SI |
| Steglich 2006 publication title contains "/H17188" | Explorer E1 | Open. A Publication bug, not an experiment one |

### 3.2 Paper limitations (flag)

| Case | Evidence |
|---|---|
| **Zero-count inflated fold changes**, no raw counts to recompute from | Capovilla 2023 MIT9303 (±21.5; only transformed counts published), Anjur-Dietrich 2026 (+24.45 / +26.13), Huang 2021 |
| **Unreplicated libraries** | Huang 2021 |
| **Arms not normalised to each other** | Aharonovich 2016 MED4: expression columns sum to 156K vs 384K, median log2FC −1.46 (−0.93 with ≥ 10 reads in both arms); 82% of values negative |
| **Sign convention not verifiable from the deposited table** | Aharonovich 2016 MED4: the published value is log2(axenic column / coculture column), while the two MIT9313 tables are coculture / axenic. The paper's text agrees with the published sign |
| **Table may not be the contrast it is labelled as** | Barreto 2022 S12–S14: probably the partner × pCO2 interaction; log2FC anticorrelates with the 400 ppm tables (r −0.61 / −0.83 / −0.97) |
| **KG significance call does not reproduce the paper's** | Pandhal 2007 (KG 6, paper ~15), Christie-Oleza 2017 S7 (32 more than the paper's own Perseus call), Huang 2021 (paper's rule gives 47 up, not 32; the KG also applies a down rule the paper doesn't define) |
| **Filtered table** | Martiny 2006 phosphate: the source itself lists only 34 genes (q < 0.05 at 48 h). Already expressed by `table_scope = filtered_subset` + detail |

### 3.3 Settled by line B, and one case needing nothing

| Case | Outcome |
|---|---|
| Aharonovich 2016, the three *Prochlorococcus* tables: pseudocount 0.01, so 1 read vs 0 reads = log2 ≈ 6.7 (562 of 1,714 MED4 genes have a zero arm; 291 have \|log2FC\| > 5 but are not significant) | **Limitation (B):** flagged, not recomputed |
| Zinser 2009 diel: Experiment with 0 expression rows | By design: only clusters are loaded (RMA values, no DE contrast). `reports_fold_change = no_fold_change` already says so. Nothing to do |

### 3.4 What the sanity statistics find on their own (trial run, 2026-11-18)

A read-only Cypher query computed the statistics of spock's `11_check_coculture2016.py` (negative fraction, |log2FC| > 5 called not significant, log2FC values shared by ≥ 10 genes, max |log2FC|) for every experiment × time point with `table_scope = all_detected_genes` and ≥ 50 rows. 32 experiments crossed loose cutoffs (negative fraction > 0.75 or < 0.25, ≥ 20 large-but-ns, ≥ 5 tied groups, max > 12). What this showed:

- **Known problems that are caught:** all five Aharonovich 2016 tables. The MED4 one stands out on every statistic (0.82 / 291 / 18 / 15.4). The two MIT9313-vs-axenic tables read 85 and 178 large-but-ns.
- **Real biology that is caught too:** Doron 2016 phage infection (WH7803, WH8102, WH8109) is 97–100% negative after 30 min, because the phage shuts down host transcription. Tolonen 2006 N deprivation reaches 0.76. **A statistic is not a verdict.**
- **Known problems that are missed:** Capovilla 2023, Anjur-Dietrich 2026 and Huang 2021 (zero-count magnitudes) don't cross these cutoffs. **Statistics can't replace curated flags.**
- **New candidates, not yet inspected:**
  - `10.1128/spectrum.03275-22` proteomics (14 experiments): negative fractions 0.09–0.25 with 2–10 large-but-ns each, max |log2FC| up to 10. Possibly imputed values.
  - `10.1038/nmicrobiol.2017.100` late coculture WH7803 proteomics: 87% negative, 40 tied groups.
  - Al-Hosani 2015 AS9601 salt: 70 tied groups. Possibly just the source's 2-decimal rounding.

Other spock scripts:
- `twin_rule.py` (identical log2FC across paralog pairs in microarray tables = one shared probe): this can be precomputed in the KG per microarray experiment. Under B, this is a platform limitation: flag it, don't change the edges.
- `02_sign_and_spot_check.py` and `01_build_lists.py` (`gene_calls_from_background`) depend on a gene set. They are explorer-side (coverage and the call matrix), not KG precomputation.

## 4. What the KG precomputes

1. **Curated flags** for §3.2. **Decided (2026-11-18):**
   - `Experiment.data_caveats` (`str[]`): dense on every Experiment with a DE contrast; `[]` = no known problem, re-materialised post-import like `background_factors`. A closed vocabulary in `config/controlled_vocabularies.yaml`, each value with a researcher-facing `value_descriptions` line that the MCP shows verbatim.
   - `Experiment.data_caveat_note` (str, sparse): one or two lines of evidence. It names the time points or tables affected when a caveat covers only part of the experiment.
   - Curated in the paperconfig under the experiment (`data_caveats:` / `data_caveat_note:`). The validator checks values against the vocabulary and requires a note whenever the list is non-empty.
   - `Experiment.data_caveats_reviewed_on` (str, sparse, ISO date): when the experiment was last reviewed against this vocabulary (explorer review §6.2). Absent = not reviewed; `[]` with a date = reviewed and clean. Curated in the paperconfig next to `data_caveats`.
   - One flag set per experiment, not per time point: no known case needs per-time-point flags.
   - No severity levels: the MCP shows each caveat with its description and the analysis decides. A severity would be a verdict.
   - Labels are plain statements of what is true about the table, without abbreviations, readable without knowing the paper. A new problem that fits no value gets a new value, minted the way `treatment_type` values are.

   Starting vocabulary:

   | Value | Description (shown by the MCP) | First experiments |
   |---|---|---|
   | `fold_change_inflated_at_low_counts` | Genes with zero or very few reads in one condition get extreme fold changes. Do not rank or filter these genes by fold-change size. | Capovilla 2023 MIT9303, Anjur-Dietrich 2026, Huang 2021, Aharonovich 2016 *Prochlorococcus* ×3 |
   | `no_biological_replicates` | Each condition was sequenced once, so the p-values and significance calls rest on a model, not on measured variation. | Huang 2021 |
   | `conditions_not_normalized_to_each_other` | The two conditions differ in overall signal, which shifts most genes in one direction. Compare genes within this experiment, not across experiments. | Aharonovich 2016 MED4 |
   | `comparison_direction_unclear` | The source table does not make clear which condition is the numerator. The KG keeps the sign the paper states. | Aharonovich 2016 MED4 |
   | `table_may_be_a_different_comparison` | The published table probably reports a different comparison than its label says (see note). | Barreto 2022 S12–S14 |
   | `significance_calls_differ_from_paper` | The KG's significant / not-significant call uses a different rule from the paper's, so its counts differ from the paper's (see note). | Pandhal 2007, Christie-Oleza 2017 S7, Huang 2021 |
   | `fold_change_recomputed_by_kg` | The KG recomputed the fold change from the paper's own expression columns, because the published values were distorted (see note), so values will not match the paper's table. The original values are only in the paper. | Aharonovich 2016 HOT1A3 ×2, Beliaev 2014, Fuszard 2012 NATL2A |
   | `fold_change_sign_reversed_by_kg` | The paper's table reports the comparison the other way round (control vs treatment). The KG negated it so that positive means up in the treatment; signs are opposite to the paper's table. | Ferrer 2020, Barber 2026 ×3, Koch 2019 (the only builders that negate, checked 2026-11-19) |

   Every first assignment is still checked against the source during curation. The suspects from §3.4 and item 2 (Moreno 2023 proteomics, `nmicrobiol.2017.100` late coculture, Al-Hosani 2015, Kaur 2018 DSS-3 exoproteome, Biller 2018 NATL2A extended darkness, Koch 2019 83-1 glucose) are tagged only if that check confirms a problem.
2. **Table-sanity statistics.** **Decided (2026-11-18):**
   - Gene-agnostic, so they leak nothing. Computed post-import from the `Changes_expression_of` edges, per time point, as Experiment arrays parallel to `time_point_hours`:
     - `time_point_negative_fractions` (float[]): share of edges with log2FC < 0
     - `time_point_median_log2fc` (float[])
     - `time_point_max_abs_log2fc` (float[])
     - `time_point_large_fc_not_significant` (int[]): edges with \|log2FC\| > 5 and `expression_status = 'not_significant'`
   - **Dropped (2026-11-19): `time_point_tied_value_groups`** (log2FC values shared by ≥ 10 genes). Measured on the live graph: with fold change alone, the largest counts come from the source's rounding (Al-Hosani 2015: 70 groups / 878 genes, 2-decimal table; `nmicrobiol.2017.100` late: 40, discrete spectral-count ratios). Requiring an equal p-value as well removes those but leaves mostly trivial groups (log2FC 0 with p = 1: darkness MIT1002, DHPS proteomics) and rounded-p coincidences. Only Aharonovich MED4 shows the intended low-count defect under either rule, and `time_point_large_fc_not_significant` (291 genes), `time_point_extreme_fc_count` and the low-signal count already catch it. Duplicate rows are caught by `shares_probe_with` and the adapter's row-collision log.
   - **Changed (2026-11-19): on every DE experiment, whatever its `table_scope`** (was `all_detected_genes` only; this overrides the explorer's §6.1 endorsement of the split). The KG publishes the numbers; the MCP interprets them with `table_scope`. In a significant-only, filtered or top-n table they describe the published subset, not the measured genome (e.g. `time_point_large_fc_not_significant` is 0 by construction in a significant-only table), and each property description says so. Still gene-agnostic, so still leak-free.
   - **Decided (2026-11-18):** `time_point_extreme_fc_count` (int[]): edges with \|log2FC\| > 10 (> 1,000-fold), on **every** table scope. It is gene-agnostic and needs no comparison with the rest of the table, so it is valid on significant-only tables too. Share = count / `time_point_totals`. Why: it is the only low-count signal for tables that carry no level columns (Huang, Capovilla), and \|log2FC\| > 5 not significant misses DESeq2 tables, where a zero condition is *called significant* (Capovilla, Anjur). A whole-graph run found 16 experiments with any: every known case (Aharonovich ×4, Capovilla 2023 MIT9303, Anjur-Dietrich 2026) plus new suspects (Kaur 2018 exoproteome `rpom_dss3_exoproteome_asw` 15, Biller 2018 NATL2A extended darkness 3, Koch 2019 83-1 glucose 2, single values elsewhere). For reference, the largest known real response, Martiny phoA, is 162-fold (log2 7.3).
   - **Decided (2026-11-18):** `time_point_low_signal_count` (int[]): large fold changes that rest on little or no signal, measured from the source table's level columns (§4.4). Computed **by the adapter at build time** (it alone sees the table), only where the paperconfig declares level columns, on every table scope. Two rules, depending on what the table carries:
     - *Level columns per condition* (one average or several replicates per arm, see §4.4) (Aharonovich 2016, Thompson 2016, Fuszard 2012, Al-Hosani 2015, Moreno 2023, Koch 2019, Beliaev 2014): genes where either condition equals that column's **minimum**, i.e. no signal. The column's own minimum makes the rule unit-free: Aharonovich's 0.01 pseudocount, Beliaev's RPKM floor and VST's 9.85 all count as "no signal". No \|log2FC\| cutoff: a zero condition is the defect whatever the fold change.
     - *One mean column* (DESeq2 `baseMean`: Hackl 2023; edgeR `logCPM`: Anjur-Dietrich 2026): genes with \|log2FC\| > 5 whose mean is in the table's bottom 10%. Check: Anjur's two extreme genes (+26.1, +24.5) have logCPM 2.74 and 1.06; the table's 10th percentile is 4.47, median 6.94.
     - Plus `low_signal_basis` (str, vocabulary): `condition_at_column_minimum` | `mean_in_bottom_decile`, so the reader knows which rule produced the number. **Empty level cells (builder's call, 2026-11-19):** an empty cell in a declared level column counts as no signal, i.e. at the column's minimum. Proteomics tables (Koch 2019 NSAF; to check: Moreno, Fuszard) leave undetected proteins blank instead of writing 0, so ignoring blanks would flag the wrong rows. One value per experiment: the build **fails** if an experiment's analyses declare level columns of different kinds (§6.3.2; no such case is known).
     - Named `low_signal_count`, not `…large_fc…`: rule 1 has no fold-change cutoff on purpose (§6.3.1). Rule 2 needs one (\|log2FC\| > 5), because "bottom 10% of the mean" alone is 10% of every table by construction; the `mean_in_bottom_decile` value description says so.
     - The levels themselves are **not** stored (§4.4).
   - How the three low-count signals divide the work:

     | Table has | Statistic | Catches |
     |---|---|---|
     | Per-condition levels | `time_point_low_signal_count` (condition at column minimum) | Aharonovich-type zero counts |
     | One mean column | `time_point_low_signal_count` (\|log2FC\| > 5, mean in bottom 10%) | Anjur, Hackl-type |
     | No levels | `time_point_extreme_fc_count` (\|log2FC\| > 10, every table) | Capovilla, Huang, unknown cases |
   - Raw numbers, no verdict, no build warnings for now. The §3.4 trial showed any cutoff catches real biology (Doron 2016 phage shutdown) and misses known problems (Capovilla, Anjur, Huang).
   - Each post-import statistic is declared in `schema_config.yaml` with a `# post-import:` comment (house rule); the adapter-computed `time_point_low_signal_count` / `low_signal_basis` are declared as ordinary adapter properties. All covered by kg-validity tests.
2b. **Edge-level flags. Decided (2026-11-18; requested by the explorer review §6.1):**
   - `Changes_expression_of.low_signal_reason` (str, sparse, vocabulary): `condition_at_column_minimum` | `mean_in_bottom_decile`, set by the adapter on each row that meets the rule above. Absent = the row passed, or the table has no level columns; the Experiment's `low_signal_basis` tells the two apart (basis present ⇒ absence means passed). A named categorical, not a `'true'`-only presence marker (R5 allows exactly two of those; this would be a third).
   - `time_point_low_signal_count` then becomes a **post-import rollup** over flagged edges (`significant_up_count` pattern), so the rule lives once, in the adapter, and the count can't drift from the edges.
   - No edge flag for \|log2FC\| > 10: readable from `log2_fold_change` itself (R3, no threshold materialised over a value the edge carries). The Experiment count stays.
   - ~~`Changes_expression_of.shares_probe_with`~~ **Removed 2026-11-19** (§10).
   - Leak-safe: edge flags travel with the gene's values, so a consumer sees them only when it reads values it is already allowed to see.
3. **Coverage support. Decided (2026-11-18, after explorer review §6.4):**
   - The coverage check must be leak-safe: for significant-only tables it may say "on the platform, reported only if significant", never "present in the table".
   - ~~`Experiment.array_design`~~ **Removed 2026-11-19** (§10). Its only remaining use was pooling tables of one design (twins, per-design gene universe), both deferred. A missing gene in a microarray table is read with `omics_type = 'MICROARRAY'`: "not significant or not on the array". The GEO findings are kept in the backlog item.

### 4.4 Per-arm expression and abundance values (raised 2026-11-18)

Many source tables carry the levels behind each fold change, and the KG drops them. A header survey (regex, approximate) of the 170 registered DE `csv` tables found such columns in **~43 tables in 11 papers**:
- Rockhopper expression per arm: Aharonovich 2016
- DESeq2 `baseMean`: Hackl 2023 (one value per gene); Al-Hosani 2015 (also per arm: `baseMeanA` / `baseMeanB`)
- RPKM per arm: Thompson 2016 (18 tables)
- Proteomics means per condition: Fuszard 2012
- Spectral-count averages per condition, plus peptide count: Moreno 2023 (10 tables; this is `10.1128/spectrum.03275-22`, the new suspects in §3.4)
- Highest / lowest mean condition + peptide count: Domínguez 2017
- NSAF per replicate: Koch 2019

Units differ in every paper, so the levels are never comparable across experiments. Their value is **within an experiment**: they show, per gene, when a fold change rests on near-zero signal or a single peptide. Under B, that is how a consumer discounts an unreliable fold change without the KG recomputing it.

Options:
1. **Use them only at curation time.** We read them to diagnose a table and write its flag. Nothing is stored.
2. **Store them on the edge.** For example `treatment_level` / `control_level` (or `base_mean`) and `peptide_count`, plus an Experiment `level_kind` vocabulary (`rockhopper_expression`, `rpkm`, `deseq2_base_mean`, `spectral_count_mean`, `nsaf`, …). Post-import could then derive a per-edge low-signal marker and a per-experiment count. Covers ~25% of tables.
3. Both: 1 now, 2 as a later phase of this spec.

**Decided (2026-11-18):** levels are **not stored** on edges (option 2 deferred). They are used twice:
- at curation time, to diagnose a table and write its flag (option 1);
- **at build time**, by the adapter, to compute one count per time point, `time_point_low_signal_count` (§4 item 2). This revises the first answer to Q7 ("curation only").

For this, the paperconfig analysis gains optional level-column fields (names for the builder to settle). **Shape (decided 2026-11-19):** each arm is a **list** of columns, because tables publish levels either as one average per arm or as one column per replicate:
- `treatment_level_cols` / `control_level_cols` (lists): one column = the arm's average; several = its replicates. Rule 1 per arm: the arm has **no signal** when every one of its columns is empty or at that column's own minimum. With replicates this means "detected in no replicate"; a protein seen in one replicate of two is not flagged (that is a weaker, different signal, out of scope).
- `mean_level_cols` (list, usually one column, e.g. `baseMean`): rule 2. Several columns are averaged per row first.
- **Peptide-count columns are not declared either (checked 2026-11-19).** Moreno 2023 (15 tables), Domínguez 2017 and Fuszard 2012 publish them. Every one of those tables already applies the standard ≥ 2-peptide identification cutoff: 0 proteins with a peptide count of 1 (Domínguez: 2 with one *unique* peptide). A single-peptide rule would flag nothing, and flagging 2-peptide proteins would second-guess the paper's own accepted cutoff with a threshold of ours. In proteomics, low signal comes from the per-condition levels (spectral counts, NSAF, means), which rule 1 already uses.
- SD / SE / CV columns are **not declared**: no rule uses them. Low-signal here means "no or almost no signal", which the levels alone answer; spread would answer "noisy", a different question. Kept out unless a later rule needs it.
- Per-arm level columns that belong to several contrasts (Koch 2019 `nsaf_g1 … nsaf_m2`, one column per condition, shared by four comparisons) are declared per analysis with the columns of that analysis's two arms. The validator checks the columns exist. Declaring them for the ~11 papers above is part of the curation work.

## 5. Open questions

1. ~~Where is the line between fix and flag?~~ **Decided:** B (§1).
2. ~~Flag shape?~~ **Decided:** §4 item 1.
3. ~~Which sanity statistics?~~ **Decided:** §4 item 2.
4. ~~Coverage?~~ **Decided:** leak-safe rule only (§4 item 3); `array_design`, twins and per-design gene universes removed to the backlog (§10).
5. ~~Process?~~ **Decided:** §7.4 (explorer proposal §6.5, accepted).
6. ~~Existing recomputations under B?~~ **Decided:** keep them, with a flag (§1).
7. ~~Per-arm levels?~~ **Decided:** not stored; used at curation time and, at build time, for one low-signal count (§4.4, revised 2026-11-18).
8. ~~F1 `fold_change_type`?~~ **Decided:** in this spec (§7.5).
9. ~~Which experiments must state caveats + review date?~~ **Decided:** allowlist (§7.4).
10. ~~First-delivery scope?~~ **Decided:** mechanism + listed cases (§7.6).
11. ~~Paperconfig fields, validator report, `/paperconfig` QA step?~~ **Proposed 2026-11-19:** §7.7.

## 6. Explorer review (2026-11-18)

**From:** `multiomics_explorer`. We reviewed the spec as it stands, with every decision through §4.4. Overall it gives the explorer what backlog 2d.5 (the DE-table QC check) and 2d.6 (the call matrix) need. Below: one gap in the data model, several naming and semantics points, answers to Q4 and Q5, and the explorer work it implies.

### 6.1 Agreed, no change requested

- **Line B (§1).** We accept the cost: published inflated values stay in the graph, and `differential_expression_by_gene` sorts by \|log2FC\|. Aharonovich MED4's 1-vs-0 genes (up to ±15) will therefore sit at the top of MED4 results, and they skew `rank_by_effect` and `gene_response_profile`. The explorer's warnings, plus the edge flag in 2b, are what make B safe, so the explorer surface ships with the first KG build that carries `data_caveats`, not later.
- **`data_caveats` shape (§4.1):** a closed list per experiment with a sparse note, no severity levels, and plain labels. The explorer already passes `ControlledVocabulary.value_descriptions` through to users, so the descriptions will reach them verbatim.
- **The statistics (§4.2):** raw numbers, no verdict, no build warnings. The split between scope-dependent statistics (`all_detected_genes` only) and scope-independent ones (`time_point_extreme_fc_count`, the low-signal count) is right.
- **2b: please approve it.** `low_signal_reason` on the edge is what lets B work row by row. Without it the explorer could only say "this experiment has N low-signal genes", never which genes. Making the experiment count a rollup of the edge flags, so the rule lives in one place, is also right.

### 6.2 Gap: `[]` means both "checked, clean" and "never checked"

§4.1 says `[]` = "no known problem". But most experiments will never be reviewed against this vocabulary, so `[]` covers two different states:
- someone checked the table and found nothing;
- nobody has looked yet.

The explorer can't tell a user "this table was reviewed" apart from "this table wasn't flagged". This is the same ambiguity we raised for F6. Proposal:
- Add `Experiment.data_caveats_reviewed` (str, sparse): the date or curation-pass id of the last review against the caveat vocabulary.
- Absent means not reviewed. `[]` together with a review date means reviewed and clean.
- It's a string rather than a bool, so it satisfies R5.
- The validator requires it on new paperconfigs (see Q5).

### 6.3 Naming and semantics

1. **`time_point_low_signal_large_fc_count` doesn't match rule 1.** The per-condition rule deliberately has no fold-change cutoff ("a zero condition is the defect whatever the fold change"). The name promises "large fc", though. After the HOT1A3 recompute, a 1-vs-0 gene has log2FC ≈ 1 and still counts. Rename it `time_point_low_signal_count`, or add the same \|log2FC\| > 5 cutoff to rule 1. We prefer the rename, since the rule as written is the useful one.
2. **`low_signal_basis` is one value per experiment, but the rule is chosen per table.** If an experiment's time points come from tables with different level columns, the single value is wrong. Either fail the build on a mix (as F1 does for scale), or make it a per-time-point array. Failing is enough if no such case exists today.
3. **`fold_change_recomputed_by_kg` records a fix, not a limitation.** Keeping it in `data_caveats` is fine, because the warning is useful: a user checking a value against the paper will see a different number. But the description should say so: "…so values will not match the paper's table." It should also say whether the original values can be found anywhere. Under §2 there is no `source_fold_change`, so the answer today is "only in the paper".
4. **`significance_calls_differ_from_paper`: "a simpler rule"** doesn't fit Christie-Oleza 2017 S7, where the KG calls 32 *more* genes than the paper. "A different rule" covers every case.
5. **`shares_probe_with` detection.** "Identical log2FC in ≥ N cells" only works where both twins have edges. In a significant-only array table, the unlisted twin has no edge, and that's exactly the case the twin rule exists for. The flag on the *listed* twin's edge is what tells the explorer that the other one is undefined. So the pairs must come from a probe-to-gene map, or from all-genes tables on the same platform, applied to every table on that platform. Please make the builder rule platform-wide, not per table.
6. **The `paralogs_share_one_array_probe` caveat vs the edge flag.** If `shares_probe_with` is on the edges, the caveat should be set by rollup ("any edge in this experiment has `shares_probe_with`"), not curated separately. Otherwise the two can drift apart, which is the same reasoning 2b applies to the low-signal count.

### 6.4 Q4, coverage: what the explorer needs from the KG

The call matrix (2d.6) must say, for a gene with no edge in an experiment, one of: "not significant", "not reported (filtered)" or "undefined (not measured)". `table_scope` already settles most cases. What remains is the measured universe of significant-only and top-n tables:

| Platform | Absent from a significant-only table means | What the KG needs |
|---|---|---|
| RNA-seq | not significant (every gene was measured) | nothing new, if the platform is known |
| Microarray | not significant **if on the array**, otherwise undefined | which genes the array covers |
| Proteomics | undefined (detection is not exhaustive) | nothing |

Minimum request:
- `Experiment.measurement_platform` (str, vocabulary), e.g. `rnaseq`, `affymetrix_md4_9313`, `agilent_…`, `itraq`. `omics_type` alone doesn't separate array designs.
- Later (a separate item): the gene universe per array design. It could be derived from the all-genes tables on the same design, the same data the twin detection would use (6.3 item 5).

Until then, the explorer will report microarray absences as "not significant or not on the array", which is honest.

### 6.5 Q5, process

- **At integration, the validator prints** the §4.2 statistics for each new analysis. This is a report, not a warning, consistent with §4.2.
- **The paperconfig must state `data_caveats`** (`[]` allowed) and `data_caveats_reviewed` (6.2) for every new experiment, so a new paper can't land unreviewed.
- **Old experiments are reviewed in batches.** The review date makes progress visible, and the explorer can report coverage, e.g. "112 / 340 experiments reviewed".
- **Explorer side:** when absorbing a rebuild (the `absorb-kg-rebuild` skill), list experiments whose statistics stand out but that have no caveat and no review date. Those go into a list for the KG, not into user warnings.

### 6.6 §2 non-goals: one exception

We'd still like **F1 `fold_change_type`** soon. It was the question that started all this ("there are several fold-change types; does the MCP support them?"), and it costs a single stamp of a field the paperconfig already declares, plus a vocabulary entry. It isn't a limitation, so it doesn't belong in `data_caveats`. If you'd rather keep it in the later provenance spec, we'll close it on our side as deferred.

### 6.7 What the explorer will build

1. **`list_experiments`:**
   - Rows get `data_caveats` and `data_caveat_note` (the note in verbose only) and `data_caveats_reviewed`.
   - The envelope gets `by_data_caveat`.
   - Filters: `data_caveats` (any of) and `exclude_data_caveats`.
   - Verbose rows also carry the §4.2 statistics.
2. **DE tools, enrichment tools and `gene_response_profile`:**
   - A warning naming each experiment in scope that has caveats, with the caveat descriptions.
   - `low_signal_reason` on DE rows (sparse).
   - The envelope counts low-signal rows.
   - Possibly an `exclude_low_signal` filter. The default stays off; the analysis decides.
3. **`de_table_qc` (backlog 2d.5)** reads the caveats, the statistics, the review date and target-gene presence. It contains no thresholds of its own.
4. **The call matrix (backlog 2d.6)** uses `table_scope`, `measurement_platform` once it exists, and `shares_probe_with`.
5. **`docs://analysis/expression`** gets a "Before you trust a DE table" section.
6. **Testing:** a `ControlledVocabulary` coverage test for `data_caveats`, `low_signal_reason` and `low_signal_basis`.

### 6.8 Effect on the explorer ask

Once this spec is final:
- F7 is replaced by §4.2 and F8 by §4.1.
- F2–F6 are deferred to the possible provenance spec.
- The Aharonovich recompute in our ask's §5.4 is withdrawn (line B).
- F1 depends on 6.6.

We'll update the ask and backlog item F to point here.

## 7. KG response to the explorer review (2026-11-18)

### 7.1 Accepted as written
- **§6.1:** 2b is now **Decided**.
- **§6.2:** reviewed marker added to §4 item 1, named `data_caveats_reviewed_on` (an ISO date) so the label says what the value is.
- **§6.3.1:** renamed to `time_point_low_signal_count`; rule 1 stays without a fold-change cutoff. Rule 2 keeps its \|log2FC\| > 5 cutoff, since without it the bottom-10% rule counts 10% of every table.
- **§6.3.2:** the build fails on mixed level-column kinds within one experiment.
- **§6.3.3, §6.3.4:** descriptions changed ("will not match the paper's table … the original values are only in the paper"; "a different rule").
- **§6.3.5 (superseded by §10):** twin detection is per array design and stamped on every table of that design. The KG has no probe-to-gene maps, so pairs come from the all-genes tables on the same design. That needs the design recorded: see 7.2.
- **§6.3.6 (superseded by §10):** `paralogs_share_one_array_probe` is set post-import from the edge flags and never curated.

### 7.2 Accepted with a change
- **(superseded by §10)** **§6.4 `measurement_platform` → `array_design`**, sparse, on microarray experiments only. Arrays are the only case where the platform changes what an absence means; for RNA-seq vs proteomics, `omics_type` already decides. One field fewer to curate on ~190 experiments, and no second field that repeats `omics_type`.

### 7.3 Effect on the work
The curation pass grows by two fields per experiment (`data_caveats`, `data_caveats_reviewed_on`), plus level-column declarations for ~11 papers. (`array_design` removed, §10.)

### 7.4 Process (Q5), accepted from §6.5
- At integration the validator prints the §4 item 2 statistics for each new analysis, computed from the resolved table (the same rules as post-import). A report, not a warning.
- **Decided (2026-11-19): allowlist.** A committed file (e.g. `config/data_caveats_unreviewed.yaml`, name left to the builder) lists the existing DE experiments not yet reviewed. The validator requires `data_caveats` (`[]` allowed) and `data_caveats_reviewed_on` on every DE experiment **not** on the list, so every new paper is covered from day one and today's paperconfigs still validate. It also errors when a listed experiment already carries a review date (stale entry), so the list only shrinks. Each review batch removes its experiments; the backfill is tracked in `plans/backlog.md` and is done when the list is empty, after which the file and its check are deleted.
- Existing experiments are reviewed in batches; `data_caveats_reviewed_on` shows progress.
- The explorer's rebuild-absorb list of "stands out, not reviewed" experiments comes back to the KG as a review queue.

### 7.5 F1 `fold_change_type` (§6.6): **Decided (2026-11-19): in this spec**
- `Experiment.fold_change_type`: categorical, one of `log2` | `linear` | `signed_linear` (the paperconfig's existing values, defaulting to `log2`). Dense on every experiment with a DE contrast, absent otherwise. Adapter-set by stamping the analyses' declared value; the build fails if one experiment's analyses declare different scales (none do today).
- **Meaning:** the scale of the table the KG read, before the adapter converted it to `log2_fold_change`. For `_modified.csv` papers (Fuszard, Ma 2022, Koch 2019) that is the KG's derived table, not the paper's; the vocabulary description says so. The stored `log2_fold_change` is always log2 whatever this says; the field tells a reader how to compare it with the paper.
- Registered in `config/controlled_vocabularies.yaml` with `value_descriptions` and declared in `schema_config.yaml`. Not a caveat: it never goes in `data_caveats`.

### 7.6 First-delivery scope: **Decided (2026-11-19): (a)**
The first delivery is the mechanism plus the cases this spec names:
- Code: vocabularies, `schema_config.yaml`, adapter (`fold_change_type`, `low_signal_reason`, low-signal count, `data_caveats` stamping), post-import (statistics, rollups, `[]` re-materialisation), validator (vocabulary, note, allowlist, level-column existence, and the §7.7 QA report: statistics, candidate level columns, caveats to consider), the `/paperconfig` skill's QA step (§7.7), unit + KG tests, docs.
- Curation: each experiment named in the §4 item 1 table, checked against its paper before tagging, with a review date.
- **Not** in it: level-column declarations for the ~11 papers in §4.4 (so `time_point_low_signal_count` / `low_signal_reason` are empty until a later batch declares them), the §3.4 suspects, the backfill of the allowlist, and anything blocked on the missing Aharonovich supplement. All are follow-ups in `plans/backlog.md`.
- The explorer surface ships with the first build that carries `data_caveats` (§6.1); it will show fewer flags than the full curation would, not a different shape.

### 7.7 Paperconfig fields and the `/paperconfig` skill (added 2026-11-19)
**New paperconfig fields** (all optional in the schema; the allowlist in §7.4 decides when the first three are required):
- Experiment: `data_caveats` (list), `data_caveat_note`, `data_caveats_reviewed_on`.
- Statistical analysis: the level-column lists of §4.4 (`treatment_level_cols` / `control_level_cols` or `mean_level_cols`; names left to the builder).

**`scripts/validate_paperconfig.py` gains a QA report** for each DE analysis. It is printed, never a warning or error, and it never fails validation:
- The §4 item 2 statistics per time point, computed with the shared per-row functions below, the same code the adapter runs. It reads `<stem>_resolved.csv` when it exists, else the source CSV, so it already works before the gene-ID steps run: the statistics need no gene IDs. Unresolved rows make the numbers differ slightly from the build's.
- When level columns are declared, the low-signal count and the column-existence check. When they are not, it lists header columns that look like levels (`baseMean*`, `*RPKM*`, `*NSAF*`, `*CPM*`, `*count*`, `*mean*`, `*expression*`, `*peptide*`) as candidates to declare.
- **Caveats to consider** (prompts for the curator, not verdicts): extreme fold changes or a low-signal count above 0 → `fold_change_inflated_at_low_counts`; a negative fraction above 0.75 or below 0.25, or a median far from 0, on a full table → `conditions_not_normalized_to_each_other`; a `_modified.csv` builder for the table → `fold_change_recomputed_by_kg` / `fold_change_sign_reversed_by_kg`; a significant count far from the paper's stated number (from `table_scope_detail`) → `significance_calls_differ_from_paper`. The thresholds are prompt cutoffs only and are never stored (R3).

**One implementation of the per-row calls, shared by adapter and validator (decided 2026-11-19).** A new module (e.g. `multiomics_kg/utils/de_calls.py`, name left to the builder) holds, as pure functions with no adapter state:
- reading a row's fold change, p-value and asterisk marks from the declared columns (including the empty-p-value rule: untested is not significant);
- the scale conversion to log2 (today `_convert_fold_change` in `omics_adapter.py`);
- the significance call (today `OMICSAdapter._check_significance`; its defaults become arguments);
- the low-signal rule of §4 item 2 (`condition_at_column_minimum` / `mean_in_bottom_decile`, empty cell = no signal);
- the per-time-point table statistics, computed from those per-row results.

`omics_adapter.py` and `validate_paperconfig.py` both call it, so the report cannot drift from the build. The move is behaviour-preserving for the adapter: `/omics-edge-snapshot` before and after must match, and the existing adapter tests keep passing unchanged. The post-import Cypher stays authoritative for the stored statistics; a KG test checks that the Python statistics equal the stored arrays on a sample of experiments, so the two implementations can't diverge silently.

**`/paperconfig` skill** (`.claude/skills/paperconfig/SKILL.md`) gains a QA step after validation:
1. For each DE table, propose the level columns from the header candidates; the curator confirms against the paper's Methods.
2. Run the validator and read its QA report.
3. Walk the caveat vocabulary with the paper open, starting from the prompts. Set `data_caveats` (`[]` allowed), a `data_caveat_note` when non-empty, and `data_caveats_reviewed_on`.
The same steps are the checklist for the allowlist backfill (§7.4).

## 8. Explorer follow-up on §7 (2026-11-19)

**All of §7 is accepted,** including the change from `measurement_platform` to `array_design`. That narrowing is right: `omics_type` already settles RNA-seq vs proteomics, and the explorer's call matrix will take its array-only path when `array_design` is present. F1 in this spec closes the question that started all this. One gap remains, plus two housekeeping notes.

### 8.1 Gap: sign flips made by a builder are flagged nowhere

F1 now describes "the table the KG read", so for `_modified.csv` papers it describes the KG's derived table. `fold_change_recomputed_by_kg` covers builders that recompute values. Builders that only **negate** the paper's sign are covered by neither:
- **Ferrer 2020:** "Table sign was control minus co-culture and is negated in the builder".
- **Barber 2026 ×3:** "Perseus Difference was glucose minus substrate and is negated".
- **Koch 2019:** its builder "negates the paper's 'X vs Y' log2 fold changes".

The flip is in `design_notes`, which is free text. A user checking a KG value against the paper's table sees the opposite sign and finds no structured reason for it. That's the same check `fold_change_recomputed_by_kg` exists to support. Proposal:
- Add one more value, `fold_change_sign_reversed_by_kg`, with the description: "The paper's table reports the comparison the other way round (control vs treatment). The KG negated it so that positive means up in the treatment; signs are opposite to the paper's table."
- Curate it on the papers above, and on any other builder that negates.
- Alternatively, widen `fold_change_recomputed_by_kg` to "the KG changed the published values (recomputed or sign-reversed; see note)". We prefer the separate value, because the description can then say exactly what changed.

### 8.2 Housekeeping

- **Status labels are out of date.** The header says "Draft, in dialog", §2 says "Non-goals (proposed)" and §4 says "(outline, open)", but every open question is now decided. Please mark the spec **Decided** / frozen once 8.1 is settled, so the explorer can freeze its own spec against it.
- **§5 numbering:** item 8 sits between 5 and 6.

### 8.3 Explorer next steps

- We'll update our ask (`2026-11-18-kg-side-de-value-provenance-ask.md`) and backlog item F to point here. The ask becomes: F1, F7 and F8 are done by this spec, F2–F6 are deferred, and the recompute in its §5.4 is withdrawn.
- The explorer surface (§6.7) waits for the first KG build that carries `data_caveats`, and ships with it, as agreed in §6.1.

## 9. KG response to §8 (2026-11-19)

- **§8.1 accepted as proposed:** separate value `fold_change_sign_reversed_by_kg`, added to the §4 item 1 table with the explorer's description. A grep of `scripts/build_modified_csv/` confirms only Ferrer 2020, Barber 2026 and Koch 2019 negate; Domínguez 2017 assigns a sign from "highest mean condition", which matches the paper's direction, so it is not flagged. A new builder that negates must add the value (checked at review, not by code).
- **§8.2 done:** status set to Decided, section titles updated, §5 renumbered.
- **§8.3:** noted. The KG side of the deferred items (F2–F6, per-design gene universes, the backfill) is in `plans/backlog.md`.

## 10. Amendment: microarray probe features removed (2026-11-19)

`Experiment.array_design`, edge `shares_probe_with` and the caveat `paralogs_share_one_array_probe` are **removed from this spec**. This supersedes §4 items 2b (twin part) and 3 (`array_design`), §6.3.5–6.3.6 and §7.1–7.2 where they concern these.

**Why.** Inferring shared probes from published fold changes needs every table mapped to a design, enough full tables per design and a matching rule, all for a handful of paralog pairs on old arrays. The reliable source is a probe list (the GEO platform table, e.g. GPL5471's 9,947 rows, or the probe lists some supplements carry), not inference. `array_design` alone only served that pooling. A GEO lookup also showed that a GEO platform id does not identify an array (MD4-9313 is GPL5471, GPL11412 and GPL19747), so even the field's value needs curation.

**For the explorer.** Key the array path of the call matrix (2d.6) on `omics_type = 'MICROARRAY'`, not on `array_design`: a gene missing from a microarray table reads "not significant or not on the array". Shared-probe twins are not flagged; a curator who notices them in a paper may say so in `data_caveat_note`.

**Kept for later:** backlog item "microarray probe-level features" in `plans/backlog.md`, with the GEO findings.

## 11. Explorer acknowledgment of §7.7 and §10 (2026-11-19)

- **§10 accepted.** Inferring shared probes is the weak link, and a GEO id that doesn't identify an array settles it. The explorer's call matrix (2d.6) keys its array path on `omics_type = 'MICROARRAY'` ("not significant or not on the array").
- **Twin detection stays in the explorer, as an inference.** Backlog 2d.6 already planned an identical-log2FC pair check within each organism's microarray tables. It will be labelled in its output as inferred from values, not from a probe map, and kept as a sensitivity option the analysis can turn on, never as a default call. No KG ask follows. If the probe-level backlog item lands later, the explorer switches to it.
- **§7.7: no explorer impact.** Sharing the per-row module between adapter and validator, and testing it against the stored arrays, is the right way to keep the report honest. One small note: the "far from the paper's stated number" prompt parses `table_scope_detail` free text. That's fine as a curator prompt, and it's also one more reason for F5 (`paper_reported_significant_count`) when the deferred provenance spec comes up.

## 12. KG note on §11 (2026-11-19)

Agreed on all three points; nothing changes in this spec. The F5 point is added to the F2–F6 item in `plans/backlog.md`: `paper_reported_significant_count` would replace the validator's parse of `table_scope_detail` free text for the significance-count prompt.
