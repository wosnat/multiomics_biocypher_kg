# Signed linear fold changes stored as log2 — scope

Backlog item: "Signed linear fold changes are stored as log2" (found 2026-10-18).
Status: DONE 2026-10-31 on branch `fix/fold-change-conventions` — commit 1 (cd3b7d43) scale fixes; commit 2
significance overhaul (all Class P / D / Z decisions below implemented, see "Commit 2"). Docker rebuild +
`/omics-edge-snapshot` + `pytest -m kg` pending in the other clone.

## Problem

Four papers publish signed linear ratios ("+3 = 3-fold up, -3 = 3-fold down"). Their
paperconfigs set no `fold_change_type`, so `Changes_expression_of.log2_fold_change` holds the
linear number. Two consequences:

1. **Magnitudes are wrong** (Kaur 2018 reaches -8,463.9; a "log2" of -8,464 is a 2^8464-fold change).
2. **Significance ignores fold change.** The adapter's default `logfc_threshold` is 1.0, applied
   as `|value| >= 1`. On signed linear data every value has `|x| >= 1`, so the threshold is a
   no-op and significance is decided by the p/q-value alone. After conversion the same
   threshold means >= 2-fold, which changes `expression_status` on many edges.

Live KG (queried 2026-10-31):

| Paper | DOI | Edges | Significant now | Significant but < 2-fold (would flip at threshold 1.0) | Range |
|---|---|---|---|---|---|
| Kaur 2018 | 10.1111/1462-2920.14012 | 5,698 | 1,658 | 12 | -8,463.9 … 353.9 |
| Oleza 2017 | 10.1038/nmicrobiol.2017.100 | 3,338 | 760 | 375 | -337.3 … 47.3 |
| Martiny 2006 | 10.1073/pnas.0601301103 | 878 | 236 | 171 | -3.06 … 162.3 |
| Steglich 2006 | 10.1128/JB.01097-06 | 1,188 | 386 | 0 | -11.3 … 61.2 |

None of the four has a value in (-1, 1), as expected for signed linear data.

## In scope

1. **New convention `fold_change_type: signed_linear`**: log2 = sign(x) · log2|x|.
   Exactly ±1 → 0. **Values in (-1, 1) fail loudly** (decided 2026-10-31). The convention cannot
   hold them (no change is ±1), and converting one would invert its direction (0.5 → -1,
   -0.5 → +1). So one such value means the table is mis-declared, not that a row is bad. The
   validator reports an error after reading the whole column, and the adapter logs an error and
   builds no edges for that analysis. Rows are never silently dropped. The shared check is
   `multiomics_kg/utils/fold_change.signed_linear_violations`. None of the four papers has such a
   value.
   - `omics_adapter._convert_fold_change`: add the branch.
   - `scripts/validate_paperconfig.py`: allowed values become `log2 | linear | signed_linear`.
   - `_validate_fc_range`: for `signed_linear`, warn if any value is in (-1, 1). For `log2`, add a
     heuristic warning when no value is in (-1, 1) **and** max|x| > 20. That pattern catches all four
     papers; it is what would have caught this.
2. **Fix the falsy-threshold bug**: `_check_significance` uses `analysis.get('logfc_threshold') or
   self.default_logfc_threshold`, so `logfc_threshold: 0` silently becomes 1.0. Same for
   `pvalue_threshold: 0`. Switch to `is not None`. This is a prerequisite for item 3.
3. **Per-paper significance criterion.** Read each paper's methods and set `logfc_threshold` to the
   criterion the paper actually used, rather than inheriting the 1.0 default:
   - Steglich 2006: legend says ">= 2-fold", already `logfc_threshold: 1.0` → 0 flips. Correct as is.
   - Kaur 2018, Oleza 2017, Martiny 2006: check the methods. If the paper calls a gene significant
     on q/p alone, set `logfc_threshold: 0`. That keeps today's significant set and fixes only the
     magnitudes. If the paper also requires a fold cut, set its log2 equivalent (1.5-fold → 0.585).
4. **Set `fold_change_type: signed_linear`** in the 4 paperconfigs, rewriting the
   "no fold_change_type — signed linear" comments (Oleza 2017 has four of them).
5. **Audit for other cases.** Run the new heuristic over every paperconfig and look at whatever it
   flags. Two papers are already worth a look: Ostrowski 2010 PtrA (10.1038/ismej.2010.24, range
   -3.39 … 2.63, no value in (-1, 1), all 160 significant) and the two `significant_only` tables
   with no value in (-1, 1) (mSystems.00070-20, s43705-022-00197-2). They are probably log2 tables
   pre-filtered at |log2FC| >= 1, but that needs checking against the paper. Record the verdict
   either way.
6. **Ma 2022 "Mean Ratio"** (10 proteomics values: up 1.21–1.31, down -0.55 … -0.83). Here the
   minus sign marks a ratio < 1; it is not a signed fold. Recommendation: **no new convention**.
   Write a `_modified.csv` with |x| and set `fold_change_type: linear`, so the down values become
   log2(0.55…0.83) = -0.86 … -0.27. Confirm against the paper's table legend first.

7. **Aharonovich 2016 HOT1A3: ratios inflated by a 0.01 pseudocount.** This is a different
   defect, added 2026-10-31. The two HOT1A3 tables (`HOT1A3 coculture with 9313 ver coculture
   with 9313 diluted.csv`, `HOT1A3 coculture with MED4 vs coculture with MIT9313.csv`) are
   correctly declared `linear` and are unsigned (0 negatives). The problem is how the authors
   computed the ratio: `(a + 0.01) / (b + 0.01)` on raw counts, with the +0.01 columns included in
   the CSV. When one side is zero the ratio becomes 100 × count + 1, e.g. 20 vs 0 → 2001 and
   1 vs 0 → 101 (log2 6.7). The repeated values 101, 201, 301 … are exactly these rows.

   | Table | Rows | One side zero | Both zero | Significant (q ≤ 0.01 and ≥ 2-fold) | Of those, one side zero | log2 range now | log2 with +1 pseudocount |
   |---|---|---|---|---|---|---|---|
   | 9313 vs 9313 diluted | 1,966 | 721 | 38 | 386 | 88 | -9.0 … 18.1 | -5.2 … 11.5 |
   | 9313 vs MED4 | 2,381 | 224 | 0 | 118 | 9 | -12.0 … 12.2 | -5.4 … 5.6 |

   Zero-count genes fill all of the top 100 |FC| rows in both tables (median non-zero count 9
   and 4 reads), so `rank_by_effect` for these two experiments ranks sequencing depth rather
   than effect size.
   Recommendation: a builder script (the `scripts/build_modified_csv/` precedent) writes
   `_modified.csv` with log2((a + 1) / (b + 1)) from the expression columns. The conventional
   +1 pseudocount bounds a 1-vs-0 gene at log2 1. Both-zero rows are dropped (no measurement),
   and the table is declared `log2`. q-values and significance stay the paper's.
   The alternative is to keep the paper's ratios and accept a meaningless top-N. Either way, record
   the before/after edge count and significant count, and add a CHANGELOG `### Data` bullet. The
   three Aharonovich tables already in log2 (MED4 / MIT9313 DE genes) are unaffected.

## Not in scope

- `rank_by_effect`: sign(x)·log2|x| keeps the |x| order within each experiment × time point, so
  ranks don't change. `rank_up` / `rank_down` move only where `expression_status` moves (item 3).
  The post-import Cypher is untouched.
- Treating 0 as "down" in the existing `expression_direction` rule (`fc > 0 → up`). This is a
  separate, pre-existing edge case. Note it in the backlog if the conversion produces any zeros.
- MCP / explorer changes. The property keeps its name and meaning (log2), and these values only
  become correct.
- Re-deriving p-values or re-running any DE.

## Acceptance criteria

- Unit tests: `_convert_fold_change` signed_linear (±3 → ±1.585, ±1 → 0); a value in (-1, 1) refuses the analysis (adapter) and is an error (validator); the
  validator accepts `signed_linear`; `logfc_threshold: 0` is honoured; the log2 heuristic fires
  on a Kaur-like fixture and not on a normal log2 fixture.
- `/omics-edge-snapshot` before/after. Edge counts are unchanged for all papers. Only
  `expression_status` changes, and only on the papers where item 3 deliberately sets a fold
  threshold. The counts are recorded here.
- Live check: max |log2_fold_change| corpus-wide drops from 8,463.9 to about 13 (log2 8,464 ≈ 13.05).
- `pytest -m "not slow and not kg"` and `pytest -m kg` both green. Regenerate `snapshot_data.json`
  if it samples edges from these papers.
- CHANGELOG `### Data`: one correction bullet per paper saying what the KG asserted before and
  what it asserts now.

## Decisions to confirm

1. Convention name `signed_linear`. Recommended; it matches the paperconfig comments.
2. Significance follows each paper's own criterion (item 3), not a uniform 2-fold default.
   Recommended. The alternative (keep 1.0 → ≥ 2-fold everywhere) would turn 558 of the
   currently-significant edges in Oleza 2017 and Martiny 2006 to not significant, with no basis
   in either paper.
3. Ma 2022 handled by a modified CSV, not by a fifth convention.
4. Aharonovich HOT1A3: recompute log2 with a +1 pseudocount and drop both-zero rows (recommended),
   or keep the paper's 0.01-pseudocount ratios.

## Progress (2026-10-31)

Offline adapter dry run, `HEAD` paperconfig vs working tree (`OMICSAdapter.get_edges`):

| Paper | Edges | Significant | max \|log2FC\| | Change |
|---|---|---|---|---|
| Kaur 2018 | 5,698 → 5,698 | 1,658 → 1,658 | 8,463.90 → 13.05 | `signed_linear` + `logfc_threshold: 0` (Methods: t-test P = 0.05, permutation FDR q = 0.05, no fold cutoff) |
| Christie-Oleza 2017 | 3,338 → 3,338 | 760 → 760 | 337.30 → 8.40 | `signed_linear` + `logfc_threshold: 0` (Fig. 2: t-test P < 0.05). Sign checked against NSAF/NAF means: positive = higher in coculture |
| Martiny 2006 | 878 → 878 | 236 → 236 | 162.29 → 7.34 | `signed_linear` + `logfc_threshold: 0` (legends: q < 0.05) |
| Steglich 2006 | 1,188 → 1,188 | 386 → 386 | 61.20 → 5.94 | `signed_linear`. Significance comes from the asterisks, so `logfc_threshold` is moot |
| Ma 2022 | 168 → 168 | 158 → 168 | (S4 4.75 → 2.25) | S4 `_modified.csv` (log2\|x\|) + `prefiltered: true`. The 10 down proteins were "not significant" because the default 2-fold was applied to \|x\| < 1 |
| Aharonovich 2016 | 10,423 → 10,392 | 1,123 → 1,123 | HOT1A3 18.11 → 11.46, 12.23 → 5.61 | HOT1A3 `_modified.csv` log2((a+1)/(b+1)); 38 both-zero rows blanked (31 of them had resolved to a gene) |

Christie-Oleza 2017 S7: the table's own "Significant" (+) column is Perseus's FDR + s0 call (263 rows),
a strict subset of q ≤ 0.05 (295). The 32 extra rows are documented in the paperconfig, not fixed.

New builders: `scripts/build_modified_csv/build_ma2022_modified_csv.py`,
`build_aharonovich2016_hot1a3_modified_csv.py`. Both re-resolved with
`resolve_paper_ids --papers ... --force`. At commit: `git rm` the now-unused tracked `_resolved.csv` +
`_resolved_report.txt` of the three original tables (Ma S4 and the two HOT1A3 tables), and add the
builders to `scripts/build_modified_csv/README.md`.

### Corpus audit (item 5)
Six report-only agents, 38 papers (the six above excluded). Brief and reports are in the session scratchpad
`fc_audit/`. Findings are verified here before any fix:
Two problem classes beyond the original four papers:

**Class S (scale): wrong `fold_change_type`.** Fixed when verified.
- Lindell 2007: **VERIFIED + FIXED.** The Supp. Table 3 legend (supp PDF) reads "Fold change (infected/control) …
  Positive and negative values indicate an increase and decline". Now `signed_linear` on all 8 analyses.
  Dry run: edges 360 → 360, significant 117 → 117 (asterisks), max 14.00 → 3.81. The warning missed it
  (max 14 < 20), so the audit, not the warning, is the safety net. Four values of exactly +1.0 → 0 →
  `down` (the tie rule, see Not in scope).

**Class P (significance): the config does not reproduce the paper's criterion.** Not fixed yet:
these go to the user as one decision, because they change `expression_status`, not magnitudes.
- Hennon 2017 Supp. Table 2 (MIT9312): **VERIFIED.** The paper: "only 35 genes were differentially
  expressed (P < 0.05, fold change > 2)" (p. 6). The table is the 103 P < 0.05 rows, but it is marked
  `prefiltered: true`. Proposed: `prefiltered: false` + `pvalue_threshold: 0.05` + `logfc_threshold: 1.0` → 26
  significant. The paper's 35 is not reproducible from the table.
- Biller 2016 MIT1002 Table S3 (batch B): only up genes are published. The table lacks `prefiltered`, so 282 + 285
  edges read not significant. Proposed: `prefiltered: true` + table_scope `filtered_subset`.
- Tai 2009 (batch B): significant-only SAM tables, no fold cutoff, no `prefiltered` → 475 of 574
  not significant. Proposed: `prefiltered: true`.
- Capovilla 2023 (batch B): significant-only DESeq2 tables (padj < 0.1), no `prefiltered` → 11/25 + 28/72
  not significant. Proposed: `prefiltered: true`.

- Tolonen 2006 + Thompson 2011 (batch E): asterisk levels. The adapter treats any `*` as significant
  and returns before thresholds. Tolonen's `*` = 0.01 < q < 0.05 but the paper's DE criterion is q < 0.01
  (`**`); `**` alone reproduces the paper's 81 MED4 / 69 MIT9313 genes (the KG has 113 / 110) → 658
  over-called. Thompson's `*` = q < 0.01 with fold < 2, the paper requires both (`**`); its declared
  `logfc_threshold: 1.0` is ignored on the asterisk path → 152 over-called. Proposed: an adapter key
  `asterisk_min_level: 2` (count the stars; fewer than N = not significant). One key fixes both, with no
  modified CSVs. Tolonen PMM1118 (hli4): two probesets interleaved in the original xls, so every
  column but t = 0h is wrong. Proposed: drop the row via a modified CSV, or NA it.

- **Default 2-fold applied to papers with no fold cutoff (batches D, F)**, the most common case. q/p-only
  papers where the KG's `default_logfc_threshold = 1.0` silently removes paper-DE genes: Weissberg 2025
  (34 analyses; authors' code `is_de = padj < 0.05`; the KG keeps ~30–70% of the paper's set), Hackl 2023
  (UV vs ctrl 886 → 117), Johnson 2026b (~half; 897/121 reproduced at `logfc_threshold: 0`), Ferrer 2020
  (−34 of 776), Read 2017 (BH q < 0.05; q rounded to 2 decimals, so `pvalue_threshold: 0.049`).
  Al-Hosani 2015 is the reverse: `prefiltered: true` on a whole-genome table marks all 1,935 significant.
  The paper uses raw p < 0.05 with no fold cutoff → 627 with `pvalue_threshold: 0.049` + `logfc_threshold: 0`
  (raw p in the adjusted field; `padj` ≤ 0.05 would give 483).
- Batch C, the same no-fold-cutoff class: Moreno 2023 (P ≤ 0.05 and ≥ 2 unique peptides, no fold cutoff; 1,063 of 2,760
  paper-significant rows read not significant), Lin 2015 (DESeq padj < 0.1, no fold cutoff → 215 vs 163; S4
  is a curated 34-gene P set incl. padj = 1 rows, so `table_scope: filtered_subset`), Domínguez-Martin 2017
  (all 408 rows are the paper's set → `prefiltered: true`; the 2-fold default keeps 155).
- **Adapter bug, empty p-value (batches D, F):** with `adjusted_p_value_col` declared, an empty cell (DESeq2
  NA padj) gives `pval = None` and `_check_significance` skips the p test. The row is then significant on fold
  alone: 200 edges in Weissberg 2025 (59 vs 25 real in MED4 day 18 vs 11) and 3 in Koch 2019. With
  `logfc_threshold: 0` all 1,543 Weissberg NA rows would be significant, so **this fix must ship before any
  `logfc_threshold: 0` on DESeq2 tables.** Proposed: declared column + empty cell → `not significant`
  (tables with no p column keep the fold-only rule). Checked 2026-10-31: none of the papers fixed above
  (Kaur, Oleza 2017, Martiny, Steglich, Aharonovich, Lindell) has an edge without a p-value; Ma 2022 S4 is
  `prefiltered`.
- **Policy question raised by the above:** should `default_logfc_threshold` stay 1.0? At least 9 papers used
  no fold cutoff, and at least 4 others needed an explicit 1.0. The alternative is no default: every analysis
  states its fold cutoff (validator error when missing).

**Class D (source-data errors):**
- Fuszard 2012 NATL2A: **VERIFIED.** The supplement's "Ratio of means" = Deplete Mean / `# peptides` on 63 of 63 rows
  (= Deplete / Control on 0). The KG stores 61 false down-regulations (published < 1: 62 rows; true: 33).
  The paper's 6 up / 1 down are reproduced by Deplete / Control. Fix: `build_fuszard2012_modified_csv.py` computes
  log2(Deplete Mean / Control Mean) for NATL2A. MIT9312 and SS120 are correct. Also: these tables list every
  quantified protein, so `table_scope: all_detected_genes`.
- Pandhal 2007: 18 single-peptide proteins have no p-value, so the empty-p bug calls them significant on fold.
  The paper requires > 1 peptide (p. 1000). Fix: the empty-p adapter fix, or drop them in the builder.
- Tolonen 2006 PMM1118 (hli4), see above.

**Class Z (zero-count artifacts), the Aharonovich kind:** Beliaev 2014 (zero-read genes floored at RPKM ≈ 0.12 drive
all top-50 effects in S3/S6/S7/S8, up to 498-fold; 97 genes rounded to 0.00 are skipped), Huang 2021 (batch F: unreplicated libraries, zero-control genes top the effect ranks; the paper's rule gives 47 up, not its 32 — UNSURE; the KG also applies a mirrored down rule the paper never defines), Capovilla 2023 MIT9303 ±19.7–21.6 (unshrunk
DESeq2 against a zero arm, Dataset S5 at the normalized-count floor), Anjur-Dietrich 2026 +24.45/+26.13
(no count table to confirm). Undecided: document, or cap in a `_modified.csv`.

OK: Ostrowski 2010, Barreto 2022, Durham 2015, Manck 2020 (A); Bagby 2015, Thompson 2016 (B); Doron 2016
(caveat already documented, D7(a)), Biller 2018, Tetu 2009, Tetu 2019 (E); Coe 2024 (minor: S4 table_scope,
README P<0.1), Oleza 2015b, Fang 2019 (D); Barber 2026, He 2022, Koch 2019 (3 NA-padj edges), Kratzl 2024 (F).
Reports: scratchpad `fc_audit/batch_{A..F}.md`. All 38 papers audited.

## Commit 2 — significance overhaul (2026-10-31, user-approved)

Implemented: (1) adapter no default `logfc_threshold` (validator error / analysis refused); 180 analyses
migrated to the paper's value; (2) empty p-value cell → not significant; (3) `asterisk_min_level`
(Tolonen, Thompson 2011 = 2); (4) per-paper `prefiltered` / threshold / scope fixes (Hennon S2, Al-Hosani,
Read, Lin, Biller MIT1002, Tai, Capovilla, Domínguez, Coe S4 scope); (5) Fuszard builder recomputes
Deplete / Control; Tolonen PMM1118 dropped; (6) Beliaev re-derived from RPKM with a one-read pseudocount;
(7) log2FC = 0 → no `expression_direction`. Capovilla / Anjur / Huang zero-count magnitudes documented
only: Capovilla's S4/S5 counts are VST-transformed (zero = 9.85), so there are no raw counts to recompute.

Corpus-wide offline adapter run, commit 1 vs commit 2: edges 358,915 → 358,998 (Beliaev +91, Tolonen −8);
significant 102,557 → 113,725 (up 33,971 → 39,689, down 68,586 → 74,036); 1,510 edges without direction.
Per-paper numbers are in CHANGELOG `### Data`. All 61 paperconfigs pass `validate_paperconfig.py`.
