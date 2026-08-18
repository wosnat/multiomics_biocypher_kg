# Scope: paperconfig updates from the GEO processed-file drop

Status: **scoped — decisions taken, implementation not started**
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

### 2.1 johnson 2026 a — GSE314951 iModulon paper

Richest of the set, all native IDs (`Refseq Locus` = `TX50_RS#####` → tier-1,
1,905 of 1,976 MED4 genes already in `specific_lookup`). Three separable pieces:

| piece | file | shape | entry type |
|---|---|---|---|
| A. DE | `media-2` sheets `Dark vs Light 2-4h (DESeq2)`, `Light vs Dark 14-16h (DESeq2)` | 1,920 rows; each sheet holds 2 timepoints as paired log2FC/p-adj columns → **4 contrasts** | `csv` x2 |
| B. periodicity | `media-3` `all_cosinor_fit` (1,872 rows: acrophase, amplitude, p_val, p_adj) + `fourier_padj_1e-2` | per-gene scalars + a peak label | `derived_metrics_table` |
| C. iModulons | `media-6`/`media-8` **M** matrix (1,872 genes x ~30 modules) + **A** activity matrix | module membership with weights | `gene_clusters` (`score_col`) |

- A is the standard path. B mirrors Waldbauer 2012 / Biller 2018. C mirrors the
  existing `gene_clusters` analyses but the M matrix is **wide** (one column per
  module) where `gene_clusters` expects a long `cluster_col` — needs a reshape
  step or a schema decision.
- `media-2` also carries **DIMA** sheets (26 rows = iModulon-level activity
  differences). No node type represents an iModulon's *activity change*; out of
  scope unless C lands first.
- Directory rename (`johnson 2026 a` to something meaningful) and stripping the
  `media-N(3).xlsx` browser suffixes should happen in the same commit.

**Acceptance:** decided per piece below; at minimum A gives 1 Publication,
1-2 Experiments, 4 analyses, ~7.7K edges on tier-1 ids.

---

## Tier 3 — blocked or needs its own plan

| item | blocker |
|---|---|
| **Hackl 2023** (GSE171435, MIT0604) | Strain **not in the KG**. Needs `/add-a-strain` (assembly, prepare_data, all per-strain tools) **plus** an ID bridge for the authors' custom `MIT0604_#####` annotation. The DE itself is clean and unfiltered (6 tables incl. antisense). Own plan file. |
| **Doron 2016** (GSE63690/GSE74921) | Microarray data is **probe-level** (`Probe.ID`, `alias`) so it needs probe-to-gene mapping. **WH8109 not in the KG.** RNA-seq processed files are **phage-only**. Host DE would come from the journal supplement (3 spreadsheets, not yet inspected). |
| **Steglich 2010** (GSE17075) | The per-gene half-life table is a **PowerPoint file**. Extraction is the whole job; GEO has CEL only. Would land as `derived_metrics_table`. |
| **munoz 2022** (GSE154594) | Environmental samples, not a strain. Needs probe-to-gene mapping off `GPL28884` and a reference-proteome-match-style organism. 129 MB RAW.tar not staged. |
| **Voigt 2014** (GSE53065) | TSS / operons / UTRs / ncRNAs have **no schema slot**. Skip unless a TSS layer is wanted. |

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

## Decisions taken

- **D1 — compute DE from count matrices? NO.** See Tier 4. Six datasets stay
  reference-only.
- **D2 — Johnson scope: A + B + C.** All three pieces, including iModulons as
  `gene_clusters`. C is the only genuinely new mechanic in this pass; see 2.1a.
- **D3 — Synechococcus registration: yes, in the Synechococcus registry.**
  Correction to an earlier note in this file: there are **two** registries and
  `create_knowledge_graph.py` loads both —
  `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` and
  `data/Synechococcus/papers_and_supp/paperconfig_files.txt` (the latter already
  registers Tal 2009, Beliaev 2014, Kratzl 2024, Ma 2022, Oleza 2015, Oleza 2017,
  kaur 2018, Bernstein 2017). Huang 2020 and Doron 2016 go in the Synechococcus
  list, not the Prochlorococcus one. `multiomics_kg/utils/paperconfig_utils.py`
  holds the canonical pair in `PAPERCONFIG_LISTS`.
  Side note: CLAUDE.md's "Reads three registry files" section names only the
  Prochlorococcus list — worth correcting when we next touch that doc.

### 2.1a iModulons as gene_clusters (the one new mechanic)

`gene_clusters` expects a long table: one row per gene, with a `cluster_col`
naming that gene's cluster. The iModulon **M** matrix is wide — 1,872 genes x
~30 module columns, each cell a gene weight — and iModulon membership is
**not mutually exclusive** (a gene can carry weight in several modules).

Two viable shapes, to settle before writing the config:
1. **Melt to long** with a threshold on the weight (standard iModulon practice
   uses the M-matrix threshold from the ICA fit) and emit one
   `Gene_in_gene_cluster` edge per (gene, module) above threshold, carrying the
   weight as `score_col`. Preserves multi-membership; needs a derived CSV.
2. **One analysis per module** — 30 ClusteringAnalysis nodes. Faithful but
   inflates the node count and misrepresents a single ICA fit as 30 analyses.

Prefer (1). Note this makes `Gene_in_gene_cluster` many-per-gene for this
analysis, which existing consumers may assume is 1:1 — verify before building.
`media-5` `tables_key` and `media-6`/`media-8` `table_key` sheets document the
matrix conventions and the threshold; read them first.

## Explicitly NOT in scope

- No schema changes, no new node/edge types.
- No new strains in this pass (MIT0604, WH8109 deferred to Tier 3 plans).
- No post-import Cypher changes.
- No re-running `prepare_data.sh` steps 0-2 (no new genomes).
- No DE computed by us (D1) — Tier 4 stays reference-only.
- Steps 3/4 (`build_gene_id_mapping`, `resolve_paper_ids`) **will** need re-running
  for every paper touched.

## Working notes

- Files are committed **gzipped**. `pd.read_csv` reads `.gz` transparently, but
  `resolve_paper_ids.py` writes `<stem>_resolved.csv` — a `.txt.gz` source yields
  an awkward stem. **Decompress the specific files we wire up**; leave the rest
  gzipped.
- `GSE118155` count CSVs use **CR-only** line endings.
- `GSE195946` and the Hackl tables carry **URL-encoded** text (`%2C`); decode
  before it reaches a node property (and note `clean_text()` silently rewrites
  the pipe and apostrophe characters).
