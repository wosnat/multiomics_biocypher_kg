# Paper README template — KG data inventory

Used by subagents to produce one `README.md` per paper directory under `data/Prochlorococcus/papers_and_supp/<paper>/`. The README inventories every file in the directory and judges suitability for the knowledge graph.

---

## Goal

For each paper, answer: **what data is here, what is already in the KG, what else should be added, and what should move from one node type to another?**

The output is a file named `README.md` in the paper's directory. It overwrites any existing README (but should preserve accurate content from the old README).

---

## Template (write exactly this structure)

````markdown
# <Paper Name> (<Year>)

**Citation:** <one-line citation if known (authors, year, title, journal); leave "(unknown — infer from PDF)" if not extractable>
**DOI:** <doi or "(not found)">
**Organism(s):** <strains studied>
**Topic:** <2–4 sentences: what the paper does, what was measured, main experimental design>

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `<filename>` | <PDF/CSV/XLSX/DOCX/TSV/GFF/FASTA/other> | <one-line description of what's in the file> | <status> | <what to do; "—" if no action> |
| ... | | | | |

### Status values (use exactly one per row)
- **already in** — data is integrated via the current `paperconfig.yaml`
- **add** — data is suitable for the KG but not yet integrated
- **change upload** — already integrated, but the node/edge type should change (e.g. `gene_clusters` → `derived_metrics_table`, or DE CSV should be split into multiple Experiments)
- **skip** — not suitable for the KG (sample-level data, raw read counts, figures, metadata without gene-level evidence, duplicate of another file already in)
- **reference** — PDF, figure, or legend; not data but kept for provenance

## Current paperconfig summary

<If paperconfig.yaml does not exist: write "No paperconfig.yaml — paper is not integrated.">
<If it exists, summarize in 4–8 bullets:>
- Experiments defined: <count + short list of keys>
- Statistical analyses (DE edges): <count; total rows / significant rows if visible>
- Supplementary materials entry types: <e.g., `csv`, `gene_clusters`, `derived_metrics_table`, `id_translation`, `annotation_gff`>
- Organisms covered: <list>
- Table scope(s): <e.g., `all_detected_genes`, `significant_only`, `significant_any_timepoint`>
- Non-DE evidence: <DerivedMetric entries, ClusteringAnalysis entries — counts>
- ID resolution: <native locus_tag? JGI? probeset? diamond-mapped? list bridge files used>

## Recommended actions

Ordered list of concrete next steps. Examples of the kind of item to write:

1. **Add** — Integrate `table_s5_periodicity.csv` as a `derived_metrics_table` entry with 4 boolean metrics (Y/N per condition). Attach to Experiment `<existing_key>`.
2. **Change upload** — Migrate `clusters.csv` from `gene_clusters` to `derived_metrics_table` (`value_kind: categorical`) because the table is a single composite classification column, not a soft-clustering assignment.
3. **Skip** — `raw_counts.xlsx` is sample × gene raw read counts; out of scope (gene-level evidence only).
4. **Add organism** — MIT9303 is referenced in the DE data but not yet deployed; add to `cyanobacteria_genomes.csv` and rebuild mapping.

Use "—" for items with no action ("PDF is reference-only, no action").

## Notes

- Gene ID format(s) used in the source tables and how they resolve (e.g., "PMM#### direct, PMED4_##### via JGI bridge, nan for ~5 rows")
- Strain(s) in/out of KG: call out missing genomes or reference-proteome-only strains
- Any anomalies: split-sheet CSVs, embedded asterisks marking significance, mixed units, etc.
- Links to existing related documentation (e.g., "See `docs/kg-changes/non-de-evidence-extension.md` for DerivedMetric rules")
````

---

## How to decide status values

### "already in" — integrated today
The file is listed as a supplementary_material entry in `paperconfig.yaml` (`csv`, `gene_clusters`, `derived_metrics_table`, `id_translation`, or `annotation_gff`). Accept this even if match rates are imperfect — the integration is present.

### "add" — missing but suitable
The file contains per-gene evidence the KG supports but no paperconfig entry references it. Common candidates:

- **DE table not yet wired up** → add as `csv` with `statistical_analyses`
- **Periodicity Y/N flags** → add as `derived_metrics_table` with `value_kind: boolean`
- **Phase, amplitude, lag, damping, timing scores per gene** → add as `derived_metrics_table` with `value_kind: numeric` (set `rankable: true` if directional)
- **Categorical per-gene tag (e.g. "early/mid/late responder", "essential/non-essential")** → `derived_metrics_table` with `value_kind: categorical` and `allowed_categories`
- **K-means or Mfuzz cluster assignments (one cluster per gene, soft-clustering membership scores)** → `gene_clusters`
- **New strain/genome** → call out; actual integration is genome deployment, not paperconfig

### "change upload" — rehome an existing integration
Most common case: a `gene_clusters` entry actually encodes classification flags or a composite label rather than a soft-clustering assignment. These should move to `derived_metrics_table`. Rule of thumb:

- **Use `derived_metrics_table`** when the column is: boolean Y/N, numeric score, or a categorical label from a small fixed vocabulary
- **Keep `gene_clusters`** when the column is: soft-clustering membership (each gene assigned to a single cluster of many, with a membership score), AND cluster descriptions exist (per-cluster name, functional_description, temporal_pattern)

Other change-upload scenarios:
- A single `csv` entry should be split into multiple Experiments (e.g., one per strain, one per timepoint family)
- Gene IDs are routed through `alternative_locus_tag` but should be `old_locus_tag` (or vice versa) based on what is in the CSV

### "skip" — out of scope
- Sample × gene raw counts or normalized expression matrices
- Sample metadata tables (sample IDs, collection dates, environmental readings)
- Figures, figure data, supporting text
- Duplicate formats of a file already listed (the Excel version of a CSV that is already in)
- Tables that are pathway-only / summary-only without per-gene evidence
- Metabolite tables (no gene target); flag if metabolite-to-gene mapping exists elsewhere

### "reference" — provenance only
PDFs, legends, figures. Mark as `reference` with action `—`.

---

## Process for each paper

1. `ls` the paper directory — enumerate every file (including nested dirs if present).
2. Read the main PDF briefly (title page + abstract + methods summary — usually first 3–5 pages is enough to identify organism, treatment, and data types produced). Skip reading figure PDFs.
3. Read `paperconfig.yaml` top-to-bottom if it exists. Note every `supplementary_materials` entry's `type` and which files it references.
4. For each file in the directory:
   - If it is referenced in the paperconfig → status `already in` (or `change upload` if the type is wrong)
   - If it is `_resolved.csv` / `_resolved_report.txt` → skip those rows (they are auto-generated; mention them once in Notes)
   - If it is a PDF, figure, legend → `reference`
   - Else decide between `add` and `skip` using the criteria above.
5. Keep the table rows in roughly the order: main PDF, DE tables, derived metrics / cluster tables, ID bridges, figures/legends, config files.
6. Fill the "Recommended actions" list with concrete verbs: **Add**, **Change upload**, **Skip**, **Add organism**, **No action**.

---

## Project conventions to carry forward

- File paths in the README table are relative to the paper directory (e.g., `table_S2.csv`, not the full repo-relative path).
- Keep the table compact — one row per file, content column is one line.
- Do not invent data that isn't there. If the PDF is the only file, say so.
- If a paperconfig exists, its `supplementary_materials` keys are the authoritative list of what's integrated. Use those, not file names, to judge "already in".
- Do not write TODO lists with checkboxes — use the "Recommended actions" numbered list.
- Keep `README.md` under ~200 lines total.

---

## Known status shortcuts (as of 2026-04-23)

- `Biller 2018` — already has detailed README; only rewrite if the inventory misses files.
- `Waldbauer  2012` — already has README; it is the canonical example of numeric DerivedMetric integration.
- `Capovilla 2023`, `zinser 2009`, `alonso 2023`, `coe 2024`, `Wang 2014` — have existing READMEs; preserve accurate content, rewrite in this template's format.
- `Labban 2022` — paperconfig commented out in `paperconfig_files.txt`. Status for its DE CSV is "add (blocked: waiting on annotation GFF)".
- `Kujawinski 2023`, `McDonagh 2012`, `Pandhal 2007`, `Szul 2019`, `biller 2022` — no paperconfig today. Most files become `add` candidates or `skip` (if metabolomics-only or sample-level).
- `MIT9313_resources` — reference resources directory, not a paper; produce a short README noting it holds a paperconfig that registers genbank-based ID translations.
- `QPCR` — one loose PDF; short README marking it as unparsed reference material.
