# Gene ID Mapping Utility — Implementation Plan

## Context

Gene IDs in supplementary paper CSVs are heterogeneous: the same gene appears as a locus tag, NCBI RS-tag, Cyanorak ID, protein accession, gene symbol, microarray probeset, or paper-specific reannotation ID. This causes ~dangling edges in the knowledge graph and makes LLM/MCP searches unreliable.

The current system partially addresses this via `gene_mapping_supp.csv` (heuristic alt-ID extraction from CSVs) and `gene_id_utils.py`. This plan extends that foundation into a comprehensive, structured gene ID mapping system with three goals:

1. **Quality**: omics adapter emits only valid locus_tag-resolved edges (no dangling)
2. **Coverage**: explicit column declarations in paperconfig capture all known alt-IDs, including from GFF/GTF reannotation files
3. **Searchability**: a per-strain `gene_id_mapping.json` enables fast forward/reverse lookup for MCP/LLM agents

Design principles:
- **All new paperconfig fields are optional** — existing configs keep working via heuristics
- **Incremental**: adding `id_columns` to any paper immediately enriches that strain's mapping
- **Backward compatible**: `gene_mapping_supp.csv` continues to be generated so check/fix-gene-ids still work

---

## Phase 1: Extend paperconfig.yaml Schema

### New optional fields per supplementary table

Each CSV contains data for a single organism. The organism is either declared directly on the table via `organism`, or inferred from the `statistical_analyses` entries that reference this table.

Two new table `type` values are supported alongside the existing `csv`:

```yaml
supplementary_materials:
  # ── Existing type: DE/expression data CSV with optional new fields ──
  table_key:
    type: csv
    filename: "path/to/data.csv"
    sep: ","         # NEW — optional; column delimiter, default "," (use "\t" for TSV)
    skip_rows: 0
    organism: "Prochlorococcus MED4"   # NEW — optional; inferred from analyses if absent
    original_filename: "data/.../original.csv"   # NEW — optional; pre-fix-gene-ids CSV

    # NEW — optional; if absent, heuristic column detection is used (backward compat)
    id_columns:                        # gene identifier columns in this CSV
      - column: "NCBI ID_3"
        id_type: locus_tag_ncbi        # one of: locus_tag, locus_tag_ncbi, locus_tag_cyanorak,
      - column: "Gene Name"            #   gene_name, protein_id_refseq, uniprot_accession,
        id_type: gene_name             #   old_locus_tag, alternative_locus_tag, gene_synonym,
                                       #   uniprot_entry_name, jgi_id, probeset,
                                       #   rast_id, annotation_specific, other

    product_columns:                   # columns with functional gene descriptions
      - column: "Genbank Annotation"
      - column: "Gene description"

    statistical_analyses:              # unchanged
      - id: "unique_id"
        organism: "Prochlorococcus MED4"
        name_col: "NCBI ID_3"
        ...

  # ── NEW type: pure ID translation table (no DE data, no statistical_analyses) ──
  gene_id_translation:
    type: id_translation
    filename: "data/.../gene_id_table.csv"
    sep: ","         # optional; default ","
    skip_rows: 0
    organism: "Prochlorococcus MED4"   # required for this type
    id_columns:                        # required for this type
      - column: "NCBI locus tag"
        id_type: locus_tag_ncbi
      - column: "MED4 gene name"
        id_type: gene_name
      - column: "Protein accession"
        id_type: protein_id_refseq
    product_columns:                   # optional
      - column: "Gene product"

  # ── NEW type: paper-specific GFF/GTF reannotation (ID bridging only) ──
  reannotation_gff:
    type: annotation_gff
    filename: "data/.../reannotation.gff"
    organism: "Prochlorococcus MED4"   # required for this type
```

The `annotation_gff` type:
- Has no `statistical_analyses`, `id_columns`, or `product_columns`
- Requires `organism` and `filename` to be declared
- Is processed by `build_gene_id_mapping.py` to bridge paper-specific gene IDs to canonical locus tags
- Is NOT processed by omics_adapter (no edges emitted directly)
- Skipped by fix-gene-ids and check-gene-ids

The `id_translation` type:
- Has no `statistical_analyses` (it's purely a mapping resource)
- Requires `organism` and `id_columns` to be declared
- Is processed by `build_gene_id_mapping.py` the same way as id_columns in a `csv` table, but is NOT processed by omics_adapter (no edges emitted)
- Skipped by fix-gene-ids and check-gene-ids (not an expression data source)

### Strain-level shared resource paperconfigs

Some ID translation files are not associated with any specific paper — they apply to all papers for a given strain (e.g., a strain-wide locus tag cross-reference table). These belong in a **strain resource paperconfig**: a `paperconfig.yaml` with no `publication` block, containing only `id_translation` and/or `annotation_gff` entries.

**Placement:** `data/Prochlorococcus/papers_and_supp/<Strain>_resources/paperconfig.yaml`
**Registration:** listed in `paperconfig_files.txt` like any other paperconfig.

Example for MIT9313:
```yaml
# Strain-level shared gene ID resources for Prochlorococcus MIT9313
# No publication block — this is not a paper
supplementary_materials:
  mit9313_id_translation:
    type: id_translation
    filename: "data/Prochlorococcus/papers_and_supp/MIT9313_resources/MIT9313_genbank.tsv"
    sep: "\t"
    organism: "Prochlorococcus MIT9313"
    id_columns:
      - column: "NCBI locus tag"
        id_type: locus_tag_ncbi
      - column: "Old locus tag"
        id_type: old_locus_tag
      - column: "Gene name"
        id_type: gene_name
```

`build_gene_id_mapping.py` processes these the same as paper-embedded `id_translation` entries; the `paper` provenance field in `paper_ids` is set to the filename (e.g., `"MIT9313_resources"`) rather than a paper name.

omics_adapter ignores strain resource paperconfigs entirely (no `publication` block → no publication node emitted; no `statistical_analyses` → no edges emitted).

**Files to modify:**
- Documenting in: `CLAUDE.md` (Adding Omics Data section) and `/skill:paperconfig` SKILL.md
- No code changes required for schema itself (it's YAML, consumed by new build script)

---

## Implementation Progress

| Paper | Status | Notes |
|-------|--------|-------|
| Anjur 2025 | ✅ paperconfig updated, 🔄 script in progress | `id_translation` for `9301_annotated_genome.csv` (JGI IDs via `uniprot_gene_name` whitespace-split); `annotation_gff` for NCBI GFF |
| barreto 2022 | pending | `annotation_gff` for `EZ55.exon.fixed2.gtf` (GTF format); `id_translation` for annotation tables |
| all others | pending | phase 2 curation |

### Anjur 2025 — JGI catalog IDs

**Problem:** DE CSV (`de_genes_freeliving_vs_biofilms.csv`) uses JGI catalog IDs (integer strings like `2626311769`) as `name_col: ID`. These are organism-specific IDs from the JGI IMG database and not standard locus tags.

**Mapping chain:** `9301_annotated_genome.csv` contains:
- `ID` column: JGI catalog ID (integer) — `id_type: jgi_id`
- `uniprot_gene_name` column: `"dnaA P9301_05911"` or `"P9301_05921"` (gene symbol + locus_tag, or locus_tag alone)
- `uniprot_entry_name` column: `"DNAA_PROM0"` (strip `_PROM0` → try as accession or gene_name)
- `em_Preferred_name` column: `"dnaA"` — `id_type: gene_name`

**Resolution strategy:** Declare `uniprot_gene_name` as primary anchor column. Build script tries full value first → fails → splits on whitespace → last token is locus_tag (e.g., `P9301_05911`). This token resolves directly via `build_id_lookup()` (it is in `old_locus_tags` for MIT9301). All 242 DE genes are in the annotation CSV, so 100% resolution is expected.

**Verification:** `P9301_05911 in lookup? P9301_05911` ✓ (confirmed via gene_id_utils.py).

### barreto 2022 — GTF + annotation tables

**Problem:** Multiple organisms; EZ55 DE tables use locus-tag-like `symbol` values (e.g., `EZ55_00798`) that already resolve. The main addition is:
1. `annotation_gff` for `EZ55.exon.fixed2.gtf` (GTF format, `gene_id "EZ55_00001"` attributes)
2. `id_translation` for three organism annotation tables providing UniProt accessions and product synonyms

**Annotation table structure:**
- `pro_9312_anot.csv`: `GID` (PMT9312_XXXX), `SYMBOL`, `PRODUCT`, `uniprot_acc`, `KO`
- `syn_9311.csv`: `GID` (sync_XXXX), `genename`, `definition`, `uniprot_acc`, `KO`
- `syn_8102_anot.csv`: `GID` (SYNW_XXXX), `SYMBOL`, `PRODUCT`, `uniprot_acc`, `KO`

GTF gene_id values (`EZ55_00001`) are already canonical locus tags — GTF adds no novel IDs but confirms the annotation source.

---

## Phase 2: Curate All 22 Existing Paperconfigs

For each paper in `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt`, inspect all supplementary CSVs and add `id_columns` and `product_columns` declarations.

### Handling pre-fix vs. post-fix state

Many papers have already been through fix-gene-ids, so their paperconfig now shows:
```yaml
filename: "results_with_locus_tag.csv"   # updated by fix-gene-ids
name_col: "locus_tag"                    # updated by fix-gene-ids
```

The `_with_locus_tag.csv` retains all original columns — nothing is removed. But the original `name_col` (e.g., `"NCBI ID_3"`) is no longer recorded in the config. Curation strategy:

1. **For post-fix papers**: inspect the `_with_locus_tag.csv` to discover all original columns (it still has them). The original primary ID column should be declared in `id_columns` with its `id_type` — even though `name_col` is now `"locus_tag"`. The build script uses `id_columns` as its source for alt-IDs regardless.

2. **For pre-fix papers** (name_col is not "locus_tag"): the CSV is the original, declare `id_columns` normally.

3. **Optional `original_filename` field**: if a paper still has its original CSV available (always true — fix-gene-ids never deletes the original), we can add:
   ```yaml
   original_filename: "data/.../results.csv"   # pre-fix CSV, for reference
   ```
   The build script can read this to discover original column structure when `_with_locus_tag.csv` is the current filename.

**Strain resource paperconfigs to create** (new files, not existing papers):
- `MIT9313_resources` — strain-wide locus tag translation CSV valid for all MIT9313 papers; create `data/Prochlorococcus/papers_and_supp/MIT9313_resources/paperconfig.yaml` with `id_translation` entry, add to `paperconfig_files.txt`

**Papers with known multiple ID columns** (priority):
- `Biller 2018` — NCBI ID, NCBI ID_2, NCBI ID_3, Gene Name, RAST annotation (scope: annotation_specific)
- `coe 2024` — Gene ID, NCBI ID, NCBI ID_2, NCBI ID_3, Gene Name, product
- `bagby 2015` — Locus tag2, Alternative locus tag, Probeset, Gene description
- `Fang 2019` — Gene ID (tuple format)
- `barreto 2022` — add a new `annotation_gff` entry: `filename: "data/Prochlorococcus/papers_and_supp/barreto 2022/EZ55.exon.fixed2.gtf"`, `organism: "Alteromonas EZ55"`
- `Anjur 2025` — two entries needed:
  1. `id_translation` for `9301_annotated_genome.csv` (organism: MIT9301):
     - `uniprot_entry_name` column: UniProt entry names (id_type: `uniprot_entry_name`);
       build script strips `_PROM0` suffix → resolves as uniprot_accession (for unnamed proteins
       like "A3PBU0_PROM0" → "A3PBU0") or gene_name (for named ones like "DNAA_PROM0" → "DNAA")
     - `ID` column: JGI gene catalog IDs (id_type: `jgi_id`, scope: organism_specific) — becomes
       alt-ID once the row is resolved via uniprot_entry_name
     - `em_Preferred_name`: gene name (id_type: `gene_name`)
     - Process BEFORE the DE CSV so JGI IDs are in the lookup
  2. `annotation_gff` for `GCF_000015965.1_ASM1596v1_genomic.gff` (organism: as declared)
  - DE CSV `de_genes_freeliving_vs_biofilms.csv`: name_col is `ID` (JGI IDs);
    resolves after id_translation is processed

All 22 paperconfigs get curated. Papers with single ID column still get an `id_columns` entry for the primary `name_col` to document it explicitly.

---

## Phase 3: Create build_gene_id_mapping.py

**New file:** `multiomics_kg/download/build_gene_id_mapping.py`

### What it produces

`cache/data/<Organism>/genomes/<Strain>/gene_id_mapping.json`:

```json
{
  "PMM0001": {
    "locus_tag": "PMM0001",
    "alt_ids": {
      "reference": [
        {"id": "TX50_RS00020",     "id_type": "locus_tag_ncbi",        "scope": "organism_specific"},
        {"id": "CK_Pro_MED4_00001","id_type": "locus_tag_cyanorak",    "scope": "organism_specific"},
        {"id": "WP_011131639.1",   "id_type": "protein_id_refseq",     "scope": "organism_specific"},
        {"id": "Q7V6L1",           "id_type": "uniprot_accession",     "scope": "organism_specific"},
        {"id": "dnaN",             "id_type": "gene_name",             "scope": "generic"},
        {"id": "beta-clamp",       "id_type": "gene_synonym",          "scope": "generic"},
        {"id": "PMM0001",          "id_type": "old_locus_tag",         "scope": "organism_specific"},
        {"id": "P9313_RS00020",    "id_type": "alternative_locus_tag", "scope": "organism_specific"}
      ],
      "paper_ids": [
        {"id": "PMED4_00071",        "id_type": "old_locus_tag",  "scope": "organism_specific", "source_col": "Alternative locus tag", "source_csv": "de_genes_gas_shock_vs_air.csv", "paper": "Bagby and Chisholm 2015"},
        {"id": "MED4_ARR_0008_x_at", "id_type": "probeset",       "scope": "organism_specific", "source_col": "Probeset",             "source_csv": "de_genes_gas_shock_vs_air.csv", "paper": "Bagby and Chisholm 2015"},
        {"id": "fig|59919.17.peg.1", "id_type": "rast_id",        "scope": "annotation_specific", "source_col": "RAST annotation", "source_csv": "...", "paper": "Biller 2018"}
      ]
    },
    "product_synonyms": [
      {"text": "DNA polymerase III beta subunit", "source": "ncbi"},
      {"text": "beta clamp",                      "source_col": "Gene description", "paper": "Bagby and Chisholm 2015"}
    ]
  }
}
```

### ID scope classification

Each alt-ID carries a `scope` field to guide LLM/MCP search disambiguation:

| `scope` | Meaning | Examples |
|---------|---------|---------|
| `organism_specific` | Uniquely identifies one gene in one genome | locus_tag, locus_tag_ncbi, old_locus_tag, protein_id, UniProt accession, probeset |
| `generic` | Shared across organisms; search returns multiple hits | gene_name ("dnaN" exists in many bacteria) |
| `annotation_specific` | Unique within a reannotation system but not standard | RAST fig\|...\|peg.N IDs |

This allows LLM agents to know: "PMM0001 → exactly one gene" vs. "dnaN → search across all organisms".

Each entry is keyed by canonical `locus_tag`. This structure supports:
- Reverse lookup (alt_id → locus_tag) for the adapter/fix-gene-ids tools
- Forward lookup (locus_tag → all known IDs and names) for MCP/LLM "find gene" queries
- Scope-aware disambiguation for LLM search ("is this ID organism-specific or generic?")

Also outputs the backward-compat `gene_mapping_supp.csv` (same format as today) from the same data.

### Algorithm

```
For each strain:
  1. Load gene_annotations_merged.json → base entries per locus_tag
     → Seed alt_ids.reference from: locus_tag_ncbi, locus_tag_cyanorak, protein_id (RefSeq WP_),
       uniprot_accession, gene_name, gene_synonyms, old_locus_tags, alternative_locus_tags
     → Apply scope: gene_name / gene_synonyms → "generic"; all others → "organism_specific"
     (gene_annotations_merged.json contains all available UniProt accessions; no separate UniProt join needed)

  2. For each paperconfig that has entries for this organism:
     Process `id_translation` tables BEFORE `csv` tables so the enriched lookup is available
     when resolving the DE data CSVs.

     a. For each supplementary table of type `csv` or `id_translation`:
        i.  Determine which CSV to read:
            - If `original_filename` declared → read that for original alt-ID columns
            - Else read `filename` (may be `_with_locus_tag.csv`, which still has all original columns)
            - Apply `skip_rows` and `sep` (default ",")
        ii. Determine resolution strategy per table type:
            For `csv`: name_col drives resolution
            - If name_col == "locus_tag" (post-fix state) → rows already resolved; use locus_tag column directly
            - Else → resolve name_col values via build_id_lookup() (with whitespace-token splitting; see below)
            For `id_translation`: no name_col; try each declared id_column in order
            - For each row, try each id_column's value via build_id_lookup(); first column that resolves
              provides the locus_tag anchor; remaining column values become paper_ids for that locus_tag
            - After processing all id_translation tables, rebuild the lookup to include the newly added
              paper_ids before processing csv tables
        iii. Special id_type handling during resolution:
            - `uniprot_entry_name`: strip the trailing `_ORGANISM` suffix (last underscore + rest)
              to obtain the accession-or-gene-name prefix, then resolve via build_id_lookup() as
              uniprot_accession (for unnamed proteins, e.g. "A3PBU0_PROM0" → "A3PBU0") or
              gene_name (for named proteins, e.g. "DNAA_PROM0" → "DNAA"); try accession first
            - Whitespace-split fallback: if a column value contains spaces and the full value doesn't
              resolve, split on whitespace and try each token individually (e.g., "dnaA P9301_05911"
              → tries "dnaA" then "P9301_05911"); use the first token that resolves
            - This fallback handles gene_name columns that concatenate symbol + locus tag
        iv. Determine ID columns for alt-ID extraction:
            - If `id_columns` declared in paperconfig → use those (with declared id_type and scope)
            - Else → heuristic fallback (is_id_like_column() from gene_id_utils.py), exclude "locus_tag" column
            - Infer scope: RAST/fig|... → "annotation_specific"; gene names (dnaN, rpsA) → "generic"; else "organism_specific"
        v. For each row: collect alt-ID values from id_columns → append to paper_ids with provenance
           - For multi-token columns, record each whitespace-split token as a separate paper_id entry
        vi. Collect product_columns values for each resolved gene → append to product_synonyms
     b. For each supplementary table of type `annotation_gff`:
        i.  Organism is taken directly from the entry's required `organism` field
        ii. Parse GFF/GTF; extract old_locus_tag + Name + ID + gene_id attributes
            → Bridge: for each GFF feature, look up any extracted attribute via build_id_lookup()
            → Record novel (not already in reference alt_ids) ID → locus_tag pairs as paper_ids

  3. Deduplicate: collapse identical (id, id_type, source_col, source_csv, paper) tuples
  4. Write gene_id_mapping.json
  5. Write gene_mapping_supp.csv (alt_id, locus_tag, source_col, source_csv, paper) for backward compat
```

### GFF/GTF ID bridging

When a supplementary table of type `annotation_gff` is present:

**Parser:** Use Biopython's `BCBio.GFF` (already a project dependency) to parse GFF3, and Biopython's `SeqIO` for GTF. Follow the same merge-CDS-into-gene logic used in `cyanorak_ncbi_adapter.py`: read `gene` features first, then merge `CDS` attributes into the parent gene by locus_tag.

**Target attributes** (mirroring the NCBI adapter):
- `locus_tag` — primary canonical ID (identity match)
- `old_locus_tag` — legacy locus tags (NCBI GFF maps RS-format → legacy PMM-style via this)
- `Name` — gene symbol
- `ID`, `gene_id`, `gene_name` — other GFF3/GTF identifier fields
- `protein_id` — CDS-level cross-reference (from merged CDS features)

**Both GFF3 and GTF are supported.** Detect format by file extension (`.gff`, `.gff3` → GFF3; `.gtf` → GTF) or by inspecting first non-comment line.

**Algorithm:**
1. Parse file → merged gene+CDS records (locus_tag keyed)
2. For each gene record: attempt to resolve `locus_tag` attribute → canonical locus_tag via `build_id_lookup()`
3. Collect all other extracted attributes as alt-IDs; mark `old_locus_tag` as `id_type: old_locus_tag, scope: organism_specific`
4. Record novel pairs (not already in reference alt_ids) as `paper_ids` entries

---

## Phase 4: Update gene_id_utils.py

**File:** `multiomics_kg/utils/gene_id_utils.py`

Add three new functions:

```python
def load_gene_id_mapping(genome_dir):
    """Load gene_id_mapping.json. Falls back to gene_annotations_merged.json if absent."""

def build_id_lookup_from_mapping(gene_id_mapping):
    """Build alt_id → locus_tag reverse index from gene_id_mapping.json.
    Covers all alt_id types: reference fields + paper_ids."""

def search_genes_by_name(genome_dir, query, fields=None):
    """Full-text search across gene names, synonyms, and product text.
    Returns list of {locus_tag, match_field, match_value} dicts.
    Used by MCP/LLM agents."""
```

Update `build_id_lookup()` to call `load_gene_id_mapping()` preferentially, falling back to the current approach when the JSON doesn't exist yet.

---

## Phase 5: Integrate into prepare_data.sh

**File:** `scripts/prepare_data.sh`

Add a new step (after Step 2 — build_gene_annotations):

```
Step 3: build_gene_id_mapping.py — builds gene_id_mapping.json per strain
```

Command:
```bash
uv run python multiomics_kg/download/build_gene_id_mapping.py
```

With `--force` / `--strains` flags matching existing step conventions. Log to `logs/prepare_data_step3.log`.

Also add to CLAUDE.md's prepare_data.sh documentation.

---

## Phase 6: Pre-resolve paper CSVs + omics_adapter.py probing

**Approach:** 2-step (transparent, auditable — like the existing fix-gene-ids pattern).

### Step A — `multiomics_kg/download/resolve_paper_ids.py` (new script, prepare_data step 4)

For each `csv` supplementary table whose `name_col` is not already `locus_tag`:
- Loads the source CSV with the declared `sep` (default `,`)
- Resolves each `name_col` value via `build_id_lookup()` (uses `gene_id_mapping.json`)
- Writes `<stem>_resolved.csv` alongside the original, adding two columns:
  - `locus_tag` — resolved canonical ID (NaN when unresolved)
  - `resolution_method` — how it was resolved (`direct`, `lookup`, `supp`, `repadded`, `unresolved`, `empty`)
- Skips tables already using `name_col: locus_tag`; skips if `_resolved.csv` is newer than source (unless `--force`)
- Prints a resolution report grouped by **publication / organism** at the end

```
======================================================================
Gene ID Resolution Report
Overall: 235/242 resolved (97.1%)
======================================================================

  Publication: Anjur 2025
    Organism: Prochlorococcus MIT9301
      [supp_table_1]  name_col=ID
        Resolved:  235/242 (97.1%)
        Output:    data/.../de_genes_freeliving_vs_biofilms_resolved.csv
        Unresolved (7):
          row 5: 2626311769
          ...
```

Usage:
```bash
uv run python -m multiomics_kg.download.resolve_paper_ids [--force] [--papers "Anjur 2025"]
```

### Step B — `omics_adapter.py` (minimal changes)

In `get_edges()`:
- Extract `sep` and `organism` from `table_data`; inject into `stat_analysis` (like `skip_rows`)
- Pass `sep` to `_load_and_create_edges(filename, analysis, sep=sep)`

In `_load_and_create_edges(filename, analysis, sep=',')`:
- Probe for `<stem>_resolved.csv` next to `filename`
- If found: read it with `sep`, set `use_locus_tag_col = 'locus_tag' in df.columns`
- In the row loop: if `use_locus_tag_col`, read `row['locus_tag']` instead of `row[name_col]`
  - NaN locus_tag rows are already skipped by the existing null-check — no extra logic needed
- Fallback (no `_resolved.csv`): unchanged behaviour (use `name_col` directly)

The adapter emits no edges for unresolved rows (they have NaN locus_tag in `_resolved.csv`).
The resolution audit trail is in the `_resolved.csv` files, not in the adapter.

### When to run step 4

Step 4 must run after step 3 (so `gene_id_mapping.json` is up to date with paper-derived IDs).
Default `STEPS` in `prepare_data.sh` is now `0 1 2 3 4`.

```bash
# Re-resolve one paper after editing its paperconfig
bash scripts/prepare_data.sh --steps 4 --force   # re-resolves all papers
uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Anjur 2025" --force
```

---

## Phase 7: Update paperconfig Skill

**File:** `.claude/skills/paperconfig/SKILL.md` (and `paperconfig.py` if it exists)

Add documentation and interactive prompts for:
- `id_columns` section: for each CSV, ask which columns contain gene IDs and what type
- `product_columns` section: which columns contain functional descriptions
- `annotation_gff` entry type: is there a GFF/GTF from paper-specific genome reannotation? If so, add it as a separate `supplementary_materials` entry with `type: annotation_gff` and `organism`

---

## Phase 8: Update check-gene-ids and fix-gene-ids

Both skills already use `build_id_lookup()` from `gene_id_utils.py`. Once Phase 4 is done (build_id_lookup prefers gene_id_mapping.json), they automatically benefit.

Add to SKILL.md: mention that running `/build-gene-mapping-supp` (or the new `build_gene_id_mapping.py`) refreshes the mapping before running these skills.

---

## File Summary

| File | Action |
|------|--------|
| `multiomics_kg/download/build_gene_id_mapping.py` | **CREATE** — core build script |
| `multiomics_kg/utils/gene_id_utils.py` | **EXTEND** — 3 new functions |
| `multiomics_kg/adapters/omics_adapter.py` | **EXTEND** — pre-resolve gene IDs |
| `scripts/prepare_data.sh` | **EXTEND** — add Step 3 |
| `data/Prochlorococcus/papers_and_supp/*/paperconfig.yaml` (×22) | **CURATE** — add id_columns/product_columns |
| `.claude/skills/paperconfig/SKILL.md` | **EXTEND** — document new fields |
| `.claude/skills/build-gene-mapping-supp/SKILL.md` | **UPDATE** — mention superseded by new build script |
| `CLAUDE.md` | **EXTEND** — new schema docs + Step 3 |

---

## Verification

```bash
# 1. Build mapping for one strain
uv run python multiomics_kg/download/build_gene_id_mapping.py --organism "Prochlorococcus MED4"
# Check output
python -c "import json; d=json.load(open('cache/data/Prochlorococcus/genomes/MED4/gene_id_mapping.json')); print(len(d), 'genes'); print(list(d.items())[:2])"

# 2. Verify backward-compat supp CSV still generated
head cache/data/Prochlorococcus/genomes/MED4/gene_mapping_supp.csv

# 3. Run check-gene-ids on a previously problematic paper
uv run python .claude/skills/check-gene-ids/check_gene_ids.py --paperconfig "data/Prochlorococcus/papers_and_supp/bagby 2015/paperconfig.yaml"

# 4. Rebuild graph and confirm match rate improvement
uv run python create_knowledge_graph.py --test
# Check logs for "gene_id resolved via gene_id_mapping" counts

# 5. KG validity tests (requires running Docker)
pytest tests/kg_validity/ -v
```
