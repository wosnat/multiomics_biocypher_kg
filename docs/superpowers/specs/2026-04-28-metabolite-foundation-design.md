# Phase 1.1 — Metabolite foundation (data + resolver + step 2 enrichment)

**Date:** 2026-04-28
**Status:** Spec ready for implementation plan
**Parent overview:** [2026-04-28-metabolite-scaffold-design.md](2026-04-28-metabolite-scaffold-design.md)

## Goal

Set up the data and resolver infrastructure that downstream Phase 1 specs (Reactions Scaffold, Transport+CAZy) consume. **No KG-visible changes ship from this spec** — the deployed graph is unchanged after Spec 1.1 lands. What changes:

- Two new sub-steps inside the existing `prepare_data.sh` step 0 (sub-step 6: reference download; sub-step 7: build resolver + hierarchies). No top-level step is renumbered.
- An extension to the existing step 2 (`build_gene_annotations.py`) that makes it MNX-aware. Same step number, same script, same output file — new fields per gene.
- A new SQLite metabolite/reaction resolver (`cache/data/mnx/metabolite_resolver.db`).
- Two new JSON hierarchy caches (TCDB, CAZy).
- New / refined fields in `gene_annotations_merged.json` via the existing YAML transforms framework: `kegg_reactions` (MNXR-resolved, replaces raw `kegg_reaction`), `transporter_classification` (TCDB-validated), `cazy_ids` (CAZy-validated).

The success of this spec is verified entirely by inspecting the cache files and the per-strain metabolism report — not by querying Neo4j.

> **Schemas confirmed against MNXref 4.5 (2025-08-13) + TCDB live + eggNOG-mapper v2.x output.** SQLite tables (`compounds`, `compound_aliases`, `compound_names`, `reactions`, `reaction_aliases`) match real MNX file columns; TCDB hierarchy is built from three TCDB TSVs; CAZy is bootstrapped from eggNOG without a bulk download. Sub-step 6 + the integration test pass against real downloads. One open caveat: the dual short/long alias-source prefix forms in MNX (`chebi:`/`CHEBI:`, `kegg.compound:`/`keggC:`) need a normalization policy in Phase 1.1B (proposed: a hardcoded canonicalization map preferring the long bioregistry-style form). Reaction direction is **not** stored in MNX `reac_prop.tsv` — all 83,796 reactions use the ` = ` separator, so the SQLite resolver treats everything as reversible and the `direction_source` field is dropped. Spec 1.2's permissive direction-handling rollup is the default.

## Pipeline changes

No top-level step is renumbered. All additions land as new sub-steps inside the existing step 0 (which already mixes raw downloads with derivation — sub-step 5 builds `gene_mapping.csv`) plus an internal extension to the existing step 2.

### Step 0 — sub-step additions

| Sub-step | Existing/New | Purpose | Module | Output |
|---|---|---|---|---|
| 1 | existing | NCBI genome (GFF, FASTA, GBFF) | `download_genome_data.py` | `cache/data/<org>/genomes/<strain>/...` |
| 2 | existing | Cyanorak GFF/GBK | `download_genome_data.py` | as today |
| 3 | existing | UniProt per unique taxid | `download_uniprot.py` | as today |
| 4 | existing | eggNOG-mapper (skipped by default) | external | as today |
| 5 | existing | Build `gene_mapping.csv` | `build_gene_mapping.py` | `cache/data/.../gene_mapping.csv` |
| **6** | **new** | Download MNX, TCDB, CAZy reference data | `download_metabolism_reference.py` | `cache/data/mnx/{chem,reac}_{prop,xref}.tsv`, `cache/data/tcdb/families.tsv`, `cache/data/cazy/families.tsv` |
| **7** | **new** | Build metabolite/reaction resolver + hierarchy caches | `build_metabolite_resolver.py` | `cache/data/mnx/metabolite_resolver.db`, `cache/data/tcdb/tcdb_hierarchy.json`, `cache/data/cazy/cazy_hierarchy.json`, `metabolite_id_mapping_report.json` |

Sub-steps 6 + 7 follow the existing pattern: skip-if-cached unless `--force`. The dispatch and re-run mechanics:

**`download_genome_data.py` extension**

- `--steps` argparse `choices` extends from `[1, 2, 3, 4, 5]` to `[1, 2, 3, 4, 5, 6, 7]`; default extends accordingly.
- Two new top-of-`main` dispatch lines added next to the existing ones:
  ```python
  if 6 in steps:
      step6_metabolism_reference(force=args.force)
  if 7 in steps:
      step7_metabolite_resolver(force=args.force)
  ```
- The two `step6_*` / `step7_*` functions are thin wrappers that import and call into `download_metabolism_reference.py` and `build_metabolite_resolver.py`. Keeps the sub-step machinery in one place (the existing entry point) while letting the work live in dedicated files.
- Help epilog updated to list 6 + 7.

**`prepare_data.sh` extension**

The step-0 case extends the default `DOWNLOAD_SUBSTEPS` lists:

```bash
if [[ $SKIP_CYANORAK -eq 1 ]]; then
    DOWNLOAD_SUBSTEPS="1 3 5 6 7"
else
    DOWNLOAD_SUBSTEPS="1 2 3 5 6 7"
fi
```

Step 4 (eggNOG-mapper) stays excluded by default per existing convention.

**Granular re-run (avoids re-running NCBI / UniProt / Cyanorak)**

Step 0 is heavy. To rebuild only the metabolism caches without re-running the rest:

```bash
# Re-download MNX/TCDB/CAZy raw files and rebuild resolver:
uv run python multiomics_kg/download/download_genome_data.py --steps 6 7 --force

# Already-downloaded raw files; rebuild resolver only:
uv run python multiomics_kg/download/download_genome_data.py --steps 7 --force
```

Both bypass `prepare_data.sh` and call `download_genome_data.py` directly with sub-step selection. This is the same escape hatch already used today for refreshing only NCBI (`--steps 1 --force`) or only UniProt (`--steps 3 --force`).

### Top-level steps — unchanged numbering

| Step | Status | Notes |
|---|---|---|
| 0 | extended | Adds sub-steps 6 + 7 above. |
| 1 | unchanged | `build_protein_annotations.py` |
| **2** | **extended internally** | `build_gene_annotations.py` becomes MNX-aware: imports the resolver + hierarchy accessors, writes new fields per gene. Same step number, same script entry point, same output file (`gene_annotations_merged.json`). |
| 3 | unchanged | `build_gene_id_mapping.py` |
| 4 | unchanged | `resolve_paper_ids.py` |
| 5 | unchanged | `build_og_descriptions.py` |

A new top-level step for the pruned scaffold cache is introduced in Spec 1.2 (next available number).

**CLI stability:** existing `--steps 1 2`, `--steps 1 2 --strains MED4 --force`, etc. invocations keep working without changes. The only new surface is the two sub-steps inside step 0, reachable via the existing step-0 sub-step filter.

## Step 0 sub-step 6 — Reference data download

`multiomics_kg/download/download_metabolism_reference.py` (new module).

| Key | URL | Output | Actual size |
|---|---|---|---|
| `mnx_chem_prop` | `https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv` | `cache/data/mnx/chem_prop.tsv` | 773 MB |
| `mnx_chem_xref` | `https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv` | `cache/data/mnx/chem_xref.tsv` | 648 MB |
| `mnx_reac_prop` | `https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv` | `cache/data/mnx/reac_prop.tsv` | 9.8 MB |
| `mnx_reac_xref` | `https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv` | `cache/data/mnx/reac_xref.tsv` | 78 MB |
| `tcdb_families` | `https://www.tcdb.org/cgi-bin/projectv/public/families.py` | `cache/data/tcdb/families.tsv` | 130 KB |
| `tcdb_substrates` | `https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py` | `cache/data/tcdb/substrates.tsv` | 436 KB |
| `tcdb_superfamilies` | `https://www.tcdb.org/cgi-bin/substrates/listSuperfamilies.py` | `cache/data/tcdb/superfamilies.tsv` | 778 KB |

Total: ~1.5 GB (MNX dominates). Files cached and skipped on re-run unless `--force`.

**No CAZy bulk download.** CAZy's published `cazy_data.zip` (~51 MB zipped, 363 MB unzipped) is per-protein membership rows — useful for validation but adds ~550 families our 27 strains never reference, with no validation benefit (eggNOG uses dbCAN HMMs which ARE CAZy). Phase 1.1B bootstraps the CAZy hierarchy from observed eggNOG `CAZy` columns + mechanical ID parsing (`GH13_1` → family `GH13` → class `GH`). CAZy clans are deferred — they cross-cut family IDs and aren't mechanically derivable.

**Versioning:** the MNX release date is parsed out of the comment block at the top of `chem_prop.tsv` (`#VERSION: 4.5`, `#DATE: 2025/08/13`) and recorded in `metabolite_id_mapping_report.json` so future rebuilds detect drift.

**Gitignore (committed to repo):** `cache/data/{mnx,tcdb,cazy}/` are git-ignored — MNX is 1.5 GB, refresh on demand via:
```
uv run python -m multiomics_kg.download.download_genome_data --steps 6 --force
```

**Granular sub-step 7 only:** `--steps 7 --force` rebuilds the resolver from already-downloaded raw files without re-fetching anything.

## Step 0 sub-step 7 — Build resolver + hierarchy caches

`multiomics_kg/download/build_metabolite_resolver.py` (new module). One pipeline step, three outputs: SQLite resolver DB + two JSON hierarchies. Internal helpers (TSV parsers, hierarchy parsers) are private functions in this file — no separate helper module needed at this size.

### SQLite schema (`cache/data/mnx/metabolite_resolver.db`)

`compounds` and `reactions` tables mirror the columns of `chem_prop.tsv` and `reac_prop.tsv` (9 and 6 columns respectively in MNXref 4.5). `*_aliases` tables are normalized from the 3-column `chem_xref.tsv` / `reac_xref.tsv` — col 1 in those files is a combined `<source>:<value>` string, which the builder splits on the first `:` to populate the `source` and `value` columns.

```
compounds (
  mnxm_id     TEXT PRIMARY KEY,   -- "MNXM41"; col 1 of chem_prop.tsv
  name        TEXT,                -- col 2
  reference   TEXT,                -- col 3, e.g. "mnx:PROTON" or "chebi:16234" — provenance source
  formula     TEXT,                -- col 4
  charge      INTEGER,             -- col 5
  mass        REAL,                -- col 6
  inchi       TEXT,                -- col 7, full InChI string
  inchikey    TEXT,                -- col 8
  smiles      TEXT                 -- col 9
)

compound_aliases (
  source      TEXT,                -- normalized source prefix (see "Source vocabulary" below)
  value       TEXT,                -- e.g. "C00031" (KEGG), "16234" (ChEBI without prefix)
  mnxm_id     TEXT REFERENCES compounds(mnxm_id)
)
CREATE INDEX idx_compound_aliases_value ON compound_aliases(value);
CREATE INDEX idx_compound_aliases_source_value ON compound_aliases(source, value);

compound_names (
  name_normalized  TEXT,           -- lowercase, whitespace collapsed, punctuation stripped (except +/-)
  mnxm_id          TEXT REFERENCES compounds(mnxm_id)
)
CREATE INDEX idx_compound_names ON compound_names(name_normalized);

reactions (
  mnxr_id        TEXT PRIMARY KEY,  -- "MNXR101234"; col 1 of reac_prop.tsv
  mnx_equation   TEXT,              -- col 2, e.g. "1 MNXM01@MNXD1 = 1 MNXM1@MNXD1"
  reference      TEXT,              -- col 3, source provenance
  classifs       TEXT,              -- col 4, semicolon-separated classification IDs (mostly EC numbers)
  is_balanced    TEXT,              -- col 5, "B" if balanced, NULL/empty otherwise
  is_transport   TEXT               -- col 6, "T" if transport, NULL/empty otherwise
)
```

**Reaction direction is not tracked.** Empirical check on MNXref 4.5: all 83,796 entries in `reac_prop.tsv` use the ` = ` separator (reversible). MNX does not encode forward/reverse-only reactions in this dump. Phase 1.1B's parser splits the equation on ` = ` and treats every reaction as reversible — which matches Spec 1.2's permissive direction-handling rollup default.

```
reaction_aliases (
  source      TEXT,                -- "kegg.reaction" | "rhea" | "bigg.reaction" | ...
  value       TEXT,                -- e.g. "R00200"
  mnxr_id     TEXT REFERENCES reactions(mnxr_id)
)
CREATE INDEX idx_reaction_aliases_value ON reaction_aliases(value);
CREATE INDEX idx_reaction_aliases_source_value ON reaction_aliases(source, value);
```

#### Source vocabulary normalization (open caveat)

MNX's xref files use **dual short/long prefix forms** for the same source — `chebi:` and `CHEBI:`, `kegg.compound:` and `keggC:`, `bigg.metabolite:` and `biggM:`, etc. About 30 distinct prefixes total. Phase 1.1B's resolver builder normalizes to a single canonical form per source (preferring the long bioregistry-style spelling). The exact mapping table is part of the implementation and TBD until we observe how Phase 2 paper extraction queries the resolver. Documented in the diagnostic report.

Estimated row counts (sized against MNXref 4.5; refine after first build):
- compounds ~340K (chem_prop.tsv has ~1.5M lines including comments + ~~340K~~ TBD data rows pending parse)
- compound_aliases: large, dominated by SwissLipids (~1.5M) and HMDB / ChEBI / Reactome (each 100K-500K)
- compound_names: ~1M after extracting from the `name` and `description` columns
- reactions: ~80K
- reaction_aliases: ~300K, distributed across rhea / bigg.reaction / kegg.reaction / metacyc.reaction / etc.

### TCDB hierarchy cache (`cache/data/tcdb/tcdb_hierarchy.json`)

Built from a join of three TCDB sources downloaded in sub-step 6:

| TCDB file | What it contributes |
|---|---|
| `families.tsv` (130 KB, 2-col `<TC_ID>\t<description>`) | Family-level (3-component) descriptions — populates `name` field on family nodes. ~2K rows. |
| `superfamilies.tsv` (778 KB, 5-col `TCID Subfamily Family Fam_abbreviation Superfamily`) | Full 5-component specificity-level enumeration (11,465 rows) + every TCID's complete ancestor chain + superfamily group + family abbreviation. |
| `substrates.tsv` (436 KB, 2-col `<TC_ID>\t<pipe-separated CHEBI:N;name>`) | `substrate_classes` field at specificity level. TCIDs are 5-component. |

Shape:

```json
{
  "1":         {"name": "Channels and Pores",  "level": 0, "level_kind": "tc_class",        "parent": null},
  "1.A":       {"name": "α-Type Channels",     "level": 1, "level_kind": "tc_subclass",     "parent": "1"},
  "1.A.1":     {"name": "VIC Family",          "level": 2, "level_kind": "tc_family",       "parent": "1.A",       "abbreviation": "VIC"},
  "1.A.1.1":   {"name": "Kv subfamily",        "level": 3, "level_kind": "tc_subfamily",    "parent": "1.A.1"},
  "1.A.1.1.1": {"name": "Shaker channel",      "level": 4, "level_kind": "tc_specificity",  "parent": "1.A.1.1",   "substrate_classes": ["potassium ion"], "superfamily": "VIC Superfamily"}
}
```

- Class (level 0) and subclass (level 1) entries are **not present** in any TCDB bulk file — they're synthesized by the builder from the prefix of TCIDs in `superfamilies.tsv` (`1` and `1.A` derived from `1.A.1.1.1`). Their `name` values are hardcoded from TCDB's class/subclass naming convention (Channels/Pores, α-Type Channels, etc.) using a small lookup table maintained inline in `build_metabolite_resolver.py`.
- Family (level 2) `name` comes from `families.tsv`. `abbreviation` and `superfamily` come from `superfamilies.tsv`.
- Subfamily (level 3) entries are derived from `superfamilies.tsv` col 2 distinct values. `name` is empty unless TCDB later exposes subfamily names.
- Specificity (level 4) entries come from `superfamilies.tsv` col 1; `substrate_classes` joined from `substrates.tsv`; `superfamily` joined from `superfamilies.tsv` col 5.
- `superfamily` is treated as cross-cutting metadata (analogous to CAZy clans), not a hierarchy parent.

### CAZy hierarchy cache (`cache/data/cazy/cazy_hierarchy.json`)

**Bootstrapped from observed eggNOG `CAZy` columns** — no bulk download. Phase 1.1B's resolver builder iterates over all strain `gene_annotations_merged.json` files (after step 2), collects the union of CAZy IDs, and synthesizes the 3-level hierarchy mechanically:

```json
{
  "GH":      {"name": "Glycoside Hydrolases",         "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "GH"},
  "GT":      {"name": "GlycosylTransferases",         "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "GT"},
  "PL":      {"name": "Polysaccharide Lyases",        "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "PL"},
  "CE":      {"name": "Carbohydrate Esterases",       "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "CE"},
  "AA":      {"name": "Auxiliary Activities",         "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "AA"},
  "CBM":     {"name": "Carbohydrate-Binding Modules", "level": 0, "level_kind": "cazy_class",     "parent": null,    "class": "CBM"},
  "GH13":    {"name": "",                             "level": 1, "level_kind": "cazy_family",    "parent": "GH",    "class": "GH"},
  "GH13_1":  {"name": "",                             "level": 2, "level_kind": "cazy_subfamily", "parent": "GH13",  "class": "GH"}
}
```

Class names (level 0) are hardcoded from CAZy's class definitions (6 classes: GH, GT, PL, CE, AA, CBM). Family and subfamily names are intentionally empty — neither bulk download nor eggNOG provides per-family human descriptions. Phase 1.3 may add them via a separate scrape if needed; for the scaffold we get by without.

Hierarchy parents are derived from the family ID syntax:
- `GH13_1` → parent `GH13` (strip `_<digits>$`)
- `GH13` → parent `GH` (strip `<digits>+$`)
- `GH` → root (no parent)

The bootstrap means: **only families/subfamilies our 27 strains actually reference will appear in the hierarchy** (estimated ~80–150 entries vs. ~700 in CAZy total). This is intentional — the KG only ever queries CAZy IDs that exist on its genes.

### Diagnostic report (`metabolite_id_mapping_report.json`)

```json
{
  "mnx_release": "MNXref 4.5 (2025-08-13)",
  "compound_count": 343219,
  "reaction_count": 81204,
  "alias_counts_by_source": {"kegg.compound": 19373, "chebi": 427380, "hmdb": 298876, "rhea": 69008, ...},
  "source_normalization_summary": {
      "chebi": "merged: chebi (223004) + CHEBI (204376)",
      "kegg.compound": "merged: kegg.compound (19373) + keggC (38746)",
      "...": "..."
  },
  "tcdb_class_count": 9,
  "tcdb_specificity_count": 11465,
  "cazy_family_count_observed_in_eggnog": 87
}
```

## Accessor module layout (`utils/`)

Three flat files in `multiomics_kg/utils/`, matching the existing `pfam_utils.py` / `brite_utils.py` convention:

| File | API surface | Used by |
|---|---|---|
| `utils/metabolite_utils.py` | `open_resolver()`, `resolve_metabolite()`, `resolve_reaction()`, path helpers | step 2 (gene annotations), Spec 1.2 scaffold builder, Phase 2 extraction, tests; Spec 1.2 adds `load_scaffold()` and `COFACTOR_MNXM_IDS` |
| `utils/tcdb_utils.py` | `load_tcdb()` (cached at module level), `tcdb_ancestors(tc_id)`, `is_valid_tcdb(tc_id)` | step 2 (validate), Spec 1.2 scaffold builder (ancestor closure), Spec 1.3 adapter |
| `utils/cazy_utils.py` | `load_cazy()` (cached), `cazy_ancestors(cazy_id)`, `is_valid_cazy(cazy_id)` | step 2 (validate), Spec 1.2 scaffold builder (ancestor closure), Spec 1.3 adapter; Spec 1.3 adds `CAZY_SUBSTRATE_HINTS` |

TCDB and CAZy live in separate files because they're conceptually unrelated reference systems (transporter classification vs. carbohydrate-active enzymes), even though the access API shape is identical.

**Boundary rules:**

- Adapters never import from `download/` — the contract surface is the cache files (`gene_annotations_merged.json`, `scaffold_pruned.json`).
- `download/build_gene_annotations.py` (step 2) imports the accessors (`metabolite_utils.resolve_reaction`, `tcdb_utils.is_valid_tcdb`, `cazy_utils.is_valid_cazy`) — it's a consumer of the resolver, not a builder of it.
- `download/build_metabolism_scaffold.py` (Spec 1.2's new top-level step) similarly consumes accessors plus its own raw-MNX-TSV loader for substrate/product extraction.
- Hierarchy JSON loaders (`load_tcdb()`, `load_cazy()`) cache in module-level globals — consumers don't deal with file I/O.

### Resolver Python API (`utils/metabolite_utils.py`)

```python
def open_resolver(path: Path | None = None) -> sqlite3.Connection: ...

def resolve_metabolite(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """
    Returns (mnxm_id, method).
    method ∈ {
        "xref:exact",          # unique alias hit
        "xref:ambiguous",      # alias maps to >1 mnxm; first returned
        "name:normalized",     # name match after normalization
        "ambiguous",           # name match maps to >1 mnxm
        "unresolved",          # nothing matched
    }
    """

def resolve_reaction(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """Same shape as resolve_metabolite, against the reactions tables."""
```

Lookup order in `resolve_metabolite`:

1. If `value` looks like an MNXM directly (`MNXM` prefix), return it after existence check.
2. `compound_aliases` exact match by value (across all sources). If unique, return `("MNXM...", "xref:exact")`. If multiple, return the first by `mnxm_id` ordering with `"xref:ambiguous"`.
3. `compound_names` after name normalization (lowercase, collapse whitespace, strip punctuation outside of `+/-`). Single hit → `"name:normalized"`; multi → `"ambiguous"`.
4. Else `"unresolved"`.

`resolve_reaction` runs the same alias-matching logic against `reactions` (direct MNXR-prefix check) and `reaction_aliases`, but **omits the name-normalization step** — reactions don't have a `reaction_names` table because the consumers (eggNOG `KEGG_Reaction` column, paper extraction) always cite reactions by structured ID, never by free-text name. Possible methods for `resolve_reaction`: `"xref:exact"`, `"xref:ambiguous"`, `"unresolved"`.

The same resolver is reused by Phase 2 paper-measurement extraction.

### Hierarchy accessors (`utils/tcdb_utils.py`, `utils/cazy_utils.py`)

```python
# tcdb_utils.py
def load_tcdb() -> dict[str, dict]: ...           # cached at module level
def tcdb_ancestors(tc_id: str) -> list[str]: ...  # ["1", "1.A", "1.A.1", "1.A.1.1"] for "1.A.1.1.1"
def is_valid_tcdb(tc_id: str) -> bool: ...

# cazy_utils.py — same shape
def load_cazy() -> dict[str, dict]: ...
def cazy_ancestors(cazy_id: str) -> list[str]: ...
def is_valid_cazy(cazy_id: str) -> bool: ...
```

## Step 2 — MNX-aware gene annotation merge (extended)

`multiomics_kg/download/build_gene_annotations.py` (existing 956-LOC file). Step number, script path, output file (`gene_annotations_merged.json`) all unchanged. **No code added to `build_gene_annotations.py` itself** — the work is done through the existing config-driven transforms framework.

### Existing config already extracts the raw eggNOG columns

`config/gene_annotations_config.yaml` already declares these fields (raw, no resolution / validation):

```yaml
kegg_reaction:                # currently raw R-numbers
  type: union
  sources: [{source: eggnog, field: KEGG_Reaction, delimiter: ","}]

transporter_classification:   # currently raw TC IDs (some may be invalid)
  type: union
  sources: [{source: eggnog, field: KEGG_TC, delimiter: ","}]

cazy_ids:                     # currently raw CAZy IDs (some may be invalid)
  type: union
  sources: [{source: eggnog, field: CAZy, delimiter: ","}]
```

`KEGG_ko`, `EC`, `BRITE`, `PFAMs`, `eggNOG_OGs`, `COG_category`, `GOs`, `Description`, `Preferred_name` continue to be ingested as before. `BiGG_Reaction`, `KEGG_Pathway`, `KEGG_Module`, `KEGG_rclass` remain unused — see overview "Out of scope".

### Three new per-token transforms

Added to `multiomics_kg/download/utils/annotation_transforms.py`, matching the existing `_tx_*(value: str) -> str | list[str] | None` signature. Returning `None` drops the token — the same mechanism used by `_tx_normalize_ec`. Internally call into the accessor modules from `utils/`.

```python
# annotation_transforms.py

from multiomics_kg.utils.metabolite_utils import open_resolver, resolve_reaction
from multiomics_kg.utils.tcdb_utils import is_valid_tcdb
from multiomics_kg.utils.cazy_utils import is_valid_cazy

_RESOLVER_CONN = None

def _get_resolver_conn():
    global _RESOLVER_CONN
    if _RESOLVER_CONN is None:
        _RESOLVER_CONN = open_resolver()
    return _RESOLVER_CONN

def _tx_resolve_kegg_reaction_to_mnxr(value: str) -> str | None:
    """KEGG R-number → MNXR ID. Returns None if unresolved (drops the token)."""
    mnxr, _method = resolve_reaction(value, _get_resolver_conn())
    return mnxr

def _tx_validate_tcdb(value: str) -> str | None:
    """Drop TC IDs not present in tcdb_hierarchy.json."""
    return value if is_valid_tcdb(value) else None

def _tx_validate_cazy(value: str) -> str | None:
    """Drop CAZy IDs not present in cazy_hierarchy.json."""
    return value if is_valid_cazy(value) else None
```

Module-level cached resolver connection (opened lazily on first call) — same pattern as `load_tcdb()` and `load_cazy()` caching in their respective accessor modules. `build_gene_annotations.py` is a one-shot CLI build; the connection lifecycle aligns with process lifetime.

### Config-driven wiring

The three existing field configs gain a `transforms` line (and the KEGG-reaction field renames to plural to reflect that values are now MNXR IDs, not raw R-numbers):

```yaml
kegg_reactions:                                # renamed from kegg_reaction (plural; values are MNXRs)
  type: union
  transforms: [resolve_kegg_reaction_to_mnxr]
  sources: [{source: eggnog, field: KEGG_Reaction, delimiter: ","}]

transporter_classification:                    # name unchanged (consumed by cyanorak_ncbi_adapter)
  type: union
  transforms: [validate_tcdb]
  sources: [{source: eggnog, field: KEGG_TC, delimiter: ","}]

cazy_ids:                                      # name unchanged (consumed by cyanorak_ncbi_adapter)
  type: union
  transforms: [validate_cazy]
  sources: [{source: eggnog, field: CAZy, delimiter: ","}]
```

The `union` type already runs each transform per-token and drops `None` returns; no framework change needed.

### Field-naming notes

- `kegg_reactions` (plural) is a new field name. The previous singular `kegg_reaction` is removed — no current consumer uses it (verified by grep).
- `transporter_classification` and `cazy_ids` keep their existing names because `cyanorak_ncbi_adapter.py` already reads them as Gene properties. Adding validation transforms only narrows the value set (drops invalid IDs); existing consumers remain compatible.

### What the transforms guarantee

Empty lists are still written for genes with no eggNOG hit (existing behaviour of the framework). Per-gene values after transforms:

- `kegg_reactions: list[MNXR]` — MNXR IDs only; unresolved R-numbers dropped.
- `transporter_classification: list[str]` — only TC IDs that exist in `tcdb_hierarchy.json`.
- `cazy_ids: list[str]` — only CAZy IDs that exist in `cazy_hierarchy.json`.

### Per-strain metabolism report

Written to `cache/data/<organism>/genomes/<strain>/step2_metabolism_report.json`:

```json
{
  "strain": "MED4",
  "gene_count": 1990,
  "kegg_reactions": {
    "raw_total":              1532,
    "raw_unique":              892,
    "resolved_total":         1488,
    "resolved_unique_mnxr":    843,
    "unresolved_unique":        49,
    "unresolved_examples":   ["R12345", "R87654"]
  },
  "transporter_classification": {
    "raw_total":           420,
    "raw_unique":           97,
    "validated_total":     412,
    "validated_unique":     94,
    "invalid_examples":   ["1.X.99"]
  },
  "cazy_ids": {
    "raw_total":          88,
    "raw_unique":         44,
    "validated_total":    88,
    "validated_unique":   44,
    "invalid_examples":   []
  }
}
```

The transforms framework itself drops invalid tokens silently. The report is generated by a **separate** post-merge sweep over the eggNOG TSV (counting raw values) and the merged JSON (counting validated/resolved values). This sweep lives next to `build_gene_annotations.py` (e.g. invoked at the end of the build) — same module-level cached resolver connection.

Aggregated cross-strain summary written to `logs/prepare_data_step2_metabolism_summary.txt`.

## Tests

### Unit tests (`tests/test_build_metabolite_resolver.py`)

- Load a small synthetic MNX TSV bundle into a temp resolver DB; verify table row counts.
- `resolve_metabolite("C00031", conn)` → `("MNXM41", "xref:exact")` (KEGG glucose example).
- `resolve_metabolite("D-glucose", conn)` → name-based hit.
- `resolve_metabolite("unicorn-mol", conn)` → `(None, "unresolved")`.
- Ambiguous alias case → `"xref:ambiguous"`.
- TCDB hierarchy: round-trip parse → JSON → re-parse, count check.
- CAZy hierarchy: same.

### Unit tests (`tests/test_build_gene_annotations_metabolism.py`)

- Build a tiny eggNOG annotations file with known KEGG_Reaction / KEGG_TC / CAZy values.
- Run step-2 merge against a small synthetic resolver.
- Assert that the resulting `gene_annotations_merged.json` contains the expected `kegg_reactions`, `transporter_classification`, `cazy_ids` lists.
- Assert empty-list defaults for genes with no eggNOG hit.
- Invalid TC / CAZy values are dropped, not retained.

### Integration test (`tests/test_prepare_data_step2_metabolism_smoke.py`)

- Marker `@pytest.mark.slow`.
- Assumes NCBI / UniProt / `gene_mapping.csv` for MED4 are already cached locally (the fixture from prior runs). Avoids the cost of re-running step 0 sub-steps 1–5.
- Runs only the new metabolism sub-steps + the gene-annotation merge:
  ```bash
  uv run python -m multiomics_kg.download.download_genome_data --steps 6 7 --strains MED4 --force
  uv run python -m multiomics_kg.download.build_gene_annotations --strains MED4 --force
  ```
- Assert: report file exists, `kegg_reactions.resolved_unique_mnxr ≥ 400` (MED4 has 526 KEGG_Reaction-annotated genes; resolution rate is the unknown), `transporter_classification.validated_unique ≥ 50` (MED4 has 121 KEGG_TC-annotated genes).

## Acceptance criteria

- `cache/data/mnx/metabolite_resolver.db` exists with non-empty `compounds`, `compound_aliases`, `compound_names`, `reactions`, `reaction_aliases` tables. `compound_count` field in the diagnostic report is ≥ 100K.
- `cache/data/tcdb/tcdb_hierarchy.json` contains ≥10K total entries (~11.5K specificity-level + ancestors) and a connected hierarchy (every non-root entry has a parent that exists in the same JSON).
- `cache/data/cazy/cazy_hierarchy.json` contains all six top-level classes and at least 50 family entries (bootstrapped from observed eggNOG `CAZy` columns; exact count depends on what eggNOG annotated).
- The seven sub-step-6 source files are all present at the paths in the SOURCES table; sizes within ±50% of the values listed (MNX is ~1.5 GB; TCDB tables are < 1 MB each).
- For every strain in `cyanobacteria_genomes.csv`, `gene_annotations_merged.json` contains the three new / refined fields (`kegg_reactions`, `transporter_classification`, `cazy_ids`) on every gene; no field is ever `null` (always at minimum an empty list).
- For MED4 specifically: `kegg_reactions.resolved_unique_mnxr ≥ 400` (lower bound; depends on KEGG R-number resolution rate), `transporter_classification.validated_unique ≥ 50`, `cazy_ids.validated_unique ≥ 5` (low because CAZy coverage on MED4 is 1.3%).
- `pytest -m "not slow and not kg"` passes.
- `bash scripts/prepare_data.sh` end-to-end succeeds with no top-level step renumbering. Existing `--steps N M` invocations still work.
- CLAUDE.md updated with the new step 0 sub-steps and the step-2 internal extension.
- **No KG validity tests change**, no schema_config.yaml change, no adapter change. The deployed graph is byte-identical to before this spec.

## Out of scope (handled in later sub-specs)

- Pruned scaffold cache (new top-level step in Spec 1.2).
- New schema nodes / edges → Specs 1.2 and 1.3.
- New adapters → Specs 1.2 and 1.3.
- Post-import Cypher → Specs 1.2 and 1.3.
- KG validity tests for new content → Specs 1.2 and 1.3.
