# Phase 1.2 — Reactions scaffold (revised, KEGG-native)

**Date:** 2026-04-30
**Status:** Ready for implementation plan
**Supersedes:** [2026-04-28-metabolite-reactions-scaffold-design.md](2026-04-28-metabolite-reactions-scaffold-design.md)
**Parent overview:** [2026-04-28-metabolite-scaffold-design.md](2026-04-28-metabolite-scaffold-design.md)
**Depends on:** Phase 1.1B — Resolver + accessors (DONE on branch `metabolite-scaffold`, not yet merged)

## Goal

Add the chemistry layer to the KG: `Reaction` and `Metabolite` nodes plus the edges that connect them to genes, pathways, and organisms. The deliverable enables two query shapes:

1. **Nutrient-exchange reasoning** — "MED4 differentially expresses genes X, Y, Z under N-limitation; what metabolite repertoire / nutrient classes does that touch?"
2. **Paper integration substrate** — Phase 2 metabolomics papers can MERGE onto the canonical Metabolite nodes via the MNX resolver.

**Explicit non-goals:**
- Flux balance / GSMM analysis. We do not build or run stoichiometric models.
- Reaction direction. Substrate/product is not modelled.
- Subcellular compartments. MNX `@MNXD<n>` markers are stripped.
- Stoichiometric coefficients. Compounds are tracked as sets.

## Schema additions

### Nodes (2 new)

**`Reaction`**

| Property | Type | Source | Notes |
|---|---|---|---|
| `id` | str | KEGG-native | `kegg.reaction:R00200` |
| `kegg_reaction_id` | str | KEGG `/list/reaction` | `R00200` (denormalized for direct query) |
| `name` | str | KEGG `/list/reaction` | "pyruvate kinase reaction" |
| `ec_numbers` | str[] | MNX `reactions.classifs` (or KEGG `/get/rn:`) | `["2.7.1.40"]` |
| `kegg_pathway_ids` | str[] | KEGG `/link/reaction/pathway` | `["ko00010", "ko00710"]` (also via edge) |
| `mnxr_id` | str \| null | MNX resolver (xref enrichment) | `"MNXR101234"` or null when not in MNX |
| `rhea_ids` | str[] | MNX `reaction_aliases` source=rhea | `["10828"]` |
| `reaction_class` | str | MNX `reactions.is_transport` | `"transport"` (T flag set) \| `"chemical"` (default) |
| `mass_balance` | str | MNX `reactions.is_balanced` | `"balanced"` (B flag set) \| `"unbalanced"`. Informational only. |
| computed: `organisms` | str[] | post-import | distinct organism names with a gene catalyzing |
| computed: `organism_count` | int | post-import |  |
| computed: `gene_count` | int | post-import | distinct genes catalyzing |

**`Metabolite`**

| Property | Type | Source | Notes |
|---|---|---|---|
| `id` | str | KEGG-native (with fallback) | `kegg.compound:C00031`; falls back to `chebi:<id>` then `mnx:<id>` for compounds without KEGG entry |
| `kegg_compound_id` | str \| null | KEGG `/list/compound` | `C00031` or null |
| `name` | str | KEGG `/list/compound` | "D-glucose" |
| `formula` | str | MNX `compounds.formula` | `"C6H12O6"` |
| `mass` | float | MNX `compounds.mass` | `180.063` |
| `inchikey` | str | MNX `compounds.inchikey` | `"WQZGKKKJIJFFOK-..."` |
| `smiles` | str | MNX `compounds.smiles` |  |
| `mnxm_id` | str \| null | MNX resolver (xref enrichment) | `"MNXM41"` or null |
| `chebi_id` | str \| null | MNX `compound_aliases` source=chebi | `"17234"` |
| `hmdb_id` | str \| null | MNX `compound_aliases` source=hmdb | `"HMDB0000122"` |
| computed: `organism_count` | int | post-import |  |
| computed: `gene_count` | int | post-import |  |

### Edges (4 new)

| Edge | Direction | Source | Notes |
|---|---|---|---|
| `Gene_catalyzes_reaction` | Gene → Reaction | `gene["kegg_reactions"]` (per strain) | One per gene-R-number pair |
| `Reaction_has_metabolite` | Reaction → Metabolite | KEGG `/link/compound/reaction` | **Direction-agnostic.** No substrate/product split. |
| `Reaction_in_kegg_pathway` | Reaction → KeggTerm | KEGG `/link/reaction/pathway` | Reuses existing KeggTerm pathway nodes (372 of them) |
| `Organism_has_metabolite` | Organism → Metabolite | post-import: derived from Gene→Reaction→Metabolite | Materialized 2-hop save; multi-evidence-source ready (Phase 2 will add paper-derived evidence) |

**Not added:** `Reaction_in_organism` (use computed `Reaction.organisms[]` property instead — 1-hop derivable, matches existing OrthologGroup convention).

## Pipeline changes

### `prepare_data.sh` — new step 6

Slot a new top-level step analogous to step 5 (`build_og_descriptions.py` — prune-then-enrich pattern):

```
0  download_genome_data           (incl. MNX/TCDB/CAZy from Phase 1.1A/B)
1  build_protein_annotations
2  build_gene_annotations          ← validates kegg_reactions / TC / CAZy
3  build_gene_id_mapping
4  resolve_paper_ids
5  build_og_descriptions           ← existing prune-then-enrich
6  build_kegg_metabolism_xrefs     ← NEW: same pattern, scoped to gene-reachable KEGG IDs
```

**Module:** `multiomics_kg/download/build_kegg_metabolism_xrefs.py`

**CLI:** `uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]`

**Pruning (the load-bearing part — like step 5):**
1. Walk every strain's `gene_annotations_merged.json`. Collect union of `gene["kegg_reactions"]` (raw KEGG R-numbers).
2. Look up each R-number's compounds via the (cached) KEGG `/link/compound/reaction` data. Collect union of C-numbers.
3. Result: a "needed" set of ~3-4K R-numbers and ~1-2K C-numbers (out of KEGG's ~14K and ~22K). Most central metabolism is shared across strains.

**Enrichment (only for the needed set):**
- For each R-number: KEGG name, EC, pathway IDs (from `kegg_data.json`); MNXR + Rhea + is_balanced + is_transport (from MNX resolver).
- For each C-number: KEGG name (from `kegg_data.json`); MNXM + ChEBI + HMDB + InChIKey + formula + mass + SMILES (from MNX resolver).

**Output:** `cache/data/kegg/kegg_metabolism_xrefs.json` (~1 MB, only the needed subset)

```json
{
  "reactions": {
    "R00200": {
      "name": "pyruvate kinase reaction",
      "ec_numbers": ["2.7.1.40"],
      "kegg_pathway_ids": ["ko00010", "ko00710"],
      "compound_ids": ["C00074", "C00008", "C00022", "C00002"],
      "mnxr_id": "MNXR101234",
      "rhea_ids": ["10828"],
      "mass_balance": "balanced",
      "reaction_class": "chemical"
    }
  },
  "compounds": {
    "C00031": {
      "name": "D-glucose",
      "formula": "C6H12O6",
      "mass": 180.063,
      "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
      "smiles": "OC[C@H]1OC(O)...",
      "mnxm_id": "MNXM41",
      "chebi_id": "17234",
      "hmdb_id": "HMDB0000122"
    }
  }
}
```

### `kegg_utils.py` extension

Extend [multiomics_kg/utils/kegg_utils.py](multiomics_kg/utils/kegg_utils.py) — currently downloads 4 endpoints (KO list, KO→pathway link, BRITE hierarchy, pathway list) into `kegg_data.json`. Add 5 endpoints:

| Endpoint | Cache key |
|---|---|
| `/list/reaction` | `reaction_names` (~14K R# → name) |
| `/list/compound` | `compound_names` (~22K C# → name) |
| `/link/compound/reaction` | `reaction_to_compounds` (R# → [C#, ...]) |
| `/link/pathway/reaction` | `reaction_to_pathways` (R# → [ko#, ...]) |
| `/link/pathway/compound` | `compound_to_pathways` (C# → [ko#, ...]) — optional, useful for compound enrichment lens |

Same `_fetch_text` / parser pattern as existing 4. Skip-if-cached.

### YAML config edit

In [config/gene_annotations_config.yaml](config/gene_annotations_config.yaml), drop the `transform: resolve_kegg_reaction_to_mnxr` line on the `kegg_reactions` field. Resulting per-gene field becomes raw KEGG R-numbers (recovers the 6.5% currently lost through MNX). The existing `validate_tcdb` / `validate_cazy` transforms stay (they read tiny static JSON, not SQLite).

### Phase 1.1B follow-ups absorbed into Spec 1.2

Three deferred tasks land here:

1. **`mnxm_to_primary_id()` and `mnxr_to_primary_id()` helpers** in `metabolite_utils.py`. Used by the new `build_kegg_metabolism_xrefs.py` to set primary IDs deterministically. KEGG > ChEBI > MNX fallback for compounds; KEGG > Rhea > MNX for reactions.
2. **Indexes on resolver DB:** `idx_compound_aliases_mnxm` and `idx_reaction_aliases_mnxr`. Needed because the new helpers query by the MNX-side ID. Add to the table DDL in `build_metabolite_resolver.py`.
3. **`source_normalization_summary` in diagnostic report.** Surface unmapped MNX prefixes (`sabiork.compound`, `vmhmetabolite`, `rheaG`, etc.). Lower priority.

### `schema_config.yaml`

Add the 2 nodes and 4 edges. Pattern matches existing `OrthologGroup` / `Pfam` entries (canonical ID, full property list, BiocCypher MERGE semantics).

### Adapter (`metabolism_adapter.py`)

New module under `multiomics_kg/adapters/`. Constructor reads:

- `cache/data/kegg/kegg_metabolism_xrefs.json` (the pruned + enriched cache)
- Each strain's `gene_annotations_merged.json` (for `gene["kegg_reactions"]`)
- `genome_config_file` (to walk strains)

`get_nodes()`:
- For each R-num in `xrefs["reactions"]`: emit `Reaction` with the enriched property dict
- For each C-num in `xrefs["compounds"]`: emit `Metabolite` with the enriched property dict

`get_edges()`:
- For each (gene, R-num) pair: emit `Gene_catalyzes_reaction`
- For each R-num: emit `Reaction_has_metabolite` (one per `compound_ids` entry)
- For each R-num: emit `Reaction_in_kegg_pathway` (one per `kegg_pathway_ids` entry)

**No SQLite at adapter time. No KEGG REST at adapter time. Pure file reads.**

The adapter follows the Multi-X pattern: a `MultiMetabolismAdapter` that wraps per-strain or one flat list (TBD which is cleaner; the per-gene edges are the only strain-scoped output, so a single flat adapter may be simpler).

### Post-import Cypher

Add to `scripts/post-import.sh` (and reference `post-import.cypher`):

```cypher
// Reaction.gene_count, organism_count, organisms[]
MATCH (r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
WITH r, count(DISTINCT g) AS gene_count, collect(DISTINCT g.organism_name) AS organisms
SET r.gene_count = gene_count,
    r.organism_count = size(organisms),
    r.organisms = organisms;

// Metabolite.gene_count, organism_count
MATCH (m:Metabolite)<-[:Reaction_has_metabolite]-(r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
WITH m, count(DISTINCT g) AS gene_count, count(DISTINCT g.organism_name) AS organism_count
SET m.gene_count = gene_count, m.organism_count = organism_count;

// Organism_has_metabolite edges (materialized 2-hop)
MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
      -[:Gene_catalyzes_reaction]->(r:Reaction)
      -[:Reaction_has_metabolite]->(m:Metabolite)
WITH DISTINCT o, m
MERGE (o)-[:Organism_has_metabolite]->(m);

// Organism rollup props
MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
WITH o, count(DISTINCT m) AS metabolite_count
SET o.metabolite_count = metabolite_count;

MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)-[:Gene_catalyzes_reaction]->(r:Reaction)
WITH o, count(DISTINCT r) AS reaction_count
SET o.reaction_count = reaction_count;
```

### Indexes (post-import)

```cypher
CREATE INDEX reaction_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.id);
CREATE INDEX reaction_kegg_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.kegg_reaction_id);
CREATE INDEX reaction_mnxr_idx IF NOT EXISTS FOR (r:Reaction) ON (r.mnxr_id);
CREATE INDEX metabolite_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.id);
CREATE INDEX metabolite_kegg_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.kegg_compound_id);
CREATE INDEX metabolite_mnxm_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.mnxm_id);
CREATE INDEX metabolite_chebi_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.chebi_id);
CREATE FULLTEXT INDEX reactionFullText IF NOT EXISTS FOR (r:Reaction) ON EACH [r.name];
CREATE FULLTEXT INDEX metaboliteFullText IF NOT EXISTS FOR (m:Metabolite) ON EACH [m.name];
```

## Tests

### Unit tests
- `tests/test_build_kegg_metabolism_xrefs.py` — synthetic gene_annotations + synthetic kegg_data + tiny resolver DB; assert pruning rule (only-gene-reachable) and enrichment shape.
- `tests/test_metabolism_adapter.py` — synthetic xrefs JSON + synthetic gene annotations; assert correct node/edge yield.

### KG validity tests
- `tests/kg_validity/test_metabolism.py`:
  - Reaction count between 2K-5K (gene-reachable subset)
  - Metabolite count between 1K-3K
  - Every Reaction has at least 1 `Reaction_has_metabolite` edge (or document the exceptions)
  - Every Reaction's `kegg_reaction_id` matches the primary ID's suffix
  - Every Metabolite's `kegg_compound_id` matches when present
  - `Organism_has_metabolite` edges: > 1K total
  - Spot check: `kegg.reaction:R00200` (pyruvate kinase) exists, has 4 metabolites, in glycolysis pathway

### Integration smoke test
- `tests/test_metabolism_smoke.py` (slow): full prepare-data step 6 + build → assert MED4 has expected reaction/metabolite count rollups

## Acceptance criteria

- `cache/data/kegg/kegg_metabolism_xrefs.json` exists, ≤2 MB, contains the gene-reachable subset
- `Reaction` nodes ≥ 2,000 with primary IDs all `kegg.reaction:R*`
- `Metabolite` nodes ≥ 1,000 with primary IDs predominantly `kegg.compound:C*`
- `Gene_catalyzes_reaction` edges: at least one per gene with non-empty `kegg_reactions`
- `Reaction_has_metabolite` edges: ~4 per reaction average
- `Reaction_in_kegg_pathway` edges: present and reusing existing `KeggTerm` pathway nodes
- `Organism_has_metabolite` edges: materialized post-import for all 25 genome organisms
- `pytest -m "not slow and not kg"` passes
- `pytest -m kg` passes (KG-validity tests)
- Build doesn't open `metabolite_resolver.db` (confirm: `lsof` or grep adapter source)
- Docker build container size doesn't grow significantly (the resolver DB is not bundled)

## Out of scope (Spec 1.3 / Phase 2)

- TCDB transporter nodes + edges → Spec 1.3
- CAZy family nodes + edges → Spec 1.3
- Paper-derived metabolite measurements → Phase 2 (uses MNX resolver via `_resolved.csv` per-paper)
- OrthologGroup-mediated reaction inheritance → Phase 3
- Compound names from `/get/cpd:` per-compound (richer chemistry metadata) → out of scope unless needed
- Reaction direction inference from pathway maps (KGML) → out of scope
