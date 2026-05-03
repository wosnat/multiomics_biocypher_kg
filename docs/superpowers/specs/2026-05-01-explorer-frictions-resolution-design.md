# Explorer-side frictions — KG-side resolution design

**Date:** 2026-05-01
**Status:** Ready for implementation plan
**Conversation seed:** [`multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-frictions-reframed.md`](../../../../multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-frictions-reframed.md) — friction descriptions, resolution-space, joint-resolutions, explorer review, KG response

## Goal

Resolve four explorer-side frictions surfaced in the 2026-04-29 trigger analysis (`axenic_up_hypotheticals_med4`):

- **F1** — term-level content informativeness is invisible (Gene-level rollup conflates catch-all terms with informative ones)
- **F2** — empty per-gene tool results conflate "no hit" with "out of scope"
- **F3** — expression-bin clustering analyses look identical to functional-cluster analyses
- **F4** — floor-case (and broader sparsely-annotated) genes have no characterization surface

The four resolutions ship as a coordinated KG release. F1 and F2 introduce new schema; F3 is a vocabulary rename; F4 adds three Gene properties (one of which requires upstream NCBI-download plumbing).

## Out of scope

- **Ortholog-mediated annotation transfer** (F4 option c). Defer to its own spec when the broader sparsely-annotated population becomes the bottleneck. Wrong-scoped under F4 because its leverage is on a different population than F4's strict floor cases.
- **Per-pipeline run version / date metadata** (F2 anti-scope).
- **Per-experiment informativeness** of ontology terms (F1 anti-scope).
- **InterPro / AlphaFold layers** (F4 anti-scope, confirms existing deferral).
- **`KGRelease` versioning node**. The per-`DataSource.version` field carries enough release context for v1.
- **`data_floor: bool` Gene flag**. Derivable from F1's `annotation_state` + F2's `contributing_sources`; no parallel field.
- **Re-curation of NATL2A diel ClusteringAnalysis** (F3 anti-scope). Different problem.
- **MCP tool `gene_neighbors`** and its window-default. Lives on the explorer/MCP roadmap (Pass B); KG side just plumbs the contig field.

---

## F1 — Term informativeness + annotation rollup refinement

Four ships, one release:

### Ship 1.1 — `is_uninformative` flag on ontology term nodes

Add `is_uninformative: str` (sentinel value `"true"` when flagged; **absent** otherwise) to nodes of the following types — vocabulary derived from a 2026-05-02 KG audit (live counts):

| Node type | Flagged count | Mechanism |
|---|---|---|
| `BiologicalProcess` | 1 | ID list (GO root only) |
| `MolecularFunction` | 1 | ID list (GO root only) |
| `CellularComponent` | 1 | ID list (GO root only) |
| `CogFunctionalCategory` | 1 | ID list (S "Function unknown" only; **R "General function prediction only" stays un-flagged** — borderline informative, same reasoning as DUF/UPF) |
| `CyanorakRole` | 5 | ID list (`R` "Other" parent + R.1, R.2, R.4, R.5 hypothetical-class sub-roles) |
| `TigrRole` | 5 | ID list (Hypothetical proteins, Not Found, Unclassified, Unknown function/General) |
| `KeggTerm` | ~210 | **Name pattern**: `^K\d+;\s+uncharacterized protein\b` — too many KO IDs to hand-list, but the name pattern is unambiguous |

**Guiding principle for the catch-all vocabulary:** if the term tells the consumer the broad **class** (transporter, enzyme, peptidase, two-component system) — even if a sub-class dimension is unknown — treat the term as **informative**. Only flag terms that convey no class signal at all ("Hypothetical protein", "Not Found", "Function unknown", GO roots).

Excluded entirely:
- **`EcNumber`** — every EC number is functionally specific by design; the level=0 entries (Hydrolases, Isomerases, etc.) are broad but informative.
- **`BriteCategory`** — every BRITE entry sits inside a class-bearing tree (Transporters, Peptidases, Two-component, etc.); "Other transporters", "Peptidases of unknown catalytic type", and similar all still convey the broad class. No catch-alls remain.
- **`Pfam` / `PfamClan`** — DUF/UPF and PfamClans are structural signals, treated as informative per the F1 anti-scope.

Examples of terms intentionally *kept* informative under this principle (would have been candidates under naive name matching, but their broad class is clearly stated): `cyanorak.role:Q.9` "Transport and binding proteins > Unknown substrate"; `cyanorak.role:R.3` "Other > Enzymes of unknown specificity"; `tigr.role:141` "Transport and binding proteins / Unknown substrate"; `tigr.role:703` "Unknown function / Enzymes of unknown specificity"; all BRITE "Other X" sub-categories.

**Type rationale.** Project convention is `str` `"true"`/absence for boolean-shaped properties — see `rankable`, `has_p_value`, `significant` on DerivedMetric. BioCypher has a known bug serializing `bool` properties, so even though `is_uninformative` is set entirely in post-import Cypher (never through the BioCypher CSV writer), we use the str-sentinel pattern for type-contract consistency with consumers and to avoid schema validation surprises.

`is_uninformative` is set to `"true"` only — informative terms have the property absent (`IS NULL`), not `"false"`. This keeps Cypher checks one-sided: `t.is_uninformative IS NULL` means informative; `t.is_uninformative = 'true'` means flagged.

Vocabulary source: hardcoded YAML at [`config/uninformative_terms.yaml`](../../../config/uninformative_terms.yaml) (new file). Each section is a node-type binding; supports both `ids:` (hand-listed CURIEs) and `name_patterns:` (regex against `t.name`):

```yaml
# Catch-all / no-content-bearing terms across ontology sources.
# Pfam DUF/UPF stay UN-flagged (structural domain = informative).
# COG R ("General function prediction only") stays UN-flagged (borderline).

biological_process:
  ids:
    - go:0008150

molecular_function:
  ids:
    - go:0003674

cellular_component:
  ids:
    - go:0005575

cog_category:
  ids:
    - cog.category:S        # Function unknown

cyanorak_role:
  ids:
    - cyanorak.role:R       # Other (parent — used standalone = no class info)
    - cyanorak.role:R.1     # Other > Conserved hypothetical domains
    - cyanorak.role:R.2     # Other > Conserved hypothetical proteins
    - cyanorak.role:R.4     # Other > Hypothetical proteins
    - cyanorak.role:R.5     # Other > Other
  # Q.9 (Transport > Unknown substrate) and R.3 (Other > Enzymes of unknown
  # specificity) intentionally NOT flagged — class (transporter, enzyme) is
  # known; only sub-dimension unknown.

tigr_role:
  ids:
    - tigr.role:156         # Hypothetical proteins / Conserved
    - tigr.role:704         # Hypothetical proteins / Domain
    - tigr.role:856         # Not Found
    - tigr.role:185         # Unclassified / Role category not yet assigned
    - tigr.role:157         # Unknown function / General
  # 141 (Transport / Unknown substrate) and 703 (Unknown function / Enzymes
  # of unknown specificity) intentionally NOT flagged — class is known.

kegg_term:
  # ~210 KOs match this pattern; ID list omitted (too many; name pattern stable)
  name_patterns:
    - '^K\d+;\s+uncharacterized protein\b'

# ec_number: section omitted — every EC is functionally specific.
# brite_category: section omitted — every BRITE entry sits inside a
#   class-bearing tree (Transporters, Peptidases, Two-component, etc.);
#   no entries pass the "no class signal" bar.
# pfam: section omitted — DUF/UPF stay informative-shaped per anti-scope.
# pfam_clan: section omitted — clans are structural superfamilies.
```

Set in post-import Cypher. For each section:
- For each `ids:` entry: `MATCH (t:<Label> {id: $id}) SET t.is_uninformative = 'true'`
- For each `name_patterns:` regex: `MATCH (t:<Label>) WHERE t.name =~ $pattern SET t.is_uninformative = 'true'`

The label binding is hardcoded per section name (`biological_process` → `BiologicalProcess`, etc.). Loader rejects unknown section names rather than silently dropping entries.

Vocabulary is small (~30 ID entries + 1 KEGG pattern matching ~210 nodes). Validation: a startup test loads the YAML and asserts each `ids:` entry resolves to an existing node post-build; pattern matches are reported as a count (not asserted to a specific number, since KEGG ref grows over time).

### Ship 1.2 — Redefine `Gene.annotation_quality` as numeric encoding of `annotation_state` (move to post-import)

Currently `annotation_quality` is computed at build-time in [`multiomics_kg/download/build_gene_annotations.py:209`](../../../multiomics_kg/download/build_gene_annotations.py) by counting presence of `gene["go_terms"]` / `kegg_ko` / `ec_numbers` / `pfam_ids` — catch-all-blind by construction, and folding in a heuristic `is_hypothetical_product_name` regex.

**Drop the product-naming heuristic. Make `annotation_quality` a 1:1 numeric encoding of `annotation_state`** (Ship 1.3). Both fields share a single derivation; consumers pick whichever shape (numeric for sort/threshold, categorical for display/exact-match) suits the query.

The shared input is **`informative_source_count`** — the number of distinct source buckets where the gene reaches at least one informative term:

| Bucket | Has-informative criterion | Status |
|---|---|---|
| `go` | `EXISTS { (g)-[:Gene_involved_in_biological_process\|Gene_enables_molecular_function\|Gene_located_in_cellular_component]->(t) WHERE t.is_uninformative IS NULL }` | live |
| `kegg` | `EXISTS { (g)-[:Gene_has_kegg_ko]->(t) WHERE t.is_uninformative IS NULL }` | live |
| `pfam` | `EXISTS { (g)-[:Gene_has_pfam]->() }` (no flag check; Pfam never flagged) | live |
| `ec` | `EXISTS { (g)-[:Gene_catalyzes_ec_number]->() }` (no flag check; EC never flagged) | live |
| `role` | `g.gene_category <> 'Unknown'` (coalesces cog_category / cyanorak_role / tigr_role; `gene_category` is already informativeness-aware by construction — see "Drift risk" below) | live |
| `reaction` | `EXISTS { (g)-[:Gene_catalyzes_reaction]->() }` (no flag check; every Reaction is KEGG-native and specific) | live (metabolism scaffold, ~21K genes have edges) |
| `transporter` | `EXISTS { (g)-[:Gene_has_tcdb_family]->() }` (no flag check; every TCDB family is functionally specific) | **live** as of 2026-05-02 (TCDB ontology landed in main; node label `TcdbFamily`) |
| `cazy` | `EXISTS { (g)-[:Gene_has_cazy_family]->() }` (no flag check; every CAZy family is functionally specific) | **live** as of 2026-05-02 (CAZy ontology landed in main; node label `CazyFamily`) |

GO BP/MF/CC coalesce into one bucket. The three role sources (cog_category, cyanorak_role, tigr_role) coalesce via `gene_category`. **All 8 buckets are now live** as of the 2026-05-02 main-merge — `transporter` and `cazy` were anticipated as forthcoming when this spec was first written (2026-05-01), and shipped the next day under the [`2026-05-01-tcdb-cazy-ontologies-design.md`](2026-05-01-tcdb-cazy-ontologies-design.md) spec.

Mapping:

| `informative_source_count` | edges exist? | `annotation_quality` | `annotation_state` |
|---|---|---|---|
| 0 | no | 0 | `no_evidence` |
| 0 | yes (all catch-all) | 1 | `catch_all_only` |
| 1 | — | 2 | `informative_single` |
| ≥ 2 | — | 3 | `informative_multi` |

Cypher (post-import, both fields set in one block):

```cypher
MATCH (g:Gene)
CALL {
  WITH g,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component]->(t) WHERE t.is_uninformative IS NULL } AS has_go,
       EXISTS { (g)-[:Gene_has_kegg_ko]->(t) WHERE t.is_uninformative IS NULL } AS has_kegg,
       EXISTS { (g)-[:Gene_has_pfam]->() } AS has_pfam,
       EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } AS has_ec,
       (g.gene_category IS NOT NULL AND g.gene_category <> 'Unknown') AS has_role,
       EXISTS { (g)-[:Gene_catalyzes_reaction]->() } AS has_reaction,
       EXISTS { (g)-[:Gene_has_tcdb_family]->() } AS has_transporter,
       EXISTS { (g)-[:Gene_has_cazy_family]->() } AS has_cazy,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component
                     |Gene_has_kegg_ko|Gene_has_pfam|Gene_catalyzes_ec_number
                     |Gene_in_cog_category|Gene_has_cyanorak_role|Gene_has_tigr_role
                     |Gene_catalyzes_reaction|Gene_has_tcdb_family|Gene_has_cazy_family]->() } AS has_any_edge
  WITH g,
       (CASE WHEN has_go THEN 1 ELSE 0 END
        + CASE WHEN has_kegg THEN 1 ELSE 0 END
        + CASE WHEN has_pfam THEN 1 ELSE 0 END
        + CASE WHEN has_ec THEN 1 ELSE 0 END
        + CASE WHEN has_role THEN 1 ELSE 0 END
        + CASE WHEN has_reaction THEN 1 ELSE 0 END
        + CASE WHEN has_transporter THEN 1 ELSE 0 END
        + CASE WHEN has_cazy THEN 1 ELSE 0 END) AS informative_count,
       has_any_edge
  SET g.annotation_state =
        CASE
          WHEN informative_count >= 2 THEN 'informative_multi'
          WHEN informative_count = 1 THEN 'informative_single'
          WHEN has_any_edge THEN 'catch_all_only'
          ELSE 'no_evidence'
        END,
      g.annotation_quality =
        CASE
          WHEN informative_count >= 2 THEN 3
          WHEN informative_count = 1 THEN 2
          WHEN has_any_edge THEN 1
          ELSE 0
        END
} IN TRANSACTIONS OF 500 ROWS;
```

The Cypher file must carry a clear comment-marked maintenance block (see "Source bucket maintenance" below) so future contributors adding a new functional Gene-edge type know where to update.

Phase consolidation: All 4 F1 ships now run in post-import, sharing `term.is_uninformative` as the single primitive. The catch-all YAML is loaded once (post-import only). `_compute_annotation_quality` deletes from `build_gene_annotations.py`; `_compute_gene_category` stays at build-time (its output `gene_category` is the input to the role bucket).

**Documented semantic shift.** `annotation_quality` semantics change:
- *Before*: 0–3 mixing product-name regex (hypothetical-named?) with arbitrary 4-source structured-count (GO/KEGG/EC/Pfam, presence-only).
- *After*: 0–3 = numeric encoding of `annotation_state`; pure informative-evidence richness over 5 source buckets, no product-name dependency.
- Existing `WHERE g.annotation_quality <= 1` queries silently shift meaning (today they catch hypothetical-named genes; refined they catch low-evidence genes — the spirit of the original filter).

**Drift risk: `gene_category` ↔ `uninformative_terms.yaml`.** Both encode catch-all role logic; `gene_category` does it via `COG_TO_CATEGORY` / `CYANORAK_TO_CATEGORY` / `TIGR_TO_CATEGORY` mappings in `build_gene_annotations.py`, the YAML does it as the post-import flag set. If the two drift (e.g., a new uninformative TIGR role added to the YAML but not mapped to "Unknown" in `TIGR_TO_CATEGORY`), the role bucket's signal goes wrong.

Mitigation (test in `test_uninformative_terms.py`): for every catch-all role entry in the YAML, assert that the corresponding entry in the role-mapping table maps to `"Unknown"`. Build-time check; failure means drift.

### Ship 1.3 — `Gene.annotation_state` (new categorical enum, 1:1 with `annotation_quality`)

New Gene property: `annotation_state: str` ∈ `{no_evidence, catch_all_only, informative_single, informative_multi}`. Set in the same post-import Cypher block as `annotation_quality` (Ship 1.2) — both fields share `informative_source_count` as their input. See the table and Cypher in Ship 1.2.

| State | Numeric `annotation_quality` | Plain-English criterion |
|---|---|---|
| `no_evidence` | 0 | no ontology/functional edges at all |
| `catch_all_only` | 1 | edges exist; every connected term is flagged `is_uninformative='true'` |
| `informative_single` | 2 | exactly 1 source bucket has informative terms |
| `informative_multi` | 3 | ≥ 2 source buckets have informative terms |

Each enum name maps directly to its criterion in plain English. Reading `annotation_state = 'catch_all_only'` tells the consumer the gene's evidence shape without consulting docs.

This explicitly carves out the `catch_all_only` bucket that F1's trigger analysis tried to build manually (and got wrong using source-presence as a proxy).

### Ship 1.4 — `Gene.informative_annotation_types: list[str]` (parallel field)

Added next to existing `Gene.annotation_types`. Domain (13 sources, all live as of 2026-05-02):

`go_bp`, `go_mf`, `go_cc`, `pfam`, `cog_category`, `kegg`, `brite`, `ec`, `cyanorak_role`, `tigr_role`, `reaction`, `transporter`, `cazy`.

Inclusion rule: a source appears in `informative_annotation_types` if at least one connected term in that source has `is_uninformative IS NULL` (or for sources that have no flagging — `pfam`, `ec`, `reaction`, `transporter`, `cazy` — at least one connected node exists).

Existing `Gene.annotation_types` is **NOT changed by this ship** — the F-AUDIT-2 layer 1 fix (2026-04-30) clarified its description as "presence-by-source — does NOT indicate content informativeness", and tightening it days later would invalidate that fix. (A separate decision to *extend* `annotation_types` with `reaction` / `transporter` / `cazy` for symmetry is left for a follow-up; it's an additive change, not a meaning-tightening change, but isn't strictly required for F1.)

Granularity note: this field keeps the 10-source domain (one entry per ontology-edge type) — *not* the 5-bucket coalescence used by `annotation_quality` / `annotation_state`. The two scales serve different purposes:

- `annotation_quality` / `annotation_state` are **coarse evidence-richness signals** (5 buckets, role-sources coalesced via `gene_category`).
- `informative_annotation_types` is a **granular routing field**, the parallel of `annotation_types`. A gene with informative `tigr_role` but flagged `cyanorak_role` should show `tigr_role` in the list, not "role".

Post-import Cypher pattern:

```cypher
MATCH (g:Gene)
CALL {
  WITH g
  SET g.informative_annotation_types =
    CASE WHEN EXISTS { (g)-[:Gene_involved_in_biological_process]->(t:BiologicalProcess)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_bp'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_enables_molecular_function]->(t:MolecularFunction)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_mf'] ELSE [] END +
    -- ... one CASE per source, mirroring existing annotation_types Cypher;
    -- pfam, ec, brite, pfam_clan have no flag check (never flagged).
} IN TRANSACTIONS OF 1000 ROWS;
```

Adjacent to the existing `annotation_types` block in [`scripts/post-import.cypher:498`](../../../scripts/post-import.cypher) and `scripts/post-import.sh`.

### Source bucket maintenance

The source-bucket list is **explicitly enumerated**, not auto-discovered from the schema. New functional Gene-edge types are NOT picked up automatically — adding a bucket requires:

1. Append a row to the bucket table in this spec (and any successor spec).
2. Add a `has_<bucket>` line to the post-import Cypher in [`scripts/post-import.cypher`](../../../scripts/post-import.cypher) and `scripts/post-import.sh`.
3. Include `has_<bucket>` in the `informative_count` sum.
4. Add the new edge type(s) to the `has_any_edge` predicate.
5. Update the CLAUDE.md "Gene annotation score" subsection (see "Documentation deliverables" below) to list the new bucket.
6. Bump the bucket-count test in `test_annotation_state.py` (assert that the live count of buckets matches the spec).

The Cypher block carries a fenced comment marker (`// SOURCE_BUCKETS:start` … `// SOURCE_BUCKETS:end`) so the maintenance points are greppable.

**Inclusion criterion for a new bucket.** A new Gene-edge type belongs in the score iff it carries **per-gene functional evidence** about what the gene's protein *does*. Examples that qualify: catalysis (`Gene_catalyzes_reaction`), enzyme classification (`Gene_catalyzes_ec_number`), domain membership (`Gene_has_pfam`), transporter classification (`Gene_has_tcdb_family`), CAZy classification (`Gene_has_cazy_family`). Examples that do NOT qualify: orthology grouping (`Gene_in_ortholog_group` — about evolutionary clustering, not function), study-specific clustering (`Gene_in_gene_cluster`), expression evidence (`Changes_expression_of`), genomic context (`start`/`end`/`strand`/`contig`).

If a new bucket's terms include catch-all entries (like KEGG's "uncharacterized protein" KOs), add an `is_uninformative` rule to `config/uninformative_terms.yaml` and a `WHERE t.is_uninformative IS NULL` clause to the Cypher; otherwise (every term is functionally specific) the bucket criterion is just `EXISTS { (g)-[:edge_type]->() }` like Pfam / EC / Reaction.

### Documentation deliverables (F1 release)

- `gene_overview` MCP tool about-content updated to describe both `annotation_quality` (refined) and `annotation_state` (new) + cross-reference between `annotation_types` and `informative_annotation_types`.
- `kg_schema` field descriptions updated for all four properties.
- **CLAUDE.md** updates:
  - "Key graph facts" — replace today's `annotation_quality` description with: numeric encoding of `annotation_state` (0=`no_evidence`, 1=`catch_all_only`, 2=`informative_single`, 3=`informative_multi`); list the 8 live source buckets (`go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter`, `cazy`); pointer to this spec.
  - "Gene properties" — list `annotation_state`, `informative_annotation_types`, `contributing_sources`, `contig`, `seed_ortholog`, `seed_ortholog_evalue` alongside the existing entries.
  - Cross-reference the "Source bucket maintenance" section above so future contributors editing the score know what to update.
- Release notes call out the `annotation_quality` semantic shift with a worked before/after example.

---

## F2 — Data source surfacing

### Ship 2.1 — `DataSource` node type

New node type. Schema additions to [`config/schema_config.yaml`](../../../config/schema_config.yaml):

```yaml
data source:
  represented_as: node
  preferred_id: data_source
  label_in_input: data_source
  properties:
    id: str                       # join key, e.g. 'eggnog'
    name: str                     # e.g. 'EggNOG-mapper'
    description: str
    version: str                  # optional v1; '' when not derivable
    scope: str                    # 'gene_level' | 'organism_restricted'
    provenance: str               # 'download' | 'tool_run'
    applies_to_organisms: str[]   # populated only for organism_restricted
    info_types: str[]             # auto-generated from gene_annotations_config.yaml
```

Initial node set (4 rows):

| id | name | scope | provenance | applies_to_organisms |
|---|---|---|---|---|
| `ncbi` | NCBI RefSeq | gene_level | download | [] |
| `cyanorak` | Cyanorak | organism_restricted | download | (Pro/Syn organism preferred_names) |
| `uniprot` | UniProt | gene_level | download | [] |
| `eggnog` | EggNOG-mapper | gene_level | tool_run | [] |

### Ship 2.2 — `gene_annotations_config.yaml` extension

One new per-source metadata key (`logical_sources`) under each entry in the `sources:` block of [`config/gene_annotations_config.yaml`](../../../config/gene_annotations_config.yaml). Always a list (even for single-logical-source files) for shape uniformity:

```yaml
sources:
  gene_mapping:
    type: csv
    path_pattern: "{data_dir}/gene_mapping.csv"
    description: "NCBI + Cyanorak merged gene annotations"
    # NEW:
    logical_sources:                  # gene_mapping is a merge of two logical sources
      - id: ncbi
        scope: gene_level
        provenance: download
      - id: cyanorak
        scope: organism_restricted
        provenance: download
        applies_to_organisms: [Prochlorococcus*, Synechococcus*, Parasynechococcus*, Thermosynechococcus*]
  eggnog:
    # ...existing...
    # NEW:
    logical_sources:
      - id: eggnog
        scope: gene_level
        provenance: tool_run
  uniprot:
    # ...existing...
    logical_sources:
      - id: uniprot
        scope: gene_level
        provenance: download
```

The `gene_mapping` case is the reason for the list shape: a single CSV carries two logical sources (NCBI columns + Cyanorak columns merged at download time). Other source files (eggnog, uniprot) declare a one-element list for uniformity. Per-field source attribution still uses the existing `source_label` mechanism — `logical_sources` is metadata for the new `DataSourceAdapter`, not field-merge logic.

### Ship 2.3 — `DataSourceAdapter`

New adapter (e.g., `multiomics_kg/adapters/data_source_adapter.py`) that:

1. Reads `config/gene_annotations_config.yaml`.
2. Inverts the `fields:` block to derive `info_types` per source (each source's `info_types` = list of all field names where that source label appears as `source` or `source_label`).
3. Emits one DataSource node per logical source.

Hooked into [`create_knowledge_graph.py`](../../../create_knowledge_graph.py) like other top-level adapters. Test mode: emits all 4 nodes regardless of `--test` flag.

### Ship 2.4 — `Gene.contributing_sources: list[str]`

Computed at build time in [`multiomics_kg/download/build_gene_annotations.py`](../../../multiomics_kg/download/build_gene_annotations.py) per gene:

```python
def _compute_contributing_sources(gene: dict) -> list[str]:
    """Determine which data sources contributed at least one field to this gene."""
    sources = set()
    # NCBI: every gene has it by definition (Gene exists ⇒ NCBI ran)
    sources.add('ncbi')
    # Cyanorak: any cyanorak-tagged field non-null
    if (gene.get('locus_tag_cyanorak')
        or _has_source_label(gene, 'cyanorak')):
        sources.add('cyanorak')
    # UniProt: uniprot_accession non-null OR any uniprot-tagged field non-null
    if (gene.get('uniprot_accession')
        or _has_source_label(gene, 'uniprot')):
        sources.add('uniprot')
    # EggNOG: any eggnog-tagged field non-null OR seed_ortholog non-null
    if (gene.get('seed_ortholog')
        or _has_source_label(gene, 'eggnog')
        or gene.get('eggnog_ogs')):
        sources.add('eggnog')
    return sorted(sources)
```

Where `_has_source_label(gene, label)` walks `gene["*_source"]` track fields and `[label]` prefixes in `alternate_functional_descriptions`.

Then add `contributing_sources` to:

- The Gene field-merge config (or set it as a computed field in `build_gene_annotations.py`).
- `GeneNodeField` enum + `properties` block in [`config/schema_config.yaml`](../../../config/schema_config.yaml) (`contributing_sources: str[]`).
- The cyanorak_ncbi adapter passthrough.

### Ship 2.5 — Eligibility note (no code)

Document in the post-import Cypher comment block + the new MCP data-source tool's docstring:

> For `provenance='tool_run'` sources, "absent from `Gene.contributing_sources`" means "the tool ran on this gene's protein and returned no hit". We do not currently distinguish "tool ran but failed" from "tool ran and returned no hit"; if that distinction becomes load-bearing, it's a separate property, not a redefinition of `contributing_sources`.

---

## F3 — `cluster_type` vocabulary refinement

Single rename + convention text. No schema work.

### Ship 3.1 — Rename `classification` → `expression_bin`

Touch points:

- `.claude/skills/paperconfig/validate_paperconfig.py`: update `VALID_CLUSTER_TYPES` set.
- `.claude/skills/paperconfig/SKILL.md` template: replace `classification` with `expression_bin` in examples.
- Tool docstrings for `list_clustering_analyses` and `gene_clusters_by_gene` (in the MCP repo) — surfaced as the convention:
  > `cluster_type=expression_bin` analyses do **not** carry per-cluster `functional_description`; cluster `name` is the metric label (e.g., VEG, HEG). All other `cluster_type` values carry per-cluster functional descriptions where curated.

Migration cost: zero. No current paperconfig uses `classification`. Validator vocab and skill template are the only touches.

### Future-ambiguity rule (documented)

If a future cluster type has unclear intent: **rename or split the `cluster_type` value** rather than adding a parallel intent field. Document in the paperconfig SKILL.md.

---

## F4 — Gene properties for sparsely-annotated genes

Three ships. Honest-framing note in release notes: these primarily benefit the broader sparsely-annotated population, not the strict floor (~14 TX50_RS-style genes which the framing originally implied).

### Ship 4.1 — `Gene.contig: str`

Required for genomic-neighbor lookups to be meaningful (current Gene has `start`/`end`/`strand` but no contig discriminator).

Plumbing path:

1. **NCBI download script** ([`multiomics_kg/download/download_genome_data.py`](../../../multiomics_kg/download/download_genome_data.py), step 5 `gene_mapping.csv` writer): preserve the GFF `seqid` column (currently dropped). Add a `seqid` column to the output CSV. (The CSV column is named `seqid` — matches GFF; the Gene-node property is named `contig` — matches biological vocabulary. Renaming happens in the merge config.)
2. **`gene_annotations_config.yaml`**: add a passthrough field:
   ```yaml
   contig:
     type: passthrough
     source: gene_mapping
     field: seqid
   ```
3. **Schema** ([`config/schema_config.yaml`](../../../config/schema_config.yaml)): `contig: str` under Gene properties.
4. **Adapter** ([`multiomics_kg/adapters/cyanorak_ncbi_adapter.py`](../../../multiomics_kg/adapters/cyanorak_ncbi_adapter.py)): add `CONTIG = 'contig'` to `GeneNodeField` enum.

### Ship 4.2 — `Gene.seed_ortholog: str` + `Gene.seed_ortholog_evalue: float`

Already in `gene_annotations_merged.json` (eggnog passthrough). Surface as Gene properties.

Touch points:

- [`config/schema_config.yaml`](../../../config/schema_config.yaml): add both fields under Gene properties.
- [`multiomics_kg/adapters/cyanorak_ncbi_adapter.py`](../../../multiomics_kg/adapters/cyanorak_ncbi_adapter.py): add `SEED_ORTHOLOG = 'seed_ortholog'` and `SEED_ORTHOLOG_EVALUE = 'seed_ortholog_evalue'` to `GeneNodeField`. Add `seed_ortholog_evalue` to the `float_fields` set in `_get_gene_nodes()`.

Field description (for `kg_schema` output and `gene_overview` about-content):

> `seed_ortholog`: The eggNOG-mapper seed ortholog — the best-scoring protein hit in eggNOG's reference DB for this gene's protein. Format: `<taxid>.<source_identifier>`, where `<source_identifier>` varies by organism (most commonly a locus tag, but `<contig>_<gene_n>` for draft assemblies, or RefSeq WP_ depending on what's in eggNOG's reference for that taxid). Examples: `160488.PP_0894` (Pseudomonas putida KT2440 locus); `225937.HP15_1417` (Marinobacter HP15 self-match); `1122222.AXWR01000046_gene1504` (Meiothermus reference draft assembly). The taxid prefix is an NCBI Taxonomy ID; the suffix is a stable identifier in that organism's eggNOG reference. Useful as a "this gene's protein most resembles X (E=Y)" pointer when functional annotations are absent or sparse.

### Ship 4.3 — Genomic-neighbor MCP tool (KG-side: nothing)

KG side has no work here once `contig` lands. The Cypher pattern for the explorer's eventual `gene_neighbors` MCP tool:

```cypher
MATCH (g:Gene {locus_tag: $lt})-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
MATCH (n:Gene)-[:Gene_belongs_to_organism]->(o)
WHERE n.contig = g.contig
  AND abs(n.start - g.start) <= $window_bp     // OR: ±N flanking via subquery on ORDER BY
RETURN n ORDER BY n.start
```

No materialized `(Gene)-[:NEIGHBOR_OF]->(Gene)` edges. Window-default decision (gene-count vs bp) is the explorer's call.

---

## Schema deltas summary

### New nodes

| Node type | Properties | Source ship |
|---|---|---|
| `DataSource` | `id`, `name`, `description`, `version`, `scope`, `provenance`, `applies_to_organisms`, `info_types` | 2.1 |

### New edges

None. (`DataSource` ↔ `Gene` is soft-joined via property values, not edges, per the F2 joint resolution.)

### Modified node properties

| Node | Property | Action | Ship |
|---|---|---|---|
| `Gene` | `annotation_quality: int` | semantic refinement (catch-all-aware); computation moves from build-time to post-import | 1.2 |
| `Gene` | `annotation_state: str` | NEW (post-import enum) | 1.3 |
| `Gene` | `informative_annotation_types: str[]` | NEW (post-import parallel) | 1.4 |
| `Gene` | `contributing_sources: str[]` | NEW (build-time) | 2.4 |
| `Gene` | `contig: str` | NEW (build-time, requires NCBI plumbing) | 4.1 |
| `Gene` | `seed_ortholog: str` | NEW (build-time passthrough) | 4.2 |
| `Gene` | `seed_ortholog_evalue: float` | NEW (build-time passthrough) | 4.2 |
| `BiologicalProcess`, `MolecularFunction`, `CellularComponent`, `CogFunctionalCategory`, `CyanorakRole`, `TigrRole`, `KeggTerm` | `is_uninformative: str` (sentinel `"true"` / absent) | NEW (post-import flag) | 1.1 |

### New config files

| Path | Purpose | Ship |
|---|---|---|
| `config/uninformative_terms.yaml` | Catch-all term vocabulary | 1.1 |

### Modified config files

| Path | Change | Ship |
|---|---|---|
| `config/gene_annotations_config.yaml` | Per-source `logical_sources:` list (id, scope, provenance, applies_to_organisms) | 2.2 |
| `config/schema_config.yaml` | Gene property additions + new `data source` node type | multiple |

### New adapters

| Path | Purpose | Ship |
|---|---|---|
| `multiomics_kg/adapters/data_source_adapter.py` | Emit `DataSource` nodes from gene_annotations_config.yaml | 2.3 |

### Modified adapters / scripts

| Path | Change | Ship |
|---|---|---|
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | `GeneNodeField` enum + passthrough for new properties | 2.4, 4.1, 4.2 |
| `multiomics_kg/download/download_genome_data.py` | Preserve GFF `seqid` in `gene_mapping.csv` | 4.1 |
| `multiomics_kg/download/build_gene_annotations.py` | Delete `_compute_annotation_quality` (moves to post-import); new `_compute_contributing_sources` | 1.2, 2.4 |
| `scripts/post-import.cypher` and `scripts/post-import.sh` | New Cypher blocks for `is_uninformative`, `annotation_quality` (refined + relocated), `annotation_state`, `informative_annotation_types` | 1.1, 1.2, 1.3, 1.4 |
| `.claude/skills/paperconfig/validate_paperconfig.py` | `VALID_CLUSTER_TYPES` rename | 3.1 |
| `.claude/skills/paperconfig/SKILL.md` | Vocab + convention text | 3.1 |

---

## Sequencing / release plan

The four resolutions land as a **single coordinated KG release**. F1 in particular MUST land all four ships together (per the explorer review's "calibrate once" request); shipping ships 1.1–1.4 separately would force consumers to rebaseline twice.

Recommended commit order (each commit independently reviewable; ship order ≠ commit order):

1. **F3 first** — vocab rename (lowest risk; sets up paperconfig skill template).
2. **F4 first half** — NCBI download script + gene_mapping.csv `seqid` plumbing (upstream; needs prepare_data step 0 rerun with `--force`).
3. **F4 second half** — schema + adapter additions for `contig` / `seed_ortholog` / `seed_ortholog_evalue`.
4. **F2** — `gene_annotations_config.yaml` `logical_sources` extension + `DataSourceAdapter` + `Gene.contributing_sources`.
5. **F1.1** — `config/uninformative_terms.yaml` + post-import Cypher to set `is_uninformative` on term nodes. Includes the drift test against `*_TO_CATEGORY` mappings.
6. **F1.2 + F1.3** — single post-import Cypher block setting `annotation_quality` and `annotation_state` from the shared `informative_source_count`. Deletes `_compute_annotation_quality` from `build_gene_annotations.py` (function relocates entirely to Cypher).
7. **F1.4** — separate post-import Cypher block for `informative_annotation_types` (granular 10-source domain).
8. **Doc updates** — CLAUDE.md, release notes, MCP tool docstrings.

Acceptance gates:

- KG validity tests (`pytest -m kg`) green after rebuild.
- New tests:
  - `test_uninformative_terms.py` — (a) every YAML `ids:` entry resolves to an existing node post-build with `is_uninformative='true'`; (b) `name_patterns:` regexes flag at least one node each; (c) **drift check**: every catch-all role entry in the YAML's `cog_category` / `cyanorak_role` / `tigr_role` sections, when fed as a gene's only role, must yield `gene_category='Unknown'` from the `*_TO_CATEGORY` mappings in `build_gene_annotations.py`.
  - `test_annotation_state.py` — every Gene has exactly one of the 4 enum values; `annotation_quality` numeric matches the expected 0–3 mapping; sample assertions per state on known genes.
  - `test_data_source.py` — 4 DataSource nodes present; `info_types` non-empty per node; `Gene.contributing_sources` always contains `'ncbi'`; cyanorak `applies_to_organisms` populated only for `organism_restricted`.
  - `test_gene_contig.py` — every Gene has a non-null `contig`; same `contig` for genes from the same NCBI accession.
- Snapshot regeneration (`uv run python tests/kg_validity/generate_snapshot.py`) — expected churn on Gene properties.
- Spot-check Cypher for floor-case TX50_RS09500 (verify outcome post-rebuild — expected to land in `no_evidence` or `catch_all_only` depending on whether eggNOG returned any catch-all hits):
  ```cypher
  MATCH (g:Gene {locus_tag: 'TX50_RS09500'})
  RETURN g.contributing_sources, g.annotation_quality, g.annotation_state,
         g.informative_annotation_types
  // expected (typical): contributing_sources=['ncbi'], annotation_quality ∈ {0,1},
  //                     annotation_state ∈ {'no_evidence', 'catch_all_only'},
  //                     informative_annotation_types=[]
  ```

---

## Acceptance criteria

- All 4 frictions' KG-side ships land in one release.
- Existing `pytest -m kg` suite still green.
- New tests added (per above).
- `omics-edge-snapshot` skill shows no expression-edge regressions.
- KG release notes include:
  - F1 worked example showing `annotation_quality` + `annotation_state` before/after. Realistic case: gene with `pfam_ids=[PF00712]` + `go_terms=[GO:0003674]` (MF root) + `gene_category='Unknown'`, and no other functional-source edges (no kegg_ko, ec_number, cyanorak_role, tigr_role, reaction, etc.). *Before*: `structured_count = 2` (`go_terms` + `pfam_ids` both non-null) and product is real → score `3`. *After*: GO MF root is flagged uninformative (`is_uninformative = 'true'`), so `has_go = false`; only `has_pfam` is true among the 8 live buckets → `informative_source_count = 1` → `annotation_quality = 2`, `annotation_state = 'informative_single'`. Population most affected: genes whose ontology hits are dominated by GO roots / COG S / TIGR catch-alls.
  - F2 worked example showing `contributing_sources` for TX50_RS09500 vs MED4 dnaA.
  - F3 vocab rename callout.
  - F4 honest-framing note (broader-population leverage, not strict floor).

## Risks

- **Term-node label binding for `is_uninformative` flagging.** The YAML keys (`biological_process`, `molecular_function`, etc.) must map exactly to BioCypher's PascalCase node labels. Validation test catches this.
- **`Gene.annotation_quality` baseline shift.** Documented in release notes; explorer-side regression baseline rebuild expected and accepted (per joint resolution).
- **`gene_annotations_config.yaml` schema extension.** The new `logical_sources` key is read only by the new `DataSourceAdapter`; it doesn't affect the existing field-merge consumers in `build_gene_annotations.py`. Adapter should fail loudly if `logical_sources` is missing or malformed rather than fall back.
- **NCBI download cache invalidation for `seqid`.** Step 5 (`gene_mapping.csv`) needs to be re-run with `--force` after the script change. Document in the spec's implementation plan.

## Future work (deferred specs)

- **Ortholog-mediated annotation transfer** (F4 option c) — its own design conversation when broader sparsely-annotated population becomes the analytical bottleneck. Higher leverage than F4 floor cases per se.
- **`KGRelease` versioning node** — when a use case actually needs cross-source release context (current per-DataSource `version` carries enough for v1).
- **Quality / score signal per (gene, pipeline)** — per F2 anti-scope; no current demand.
- **Tool-failure vs no-hit distinction** for `provenance='tool_run'` sources — separate property if it ever earns its keep.
