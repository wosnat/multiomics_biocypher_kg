# TigrRole hierarchy + NcbifamFamily → TigrRole bridge — design

**Date:** 2026-08-28
**Status:** draft, awaiting review
**Scope class:** architectural (new edge type, new node level, one merge rule)

## 1. Problem

`TigrRole` is the KG's only hierarchical ontology with no `is_a` edges: 114
flat nodes (`level = 0` everywhere) whose two-level JCVI mainrole/subrole
scheme is embedded in compound names (`"Energy metabolism / TCA cycle"`).
Its only source is the Cyanorak GFF `tIGR_Role` field, so only the 19
Cyanorak-annotated Pro/Syn strains carry `Gene_has_tigr_role` edges; the 24
other strains (all Alteromonas, heterotrophs, PCC/UTEX/BP1, MIT1314, RSP50)
have none.

Meanwhile every strain carries `Gene_has_ncbifam_family` edges to 2,204
observed `TIGR*` `NcbifamFamily` nodes. When NCBI absorbed TIGRFAMs into
NCBIfam it kept the accessions (only 2 of the 2,963 role-bearing TIGRFAMs 15.0
families are retired today) but dropped the role metadata — `hmm_PGAP.tsv`
has no role column. The role assignments survive only in JCVI's frozen
release-15.0 archive on the NCBI FTP.

## 2. Measurements (2026-08-28, all 43 strains / 127,458 genes)

Archive: `TIGRFAMS_ROLE_LINK` (2,963 families → 117 role ids) +
`TIGR_ROLE_NAMES` (116 named roles; role `719` is unnamed).

| Measure | Value |
|---|---|
| Observed `TIGR*` NcbifamFamily nodes with an archive role | **1,720 / 2,204** (78%) |
| Observed `NF*` nodes (can never link — post-date the role system) | 2,753 |
| Distinct roles reached | 108 — 104 already `TigrRole` nodes, 4 new (`110`, `186`, `187`, `719`) |
| Main-role name concordance, Cyanorak `tIGR_Role_description` vs archive, shared codes | **110 / 110** identical |
| Cyanorak-only codes not in archive | `128`, `270`, `701`, `856` |
| Non-Cyanorak genes reaching a role via 2-hop | **19,115 / 75,996** (25%) — NCBIfam-hit-bounded |
| Non-Cyanorak `gene_category = Unknown` → category via role | **1,032** (+40 more if Cyanorak `protein_domains` TIGR tokens were also read) |
| Cyanorak-strain `Unknown` → category | 101 |
| Genes whose `annotation_state` would change | **12** — the bridge only reaches genes that already have an NCBIfam edge, which is already an informative bucket |
| Non-Cyanorak genes whose TIGR-derived category differs from the current COG-derived one | 7,837 of 19,115 (41%) — schemes carve biology differently; NOT an error rate |

Dropped sources (checked because the merge discards them): Cyanorak
`protein_domains` TIGR tokens are 93% already in the gene's InterProScan hits;
PGAP `inference` HMM ids are 84% `NF*` `domain`/`PfamEq`-type models that
InterPro's NCBIfam member DB excludes **by design** (0 / 1,860 `PfamEq` and
108 / 14,234 `domain` families ever observed across 43 genomes; the ids are
numerically interleaved with observed ones, so not a version gap). Neither
moves the numbers above. Multi-source `Gene_has_ncbifam_family` is **out of
scope** and recorded in `plans/backlog.md`.

## 3. Decisions

1. **Hierarchy**: `TigrRole` becomes two-level. Existing subrole nodes keep
   their ids (`tigr.role:<id>`) and their compound `name` (non-breaking for
   name consumers); they get `level = 1`, `level_kind = 'tigr_subrole'`. New
   mainrole nodes `tigr.mainrole:<slug>` (`level = 0`,
   `level_kind = 'tigr_mainrole'`, `name` = mainrole text, `code` = slug) and
   `Tigr_role_is_a_tigr_role` (subrole → mainrole). Cyanorak-only codes whose
   description has no `" / "` (e.g. `856` "Not Found") stay level-0 roots
   with no parent.
2. **Bridge edge** `Ncbifam_family_has_tigr_role` (NcbifamFamily → TigrRole),
   ~1,720 edges, `TIGR*` only, **no properties** (vocabulary-contract R3/R5:
   the archive is a single frozen source; provenance is documented on the
   edge type, not repeated per edge). Verb `has`, same composition semantics
   as `Tcdb_family_has_pfam_domain`: "JCVI assigned this family this role."
   It is a **router** — read outward from a gene's NCBIfam family it suggests
   a role; it never asserts a role on a gene. No `Gene_has_tigr_role` edges
   are created from it.
3. **Role `719`** (unnamed in the archive, 43 observed families) is dropped —
   an unnamed role node is noise. Logged at build time.
4. **New role nodes** (`110` Energy metabolism / Anaerobic, `186` Mobile … /
   Plasmid functions, `187` Cellular processes / Pathogenesis) are emitted
   only because a bridge edge reaches them (observed-only, TCDB pattern).
5. **`gene_category` fill-only rule** (build-time, step 2): after the
   existing Cyanorak → TIGR → COG chain, if the result is still `Unknown`,
   map the gene's `ncbifam_ids` → archive role → mainrole →
   `TIGR_TO_CATEGORY`; take the most frequent non-`Unknown` category (ties →
   alphabetical first). It can only turn `Unknown` into a category, never
   change an existing one — zero churn on the 7,837 disagreements. Expected
   +~1,130 genes. Priority order becomes Cyanorak role → Cyanorak TIGR role →
   COG → **NCBIfam-bridged TIGR role**. (Placed after COG deliberately: COG is
   per-gene orthology; the bridge is family-level inference.)
6. **Uninformative flags**: the 5 existing subrole flags stay; the mainrole
   nodes for "Hypothetical proteins", "Unknown function", "Unclassified",
   "Not Found" are added to `config/uninformative_terms.yaml` and the
   post-import F1.1 list. `141`/`703` stay unflagged (class known).
7. **Not done**: `annotation_types` / `informative_annotation_types` /
   `annotation_quality` are untouched (bridge is ontology→ontology; measured
   effect 12 genes). Enrichment backgrounds over `TigrRole` for heterotrophs
   must be documented as coverage-biased (25% vs ~90%).

## 4. Components

### 4.1 Reference data — prepare_data step 9 (`build_ncbifam_reference.py`)

- Downloads (with `--refetch-raw`) into gitignored `cache/data/ncbifam/raw/`:
  `https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMS_ROLE_LINK`
  and `.../TIGR_ROLE_NAMES` (45 KB + 11 KB).
- Writes a **separate** committed file `cache/data/ncbifam/tigr_roles.json`
  (keeping `ncbifam_reference.json`'s flat `{acc: {...}}` shape untouched so
  no consumer iterating entries sees a foreign key):
  ```json
  {"release": "TIGRFAMs 15.0 (frozen 2018)",
   "roles": {"132": {"mainrole": "DNA metabolism",
                     "sub1role": "DNA replication, recombination, and repair"}},
   "family_role": {"TIGR00001": "158", ...}}
  ```
  Roles with no `mainrole` name (`719`) are excluded from both maps.
- Pure parser in `multiomics_kg/utils/ncbifam.py` (`parse_tigr_role_link`,
  `parse_tigr_role_names`), fail-loud on zero rows.
- Outage tolerance: if the FTP download fails and `tigr_roles.json` exists,
  warn and reuse it (TCDB precedent); only a missing file is fatal.

### 4.2 Nodes + hierarchy — `functional_annotation_adapter.MultiCogRoleAnnotationAdapter`

- New constructor arg `extra_tigr_roles: dict[str, str] | None` — `{code:
  compound_name}` for roles reached by the bridge but absent from Cyanorak
  data (computed by `create_knowledge_graph`, §4.4).
- `get_nodes()` emits, for the union of observed Cyanorak codes ∪
  `extra_tigr_roles`: the subrole node (`level 1`, `level_kind
  'tigr_subrole'`, name unchanged) and, from the `" / "` split of the
  compound name, a deduplicated mainrole node (`tigr.mainrole:<slug>`, `level
  0`, `level_kind 'tigr_mainrole'`). Codes without `" / "` → level 0,
  `level_kind 'tigr_mainrole'`, no parent.
- `get_edges()` adds `tigr_role_is_a_tigr_role` (subrole → mainrole).
- Exposes `tigr_role_node_ids() -> set[str]` for the bridge's dangling guard.
- Slug: lowercase, non-alphanumerics → `_`, collapsed (`"Purines, pyrimidines,
  nucleosides, and nucleotides"` → `purines_pyrimidines_nucleosides_and_nucleotides`).

### 4.3 Bridge edges — `ncbifam_adapter.MultiNcbifamAdapter`

- New constructor arg `tigr_role_node_ids: set[str] | None` (TCDB
  `pfam_node_ids` contract: `None` → emit no bridge edges).
- `download_data()` also loads `tigr_roles.json`.
- `get_edges()` emits `ncbifam_family_has_tigr_role` for every observed
  `TIGR*` accession with a `family_role` whose target id is in
  `tigr_role_node_ids`; edge id `{acc}-tigrrole-{role}`; no properties.
  Logs emitted / skipped-unnamed / skipped-no-node counts.

### 4.4 Orchestration — `create_knowledge_graph.py`

Order today: `cog_role_adapter` (line ~256) runs before `ncbifam_adapter`
(~345). To keep that order:

1. Load `tigr_roles.json` once.
2. Compute the observed NCBIfam id set early (`MultiNcbifamAdapter` is
   constructed before `cog_role_adapter.get_nodes()` — its `_observed_ids()`
   reads only calls.json, no InterPro dependency; the InterPro kept-id set is
   still passed for the existing bridge).
3. `extra_tigr_roles = {role: f"{main} / {sub}" for observed TIGR accs' roles}
   − Cyanorak-observed codes` → pass to `MultiCogRoleAnnotationAdapter`.
4. After `cog_role_adapter.get_nodes()`, pass
   `cog_role_adapter.tigr_role_node_ids()` to the NCBIfam adapter.

### 4.5 `gene_category` fill — `build_gene_annotations.py`

- `_compute_gene_category(result, tigr_roles: dict | None = None)`: new
  priority 4 as in §3.5. `tigr_roles` is loaded lazily next to `ncbifam_ref`
  (step 9 artefact; absent → rule is a no-op with one warning).
- `TIGR_TO_CATEGORY` unchanged (already covers every archive mainrole; a
  build-time assertion checks this so an archive refresh can't introduce an
  unmapped mainrole silently).

### 4.6 Post-import (`scripts/post-import.sh` + `.cypher`, identical logic)

- `TigrRole` rollup switches to the CyanorakRole subtree form
  (`*0..` over `Tigr_role_is_a_tigr_role`): `gene_count` (subtree),
  `direct_gene_count`, `organism_count`.
- New `TigrRole.ncbifam_family_count` (int, subtree count of incoming bridge
  edges; 0 when none).
- Indexes `tigr_role_level_idx`, `tigr_role_level_kind_idx`;
  `tigrRoleFullText` unchanged.
- F1.1 uninformative list extended with the 4 mainrole ids (§3.6).

### 4.7 Schema + vocabulary

- `config/schema_config.yaml`: `tigr role` gains `level_kind: str`; new
  `tigr role hierarchical association` (`Tigr_role_is_a_tigr_role`) and
  `ncbifam family to tigr role association` (`Ncbifam_family_has_tigr_role`,
  no properties).
- `config/controlled_vocabularies.yaml`: `TigrRole.level_kind`
  (`tigr_mainrole | tigr_subrole`, closed).
- `config/uninformative_terms.yaml`: 4 mainrole ids.

## 5. Data flow

```
NCBI FTP archive ──step 9──▶ cache/data/ncbifam/tigr_roles.json (committed)
                                   │                      │
                    step 2 merge ◀─┘                      └─▶ create_knowledge_graph
                    gene_category fill                         ├─ CogRole adapter: TigrRole nodes (2 levels) + is_a
                    (Unknown only)                              └─ Ncbifam adapter: Ncbifam_family_has_tigr_role
                                                          post-import: subtree rollups, ncbifam_family_count, flags
```

## 6. Error handling

- Archive download failure: reuse committed `tigr_roles.json`, warn.
- Archive role with no name: excluded at build (logged count).
- Bridge target with no TigrRole node: skipped (dangling guard), counted.
- Mainrole not in `TIGR_TO_CATEGORY`: assertion at step 2 (fail loud).
- `tigr_roles.json` absent at step 2: warn once, fill rule disabled.

## 7. Testing

- Unit: archive parsers (fixture of 5 lines each, incl. an unnamed role);
  `MultiCogRoleAnnotationAdapter` node/edge emission with `extra_tigr_roles`
  (mainrole dedup, no-`" / "` root case, slug); `MultiNcbifamAdapter` bridge
  (None → 0 edges, dangling target skipped, NF* ignored);
  `_compute_gene_category` fill (Unknown → category, existing category
  untouched, tie rule, missing `tigr_roles`).
- Static: `tests/test_controlled_vocab.py` picks up the new entry;
  `tests/test_annotation_quality_buckets.py` unchanged (no bucket change).
- KG validity: `test_ontology_level.py` — `TigrRole` leaves the flat list;
  new assertions: ≥ 1,600 bridge edges, every bridge source is `TIGR*`, every
  level-1 node has exactly one parent, `direct_gene_count <= gene_count`,
  mainrole `gene_count` = union of children. Regenerate `snapshot_data.json`.
  `capture_annotation_state.py --save/--compare` (expect ≤ 12 moves).
- Edge regression: `/omics-edge-snapshot` before/after (no expression change
  expected).

## 8. Breaking / documentation

- **Breaking**: `TigrRole` count 114 → ~140 (3 new subroles + ~22 mainroles);
  `level` no longer all 0; `gene_count` on subroles unchanged, but the new
  mainrole nodes carry subtree counts. `gene_category` changes for ~1,130
  genes (`Unknown` → category only).
- CHANGELOG `### Breaking` + `### Added`; CLAUDE.md (TigrRole bullet, label
  lists, `gene_category` priority chain, step 9 outputs, Data Locations);
  `docs/kg-changes/tigr-role-bridge.md` (router semantics, coverage caveat,
  archive provenance); `plans/backlog.md` — amend the 2026-08-18 rejection
  (gene-level rejected; ontology-level bridge shipped) and add the
  multi-source NCBIfam measurement; `docs/kg-changes/vocabulary-contract.md`.
- MCP/explorer surfacing (2-hop role lookup for heterotrophs) is a separate
  follow-up.
