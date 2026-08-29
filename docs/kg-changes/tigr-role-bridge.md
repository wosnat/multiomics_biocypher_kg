# TigrRole hierarchy + NcbifamFamily → TigrRole bridge + inferred gene roles

**Date:** 2026-08-29
**Spec:** [`docs/superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md`](../superpowers/specs/2026-08-28-tigrrole-hierarchy-ncbifam-bridge-design.md)
**Track:** ontology hygiene (hierarchy) · ontology→ontology bridge · gated gene-edge inference

## What changed

Until now `TigrRole` was the KG's only hierarchical ontology with no `is_a`
edges — 114 flat nodes (`level = 0` everywhere) whose two-level JCVI
mainrole/subrole scheme was buried inside compound names
(`"Energy metabolism / TCA cycle"`). And its only source was the Cyanorak GFF
`tIGR_Role` field, so only the 22 Cyanorak-annotated *Prochlorococcus* /
*Synechococcus* strains carried `Gene_has_tigr_role` edges at all; the other 21
strains (every *Alteromonas*, every heterotroph, PCC/UTEX/BP1, MIT1314, RSP50)
had none.

Both are fixed here, from one source: JCVI's frozen TIGRFAMs 15.0 role archive.

### 1. `TigrRole` is two-level

- **Subroles** keep their ids (`tigr.role:<numeric code>`) and their compound
  `name` — non-breaking for name consumers — and gain `level = 1`,
  `level_kind = 'tigr_subrole'`.
- **Mainroles** are new nodes: `tigr.role:<slug>` (e.g.
  `tigr.role:energy_metabolism`), `level = 0`,
  `level_kind = 'tigr_mainrole'`, `name` = the mainrole text,
  `code` = the slug.
- **`Tigr_role_is_a_tigr_role`** (subrole → mainrole), one parent per subrole.
- Cyanorak-only codes whose description has no `" / "` separator (e.g. `856`
  "Not Found") stay level-0 roots with `level_kind = 'tigr_mainrole'` and no
  parent.

Mainroles reuse the single existing `tigr.role` prefix rather than minting a
second one (the CyanorakRole precedent — `cyanorak.role:R` / `R.1`). `tigr.role`
is house-minted: it is not a bioregistry prefix, and the registered `tigrfam`
prefix is the *family* namespace (`^TIGR\d+$`), not roles. No id collision is
possible (subrole ids are numeric, slug ids alphabetic), but two Cyanorak-only
level-0 roots keep numeric ids (`856`, `270`) — `level_kind` is the only
reliable level discriminator.

### 2. `Ncbifam_family_has_tigr_role` — the ontology bridge

New edge `NcbifamFamily → TigrRole`, **no properties**, one per (family, role),
`TIGR*` families only. Edge id `{accession}-has_tigr_role-{role_code}`.
Dangling-proof: `create_knowledge_graph` injects
`MultiCogRoleAnnotationAdapter.tigr_role_node_ids()` into the NCBIfam adapter
(the TCDB `pfam_node_ids` contract — `None` → no bridge edges at all).

Why the archive is the only source: when NCBI absorbed TIGRFAMs into NCBIfam it
kept the **accessions** (only 2 of the 2,963 role-bearing TIGRFAMs 15.0 families
are retired today) but **dropped the role metadata** — `hmm_PGAP.tsv` has no
role column. The assignments survive only in JCVI's frozen release-15.0 archive
on the NCBI FTP, which step 9 now downloads and freezes into a committed
artifact.

The verb is `has`, with the same composition semantics as
`Tcdb_family_has_pfam_domain`: *"JCVI assigned this family this role."* No
properties, per vocabulary-contract R3/R5 — the archive is a single frozen
source, so provenance belongs on the edge type, not repeated on every edge.

### 3. Equivalog-gated inferred `Gene_has_tigr_role`

For every gene, for every `ncbifam_ids` accession whose NCBIfam
`family_type` is **`equivalog`** and that carries an archive role, the graph
emits `Gene_has_tigr_role` with `sources: ['interproscan']`,
`evidence: 'family_inferred'`.

- Edge id `{locus_tag}-tigrrole-{code}` — **the same id the curated Cyanorak
  edge uses.** When both name the same role they are **merged adapter-side into
  ONE edge** with `sources: ['cyanorak', 'interproscan']` and
  `evidence: 'curated'` (curated wins the ladder; the TCDB eggNOG+diamond
  precedent). The merge must happen in the adapter, not by relying on BioCypher
  dedup, which keeps the *first* property set.
- When the two disagree the gene carries **two** edges and the disagreement is
  visible in the graph (see Semantics).
- Non-equivalog families (`subfamily`, `equivalog_domain`, `hypoth_equivalog`,
  `domain`, …) **never** produce a gene edge. The bridge (§2) is their only
  route.
- Junk roles (156 "Unknown function", 157, …) are emitted like any other role —
  the nodes are already flagged `is_uninformative`, and dropping the edges would
  hide the fact that a family *is* "hypothetical".

Result: `Gene_has_tigr_role` now spans **all 43 organisms**, not just the 22 Cyanorak-curated ones.

### 4. `gene_category` fill-only + `[tigr_role_inferred]` description lines

Both in the step-2 merge (`build_gene_annotations.apply_tigr_role_inference`),
behind the same equivalog gate.

- **`gene_category` priority chain becomes** Cyanorak role → Cyanorak TIGR role
  → COG → **NCBIfam-bridged TIGR role** (new priority 4). It is **fill-only**:
  it can only turn `Unknown` into a category, never change an existing one, so
  the 7,837 measured TIGR-vs-COG disagreements produce zero churn. Ties among
  the gene's roles break alphabetically. Placed after COG deliberately — COG is
  per-gene orthology, the bridge is family-level inference.
- **`[tigr_role_inferred] <Main> / <Sub>`** lines are appended to
  `alternate_functional_descriptions`, one per distinct inferred role,
  deduplicated against the curated `[tigr_role]` text so a Cyanorak gene whose
  archive role matches gets no duplicate. The label is deliberately distinct
  from `[tigr_role]` so curated-vs-inferred stays readable in the text, mirroring
  the `evidence` ladder. Same pattern as the existing `[interpro]` / `[pfam]`
  family-level lines. Effect: ~13.7K heterotroph genes become full-text
  searchable by role via `geneFullText` — no graph-structure change.
- `_check_tigr_roles_mapped` asserts at build time that every archive mainrole
  is present in `TIGR_TO_CATEGORY`, so an archive refresh cannot silently
  introduce an unmapped mainrole that lands as `Unknown`.

### 5. Reference artifact — prepare_data step 9

`build_ncbifam_reference.py::build_tigr_roles` downloads
`TIGRFAMS_ROLE_LINK` + `TIGR_ROLE_NAMES` from
`https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/` (45 KB + 11 KB, into
gitignored `cache/data/ncbifam/raw/`) and writes the **committed**
`cache/data/ncbifam/tigr_roles.json`:

```json
{"release": "TIGRFAMs 15.0 (frozen 2018)",
 "roles": {"132": {"mainrole": "DNA metabolism",
                   "sub1role": "DNA replication, recombination, and repair"}},
 "family_role": {"TIGR00001": ["158"], "TIGR00202": ["141", "271"]}}
```

It is a **separate file** from `ncbifam_reference.json`, whose flat
`{accession: {...}}` shape stays untouched so no consumer iterating entries sees
a foreign key. `family_role` values are **lists** — 294 of the 2,862 families
carry two or three roles (e.g. CsrA = *Glycolysis* + *RNA interactions*) — and
both/all are honoured everywhere: one bridge edge per (family, role), and gene
inference fans out. Role `719` is unnamed in the archive (102 link lines) and is
excluded from both maps; an unnamed role node is noise.

Outage tolerance follows the TCDB precedent: if the FTP download fails and
`tigr_roles.json` exists, step 9 warns and reuses the committed file; only a
*missing* file is fatal. If the file is absent at step 2, the merge warns once
and both the `gene_category` fill and the `[tigr_role_inferred]` lines become
no-ops.

### 6. Post-import

- `TigrRole` rollups switch to the CyanorakRole subtree form (`*0..` over
  `Tigr_role_is_a_tigr_role`): **`gene_count`** (subtree),
  **`direct_gene_count`**, **`organism_count`**.
- New **`TigrRole.ncbifam_family_count`** (int, subtree count of incoming
  bridge edges; 0 when none).
- New indexes `tigr_role_level_idx`, `tigr_role_level_kind_idx`;
  `tigrRoleFullText` unchanged.
- F1.1 uninformative flags gain the three junk **mainrole** nodes
  `tigr.role:hypothetical_proteins`, `tigr.role:unknown_function`,
  `tigr.role:unclassified`, alongside the 5 existing subrole flags (`156`,
  `157`, `185`, `704`, `856`). `141` / `703` stay unflagged — the class is
  known there.
- `annotation_types` / `informative_annotation_types` `tigr_role` rules are
  **unchanged**: they key off `Gene_has_tigr_role` existence plus
  `is_uninformative`, so the inferred edges join automatically with no new rule.

### 7. Vocabulary

- **New**: `TigrRole.level_kind` (`tigr_mainrole | tigr_subrole`, closed).
- **Widened**: `Gene_has_tigr_role.sources` → `[cyanorak, interproscan]`
  (closed, R2-joinable — `data_source:interproscan` exists);
  `Gene_has_tigr_role.evidence` → `[curated, family_inferred]`, two rungs of the
  shared KG-SYNC-005 ladder `curated > signature > homology > family_inferred >
  domain_inferred`.
- `Ncbifam_family_has_tigr_role` carries **no properties** and therefore
  declares no vocabulary entries.

## Semantics — read this before querying

**The bridge is a router, not an annotation.** Read *outward* from a gene's
known NCBIfam family, `Ncbifam_family_has_tigr_role` suggests a role. Read
*backward* — "this family maps to role X, therefore any gene hitting it does
role X" — it is only safe for the equivalog subset, which is exactly why the
gene edges are gated and the bridge is not. For the ~1,000 non-equivalog
role-bearing families the graph **never asserts a gene role**; the bridge is the
only route, and it is the ontology-level record of *why* an inferred gene role
exists, not a second one.

**Why the equivalog gate.** `equivalog` means "same function in every member" —
the only NCBIfam family type JCVI meant for unconditional role transfer. Gating
on it costs about a third of heterotroph reach and buys a qualitatively
different edge:

| | all role-bearing types | **equivalog only** |
|---|---|---|
| Non-Cyanorak genes gaining a role edge | 20,869 | **13,667** (18%) |
| Genes with >1 informative role / conflicting mainroles | 375 / 308 | **15 / 14** |
| Cyanorak genes carrying curated **and** inferred roles | 13,484 | 10,369 |
| … exact agreement / same mainrole / contradiction / junk-only | 83 / 3 / 7 / 7 % | **88 / 3 / 5.4 / 3 %** |
| `Unknown` → category fill (non-Cyanorak + Cyanorak) | 1,032 + 101 | **704 + 17** |

**The 557 curated/inferred contradictions are facet choices, not errors.** They
were sampled: FtsH1 is Cyanorak *Cell division* vs TIGR *Protein degradation*;
RsgA is *Thiamine* vs *Translation factors*; polyphosphate kinase 2 gets TIGR's
*Phosphorus compounds* over Cyanorak's *Enzymes of unknown specificity*; 44
photosynthesis genes TIGR files under *One-carbon metabolism*. Both views are
kept, as **two separate edges** distinguishable by `evidence` — the graph does
not pick a winner on a multi-role protein.

**Coverage is biased — ORA over `TigrRole` needs a corrected background.**
Inferred coverage is NCBIfam-hit-bounded and equivalog-bounded: roughly
**700–860 genes per heterotroph genome** (KT2440 860, HOT1A3 700) against
**~1,800–1,976 on MED4**, where curated Cyanorak roles cover ~90% of genes. An
enrichment run whose background is "all genes of the organism" will report
role-category enrichment that is really an annotation-density artifact. Use a
**per-organism background of genes with ≥ 1 `Gene_has_tigr_role` edge**.

**`annotation_quality` does not move.** Measured: **0** genes change
`annotation_state` / `annotation_quality`. Every gene gaining a first role edge
already has an NCBIfam edge, which is already an informative bucket, and every
one of the 13,255 measured is already `informative_multi`. No new quality
bucket, no gating, no bucket-count change (still 9). This layer buys
**findability and cross-genus comparability**, explicitly not quality.

**`gene_count` on the new mainrole nodes is a subtree count.** Subrole
`gene_count` is unchanged in meaning; a mainrole's is the union over its
children. Use `direct_gene_count` when you need the non-recursive number, and do
not mix levels in one ORA without collapsing them first.

**Filter `sources` by membership, never by list equality.** A curated-only edge
is `['cyanorak']`, an inferred-only edge is `['interproscan']`, and an agreeing
edge is `['cyanorak','interproscan']`. `r.sources = ['cyanorak']` silently drops
every agreeing edge; `'cyanorak' IN r.sources` does not.

## Query cookbook

**Heterotroph genes by functional role** (the capability this layer adds):

```cypher
MATCH (g:Gene)-[r:Gene_has_tigr_role]->(t:TigrRole)
WHERE g.organism_name = 'Alteromonas macleodii HOT1A3'
  AND 'interproscan' IN r.sources
  AND t.is_uninformative IS NULL
RETURN t.name, t.level_kind, count(DISTINCT g) AS genes
ORDER BY genes DESC
```

**Roll a genome up to mainroles** (the level the two-level hierarchy makes
possible):

```cypher
MATCH (g:Gene)-[:Gene_has_tigr_role]->(:TigrRole)
      -[:Tigr_role_is_a_tigr_role*0..]->(m:TigrRole {level_kind: 'tigr_mainrole'})
WHERE g.organism_name = 'Pseudomonas putida KT2440'
RETURN m.name, count(DISTINCT g) AS genes ORDER BY genes DESC
```

**A correct ORA background** (genes with any role edge, per organism):

```cypher
MATCH (g:Gene)-[:Gene_has_tigr_role]->()
RETURN g.organism_name, count(DISTINCT g) AS background_size
ORDER BY background_size DESC
```

**Curated vs inferred disagreements** (the QC cross-check):

```cypher
MATCH (g:Gene)-[a:Gene_has_tigr_role]->(x:TigrRole),
      (g)-[b:Gene_has_tigr_role]->(y:TigrRole)
WHERE a.evidence = 'curated' AND b.evidence = 'family_inferred' AND x <> y
RETURN g.locus_tag, g.product, x.name AS cyanorak_role, y.name AS tigr_role
LIMIT 50
```

**Agreement between the two sources** (both on one edge):

```cypher
MATCH ()-[r:Gene_has_tigr_role]->()
WHERE 'cyanorak' IN r.sources AND 'interproscan' IN r.sources
RETURN count(r) AS agreeing_edges
```

**Which NCBIfam families back a role** (the bridge, read outward):

```cypher
MATCH (f:NcbifamFamily)-[:Ncbifam_family_has_tigr_role]->(t:TigrRole)
WHERE t.name STARTS WITH 'Transport and binding proteins'
RETURN t.name, collect(f.ncbifam_id)[..20] AS families, count(f) AS n
ORDER BY n DESC
```

**How much of a role is family-backed at all:**

```cypher
MATCH (t:TigrRole)
RETURN t.name, t.level_kind, t.gene_count, t.direct_gene_count,
       t.ncbifam_family_count
ORDER BY t.gene_count DESC LIMIT 25
```

## Refresh procedure

The archive is frozen (2018), so a refresh is only needed when the NCBIfam
reference itself moves or a new strain's InterProScan batch lands.

```bash
# 1. Rebuild the reference caches (writes cache/data/ncbifam/tigr_roles.json).
#    Add --refetch-raw only to re-pull the FTP archive.
bash scripts/prepare_data.sh --steps 9 --force

# 2. Re-merge gene annotations so gene_category fill + [tigr_role_inferred]
#    lines pick the roles up. Step 2 REQUIRES step 9 to have run first.
bash scripts/prepare_data.sh --steps 2 --force

# 3. Rebuild + reimport the graph (nodes, hierarchy, bridge, gene edges and
#    the post-import rollups all come from the build).
docker compose up -d
```

Verification after the rebuild:

```bash
pytest -m kg -v -k "tigr or ncbifam"
uv run python tests/kg_validity/capture_annotation_state.py --compare   # expect 0 moves
```

## Numbers measured 2026-08-29

**Shipped artifact** (`cache/data/ncbifam/tigr_roles.json`, TIGRFAMs 15.0):

| Measure | Value |
|---|---|
| Named roles | 116 (role `719`, unnamed, 102 link lines, excluded) |
| Distinct mainroles | 19 |
| Families with ≥ 1 role | 2,862 |
| … carrying 2 roles / 3 roles | 281 / 13 (294 total) |
| Role-bearing families by `family_type` | equivalog 1,870 · subfamily 523 · hypoth_equivalog 184 · equivalog_domain 118 · subfamily_domain 55 · domain 41 · exception 25 · superfamily 20 |

These counts are over the **shipped** 2,862 families. Spec §2 quotes the same
breakdown over the raw archive's 2,963 (equivalog 1,931 · subfamily 526 ·
hypoth_equivalog 185 · equivalog_domain 123 · other 198) — the difference is the
~101 families whose only role was the unnamed `719`, which the artifact drops.

**Corpus** (43 strains / 127,458 genes, spec §2):

| Measure | Value |
|---|---|
| Observed `TIGR*` NcbifamFamily nodes with an archive role | 1,720 / 2,204 (78%) |
| Observed `NF*` nodes (can never link — post-date the role system) | 2,753 |
| Distinct roles reached | 108 — 104 already `TigrRole` nodes, 4 new (`110`, `186`, `187`, `719`; `719` then dropped as unnamed) |
| Mainrole name concordance, Cyanorak vs archive, shared codes | **110 / 110** identical |
| Cyanorak-only codes not in the archive | `128`, `270`, `701`, `856` |
| Non-Cyanorak genes gaining a role edge (equivalog gate) | **13,667** (18% of 75,996) |
| Genes with > 1 inferred role (list fan-out) | ~1,751 of 24,052 equivalog-role genes; 1,429 with > 1 mainrole |
| Curated + inferred contradictions | 557 (5.4%) — facet choices, kept as two edges |
| `gene_category` `Unknown` → category | ~720 (704 non-Cyanorak + 17 Cyanorak) |
| `annotation_state` / `annotation_quality` movement | **0** |
| Non-Cyanorak genes whose TIGR category differs from the COG one | 7,837 of 19,115 (41%) — schemes carve biology differently, **not** an error rate; fill-only means zero churn |

**Graph shape** (measured full build, 2026-08-29):

| Object | Value |
|---|---|
| `TigrRole` nodes | 136 (115 subroles, 21 level-0: 19 mainroles + 2 Cyanorak-only numeric roots `tigr.role:856` / `tigr.role:270`) — was 114, all flat |
| `Tigr_role_is_a_tigr_role` | 115 |
| `Ncbifam_family_has_tigr_role` | 1,847 (0 skipped as dangling) |
| `Gene_has_tigr_role` | 49,213 → 65,542 (16,329 inferred-only + 9,640 corroborated merges), across 43 organisms (was 22 Cyanorak-curated only) |

## Dropped sources (checked, not shipped)

Both were measured because the merge discards them, and both are recorded in
`plans/backlog.md`:

- **Cyanorak `protein_domains` TIGR tokens** — 93% are already in the gene's
  InterProScan hits (1,244 extra ids over 903 genes, Pro/Syn only).
- **PGAP `inference` HMM ids** — 84% are `NF*` `domain` / `PfamEq` models that
  InterPro's NCBIfam member DB excludes **by design** (0 / 1,860 `PfamEq` and
  108 / 14,234 `domain` families ever observed across 43 genomes; the ids
  interleave numerically with observed ones, so this is not a version gap).

Neither moves any number above. Multi-source `Gene_has_ncbifam_family` is out of
scope here.

## See also

- [`interpro-multi-ontology.md`](interpro-multi-ontology.md) — the NCBIfam
  ontology this bridge starts from
- [`annotation-trust-surface.md`](annotation-trust-surface.md) — the shared
  `sources` / `evidence` ladder these edges join
- [`vocabulary-contract.md`](vocabulary-contract.md) — R2/R3/R5 house rules
  cited above
- [`ontology-level.md`](ontology-level.md) — the unified `level` convention
  `TigrRole` now satisfies
