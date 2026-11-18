# DataSource citations — design (short)

**Date:** 2026-11-18 · **Size:** small change · **Revised** after explorer review (same day)

## Intent
A new explorer MCP tool should answer "what do I cite for this KG?" by returning
`{data_source_id: [citations]}` for every tool and database the KG was built from, plus the
build framework. The KG side supplies the data; the MCP tool (multiomics_explorer) is a
separate follow-up — nothing in the explorer queries DataSource yet, so the follow-up is a
new tool + a schema-baseline refresh.

## Goals / non-goals
- Goal: every DataSource node carries the citation(s) that its provenance unit needs.
- Goal: one query serves the MCP tool:
  `MATCH (d:DataSource) RETURN d.id, d.name, d.citations, d.dois, d.url, d.version, d.license`
  `+ MATCH (s:Schema_info) RETURN s.build_citations`.
- Non-goal: per-query scoping (which sources back which labels/edge types) — see Approaches.
- Non-goal: the paper DOIs behind the expression/assay data — already `Publication` nodes.
- Anti-goal: wrong or invented citations. Every entry is taken from the tool's/database's own
  "how to cite" page.
- Anti-goal: muddying DataSource. It stays a **provenance target** (R2: edge `sources` values
  join to it); build software is not a DataSource.

## Approaches
- **A — whole-KG list on DataSource (chosen).** Citation lists on the existing nodes, a few
  reference databases added; the MCP returns the full map.
- B — also store `applies_to: str[]` (labels/edge types each source backs) so the MCP
  returns only what an answer touched. Rejected: hand-curated mapping over ~90 edge types +
  ~35 labels, kept in sync with every new source; user: "only if really easy".
- C — graph edges DataSource → label. Rejected: edges to labels are not modelable.

## Design
1. **New committed `config/data_source_citations.yaml`** — the single source of per-source
   metadata, keyed by DataSource id: `name`, `description`, `citations` (list of formatted
   references), `dois` (list, bare DOIs, parallel to `citations`, `""` where none), `url`,
   `license` (optional), and `version` (literal) or `version_from` (named reader). Plus a
   top-level `build_citations` list (BioCypher, Lobentanzer 2023). The hard-coded
   `_name_for` / `_description_for` dicts in `data_source_adapter.py` are **deleted**
   (their text is already stale, e.g. "InterProScan 5").
2. **One provenance unit, several citations** (`citations: str[]`, not `citation: str`).
   Examples (the builder confirms each against the upstream citation page):
   - `eggnog`: eggNOG-mapper v2 + eggNOG 5 + COG (Galperin; COG categories arrive via eggNOG)
   - `tcdb_diamond`: TCDB + DIAMOND (the TCDB ontology/substrate data cites under it too)
   - `merops_diamond`: MEROPS + DIAMOND
   - `interproscan`: InterProScan + InterPro + NCBIfam + TIGRFAMs (roles archive)
   - `kegg`: KEGG's requested papers; **description widened** to KO / pathways / BRITE /
     COMPOUND (it currently says COMPOUND only)
   So DIAMOND, InterPro, NCBIfam and TIGRFAMs get **no nodes of their own**.
3. **New reference sources** (in `reference_sources` of `gene_annotations_config.yaml`,
   `scope: global`, `provenance: download`) — databases whose own downloads become KG
   content and that no existing node covers: `go`, `pfam`, `cazy`, `expasy_enzyme`
   (`cache/data/ec/ec_data.json`), `metanetx`, `ncbi_taxonomy` (`cache/data/taxonomy/`).
   PubChem / Crossref (referenced in `metabolite_identity/pipeline.py`): builder checks
   whether their data lands in the graph; include only if it does.
4. **Build framework on `Schema_info`, not DataSource.** `Schema_info.build_citations: str[]`,
   written by `create_knowledge_graph.py` to a `build_citations.json` sidecar next to the
   BioCypher output and stamped by post-import — the existing `chebi_release.txt` pattern
   (both `post-import.sh` and `post-import.cypher`; declared in `schema_config.yaml`).
   Neo4j gets nothing (no citation to give).
5. **Schema** (`data source` in `schema_config.yaml`): add `citations: str[]`, `dois: str[]`,
   `url: str`, `license: str` (sparse); fix the `scope` comment to include `global`.
6. **Versions** — `version_from` readers where a committed cache records one:
   `chebi` (existing), `go` (`go_terms/go_release.txt`). TIGRFAMs has no node, so its
   release ("TIGRFAMs 15.0", from `ncbifam/tigr_roles.json`) is named in its citation text
   instead. **Literals** in the yaml for the rest: `interproscan` "111.0"
   (InterPro release), `merops_diamond` "12.5", `eggnog` (eggNOG DB version), others as known
   or `""`. Neither `interpro_reference.json` nor `merops_reference.json` records a release,
   so literals are unavoidable there. The unit test pins the InterPro literal to
   `INTERPRO_RELEASE` in `.claude/skills/interproscan-run/run_interproscan.py` (parsed, not
   imported), so a rescan that bumps the pin fails the test until the yaml follows.
7. **Adapter** reads the yaml and emits the new properties; `_clean_str` applies (it maps
   `'` → `^` and **deletes** `|`), so curated text avoids both.
8. **Tests** (`tests/test_data_source_adapter.py`):
   - update `test_emits_eleven_nodes` (exact id set → new set) and the
     `nodes["kegg"]["version"] == ""` assertion;
   - every DataSource id has a yaml entry with non-empty `name`, `citations`, `url`, and
     `len(dois) == len(citations)`; no yaml entry without a DataSource;
   - no `'` or `|` in curated text; InterPro literal == `INTERPRO_RELEASE`.
   (The R2 vocabulary tests are unaffected — adding ids only enlarges `data_source_ids`.)
9. Docs: CLAUDE.md DataSource bullet (count, new properties, Schema_info.build_citations),
   CHANGELOG `### Added`.

## Left to the builder
Exact citation text and DOIs per source; the eggNOG DB version literal; whether PubChem /
Crossref qualify (step 3).
