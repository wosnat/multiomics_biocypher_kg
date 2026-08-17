# InterPro multi-ontology redesign — NCBIfam split out, faceted calls.json, provenance cleanup

**Date:** 2026-08-17
**Spec:** [`docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`](../superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md)
**Plan:** [`docs/superpowers/plans/2026-08-17-interpro-multi-ontology-redesign.md`](../superpowers/plans/2026-08-17-interpro-multi-ontology-redesign.md)
**Branch:** `interpro-multi-ontology` (21 tasks, commits `da9c8617`..`514f4a3a`)
**Supersedes:** [`interproscan-extension.md`](interproscan-extension.md) (2026-07-26 integration) and
[`interpro-two-layer.md`](interpro-two-layer.md) (2026-08-10 Layer A/B provenance) — both kept for history,
now headed "SUPERSEDED".

> **For the MCP / explorer:** this doc is the current integration contract for both `InterproEntry` and the
> new `NcbifamFamily` ontology — node/edge labels, property lists, id conventions, indexes, and the
> `calls.json` artifact shape. It supersedes the two docs above; read this one first.

## Why this shipped

InterProScan bundles NCBIfam/TIGRFAM as one of ~15 member DBs feeding the integrated `InterproEntry` layer.
NCBIfam is itself a curated, functionally-specific ontology (TIGRFAM equivalog-level calls in particular are
near-EC-grade precision) — folding it entirely into the broad, conduit-like InterPro FAMILY/DOMAIN layer threw
away that precision. This redesign splits NCBIfam into its own first-class ontology (mirroring Pfam/TCDB/CAZy),
tightens `Gene_has_interpro_entry`'s edge properties (no more synthesized cross-library `score`), gives both
`InterproEntry` and the new `NcbifamFamily` a short curated `description`, adds naming-recovery for genes whose
`product` was previously uninformative, and closes a bucket-accounting gap in `annotation_quality` where
interpro/ncbifam-only genes were being scored as `no_evidence` instead of `catch_all_only`.

## New node label: `NcbifamFamily`

- **Count:** 4,957 nodes (2,204 `TIGR*` + 2,753 `NF*`) — observed-only, flat (no ancestor walk; the source TSV
  carries no parent/child column, unlike InterPro's ParentChildTreeFile).
- **ID prefix:** `ncbifam_<accession>` (underscore form — `ncbifam` is not a registered bioregistry prefix,
  same convention as `psortb_*` / `signalp_*`).
- **Properties (adapter-emitted):** `name`, `ncbifam_id` (`TIGR*` | `NF*`), `family_type` (NCBIfam's own
  classification — **equivalog-dominant** in the observed set, i.e. most calls are precise
  "this exact enzymatic/functional role" assertions, not broad homology families),
  `gene_symbol` (sparse), `description` (sparse — short curated text from `ncbifam_reference.json`).
  Retired accessions (present on a gene but absent from the reference — cross-strain observed via calls.json)
  still emit a minimal node (`family_type` omitted) so gene edges never dangle.
- **Computed (post-import):** `gene_count`, `organism_count`.

## New edge types

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_ncbifam_family` | Gene → NcbifamFamily | 67,459 | `start:int`, `end:int`, `evalue:float`, `score:float`, `match_count:int` — **keeps both `evalue` AND `score`**, unlike `Gene_has_interpro_entry` (see below), because NCBIfam is a single homogeneous HMMER-scale library, not a cross-library rollup |
| `Ncbifam_family_in_interpro_entry` | NcbifamFamily → InterproEntry | 2,630 | — (bridge; overlap link between the two ontologies, not a merge) |

`Gene_has_interpro_entry` stays at **397,342 edges** (unchanged — NCBIfam calls were already folding into
InterPro entries before this redesign; splitting them into a parallel ontology doesn't remove the InterPro-side
edge, it adds a second, more specific one). 47,324 genes carry ≥1 `ncbifam` annotation (46,159 of those pass the
informativeness filter — see the bucket section below).

## Edge property change: `Gene_has_interpro_entry` drops `score`, gains `evalue_library`

**Before:** the edge carried both `evalue` (best/min across all contributing library rows) and a synthesized
`score` (best/max), silently blending values from heterogeneous scoring scales (Pfam bit scores, PROSITE
pattern scores, HAMAP profile scores — not comparable to each other).

**After ("count, don't combine"):** the edge keeps `evalue` (min across contributing rows) and adds
`evalue_library` (str — which library reported that minimum), and drops `score` entirely. Cross-library
corroboration is read from `size(libraries)` (independent model families agreeing) vs. `match_count` (total
rows, including within-library repeats) — never from a blended number. `Gene_has_ncbifam_family` is exempt from
this change because it is a single-library edge (see above).

## `description` field (new on both `InterproEntry` and `NcbifamFamily`)

Both node types now carry a sparse curated `description` string (short, 400-char-capped abstract for InterPro;
short curated text for NCBIfam). Both full-text indexes (`interproEntryFullText`, new `ncbifamFamilyFullText`)
include `description`.

**InterPro descriptions shipped in OBSERVED-ONLY mode, not full-corpus.** The full description set (all
54,068 abstract-bearing InterPro entries) serializes to ~27.1 MB, over the reference cache's ~25 MB size gate.
The documented fallback fired automatically: descriptions are kept only for the ~13K InterPro accessions
actually observed on some strain's `calls.json` **plus their is-a ancestors** (12,955 of 12,999 kept nodes carry
a non-empty description; 44 have no InterPro abstract at all). This makes the InterPro reference **corpus-coupled**
in a new way it wasn't before (name/type/hierarchy/GO/pathway/EC/CAZy fields were already snapshot-not-live, but
this is the first field gated on which strains have been scanned): **onboarding a new strain with genuinely novel
domain content requires re-running `prepare_data.sh --steps 9 --force` after that strain's InterProScan batch
lands, or the new strain's InterPro nodes will fail-soft with a missing `description`** (every other field stays
correct — this is additive-only degradation, not a data-loss risk). See Data Locations in `CLAUDE.md`.

NCBIfam's reference (`ncbifam_reference.json`, 38,394 entries) did not need an observed-only fallback — the raw
descriptions are short enough that the full set fits under the size gate comfortably.

## `calls.json` faceted format (supersedes the old flat `matches[]`/`interpro_entries[]` shape)

`<strain>.interproscan.calls.json` is now a WP_-keyed dict of per-protein records with three sparse facets:

```jsonc
{
  "<WP_accession>": {
    "md5": "<str>", "match_count": 0,
    "libraries": {"<LIBRARY>": [{"accession","name","ipr","start","end","evalue","score"}]},
    "interpro_entries": {"<IPR*>": {"type","libraries","match_count","start","end","evalue","evalue_library"}},
    "go_terms": {"<GO:NNNNNNN>": ["<IPR*>", "..."]}
  }
}
```

- `libraries` is sparse per member DB; each row's `ipr` is `null` if that signature isn't integrated into an
  InterPro entry (common for PANTHER/PIRSF/MobiDBLite and some NCBIfam rows).
- `interpro_entries` carries **no cross-library `score`** — only min `evalue` + `evalue_library` (the "count,
  don't combine" rule above, now visible at the artifact level, not just the edge level).
- `go_terms` maps a GO id to the sorted list of `IPR*` entries that donated it — attribution, not a flat list.
- **No `pathways` field anywhere** — InterPro carries no KEGG pathway xrefs at all (verified: `db="KEGG"` occurs
  0 times in the full release), and the Reactome/MetaCyc xrefs that do exist are entry-level and fully derivable
  from the central `interpro_reference.json`; duplicating them per-strain per-protein was pure bloat.
- **`<strain>.interproscan.entry_xrefs.json` no longer exists.** GO provenance moved inline into `go_terms`;
  entry-level lookups (names, types, hierarchy, GO, pathways, EC, CAZy) live only in the two central reference
  caches (`interpro_reference.json`, `ncbifam_reference.json`).
- New `--normalize` mode on `/interproscan-run`: re-parses a strain's already-cached `<strain>.interproscan.raw.json`
  through the current parser and rewrites `calls.json` in place — no Docker re-run needed after a parser change.
  Used to migrate all 42 strains to the faceted format in this redesign without re-scanning.

Full artifact contract: `.claude/skills/interproscan-run/SKILL.md` (Output Schema section).

## Naming recovery

Genes whose `product` field was previously uninformative (or only a raw domain name) now get better names from
two new sources, both gated to avoid restating what's already known:

- `[hamap] <description>` — from HAMAP entries' descriptions, skipped when the text case-insensitively equals
  `gene["product"]` (i.e. it would just echo the existing product string).
- `[ncbifam] <product name>` — from the NCBIfam reference's product-name field for the gene's `ncbifam_ids`,
  same echo-skip rule. On MED4 alone, 325 of 1,079 considered `[ncbifam]` name tokens were skipped as echoes
  (out of the full corpus this ratio holds broadly) — the gate is doing real work, not a no-op.
- `gene_name` fill-if-empty: when a gene has no `gene_name` yet, the first sorted `gene_symbol` among its
  `ncbifam_ids` reference entries fills it (`gene_name_source = "ncbifam"`), never overwriting an existing name.
  Verified examples on MED4: `yvcK`, `trhO`.

Both naming-recovery strings flow through the existing exact-string `afd_set` dedup on
`alternate_functional_descriptions`, so a name recovered this way never duplicates an identical string already
contributed by another source.

## `is_uninformative` rules (new, both ontologies)

- **InterProEntry:** name-pattern rule flags entries whose `name` matches
  `^Protein of unknown function.*` / `^Domain of unknown function.*` / `^Uncharacteri[sz]ed protein family.*`
  (case-sensitive, exact InterPro naming convention).
- **NcbifamFamily:** a typed rule (`family_type IN ['hypoth_equivalog', 'hypoth_equivalog_domain']` — NCBIfam's
  own "this looks like an equivalog but function isn't confirmed" categories) plus a name-pattern fallback
  (`(?i).*hypothetical.*`, `(?i).*uncharacterized.*`, `.*DUF\d.*`). 114 NcbifamFamily nodes flagged this way.

Both rules exist so ORA / routing signals can down-weight "we found a domain but don't know what it does" hits
without hand-filtering by name at query time — see the 9-bucket `annotation_quality` section, which uses these
flags via the same `t.is_uninformative IS NULL` predicate pattern established for TCDB/BRITE/CAZy.

## `annotation_quality` bucket resolution: 9 buckets, `has_any_edge` fix

**New 9th bucket: `ncbifam`.** The live source-bucket list (see `CLAUDE.md` "Source bucket maintenance") is now
`go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter`, `cazy`, `ncbifam` — **NOT** `interpro`. InterPro
stays deliberately out of the bucket count: it is a conduit/routing signal (folded into `annotation_types` and
`Gene.interpro_entry_count`, but not into `informative_annotation_types` / `annotation_quality`) because it is
largely redundant with the `pfam`/`go`/`ec` buckets it feeds — counting it again would inflate quality scores
without adding independent information. NCBIfam, by contrast, is curated and functionally specific enough
(equivalog-dominant, per above) to earn its own bucket.

**`has_any_edge` fix.** The `annotation_state = 'no_evidence'` classification checks whether a gene has *any*
ontology edge at all, via an explicit relationship-type list. That list was missing
`Gene_has_interpro_entry`/`Gene_has_ncbifam_family` — so a gene whose *only* annotation was an InterPro or
NCBIfam call was being scored `no_evidence` even though it visibly has edges in the graph. Fixed by extending the
list. Both fixes landed in the same post-import Cypher change (`scripts/post-import.cypher` +
`scripts/post-import.sh`, kept byte-identical in logic per the usual invariant).

**Measured movement on the live rebuild** (whole-KG `annotation_state`):

| state | before | after | delta | driver |
|---|---:|---:|---:|---|
| `no_evidence` | 12,481 | 12,061 | **−420** | `has_any_edge` fix picks up interpro/ncbifam-only genes |
| `catch_all_only` | 5,752 | 5,988 | **+236** | net of 420 movers in − ~184 climbers out (climbers moved on to `informative_single`/`informative_multi` via the new `ncbifam` bucket) |

Both deltas independently matched the spec's predictions (~367–420 movers, ~185 climbers) almost exactly — see
`.superpowers/sdd/2026-08-17-interpro-multi-ontology-redesign/task-19-report.md` for the full verification
(including the caveat that the pre-rebuild baseline was supplied as prior context, not captured live in the same
session — internally consistent with two independent predictions, treated as reliable).

`annotation_types` also gains `'ncbifam'` alongside the pre-existing `'interpro'` (interpro was already folded
there from the original 2026-07-26 integration; ncbifam is new).

## The `IPR000362` / `IPR009049` multi-EC-gate nuance

The InterPro→gene EC propagation gate (from the earlier `interpro-two-layer.md` integration, unchanged by this
redesign but re-verified during Task 18's spot checks) only stamps an EC onto a gene from a **single-EC FAMILY**
entry — a **multi-EC FAMILY** entry (like `IPR000362`, argininosuccinate lyase family, which carries 5 ECs:
`4.2.1.2`/`4.3.1.1`/`4.3.2.1`/`4.3.2.2`/`5.5.1.2`) contributes none of its ECs to `ec_numbers` (they go to the
Layer-A router edge instead — see `interpro-two-layer.md`).

The nuance: on MED4's `PMM0012` (argininosuccinate lyase, `WP_011131650.1`), `ec_numbers` legitimately **does**
contain `4.3.2.1` with `interpro` in its source list — but not because `IPR000362` leaked it. A *second,
independent* InterPro entry on the same gene, `IPR009049`, is itself **single-EC** (its sole EC is `4.3.2.1`),
so it passes the single-EC gate on its own merits. The multi-EC gate on `IPR000362` still held perfectly — none
of the other 4 ECs it carries (`4.2.1.2`, `4.3.1.1`, `4.3.2.2`, `5.5.1.2`) leaked in. This is worth documenting
because it's easy to misread as a gate failure from the outside: **check per-entry, not per-gene**, when
verifying the multi-EC refuse rule — a gene can legitimately carry an EC from InterPro while also carrying a
different, un-donated multi-EC entry. Full detail:
`.superpowers/sdd/2026-08-17-interpro-multi-ontology-redesign/task-18-report.md`.

## `prepare_data.sh` step reordering

**Step 9** now builds *two* central reference caches (`build_interpro_reference.py` then
`build_ncbifam_reference.py`, both logged into a shared `logs/prepare_data_step9.log` via append mode), not
just the InterPro one.

**Default `STEPS` order changed to `"0 9 1 2 3 4 5 6 7 8"`** (step 9 moved before 1/2, not run last as its
number would suggest) — this is dependency-ordering, not a renumber: step 1 (`build_protein_annotations.py`) and
especially step 2 (`build_gene_annotations.py`, which lazily loads both `interpro_reference.json` and
`ncbifam_reference.json` to enrich the merged JSON) need both reference caches to exist first. Step *numbers*
are unchanged (step 9 is still called "step 9" everywhere — flags, logs, CLI); only the default execution order
within a full `prepare_data.sh` run moved. Running `--steps 9` alone still works exactly as before.

## Consistency guards (new tests)

- **`tests/test_interproscan_consistency.py`** — asserts exact equality between what a strain's
  `<strain>.interproscan.calls.json` says (`interpro_entries` keys, `ncbifam_ids` from the NCBIFAM library facet)
  and what `gene_annotations_merged.json` actually merged in, per strain. Was already partially in place;
  this redesign extended it to also check `ncbifam_ids` and confirmed 42/42 strains pass post-re-merge (Task 18).
- **`tests/test_annotation_quality_buckets.py`** (**NEW file**) — a fast, no-Neo4j static gate that parses
  `scripts/post-import.cypher` / `scripts/post-import.sh` and asserts the declared `SOURCE_BUCKETS` comment and
  the actual `has_<bucket>` terms summed into `informative_count` both equal the 9-bucket set, and that the two
  files agree with each other. **CLAUDE.md's "Source bucket maintenance" section previously referenced "the
  bucket-count test" as an existing thing — it did not exist before this task** (verified by search across
  `tests/` and git history); this file is that test, created for real.

## Process notes

- **`docker compose up <service>` re-cascades the full chain every invocation.** This compose file has no
  cross-invocation "already completed" memory — calling `docker compose up import` (say) after a prior, separate
  `docker compose up build` that already exited 0 re-runs `build` again before `import`, even though nothing
  changed. Observed cost: ~5 min per redundant cascade against warm caches during this branch's verification
  (three full cascades run back-to-back). Not a correctness issue, just a planning note for anyone scripting
  repeated `docker compose up <service>` calls expecting partial-cascade idempotency — it doesn't have any.
- **Known Issues update:** the CLAUDE.md-documented `test_no_orphan_proteins` / `test_no_orphan_proteins_without_gene`
  failures **did not reproduce** on the 2026-08-17 full rebuild verified for this redesign (both passed clean,
  `pytest -m kg` → 1098 passed / 4 skipped / 0 failed). Not investigated further as part of this work — the
  Known Issues entry is softened to "did not reproduce on the 2026-08-17 rebuild; verify on next rebuild before
  closing" rather than removed, since a single clean run isn't proof the underlying gap is gone.

## Full actuals reference

| Metric | Value |
|---|---|
| `NcbifamFamily` nodes | 4,957 (2,204 `TIGR*`, 2,753 `NF*`) |
| `Gene_has_ncbifam_family` edges | 67,459 |
| `Ncbifam_family_in_interpro_entry` edges | 2,630 |
| `Gene_has_interpro_entry` edges | 397,342 (unchanged) |
| Genes with ≥1 ncbifam annotation | 47,324 (46,159 pass informativeness) |
| `annotation_state` bucket count | 9 (was 8) |
| `no_evidence` movement | 12,481 → 12,061 (−420) |
| `catch_all_only` movement | 5,752 → 5,988 (+236 net) |
| Import | exit 0, zero skipped relationships |
| `/omics-edge-snapshot` before/after | zero regressions (expression, metabolism, DerivedMetric edges all +0) |
| `pytest -m kg` | 1098 passed, 4 skipped, 0 failed |

Full task-by-task evidence trail: `.superpowers/sdd/2026-08-17-interpro-multi-ontology-redesign/` (21 task
briefs + reports, `progress.md` ledger, `task-19-report.md` for the live-rebuild verification this doc's
actuals are drawn from).
