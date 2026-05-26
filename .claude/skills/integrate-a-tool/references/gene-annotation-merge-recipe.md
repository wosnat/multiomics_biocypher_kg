# Gene-annotation merge recipe — exact edit points (Step 2)

Both tracks flow through the **same front door**: the tool's calls.json becomes a
new *source* in the gene-annotation merge, producing a field in
`gene_annotations_merged.json`. Merging first is **mandatory** — it is what makes
`contributing_sources` and the auto-generated `DataSource` node pick the tool up
for free ([data_source_adapter.py](../../../../multiomics_kg/adapters/data_source_adapter.py)
derives everything from the `logical_sources` block, and **fails loudly** if it's missing).

> The merge is **hand-wired**, not config-driven — `_get_raw` dispatches on the
> source *name* and the source dicts are threaded as **positional args** (`gm`, `eg`,
> `up`). Adding a source means editing every signature + call site in that chain.
> Generalizing the builder is an explicit non-goal. Follow the
> hand-wired pattern; don't refactor it.

> ⚠️ **Add the new source arg as an OPTIONAL default — `tool: dict | None = None`
> (coerce to `{}` in the body) — NOT a required positional like `eg`/`up`.**
> `tests/test_build_gene_annotations.py` calls `build_wide(...)`, `build_merged(...)` and
> `_resolve_union(...)` **positionally ~120 times** (117 `build_*` + 3 `_resolve_union`); a required
> 4th positional detonates the whole suite with zero error message pointing at the cause. `eg`/`up`
> are required positionals only because they predate that test suite — do **not** copy that detail.
> (Internal resolvers the tests never call directly can take a required positional — see point 3.)

Line anchors verified 2026-05-26 — they drift; re-`grep` the symbol.

## A. `config/gene_annotations_config.yaml`

Two paste targets. See `assets/gene_annotations_source_snippet.yaml` for the full stub.

1. **Under `sources:`** ([~L24](../../../../config/gene_annotations_config.yaml#L24)) — add a `<tool>:` block.
   **Only `logical_sources` is read by code** (by `DataSourceAdapter`). `type` / `path_pattern` /
   `join_to` are **documentary only** — no code dispatches on them (`grep join_to` → YAML only; the
   loader + join are hand-wired in `process_strain`, section B). Write them for the reader, but don't
   expect them to wire anything:
   - `type:` / `path_pattern:` / `join_to:` — documentary; describe the loader + join you hand-wire below.
   - `description:`
   - `logical_sources:` — list with `id` / `scope` / `provenance` — **the part that actually does
     something** (emits a `DataSource` node + a `contributing_sources` value). **Valid values in
     the live config:** `scope: gene_level | organism_restricted`; `provenance:
     tool_run | download`. Use `scope: gene_level`, `provenance: tool_run` for a tool.
     (`organism_restricted` additionally takes `applies_to_organisms: [...]`.)

2. **Under `fields:`** ([~L78](../../../../config/gene_annotations_config.yaml#L78)) — add field rule(s).
   Valid rule `type:` values in use: `passthrough`, `passthrough_list`, `union`,
   `single`, `integer`, `float`, `computed`, `extract_first_match`. Shapes:
   - 1:1 scalar (psortb/signalp): two `passthrough` rules (`<tool>_call`, `<tool>_score`).
   - multi-call scored: one structured list-of-dicts via `passthrough` of the loader's list.
   - bare multi-call (cazy-style `str[]`): `union` with a delimiter (model `cazy_ids` at [~L540](../../../../config/gene_annotations_config.yaml#L540)).

## B. `multiomics_kg/download/build_gene_annotations.py`

Six edit points. Model every one on **eggNOG** — `grep -n eggnog` traces the whole thread.

1. **`load_<tool>(data_dir, strain_name) -> dict[str, dict]`** — model
   [`load_eggnog`](../../../../multiomics_kg/download/build_gene_annotations.py#L334)
   ([~L334](../../../../multiomics_kg/download/build_gene_annotations.py#L334)). Returns
   `{join_key: {field: value, …}}`. Keep parsing pure where possible.

2. **`_get_raw`** ([~L411](../../../../multiomics_kg/download/build_gene_annotations.py#L411)) —
   signature today is `(self, src_cfg, gm, eg, up)`; the dispatch block (~L419–428) is:
   ```python
   if source == "gene_mapping":   raw = gm.get(field)
   elif source == "eggnog":       raw = eg.get(field)
   elif source == "uniprot":      raw = up.get(field)
   else:                          return None
   ```
   Add `tool` to the signature **and** `elif source == "<tool>": raw = tool.get(field)`.

3. **Thread the new `tool` dict through ALL SIX resolvers + a second dispatch site.**
   `grep -n "def _resolve"` — there are **six**, not three: `_resolve_passthrough`,
   `_resolve_passthrough_list`, `_resolve_single`, `_resolve_union`, `_resolve_integer`,
   `_resolve_float`. Each calls `self._get_raw(fconf, gm, eg, up)`; thread `tool` through to
   `_get_raw`. Two snags:
   - `_resolve_single`'s signature is `(fconf, gm, eg, up, source_tracking, locus_tag="")` — the new
     arg slots **before** the non-defaulted `source_tracking`, so it's a required positional there
     (fine — tests don't call `_resolve_single` directly).
   - `_resolve_union` **is** called positionally by the test suite — give the new arg a default there
     (and on `build_wide`/`build_merged`), per the ⚠️ callout above.
   - **Second dispatch site, outside the resolvers:**
     [`extract_first_match_in_sources`](../../../../multiomics_kg/download/utils/annotation_helpers.py#L39)
     in `annotation_helpers.py` (called from `build_merged` for `type: extract_first_match` fields) has
     its own `if source == "gene_mapping"/"eggnog"/"uniprot"` block (~L58–62). Add a `tool` param +
     branch there **iff** any of your fields use `extract_first_match` — otherwise that source silently
     won't resolve, with no error.

4. **`build_wide`** ([~L613](../../../../multiomics_kg/download/build_gene_annotations.py#L613)) —
   signature `(self, gm, eg, up)`; the body source-prefixes each dict
   (`wide[f"eggnog_{k}"] = v` at ~L626). Add `tool` param + a
   `for k, v in tool.items(): wide[f"<tool>_{k}"] = v` block.

5. **`build_merged`** ([~L634](../../../../multiomics_kg/download/build_gene_annotations.py#L634)) —
   signature `(self, gm, eg, up, organism_name=None)`. Add the `tool` param; it calls
   the resolvers (point 3) which call `_get_raw`.

6. **`process_strain`** ([~L866](../../../../multiomics_kg/download/build_gene_annotations.py#L866)) —
   **two edits, and the second is the actual join.**
   - (a) Load the source: `tool_data = load_<tool>(data_dir, strain_name)` (model `eg_data = load_eggnog(...)` ~L900).
   - (b) **The row-level join is what makes `join_to` real:** inside the per-gene `for locus_tag, gm_row …`
     loop, fetch this gene's record by the join key — `tool_row = tool_data.get(protein_id, {})` — and pass
     `tool_row` into the `build_wide`/`build_merged` calls (mirror how the `eg`/`up` rows are fetched +
     passed). Without this `.get(<join_key>)` lookup the field never lands, even with everything else wired.

7. **`_compute_contributing_sources`** ([~L292](../../../../multiomics_kg/download/build_gene_annotations.py#L292)) —
   add a presence branch (model the eggnog block at ~L307–309: a field-presence OR
   `_has_source_label(gene, "<tool>")` test → `sources.add("<tool>")`).

## C. `multiomics_kg/download/utils/annotation_transforms.py` (only if reshaping)

If a raw value needs prefix/split/normalize/parse, add `_tx_<name>` and register it in
`_TRANSFORMS` ([~L263](../../../../multiomics_kg/download/utils/annotation_transforms.py#L263)),
reference it via `transform: <name>` in the field rule, and unit-test it in
`tests/test_annotation_transforms.py`. Reusable today: `_tx_first_token_space`
([~L27](../../../../multiomics_kg/download/utils/annotation_transforms.py#L27)) — WP_
accession from a FASTA-header key (signalp needs this).

## Two consequences of adding a source — fix these or tests break

- **The `DataSource` node gets a junk name unless you edit a file the rest of this recipe never
  touches.** [`data_source_adapter.py`](../../../../multiomics_kg/adapters/data_source_adapter.py)
  has hardcoded `_name_for`/`_description_for` dicts (~L110/L120); an unlisted source falls back to
  `source_id.title()` (→ "Psortb") with an empty description. Add your tool to **both** dicts.
- **`tests/test_data_source_adapter.py` asserts an exact node count** (`test_emits_four_nodes`).
  A new `logical_source` emits one more node → that test fails. Update the count and add a
  provenance assertion for your source.

## Verify (before forking to a track)

```bash
bash scripts/prepare_data.sh --steps 2 --strains MED4 --force
```
Then confirm:
- the new field landed in `cache/data/<org>/genomes/MED4/gene_annotations_merged.json`
  (`jq '.[] | select(.<field>) | .<field>' …`);
- `contributing_sources` lists `<tool>` on a gene that has the annotation;
- a `DataSource` node emits with a *real* name (not "Psortb") — confirms the `_name_for` edit;
- `pytest tests/test_build_gene_annotations.py tests/test_data_source_adapter.py -q` is green
  (catches the positional-arg + node-count breakage before you move on).
