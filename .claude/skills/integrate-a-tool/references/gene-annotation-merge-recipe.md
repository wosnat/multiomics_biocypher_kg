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

Line anchors verified 2026-05-26 — they drift; re-`grep` the symbol.

## A. `config/gene_annotations_config.yaml`

Two paste targets. See `assets/gene_annotations_source_snippet.yaml` for the full stub.

1. **Under `sources:`** ([~L24](../../../../config/gene_annotations_config.yaml#L24)) — add a `<tool>:` block:
   - `type:` — must match a `load_<tool>()` you add below.
   - `path_pattern: "{data_dir}/<tool>/{strain_name}.<tool>.calls.json"`
   - `join_to: protein_id` — calls.json keyed by RefSeq WP_ accession (psortb, signalp, eggnog all do).
   - `description:`
   - `logical_sources:` — list with `id` / `scope` / `provenance`. **Valid values in
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

3. **Thread the new `tool` dict through the resolvers.** Every `_resolve_*` method
   (`_resolve_passthrough` ~L455, `_resolve_passthrough_list` ~L466, `_resolve_union`,
   `_resolve_single`, …) calls `self._get_raw(fconf, gm, eg, up)`. `grep -n "_get_raw(" `
   and `grep -n "def _resolve"` to find them all; add the `tool` positional everywhere
   so it reaches `_get_raw`.

4. **`build_wide`** ([~L613](../../../../multiomics_kg/download/build_gene_annotations.py#L613)) —
   signature `(self, gm, eg, up)`; the body source-prefixes each dict
   (`wide[f"eggnog_{k}"] = v` at ~L626). Add `tool` param + a
   `for k, v in tool.items(): wide[f"<tool>_{k}"] = v` block.

5. **`build_merged`** ([~L634](../../../../multiomics_kg/download/build_gene_annotations.py#L634)) —
   signature `(self, gm, eg, up, organism_name=None)`. Add the `tool` param; it calls
   the resolvers (point 3) which call `_get_raw`.

6. **`process_strain`** ([~L866](../../../../multiomics_kg/download/build_gene_annotations.py#L866)) —
   loads each source (`eg_data = load_eggnog(data_dir, strain_name)` ~L900) and passes
   them into `build_wide`/`build_merged`. Add `tool_data = load_<tool>(data_dir, strain_name)`
   and thread it into both calls.

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

## Verify (before forking to a track)

```bash
bash scripts/prepare_data.sh --steps 2 --strains MED4 --force
```
Then confirm:
- the new field landed in `cache/data/<org>/genomes/MED4/gene_annotations_merged.json`
  (`jq '.[] | select(.<field>) | .<field>' …`);
- `contributing_sources` lists `<tool>` on a gene that has the annotation;
- a `DataSource` node will emit — `DataSourceAdapter` reads the `logical_sources` you added.
