# InterPro integration review — KG correctness

**Date:** 2026-08-13
**Branch reviewed:** `tcdb-rollup-depth-and-uniprot-sanitisation` (InterPro work merged via `e4a3f717`, `88b20f24`, `b30356cd`, `6c962f50`, `4a609b77`)
**Build audited:** `biocypher-out/20260812165914/` + live Docker graph (`deploy` container)
**Focus:** KG correctness; not overrepresenting uncertain evidence; plus standard code-review items.

Every finding below was verified against the raw InterPro release, the shipped CSVs, and/or the
live Neo4j graph. Repro commands are included so each can be re-checked independently.

---

## Summary

| # | Severity | Finding | Blast radius |
|---|---|---|---|
| 1 | **Critical** | `ambiguous` boolean destroyed at import — all Layer A router edges land `false` | 6,854 EC + 122 CAZy edges |
| 2 | **Critical** | Narrative `<abstract>` prose harvested as curated EC/CAZy xrefs | 979 gene→EC + 175 gene→CAZy edges |
| 3 | Medium | `evidence='curated'` is the default for anything InterPro didn't touch | 429,287 GO + 11,537 EC edges |
| 4 | Medium | Layer A over-prunes exactly the edges it exists for | 757 candidates dropped (755 ambiguous) |
| 5 | Low | Pfam evidence ordering — `signature` where docstring promises `curated` | 59,363 edges (label only) |
| 6 | Low | `_normalize_ec_token` missing `_VALID_EC_RE` (parity claim false) | masked by EC-node prune |
| 7 | Low | `_CURATED_SOURCES` dead code in `annotation_provenance.py` | — |
| 8 | Low | `evalue=min` / `score=max` folded across incomparable member DBs | 397,342 edges |
| 9 | Low | Latent coupling: edge targets from calls.json, node set from merged JSON | risk only, 0 actual |
| 10 | Low | Doc drift: `Pfam_in_interpro_entry` count | ~5,236 documented vs 5,972 actual |

---

## 1. CRITICAL — `ambiguous` flag silently destroyed at import

**Where:** `multiomics_kg/adapters/interpro_adapter.py:399` and `:417`

Both Layer A edge types emit a Python `bool`. BioCypher writes `True`/`False` (Python
capitalisation); `neo4j-admin import` parses only lowercase `true`, so **every `True` becomes
`false`**.

| Edge type | CSV | Live graph |
|---|---|---|
| `Interpro_entry_related_to_ec_number` | 4,970 True / 1,884 False | **6,854 all FALSE** |
| `Interpro_entry_related_to_cazy_family` | 84 True / 38 False | **122 all FALSE** |

Two edges confirmed individually — `IPR000023-related_ec-2.7.1.11` and
`IPR000043-related_ec-3.13.2.1` are `True` in the CSV and `FALSE` in Neo4j.

### Why this matters

`ambiguous` is the **single safeguard** marking a Layer A edge as a one-of-several candidate
rather than a claim. With it uniformly false, every router edge in the graph reads as an
unambiguous assertion. This is the highest-impact "overrepresenting uncertainty" defect found.

### This lesson was already learned

CLAUDE.md records it verbatim for TCDB `substrate_depth`:

> Categorical **string**, not bool — BioCypher mishandles boolean properties.

It was not carried over to Layer A.

### Scope check

`ambiguous` is the **only** CSV-imported boolean in the entire build. `is_promiscuous`
(InterproEntry + TcdbFamily) and the three `Gene_has_tcdb_family` booleans are all written
post-import via Cypher and are correct (verified: 22 InterproEntry nodes carry `TRUE`).

### Proposed fix

Convert to a categorical string, mirroring `substrate_depth`:

```python
# interpro_adapter.py:399 / :417
{"ambiguous": "ambiguous" if amb else "specific", "source_db": "interpro.xml"}
```

Then update `config/schema_config.yaml:1532` and `:1543` (`ambiguous: bool` → `str`), plus any
MCP/explorer consumer that reads the flag.

### Why it escaped

`tests/test_interpro_adapter.py:137,139` assert `is False` / `is True` on the **adapter's Python
dict**, which passes. Nothing covers the CSV → Neo4j round trip. Consider a KG-validity assertion
in `tests/kg_validity/` that the flag has both values present.

### Repro

```bash
# CSV ground truth
awk -F'\t' '{print $3}' biocypher-out/20260812165914/Interpro_entry_related_to_ec_number-part000.csv \
  | sort | uniq -c

# Live graph
docker exec deploy cypher-shell --format plain "
MATCH ()-[r:Interpro_entry_related_to_ec_number]->()
RETURN apoc.meta.cypher.type(r.ambiguous) AS type, r.ambiguous AS val, count(*) AS n;"

# Scan every CSV-imported boolean in a build
cd biocypher-out/20260812165914
for h in *-header.csv; do b=${h%-header.csv}; cols=$(cat "$h")
  idx=$(echo "$cols" | tr '\t' '\n' | grep -n ':boolean$' | cut -d: -f1)
  [ -z "$idx" ] && continue
  for i in $idx; do nm=$(echo "$cols" | cut -f"$i"); [ -f "$b-part000.csv" ] || continue
    printf "  %-42s %-18s %s\n" "$b" "$nm" \
      "$(awk -F'\t' -v c="$i" '{print $c}' "$b-part000.csv" | sort | uniq -c | tr '\n' ' ')"
  done
done
```

---

## 2. CRITICAL — narrative prose harvested as curated cross-references

**Where:** `multiomics_kg/utils/interpro_reference.py:197` (`parse_entry_db_xrefs`)

The parser scans **every line** between `<interpro id=…>` and `</interpro>`, which includes the
`<abstract>` block. InterPro renders inline EC/CAZy mentions in narrative text as `<db_xref>`
elements, so prose mentions become entry-level annotations indistinguishable from curated ones.

### Measured against the release (2026-08-06)

| DB | `<external_doc_list>` (curated) | `<abstract>` (prose) |
|---|---|---|
| EC | 12,779 xrefs / 10,191 entries | **5,521 xrefs / 3,032 entries** |
| CAZY | 116 / 116 | **223 / 207** — prose is the majority |
| METACYC | 79,683 / 5,091 | **0** |

- **MetaCyc pathways are unaffected** — zero abstract xrefs.
- **GO is unaffected** — it comes from the separate `interpro2go` file, not the XML.
- The defect is confined to **EC and CAZy**.

### Reference-cache contamination

- 658 of 10,849 entries have wholly prose-derived `ec_numbers`; **211 pass Layer B's
  FAMILY + single-EC gate** (5.2% of the 4,094 propagating entries).
- **111 of 227** entries with `cazy_ids` are wholly prose-derived (49%).

### Blast radius in the shipped graph

- **979 gene→EC edges** traceable to prose (265 with InterPro as the *only* source — i.e. net-new,
  uncorroborated catalysis claims).
- **175 gene→CAZy edges** (105 sole-source).

### Two distinct failure modes

**(a) Explicitly hedged prose becomes a hard claim.**
IPR001602 "UPF0047 protein YjbQ-like" — the abstract states members "have been shown to have
*sufficient* thiamine phosphate synthase activity (EC 2.5.1.3) to complement thiE mutants.
**However, it is presumed that this is a secondary activity, and the primary function of the YjbQ
family enzyme remains unknown.**" That EC is now a `Gene_catalyzes_ec_number` edge on **43 genes**.

IPR001128 "Cytochrome P450" picks up EC 1.14.14.95 from a sentence about germacrene A hydroxylase
in **lettuce** — applied to marine bacterial P450s.

**(b) Complex-level EC stamped on structural subunits.**
EC 7.1.2.2 (ATP synthase holoenzyme) is propagated to individual subunit families — F0 subunit
b/b′, subunit c, subunit a, F1 delta/epsilon, gamma, OSCP/delta — ~40 genes each (461 total
EC-7.1.2.2 edges). Subunit c is a proteolipid; it catalyses nothing on its own.

Top contributing entries:

```
IPR006062 FAMILY EC=5.3.1.16  n=94  Histidine biosynthesis protein
IPR002146 FAMILY EC=7.1.2.2   n=69  ATP synthase, F0 complex, subunit b/b'
IPR003445 FAMILY EC=7.1.2.2   n=54  Cation transporter
IPR002698 FAMILY EC=6.3.3.2   n=44  5-formyltetrahydrofolate cyclo-ligase
IPR001602 FAMILY EC=2.5.1.3   n=43  UPF0047 protein YjbQ-like
IPR001469 FAMILY EC=7.1.2.2   n=42  ATP synthase, F1 complex, delta/epsilon subunit
IPR000131 FAMILY EC=7.1.2.2   n=42  ATP synthase, F1 complex, gamma subunit
IPR001249 FAMILY EC=6.4.1.2   n=41  Acetyl-CoA biotin carboxyl carrier
IPR000454 FAMILY EC=7.1.2.2   n=41  ATP synthase, F0 complex, subunit C
IPR000568 FAMILY EC=7.1.2.2   n=41  ATP synthase, F0 complex, subunit A
IPR005953 FAMILY EC=7.1.2.2   n=40  ATP synthase, F0 complex, subunit C, bacterial/chloroplast
IPR000711 FAMILY EC=7.1.2.2   n=39  ATPase, OSCP/delta subunit
```

### The single-EC gate inverts

`build_gene_annotations.py:387` gates EC propagation on `len(ecs) == 1`, intended as the noise
filter ("a multi-EC family is a candidate set, not a claim"). But abstracts typically mention
**exactly one** EC while curated entries carry **several** — so the gate systematically *admits*
prose mentions and *rejects* curated multi-EC sets. Fixing the parser also repairs the gate's
intent.

### Proposed fix

Restrict `parse_entry_db_xrefs` to the curated section. Track the enclosing element and only keep
xrefs inside `<external_doc_list>` (and `<member_list>` if ever needed), never `<abstract>`,
`<pub_list>` or `<structure_db_links>`.

Suggested shape — add a section filter parameter so the pathway wrapper keeps working unchanged:

```python
_SEC_OPEN = re.compile(r'<(abstract|external_doc_list|member_list|class_list'
                       r'|pub_list|structure_db_links)\b')
_SEC_CLOSE = re.compile(r'</(abstract|external_doc_list|member_list|class_list'
                        r'|pub_list|structure_db_links)>')

def parse_entry_db_xrefs(lines, include_dbs, sections=("external_doc_list",)):
    ...
    # track `sect`; only collect when `sect in sections`
```

Then rebuild: `bash scripts/prepare_data.sh --steps 9 --force` (re-parses cached raw files, no
network), followed by steps 2 + a KG rebuild.

**Optional refinement:** rather than discarding prose xrefs entirely, they could be routed to
Layer A (the recall-biased router) with `ambiguous` set — that preserves recall while keeping them
out of `Gene_catalyzes_ec_number`. Requires finding 1 to be fixed first, or the flag is inert.

### Expected effect

- ~211 entries stop propagating EC → ~979 gene→EC edges removed (265 of them sole-source).
- ~76 propagating entries stop contributing CAZy → ~175 gene→CAZy edges removed (105 sole-source).
- `Gene_catalyzes_ec_number` drops from 69,026 to roughly 68,047.
- MetaCyc pathways and all GO unchanged.

### Repro

```bash
python3 - <<'EOF'
import gzip, re
from collections import Counter
open_re = re.compile(r'<interpro\s+id="(IPR\d{6})"')
xref_re = re.compile(r'<db_xref[^>]*\bdb="([A-Z]+)"[^>]*\bdbkey="([^"]+)"')
so = re.compile(r'<(abstract|external_doc_list|member_list|class_list|pub_list|structure_db_links)\b')
sc = re.compile(r'</(abstract|external_doc_list|member_list|class_list|pub_list|structure_db_links)>')
ctx = Counter(); ents = {}
cur = sect = None
with gzip.open("cache/data/interpro/raw/interpro.xml.gz", "rt", errors="replace") as fh:
    for line in fh:
        m = open_re.search(line)
        if m: cur = m.group(1)
        s = so.search(line)
        if s: sect = s.group(1)
        if cur:
            for db, key in xref_re.findall(line):
                if db in ("EC", "CAZY", "METACYC"):
                    k = (db, sect or "TOP"); ctx[k] += 1
                    ents.setdefault(k, set()).add(cur)
        if sc.search(line): sect = None
        if "</interpro>" in line: cur = sect = None
for k in sorted(ctx, key=lambda x: (x[0], -ctx[x])):
    print(f"{k[0]:8s} in <{k[1]:20s}>  xrefs={ctx[k]:7d}  entries={len(ents[k]):6d}")
EOF
```

---

## 3. MEDIUM — `evidence='curated'` is the default for anything InterPro didn't touch

**Where:** `multiomics_kg/utils/annotation_provenance.py:34`

```python
evidence = ev_map.get(token) or "curated"
```

The `<field>_evidence` map is **sparse and InterPro-centric** — only tokens InterPro touched get an
explicit label (set in `_fold_interpro_field`, `build_gene_annotations.py:337-345`). Every other
token falls through to the strongest value on the ladder.

Consequence, measured in the shipped build:

```
429,287  Gene_involved_in_biological_process   evidence='curated'  sources='eggnog'
 11,537  Gene_catalyzes_ec_number              evidence='curated'  sources='eggnog'
```

eggNOG GO/EC terms are **ortholog transfer** — inference by definition. CLAUDE.md documents
`evidence` as "`curated`>`signature`>`family_inferred`>`domain_inferred` — the field to read for
curated-vs-inferred", so a consumer filtering `evidence = 'curated'` to get trustworthy
annotations gets 429K ortholog-transferred GO terms.

Note `_CURATED_SOURCES` includes `eggnog`, so this is deliberate at the source-bucket level — but
"came from a curated database" and "is a curated annotation" are different claims, and the field
name asserts the latter.

### Options

- **(a) Rename the semantics.** Treat `curated` as `from_curated_source` and document it that way.
  Cheapest; no rebuild. But the ladder still can't distinguish eggNOG transfer from a UniProt
  experimental annotation.
- **(b) Give eggNOG its own rung** — e.g. `ortholog_transfer`, slotting below `signature`. Requires
  populating the evidence map for eggNOG tokens (currently untouched), a step-2 rebuild, and a
  post-import/MCP audit.
- **(c) Drop the default** — leave `evidence` absent when unknown rather than defaulting to the
  strongest value, forcing consumers to handle the null.

(b) is the honest fix; (a) is the cheap one. Worth deciding explicitly since MCP filters are
planned on this field (spec §7 follow-up).

---

## 4. MEDIUM — Layer A over-prunes exactly the edges it exists for

**Where:** `multiomics_kg/adapters/interpro_adapter.py:379`

```python
observed_ec = set() if self.ec_node_ids is None else observed_ec & self.ec_node_ids
```

`observed_ec` is built from genes' `ec_numbers`; `ec_node_ids` is the injected Expasy node universe.
The **`ec_node_ids` prune alone is sufficient** for dangling-safety — EcNumber nodes are the full
unpruned Expasy hierarchy (7,337 nodes vs only 2,217 gene-referenced ECs). The extra intersection
with `observed_ec` restricts Layer A targets to ECs some gene *already carries*.

Measured impact:

```
Layer A candidates pruned to EcNumber nodes ONLY : 7,714
  surviving observed_ec intersection (shipped)   : 6,957  (→ 6,854 after id dedup)
  DROPPED by the extra observed_ec intersection  :   757
     ambiguous=true (multi-EC / non-FAMILY)      :   755
     ambiguous=false (single-EC FAMILY)          :     2
```

**755 of 757 dropped edges are `ambiguous=true`** — precisely the multi-EC/DOMAIN router edges
Layer A was built to home. Layer A's stated purpose is to carry the EC candidates Layer B
*refuses*; intersecting with what Layer B already produced defeats it.

The 6,957 → 6,854 gap is the ~103 duplicate edge ids CLAUDE.md already documents (two raw tokens
normalising to one EC); they dedup harmlessly and the shipped CSV has 6,854 rows / 6,854 distinct
ids.

### Proposed fix

Drop the `observed_ec` intersection and prune EC targets to `ec_node_ids` only. Consider the same
review for `observed_cazy` (CazyFamily nodes are observed-only, so that one likely *is* load-bearing
— verify before changing).

---

## 5. LOW — Pfam evidence ordering

`enrich_interpro_fields` runs at `build_gene_annotations.py:1297`; `enrich_pfam_fields` (which
re-keys eggNOG **shortnames** → `PF*` accessions) runs later at `:1345`. So when
`_fold_interpro_field` evaluates `srcs & _CURATED_SOURCES` at `:337`, an eggNOG-asserted domain is
still keyed as e.g. `YjbQ`, not `PF01894` — the overlap isn't seen.

Result: **59,363 edges** ship as `evidence='signature'` with `sources='eggnog|interpro'`, where the
docstring at `:319-322` promises `curated` when a curated source also asserts the token.

**Impact is label-only** — `evidence_score` is identical either way (`signature` and `curated` both
earn +1 on rung 2 and +1 on rung 3; the eggnog+interpro pair collapses to one effective source on
rung 1). Arguably `signature` is the more *accurate* label for a direct HMM hit, so the cheapest
fix may be correcting the docstring rather than the ordering.

---

## 6. LOW — `_normalize_ec_token` parity claim is false

`interpro_adapter.py:72-82` — the docstring says it normalises "the SAME way Layer B does", but it
omits Layer B's `_VALID_EC_RE` filter (`build_gene_annotations.py:295,311`). Currently harmless:
invalid ECs have no EcNumber node, so the `ec_node_ids` prune drops them. Fix the filter or the
docstring so the next change doesn't rely on that accident.

---

## 7. LOW — dead code

`multiomics_kg/utils/annotation_provenance.py:22` defines `_CURATED_SOURCES`, never used in that
module. The live copy is `build_gene_annotations.py:290`. Two definitions of the same constant that
can silently drift; delete the dead one.

---

## 8. LOW — evidence folded across incomparable member DBs

`interpro_adapter.py:104-107` — `evalue = min(...)`, `score = max(...)` across all matches to one
entry, but those come from different member databases (PANTHER, Gene3D, PROSITE, NCBIfam…) whose
e-values and bit-scores are not on a comparable scale. `libraries` is retained so a consumer *can*
tell, but the folded scalars invite cross-DB comparison. Consider documenting them as
"best-of-any-library, not comparable across entries", or storing per-library values.

Note this is evidence-only and never a filter (no e-value cutoff by design), so impact is limited
to display/ranking.

---

## 9. LOW — latent coupling between calls.json and the merged JSON

`InterproAnnotationAdapter.get_edges` (`:205`) builds edge targets from **calls.json**
`matches[].interpro_accession`, while the node set comes from `get_all_interpro_ids()` reading the
**merged JSON**'s `interpro_entries`. Nothing intersects the two.

Verified consistent today — 0 mismatches over 24,039 proteins across 6 strains — because both
derive from the same artifact via `load_interproscan`. But a `/interproscan-run --force` re-scan
without re-running `prepare_data.sh --steps 2` would desynchronise them and silently produce
skipped edges at import. Cheap guard: filter `by_ipr` to `set(gene["interpro_entries"])`, or assert
equality and log.

---

## 10. LOW — documentation drift

CLAUDE.md states `Pfam_in_interpro_entry` (Pfam → node; **~5,236**). Actual in the live graph:
**5,972**. Other documented counts check out (12,999 nodes; 397,342 gene edges; 1,569 hierarchy;
102,895 genes covered ≈ 85%).

Level distribution is 11,430 / 1,490 / 79 for levels 0/1/2 — 87.9% at level 0, vs "~86%"
documented. Close enough, but refresh alongside the bridge count.

---

## What is solid (verified, no action)

- **Zero dangling endpoints** across all six InterPro edge types plus `Gene_catalyzes_ec_number`
  and `Gene_has_cazy_family`. The `4a609b77` EC prune works as intended.
  *(An initial "34,912 dangling" reading during this review was a column-index error on my part —
  `:END_ID` is column 9 in `Gene_has_interpro_entry`, not column 6. Re-verified correctly.)*
- **`InterproEntry.gene_count` is genuinely DIRECT**, not a subtree sum (confirmed: parent entries
  have `gene_count` independent of the sum of their children). This is the ORA-correct semantics
  CLAUDE.md claims, and it correctly differs from `TcdbFamily.gene_count`, which is a subtree count.
- **No duplicate edge ids** in any shipped InterPro CSV.
- **Post-import booleans work** — `is_promiscuous` carries 22 TRUE on InterproEntry; only the
  CSV-imported boolean path is broken.
- **All 2,174 unit tests pass** (`pytest -m "not slow and not kg"`, 13s).

---

## Suggested order of work

1. **Finding 1** (`ambiguous` → categorical string) — self-contained adapter + schema change, plus a
   KG-validity test asserting both values are present. No prepare_data re-run needed.
2. **Finding 2** (restrict parser to `<external_doc_list>`) — `interpro_reference.py` change, then
   `prepare_data.sh --steps 9 --force` (offline, uses cached raw), then step 2, then rebuild.
   Add a unit test with an abstract-embedded `<db_xref>` fixture that must **not** be harvested.
3. **Finding 4** (drop the `observed_ec` intersection) — one line; rides along with the same rebuild.
4. **Finding 3** (`evidence` semantics) — decide (a) / (b) / (c) before MCP filters ship on this
   field; (b) needs a step-2 rebuild so it batches with 2.
5. Findings 5–10 — cleanup, batch freely.

Findings 2, 3 and 4 all change edge counts, so run `/omics-edge-snapshot` before/after and refresh
`tests/kg_validity/snapshot_data.json` once at the end.
