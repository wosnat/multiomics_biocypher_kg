---
name: extract-discussed-topics
description: Extract the genes and KEGG pathways a publication discusses in its narrative (from the full PDF) for the Publication_discusses_gene / Publication_discusses_kegg_pathway edges. Stage 1 of the discuss-edges pipeline. Use when adding narrative literature links for a paper.
argument-hint: [paper-name]
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(uv *), Bash(python *)
---

# Extract Discussed Topics (Stage 1)

Extracts, from a paper's **full PDF**, the genes and KEGG pathways it substantively
**discusses in its narrative** — regulators named in the discussion, genes mentioned but
not differentially expressed, pathways referred to by name — which the supplementary DE
tables do not capture. Output feeds the `Publication_discusses_gene` /
`Publication_discusses_kegg_pathway` edges.

This is **Stage 1** of the three-stage pipeline (extract → resolve → adapter). See
`plans/publication_discusses_edges.md`. It only writes a raw, human-reviewable JSON;
resolution to canonical Gene / KeggTerm nodes happens in Stage 2.

## What it produces

All artifacts live in a `publication_topics/` subfolder of the paper directory:
- `topics.json` — raw extraction (this skill, Stage 1)
- `topics_resolved.json` — resolved to canonical nodes (Stage 2, `resolve_paper_topics`)
- `resolution_report.txt` — per-paper stats + unresolved-mention diagnostics (Stage 2)

`publication_topics/topics.json`:

```json
{
  "metadata": { "paper": "...", "doi": "...", "strains": ["Prochlorococcus MED4", ...],
                "model": "gpt-4.1", "self_assessment": "high", ... },
  "genes": [
    { "surface_form": "NtcA", "gene_name": "ntcA", "identifiers": ["PMM0552"],
      "strain": "Prochlorococcus MED4", "prominence": "central",
      "evidence": "‘ntcA was strongly induced under N stress’ (Results, p4)" }
  ],
  "pathways": [
    { "surface_form": "nitrogen metabolism", "kegg_pathway_id": "ko00910",
      "prominence": "central", "evidence": "..." }
  ]
}
```

Key fields:
- **identifiers** — every locus tag / ORF / protein id the paper prints, verbatim. These
  resolve 1:1 via Tier-1 `specific_lookup` in Stage 2, and the locus-tag prefix doubles as
  the strain fingerprint. Highest-value field.
- **strain** — one of the paper's candidate strains (from `experiments[].organism`), or
  `all` (conserved/generic across the paper's strains) or `unspecified`.
- **prominence** — `central` (a focus of the paper) or `peripheral` (passing mention).

## Usage

```bash
# See what would run (no API calls) — lists papers + strain counts + EXISTS/NEW
uv run python -m multiomics_kg.extraction.topics.extract --dry-run

# Extract one paper (recommended first, then review the JSON)
uv run python -m multiomics_kg.extraction.topics.extract --paper "Aharonovich 2016"

# Re-extract (overwrite existing)
uv run python -m multiomics_kg.extraction.topics.extract --paper "Aharonovich 2016" --force

# Extract every paper that has a PDF and no existing extraction
uv run python -m multiomics_kg.extraction.topics.extract

# Diff-friendly markdown report from existing JSONs (no API calls)
uv run python -m multiomics_kg.extraction.topics.extract --report
```

Model defaults to `gpt-4.1` (override with `--model` or `TOPICS_EXTRACTION_MODEL`).
Requires `OPENAI_API_KEY` in `.env`.

## Resolution (Stage 2)

After extraction, resolve mentions to canonical Gene / KEGG-pathway nodes:

```bash
uv run python -m multiomics_kg.download.resolve_paper_topics --papers "Aharonovich 2016" --force
```

This writes `publication_topics/topics_resolved.json` + `resolution_report.txt`. The report
gives per-paper resolution rates, a method breakdown, a truncated-identifier count, and an
`unresolved_reasons` tally (`descriptive_only` = function/enzyme phrase with no symbol/id;
`truncated_id_only` = bare locus number; `lookup_miss` = had a name/id but no mapping match).

## Review checklist

After running on a paper, open `publication_topics/topics.json` and verify:
1. **Identifiers captured** — every gene the paper prints a locus tag / ORF / protein id
   for has it in `identifiers`. This is the explicit requirement; gaps here hurt Stage 2.
2. **Strain attribution** — multi-strain papers: each gene points at the right strain, or
   `all`/`unspecified` only when genuinely ambiguous. Cross-check against `identifiers`
   prefixes (PMM*→MED4, PMT*→MIT9313, P9301_*→MIT9301, SYNW*→WH8102, …).
3. **Prominence** — `central` genes/pathways really are foci of the paper.
4. **No hallucination** — every entry has an `evidence` quote that actually appears in the
   paper; nothing invented from outside knowledge.
5. **Sanitization** — no `|` or stray single quotes used as delimiters in any field.

Then proceed to Stage 2 (resolution).
