# Plan: Quarto Manuscript for Multi-Omics KG Paper

## Context

Write a bioinformatics resource paper about the marine cyanobacteria multi-omics knowledge graph. The paper will be authored collaboratively with Claude Code using Quarto with Typst PDF backend for fast iteration and LaTeX fallback for journal submission.

**Target journals:** Bioinformatics (Oxford), GigaScience, Database (Oxford), PLOS Computational Biology

**Paper angle:** The KG as an LLM-queryable bioinformatics resource — architecture, gene ID harmonization, data integration pipeline, MCP-based LLM access layer, and validation use cases.

### Scope: Two Repositories, One Paper

This paper covers two repositories that together form the complete resource:

| Repo | Role | Status |
|------|------|--------|
| `second_multiomics` (this repo) | KG construction: data integration, gene ID harmonization, pipeline, Docker deployment | Complete |
| `multiomics_explorer` (`/home/osnat/github/multiomics_explorer`) | LLM access layer: MCP server (7 tools), NL→Cypher agent, CLI, evaluation framework | Core working; Streamlit UI and multi-hop reasoning are stubs |

A separate **biology paper** (future, third repo) will use the KG+explorer for biological discovery and cite this resource paper for both the KG and the tooling.

### What's in multiomics_explorer (for paper purposes)

**Working and paper-ready:**
- MCP server (`mcp_server/`) — 7 specialized tools: `get_schema`, `search_genes`, `get_gene_details`, `query_expression`, `compare_conditions`, `get_homologs`, `run_cypher`
- NL→Cypher agent (`agents/cypher_agent.py`) — LangChain GraphCypherQAChain, provider-agnostic LLM
- CLI (`cli/main.py`) — 5 commands: schema, cypher, stats, query, interactive
- Schema introspection (`kg/schema.py`) — live from Neo4j, no static files
- Curated queries + few-shot examples (`kg/queries.py`)
- Evaluation data (`evaluation_data/`) — 3-stage YAML test cases (KG validation, NL→Cypher pairs, reasoning cases)
- Configuration (`config/prompts.yaml`) — system prompts with biological context

**Draft/future work (mention in Discussion):**
- Streamlit UI (stub)
- Multi-hop reasoning agent via LangGraph (stub)
- Programmatic API (stub)

---

## Why Quarto with Typst Backend

| Criterion | Decision |
|-----------|----------|
| Code execution | Quarto runs Python/Jupyter natively — figures generated inline from Neo4j queries |
| PDF rendering | Typst backend: sub-second compilation vs 10-90s for LaTeX |
| Multi-format | One source → Typst PDF (fast), LaTeX PDF (submission), HTML (sharing), DOCX (collaborators) |
| Journal templates | PLOS, Elsevier, Springer Nature, bioRxiv available; Oxford accepts rendered LaTeX/DOCX |
| Reproducibility | `freeze: auto` caches computation; paper compiles without running Neo4j |
| Collaboration | Plain-text `.qmd` files; git-friendly; Claude Code reads/writes sections directly |
| Maturity | Stable (Posit-backed), large academic user base, successor to R Markdown |

**Standalone Typst was rejected** because it has no meaningful Python code execution (pyrunner can't import pandas/matplotlib) and very few biology journal templates. Quarto-with-Typst gives both Quarto's computation features and Typst's fast rendering.

---

## Project Structure

```
paper/
  _quarto.yml              # Manuscript project config
  index.qmd                # Front matter + section includes (thin shell)
  supp_index.qmd           # Supplementary materials shell (includes from supplementary/)
  references.bib           # BibTeX bibliography
  sections/                # One .qmd per section — pure content, no front matter
    introduction.qmd
    system-overview.qmd       # Data sources, schema, pipeline, deployment, QA
    gene-id-harmonization.qmd # Three-tier system + position merge + resolution rates
    expression-integration.qmd # Omics adapter + coverage + homology propagation
    llm-integration.qmd      # MCP server, NL→Cypher agent, tool design
    llm-evaluation.qmd       # NL→Cypher accuracy, tool selection, gene finder, E2E QA, hallucination
    use-cases.qmd             # LLM-driven validation queries
    discussion.qmd
    data-availability.qmd
  notebooks/               # Computational notebooks for figures
    helpers.py             # Shared Neo4j connection + plot styling utilities
    fig-schema.qmd         # Figure: Graph schema visualization
    fig-pipeline.qmd       # Figure: Pipeline architecture
    fig-gene-id-resolution.qmd  # Figure: ID resolution rates per strain
    fig-expression-coverage.qmd # Figure: Expression heatmap (pubs × organisms)
    fig-organism-stats.qmd      # Figure: Gene/protein counts per organism
    fig-query-examples.qmd      # Figure: Biological use-case queries
  static/                  # Hand-drawn diagrams (SVG from draw.io)
    pipeline-architecture.svg
  supplementary/           # Supplementary materials
    supp-text-position-merge.qmd    # S1: Full position-based merge method
    supp-text-id-tier-rules.qmd     # S2: Detailed tier classification rules
    supp-text-mcp-tool-specs.qmd    # S3: MCP tool input/output schemas
    supp-text-eval-pairs.qmd        # S4: NL→Cypher evaluation pairs
    supp-text-prompts.qmd           # S5: System prompts with biological context
    supp-text-eval-results.qmd     # S6: LLM evaluation detailed results
    supp-table-publications.qmd     # Table S1: All integrated publications
    supp-table-gene-id-rates.qmd    # Table S2: Per-strain ID resolution detail
    supp-table-paper-match.qmd      # Table S3: Per-paper match rates
    supp-table-eval-summary.qmd    # Table S4: LLM evaluation summary
    supp-fig-resolution-detail.qmd  # Fig S1: Per-strain resolution breakdown
```

### Multi-file design

- **`index.qmd`** is a thin shell: YAML front matter + section headers + `{{< include >}}` shortcodes. This keeps the full paper outline visible in one file.
- **Section files** contain pure content (prose, code blocks, cross-references) — no YAML front matter.
- **Cross-references work across includes** — `@fig-schema` in a results section can reference a figure defined in a notebook, since Quarto resolves everything after assembly.
- **Benefits**: each section can be edited independently without touching other parts; git diffs are isolated to the section changed; Claude Code can draft/revise one section at a time.

---

## Quarto Configuration (`paper/_quarto.yml`)

```yaml
project:
  type: manuscript
  output-dir: _output

manuscript:
  article: index.qmd
  notebooks:
    - notebook: supp_index.qmd
      title: "Supplementary Materials"
    - notebook: notebooks/fig-schema.qmd
      title: "Graph Schema Visualization"
    - notebook: notebooks/fig-gene-id-resolution.qmd
      title: "Gene ID Resolution Performance"
    - notebook: notebooks/fig-expression-coverage.qmd
      title: "Expression Data Coverage"
    - notebook: notebooks/fig-organism-stats.qmd
      title: "Organism Statistics"
    - notebook: notebooks/fig-query-examples.qmd
      title: "Query Use Cases"
  code-links:
    - text: "GitHub Repository"
      icon: github
      href: https://github.com/wosnat/second_multiomics

format:
  typst:
    keep-typ: true
  pdf:
    documentclass: article
    keep-tex: true
  html:
    theme: cosmo
  docx: default

execute:
  freeze: auto
  echo: false
  warning: false

bibliography: references.bib
csl: bioinformatics.csl
```

Key settings:
- `freeze: auto` — notebooks only re-execute when source changes; paper compiles without Neo4j after first render
- `echo: false` — hide code blocks in the rendered article
- Four output formats available; default to Typst PDF for daily iteration

---

## `index.qmd` Structure (thin shell with includes)

```qmd
---
title: "A Multi-Omics Knowledge Graph for Marine Cyanobacteria..."
author:
  - name: ...
    affiliations: ...
abstract: |
  ...
keywords: [knowledge graph, cyanobacteria, multi-omics, gene ID harmonization]
---

# Introduction
{{< include sections/introduction.qmd >}}

# System Overview
{{< include sections/system-overview.qmd >}}

# Gene ID Harmonization
{{< include sections/gene-id-harmonization.qmd >}}

# Expression Data Integration
{{< include sections/expression-integration.qmd >}}

# LLM Integration
{{< include sections/llm-integration.qmd >}}

# LLM Evaluation
{{< include sections/llm-evaluation.qmd >}}

# Use Cases
{{< include sections/use-cases.qmd >}}

# Discussion
{{< include sections/discussion.qmd >}}

# Data Availability
{{< include sections/data-availability.qmd >}}

# References
```

---

## `supp_index.qmd` Structure (supplementary shell)

```qmd
---
title: "Supplementary Materials"
---

# Supplementary Text S1: Position-Based Annotation Merge
{{< include supplementary/supp-text-position-merge.qmd >}}

# Supplementary Text S2: Gene ID Tier Classification Rules
{{< include supplementary/supp-text-id-tier-rules.qmd >}}

# Supplementary Text S3: MCP Tool Specifications
{{< include supplementary/supp-text-mcp-tool-specs.qmd >}}

# Supplementary Text S4: NL→Cypher Evaluation Pairs
{{< include supplementary/supp-text-eval-pairs.qmd >}}

# Supplementary Text S5: System Prompts
{{< include supplementary/supp-text-prompts.qmd >}}

# Supplementary Text S6: LLM Evaluation Detailed Results
{{< include supplementary/supp-text-eval-results.qmd >}}

# Supplementary Tables
{{< include supplementary/supp-table-publications.qmd >}}
{{< include supplementary/supp-table-gene-id-rates.qmd >}}
{{< include supplementary/supp-table-paper-match.qmd >}}
{{< include supplementary/supp-table-eval-summary.qmd >}}

# Supplementary Figures
{{< include supplementary/supp-fig-resolution-detail.qmd >}}
```

---

## Writing Guidelines

### Tone
Academic throughout. Formal scientific register appropriate for bioinformatics journals (Bioinformatics, Database, GigaScience). Avoid blog-style explanations, colloquialisms, and first-person narrative except where convention allows ("we present", "we developed"). Claims must be precise and supported by data or citations.

### Specifics as Examples
The main text presents **general methods and principles**. Specific organisms, strains, papers, or ID types appear only as **illustrative examples** — never as exhaustive lists. Comprehensive details (per-strain breakdowns, full ID type tables, complete publication lists, tool schemas, evaluation pair-by-pair results) belong in **supplementary materials**.

**Example pattern for main text:**
> "The three-tier system classifies identifier types by their expected uniqueness. Tier 1 identifiers (e.g., locus tags and old locus tags) map unambiguously to a single gene, while Tier 2 identifiers (e.g., RefSeq protein IDs) may map to paralogs and are accepted only when the mapping is unambiguous for the organism in question (Supplementary Text S2)."

**Not:**
> "Tier 1 includes locus_tag, locus_tag_ncbi, locus_tag_cyanorak, old_locus_tag, alternative_locus_tag, jgi_id, probeset_id, and uniprot_entry_name. These are stored in specific_lookup as 1:1 mappings..."

This keeps the main text readable and journal-length while preserving rigor through the supplements.

---

## Paper Structure

### Title Candidates

1. "A Multi-Omics Knowledge Graph for Marine Cyanobacteria: Integrating Heterogeneous Gene Identifiers and Differential Expression Data for LLM-Assisted Biological Discovery"
2. "ProCyanoGraph: A BioCypher Knowledge Graph Integrating Transcriptomic and Proteomic Data Across 13 Marine Cyanobacterial Strains"
3. "From Spreadsheets to Subgraph Queries: A Knowledge Graph for Cross-Study Multi-Omics Integration in Marine Microbiology"

### Section Outline (feature-driven — method + results unified per topic)

**1. Introduction**
- Marine cyanobacteria as model systems; the cross-study data integration challenge
- Why knowledge graphs suit multi-omics integration (contrast with RAG, fine-tuning)
- Contribution statement: KG + LLM access layer as an integrated resource
- Sources: `doc/blog posts/1 why omics and llms.md`, `doc/blog posts/3 knowledge graph.md` — formalize from blog tone to academic register

**2. System Overview** — architecture and deployment
- Data sources and organism selection — describe the rationale (ecotype coverage, data availability); present representative strains as examples, full strain table as Table 1
- Graph schema and ontology alignment (Biolink via BioCypher) — describe the design principles; reference `config/schema_config.yaml`
- Pipeline architecture — adapter pattern described generically with one adapter as a worked example; Docker deployment
- Graph statistics — summary counts (Figure: organism stats); per-type breakdowns in text only where they illustrate a point
- Quality assurance — describe the testing philosophy (regression snapshots, biological invariants); cite specific tests as examples, not an exhaustive catalog

**3. Gene ID Harmonization** — the core technical contribution
- The gene ID problem — frame generally (heterogeneous identifiers across sequencing eras and databases), with one or two concrete examples (e.g., a gene known by four different IDs across three studies)
- Three-tier resolution system — describe the principle (uniqueness-based classification) with representative ID types as examples; full tier rules → **Supp S2**
  - **Source: `docs/methods_gene_id_mapping.md` — condense from 2,800 words to ~1,200 for main text**
- Position-based annotation merge — summarize the approach (~400 words); MIT9313 as the illustrative case; full method → **Supp S1**
  - **Source: `docs/methods_position_fallback_merge.md` — condense from 2,400 words**
- Resolution performance — overall rates and a few notable strains as examples (Figure: resolution rates); per-strain detail → **Supp Table S2**

**4. Expression Data Integration** — omics data in the graph
- Omics adapter — describe the paperconfig-driven design generically; one publication as a worked example
- Expression coverage — summary statistics; publications × organisms heatmap (Figure); full publication list → **Supp Table S1**
- Homology propagation — explain the principle (Cyanorak clusters → cross-strain edges) with one gene family as an example
- Homology network properties — aggregate statistics

**5. LLM Integration** — making the graph queryable by AI agents
- Design goals: structured tool access vs. free-form Cypher; write-blocking for safety
- Gene finder: multi-strategy gene resolution — describe the approach (exact lookup → synonym resolution → text search); connect to the three-tier system from §3; one worked example showing the resolution chain
- MCP server architecture — describe the tool design philosophy (specialized tools vs. monolithic query); list the 7 tools with one-line descriptions; full input/output schemas → **Supp S3**
  - **Source: `multiomics_explorer/mcp_server/tools.py`**
- NL→Cypher translation agent — describe the approach (few-shot prompting, schema injection); one example query pair; full prompt text → **Supp S5**, evaluation pairs → **Supp S4**
  - **Source: `multiomics_explorer/agents/cypher_agent.py`, `multiomics_explorer/config/prompts.yaml`**

**6. LLM Evaluation** — systematic assessment of the LLM access layer
- Describe the evaluation framework: five complementary assessment types, each testing a different layer of the pipeline
- 6.1 NL→Cypher accuracy — methodology and summary metrics; per-query breakdowns → **Supp S6**
- 6.2 Tool selection accuracy — methodology and summary metrics
- 6.3 Gene finder precision/recall — methodology with example test cases; full test set → **Supp S6**
- 6.4 End-to-end QA — methodology with one representative question as a worked example
- 6.5 Hallucination detection — methodology with one illustrative probe
- Main text: summary table of accuracy metrics across all five evaluation types
- **All detailed results → Supp S6**
- **Source: `multiomics_explorer/evaluation/` (to be built) + `evaluation_data/` (existing YAML test cases)**

**7. Use Cases** — LLM-driven validation queries
- Two to three end-to-end examples: NL question → tool call → Cypher → result → biological interpretation
- Select examples that demonstrate cross-study or cross-condition comparisons
- Framed as validation of the KG+LLM pipeline, not as biological discovery (companion paper)

**8. Discussion**
- Comparison to existing resources (STRING, BioCyc, CyanoBase) — what this resource adds
- KG-as-tool paradigm: structured graph access vs. RAG vs. fine-tuning for domain-specific omics data
- Limitations — presented frankly (orphan proteins, Cyanorak server availability, taxonomic scope)
- Future work: Streamlit UI, multi-hop reasoning, exoenzyme pipeline, companion biology paper

**9. Data Availability**
- GitHub repos (KG builder + explorer), Docker deployment, Neo4j dump

---

## Content Mapping: Existing Material → Paper Sections

| Existing Content | Repo | Paper Section | Transformation |
|---|---|---|---|
| `doc/blog posts/1 why omics and llms.md` | KG | Introduction (para 1-2) | Rewrite in academic register; remove conversational framing |
| `doc/blog posts/3 knowledge graph.md` | KG | Introduction (para 3-5: KG vs RAG) | Formalize argument; add citations to KG/RAG literature |
| `doc/blog posts/2 leve of details.md` | KG | Expression Integration (opening) | Extract the "level of detail" rationale (why fold change, not raw counts); present as a design decision with brief justification |
| `data/.../cyanobacteria_genomes.csv` | KG | System Overview — Table 1 | Format as publication table; main text references representative strains |
| `config/schema_config.yaml` | KG | System Overview — schema | Describe design principles in prose; schema diagram (Figure 2) conveys structure |
| KG validity tests (`tests/kg_validity/`) | KG | System Overview — QA | Describe testing philosophy; cite specific tests as examples of biological invariants |
| `docs/methods_gene_id_mapping.md` (2,800 words) | KG | §3 (~1,200 words) + Supp S2 (full) | Condense to general principles with examples; exhaustive ID type tables → S2 |
| `docs/methods_position_fallback_merge.md` (2,400 words) | KG | §3 (~400 words) + Supp S1 (full) | Summarize approach; MIT9313 as illustrative example; full case study → S1 |
| `mcp_server/tools.py` (368 lines, 7 tools) | Explorer | §5 — MCP server | Describe design philosophy; one-line tool descriptions; full schemas → Supp S3 |
| `agents/cypher_agent.py` (102 lines) | Explorer | §5 — NL→Cypher | Describe approach with one example query; full prompt → Supp S5 |
| `config/prompts.yaml` | Explorer | §5 — prompt engineering | Describe biological context strategy; full prompts → Supp S5 |
| `evaluation_data/*.yaml` (3 stages) | Explorer | §6 — evaluation | Reference as test data; example pairs in main text; full set → Supp S4/S6 |
| `evaluation/` (to be built) | Explorer | §6 — metrics | Build: evaluation runner; main text gets summary table; detailed results → Supp S6 |
| `AGENT.md` | Explorer | §5 | Architectural overview as source material |
| `.claude/skills/cypher-queries/SKILL.md` | KG | §7 Use Cases | Select 2-3 representative queries that illustrate cross-study capability |
| `docs/paper/kg_improvements_exoenzyme_analysis.md` | KG | Discussion (future work) | Brief mention as planned extension; not a detailed description |

---

## Figure Plan

| # | Figure | Section | Data Source | Type | Notebook |
|---|--------|---------|-------------|------|----------|
| 1 | Pipeline architecture | System Overview | Hand-drawn | Static SVG | `static/pipeline-architecture.svg` |
| 2 | Graph schema | System Overview | `schema_config.yaml` | graphviz/networkx | `fig-schema.qmd` |
| 3 | Gene ID resolution rates | Gene ID Harmonization | `gene_id_mapping_report.json` files | Grouped bar chart | `fig-gene-id-resolution.qmd` |
| 4 | Expression coverage | Expression Integration | Neo4j (pubs × organisms) | Heatmap | `fig-expression-coverage.qmd` |
| 5 | Organism statistics | System Overview | Neo4j (gene/protein counts) | Stacked bars | `fig-organism-stats.qmd` |
| 6 | LLM query flow | LLM Integration | Hand-drawn | Static SVG | `static/llm-query-flow.svg` |
| 7 | LLM evaluation results | LLM Evaluation | evaluation runner output | Summary table/bar chart | `fig-eval-results.qmd` |
| 8 | Use case examples | Use Cases | Neo4j (via MCP tools) | Multi-panel | `fig-query-examples.qmd` |
| 9 | Graph neighborhood: gene context | System Overview or Use Cases | Neo4j Browser export | Neo4j graph viz | `static/neo4j-gene-neighborhood.svg` |
| 10 | Graph neighborhood: cross-strain expression | Use Cases | Neo4j Browser export | Neo4j graph viz | `static/neo4j-expression-subgraph.svg` |

**Neo4j graph visualizations** (Figures 9–10): exported from Neo4j Browser as SVG. Show real subgraphs — e.g., a gene node with its protein, organism, expression edges, homologs, and GO terms — to complement the abstract schema diagram (Figure 2). Generated via Cypher queries in the browser; saved as static SVGs.

Figures built incrementally as each section is written — not all upfront.

---

## Supplementary Materials

The main text presents general methods and results, using specific cases as illustrative examples. All exhaustive detail — per-strain breakdowns, complete ID type catalogs, full tool schemas, evaluation pair-by-pair results, and comprehensive publication lists — goes to supplements. This keeps the main text at journal length (target ~4,000–6,000 words depending on venue) while preserving full reproducibility and rigor through the supplementary materials.

### Supplementary Texts

| ID | Title | Source Material | Main Text Treatment |
|----|-------|----------------|---------------------|
| S1 | Position-based annotation merge: full method | `docs/methods_position_fallback_merge.md` (2,400 words) | §3 gets ~800 word summary; full MIT9313 case study in S1 |
| S2 | Gene ID tier classification rules | `docs/methods_gene_id_mapping.md` + `CLAUDE.md` | §3 describes the three-tier concept; S2 has exhaustive ID type tables, conflict resolution rules, iterative convergence details |
| S3 | MCP tool specifications | `multiomics_explorer/mcp_server/tools.py` | §5 describes design philosophy; S3 has full input/output schemas for all 7 tools |
| S4 | NL→Cypher evaluation pairs | `multiomics_explorer/evaluation_data/stage2_nl_cypher_pairs.yaml` | §5 summarizes evaluation framework; S4 lists all 12 NL→Cypher test cases with expected queries |
| S5 | System prompts | `multiomics_explorer/config/prompts.yaml` | §5 describes prompt strategy; S5 has full prompts with biological context and few-shot examples |
| S6 | LLM evaluation detailed results | `multiomics_explorer/evaluation/` outputs | §6 has summary accuracy table; S6 has per-question breakdowns, error analysis, hallucination examples |

### Supplementary Tables

| ID | Title | Source | Notes |
|----|-------|--------|-------|
| Table S1 | Integrated publications | Neo4j + paperconfig files | Full list of ~19 papers: DOI, organism, DE type, condition, edge count |
| Table S2 | Per-strain gene ID resolution | `gene_id_mapping_report.json` per strain | Detailed breakdown: per-ID-type stats, reclassification warnings |
| Table S3 | Per-paper match rates | check-gene-ids output | Matched / RNA-skipped / mismatched counts per paper |
| Table S4 | LLM evaluation summary | `multiomics_explorer/evaluation/` | Accuracy per evaluation type, per LLM provider |

### Supplementary Figures

| ID | Title | Source | Notes |
|----|-------|--------|-------|
| Fig S1 | Per-strain resolution breakdown | `gene_id_mapping_report.json` | Expanded version of main Figure 3, showing per-ID-type resolution for each strain |

---

## Shared Notebook Utilities (`paper/notebooks/helpers.py`)

- `run_cypher(query, **params)` → pandas DataFrame (Neo4j at localhost:7687)
- `set_paper_style()` — Arial 8pt, 300 DPI, no top/right spines
- `ORGANISM_COLORS` — consistent palette: HLI=red, HLII=orange, LLI=green, LLII=light green, LLIV=purple, Syn=blue, Alt=brown/pink
- `PROJECT_ROOT` — computed path for importing project modules

---

## Dependencies

**No new Python packages needed.** Existing `pyproject.toml` already includes: `seaborn` (brings matplotlib), `neo4j`, `pandas`, `jupyter`. Only add `networkx`/`graphviz` if schema diagram is programmatic.

**External tool:** Quarto CLI (not a Python package). Typst is bundled with Quarto.

```bash
# Install Quarto
wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.42/quarto-1.6.42-linux-amd64.deb
sudo dpkg -i quarto-1.6.42-linux-amd64.deb
quarto check

# Register project kernel
uv run python -m ipykernel install --user --name=multiomics-kg

# Download citation style
curl -o paper/bioinformatics.csl https://raw.githubusercontent.com/citation-style-language/styles/master/bioinformatics.xml
```

---

## Claude Code Skill: `/paper`

Create `.claude/skills/paper/SKILL.md` with commands:
- `render` — all formats
- `render --to typst` — fast Typst PDF
- `preview` — live hot-reload
- `wordcount` — article word count
- `status` — what exists, frozen notebooks, last render
- `cite <DOI>` — fetch BibTeX from doi.org, append to `references.bib`
- `cite search <query>` — search CrossRef/Semantic Scholar, present matches, fetch selected BibTeX
- `cite check` — scan all `.qmd` files for `@cite-key` references, report any missing from `references.bib`
- `cite unused` — find entries in `references.bib` not referenced in any `.qmd` file

Workflow instructions for iterative section drafting, figure generation, and citation management.

---

## Implementation Sequence

| Step | Action |
|------|--------|
| 1 | Install Quarto CLI |
| 2 | Create `paper/` directory structure + `_quarto.yml` |
| 3 | Create `index.qmd` with front matter + section headers |
| 4 | Create `references.bib` skeleton |
| 5 | Create `notebooks/helpers.py` |
| 6 | First render test (skeleton only) |
| 7 | Create `/paper` Claude Code skill |
| 8 | Add `.gitignore` entries (`paper/_output/`, `paper/_freeze/`, `paper/**/*.pdf`, `paper/**/*.typ`) |
| 9 | Begin writing — Methods sections first (most content exists) |
| 10 | Build figures as Results sections develop |

---

## Collaborative Workflow: User + Claude Code

See `docs/paper/quarto_collaborative_workflow.md` for the full methodology section describing how the paper is co-authored with an AI coding assistant using Quarto.

---

## Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| Quarto manuscript type + Typst compatibility gaps | Fall back to `type: default` if needed; core workflow identical |
| Neo4j must run for notebook execution | `freeze: auto` caches outputs; paper compiles without Neo4j after first render |
| Large intermediate files | `.gitignore` for `_output/`, `_freeze/`, `paper/**/*.pdf`, `paper/**/*.typ` (scoped to `paper/` so `data/` PDFs remain tracked) |
| Journal requires LaTeX source | Switch to `format: pdf` at submission; same `.qmd` source |
| Python imports from `multiomics_kg/` in notebooks | `helpers.py` computes `PROJECT_ROOT`; `uv sync` installs package |
| Cross-repo coordination | Paper lives in KG repo; explorer content referenced by path, not imported |
| Explorer repo still in draft | Core (MCP server, CypherAgent) is working; stubs mentioned as future work in Discussion |
