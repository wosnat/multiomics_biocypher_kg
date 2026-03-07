# AI-Assisted Scientific Paper Generation Using Quarto and Claude Code

## Overview

This paper was co-authored using a structured collaboration between a domain scientist and an AI coding assistant (Claude Code, Anthropic) within a Quarto manuscript project. This section describes the methodology, tooling, and workflow that enabled this collaboration, which we believe represents a reproducible model for AI-assisted scientific writing in computational biology.

## Motivation

Scientific papers in bioinformatics typically combine three elements that are traditionally authored in disconnected workflows: (1) prose describing methods and results, (2) computational analyses generating figures and statistics, and (3) bibliography management. Quarto unifies these into a single reproducible document where Python code blocks execute inline, generating figures that are automatically numbered, captioned, and cross-referenced within the text.

The addition of an AI coding assistant addresses a complementary challenge: the domain scientist holds the biological expertise and architectural knowledge of the system, but faces the time-intensive work of translating that knowledge into polished prose, generating publication-quality figures, and maintaining internal consistency across sections. The AI assistant can draft sections from structured notes, generate figures by querying the live knowledge graph, manage citations, and maintain cross-references — all within the same repository and toolchain used to build the knowledge graph itself.

## Toolchain

### Quarto with Typst Backend

We used Quarto (v1.6+), an open-source scientific publishing system built on Pandoc, configured with the Typst PDF backend for daily iteration and LaTeX as a fallback for journal submission. The choice of Typst over LaTeX as the primary backend was driven by compilation speed: Typst renders a complete manuscript in under one second, compared to 10–90 seconds for LaTeX, enabling a tight edit–render–review cycle during collaborative writing sessions.

The Quarto manuscript project type (`project: type: manuscript`) provided several features critical to our workflow:

- **Multi-format output**: A single `.qmd` source file rendered to Typst PDF (for review), HTML (for sharing), and DOCX (for collaborators preferring Word), eliminating format-conversion overhead.
- **Computational notebooks as supplements**: Analysis notebooks in `paper/notebooks/` were rendered as browsable HTML alongside the article, providing a transparent record of every figure-generating computation.
- **Freeze mechanism**: `execute: freeze: auto` cached notebook outputs after first execution, allowing the paper to compile in seconds without requiring a running Neo4j instance. Only changes to notebook source code triggered re-execution.
- **Embedded notebook outputs**: The `{{< embed >}}` shortcode pulled specific figures from computational notebooks into the article, keeping analysis code separated from prose while maintaining the link between computation and presentation.

### Claude Code as Writing Partner

Claude Code (Anthropic) is a CLI-based AI coding assistant that operates directly in the project repository with access to the file system, terminal, and version control. For paper writing, its capabilities were extended through a custom project skill (`/paper`) that provided:

- Rendering commands (`quarto render --to typst` for fast PDF, `quarto render --to html` for web preview)
- Word count tracking
- Section status overview
- Workflow instructions for iterative drafting

The AI assistant had full read access to the knowledge graph codebase, configuration files, test suites, existing methods documentation, and the deployed Neo4j instance. This meant it could:

1. **Draft sections from structured notes**: Existing methods documents (`docs/methods_gene_id_mapping.md`, `docs/methods_position_fallback_merge.md`) and architecture documentation (`CLAUDE.md`) were transformed into journal-formatted prose.
2. **Generate figures programmatically**: The assistant wrote Python notebooks that queried Neo4j via Cypher, computed statistics, and produced publication-quality matplotlib/seaborn figures — all within the Quarto framework.
3. **Maintain internal consistency**: Cross-references between figures, tables, and sections were managed through Quarto's `@fig-` / `@tbl-` / `@sec-` labeling system.
4. **Populate the bibliography**: BibTeX entries were added to `references.bib` as sections referenced new works, using DOI-based lookup.

### Repository Integration

The manuscript lived in a `paper/` directory within the same git repository as the knowledge graph code:

```
second_multiomics/
  multiomics_kg/          # KG source code
  config/                 # Graph schema, BioCypher config
  tests/                  # KG validity tests
  docs/                   # Methods documentation
  paper/                  # Quarto manuscript
    _quarto.yml           # Project configuration
    index.qmd             # Main article
    references.bib        # Bibliography
    notebooks/            # Figure-generating notebooks
      helpers.py          # Shared Neo4j + plotting utilities
    static/               # Hand-drawn diagrams
```

This co-location meant that:
- Figure notebooks could import project modules directly (e.g., `from multiomics_kg.utils import gene_id_utils`)
- The paper's computational results derived from the exact same codebase and data as the knowledge graph itself
- Version control tracked both code and manuscript changes together
- CI/CD could potentially render the manuscript alongside running KG tests

## Collaborative Workflow

### Section-by-Section Drafting

The paper was written iteratively, one section at a time. A typical writing session followed this pattern:

1. **Context loading**: The scientist identified the next section to write and pointed the AI assistant to relevant source material (existing documentation, code files, Neo4j query results).
2. **First draft**: The AI assistant produced a draft section in `index.qmd`, drawing on the source material and adhering to the target journal's style conventions.
3. **Render and review**: The scientist rendered the Typst PDF (`quarto render --to typst`, <1 second) and reviewed the formatted output.
4. **Iterative refinement**: The scientist provided feedback ("tighten paragraph 2", "add comparison to BioCyc here", "this claim needs a citation"), and the assistant revised in place.
5. **Commit**: Once the section reached an acceptable state, changes were committed to git.

### Figure Generation

Each computational figure was developed as a separate Quarto notebook in `paper/notebooks/`:

1. The AI assistant wrote a `.qmd` notebook containing Python code that queried Neo4j (via the shared `helpers.py` utilities) and generated a matplotlib/seaborn figure.
2. The notebook was rendered in isolation (`quarto render notebooks/fig-expression-coverage.qmd`) to verify the output.
3. The figure was embedded in the main article via Quarto's cross-reference system (`@fig-expression-coverage`).
4. The `freeze: auto` setting cached the rendered output, so subsequent full renders did not require Neo4j to be running.

This approach ensured that every figure in the paper was traceable to a specific Cypher query and Python script, and could be regenerated from the live knowledge graph at any time.

### Quality Control

Several mechanisms maintained quality throughout the collaboration:

- **Domain expert review**: Every AI-generated draft was reviewed by the scientist before rendering. The AI assistant proposed; the scientist disposed.
- **Render verification**: The Typst PDF provided immediate visual feedback on formatting, figure placement, and cross-reference integrity.
- **Git history**: All changes were version-controlled, enabling rollback and diff-based review of AI contributions.
- **Existing test suite**: The KG validity tests (`tests/kg_validity/`) served as ground truth — any statistics or claims in the paper could be verified against the test assertions.
- **Freeze integrity**: Cached notebook outputs in `_freeze/` ensured that paper content remained stable even as the underlying KG evolved between writing sessions.

## Reproducibility Considerations

The Quarto manuscript framework provides several layers of reproducibility:

1. **Computational reproducibility**: All figures and statistics derive from executable notebooks that query the deployed knowledge graph. The notebooks, Neo4j connection parameters, and Cypher queries are all committed to the repository.
2. **Environment reproducibility**: The Python environment is managed by `uv` with pinned dependencies in `pyproject.toml`. The Jupyter kernel is registered from this environment.
3. **Rendering reproducibility**: Quarto's `freeze` mechanism stores computed outputs alongside their source code. Collaborators can render the complete manuscript without installing Python or running Neo4j — they only need the Quarto CLI.
4. **Data provenance**: The paper repository includes the complete KG build pipeline, so readers can trace any result back through the adapter code to the original data source (NCBI, UniProt, Cyanorak, or the specific publication's supplementary table).

## Limitations

- **AI hallucination risk**: The AI assistant could generate plausible but incorrect biological claims. Every factual statement required verification by the domain expert, particularly in the Introduction and Discussion where the assistant drew on its general training data rather than project-specific files.
- **Context window constraints**: Long manuscripts can exceed the AI assistant's context window. The section-by-section approach mitigated this by focusing each session on a manageable portion of the text.
- **Figure aesthetics**: While the AI assistant produced functional figures, fine-tuning visual aesthetics (color choices, label placement, panel arrangement) often required multiple rounds of feedback.
- **Citation accuracy**: The assistant could suggest relevant references from its training data, but DOIs and publication details required verification against the actual bibliography.

## Assessment

The Quarto + Claude Code workflow proved effective for this type of paper — a bioinformatics resource paper where the majority of content describes a system the AI assistant has full access to. The assistant's ability to read the codebase, query the knowledge graph, and write both prose and code within a single framework significantly accelerated the writing process. The tight Quarto render loop (edit → render in <1s → review) made iterative refinement practical in a way that would not have been possible with slower LaTeX compilation.

We note that this workflow is best suited for papers where the AI assistant has access to the underlying data and code. For papers requiring synthesis of external literature or novel biological interpretation, the domain expert's role remains central, with the AI assistant serving primarily as a drafting and formatting accelerator.
