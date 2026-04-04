# multiomics_kg/extraction/cluster/prompts.py
"""All LLM prompts for cluster description extraction.

Separated from code so prompts can be iterated independently.
Each prompt is a module-level string constant.
"""

# ── Extraction fields (shared across all prompts) ──

EXTRACTION_FIELDS_DESCRIPTION = """\
Extract ALL of the following fields for each cluster. If a field does not apply \
or you cannot determine it with confidence, write "not described in paper".
Do NOT guess. "not described in paper" is always better than an incorrect value.

Fields:
- enrichment_category: functional category enriched in this cluster
- enrichment_pvalue: p-value or FDR for that enrichment
- enrichment_significant: true/false
- enrichment_details: additional enrichment info (GO terms, specific subcategories)
- genes: specific gene names attributed to this cluster (comma-separated)
- direction: "up" / "down" / "periodic"
- cluster_description: general description of what this cluster represents
- temporal_pattern: when and how expression changes
- peak_time: peak expression time for diel/periodic clusters
- period_description: oscillation behavior description
- light_phase: "light" / "dark" / "transition" / "not described in paper"
- supporting_quotes: ALL exact quotes that describe this cluster (no limit)
"""

# ── Path: visual ──

VISUAL_PROMPT = """\
You are extracting information about gene expression clusters from a scientific \
paper and its supplementary materials.

The paper reports {n_clusters} clusters from {cluster_method} analysis of \
{organism} under {treatment}. The cluster IDs in the data are: {cluster_keys}.

{fields_description}

For each cluster, also extract:
- supporting_quotes: exact quotes from the paper with figure/table location

IMPORTANT:
- Use the exact cluster IDs ({cluster_keys}) as JSON keys.
- Read figure legends carefully — they often have enrichment categories and p-values.
- Be precise about which cluster information refers to.
- Only extract data for {organism} (not other organisms in the same paper).
- If you cannot determine a field, use "not described in paper" — do NOT guess.

Return valid JSON:
{{
  "<cluster_id>": {{
    "enrichment_category": "...",
    "enrichment_pvalue": "...",
    "enrichment_significant": true,
    "enrichment_details": "...",
    "genes": "...",
    "direction": "up|down|periodic",
    "cluster_description": "...",
    "temporal_pattern": "...",
    "peak_time": "...",
    "period_description": "...",
    "light_phase": "...",
    "supporting_quotes": [{{"quote": "...", "location": "Figure 3 legend"}}]
  }}
}}
"""

# ── Path: semantic ──

SEMANTIC_PROMPT = """\
You are extracting information about gene cluster {cluster_key} from a \
scientific paper about {organism} under {treatment}.

Below are text passages retrieved from the paper and supplementary materials. \
Extract ONLY information that is explicitly about cluster {cluster_key}. \
If a passage discusses multiple clusters, only use parts about THIS cluster.

If you cannot determine something with confidence, use "not described in paper". \
Do NOT guess. "not described in paper" is always better than wrong.

{fields_description}

Each supporting_quote must include a relevance_score (provided with the passage).

Return valid JSON:
{{
  "enrichment_category": "...",
  "enrichment_pvalue": "...",
  "enrichment_significant": true,
  "enrichment_details": "...",
  "genes": "...",
  "direction": "up|down|periodic",
  "cluster_description": "...",
  "temporal_pattern": "...",
  "peak_time": "...",
  "period_description": "...",
  "light_phase": "...",
  "supporting_quotes": [{{"quote": "...", "location": "...", "relevance_score": 0.85}}]
}}

Passages:
{passages}
"""

# ── Stage 2: Synthesis ──

SYNTHESIS_PROMPT = """\
You are writing descriptions for gene expression clusters from a scientific paper.

For each cluster below, you are given extracted data from multiple sources \
(table=structured data files, visual=PDF figures, semantic=paper text). \
Each field is a list of extractions tagged with source and confidence level \
(very_high > high > medium > low).

Write these outputs per cluster:

1. **id** — Short snake_case identifier prefixed with organism. \
Format: {{organism}}_{{direction}}_{{theme}}. \
Examples: med4_up_n_transport, mit9313_down_translation, med4_diel_dawn_psi. \
Must be unique across all clusters in this paper.

2. **name** — Short human-readable cluster label. \
Format: "{{Organism}} cluster {{key}} ({{direction}}, {{theme}})". \
Examples: "MED4 cluster 1 (up, N transport)", "MIT9313 cluster 3 (down, translation)". \
Keep under 60 characters.

3. **functional_description** — What types of genes are in this cluster. \
2-3 sentences max. Include enrichment categories/pathways if mentioned in the paper. \
For p-values use max 3 decimal places or scientific notation (e.g., p=0.013 or p=7.9e-10). \
List only genes mentioned BY NAME in the paper text (not from CSV), max 3-5 genes. \
Do NOT include treatment/experimental conditions (those belong on the analysis, not clusters).

4. **behavioral_description** — The temporal/response pattern. \
1-2 sentences max. Include direction (up/down), timing, and magnitude if available. \
Do NOT repeat treatment conditions or experimental setup.

5. **peak_time_hours** — Peak expression time in hours (float). \
Only for diel/periodic clusters. null for non-periodic clusters.

6. **period_hours** — Oscillation period in hours (float). \
Only for periodic clusters. null for non-periodic clusters.

RULES:
- Descriptions must be SHORT: functional_description 2-3 sentences, behavioral_description 1-2 sentences.
- Never leave fields empty. Use these sentinel values:
  - "not described in paper" — paper doesn't discuss this aspect
  - "insufficient data" — sources found but too ambiguous
  - "conflicting sources: [source A says X, source B says Y]" — paths disagree
- Partial descriptions are fine; incorrect descriptions are not.
- Prefer higher-confidence sources when they agree on meaning.
- If uncertain, use the sentinel value — do NOT guess.
- IDs must be unique within this paper (organism prefix helps).
- peak_time_hours and period_hours must be null (not "null") for non-periodic clusters.
- Do NOT include treatment conditions or experimental context in cluster descriptions.
- Only cite genes that are mentioned by name in the paper text, not from data tables.
- Maximum 3-5 named genes per cluster description.

Paper: {paper_name}
Organism: {organism}
Treatment: {treatment}
Cluster method: {cluster_method}

{cluster_blocks}

Return valid JSON using the original cluster keys:
{{
  "<cluster_key>": {{
    "id": "organism_direction_theme",
    "name": "Organism cluster N (direction, theme)",
    "functional_description": "...",
    "behavioral_description": "...",
    "peak_time_hours": null,
    "period_hours": null
  }}
}}
"""

# ── Stage 3: Validation (LLM-as-judge) ──

JUDGE_PROMPT = """\
You are validating extracted gene cluster descriptions against the original \
paper and data.

You can see: the paper PDF, the cluster CSV data, and the extracted descriptions.
Be strict — flag anything not clearly supported by the paper or data.

For each cluster, check:
1. enrichment_correct — Does the paper support this enrichment category? (yes/no/uncertain)
2. genes_correct — Are the listed genes attributed to this cluster in the paper? (yes/no/uncertain)
3. direction_correct — Is the up/down/periodic direction correct? (yes/no)
4. behavioral_correct — Is the temporal description supported? (yes/no/uncertain)
5. hallucination — Any claims NOT in the paper? (yes/no — if yes, list them)
6. missing_info — Important paper info that is missing? (yes/no — if yes, what)
7. verdict — "pass" (all correct), "warn" (minor issues), "fail" (factual errors)
8. explanation — Brief explanation of any issues found

Sentinel values ("not described in paper", "insufficient data") are EXPECTED \
for clusters the paper doesn't discuss in detail. Only "warn" for missing info \
if the paper actually describes something the extraction missed. Do NOT "fail" \
for fields the paper genuinely doesn't cover.

Extracted descriptions to validate:
{descriptions}

Return valid JSON using the cluster keys:
{{
  "<cluster_key>": {{
    "enrichment_correct": "yes|no|uncertain",
    "genes_correct": "yes|no|uncertain",
    "direction_correct": "yes|no",
    "behavioral_correct": "yes|no|uncertain",
    "hallucination": "no",
    "missing_info": "no",
    "verdict": "pass|warn|fail",
    "explanation": "..."
  }}
}}
"""
