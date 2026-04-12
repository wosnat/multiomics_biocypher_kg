"""Prompt construction for the timepoint/growth-phase extractor."""
from __future__ import annotations

import json

VALID_GROWTH_PHASES_LIST = [
    "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
    "infected", "recovery", "diel", "darkness", "death", "acute_stress", "unknown",
]


SHARED_RULES = f"""\
You are extracting sampling-time and growth-phase metadata for a list of
statistical analyses from one scientific publication. The publication's PDF
is attached; additional supplementary PDFs may also be attached.

For each analysis listed in the "Extraction targets" section below, return
the fields listed in its `fields_requested` array — and no others.

FIELD VALUE RULES

1. `timepoint` (string, human-readable): the sampling time label, e.g.
   "120 min", "24h", "acclimated steady-state". Use "unknown" ONLY if the
   paper truly doesn't state a time.

2. `timepoint_hours` (number or null): numeric conversion of `timepoint`
   to hours. null ONLY if the paper truly doesn't state a time. "120 min" →
   2. "48h" → 48. "acclimated" with no specific hours → null.

3. `growth_phase` (string): the physiological state at sampling. Must be
   one of these canonical values:
     {", ".join(VALID_GROWTH_PHASES_LIST)}

   If the paper describes a phase not in this list, emit `other:<slug>`
   with a short snake_case slug, e.g. `other:heat_acclimated`,
   `other:late_decline`. Prefer `other:<slug>` over `unknown` whenever
   the paper gives any positional information.

   `unknown` means the paper gives no information about physiological
   state at sampling. `other:<slug>` means the paper describes a state
   that doesn't fit the enum.

QUALITY AND ANCHORING RULES

4. `self_assessment`: "high", "medium", or "low". Based on how directly
   the paper states the value and how unambiguous the evidence is.

5. `assessment_notes`: free-text rationale. Short, specific. Use this to
   flag edge cases, disagreements between the CSV label and the paper
   methods, or anything a reviewer should know.

6. `supporting_quotes`: list of objects {{"quote": "...", "location": "..."}}.
   Verbatim from the paper's methods or results. At least one quote is
   required unless `self_assessment` is "low" AND `assessment_notes`
   explains why no direct quote is available.

7. `source_figures`: list of figure or table numbers you consulted.
   Empty list is acceptable for purely text-derived extractions.

CROSS-CHECK RULE

8. Each analysis has a `logfc_col` — the fold-change column name from
   the source CSV. This column name often contains the timepoint
   verbatim ("24 hours", "Axenic, 36 hours", "log2FC_48h_vs_T0"). Cross-check:

   - If `logfc_col` agrees with paper methods, cite both in
     `supporting_quotes`.
   - If they disagree, PREFER paper methods but note the discrepancy
     in `assessment_notes` and lower `self_assessment` accordingly.

DEFAULT BIAS

9. Unknown is always better than wrong. If you cannot find clear evidence
   from the paper, set `growth_phase` to "unknown" and explain in
   `assessment_notes`.

OUTPUT SHAPE

Return one JSON object:

```
{{
  "analyses": [
    {{
      "analysis_id": "...",
      "timepoint": "...",        // only if in fields_requested
      "timepoint_hours": 2.0,    // only if in fields_requested; may be null
      "growth_phase": "...",     // only if in fields_requested
      "self_assessment": "high", // always
      "assessment_notes": "...", // always (may be empty string)
      "supporting_quotes": [     // always
        {{"quote": "...", "location": "..."}}
      ],
      "source_figures": []       // always (may be empty list)
    }}
  ]
}}
```

Do NOT emit `metadata`, `experiment_key`, or `fields_requested` — these
are added by the extraction tool after your response.

PER-EXPERIMENT-TYPE HINTS

- Diel studies (treatment_type contains "diel") → `diel` for
  cycling-phase samples. Do not split into light/dark sub-phases.
- Extended darkness (treatment_type contains "darkness", prolonged
  dark exposure of hours or days) → `darkness`. Do NOT confuse with
  the dark half of a normal diel cycle — that stays `diel`.
- Acclimation / chronic exposure (e.g. treatment_type contains
  "carbon" or "temperature" with language about ≥5 generations) →
  `acclimated_steady_state`.
- Nutrient starvation (treatment_type contains "nitrogen", "phosphorus",
  or "iron") → `exponential` at early timepoints before depletion is
  confirmed, `nutrient_limited` once cells are clearly nutrient-starved,
  `death` if sampling extends past viable-cell peak.
- Phage / viral infection → `exponential` at t=0 (pre-infection),
  `infected` at post-infection timepoints.
- Rescue / re-addition experiments → `recovery` for post-intervention
  timepoints.
- Short-exposure stress (≤6 h) at still-dividing cells → `acute_stress`,
  ONLY when no more-specific phase fits. `nutrient_limited`, `infected`,
  and `darkness` take precedence.
"""


def build_prompt(
    background: dict,
    targets: list[dict],
    pdf_cache_entry: dict | None = None,
) -> str:
    """Assemble the text prompt. PDFs are attached separately by extract.py."""
    sections = [SHARED_RULES]

    sections.append("\n## PAPER CONTEXT")
    sections.append(f"**Paper:** {background.get('papername', '')}")
    sections.append(f"**DOI:** {background.get('doi', '')}")

    if pdf_cache_entry:
        sections.append("\n### Publication metadata")
        for k in ("title", "abstract", "description", "study_type"):
            if v := pdf_cache_entry.get(k):
                sections.append(f"**{k.title()}:** {v}")

    sections.append("\n### Experiments block")
    sections.append("```yaml")
    sections.append(json.dumps(background.get("experiments", {}), indent=2, default=str))
    sections.append("```")

    sections.append("\n## EXTRACTION TARGETS")
    if not targets:
        sections.append("_(no analyses require extraction — all fields populated)_")
    else:
        sections.append("```json")
        sections.append(json.dumps(targets, indent=2, default=str))
        sections.append("```")

    return "\n".join(sections)
