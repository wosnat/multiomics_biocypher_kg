"""Pure-Python SignalP 6.0 signal-peptide vocabulary + parser.

SignalP 6.0 emits one of a closed 6-value set per protein. Five are real
signal-peptide types; ``OTHER`` is the no-signal sentinel (assigned when no
peptide class scores above SignalP's threshold). The KG models the five real
types as ``SignalPeptideType`` ontology nodes and skips ``OTHER`` — absence of a
``Gene_has_signal_peptide_type`` edge encodes "no signal peptide".

This module is intentionally file-I/O-free so it is unit-testable in isolation.
``signalp-run --normalize`` reads each strain's raw ``prediction_results.txt``
and calls :func:`parse_prediction_results` to produce the normalized
``<strain>.signalp.calls.json`` (the ``add-a-tool`` Phase-1 artifact convention,
which ``signalp-run`` predates).

Public API:
    SIGNALP_VOCAB  — dict of raw code → display name (real types only).
    OTHER_SENTINEL — the literal string SignalP uses for no-signal calls.
    PROB_COLUMNS   — fixed per-class probability column order in the raw TSV.
    is_kept(call)            — True when ``call`` is a real (non-sentinel) type.
    display_name(call)       — human-readable name for a real type.
    first_token(header)      — RefSeq WP_ accession = first whitespace token.
    parse_cs_pos(cs)         — "CS pos: 21-22. Pr: 0.3982" → (21, 0.3982).
    parse_prediction_results(text) — full TSV → {wp_accession: record}.
"""
from __future__ import annotations

import re

# The no-signal sentinel — kept verbatim in gene_annotations_merged.json so the
# field is round-trippable, but NOT emitted as a node or edge.
OTHER_SENTINEL = "OTHER"

# Closed vocabulary of the FIVE real SignalP 6.0 signal-peptide types.
# code → human-readable display name (used as the node ``name`` property).
SIGNALP_VOCAB: dict[str, str] = {
    "SP": "Signal peptide (Sec/SPI)",
    "LIPO": "Lipoprotein signal peptide (Sec/SPII)",
    "TAT": "TAT signal peptide (Tat/SPI)",
    "TATLIPO": "TAT lipoprotein signal peptide (Tat/SPII)",
    "PILIN": "Pilin-like signal peptide (Sec/SPIII)",
}

# Per-class probability columns in prediction_results.txt, in the fixed order
# SignalP 6.0 writes them (header:
#   # ID  Prediction  OTHER  SP(Sec/SPI)  LIPO(Sec/SPII)  TAT(Tat/SPI)
#         TATLIPO(Tat/SPII)  PILIN(Sec/SPIII)  CS Position
# ). Probabilities start at column index 2; the winning class's probability is
# PROB_COLUMNS.index(code) positions in.
PROB_COLUMNS: list[str] = ["OTHER", "SP", "LIPO", "TAT", "TATLIPO", "PILIN"]

_CS_SITE_RE = re.compile(r"pos[.:]?\s*(\d+)")
_CS_PROB_RE = re.compile(r"(?:Pr|Probability)[.:]?\s*([\d.]+)")


def is_kept(call: str | None) -> bool:
    """True when ``call`` is a real signal-peptide type the KG keeps."""
    return bool(call) and call in SIGNALP_VOCAB


def display_name(call: str) -> str:
    """Human-readable display name for a real type.

    Falls back to the raw code for anything outside the vocabulary (defensive —
    callers should gate on :func:`is_kept` first).
    """
    return SIGNALP_VOCAB.get(call, call)


def first_token(header: str) -> str:
    """RefSeq WP_ accession = first whitespace token of a FASTA header.

    The prediction_results.txt ID column carries the full FASTA header
    (e.g. ``"WP_011131644.1 MULTISPECIES: ... [Prochlorococcus]"``); the join
    key into gene_mapping.csv is the leading accession.
    """
    return header.split()[0] if header and header.strip() else ""


def parse_cs_pos(cs: str | None) -> tuple[int | None, float | None]:
    """Parse a SignalP cleavage-site string into (site, probability).

    Handles both the TSV form (``"CS pos: 21-22. Pr: 0.3982"``) and the
    verbose output.json form (``"Cleavage site between pos. 21 and 22.
    Probability 0.398243"``). ``site`` is the residue *before* the cut (the
    first position number). Returns ``(None, None)`` for an empty/absent CS.
    """
    if not cs or not cs.strip():
        return None, None
    site: int | None = None
    prob: float | None = None
    m = _CS_SITE_RE.search(cs)
    if m:
        site = int(m.group(1))
    m = _CS_PROB_RE.search(cs)
    if m:
        try:
            prob = float(m.group(1))
        except ValueError:
            prob = None
    return site, prob


def parse_prediction_results(text: str) -> dict[str, dict]:
    """Parse a SignalP 6.0 ``prediction_results.txt`` into a WP_-keyed dict.

    Each value is ``{signalp_type, probability, cleavage_site,
    cleavage_probability}``:
    - ``signalp_type`` — winning code (``SP`` / ``LIPO`` / ``TAT`` / ``TATLIPO``
      / ``PILIN`` / ``OTHER``); ``OTHER`` rows are kept verbatim (the adapter
      skips them at edge-build time, mirroring PSORTb's ``Unknown``).
    - ``probability`` — the winning class's likelihood (∈[0,1]).
    - ``cleavage_site`` / ``cleavage_probability`` — parsed from the CS column;
      ``None`` for ``OTHER`` and for kept types with no reported cleavage site.

    Pure: takes the file *text*, not a path. Malformed / short rows are skipped.
    """
    records: dict[str, dict] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        # ID, Prediction, 6 probabilities = 8 minimum; CS column optional.
        if len(parts) < 2 + len(PROB_COLUMNS):
            continue
        wp = first_token(parts[0])
        if not wp:
            continue
        code = parts[1].strip()
        # Winning class probability column.
        probability: float | None = None
        if code in PROB_COLUMNS:
            try:
                probability = float(parts[2 + PROB_COLUMNS.index(code)].strip())
            except (ValueError, IndexError):
                probability = None
        cs = parts[8] if len(parts) > 8 else ""
        site, cs_prob = parse_cs_pos(cs)
        records[wp] = {
            "signalp_type": code,
            "probability": probability,
            "cleavage_site": site,
            "cleavage_probability": cs_prob,
        }
    return records
