"""Calls.json <-> gene_annotations_merged.json consistency (spec 2026-08-17 §4.2/§6).

The two ungated passthrough fields must agree exactly per gene:
merged interpro_entries == calls interpro_entries keys, and
merged ncbifam_ids == calls libraries.NCBIFAM accessions.
A failure means calls.json was re-normalized without re-running
`prepare_data --steps 2` (or vice versa) — fix by re-running the merge.

A strain whose merged JSON predates the ncbifam_ids-era merge is SKIPPED
(not failed) with a re-run-the-merge hint; see the skip guard below.
"""
import csv
import json
from pathlib import Path

import pytest

GENOMES_CSV = Path("data/Prochlorococcus/genomes/cyanobacteria_genomes.csv")


def _strains():
    with open(GENOMES_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.lstrip().startswith("#")))
    return [(r["strain_name"], Path(r["data_dir"])) for r in rows if r.get("data_dir")]


def _calls_has_ncbifam(calls: dict) -> bool:
    """True if any protein in calls.json carries an NCBIFAM library facet."""
    for call in calls.values():
        if not isinstance(call, dict):
            continue
        if (call.get("libraries") or {}).get("NCBIFAM"):
            return True
    return False


def _merged_has_ncbifam_marker(merged: dict) -> bool:
    """True if any gene in the merged JSON has an ncbifam_ids key at all.

    Presence of the key (even as an empty list) marks a merge that ran
    after Task 9's ncbifam_ids passthrough was added. Its total absence
    across every gene marks a stale, pre-Task-9 merge.
    """
    return any("ncbifam_ids" in gene for gene in merged.values())


@pytest.mark.parametrize("strain,data_dir", _strains(), ids=lambda v: str(v))
def test_calls_merge_consistency(strain, data_dir):
    calls_p = data_dir / "interproscan" / f"{strain}.interproscan.calls.json"
    merged_p = data_dir / "gene_annotations_merged.json"
    if not calls_p.exists() or not merged_p.exists():
        pytest.skip("artifact missing")
    calls = json.loads(calls_p.read_text())
    merged = json.loads(merged_p.read_text())

    if _calls_has_ncbifam(calls) and not _merged_has_ncbifam_marker(merged):
        pytest.skip(
            "merged JSON predates ncbifam_ids — re-run prepare_data --steps 2"
        )

    for locus_tag, gene in merged.items():
        wp = (gene.get("protein_id") or "").strip()
        call = calls.get(wp)
        if not call:
            assert not gene.get("interpro_entries"), (
                f"{strain}:{locus_tag} merged has entries, calls lack protein"
            )
            continue
        assert sorted(gene.get("interpro_entries") or []) == sorted(
            call.get("interpro_entries") or {}
        ), f"{strain}:{locus_tag} interpro_entries skew — re-run prepare_data --steps 2"
        calls_ncbifam = sorted(
            {r["accession"] for r in (call.get("libraries") or {}).get("NCBIFAM", [])}
        )
        assert sorted(gene.get("ncbifam_ids") or []) == calls_ncbifam, (
            f"{strain}:{locus_tag} ncbifam_ids skew — re-run prepare_data --steps 2"
        )
