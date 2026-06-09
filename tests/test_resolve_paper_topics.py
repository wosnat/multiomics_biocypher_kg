"""Tests for Stage 2 topic resolution (deterministic functions)."""
import pytest

from multiomics_kg.utils.gene_id_utils import MappingData
from multiomics_kg.download.resolve_paper_topics import (
    build_pathway_lookup,
    normalize_kegg_pathway_id,
    resolve_gene_mention,
    resolve_pathway,
    _match_strain,
)


# ── KEGG pathway normalization + resolution ──


@pytest.mark.parametrize("raw,expected", [
    ("ko00910", "ko00910"),
    ("map00910", "ko00910"),
    ("path:map00910", "ko00910"),
    ("00910", "ko00910"),
    ("ko00010", "ko00010"),
    ("", None),
    (None, None),
    ("photosynthesis", None),     # no number
    ("123", None),                # too few digits
])
def test_normalize_kegg_pathway_id(raw, expected):
    assert normalize_kegg_pathway_id(raw) == expected


@pytest.fixture
def kegg_data():
    return {"pathways": {
        "ko00910": {"name": "Nitrogen metabolism"},
        "ko00010": {"name": "Glycolysis / Gluconeogenesis"},
    }}


def test_build_pathway_lookup(kegg_data):
    id_lookup, name_lookup = build_pathway_lookup(kegg_data)
    assert id_lookup["ko00910"] == "ko00910"
    assert name_lookup["nitrogen metabolism"] == "ko00910"
    assert name_lookup["glycolysis / gluconeogenesis"] == "ko00010"


def test_resolve_pathway_id_first(kegg_data):
    id_lookup, name_lookup = build_pathway_lookup(kegg_data)
    ko, name, method = resolve_pathway(
        {"kegg_pathway_id": "map00910", "surface_form": "N metab"}, id_lookup, name_lookup)
    assert (ko, method) == ("ko00910", "id")
    assert name == "nitrogen metabolism"


def test_resolve_pathway_name_fallback(kegg_data):
    id_lookup, name_lookup = build_pathway_lookup(kegg_data)
    ko, _name, method = resolve_pathway(
        {"kegg_pathway_id": "", "surface_form": "Nitrogen Metabolism"}, id_lookup, name_lookup)
    assert (ko, method) == ("ko00910", "name")


def test_resolve_pathway_unresolved_drops(kegg_data):
    """A pathway not in the pruned set (e.g. ko04110) yields no node -> drop."""
    id_lookup, name_lookup = build_pathway_lookup(kegg_data)
    ko, _name, method = resolve_pathway(
        {"kegg_pathway_id": "ko04110", "surface_form": "cell cycle"}, id_lookup, name_lookup)
    assert ko is None and method == "unresolved"


# ── Strain matching ──


def test_match_strain():
    strains = ["Prochlorococcus MED4", "Prochlorococcus MIT9313"]
    assert _match_strain("Prochlorococcus MED4", strains) == "Prochlorococcus MED4"
    assert _match_strain("prochlorococcus med4", strains) == "Prochlorococcus MED4"  # ci
    assert _match_strain("all", strains) is None
    assert _match_strain("", strains) is None


# ── Gene mention resolution (strain-aware) ──


@pytest.fixture
def med4_mapping():
    return MappingData(
        specific_lookup={"ntcA": "PMM0246", "PMM0552": "PMM0552", "psbA": "PMM0223"},
        multi_lookup={}, conflicts={},
        locus_tags={"PMM0246", "PMM0552", "PMM0223"}, version=2,
    )


@pytest.fixture
def mit9313_mapping():
    return MappingData(
        specific_lookup={"ntcA": "PMT1831", "PMT0956": "PMT0956"},
        multi_lookup={}, conflicts={},
        locus_tags={"PMT1831", "PMT0956"}, version=2,
    )


def test_gene_identifier_first(med4_mapping):
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "X", "gene_name": "", "identifiers": ["PMM0552"],
         "strain": "Prochlorococcus MED4", "prominence": "central"},
        ["Prochlorococcus MED4"], cache)
    assert len(recs) == 1
    assert recs[0]["locus_tag"] == "PMM0552"
    assert recs[0]["resolved_from"] == "PMM0552"


def test_gene_family_multi_identifier_fans_out(med4_mapping):
    """A family mention with N distinct member ids yields N records."""
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "FAM", "gene_name": "", "identifiers": ["PMM0552", "PMM0223"],
         "strain": "Prochlorococcus MED4", "prominence": "central"},
        ["Prochlorococcus MED4"], cache)
    assert {r["locus_tag"] for r in recs} == {"PMM0552", "PMM0223"}
    assert all(r["strain"] == "Prochlorococcus MED4" for r in recs)


def test_gene_synonym_identifiers_collapse(med4_mapping):
    """Repeated ids resolving to the same locus_tag collapse to one record."""
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "X", "gene_name": "", "identifiers": ["PMM0552", "PMM0552"],
         "strain": "Prochlorococcus MED4", "prominence": "central"},
        ["Prochlorococcus MED4"], cache)
    assert len(recs) == 1 and recs[0]["locus_tag"] == "PMM0552"


def test_gene_name_skipped_when_identifier_resolves(med4_mapping):
    """gene_name is a fallback only — not added when an identifier already resolved."""
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "X", "gene_name": "ntcA", "identifiers": ["PMM0552"],
         "strain": "Prochlorococcus MED4", "prominence": "central"},
        ["Prochlorococcus MED4"], cache)
    # ntcA would resolve to PMM0246, but the PMM0552 identifier wins; name skipped.
    assert {r["locus_tag"] for r in recs} == {"PMM0552"}


def test_gene_name_specific_strain(med4_mapping):
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "NtcA", "gene_name": "ntcA", "identifiers": [],
         "strain": "Prochlorococcus MED4", "prominence": "peripheral"},
        ["Prochlorococcus MED4"], cache)
    assert [r["locus_tag"] for r in recs] == ["PMM0246"]


def test_gene_all_expands_per_strain(med4_mapping, mit9313_mapping):
    """A gene_name with strain='all' resolves in every paper strain that has it."""
    cache = {"Prochlorococcus MED4": med4_mapping,
             "Prochlorococcus MIT9313": mit9313_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "NtcA", "gene_name": "ntcA", "identifiers": [],
         "strain": "all", "prominence": "peripheral"},
        ["Prochlorococcus MED4", "Prochlorococcus MIT9313"], cache)
    by_strain = {r["strain"]: r["locus_tag"] for r in recs}
    assert by_strain == {"Prochlorococcus MED4": "PMM0246",
                         "Prochlorococcus MIT9313": "PMT1831"}


def test_gene_unresolved_single_record(med4_mapping):
    cache = {"Prochlorococcus MED4": med4_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "Ferredoxin", "gene_name": "", "identifiers": [],
         "strain": "Prochlorococcus MED4", "prominence": "peripheral"},
        ["Prochlorococcus MED4"], cache)
    assert len(recs) == 1
    assert recs[0]["locus_tag"] is None


def test_gene_all_unresolved_yields_one_audit_record(med4_mapping, mit9313_mapping):
    cache = {"Prochlorococcus MED4": med4_mapping,
             "Prochlorococcus MIT9313": mit9313_mapping}
    recs = resolve_gene_mention(
        {"surface_form": "Mystery", "gene_name": "nope", "identifiers": [],
         "strain": "all", "prominence": "peripheral"},
        ["Prochlorococcus MED4", "Prochlorococcus MIT9313"], cache)
    assert len(recs) == 1
    assert recs[0]["locus_tag"] is None
    assert recs[0]["strain"] == "all"
