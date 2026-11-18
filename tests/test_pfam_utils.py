"""
Unit tests for multiomics_kg/utils/pfam_utils.py.

Tests cover:
- parse_clans_tsv: parsing gzipped TSV into PfamData
- Entries with and without clan membership
- JSON round-trip serialization (_pfam_data_to_json / _pfam_data_from_json)
- No network calls (all fixtures use local temp files)
"""

import gzip
import json
import pytest
from pathlib import Path

from multiomics_kg.utils.pfam_utils import (
    PfamData,
    PfamEntry,
    parse_clans_tsv,
    _pfam_data_to_json,
    _pfam_data_from_json,
    parse_dead_families,
    normalize_pfam_token,
    parse_release,
)


# ---------------------------------------------------------------------------
# Test data: a small Pfam-A.clans.tsv content
# ---------------------------------------------------------------------------

# 5 columns: accession, clan_accession, clan_name, shortname, description
# Entries with clan, entries without clan, entries with special chars
MINI_CLANS_TSV = """\
PF00712\tCL0060\tDNA_clamp\tDNA_pol3_beta\tDNA polymerase III beta subunit, N-terminal domain
PF02768\tCL0060\tDNA_clamp\tDNA_pol3_beta_2\tDNA polymerase III beta subunit, central domain
PF02767\tCL0060\tDNA_clamp\tDNA_pol3_beta_3\tDNA polymerase III beta subunit, C-terminal domain
PF00001\tCL0192\tGPCR_A\t7tm_1\t7 transmembrane receptor (rhodopsin family)
PF13927\t\t\t3H\t3H domain
PF00069\t\t\tPkinase\tProtein kinase domain
PF99999\tCL9999\tTest_clan\tTest_entry\tEntry with 'quotes' and | pipes
"""


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def clans_gz(tmp_path):
    """Write MINI_CLANS_TSV as a gzipped file and return its path."""
    gz_path = tmp_path / "Pfam-A.clans.tsv.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
        fh.write(MINI_CLANS_TSV)
    return gz_path


@pytest.fixture
def pfam_data(clans_gz):
    """Parse the mini clans TSV and return PfamData."""
    return parse_clans_tsv(clans_gz)


# ===========================================================================
# Tests for parse_clans_tsv
# ===========================================================================


class TestParseClansBasic:
    def test_returns_pfam_data(self, pfam_data):
        assert isinstance(pfam_data, PfamData)

    def test_all_accessions_present(self, pfam_data):
        expected = {"PF00712", "PF02768", "PF02767", "PF00001", "PF13927", "PF00069", "PF99999"}
        assert set(pfam_data.by_accession.keys()) == expected

    def test_entry_count(self, pfam_data):
        assert len(pfam_data.by_accession) == 7

    def test_entry_is_pfam_entry(self, pfam_data):
        entry = pfam_data.by_accession["PF00712"]
        assert isinstance(entry, PfamEntry)


class TestParseClansWithClan:
    def test_clan_accession_populated(self, pfam_data):
        entry = pfam_data.by_accession["PF00712"]
        assert entry.clan_accession == "CL0060"

    def test_clan_name_populated(self, pfam_data):
        entry = pfam_data.by_accession["PF00712"]
        assert entry.clan_name == "DNA_clamp"

    def test_shortname_populated(self, pfam_data):
        entry = pfam_data.by_accession["PF00712"]
        assert entry.shortname == "DNA_pol3_beta"

    def test_description_populated(self, pfam_data):
        entry = pfam_data.by_accession["PF00712"]
        assert entry.description == "DNA polymerase III beta subunit, N-terminal domain"

    def test_clan_in_clans_dict(self, pfam_data):
        assert "CL0060" in pfam_data.clans
        assert pfam_data.clans["CL0060"] == "DNA_clamp"

    def test_multiple_entries_same_clan(self, pfam_data):
        """All three DNA_pol3_beta domains share CL0060."""
        for pf_id in ("PF00712", "PF02768", "PF02767"):
            assert pfam_data.by_accession[pf_id].clan_accession == "CL0060"

    def test_second_clan_present(self, pfam_data):
        assert "CL0192" in pfam_data.clans
        assert pfam_data.clans["CL0192"] == "GPCR_A"


class TestParseClansWithoutClan:
    def test_no_clan_accession(self, pfam_data):
        entry = pfam_data.by_accession["PF13927"]
        assert entry.clan_accession == ""

    def test_no_clan_name(self, pfam_data):
        entry = pfam_data.by_accession["PF13927"]
        assert entry.clan_name == ""

    def test_not_in_clans_dict(self, pfam_data):
        # PF13927 and PF00069 have no clan; their empty clan_accession
        # should NOT appear in the clans dict
        assert "" not in pfam_data.clans

    def test_shortname_still_present(self, pfam_data):
        entry = pfam_data.by_accession["PF00069"]
        assert entry.shortname == "Pkinase"

    def test_description_still_present(self, pfam_data):
        entry = pfam_data.by_accession["PF00069"]
        assert entry.description == "Protein kinase domain"


class TestParseClansShortnameLookup:
    def test_shortname_to_accession(self, pfam_data):
        assert pfam_data.by_shortname["DNA_pol3_beta"] == "PF00712"

    def test_all_shortnames_present(self, pfam_data):
        expected_shortnames = {
            "DNA_pol3_beta", "DNA_pol3_beta_2", "DNA_pol3_beta_3",
            "7tm_1", "3H", "Pkinase", "Test_entry",
        }
        assert set(pfam_data.by_shortname.keys()) == expected_shortnames

    def test_shortname_count_matches_entry_count(self, pfam_data):
        assert len(pfam_data.by_shortname) == len(pfam_data.by_accession)

    def test_reverse_lookup_consistency(self, pfam_data):
        """Every shortname maps to an accession that maps back to that shortname."""
        for shortname, acc in pfam_data.by_shortname.items():
            entry = pfam_data.by_accession[acc]
            assert entry.shortname == shortname


class TestParseClansCount:
    def test_clan_count(self, pfam_data):
        # CL0060 (DNA_clamp), CL0192 (GPCR_A), CL9999 (Test_clan)
        assert len(pfam_data.clans) == 3


class TestParseClansEmptyLines:
    def test_empty_lines_skipped(self, tmp_path):
        """Empty lines in the TSV should be silently skipped."""
        content = "\nPF00712\tCL0060\tDNA_clamp\tDNA_pol3_beta\tSome description\n\n"
        gz_path = tmp_path / "test.tsv.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
            fh.write(content)
        data = parse_clans_tsv(gz_path)
        assert len(data.by_accession) == 1

    def test_short_lines_skipped(self, tmp_path):
        """Lines with fewer than 5 columns should be silently skipped."""
        content = "PF00712\tCL0060\tDNA_clamp\n"  # only 3 columns
        gz_path = tmp_path / "test.tsv.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
            fh.write(content)
        data = parse_clans_tsv(gz_path)
        assert len(data.by_accession) == 0


# ===========================================================================
# Tests for JSON round-trip
# ===========================================================================


class TestJsonRoundTrip:
    def test_roundtrip_preserves_accessions(self, pfam_data):
        json_dict = _pfam_data_to_json(pfam_data)
        restored = _pfam_data_from_json(json_dict)
        assert set(restored.by_accession.keys()) == set(pfam_data.by_accession.keys())

    def test_roundtrip_preserves_shortnames(self, pfam_data):
        json_dict = _pfam_data_to_json(pfam_data)
        restored = _pfam_data_from_json(json_dict)
        assert restored.by_shortname == pfam_data.by_shortname

    def test_roundtrip_preserves_clans(self, pfam_data):
        json_dict = _pfam_data_to_json(pfam_data)
        restored = _pfam_data_from_json(json_dict)
        assert restored.clans == pfam_data.clans

    def test_roundtrip_preserves_entry_fields(self, pfam_data):
        json_dict = _pfam_data_to_json(pfam_data)
        restored = _pfam_data_from_json(json_dict)
        for acc in pfam_data.by_accession:
            orig = pfam_data.by_accession[acc]
            rest = restored.by_accession[acc]
            assert rest.accession == orig.accession
            assert rest.shortname == orig.shortname
            assert rest.description == orig.description
            assert rest.clan_accession == orig.clan_accession
            assert rest.clan_name == orig.clan_name

    def test_roundtrip_via_json_string(self, pfam_data):
        """Serialize to JSON string and back."""
        json_str = json.dumps(_pfam_data_to_json(pfam_data))
        restored = _pfam_data_from_json(json.loads(json_str))
        assert len(restored.by_accession) == len(pfam_data.by_accession)
        assert restored.by_shortname == pfam_data.by_shortname
        assert restored.clans == pfam_data.clans

    def test_entries_are_pfam_entry_after_roundtrip(self, pfam_data):
        json_dict = _pfam_data_to_json(pfam_data)
        restored = _pfam_data_from_json(json_dict)
        for entry in restored.by_accession.values():
            assert isinstance(entry, PfamEntry)


# ===========================================================================
# Pfam 39 reference: killed families, token normalization, release
# ===========================================================================

MINI_DEAD = """\
# STOCKHOLM 1.0
#=GF ID   MR_MLE
#=GF AC   PF01188
#=GF KL   This family has been killed
#=GF FW   PF13378
#=GF CC   now covered by PF13378.
//
# STOCKHOLM 1.0
#=GF ID   Spurious
#=GF AC   PF13659
#=GF KL   This family has been killed
#=GF FW   
#=GF CC   Not what it seems; adds nothing real
//
# STOCKHOLM 1.0
#=GF ID   Old_split
#=GF AC   PF00065
#=GF KL   This family has been killed
#=GF FW   PF02931 PF02932
#=GF CC   This family has been killed and split into two
//
"""


@pytest.fixture
def dead_gz(tmp_path):
    p = tmp_path / "Pfam-A.dead.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(MINI_DEAD)
    return p


def _pfam_with_dead(pfam_data, dead):
    pfam_data.dead = dead
    for acc, sn in [("PF13378", "MR_MLE_C"), ("PF02931", "Neur_chan_LBD"), ("PF02932", "Neur_chan_memb"),
                    ("PF07969", "Amidohydro_3")]:
        pfam_data.by_accession[acc] = PfamEntry(acc, sn, sn, "", "")
        pfam_data.by_shortname[sn] = acc
    pfam_data.by_shortname["PFK"] = "PF00365"
    pfam_data.by_accession["PF00365"] = PfamEntry("PF00365", "PFK", "Phosphofructokinase", "", "")
    return pfam_data


class TestDeadFamilies:
    def test_parse(self, dead_gz):
        dead = parse_dead_families(dead_gz)
        assert set(dead) == {"PF01188", "PF13659", "PF00065"}
        assert dead["PF01188"] == {"shortname": "MR_MLE", "forward_to": ["PF13378"],
                                   "comment": "now covered by PF13378."}
        assert dead["PF13659"]["forward_to"] == []
        assert dead["PF00065"]["forward_to"] == ["PF02931", "PF02932"]

    def test_resolve_forwards_killed_family(self, pfam_data, dead_gz):
        pd_ = _pfam_with_dead(pfam_data, parse_dead_families(dead_gz))
        assert pd_.resolve_accession("PF01188") == ["PF13378"]
        assert pd_.resolve_accession("PF00065") == ["PF02931", "PF02932"]
        assert pd_.resolve_accession("PF13659") == []          # killed, no successor
        assert pd_.resolve_accession("PF00712") == ["PF00712"]  # live
        assert pd_.resolve_accession("PF55555") == ["PF55555"]  # unknown: kept

    def test_resolve_follows_chains_and_survives_cycles(self, pfam_data):
        pfam_data.dead = {"PF00010": {"forward_to": ["PF00011"]},
                          "PF00011": {"forward_to": ["PF00712"]},
                          "PF00020": {"forward_to": ["PF00021"]},
                          "PF00021": {"forward_to": ["PF00020"]}}
        assert pfam_data.resolve_accession("PF00010") == ["PF00712"]
        assert pfam_data.resolve_accession("PF00020") == []


class TestNormalizePfamToken:
    @pytest.fixture
    def pd_(self, pfam_data, dead_gz):
        return _pfam_with_dead(pfam_data, parse_dead_families(dead_gz))

    @pytest.mark.parametrize("token, expected", [
        ("PF00712", "PF00712"),
        ("PF00712.24", "PF00712"),      # versioned accession
        ("PF7969", "PF07969"),          # Cyanorak typo, padded accession exists
        ("PF1234", None),               # padded accession unknown: not guessed
        ("PFK", "PF00365"),             # a short name that starts with PF
        ("PFOR_II", None),              # unknown short name
        ("DNA_pol3_beta", "PF00712"),
        ("TIGR00663", None),
        ("PF01188", "PF01188"),         # killed family: normalized, forwarded later
    ])
    def test_tokens(self, pd_, token, expected):
        assert normalize_pfam_token(token, pd_) == expected


class TestReleaseAndRoundtrip:
    def test_parse_release(self, tmp_path):
        p = tmp_path / "Pfam.version.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("Pfam release       : 39.0\nPfam-A families    : 53275\n")
        assert parse_release(p) == "39.0"

    def test_dead_and_release_roundtrip(self, pfam_data, dead_gz):
        pfam_data.dead = parse_dead_families(dead_gz)
        pfam_data.release = "39.0"
        back = _pfam_data_from_json(json.loads(json.dumps(_pfam_data_to_json(pfam_data))))
        assert back.dead == pfam_data.dead and back.release == "39.0"

    def test_old_cache_without_dead_loads(self, pfam_data):
        data = _pfam_data_to_json(pfam_data)
        data.pop("dead"); data.pop("release")
        back = _pfam_data_from_json(data)
        assert back.dead == {} and back.release == ""
