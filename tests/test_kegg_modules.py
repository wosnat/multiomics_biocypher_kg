"""Unit tests for the KEGG module definition evaluator (multiomics_kg.utils.kegg_modules)."""
import pytest

from multiomics_kg.utils.kegg_modules import evaluate_organism, ko_roles, parse_steps

NITRATE = "(K00367,K10534,K00372+K00360) (K00366,K17877)"           # M00531
GLYC = ("K01803 (((K00134,K00150) K00927),K11389) (K01834,K15633,K15634,K15635) "
        "(K01689,K27394) (K00873,K12406)")                            # M00002


def ev(defn, kos, ec=None, extra=None):
    mods = {"M1": parse_steps(defn), **(extra or {})}
    return evaluate_organism(mods, set(kos), ec or {})["M1"]


def test_steps_and_dash_placeholder():
    assert [s.definition for s in parse_steps("K00001 -- K00002")] == ["K00001", "K00002"]
    assert [s.index for s in parse_steps("K00001 -- K00002")] == [1, 2]
    assert len(parse_steps("K00001 K00002\nK00003")) == 3                     # multi-line joined


def test_step_definition_keeps_group_text():
    assert [s.definition for s in parse_steps(NITRATE)] == [
        "(K00367,K10534,K00372+K00360)", "(K00366,K17877)"]
    assert len(parse_steps(GLYC)) == 5


def test_or_and_complex_optional():
    assert ev(NITRATE, {"K00366", "K00372"}).steps_present == 1               # complex needs K00360 too
    assert ev(NITRATE, {"K00366", "K00372", "K00360"}).status == "complete"
    assert ev("K00001+K00002-K00003", {"K00001", "K00002"}).status == "complete"   # -K optional
    assert ev("K00001+(K00002,K00003)", {"K00001", "K00003"}).status == "complete"
    assert ev("(K00001,K00002)-K00003", {"K00002"}).status == "complete"
    assert ev("-K00003+K00001", {"K00001"}).status == "complete"               # leading optional
    assert ev("K00001+K00002", {"K00001"}).status == "one_step_missing"


def test_and_inside_group():
    assert ev(GLYC, {"K01803", "K00134", "K01834", "K01689", "K00873"}).missing_step_indices == [2]
    assert ev(GLYC, {"K01803", "K11389", "K01834", "K01689", "K00873"}).status == "complete"


def test_missing_fields_and_status():
    r = ev(NITRATE, {"K00366"})
    assert (r.status, r.missing_step_indices) == ("one_step_missing", [1])
    assert r.missing_step_definitions == ["(K00367,K10534,K00372+K00360)"]
    assert r.missing_kos == ["K00360", "K00367", "K00372", "K10534"]
    assert ev("K00001 K00002 K00003", {"K00001"}).status == "two_steps_missing"
    assert ev("K00001 K00002 K00003 K00004", {"K00001"}).status == "partial"


def test_ec_fallback_fills_plain_term():
    r = ev(GLYC, {"K01803", "K00134", "K00927", "K01834", "K01689"}, {"K00873": {"PMN2A_0295"}})
    assert (r.status, r.steps_present_via_ec, r.ec_filled_step_indices, r.ec_genes) == (
        "complete", 1, [5], ["PMN2A_0295"])


def test_ec_never_fills_complex_subunit():
    r = ev(NITRATE, {"K00366"}, {"K00372": {"g1"}, "K00360": {"g1"}})
    assert r.steps_present == 1 and r.ec_genes == []


def test_ec_one_step_per_gene_and_ko_wins():
    r = ev("K00001 K00002", set(), {"K00001": {"g1"}, "K00002": {"g1"}})
    assert r.steps_present == 1
    r = ev("K00001 K00002", {"K00001"}, {"K00001": {"g1"}, "K00002": {"g1"}})
    assert (r.steps_present, r.steps_present_via_ec, r.ec_filled_step_indices) == (2, 1, [2])


def test_ec_matching_maximises_filled_steps():
    # g1 can fill step 1 or 2, g2 only step 1: a maximum matching fills both.
    r = ev("K00001 K00002", set(), {"K00001": {"g1", "g2"}, "K00002": {"g1"}})
    assert (r.steps_present_via_ec, r.ec_genes) == (2, ["g2", "g1"])


def test_module_reference_via_ec_propagates():
    sub = {"M00902": parse_steps("K00001")}
    r = ev("M00902 K00002", {"K00002"}, {"K00001": {"g1"}}, extra=sub)
    assert (r.status, r.steps_present_via_ec, r.ec_genes) == ("complete", 1, ["M00902"])
    assert ev("M00902 K00002", {"K00002"}, {}, extra=sub).missing_kos == []
    assert ev("M00902 K00002", {"K00002", "K00001"}, {}, extra=sub).steps_present_via_ec == 0


def test_missing_kos_skip_optional_and_ko_pass_ignores_ec_completed_reference():
    assert ev("K00001+K00002-K00003", set()).missing_kos == ["K00001", "K00002"]
    sub = {"M00902": parse_steps("K00001")}
    r = ev("(M00902,K00009) K00002", {"K00002"}, {"K00001": {"g1"}}, extra=sub)
    assert (r.ec_filled_step_indices, r.ec_genes) == ([1], ["M00902"])   # KO-only pass sees step 1 missing


def test_reference_errors():
    with pytest.raises(ValueError):
        evaluate_organism({"M1": parse_steps("M00909 K00001")}, set(), {})
    with pytest.raises(ValueError):
        evaluate_organism({"M00901": parse_steps("M00902"), "M00902": parse_steps("M00901")}, set(), {})


def test_ko_roles_and_conflict():
    assert ko_roles(parse_steps(NITRATE))["K00360"] == ([1], "required")
    assert ko_roles(parse_steps("K00001+K00002-K00003"))["K00003"] == ([1], "optional")
    assert ko_roles(parse_steps("K00001 (K00001,K00002)"))["K00001"] == ([1, 2], "required")
    with pytest.raises(ValueError):
        ko_roles(parse_steps("K00001-K00003 K00003"))


def test_standalone_optional_step_is_dropped():
    """A whole step written as an optional subunit (M00955 `... K21749 -K08080`) never needs a
    gene: it is dropped like `--`, so the module can be complete without it."""
    assert [s.definition for s in parse_steps("K00001 -K00002")] == ["K00001"]
    assert ev("K00001 -K00002", {"K00001"}).status == "complete"
    assert "K00002" not in ko_roles(parse_steps("K00001 -K00002"))


def test_ec_gene_must_close_the_step_on_its_own():
    # Step 1 needs K00001 AND K00002; the organism has K00001. gX stands in for K00001 (already
    # present) and cannot close step 1; gB stands in for K00002 and K00003 — one gene, one step.
    r = ev("(K00001 K00002) K00003", {"K00001"}, {"K00001": {"gX"}, "K00002": {"gB"}, "K00003": {"gB"}})
    assert (r.steps_present, r.steps_present_via_ec) == (1, 1)
    assert len(r.ec_genes) == 1 and r.ec_genes == ["gB"]


def test_step_needing_two_ec_genes_stays_missing():
    r = ev("(K00001 K00002)", set(), {"K00001": {"g1"}, "K00002": {"g2"}})
    assert (r.steps_present, r.ec_genes) == (0, [])


def test_ec_gene_named_is_the_one_that_supplies_the_missing_ko():
    # MIT9313 M00097 shape: only K15744 is missing; g_iso supplies it, g_other maps to a KO the
    # organism already has.
    r = ev("((K02293 K15744 K00514 K09835),K10027,K15745)", {"K02293", "K00514", "K09835"},
           {"K15744": {"g_iso"}, "K09835": {"a_other"}})
    assert r.ec_genes == ["g_iso"]
