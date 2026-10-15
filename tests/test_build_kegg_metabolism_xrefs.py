"""Unit tests for step 6: unified pruned kegg_data.json builder."""
from __future__ import annotations
import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def _make_resolver(tmp_path):
    """Minimal MNX resolver with synthetic compounds and reactions."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
            reference TEXT, classifs TEXT, is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'I=', 'WQZ-...', 'OC[CH]1...', 0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO reactions VALUES
            ('MNXR101234', '...', 'kegg.reaction:R00200', '2.7.1.40', 'B', NULL);
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
    """)
    conn.commit()
    return conn


@pytest.fixture
def synthetic_kegg_raw(tmp_path):
    """Return a kegg_data dict matching the parser output, suitable as input to step 6."""
    return {
        "ko_names":               {"K02338": "DNA polymerase III subunit beta",
                                   "K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00710": "Carbon fixation",
                                   "ko03030": "DNA replication",
                                   "ko00500": "Starch and sucrose metabolism"},
        "ko_to_pathways":         {"K02338": ["ko03030"], "K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00710": "09102",
                                   "ko03030": "09124", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism",
                                   "09102": "Energy metabolism",
                                   "09124": "Replication and repair"},
        "subcategory_to_category":{"09101": "09100", "09102": "09100", "09124": "09120"},
        "category_names":         {"09100": "Metabolism", "09120": "Genetic Information Processing"},
        "reaction_names":         {"R00200": "ATP:pyruvate ..."},
        "compound_names":         {"C00031": "D-glucose", "C00074": "PEP",
                                   "C99999": "compound-only-pathway-trigger"},
        "reaction_to_compounds":  {"R00200": ["C00031", "C00074"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"],
                                   "C99999": ["ko99999"]},   # ko99999 is compound-only
    }


def test_build_pruned_kegg_data(tmp_path, monkeypatch, synthetic_kegg_raw):
    """End-to-end: 1 strain, 2 genes, ensure pruned structure is correct."""
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
        "gene2": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])

    conn = _make_resolver(tmp_path)
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    data = bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, sets=sets)

    # KOs: only gene-annotated KOs survive
    assert set(data["kos"].keys()) == {"K02338", "K01006"}
    assert data["kos"]["K02338"]["name"] == "DNA polymerase III subunit beta"
    assert data["kos"]["K02338"]["pathways"] == ["ko03030"]

    # Pathways: KO ∪ Rxn-reachable. K02338→ko03030, K01006→ko00010, R00200→ko00010
    # ko00710 not present — no KO/Rxn maps there.
    assert set(data["pathways"].keys()) == {"ko03030", "ko00010"}

    # Reactions: only gene-reachable
    assert set(data["reactions"].keys()) == {"R00200"}
    rxn = data["reactions"]["R00200"]
    assert rxn["name"].startswith("ATP:pyruvate")
    assert rxn["pathways"] == ["ko00010"]
    assert rxn["compounds"] == ["C00031", "C00074"]
    assert rxn["mnxr_id"] == "MNXR101234"
    assert rxn["ec_numbers"] == ["2.7.1.40"]
    assert rxn["mass_balance"] == "balanced"
    assert rxn["reaction_class"] == "chemical"

    # Compounds: only reaction-reachable
    assert set(data["compounds"].keys()) == {"C00031", "C00074"}
    glucose = data["compounds"]["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"

    # Option B: compound.pathways pruned to KO ∪ Rxn-reachable set.
    # ko00500 is NOT in the set (no KO or Rxn maps there in the fixture).
    # ko00010 IS in the set.
    assert glucose["pathways"] == ["ko00010"]

    # Subcategories/categories present iff parent of an emitted pathway
    assert set(data["subcategories"].keys()) == {"09101", "09124"}
    assert data["subcategories"]["09101"]["category"] == "09100"
    assert set(data["categories"].keys()) == {"09100", "09120"}


def test_compound_only_pathways_dropped(tmp_path, monkeypatch, synthetic_kegg_raw):
    """The C99999→ko99999 compound-only pathway should not appear anywhere.

    Critical for Option B: even if a compound were in the pruned set, its pathways
    should be filtered to KO ∪ Rxn-reachable. C99999 isn't even in the pruned set
    (no reaction reaches it), so this is a double-defensive check.
    """
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    data = bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, sets=sets)
    assert "C99999" not in data["compounds"]
    assert "ko99999" not in data["pathways"]


def test_compound_pathway_filter_drops_unevidenced(tmp_path, monkeypatch):
    """Compound pathways are filtered to KO ∪ Rxn-reachable.

    Force a case where a compound IS in the gene-reachable set (via R00200) but has
    a pathway (ko00500) that's not reachable from any gene-KO or gene-reaction.
    Expected: compound emitted, but ko00500 NOT in compound.pathways.
    """
    raw = {
        "ko_names":               {"K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00500": "Starch and sucrose"},
        "ko_to_pathways":         {"K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism"},
        "subcategory_to_category":{"09101": "09100"},
        "category_names":         {"09100": "Metabolism"},
        "reaction_names":         {"R00200": "rxn"},
        "compound_names":         {"C00031": "D-glucose"},
        "reaction_to_compounds":  {"R00200": ["C00031"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"]},  # ko00500 unevidenced
    }
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "g1": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    sets = bx._gene_reachable_sets(raw)
    data = bx.build_pruned_kegg_data(raw, conn, sets=sets)

    # ko00500 must NOT be in pathways top-level (Option B: only KO ∪ Rxn-reachable)
    assert "ko00500" not in data["pathways"]
    # And compound's pathway list filtered
    assert data["compounds"]["C00031"]["pathways"] == ["ko00010"]


def test_bulk_enrich_reactions_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_reactions(conn, kegg_ids, allowed_pathways, raw) returns
    a dict mapping kegg_reaction_id → enrichment dict, with the same fields
    as the per-entity _enrich_reaction.
    """
    conn = _make_resolver(tmp_path)
    raw = {
        "reaction_names": {"R00200": "ATP:pyruvate ..."},
        "reaction_to_pathways": {"R00200": ["ko00010", "ko00710"]},
        "reaction_to_compounds": {"R00200": ["C00031"]},
    }
    allowed = {"ko00010", "ko00710"}

    result = bx._bulk_enrich_reactions(conn, ["R00200", "R99999"], allowed, raw)

    # Every requested ID is in the result (R99999 has no MNX entry but still gets a stub)
    assert set(result.keys()) == {"R00200", "R99999"}

    r = result["R00200"]
    assert r["name"].startswith("ATP:pyruvate")
    assert r["mnxr_id"] == "MNXR101234"
    assert r["ec_numbers"] == ["2.7.1.40"]
    assert r["mass_balance"] == "balanced"
    assert r["reaction_class"] == "chemical"
    assert r["pathways"] == ["ko00010", "ko00710"]
    assert r["compounds"] == ["C00031"]

    r99 = result["R99999"]
    assert r99["mnxr_id"] is None
    assert r99["ec_numbers"] == []
    assert r99["mass_balance"] == "unbalanced"  # default
    assert r99["reaction_class"] == "chemical"  # default


def test_step6_builds_tcdb_hierarchy(tmp_path):
    """Step 6 must produce cache/data/tcdb/tcdb_hierarchy.json (was previously
    built by build_metabolite_resolver.py)."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    # Stage minimal TCDB TSVs the way download_metabolism_reference.download_all would:
    cache_root = tmp_path / "cache" / "data"
    tcdb_dir = cache_root / "tcdb"
    raw_dir = tcdb_dir / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "families.tsv").write_text("1.A.1\tThe Voltage-gated Ion Channel\n")
    (raw_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (raw_dir / "substrates.tsv").write_text(
        "1.A.1.5.2\tCHEBI:9314;sucrose|CHEBI:3308;calcium(2+)\n"
    )
    (raw_dir / "acc2tcid.tsv").write_text("P12345\t1.A.1.5.2\n")

    out = tcdb_dir / "tcdb_hierarchy.json"
    mod._build_tcdb_hierarchy(cache_root=cache_root)

    assert out.exists()
    import json
    h = json.loads(out.read_text())
    # 5 levels for 1.A.1.5.2 → expect class+subclass+family+subfamily+specificity
    assert "1" in h
    assert "1.A" in h
    assert "1.A.1" in h
    assert "1.A.1.5" in h
    assert "1.A.1.5.2" in h
    assert h["1.A.1.5.2"]["level_kind"] == "tc_specificity"
    # substrate_classes preserves the full 'CHEBI:N;name' string so downstream
    # MNX resolution can map to a canonical Metabolite primary ID.
    assert h["1.A.1.5.2"]["substrate_classes"] == [
        "CHEBI:9314;sucrose",
        "CHEBI:3308;calcium(2+)",
    ]


def test_build_tcdb_hierarchy_seeds_3part_families_without_5part_children(tmp_path):
    """Family-level entries in families.tsv must land in the hierarchy even
    when superfamilies.tsv has no 5-part child for them. Regression for the
    1.A.11 (Ammonium Transporter Channel) class of dropped seeds."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    raw_dir = cache_root / "tcdb" / "raw"
    raw_dir.mkdir(parents=True)
    # 1.A.11 has a description but no curated 5-part children
    (raw_dir / "families.tsv").write_text(
        "1.A.1\tThe Voltage-gated Ion Channel\n"
        "1.A.11\tThe Ammonium Transporter Channel (Amt) Family\n"
    )
    (raw_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (raw_dir / "substrates.tsv").write_text("")
    (raw_dir / "acc2tcid.tsv").write_text("")

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    assert "1.A.11" in h
    assert h["1.A.11"]["level_kind"] == "tc_family"
    assert h["1.A.11"]["name"] == "The Ammonium Transporter Channel (Amt) Family"
    assert h["1.A.11"]["parent"] == "1.A"


def test_build_tcdb_hierarchy_seeds_5part_from_acc2tcid(tmp_path):
    """5-part TCIDs present only in acc2tcid.tsv (not in superfamilies.tsv)
    must land in the hierarchy with synthesized ancestors. Regression for
    1.A.12.2.2 / 1.A.30.2.1 / 2.A.23.1.1 class of dropped seeds."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    raw_dir = cache_root / "tcdb" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "families.tsv").write_text("1.A.12\tThe Intracellular Chloride Channel\n")
    (raw_dir / "superfamilies.tsv").write_text("")
    # The TCID has a substrate even though no superfamilies row exists
    (raw_dir / "substrates.tsv").write_text("1.A.12.2.2\tCHEBI:3731;chloride\n")
    (raw_dir / "acc2tcid.tsv").write_text("Q9X9X9\t1.A.12.2.2\n")

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    # Full 5-level chain synthesized
    assert {"1", "1.A", "1.A.12", "1.A.12.2", "1.A.12.2.2"} <= set(h.keys())
    assert h["1.A.12.2.2"]["level_kind"] == "tc_specificity"
    assert h["1.A.12.2.2"]["substrate_classes"] == ["CHEBI:3731;chloride"]
    assert h["1.A.12.2"]["level_kind"] == "tc_subfamily"
    assert h["1.A.12"]["name"] == "The Intracellular Chloride Channel"


def _write_tcdb_tsv_fixture(cache_root):
    """Minimal 4-TSV fixture shared by the names-artifact tests (same shape as
    test_build_tcdb_hierarchy_seeds_5part_from_acc2tcid)."""
    raw_dir = cache_root / "tcdb" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / "families.tsv").write_text("1.A.1\tThe Voltage-gated Ion Channel\n")
    (raw_dir / "superfamilies.tsv").write_text(
        "1.A.1.5.2\t1.A.1.5\t1.A.1\tVIC\tVIC Superfamily\n"
    )
    (raw_dir / "substrates.tsv").write_text("")
    (raw_dir / "acc2tcid.tsv").write_text("P12345\t1.A.1.5.2\n")


def _write_names_file(cache_root, names, scraped_families):
    tcdb_dir = cache_root / "tcdb"
    tcdb_dir.mkdir(parents=True, exist_ok=True)
    (tcdb_dir / "tcdb_names.json").write_text(json.dumps({
        "meta": {"scraped_at": "2026-08-18", "source": "test",
                 "scraped_families": scraped_families, "failed_families": []},
        "names": names,
    }))


def test_build_tcdb_hierarchy_consumes_names_artifact(tmp_path):
    """Names from the committed tcdb_names.json (T6 scraper) fill every level;
    tc_specificity entries also gain a sparse description."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    _write_tcdb_tsv_fixture(cache_root)
    _write_names_file(cache_root, {
        "1": {"name": "Channels/Pores"},
        "1.A": {"name": "\u03b1-Type Channels"},
        "1.A.1.5": {"name": "The Epsilon Subfamily"},
        "1.A.1.5.2": {"name": "Short name.",
                      "description": "Short name. Much longer curated text."},
    }, scraped_families=["1.A.1"])

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    assert h["1"]["name"] == "Channels/Pores"     # replaces _TC_CLASS_NAMES
    assert h["1.A"]["name"] == "\u03b1-Type Channels"  # was always ""
    assert h["1.A.1.5"]["name"] == "The Epsilon Subfamily"
    assert h["1.A.1.5.2"]["name"] == "Short name."
    assert h["1.A.1.5.2"]["description"] == "Short name. Much longer curated text."
    # families.tsv stays the primary family-name source
    assert h["1.A.1"]["name"] == "The Voltage-gated Ion Channel"


def test_build_tcdb_hierarchy_without_names_file_leaves_blank(tmp_path):
    """No names artifact → deep levels stay unnamed; class names come ONLY
    from the artifact now (no hardcoded fallback)."""
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    _write_tcdb_tsv_fixture(cache_root)

    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())

    assert h["1"]["name"] == ""
    assert h["1.A"]["name"] == ""
    assert "description" not in h["1.A.1.5.2"]


def test_warn_unscraped_families(tmp_path, caplog):
    """Spec S4 guard: kept 4/5-part ids whose family the scraper has not
    covered are reported loudly with a top-up command."""
    import logging as _logging
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    cache_root = tmp_path / "cache" / "data"
    _write_names_file(cache_root, {}, scraped_families=["2.A.1"])
    names_path = cache_root / "tcdb" / "tcdb_names.json"

    with caplog.at_level(_logging.WARNING):
        missing = mod._warn_unscraped_families(
            {"2.A.1.1", "3.A.1.5.9", "1.A"},  # 1.A too shallow to need scraping
            names_path)
    assert missing == ["3.A.1"]
    assert "3.A.1" in caplog.text and "scrape_tcdb_names" in caplog.text


def test_warn_unscraped_families_missing_file(tmp_path, caplog):
    import logging as _logging
    from multiomics_kg.download import build_kegg_metabolism_xrefs as mod

    with caplog.at_level(_logging.WARNING):
        missing = mod._warn_unscraped_families(
            {"2.A.1.1"}, tmp_path / "nope" / "tcdb_names.json")
    assert missing == ["2.A.1"]
    assert "MISSING" in caplog.text


def test_bulk_enrich_compounds_returns_dict_keyed_by_kegg_id(tmp_path):
    """_bulk_enrich_compounds returns a dict mapping kegg_compound_id → enrichment."""
    conn = _make_resolver(tmp_path)
    raw = {
        "compound_names": {"C00031": "D-glucose", "C99999": "obscure"},
        "compound_to_pathways": {
            "C00031": ["ko00010", "ko00500"],
            "C99999": [],
        },
    }
    allowed = {"ko00010"}

    result = bx._bulk_enrich_compounds(conn, ["C00031", "C99999"], allowed, raw)

    assert set(result.keys()) == {"C00031", "C99999"}

    glucose = result["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["pathways"] == ["ko00010"]  # ko00500 filtered out by allowed

    obscure = result["C99999"]
    assert obscure["mnxm_id"] is None
    assert obscure["chebi_id"] is None
    assert obscure["formula"] is None


# ── Step 6 part 2: TCDB substrate resolution + evidence_sources tagging ──────


def test_compounds_get_metabolism_evidence_source():
    """A compound in catalysis_cpds with no transport overlap gets
    evidence_sources=['metabolism']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {
        "compounds": {
            "C00031": {"name": "D-Glucose", "formula": "C6H12O6"},
            "C00002": {"name": "ATP", "formula": "C10H16N5O13P3"},
        },
        "reactions": {},
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds={"C00031", "C00002"},
        substrate_kegg_cpds=set(),
        substrate_primaries=set(),
        compound_props={},
    )
    for cpd in kegg_data["compounds"].values():
        assert cpd["evidence_sources"] == ["metabolism"]


def test_transport_only_compounds_land_in_additional_compounds():
    """A non-KEGG TCDB substrate primary becomes an additional_compounds entry
    with evidence_sources=['transport']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {"compounds": {}, "reactions": {}}
    substrate_primaries = {"chebi:9999"}
    compound_props = {
        "chebi:9999": {
            "name": "tetracycline", "formula": "C22H24N2O8",
            "mass": 444.43, "inchikey": None,
            "mnxm_id": "MNXM00099", "chebi_id": "9999",
        }
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds=set(),
        substrate_kegg_cpds=set(),
        substrate_primaries=substrate_primaries,
        compound_props=compound_props,
    )
    assert "chebi:9999" in kegg_data["additional_compounds"]
    entry = kegg_data["additional_compounds"]["chebi:9999"]
    assert entry["evidence_sources"] == ["transport"]
    assert entry["chebi_id"] == "9999"


def test_overlap_compound_gets_both_evidence_sources():
    """A compound in BOTH catalysis_cpds AND substrate_kegg_cpds gets
    evidence_sources=['metabolism', 'transport']."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _fold_substrates_into_kegg_data,
    )
    kegg_data = {
        "compounds": {"C00089": {"name": "Sucrose", "formula": "C12H22O11"}},
        "reactions": {},
    }
    _fold_substrates_into_kegg_data(
        kegg_data,
        catalysis_cpds={"C00089"},
        substrate_kegg_cpds={"C00089"},
        substrate_primaries={"kegg.compound:C00089"},
        compound_props={},
    )
    # KEGG-flavor substrate stays out of additional_compounds — it's already in compounds.
    assert "kegg.compound:C00089" not in kegg_data["additional_compounds"]
    assert "C00089" not in kegg_data["additional_compounds"]
    assert kegg_data["compounds"]["C00089"]["evidence_sources"] == ["metabolism", "transport"]


def test_transport_only_kegg_compound_in_compounds_with_pathways(tmp_path, monkeypatch):
    """Refactor invariant: when a kegg.compound primary is in substrate_kegg_cpds
    but NOT in catalysis_cpds, build_pruned_kegg_data (called with the extended
    cpds set) puts the entry in compounds (not additional_compounds) AND populates
    its pathways field with the gene-reachable pathways (extended pws set).

    Setup: the substrate compound C00089 participates in 2 KEGG pathways. Only
    one (ko00010) is in the extended pws; ko00500 is excluded. Verify the entry
    lands in compounds with pathways=['ko00010'].
    """
    raw = {
        "ko_names":               {"K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00500": "Starch and sucrose"},
        "ko_to_pathways":         {"K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism"},
        "subcategory_to_category": {"09101": "09100"},
        "category_names":         {"09100": "Metabolism"},
        "reaction_names":         {},
        "compound_names":         {"C00089": "Sucrose"},
        "reaction_to_compounds":  {},
        "reaction_to_pathways":   {},
        "compound_to_pathways":   {"C00089": ["ko00010", "ko00500"]},
    }
    # Catalysis side: 1 KO (ko00010), no reactions, no compounds.
    catalysis_sets = {
        "kos": {"K01006"},
        "rxns": set(),
        "cpds": set(),
        "pws": {"ko00010"},
        "tcdb_ids": set(),
    }
    # Caller-extended: C00089 is a transport-reachable kegg compound. Curate
    # extended_pws = {"ko00010"} so ko00500 stays excluded; we want to confirm
    # the compound's pathways field is filtered to that extended set.
    extended_sets = {
        **catalysis_sets,
        "cpds": {"C00089"},
        "pws": {"ko00010"},
    }

    conn = sqlite3.connect(":memory:")
    conn.executescript("""
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, formula TEXT, mass REAL,
            inchikey TEXT, smiles TEXT);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, classifs TEXT,
            is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
    """)
    data = bx.build_pruned_kegg_data(raw, conn, sets=extended_sets)
    # Entry lands in compounds (not additional_compounds — the latter is created
    # by _fold_substrates_into_kegg_data, not build_pruned_kegg_data).
    assert "C00089" in data["compounds"]
    # Pathways filtered to extended pws — ko00500 dropped.
    assert data["compounds"]["C00089"]["pathways"] == ["ko00010"]
    # ko00500 does not appear in top-level pathways either.
    assert "ko00500" not in data["pathways"]


def test_prune_tcdb_walks_up_only():
    """_prune_tcdb walks UP to tc_class and does NOT walk down to descendants.

    The downward arm used to expand 528 gene-annotated seeds into 12,902 nodes,
    94.5% of which had gene_count = 0. It existed so leaf substrate data stayed
    reachable, but `_compute_subtree_substrate_strings` already rolls substrates
    up over the FULL hierarchy pre-pruning, so every kept ancestor carries its
    descendants' substrates without the descendant nodes existing.
    """
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"name": "Channels", "level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"name": "VIC", "level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"name": "", "level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["calcium(2+)"]},
        "1.A.1.5.3": {"name": "", "level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["sodium(1+)"]},
        # An unrelated branch that should NOT be kept
        "2": {"name": "ECP", "level": 0, "level_kind": "tc_class", "parent": None},
        "2.A": {"name": "", "level": 1, "level_kind": "tc_subclass", "parent": "2"},
    }
    # Seed at family level — walk up to 1 only; the two leaves below stay out
    kept, seed_aliases = _prune_tcdb(hierarchy, {"1.A.1"})
    assert kept == {"1", "1.A", "1.A.1"}
    assert seed_aliases == {}
    # No downward expansion
    assert "1.A.1.5" not in kept
    assert "1.A.1.5.2" not in kept
    assert "1.A.1.5.3" not in kept
    # Branch 2 untouched
    assert "2" not in kept


def test_prune_tcdb_keeps_only_seeds_and_ancestors():
    """The kept set is exactly ⋃(seed ∪ ancestors(seed)) — nothing else.

    Guards the invariant that made the cleanup safe: on real data the new kept
    set was a strict subset of the old one (704 ⊂ 12,902) and lost zero
    metabolite reachability.
    """
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"level": 4, "level_kind": "tc_specificity", "parent": "1.A.1.5"},
        "1.A.2": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.2.9": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.2"},
    }
    kept, _ = _prune_tcdb(hierarchy, {"1.A.1.5.2", "1.A.2"})
    assert kept == {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2", "1.A.2"}
    # Sibling subfamily under a seeded family is NOT pulled in
    assert "1.A.2.9" not in kept


def test_prune_tcdb_seed_at_leaf():
    """Seeding at a leaf still walks up to root.

    Regression: previously, sibling specificities under the same subfamily were
    NOT in `kept`. Substrate rollup is now decoupled (computed pre-prune in
    `_compute_subtree_substrate_strings`), so this test only asserts the
    structural result — siblings still aren't kept (correct: no gene reaches
    them, so no TcdbFamily node is needed for them).
    """
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5", "substrate_classes": ["x"]},
    }
    kept, seed_aliases = _prune_tcdb(hierarchy, {"1.A.1.5.2"})
    assert kept == {"1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2"}
    assert seed_aliases == {}


def test_prune_tcdb_unknown_seed_walks_up_to_ancestor():
    """Seeds not in the hierarchy are remapped to the nearest curated ancestor
    (so retired eggNOG TCIDs collapse to family-level, not lost)."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "3": {"level": 0, "level_kind": "tc_class", "parent": None},
        "3.A": {"level": 1, "level_kind": "tc_subclass", "parent": "3"},
        "3.A.1": {"level": 2, "level_kind": "tc_family", "parent": "3.A"},
    }
    # 3.A.1.35 doesn't exist — should anchor onto 3.A.1.
    kept, seed_aliases = _prune_tcdb(hierarchy, {"3.A.1.35"})
    assert seed_aliases == {"3.A.1.35": "3.A.1"}
    assert kept == {"3", "3.A", "3.A.1"}


def test_prune_tcdb_unknown_seed_with_no_ancestor_marked_dropped():
    """When no truncation hits an ancestor, the seed is recorded with empty
    alias so the caller can warn and drop it."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _prune_tcdb
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
    }
    kept, seed_aliases = _prune_tcdb(hierarchy, {"99.X.Y.Z.Q"})
    assert kept == set()
    assert seed_aliases == {"99.X.Y.Z.Q": ""}


def test_compute_subtree_substrate_strings_rolls_up_through_full_hierarchy():
    """Bug fix: ancestor nodes capture substrates from EVERY tc_specificity
    descendant in the FULL hierarchy, not just gene-annotated ones. Done
    pre-pruning so pruning is a pure structural filter.
    """
    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _compute_subtree_substrate_strings,
    )
    hierarchy = {
        "1": {"level": 0, "level_kind": "tc_class", "parent": None},
        "1.A": {"level": 1, "level_kind": "tc_subclass", "parent": "1"},
        "1.A.1": {"level": 2, "level_kind": "tc_family", "parent": "1.A"},
        "1.A.1.5": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.5.2": {"level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5",
                      "substrate_classes": ["CHEBI:1;a", "CHEBI:2;b"]},
        "1.A.1.5.3": {"level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.5",
                      "substrate_classes": ["CHEBI:3;c"]},
        # A second subfamily under the same family — its substrate must roll up
        # to 1.A.1 even though no seed touches its leaves.
        "1.A.1.6": {"level": 3, "level_kind": "tc_subfamily", "parent": "1.A.1"},
        "1.A.1.6.1": {"level": 4, "level_kind": "tc_specificity",
                      "parent": "1.A.1.6",
                      "substrate_classes": ["CHEBI:4;d"]},
    }
    rollup = _compute_subtree_substrate_strings(hierarchy)
    # Family node sees substrates from ALL descendant specificities
    assert rollup["1.A.1"] == {"CHEBI:1;a", "CHEBI:2;b", "CHEBI:3;c", "CHEBI:4;d"}
    # Subfamily 1.A.1.5 sees only its own children
    assert rollup["1.A.1.5"] == {"CHEBI:1;a", "CHEBI:2;b", "CHEBI:3;c"}
    # Sibling subfamily 1.A.1.6 sees only its own child
    assert rollup["1.A.1.6"] == {"CHEBI:4;d"}
    # Specificity node = its own substrate_classes
    assert rollup["1.A.1.5.2"] == {"CHEBI:1;a", "CHEBI:2;b"}
    # Class node sees union of everything below
    assert rollup["1"] == {"CHEBI:1;a", "CHEBI:2;b", "CHEBI:3;c", "CHEBI:4;d"}


def test_resolve_substrate_strings_strips_chebi_prefix(tmp_path):
    """Regression: substrate strings 'CHEBI:NNNN;name' must resolve through the
    MNX SQLite, which stores chebi aliases as (source='chebi', value='<NNNN>') —
    NOT 'CHEBI:NNNN'."""
    import sqlite3

    from multiomics_kg.download.build_kegg_metabolism_xrefs import (
        _resolve_substrate_strings,
    )

    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compound_aliases (
            source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id)
        );
        CREATE TABLE compounds (
            mnxm_id TEXT PRIMARY KEY,
            name TEXT, formula TEXT, mass REAL, inchikey TEXT
        );
    """)
    cur.execute("INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000')")
    cur.execute("INSERT INTO compounds VALUES ('MNXM50000', 'sucrose', 'C12H22O11', 342.30, 'CZMRCDWAGMRECN-UGDNZRGBSA-N')")
    conn.commit()

    string_to_primary, compound_props = _resolve_substrate_strings(
        {"CHEBI:9314;sucrose"}, conn,
    )

    assert string_to_primary == {"CHEBI:9314;sucrose": "chebi:9314"}
    entry = compound_props["chebi:9314"]
    assert entry["mnxm_id"] == "MNXM50000"
    assert entry["chebi_id"] == "9314"
    assert entry["name"] == "sucrose"
    assert entry["formula"] == "C12H22O11"


def _resolve_subs_db(tmp_path):
    """Shared resolver fixture for _resolve_substrate_strings tests."""
    import sqlite3
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compound_aliases (
            source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id)
        );
        CREATE TABLE compounds (
            mnxm_id TEXT PRIMARY KEY,
            name TEXT, formula TEXT, mass REAL, inchikey TEXT
        );
    """)
    return conn


def test_resolve_substrate_strings_dedupes_compound_props(tmp_path):
    """Two strings citing the same CHEBI value yield one compound_props entry;
    sorted iteration → first-by-string-sort name_part wins."""
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.execute("INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000')")
    # No `compounds` row → mnx name is None, so fallback = name_part from input
    conn.commit()

    string_to_primary, compound_props = bx._resolve_substrate_strings(
        {"CHEBI:9314;sucrose", "CHEBI:9314;table sugar"},
        conn,
    )

    assert string_to_primary == {
        "CHEBI:9314;sucrose": "chebi:9314",
        "CHEBI:9314;table sugar": "chebi:9314",
    }
    assert list(compound_props.keys()) == ["chebi:9314"]
    # Sorted iteration: 'CHEBI:9314;sucrose' < 'CHEBI:9314;table sugar' → 'sucrose' wins
    assert compound_props["chebi:9314"]["name"] == "sucrose"


def test_resolve_substrate_strings_kegg_primary_skips_compound_props(tmp_path):
    """When CHEBI resolves to an mnxm with a kegg.compound alias, the primary is
    `kegg.compound:Cnnnnn` and compound_props has NO entry — KEGG primaries flow
    through `_bulk_enrich_compounds` via the extended cpds set."""
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.executescript("""
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', 'C6H12O6', 180.06, 'WQZ-');
    """)
    conn.commit()

    string_to_primary, compound_props = bx._resolve_substrate_strings(
        {"CHEBI:17234;glucose"},
        conn,
    )

    assert string_to_primary == {"CHEBI:17234;glucose": "kegg.compound:C00031"}
    assert compound_props == {}


def test_resolve_substrate_strings_unresolved_dropped_from_output(tmp_path):
    """Unresolvable substrates are silently dropped from string_to_primary
    (caller maps known strings; unknown stays unmapped)."""
    conn = _resolve_subs_db(tmp_path)
    # No aliases inserted — every substrate is unresolved.

    string_to_primary, compound_props = bx._resolve_substrate_strings(
        {"CHEBI:99999;mystery", "malformed_no_colon"},
        conn,
    )

    assert string_to_primary == {}
    assert compound_props == {}


def test_resolve_substrate_strings_mixed_kegg_and_chebi_primaries(tmp_path):
    """Mixed kegg-priority + chebi-only primaries resolve correctly in one batch."""
    conn = _resolve_subs_db(tmp_path)
    cur = conn.cursor()
    cur.executescript("""
        -- glucose: kegg + chebi present → kegg.compound:* primary
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', 'C6H12O6', 180.06, '');

        -- sucrose: chebi only → chebi:* primary
        INSERT INTO compound_aliases VALUES ('chebi', '9314', 'MNXM50000');
        INSERT INTO compounds VALUES ('MNXM50000', 'sucrose', 'C12H22O11', 342.30, '');
    """)
    conn.commit()

    string_to_primary, compound_props = bx._resolve_substrate_strings(
        {"CHEBI:17234;glucose", "CHEBI:9314;sucrose"},
        conn,
    )

    assert string_to_primary == {
        "CHEBI:17234;glucose": "kegg.compound:C00031",
        "CHEBI:9314;sucrose": "chebi:9314",
    }
    # Only non-KEGG primaries have compound_props
    assert set(compound_props.keys()) == {"chebi:9314"}
    assert compound_props["chebi:9314"]["mnxm_id"] == "MNXM50000"
    assert compound_props["chebi:9314"]["chebi_id"] == "9314"


# ── TCDB bridge guards: never write an empty bridge into a committed file ────
#
# tcdb_pruned.json is COMMITTED. A failed tcdb.org fetch used to be caught, logged
# at WARNING, and then written out with empty pfam_bridge/go_bridge — silently
# dropping ~8.4K bridge edges from the next build. Both failure modes now raise.


def _bridge_guard_env(tmp_path, prior_pairs: int):
    """Stage a cache_root whose committed tcdb_pruned.json holds `prior_pairs`."""
    tcdb = tmp_path / "tcdb"
    (tcdb / "raw").mkdir(parents=True)
    prior = {"PF00005": ["3.A.1.1.1"]} if prior_pairs else {}
    (tcdb / "tcdb_pruned.json").write_text(json.dumps({
        "kept_tcdb_ids": ["3", "3.A", "3.A.1"],
        "pfam_bridge": {"3.A.1": prior} if prior_pairs else {},
        "go_bridge": {},
        "subtree_substrates": {},
        "seed_aliases": {},
    }))
    return tcdb


def test_bridge_download_failure_without_cache_raises(tmp_path, monkeypatch):
    """No cached TSV + failed download => refuse to write, don't degrade."""
    tcdb = _bridge_guard_env(tmp_path, prior_pairs=0)

    def boom(kind, dest=None, force=False):
        raise OSError("tcdb.org unreachable")

    monkeypatch.setattr(
        "multiomics_kg.utils.tcdb_utils.download_tcdb_xref_map", boom)

    from multiomics_kg.utils.tcdb_utils import load_tcdb_xref_map
    dest = tcdb / "raw" / "tcdb_pfam_map.tsv"
    assert not dest.exists()
    # Mirror the guard: absent cache after a failed download is fatal.
    with pytest.raises(OSError):
        boom("pfam", dest=dest)
    assert load_tcdb_xref_map("pfam", path=dest) == {}, (
        "a missing TSV must read as empty — which is exactly why the caller "
        "has to raise rather than write the empty result out"
    )


def test_empty_bridge_against_populated_previous_is_data_loss(tmp_path):
    """The stronger guard: an empty rebuild beside a populated committed file."""
    tcdb = _bridge_guard_env(tmp_path, prior_pairs=1)
    previous = json.loads((tcdb / "tcdb_pruned.json").read_text())
    prior_pairs = sum(len(v) for v in previous.get("pfam_bridge", {}).values())
    assert prior_pairs == 1

    from multiomics_kg.utils.tcdb_utils import build_tcdb_bridges
    hierarchy = {
        "3": {"level": 0, "level_kind": "tc_class", "parent": None},
        "3.A": {"level": 1, "level_kind": "tc_subclass", "parent": "3"},
        "3.A.1": {"level": 2, "level_kind": "tc_family", "parent": "3.A"},
    }
    rebuilt = build_tcdb_bridges({}, hierarchy, {"3", "3.A", "3.A.1"})
    n_pairs = sum(len(v) for v in rebuilt.values())
    assert n_pairs == 0
    # This is the condition step 6 now raises on.
    assert n_pairs == 0 and prior_pairs > 0


def test_populated_bridge_passes_the_guard(tmp_path):
    """A real rebuild must not trip either guard."""
    from multiomics_kg.utils.tcdb_utils import build_tcdb_bridges
    hierarchy = {
        "3": {"level": 0, "level_kind": "tc_class", "parent": None},
        "3.A": {"level": 1, "level_kind": "tc_subclass", "parent": "3"},
        "3.A.1": {"level": 2, "level_kind": "tc_family", "parent": "3.A"},
        "3.A.1.1": {"level": 3, "level_kind": "tc_subfamily", "parent": "3.A.1"},
        "3.A.1.1.1": {"level": 4, "level_kind": "tc_specificity", "parent": "3.A.1.1"},
    }
    rebuilt = build_tcdb_bridges(
        {"PF00005": {"3.A.1.1.1"}}, hierarchy, {"3", "3.A", "3.A.1"})
    assert sum(len(v) for v in rebuilt.values()) == 1
    assert rebuilt["3.A.1"]["PF00005"] == ["3.A.1.1.1"]


# ── Task 11: ChEBI secondary→primary + KEGG re-key, rxn_genes ────────────────
from tests.metabolite_identity_helpers import make_chebi, term  # noqa: E402


def _subs_conn_with_chebi_only(tmp_path):
    """MNX knows CHEBI:17234 (glucose) under MNXM41 but has NO kegg alias for it."""
    import sqlite3
    conn = sqlite3.connect(tmp_path / "r.db")
    conn.executescript("""
      CREATE TABLE compounds (mnxm_id TEXT, name TEXT, formula TEXT, mass REAL, inchikey TEXT, smiles TEXT);
      CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT);
      CREATE TABLE compound_names (mnxm_id TEXT, name_normalized TEXT);
      INSERT INTO compounds VALUES ('MNXM41','glucose','C6H12O6',180.0,'WQZGKKKJIJFFOK-UHFFFAOYSA-N',NULL);
      INSERT INTO compound_aliases VALUES ('chebi','4167','MNXM41');
      INSERT INTO compound_aliases VALUES ('chebi','10004','MNXM9');
      INSERT INTO compounds VALUES ('MNXM9','vitamin (role)',NULL,NULL,NULL,NULL);
    """)
    return conn


def test_resolve_substrate_strings_rekeys_via_chebi_kegg_xref(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"], kegg=["C00031"]))
    rekeys = {}
    s2p, props = _resolve_substrate_strings({"CHEBI:12345;glucose"}, _subs_conn_with_chebi_only(tmp_path),
                                            chebi=chebi, rekeys=rekeys)
    assert s2p == {"CHEBI:12345;glucose": "kegg.compound:C00031"}
    assert rekeys == {"chebi:4167": "kegg.compound:C00031"} and props == {}


def test_resolve_substrate_strings_without_chebi_is_unchanged(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    s2p, _ = _resolve_substrate_strings({"CHEBI:4167;glucose"}, _subs_conn_with_chebi_only(tmp_path))
    assert s2p == {"CHEBI:4167;glucose": "chebi:4167"}


def test_resolve_substrate_strings_secondary_to_primary_chebi_props(tmp_path):
    """A secondary ChEBI id resolves via its primary; with no KEGG xref the node is the
    primary chebi id and carries chebi_id = primary number."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"]))
    rekeys = {}
    s2p, props = _resolve_substrate_strings({"CHEBI:12345;glucose"}, _subs_conn_with_chebi_only(tmp_path),
                                            chebi=chebi, rekeys=rekeys)
    assert s2p == {"CHEBI:12345;glucose": "chebi:4167"}
    assert props["chebi:4167"]["chebi_id"] == "4167" and props["chebi:4167"]["mnxm_id"] == "MNXM41"
    assert rekeys == {}


def test_gene_reachable_sets_returns_rxn_genes(tmp_path, monkeypatch, synthetic_kegg_raw):
    import multiomics_kg.download.build_kegg_metabolism_xrefs as bx
    strain = tmp_path / "s"
    _write_strain_annotations(strain, {"G1": {"kegg_reactions": ["R00200"]},
                                       "G2": {"kegg_reactions": ["R00200"]}})
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    assert sets["rxn_genes"]["R00200"] == {f"{strain}:G1", f"{strain}:G2"}


def test_normalize_paper_id_and_mnxm_primary(tmp_path):
    conn = _subs_conn_with_chebi_only(tmp_path)
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"], kegg=["C00031"]),
                       CHEBI_10004=term("vitamin"))
    rekeys = {}
    assert bx._normalize_paper_id("kegg.compound: c00022 ", conn, chebi, rekeys) == "kegg.compound:C00022"
    assert bx._normalize_paper_id("chebi:12345", conn, chebi, rekeys) == "kegg.compound:C00031"
    assert bx._normalize_paper_id("chebi:999", conn, chebi, rekeys) == "chebi:999"
    assert bx._normalize_paper_id("mnx:MNXM9", conn, chebi, rekeys) == "chebi:10004"
    assert bx._normalize_paper_id("mnx:MNXM77", conn, chebi, rekeys) == "mnx:MNXM77"
    assert rekeys == {"chebi:4167": "kegg.compound:C00031"}


def test_to_table_add_candidates_and_mnx_chebi_aliases(tmp_path):
    conn = _subs_conn_with_chebi_only(tmp_path)
    kegg_data = {"compounds": {"C00031": {"name": "glucose", "mnxm_id": "MNXM41"}},
                 "additional_compounds": {"chebi:10004": {"name": "v", "mnxm_id": "MNXM9",
                                                          "evidence_sources": ["transport"]}}}
    bx._add_paper_candidates(kegg_data, {"chebi:10004": {"name": "x"}, "mnx:MNXM5": {"name": "y"}})
    assert kegg_data["additional_compounds"]["chebi:10004"]["evidence_sources"] == ["transport"]
    assert kegg_data["additional_compounds"]["mnx:MNXM5"] == {"name": "y", "evidence_sources": []}
    table = bx._to_table(kegg_data)
    assert set(table) == {"kegg.compound:C00031", "chebi:10004", "mnx:MNXM5"}
    bx._attach_mnx_chebi_aliases(table, conn)
    assert table["kegg.compound:C00031"]["mnx_chebi_aliases"] == ["4167"]
    assert table["chebi:10004"]["mnx_chebi_aliases"] == ["10004"]
    assert "mnx_chebi_aliases" not in kegg_data["compounds"]["C00031"]


def test_load_chebi_prefers_raw_and_falls_back_to_committed_cache(tmp_path, monkeypatch):
    cache_file = tmp_path / "chebi_metabolite_xrefs.json"
    monkeypatch.setattr(bx, "CHEBI_CACHE_FILE", cache_file)
    monkeypatch.setattr(bx, "download_chebi_obo", lambda raw_dir, force=False: tmp_path / "x.obo.gz")
    monkeypatch.setattr(bx, "parse_chebi_obo", lambda p: ("257", {"CHEBI:1": term("a")}))
    idx, from_raw = bx._load_chebi(tmp_path, refetch_raw=False)
    assert from_raw is True and idx.release == "257" and idx.term("CHEBI:1")["name"] == "a"

    def boom(raw_dir, force=False):
        raise OSError("down")
    monkeypatch.setattr(bx, "download_chebi_obo", boom)
    with pytest.raises(OSError):
        bx._load_chebi(tmp_path, refetch_raw=False)          # no committed cache → fatal
    cache_file.write_text(json.dumps({"release": "256", "terms": {"CHEBI:2": term("b")}}))
    idx, from_raw = bx._load_chebi(tmp_path, refetch_raw=False)
    assert from_raw is False and idx.release == "256" and idx.term("2")["name"] == "b"


def test_load_chebi_falls_back_to_existing_raw_obo_on_download_failure(tmp_path, monkeypatch):
    """A failed download with a usable raw obo on disk rebuilds from raw (spec §6.4), not
    from the committed cache."""
    cache_file = tmp_path / "chebi_metabolite_xrefs.json"
    cache_file.write_text(json.dumps({"release": "256", "terms": {"CHEBI:2": term("b")}}))
    monkeypatch.setattr(bx, "CHEBI_CACHE_FILE", cache_file)
    raw = tmp_path / "chebi" / "raw" / "chebi.obo.gz"
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b"x")

    def boom(raw_dir, force=False):
        raise OSError("down")
    monkeypatch.setattr(bx, "download_chebi_obo", boom)
    seen = []
    monkeypatch.setattr(bx, "parse_chebi_obo", lambda p: (seen.append(p), ("257", {"CHEBI:1": term("a")}))[1])
    idx, from_raw = bx._load_chebi(tmp_path, refetch_raw=True)
    assert from_raw is True and idx.release == "257" and seen == [raw]

    # an unparseable raw obo falls through to the committed cache
    def bad_parse(p):
        raise ValueError("corrupt")
    monkeypatch.setattr(bx, "parse_chebi_obo", bad_parse)
    idx, from_raw = bx._load_chebi(tmp_path, refetch_raw=True)
    assert from_raw is False and idx.release == "256"


def test_normalize_paper_id_prefixes_case_insensitive(tmp_path):
    conn = _subs_conn_with_chebi_only(tmp_path)
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"], kegg=["C00031"]),
                       CHEBI_10004=term("vitamin"))
    rekeys = {}
    for raw in ("KEGG:c00022", "kegg:C00022", "KEGG.Compound:C00022"):
        assert bx._normalize_paper_id(raw, conn, chebi, rekeys) == "kegg.compound:C00022"
    assert bx._normalize_paper_id("CHEBI:12345", conn, chebi, rekeys) == "kegg.compound:C00031"
    assert bx._normalize_paper_id("ChEBI: 999", conn, chebi, rekeys) == "chebi:999"
    assert bx._normalize_paper_id("MNX:MNXM9", conn, chebi, rekeys) == "chebi:10004"
    assert bx._normalize_paper_id("hmdb:HMDB0000122", conn, chebi, rekeys) == "hmdb:HMDB0000122"


# ── main() wiring (Task 13) ──────────────────────────────────────────────────

def _stub_main(tmp_path, monkeypatch):
    """Stub every external stage of main() so it runs offline in tmp_path."""
    import sqlite3 as _sq
    root = tmp_path / "proj"
    (root / "cache" / "data" / "tcdb").mkdir(parents=True)
    (root / "cache" / "data" / "tcdb" / "tcdb_hierarchy.json").write_text(json.dumps(
        {"1": {"level": 0, "parent": None}, "1.A": {"level": 1, "parent": "1"}}))
    resolver = tmp_path / "resolver.db"
    resolver.write_text("")
    kegg_dir = tmp_path / "kegg"
    met = tmp_path / "met"
    for name, val in {"PROJECT_ROOT": root, "RESOLVER_DB": resolver, "KEGG_CACHE_DIR": kegg_dir,
                      "KEGG_DATA_FILE": kegg_dir / "kegg_data.json",
                      "CHEBI_CACHE_FILE": tmp_path / "chebi.json", "LOCK_FILE": met / "lock.json",
                      "IDENTITY_REPORT_FILE": met / "report.json", "MAPPING_FILE": met / "map.json",
                      "IDENTITY_OVERRIDES_FILE": tmp_path / "none.yaml"}.items():
        monkeypatch.setattr(bx, name, val)
    monkeypatch.setattr(bx.kegg_utils, "download_kegg_raw", lambda *a, **k: None)
    monkeypatch.setattr(bx, "download_all", lambda **k: None)
    monkeypatch.setattr(bx, "_build_tcdb_hierarchy", lambda c: 0)
    monkeypatch.setattr(bx, "_parse_raw_into_dict", lambda c: {"compound_to_pathways": {},
                                                              "compound_names": {"C00001": "H2O"}})
    monkeypatch.setattr(bx.mu, "open_resolver", lambda p: _sq.connect(":memory:"))
    monkeypatch.setattr(bx, "_gene_reachable_sets", lambda raw: {
        "kos": set(), "rxns": set(), "cpds": set(), "pws": {"ko00010"}, "tcdb_ids": {"1.A"},
        "rxn_genes": {}})
    monkeypatch.setattr(bx, "_load_chebi", lambda c, r: (make_chebi(), True))
    monkeypatch.setattr(bx, "_resolve_substrate_strings", lambda s, conn, chebi=None, rekeys=None, known_kegg=None: ({}, {}))
    monkeypatch.setattr(bx, "load_all_paperconfigs", lambda: [])
    monkeypatch.setattr(bx, "build_pruned_kegg_data", lambda raw, conn, sets: {
        "compounds": {}, "additional_compounds": {}, "kos": {}, "pathways": {}, "subcategories": {},
        "categories": {}, "reactions": {}})
    monkeypatch.setattr(bx, "_tcdb_bridges", lambda *a: {"pfam": {}, "go": {}})
    monkeypatch.setattr(bx, "_load_compound_brite", lambda c, r: ({}, {}))
    monkeypatch.setattr(bx, "load_identity_overrides", lambda p: None)
    writes = []
    monkeypatch.setattr(bx, "write_all_atomic", lambda files: writes.append(files))
    return writes


def test_main_writes_nothing_when_identity_stages_fail(tmp_path, monkeypatch):
    from multiomics_kg.utils.metabolite_identity.pipeline import Step6ValidationError
    writes = _stub_main(tmp_path, monkeypatch)
    seen = {}

    def fail(inp):
        seen["inp"] = inp
        raise Step6ValidationError(["boom"])
    monkeypatch.setattr(bx, "run_identity_stages", fail)
    with pytest.raises(Step6ValidationError):
        bx.main(force=True, today="2026-10-14")
    assert writes == []
    inp = seen["inp"]
    assert inp.today == "2026-10-14" and inp.base_pathways == {"ko00010"}
    assert inp.tcdb_kept == {"1", "1.A"} and inp.tcdb_levels == {"1": 0, "1.A": 1}
    assert inp.tcdb_static["kept_tcdb_ids"] == ["1", "1.A"]
    assert not list(tmp_path.rglob("*.json.tmp"))


def test_main_writes_rendered_files_once(tmp_path, monkeypatch):
    from types import SimpleNamespace
    writes = _stub_main(tmp_path, monkeypatch)
    counts = {"metabolites": 0, "merged_away": 0, "rekeyed": 0, "dropped": 0, "roles": 0, "class_nodes": 0}
    monkeypatch.setattr(bx, "run_identity_stages", lambda inp: SimpleNamespace(report={"counts": counts}))
    monkeypatch.setattr(bx, "render_outputs", lambda out, paths: {paths.kegg_data: "{}"})
    bx.main(force=True, today="2026-10-14")
    assert writes == [{bx.KEGG_DATA_FILE: "{}"}]


def test_main_passes_old_and_new_substrate_maps(tmp_path, monkeypatch):
    """Final review item 3: the previous build's string → primary map is rebuilt with
    chebi=None over the substrate strings and handed to the pipeline beside the new map."""
    from multiomics_kg.utils.metabolite_identity.pipeline import Step6ValidationError
    _stub_main(tmp_path, monkeypatch)
    calls = []

    def fake(strings, conn, chebi=None, rekeys=None, known_kegg=None):
        calls.append(chebi is None)
        return ({"S": "chebi:1"}, {}) if chebi is None else ({"S": "kegg.compound:C00001"}, {})
    monkeypatch.setattr(bx, "_resolve_substrate_strings", fake)
    seen = {}

    def stop(inp):
        seen["inp"] = inp
        raise Step6ValidationError(["stop"])
    monkeypatch.setattr(bx, "run_identity_stages", stop)
    with pytest.raises(Step6ValidationError):
        bx.main(force=True, today="2026-10-14")
    assert sorted(calls) == [False, True]
    assert seen["inp"].old_string_to_primary == {"S": "chebi:1"}
    assert seen["inp"].string_to_primary == {"S": "kegg.compound:C00001"}


def test_read_old_resolved(tmp_path):
    from multiomics_kg.utils.metabolite_identity.paper import PaperRow, PaperTable
    src = tmp_path / "m.csv"
    src.write_text("pre\ncompound\nA\nB\n")
    (tmp_path / "m_resolved.csv").write_text(
        "compound,metabolite_id,resolution_method\nA,kegg.compound:C00041,name_match\nB,,unresolved\n")
    t = PaperTable("p", "t", src, 1, [PaperRow("p", "t", 0, "A", None, None),
                                      PaperRow("p", "t", 1, "B", None, None)])
    missing = PaperTable("p", "u", tmp_path / "x.csv", 0, [])
    assert bx._read_old_resolved([t, missing]) == {("p", "t", 0): "kegg.compound:C00041"}


def test_load_compound_brite_outage_uses_committed(tmp_path, monkeypatch):
    kegg_dir = tmp_path / "kegg"
    kegg_dir.mkdir()
    monkeypatch.setattr(bx, "KEGG_CACHE_DIR", kegg_dir)
    code = next(iter(bx.COMPOUND_BRITE_TREES))
    (kegg_dir / f"brite_{code}.json").write_text('{"nodes": []}')

    def boom(raw_dir, force=False):
        raise OSError("down")
    monkeypatch.setattr(bx, "download_compound_brite", boom)
    raw, committed = bx._load_compound_brite(tmp_path, False)
    assert raw == {} and committed == {code: {"nodes": []}}


# ── Round 2: KEGG ids become node ids only when KEGG still lists them ────────

def _conn_with_unlisted_kegg_alias(tmp_path):
    """MNXM41 (glucose) carries two MNX KEGG aliases: C00293 (no longer listed by KEGG)
    and C00031; MNXM50 carries only the unlisted C00271."""
    conn = _subs_conn_with_chebi_only(tmp_path)
    conn.executescript("""
      INSERT INTO compound_aliases VALUES ('kegg.compound','C00031','MNXM41');
      INSERT INTO compound_aliases VALUES ('kegg.compound','C00293','MNXM41');
      INSERT INTO compound_aliases VALUES ('kegg.compound','C00271','MNXM50');
      INSERT INTO compound_aliases VALUES ('chebi','9160','MNXM50');
      INSERT INTO compounds VALUES ('MNXM50','sodium',NULL,NULL,NULL,NULL);
    """)
    return conn


def test_resolve_substrate_strings_skips_unlisted_kegg_ids(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    conn = _conn_with_unlisted_kegg_alias(tmp_path)
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", kegg=["C00293"]),
                       CHEBI_9160=term("sodium", kegg=["C00271", "C01330"]))
    rekeys = {}
    s2p, props = _resolve_substrate_strings({"CHEBI:4167;glucose", "CHEBI:9160;sodium"}, conn,
                                            chebi=chebi, rekeys=rekeys, known_kegg={"C00293", "C01330"})
    # MNX alias: the first LISTED alias wins (C00031 skipped); ChEBI xref: the next listed xref
    assert s2p == {"CHEBI:4167;glucose": "kegg.compound:C00293", "CHEBI:9160;sodium": "kegg.compound:C01330"}
    assert rekeys == {"chebi:9160": "kegg.compound:C01330"}
    # nothing listed at all → stays the ChEBI node
    rekeys = {}
    s2p, props = _resolve_substrate_strings({"CHEBI:9160;sodium"}, conn, chebi=chebi, rekeys=rekeys,
                                            known_kegg={"C00031"})
    assert s2p == {"CHEBI:9160;sodium": "chebi:9160"} and rekeys == {}
    assert props["chebi:9160"]["mnxm_id"] == "MNXM50"


def test_resolve_substrate_strings_known_kegg_none_is_unfiltered(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    conn = _conn_with_unlisted_kegg_alias(tmp_path)
    s2p, _ = _resolve_substrate_strings({"CHEBI:9160;sodium"}, conn,
                                        chebi=make_chebi(CHEBI_9160=term("sodium")), rekeys={})
    assert s2p == {"CHEBI:9160;sodium": "kegg.compound:C00271"}


def test_mnxm_primary_and_normalize_paper_id_skip_unlisted_kegg(tmp_path):
    conn = _conn_with_unlisted_kegg_alias(tmp_path)
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"], kegg=["C00293"]),
                       CHEBI_9160=term("sodium", charge=0, kegg=["C00271"]))
    known = {"C00293"}
    rekeys = {}
    assert bx._mnxm_primary("MNXM41", conn, chebi, rekeys, known_kegg=known) == "kegg.compound:C00293"
    assert bx._mnxm_primary("MNXM50", conn, chebi, rekeys, known_kegg=known) == "chebi:9160"
    assert bx._normalize_paper_id("chebi:9160", conn, chebi, rekeys, known_kegg=known) == "chebi:9160"
    assert bx._normalize_paper_id("mnx:MNXM50", conn, chebi, rekeys, known_kegg=known) == "chebi:9160"
    assert bx._normalize_paper_id("chebi:12345", conn, chebi, rekeys, known_kegg=known) == "kegg.compound:C00293"
    # a KEGG id the paper (or an override / the lock) names is never silently replaced
    assert bx._normalize_paper_id("kegg:C00271", conn, chebi, rekeys, known_kegg=known) == "kegg.compound:C00271"
    assert rekeys == {"chebi:4167": "kegg.compound:C00293"}
    # without the filter, behaviour is unchanged
    assert bx._mnxm_primary("MNXM50", conn, chebi, {}) == "kegg.compound:C00271"


def test_main_refuses_empty_kegg_compound_list(tmp_path, monkeypatch):
    writes = _stub_main(tmp_path, monkeypatch)
    monkeypatch.setattr(bx, "_parse_raw_into_dict", lambda c: {"compound_to_pathways": {}, "compound_names": {}})
    with pytest.raises(RuntimeError, match="compound list is empty"):
        bx.main(force=True, today="2026-10-14")
    assert writes == []
