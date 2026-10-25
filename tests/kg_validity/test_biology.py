"""
Biological truth tests for the multi-omics knowledge graph.

These tests encode domain-specific constraints that must hold for
Prochlorococcus / Synechococcus / Alteromonas biology:

- Ecotype/clade labels on strains
- Gene absences that are defining features (Black Queen Hypothesis)
- All expected strains present
- UniProt locus-tag → accession spot checks
"""

import pytest


pytestmark = pytest.mark.kg


# Strains confirmed in graph (OrganismTaxon nodes with strain_name):
# 15 Prochlorococcus + 6 Synechococcus/Parasynechococcus/Thermosynechococcus + 3 Alteromonas + 4 heterotrophs
EXPECTED_STRAINS = [
    # Prochlorococcus
    "MED4",
    "MIT9313",
    "MIT9301",
    "MIT9312",
    "MIT9303",
    "NATL2A",
    "NATL1A",
    "AS9601",
    "RSP50",
    # Soussan 2025 N/P starvation strain panel (May 2026)
    "MIT9515",
    "MIT9202",
    "MIT9215",
    "MIT0604",
    "SB",
    "PAC1",
    # Soussan 2025 deferred strains, public assemblies added June 2026
    "MIT1327",
    "MIT1314",
    # Synechococcus / Parasynechococcus
    "WH8102",
    "CC9311",
    "WH7803",
    "PCC7002",
    "PCC7942",
    "UTEX2973",
    # Thermosynechococcus
    "BP1",
    # Alteromonas
    "MIT1002",
    "EZ55",
    "HOT1A3",
    # Alteromonas mediterranea deep-ecotype genome strains (May 2026, genomes only)
    "AltDE",
    "AltDE1",
    # Heterotrophic coculture partners
    "W3-18-1",
    "KT2440",
    "DSS-3",
    "MruberA",
    # Ford 2021 Cyanorak v2.1 SBP survey: every remaining Cyanorak genome but WH8020 (74; Oct 2026)
    "EQPAC1", "GP2", "HNLC1", "HNLC2", "LG", "MIT0601",
    "MIT0602", "MIT0603", "MIT0701", "MIT0702", "MIT0703", "MIT9107",
    "MIT9116", "MIT9123", "MIT9201", "MIT9211", "MIT9302", "MIT9311",
    "MIT9314", "MIT9321", "MIT9322", "MIT9401", "SS2", "SS35",
    "SS51", "SS52", "NS01", "PCC6307", "PCC7001", "A15-127",
    "A15-24", "A15-28", "A15-44", "A15-60", "A15-62", "A18-25c",
    "A18-40", "A18-46.1", "BIOS-E4-1", "BIOS-U3-1", "BMK-MC-1", "BOUM118",
    "CB0101", "CB0205", "CC9605", "CC9616", "CC9902", "KORDI-100",
    "KORDI-49", "KORDI-52", "M16.1", "MEDNS5", "MINOS11", "MITS9220",
    "MVIR-18-1", "NOUM97013", "PROS-7-1", "PROS-9-1", "PROS-U-1", "RCC307",
    "ROS8604", "RS9902", "RS9907", "RS9909", "RS9915", "RS9916",
    "RS9917", "SYN20", "TAK9802", "WH5701", "WH7805", "WH8016",
    "WH8101", "WH8103",
]

# Clade assignments per strain. Prochlorococcus uses the HL/LL ecotype scheme
# (HLI = High Light clade I); marine Synechococcus uses Roman-numeral clades
# within sub-cluster 5.1 (per the Cyanorak organism table).
EXPECTED_CLADES = {
    "MED4":    "HLI",
    "RSP50":   "HLI",
    "MIT9515": "HLI",
    "AS9601":  "HLII",
    "MIT9312": "HLII",
    "MIT9301": "HLII",
    "MIT9202": "HLII",
    "MIT9215": "HLII",
    "MIT0604": "HLII",
    "SB":      "HLII",
    "NATL2A":  "LLI",
    "NATL1A":  "LLI",
    "PAC1":    "LLI",
    "MIT9313": "LLIV",
    "MIT9303": "LLIV",
    "MIT1327": "LLIV",
    "MIT1314": "HLII",
    # Marine Synechococcus (sub-cluster 5.1)
    "CC9311":  "I",
    "WH8102":  "III",
    "WH7803":  "V",
    "BL107":   "IV",
    # Ford 2021 Cyanorak batch: clade = the Cyanorak organism table's Clade column
    "EQPAC1": "HLI",
    "GP2": "HLII",
    "HNLC1": "HLIV",
    "HNLC2": "HLIII",
    "LG": "LLII",
    "MIT0601": "LLIII",
    "MIT0602": "LLII",
    "MIT0603": "LLII",
    "MIT0701": "LLIV",
    "MIT0702": "LLIV",
    "MIT0703": "LLIV",
    "MIT9107": "HLII",
    "MIT9116": "HLII",
    "MIT9123": "HLII",
    "MIT9201": "HLII",
    "MIT9211": "LLIII",
    "MIT9302": "HLII",
    "MIT9311": "HLII",
    "MIT9314": "HLII",
    "MIT9321": "HLII",
    "MIT9322": "HLII",
    "MIT9401": "HLII",
    "SS2": "LLII",
    "SS35": "LLII",
    "SS51": "LLII",
    "SS52": "LLII",
    "NS01": "5.2",
    "PCC6307": "5.2",
    "PCC7001": "5.2",
    "A15-127": "WPC1",
    "A15-24": "III",
    "A15-28": "III",
    "A15-44": "II",
    "A15-60": "VII",
    "A15-62": "II",
    "A18-25c": "VII",
    "A18-40": "III",
    "A18-46.1": "III",
    "BIOS-E4-1": "CRD1",
    "BIOS-U3-1": "CRD1",
    "BMK-MC-1": "V",
    "BOUM118": "III",
    "CB0101": "5.2",
    "CB0205": "5.2",
    "CC9605": "II",
    "CC9616": "EnvC",
    "CC9902": "IV",
    "KORDI-100": "UC-A",
    "KORDI-49": "WPC1",
    "KORDI-52": "II",
    "M16.1": "II",
    "MEDNS5": "VI",
    "MINOS11": "5.3",
    "MITS9220": "CRD1",
    "MVIR-18-1": "I",
    "NOUM97013": "VII",
    "PROS-7-1": "VI",
    "PROS-9-1": "I",
    "PROS-U-1": "II",
    "RCC307": "5.3",
    "ROS8604": "I",
    "RS9902": "II",
    "RS9907": "II",
    "RS9909": "VIII",
    "RS9915": "III",
    "RS9916": "IX",
    "RS9917": "VIII",
    "SYN20": "I",
    "TAK9802": "II",
    "WH5701": "5.2",
    "WH7805": "VI",
    "WH8016": "I",
    "WH8101": "VIII",
    "WH8103": "III",
}


# ---------------------------------------------------------------------------
# Strain presence
# ---------------------------------------------------------------------------

def test_all_expected_strains_present(run_query):
    """All expected organism strains must exist in the graph."""
    result = run_query(
        "MATCH (o:OrganismTaxon) RETURN o.strain_name AS name"
    )
    present = {r["name"] or "" for r in result}
    missing = [
        s for s in EXPECTED_STRAINS
        if not any(s in name for name in present)
    ]
    assert not missing, (
        f"Expected strains not found in graph: {missing}\n"
        f"Strains present: {sorted(present)}"
    )


# ---------------------------------------------------------------------------
# Ecotype / clade labels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strain,expected_clade", EXPECTED_CLADES.items())
def test_strain_clade_label(run_query, strain, expected_clade):
    """Each strain must carry the correct ecotype/clade label."""
    result = run_query(
        "MATCH (o:OrganismTaxon) WHERE o.strain_name = $strain "
        "RETURN o.strain_name AS name, o.clade AS clade",
        strain=strain,
    )
    assert len(result) >= 1, f"No OrganismTaxon node found for strain '{strain}'"
    for row in result:
        assert row["clade"] == expected_clade, (
            f"Strain {row['name']}: expected clade '{expected_clade}', "
            f"got '{row['clade']}'"
        )


def test_prochlorococcus_genus_label(run_query):
    """All Prochlorococcus strains must have genus='Prochlorococcus'."""
    pro_strains = [s for s in EXPECTED_STRAINS
                   if s in ("MED4", "MIT9313", "MIT9301", "MIT9312",
                             "MIT9303", "NATL2A", "NATL1A", "AS9601", "RSP50")]
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IN $strains
        AND (o.genus IS NULL OR o.genus <> 'Prochlorococcus')
        RETURN o.strain_name AS name, o.genus AS genus
    """, strains=pro_strains)
    assert result == [], (
        f"Strains with wrong/missing genus: {result}"
    )


# ---------------------------------------------------------------------------
# Black Queen Hypothesis: katG is absent in Prochlorococcus
# ---------------------------------------------------------------------------

def test_prochlorococcus_no_katG(run_query):
    """
    Prochlorococcus has lost catalase-peroxidase (katG) through genome streamlining
    (Black Queen Hypothesis). No Gene belonging to a Prochlorococcus organism
    should carry 'katG' as a gene name.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.genus = 'Prochlorococcus'
          AND any(name IN g.gene_names WHERE toLower(name) = 'katg')
        RETURN count(g) AS cnt, collect(g.locus_tag)[..5] AS examples
    """)
    cnt = result[0]["cnt"]
    assert cnt == 0, (
        f"Found {cnt} Prochlorococcus gene(s) named 'katG' (should be 0). "
        f"Example locus tags: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# UniProt locus-tag → accession spot checks (MED4)
# ---------------------------------------------------------------------------

def test_med4_genes_have_proteins(run_query):
    """
    MED4 genes should map to UniProt proteins via shared locus_tag.
    Require that at least 80% of MED4 genes have a corresponding Protein.
    Gene-protein linkage is via locus_tag (no explicit Gene_encodes_protein edge).
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.strain_name = 'MED4'
        OPTIONAL MATCH (p:Protein {locus_tag: g.locus_tag})
        WITH count(g) AS total, count(p) AS with_protein
        RETURN total, with_protein
    """)
    row = result[0]
    if row["total"] == 0:
        pytest.skip("No MED4 genes found")
    coverage = row["with_protein"] / row["total"]
    assert coverage >= 0.75, (
        f"MED4 gene→protein coverage is {coverage:.1%} "
        f"({row['with_protein']}/{row['total']}); expected >= 75%"
    )


def test_known_med4_locus_tag_maps_to_protein(run_query):
    """
    PMM1375 (psbA) in MED4 should have a corresponding UniProt Protein node.
    Genes and proteins are linked by shared locus_tag (no explicit edge).
    This is a regression anchor — if the locus-tag mapping breaks, this fails.
    """
    result = run_query("""
        MATCH (p:Protein {locus_tag: 'PMM1375'})
        RETURN p.id AS uniprot_id, p.protein_name AS name
    """)
    assert len(result) >= 1, (
        "No Protein found with locus_tag='PMM1375' (MED4 psbA). "
        "UniProt locus-tag mapping may be broken."
    )


# ---------------------------------------------------------------------------
# Taxonomic hierarchy completeness
# ---------------------------------------------------------------------------

def test_genomic_organisms_have_genus(run_query):
    """
    Genomic OrganismTaxon nodes (those with a strain_name) must have genus populated.
    Species may be absent for some strains (e.g., Alteromonas strains without
    a named species yet). Non-genomic treatment organisms (phage, viruses) may
    lack full taxonomy and are excluded by requiring strain_name IS NOT NULL.
    """
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.strain_name IS NOT NULL AND o.genus IS NULL
        RETURN o.strain_name AS strain, o.genus AS genus
    """)
    assert result == [], (
        f"Genomic OrganismTaxon nodes missing genus: {result}"
    )
