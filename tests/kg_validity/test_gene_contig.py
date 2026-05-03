"""KG validity: every Gene must have a non-null contig.

F4 ship 4.1 requires the contig (NCBI GFF seqid) on every Gene node so
genomic-neighbor lookups have a valid discriminator.
"""

import pytest


pytestmark = pytest.mark.kg


def test_every_positioned_gene_has_contig(run_query):
    """Every Gene node with positional data must have a non-null contig.

    Some Cyanorak-only orphan records lack start/end (and therefore no contig).
    Filter to genes with start populated — the contig invariant only applies
    where genomic-neighbor lookups would actually be performed.
    """
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.start IS NOT NULL
          AND (g.contig IS NULL OR g.contig = '')
        RETURN count(*) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} positioned Gene nodes have null/empty contig"


def test_every_organism_has_at_least_one_contig(run_query):
    """Every organism must have at least one distinct contig across its genes.

    Single-chromosome organisms have 1 contig; multi-replicon organisms have
    more, but every organism must have at least one.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WITH o.preferred_name AS organism, collect(DISTINCT g.contig) AS contigs
        WHERE size(contigs) = 0
        RETURN organism, contigs
    """)
    assert result == [], f"organisms with no distinct contig: {result}"


def test_med4_has_known_chromosome(run_query):
    """Spot-check: MED4 genes should land on its known NCBI accession.

    Pinning the actual value (NC_005072.1) catches both null contigs (the
    pre-rebuild failure mode) and accession-swap regressions.
    """
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.preferred_name = 'Prochlorococcus MED4'
          AND g.start IS NOT NULL
        RETURN DISTINCT g.contig AS contig
    """)
    contigs = {row["contig"] for row in result}
    # MED4 single-replicon — expect exactly the known NCBI chromosome accession
    assert contigs == {'NC_005072.1'}, f"MED4 expected {{'NC_005072.1'}}, got: {contigs}"
