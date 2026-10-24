"""
Shared fixtures for KG validity tests.

Tests in this directory connect to a live Neo4j instance (the deployed Docker graph).
They are skipped automatically if Neo4j is not reachable.

Run with:
    pytest tests/kg_validity/ -v
    pytest tests/kg_validity/ -v --neo4j-url bolt://localhost:7688
    pytest -m kg
"""

import pytest
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable, AuthError


def pytest_addoption(parser):
    parser.addoption(
        "--neo4j-url",
        default="bolt://localhost:7687",
        help="Bolt URL for the Neo4j instance to validate (default: bolt://localhost:7687)",
    )


@pytest.fixture(scope="session")
def neo4j_driver(request):
    url = request.config.getoption("--neo4j-url")
    try:
        # Auth is disabled in docker-compose (NEO4J_dbms_security_auth__enabled=false)
        driver = GraphDatabase.driver(url, auth=None)
        driver.verify_connectivity()
    except (ServiceUnavailable, AuthError, Exception) as e:
        pytest.skip(f"Neo4j not available at {url}: {e}")
        return
    yield driver
    driver.close()


@pytest.fixture(scope="session")
def run_query(neo4j_driver):
    """Helper: run a Cypher query and return results as a list of dicts."""
    def _run(cypher, **params):
        with neo4j_driver.session() as session:
            return session.run(cypher, **params).data()
    return _run


@pytest.fixture(scope="session")
def kegg_cache_counts():
    """Metabolite / pathway totals the build derives from the committed
    cache/data/kegg/kegg_data.json, so count tests track the cache instead of
    pinning literals (adding a paper's compound must not break them).

    - metabolites: compounds + additional_compounds (metabolism_adapter nodes)
    - pathways: pathway keys (MultiKeggAnnotationAdapter pathway KeggTerm nodes)
    - compound_pathway_pairs: distinct (compound, pathway) pairs whose pathway is a
      node — the Metabolite_in_pathway edge set, = sum(KeggTerm.metabolite_count)
    """
    import json
    import pathlib

    path = pathlib.Path(__file__).resolve().parents[2] / "cache/data/kegg/kegg_data.json"
    if not path.exists():
        pytest.skip(f"{path} missing; run prepare_data step 6")
    data = json.loads(path.read_text())
    pathways = set(data["pathways"])
    pairs = {(cpd, pw) for cpd, v in data["compounds"].items()
             for pw in v.get("pathways", []) if pw in pathways}
    return {
        "metabolites": len(data["compounds"]) + len(data.get("additional_compounds", {})),
        "pathways": len(pathways),
        "compound_pathway_pairs": len(pairs),
    }
