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
