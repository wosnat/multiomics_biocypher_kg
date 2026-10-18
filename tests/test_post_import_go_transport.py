"""Static checks: GO-route organism rollup + Schema_info exclusion list."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SH = (ROOT / "scripts" / "post-import.sh").read_text()
CYPHER = (ROOT / "scripts" / "post-import.cypher").read_text()
MARKER = "// GO transport organism rollup"


def _block(text: str) -> str:
    start = text.index(MARKER)
    end = text.index("IN TRANSACTIONS", start)
    return text[start:end].strip()


def test_go_transport_rollup_identical_in_both_files():
    assert MARKER in SH and MARKER in CYPHER
    assert _block(SH) == _block(CYPHER)
    assert "go_transport_organism_names" in _block(SH)


def test_schema_info_sets_excluded_targets():
    assert "s.go_transport_excluded_targets" in SH
    assert "go_transport_excluded_targets.json" in SH
    assert "go_excluded" in SH
