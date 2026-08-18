"""Build-level gates for create_knowledge_graph.py.

Two modes share one helper: `--test` (100 items/adapter) runs by default so a
broken build surfaces in `pytest -m "not slow and not kg"`; the full build stays
behind @pytest.mark.slow. Both need a populated cache/ -- neither is portable to
a fresh checkout.

Output-dir note: `create_knowledge_graph.py --output-dir` only controls where
the optional `--go`/`--ec` CSV exports land; it does NOT change where BioCypher
writes its main node/edge CSVs. Those always go to `<repo>/biocypher-out/
<timestamp>/` (BioCypher's own default output_directory, since neither
`config/biocypher_config.yaml` nor `create_knowledge_graph.py` overrides it).
So `_run_build` locates the real output by diffing `biocypher-out/` before and
after the subprocess runs, rather than trusting `--output-dir`.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from multiomics_kg.utils.controlled_vocab import load_vocabularies

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BIOCYPHER_OUT = PROJECT_ROOT / "biocypher-out"


def _run_build(tmp_path: Path, test_mode: bool) -> tuple[Path, str]:
    """Run the build once; return (biocypher_out_dir, combined stdout+stderr).

    `biocypher_out_dir` is the real BioCypher CSV output directory, discovered
    by diffing `biocypher-out/` (see module docstring) rather than trusting
    `--output-dir`, which does not govern it.
    """
    before = set(BIOCYPHER_OUT.iterdir()) if BIOCYPHER_OUT.exists() else set()

    # --output-dir only matters for --go/--ec CSV exports (not used here); it
    # still needs a writable path because the script unconditionally
    # os.makedirs()s it.
    unused_go_ec_dir = tmp_path / "go_ec_csv"
    cmd = [sys.executable, "create_knowledge_graph.py",
           "--output-dir", str(unused_go_ec_dir)]
    if test_mode:
        cmd.append("--test")
    result = subprocess.run(
        cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=3600,
    )
    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    assert result.returncode == 0, (
        f"Build exited {result.returncode}.\nstderr:\n{(result.stderr or '')[-2000:]}")

    after = set(BIOCYPHER_OUT.iterdir()) if BIOCYPHER_OUT.exists() else set()
    new_dirs = after - before
    assert new_dirs, (
        f"No new directory appeared under {BIOCYPHER_OUT} after the build; "
        f"before={sorted(p.name for p in before)}, after={sorted(p.name for p in after)}")
    out_dir = max(new_dirs, key=lambda p: p.stat().st_mtime)
    return out_dir, combined


@pytest.fixture(scope="module")
def build_test_mode(tmp_path_factory):
    return _run_build(tmp_path_factory.mktemp("build_test"), test_mode=True)


@pytest.fixture(scope="module")
def build_full(tmp_path_factory):
    return _run_build(tmp_path_factory.mktemp("build_full"), test_mode=False)


def _assert_no_error_lines(combined: str) -> None:
    errors = [l for l in combined.splitlines() if l.strip().startswith("ERROR")]
    assert errors == [], f"Found {len(errors)} ERROR line(s):\n" + "\n".join(errors)


def _observed_values(out_dir: Path) -> dict[tuple[str, str], set[str]]:
    """Scan BioCypher CSVs -> {(label_or_edge_type, property): {values}}.

    Real on-disk layout (verified against an actual --test build output, NOT
    assumed): BioCypher writes tab-delimited `<Label>-header.csv` +
    `<Label>-part000.csv` pairs for both node labels (e.g. `Gene`) and edge
    types (e.g. `Gene_has_pfam`) side by side in the same directory -- the
    field delimiter is a TAB (`config/biocypher_config.yaml` sets
    `delimiter: '\\t'`), not a semicolon. Header columns are named
    `prop` or `prop:type` (e.g. `values:string[]`, `end:long`, `:START_ID`);
    the bare property name is everything before the first `:`. String-typed
    field values are wrapped in single quotes (`'value'`); array values are a
    single quoted field with `|`-separated members (`array_delimiter: "|"`),
    not one quoted token per member. Numeric fields are unquoted. An
    empty-string property is written as a quoted empty (`''`) -- distinct
    from a fully absent field (empty `raw`) but equally "no value"; both are
    excluded from the observed set after stripping quotes.
    """
    import csv
    observed: dict[tuple[str, str], set[str]] = {}
    for header_file in out_dir.rglob("*-header.csv"):
        label = header_file.name.replace("-header.csv", "")
        with open(header_file, newline="") as hf:
            cols = [c.split(":")[0] for c in next(csv.reader(hf, delimiter="\t"))]
        for part in header_file.parent.glob(f"{label}-part*.csv"):
            with open(part, newline="") as f:
                for row in csv.reader(f, delimiter="\t"):
                    for col, raw in zip(cols, row):
                        if not raw:
                            continue
                        for v in raw.split("|"):
                            v = v.strip("'")
                            if not v:
                                # BioCypher writes an empty string property as a
                                # quoted empty (''); after stripping quotes that
                                # is indistinguishable from "field genuinely
                                # empty" and must not count as an observed value.
                                continue
                            observed.setdefault((label, col), set()).add(v)
    return observed


def _check_vocabularies(out_dir: Path, both_directions: bool) -> None:
    declared = load_vocabularies()
    observed = _observed_values(out_dir)
    problems = []
    for entry in declared.values():
        if not entry.closed:
            continue
        seen = observed.get((entry.applies_to, entry.property), set())
        undeclared = seen - set(entry.values)
        if undeclared:
            problems.append(
                f"{entry.id}: emitted undeclared value(s) {sorted(undeclared)}; "
                f"declared {sorted(entry.values)}")
        if both_directions and entry.exhaustive and not entry.expected_empty:
            unseen = set(entry.values) - seen
            if unseen and seen:
                problems.append(
                    f"{entry.id}: declared but never emitted: {sorted(unseen)}")
    assert not problems, "Vocabulary drift:\n" + "\n".join(problems)


def test_build_no_errors_test_mode(build_test_mode):
    _assert_no_error_lines(build_test_mode[1])


def test_vocabularies_in_csvs_test_mode(build_test_mode):
    """observed <= declared only: 100 items/adapter cannot contain rare values."""
    _check_vocabularies(build_test_mode[0], both_directions=False)


@pytest.mark.slow
def test_build_no_errors_full(build_full):
    _assert_no_error_lines(build_full[1])


@pytest.mark.slow
def test_vocabularies_in_csvs_full(build_full):
    _check_vocabularies(build_full[0], both_directions=True)
