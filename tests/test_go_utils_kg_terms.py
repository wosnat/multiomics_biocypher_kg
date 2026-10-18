"""kg_go_terms: the single definition of the KG's GO node set."""
import json

from multiomics_kg.utils.go_utils import kg_go_terms

GO_DATA = {
    "GO:1": {"name": "a", "namespace": "molecular_function", "parents": [["GO:0", "is_a"]]},
    "GO:0": {"name": "root", "namespace": "molecular_function", "parents": []},
    "GO:9": {"name": "b", "namespace": "biological_process", "parents": [["GO:1", "regulates"]]},
}


def _dir(tmp_path, name, terms):
    d = tmp_path / name
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(
        json.dumps({"g1": {"locus_tag": "g1", "go_terms": terms}}), encoding="utf-8"
    )
    return d


def test_kg_go_terms_seeds_and_ancestors(tmp_path):
    dirs = [_dir(tmp_path, "a", ["GO:1"]), _dir(tmp_path, "b", ["GO:9", "GO:404"])]
    assert kg_go_terms(GO_DATA, dirs) == {
        "GO:1": "molecular function",
        "GO:0": "molecular function",
        "GO:9": "biological process",
    }


def test_kg_go_terms_resolves_alt_ids(tmp_path):
    dirs = [_dir(tmp_path, "a", ["GO:ALT"])]
    assert kg_go_terms(GO_DATA, dirs) == {}
    got = kg_go_terms(GO_DATA, dirs, alt_ids={"GO:ALT": "GO:1"})
    assert got == {"GO:1": "molecular function", "GO:0": "molecular function"}
