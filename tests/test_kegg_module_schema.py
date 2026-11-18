"""KeggModule wiring: schema, controlled vocabularies, post-import (.sh and .cypher), build script."""
from pathlib import Path

import yaml

from multiomics_kg.utils.controlled_vocab import load_vocabularies

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = yaml.safe_load((ROOT / "config" / "schema_config.yaml").read_text())
EDGES = {"kegg_term_in_kegg_module": ("kegg term", "kegg module"),
         "kegg_module_is_a_kegg_module": ("kegg module", "kegg module"),
         "kegg_module_in_kegg_pathway": ("kegg module", "kegg term"),
         "organism_has_kegg_module": ("organism taxon", "kegg module")}


def test_schema_declares_kegg_module_node():
    node = SCHEMA["kegg module"]
    assert node["represented_as"] == "node" and node["preferred_id"] == "kegg.module"
    props = node["properties"]
    for p in ("name", "level", "level_kind", "module_type", "module_id", "definition", "step_count",
              "gene_count", "direct_gene_count", "organism_count", "complete_organism_count"):
        assert p in props, p


def test_schema_declares_the_four_edges():
    by_label = {v.get("label_in_input"): v for v in SCHEMA.values() if isinstance(v, dict)}
    for label, (src, tgt) in EDGES.items():
        e = by_label[label]
        assert (e["label_as_edge"], e["source"], e["target"]) == (label, src, tgt)
    ko = by_label["kegg_term_in_kegg_module"]["properties"]
    assert ko == {"steps": "int[]", "component": "str"}
    org = by_label["organism_has_kegg_module"]["properties"]
    assert set(org) == {"status", "steps_total", "steps_present", "steps_present_via_ec",
                        "missing_step_indices", "missing_step_definitions", "missing_kos",
                        "ec_filled_step_indices", "ec_genes"}


def test_post_import_declarations_carry_a_comment():
    text = (ROOT / "config" / "schema_config.yaml").read_text()
    block = text.split("\nkegg module:\n", 1)[1].split("\n\n", 1)[0]
    for p in ("gene_count", "direct_gene_count", "organism_count", "complete_organism_count"):
        line = next(ln for ln in block.splitlines() if ln.strip().startswith(p + ":"))
        assert "# post-import" in line, line


def test_vocabularies():
    entries = {f"{e.applies_to}.{e.property}": e
               for e in load_vocabularies(ROOT / "config" / "controlled_vocabularies.yaml").values()}
    assert entries["KeggModule.module_type"].values == ["pathway", "signature"]
    assert entries["KeggModule.level_kind"].values == [
        "kegg_module_type", "kegg_module_class", "kegg_module_subclass", "kegg_module"]
    assert entries["Kegg_term_in_kegg_module.component"].values == ["required", "optional"]
    status = entries["Organism_has_kegg_module.status"]
    assert status.values == ["complete", "one_step_missing", "two_steps_missing", "partial"]
    assert status.value_descriptions["complete"] == "every step has an annotated gene; not shown to run"
    assert entries["Organism_has_kegg_module.steps_total"].min_value == 1
    for p in ("steps_present", "steps_present_via_ec"):
        assert entries[f"Organism_has_kegg_module.{p}"].min_value == 0


def test_both_post_import_files_carry_the_kegg_module_block():
    for name in ("post-import.sh", "post-import.cypher"):
        text = (ROOT / "scripts" / name).read_text()
        for needle in ("kegg_module_id_idx", "kegg_module_level_idx", "kegg_module_level_kind_idx",
                       "kegg_module_type_idx", "kegg_module_module_id_idx", "keggModuleFullText",
                       "complete_organism_count", "Kegg_term_in_kegg_module", "kegg_release"):
            assert needle in text, f"{needle} missing from {name}"


def test_build_wires_the_adapter_and_release_file():
    text = (ROOT / "create_knowledge_graph.py").read_text()
    assert "MultiKeggModuleAdapter(" in text
    assert "all_pathway_ids()" in text
    assert "kegg_release.txt" in text


def test_post_import_rematerializes_empty_missing_lists():
    for name in ("post-import.sh", "post-import.cypher"):
        text = (ROOT / "scripts" / name).read_text()
        for prop in ("missing_step_indices", "missing_step_definitions", "missing_kos"):
            assert f"r.{prop} = coalesce(r.{prop}, [])" in text, f"{prop} in {name}"
