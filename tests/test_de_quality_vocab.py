from multiomics_kg.utils.controlled_vocab import load_vocabularies

DATA_CAVEATS = {
    "fold_change_inflated_at_low_counts",
    "no_biological_replicates",
    "conditions_not_normalized_to_each_other",
    "comparison_direction_unclear",
    "table_may_be_a_different_comparison",
}
FOLD_CHANGE_PROVENANCE = {
    "computed_by_kg",
    "sign_reversed_by_kg",
    "significance_rule_differs_from_paper",
}


def test_de_quality_vocabularies():
    v = load_vocabularies()
    for vid in ("Experiment.data_caveats", "Experiment.fold_change_provenance",
                "Experiment.fold_change_type",
                "Experiment.low_signal_basis", "Changes_expression_of.low_signal_reason"):
        assert vid in v, vid
    caveats = v["Experiment.data_caveats"]
    assert set(caveats.values) == DATA_CAVEATS
    assert "paralogs_share_one_array_probe" not in caveats.values
    assert caveats.value_type == "string_array"
    provenance = v["Experiment.fold_change_provenance"]
    assert set(provenance.values) == FOLD_CHANGE_PROVENANCE
    assert provenance.value_type == "string_array"
    assert not FOLD_CHANGE_PROVENANCE & DATA_CAVEATS
    assert set(v["Experiment.fold_change_type"].values) == {"log2", "linear", "signed_linear"}
    basis = {"condition_at_column_minimum", "mean_in_bottom_decile"}
    assert set(v["Experiment.low_signal_basis"].values) == basis
    assert set(v["Changes_expression_of.low_signal_reason"].values) == basis
    for vid in ("Experiment.data_caveats", "Experiment.fold_change_provenance",
                "Experiment.low_signal_basis",
                "Changes_expression_of.low_signal_reason"):
        descs = v[vid].value_descriptions
        assert len(descs) == len(v[vid].values), vid
