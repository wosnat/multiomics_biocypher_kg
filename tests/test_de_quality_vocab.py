from multiomics_kg.utils.controlled_vocab import load_vocabularies

DATA_CAVEATS = {
    "fold_change_inflated_at_low_counts",
    "no_biological_replicates",
    "conditions_not_normalized_to_each_other",
    "comparison_direction_unclear",
    "table_may_be_a_different_comparison",
    "significance_calls_differ_from_paper",
    "fold_change_computed_by_kg",
    "fold_change_sign_reversed_by_kg",
}


def test_de_quality_vocabularies():
    v = load_vocabularies()
    for vid in ("Experiment.data_caveats", "Experiment.fold_change_type",
                "Experiment.low_signal_basis", "Changes_expression_of.low_signal_reason"):
        assert vid in v, vid
    caveats = v["Experiment.data_caveats"]
    assert set(caveats.values) == DATA_CAVEATS
    assert "paralogs_share_one_array_probe" not in caveats.values
    assert caveats.value_type == "string_array"
    assert set(v["Experiment.fold_change_type"].values) == {"log2", "linear", "signed_linear"}
    basis = {"condition_at_column_minimum", "mean_in_bottom_decile"}
    assert set(v["Experiment.low_signal_basis"].values) == basis
    assert set(v["Changes_expression_of.low_signal_reason"].values) == basis
    for vid in ("Experiment.data_caveats", "Experiment.low_signal_basis",
                "Changes_expression_of.low_signal_reason"):
        descs = v[vid].value_descriptions
        assert len(descs) == len(v[vid].values), vid
