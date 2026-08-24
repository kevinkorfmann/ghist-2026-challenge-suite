from scripts.evaluate_human_background_corrected_holdout import aggregate


def test_corrected_holdout_aggregate():
    summary = aggregate(
        [
            {"exact_accuracy": 0.4, "meioses_accuracy": 0.5},
            {"exact_accuracy": 0.6, "meioses_accuracy": 0.7},
        ]
    )
    assert summary["mean_exact_accuracy"] == 0.5
    assert summary["minimum_exact_accuracy"] == 0.4
    assert summary["mean_meioses_accuracy"] == 0.6
