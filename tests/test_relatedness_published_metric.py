"""Hand-calculated published-metric examples, including incomplete predictions."""

import pytest

from scripts.score_relatedness_published_metric import evaluate_published_metric


def row(pair, label):
    return {"indID1": str(pair), "indID2": "other", "relation": label}


def test_truth_class_weights_and_false_positive_denominators():
    truth = [row(0, "a"), row(1, "a"), row(2, "a"), row(3, "b")]
    predictions = [row(0, "a"), row(1, "a"), row(2, "b"), row(3, "a")]
    score = evaluate_published_metric([
        {"cohort_id": "one", "truth": truth, "predictions": predictions},
    ], ["a", "b"])
    assert score["pooled"]["weighted_jaccard"] == pytest.approx(0.375)
    assert score["pooled"]["per_label"]["a"]["jaccard"] == 0.5
    assert score["pair_coverage_passed"]


def test_missing_abstained_and_extra_pairs_are_retained():
    truth = [row(0, "a"), row(1, "a"), row(2, "b")]
    predictions = [row(0, "a"), row(2, None), row(3, "a")]
    score = evaluate_published_metric([
        {"cohort_id": "one", "truth": truth, "predictions": predictions},
    ], ["a", "b", "absent"])
    assert score["pooled"]["weighted_jaccard"] == pytest.approx(2 / 9)
    assert score["cohorts"][0]["missing_pairs"] == 1
    assert score["cohorts"][0]["extra_pairs"] == 1
    assert score["cohorts"][0]["abstained_pairs"] == 1
    assert not score["pair_coverage_passed"]
    assert score["pooled"]["per_label"]["absent"]["weight"] == 0


def test_pooled_counts_keep_cohorts_separate_and_reweight_after_pooling():
    cohorts = [
        {"cohort_id": "one", "truth": [row(0, "a")], "predictions": [row(0, "a")]},
        {"cohort_id": "two", "truth": [row(0, "b"), row(1, "b"), row(2, "b")],
         "predictions": [row(0, "a"), row(1, "a"), row(2, "a")]},
    ]
    score = evaluate_published_metric(cohorts, ["a", "b"])
    assert score["pooled"]["weighted_jaccard"] == pytest.approx(1 / 16)
    assert score["pooled"]["truth_pairs"] == 4
    assert score["cohorts"][0]["weighted_jaccard"] == 1


def test_reversed_duplicate_pair_is_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        evaluate_published_metric([{
            "cohort_id": "one", "truth": [row(0, "a")],
            "predictions": [row(0, "a"), {"indID1": "other", "indID2": "0", "relation": "a"}],
        }], ["a"])
