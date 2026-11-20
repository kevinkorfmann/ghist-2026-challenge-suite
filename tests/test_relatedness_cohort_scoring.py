from copy import deepcopy

import pytest

from scripts.score_relatedness_cohorts import evaluate_cohorts


LABELS = ["a", "b"]


def row(index, label, probabilities=None):
    result = {"indID1": f"left_{index}", "indID2": f"right_{index}", "relation": label}
    if probabilities is not None:
        result["probabilities"] = probabilities
    return result


def cohort(name="one", correct=True, repeat=1):
    truth = [row(i, LABELS[i % 2]) for i in range(2 * repeat)]
    predictions = [row(i, LABELS[i % 2 if correct else 1 - i % 2]) for i in range(2 * repeat)]
    return {"cohort_id": name, "truth": truth, "predictions": predictions}


def test_missing_and_abstained_pairs_stay_in_denominators():
    item = cohort()
    item["truth"].append(row(2, "a"))
    item["predictions"] = [row(0, "a"), row(2, None)]
    report = evaluate_cohorts([item], LABELS)
    assert report["pooled"]["pairs"] == 3
    assert report["pooled"]["exact_accuracy"] == pytest.approx(1 / 3)
    assert report["pooled"]["balanced_accuracy"] == 0.25
    assert report["missing_pairs"] == report["abstained_pairs"] == 1
    assert not report["pair_coverage_passed"]
    assert report["calibration"] is None


def test_extra_cross_family_pairs_fail_coverage_without_inflating_denominator():
    item = cohort()
    item["predictions"].append(row(999, "a"))
    report = evaluate_cohorts([item], LABELS)
    assert report["pooled"]["pairs"] == 2
    assert report["pooled"]["exact_accuracy"] == 1
    assert report["extra_pairs"] == 1
    assert not report["pair_coverage_passed"]
    assert not report["release_authorized"]


@pytest.mark.parametrize("section", ["truth", "predictions"])
def test_reversed_duplicates_are_rejected(section):
    item = cohort()
    duplicate = deepcopy(item[section][0])
    duplicate["indID1"], duplicate["indID2"] = duplicate["indID2"], duplicate["indID1"]
    item[section].append(duplicate)
    with pytest.raises(ValueError, match="Duplicate unordered"):
        evaluate_cohorts([item], LABELS)


def test_pair_order_is_irrelevant():
    item = cohort()
    original = evaluate_cohorts([item], LABELS)
    item["predictions"].reverse()
    for prediction in item["predictions"]:
        prediction["indID1"], prediction["indID2"] = prediction["indID2"], prediction["indID1"]
    assert evaluate_cohorts([item], LABELS) == original


def test_cohort_bootstrap_does_not_treat_repeated_pairs_as_new_units():
    small = evaluate_cohorts([cohort("good"), cohort("bad", False)], LABELS)
    large = evaluate_cohorts([cohort("good", repeat=100), cohort("bad", False, 100)], LABELS)
    assert small["bootstrap"] == large["bootstrap"]
    assert small["bootstrap"]["exact_accuracy_interval"] == [0, 1]
    assert small["bootstrap"]["balanced_accuracy_interval"] == [0, 1]


def test_no_intervals_for_one_cohort():
    report = evaluate_cohorts([cohort()], LABELS)
    assert report["bootstrap"]["status"] == "insufficient_cohorts"
    assert report["bootstrap"]["balanced_accuracy_interval"] is None


def test_missing_class_bootstrap_is_not_conditioned_on_complete_resamples():
    first, second = cohort("a"), cohort("b")
    for section in ("truth", "predictions"):
        first[section] = first[section][:1]
        second[section] = second[section][1:]
    report = evaluate_cohorts([first, second], LABELS)
    assert report["pooled"]["all_labels_represented"]
    assert report["bootstrap"]["replicates_missing_labels"] > 0
    assert report["bootstrap"]["status"] == "missing_label_support"
    assert report["bootstrap"]["balanced_accuracy_interval"] is None


def test_unrepresented_labels_are_explicit():
    report = evaluate_cohorts([cohort()], LABELS + ["c"])
    assert report["pooled"]["balanced_accuracy"] is None
    assert report["pooled"]["observed_label_macro_recall"] == 1
    assert report["pooled"]["per_label"]["c"] == {"pairs": 0, "correct": 0, "recall": None}


def test_calibration_uses_probability_of_emitted_label_not_maximum():
    item = {"cohort_id": "one", "truth": [row(0, "a")],
            "predictions": [row(0, "b", {"a": 0.9, "b": 0.1})]}
    report = evaluate_cohorts([item], LABELS)
    assert report["calibration"]["ece"] == pytest.approx(0.1)
    assert report["calibration"]["bins"][1]["pairs"] == 1


def test_calibration_includes_confidence_one_and_requires_all_probabilities():
    item = cohort()
    item["predictions"][0]["probabilities"] = {"a": 1.0, "b": 0.0}
    assert not evaluate_cohorts([item], LABELS)["calibration_available"]
    item["predictions"][1]["probabilities"] = {"a": 0.0, "b": 1.0}
    report = evaluate_cohorts([item], LABELS)
    assert report["calibration"]["ece"] == 0
    assert report["calibration"]["bins"][9]["pairs"] == 2


@pytest.mark.parametrize("probabilities", [
    {"a": 1}, {"a": 1, "b": 1}, {"a": float("nan"), "b": 0},
    {"a": -0.1, "b": 1.1}, {"a": True, "b": False}, {"a": "1", "b": 0},
])
def test_invalid_probability_records_fail(probabilities):
    item = cohort()
    item["predictions"][0]["probabilities"] = probabilities
    with pytest.raises(ValueError, match="[Pp]robabilit"):
        evaluate_cohorts([item], LABELS)


def test_empty_cohort_or_reused_identifier_is_not_independence():
    with pytest.raises(ValueError, match="unique"):
        evaluate_cohorts([cohort(), cohort()], LABELS)
    with pytest.raises(ValueError, match="truth pairs"):
        evaluate_cohorts([{"cohort_id": "empty", "truth": [], "predictions": []}], LABELS)
