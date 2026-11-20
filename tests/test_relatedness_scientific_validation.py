import pytest

from ghist_relatedness.validation import validate_simulation_splits
from scripts.validate_relatedness_scientifically import probability_metrics, pedigree


def test_copied_simulations_cannot_cross_splits(tmp_path):
    training = tmp_path / "training.seg"
    renamed = tmp_path / "evaluation.seg"
    training.write_text("one generated pedigree\n")
    renamed.write_bytes(training.read_bytes())
    with pytest.raises(ValueError, match="Duplicate simulation content"):
        validate_simulation_splits({"training": [training], "evaluation": [renamed]})


def test_independent_simulation_manifest(tmp_path):
    training = tmp_path / "training.seg"
    evaluation = tmp_path / "evaluation.seg"
    training.write_text("pedigree seed 1\n")
    evaluation.write_text("pedigree seed 2\n")
    result = validate_simulation_splits(
        {"training": [training], "evaluation": [evaluation]}
    )
    assert result["training"][0]["sha256"] != result["evaluation"][0]["sha256"]


def test_probability_metrics_groups_correlated_pedigree_pairs():
    report = probability_metrics(
        ["a", "b"],
        [[0.8, 0.2], [0.8, 0.2], [0.1, 0.9]],
        ["a", "a", "b"],
        ["ped1", "ped1", "ped2"],
        {"a": 2, "b": 3},
    )
    assert report["pairs"] == 3
    assert report["independent_pedigrees"] == 2
    assert report["accuracy"] == 1
    assert report["multiclass_brier_sum"] == pytest.approx(0.06)
    assert report["ece_10_equal_width_bins"] == pytest.approx(1 / 6)
    assert report["pedigree_bootstrap_accuracy_95pct"] == [1.0, 1.0]


def test_probability_metrics_rejects_bad_probability():
    with pytest.raises(ValueError, match="sum to one"):
        probability_metrics(["a", "b"], [[0.8, 0.4]], ["a"], ["ped1"], {"a": 2, "b": 3})


def test_pedigree_group_preserves_grandparent_name():
    assert pedigree("grandparent_grandchild17_g3-b1-i1") == "grandparent_grandchild17"


def test_background_evaluation_uses_independent_reference(monkeypatch, tmp_path):
    from scripts import fit_human_empirical_background as empirical

    evaluation = [
        {"label": "siblings", "indID1": "evaluation1", "indID2": "evaluation2"}
    ]
    clean_reference = [
        {"label": "siblings", "indID1": "reference1", "indID2": "reference2"}
    ]
    seen = []
    monkeypatch.setattr(
        empirical, "read_pedsim_with_lengths", lambda path, genome: (evaluation, {})
    )
    monkeypatch.setattr(empirical, "augment_empirical", lambda *args: evaluation)

    def correct(records, reference, minimum, genome):
        seen.append(reference)
        return records

    monkeypatch.setattr(empirical, "correct_records", correct)
    empirical.regime_datasets(
        tmp_path / "evaluation.seg",
        [],
        (empirical.EmpiricalRegime(1.0, 0.4),),
        123,
        5,
        3300,
        clean_reference,
    )
    assert seen == [clean_reference]
    assert seen[0] is not evaluation


@pytest.mark.parametrize(
    "module_name,function_name,augment_name",
    [
        ("scripts.fit_human_empirical_background", "training_set", "augment_empirical"),
        (
            "scripts.fit_human_background_corrected",
            "corrected_training_set",
            "augment_records",
        ),
    ],
)
def test_background_training_reference_excludes_own_replicate(
    monkeypatch, tmp_path, module_name, function_name, augment_name
):
    from importlib import import_module

    empirical = import_module(module_name)

    inputs = {
        tmp_path / "a.seg": [{"label": "siblings", "source": "a"}],
        tmp_path / "b.seg": [{"label": "siblings", "source": "b"}],
    }
    seen = []
    monkeypatch.setattr(
        empirical, "read_pedsim_with_lengths", lambda path, genome: (inputs[path], {})
    )
    monkeypatch.setattr(empirical, augment_name, lambda records, *args: records)

    def correct(records, reference, minimum, genome):
        seen.append((records[0]["source"], {r["source"] for r in reference}))
        return records

    monkeypatch.setattr(empirical, "correct_records", correct)
    args = (
        (list(inputs), [], 1, 123, 5, 3300)
        if function_name == "training_set"
        else (list(inputs), 1, 123, 5, 3300)
    )
    getattr(empirical, function_name)(*args)
    assert seen == [("a", {"b"}), ("b", {"a"})]


def test_background_holdout_blend_uses_saved_weight(monkeypatch):
    import numpy as np
    from scripts import fit_human_background_corrected as corrected

    class Model:
        classes_ = np.asarray(["siblings", "avuncular"])

        def predict_proba(self, values):
            return np.asarray([[0.9, 0.1]])

    monkeypatch.setattr(corrected, "matrix", lambda *args: np.ones((1, 1)))
    monkeypatch.setattr(
        corrected,
        "predict_probabilities",
        lambda *args: (np.asarray(["siblings", "avuncular"]), np.asarray([[0.1, 0.9]])),
    )
    labels, values = corrected.blended_probabilities({}, Model(), [], [{}], [], 0.25)
    assert values[0, list(labels).index("siblings")] == pytest.approx(0.7)
    assert values[0, list(labels).index("avuncular")] == pytest.approx(0.3)
