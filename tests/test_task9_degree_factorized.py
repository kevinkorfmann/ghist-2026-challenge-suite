import numpy as np

from scripts.fit_task9_degree_factorized import LABELS, aligned, metrics


def test_task9_probability_alignment_and_metrics():
    probabilities = aligned(
        np.asarray(["siblings", "avuncular"]),
        np.asarray([[0.8, 0.2], [0.1, 0.9]]),
    )
    assert probabilities.shape == (2, len(LABELS))
    result = metrics(probabilities, np.asarray(["siblings", "avuncular"]))
    assert result["accuracy"] == 1.0
    assert result["meioses_accuracy"] == 1.0


def test_legacy_probabilities_are_restricted_to_current_labels():
    probabilities = aligned(np.array(['half_sibling', 'siblings']), np.array([[0.8, 0.2]]))
    assert probabilities.sum() == 1.0
    assert probabilities[0, LABELS.index('siblings')] == 1.0
    assert 'half_sibling' not in LABELS
    assert 'half_2nd_cousin' not in LABELS


def test_factorized_model_handles_single_observed_class_per_degree():
    from types import SimpleNamespace
    from scripts.fit_task9_degree_factorized import factorized_probabilities
    model = SimpleNamespace(classes_=np.array(['full_2nd_cousin']))
    degree = SimpleNamespace(classes_=np.array([6]), predict_proba=lambda x: np.ones((len(x), 1)))
    result = factorized_probabilities({'exact': model, 'degree': degree, 'conditional': {}}, np.zeros((2, 1)))
    assert np.all(result[:, LABELS.index('full_2nd_cousin')] == 1.0)


def test_frozen_model_adapter_drops_retired_labels_without_retraining():
    import pytest
    from scripts.predict_task9_current import accepted_probabilities
    classes, probabilities = accepted_probabilities(['half_sibling', 'siblings', 'avuncular'], [[0.6, 0.3, 0.1]])
    assert classes.tolist() == ['siblings', 'avuncular']
    np.testing.assert_allclose(probabilities, [[0.75, 0.25]])
    with pytest.raises(ValueError):
        accepted_probabilities(['half_sibling', 'siblings'], [[1.0, 0.0]])
