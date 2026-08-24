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
