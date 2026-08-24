import numpy as np

from scripts.fit_task11_degree_factorized import (
    degree_guarded_probabilities,
    rich_feature_names,
    rich_summary,
)
from scripts.fit_task11_dog_relationships import feature_names


def test_rich_summary_adds_ibd_state_and_chromosome_features():
    records = [
        ("1", 12.0, False),
        ("1", 8.0, True),
        ("2", 3.0, False),
    ]
    values = rich_summary(records)
    result = dict(zip(rich_feature_names(), values, strict=True))
    assert len(values) > len(feature_names())
    assert result["ibd1_segments"] == 2
    assert result["ibd2_segments"] == 1
    assert result["chromosomes_with_ibd2"] == 1
    assert np.isclose(result["chromosome_sharing_max"], 28.0)


def test_rich_summary_empty_pair_is_finite():
    values = rich_summary([])
    assert len(values) == len(rich_feature_names())
    assert np.all(np.isfinite(values))


def test_degree_guard_preserves_baseline_degree():
    candidate = np.asarray([[0.05, 0.05, 0.01, 0.01, 0.40, 0.02, 0.45, 0.01]])
    baseline = np.asarray([[0.10, 0.55, 0.02, 0.02, 0.10, 0.01, 0.19, 0.01]])
    guarded = degree_guarded_probabilities(candidate, baseline)
    assert np.argmax(guarded, axis=1).item() in (0, 1)
    assert np.isclose(guarded.sum(), 1.0)
