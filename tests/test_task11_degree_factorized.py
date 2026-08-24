import numpy as np

from scripts.fit_task11_degree_factorized import rich_feature_names, rich_summary
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
