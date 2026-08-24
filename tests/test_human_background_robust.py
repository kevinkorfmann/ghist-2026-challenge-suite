import numpy as np

from scripts.fit_human_background_robust import BackgroundRegime, augment_records
from scripts.fit_human_family_relationships import summarize_raw_lengths


def test_background_augmentation_is_reproducible_and_additive():
    record = {"indID1": "a", "indID2": "b", "label": "siblings"}
    record.update(summarize_raw_lengths([("1", 50.0)], genome_cm=100.0))
    regime = BackgroundRegime(
        mean_segments=3.0,
        excess_scale_cm=10.0,
        pair_heterogeneity=0.0,
    )
    first = augment_records(
        [record], {("a", "b"): [50.0]}, regime, np.random.default_rng(7), 100.0
    )
    second = augment_records(
        [record], {("a", "b"): [50.0]}, regime, np.random.default_rng(7), 100.0
    )
    assert first == second
    assert first[0]["segments_ge_20"] >= record["segments_ge_20"]
    assert first[0]["total_cm_ge_20"] >= record["total_cm_ge_20"]


def test_zero_background_preserves_summary():
    record = {"indID1": "a", "indID2": "b", "label": "siblings"}
    record.update(summarize_raw_lengths([("1", 12.0), ("2", 45.0)], genome_cm=100.0))
    output = augment_records(
        [record],
        {("a", "b"): [12.0, 45.0]},
        BackgroundRegime(0.0, 10.0, 0.0),
        np.random.default_rng(8),
        100.0,
    )[0]
    for key, value in record.items():
        assert output[key] == value
