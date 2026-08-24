import numpy as np

from scripts.fit_human_empirical_background import (
    EmpiricalRegime,
    augment_empirical,
    estimate_bin_profile,
    sample_bin_lengths,
)
from scripts.fit_human_family_relationships import summarize_raw_lengths


def test_empirical_profile_is_nonnegative_median_difference():
    reference = [{f"length_bin_{i}": float(i) for i in range(9)}]
    observed = [{f"length_bin_{i}": float(i + 2) for i in range(9)}]
    assert estimate_bin_profile(observed, reference) == [2.0] * 9
    assert estimate_bin_profile(reference, observed) == [0.0] * 9


def test_bounded_bin_sampler_stays_inside_bin():
    values = sample_bin_lengths(100, 5.0, 10.0, 0.7, np.random.default_rng(4))
    assert np.all(values >= 5.0)
    assert np.all(values <= 10.0)


def test_empirical_augmentation_adds_profiled_segments():
    record = {"indID1": "a", "indID2": "b", "label": "siblings"}
    record.update(summarize_raw_lengths([("1", 50.0)], genome_cm=100.0))
    profile = [0.0, 20.0] + [0.0] * 7
    result = augment_empirical(
        [record],
        {("a", "b"): [50.0]},
        profile,
        EmpiricalRegime(1.0, 0.0),
        np.random.default_rng(5),
        100.0,
    )[0]
    assert result["segments_ge_5"] > record["segments_ge_5"]
    assert result["segments_ge_10"] == record["segments_ge_10"]
