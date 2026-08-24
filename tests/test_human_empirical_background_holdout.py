from scripts.fit_human_empirical_background import FINAL_REGIMES


def test_empirical_final_regimes_are_fixed_and_observed_centered():
    assert len(FINAL_REGIMES) == 5
    assert [regime.intensity for regime in FINAL_REGIMES] == [
        0.55,
        0.80,
        1.00,
        1.20,
        1.45,
    ]
