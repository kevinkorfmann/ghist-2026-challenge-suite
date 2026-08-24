import pytest

from scripts.evaluate_human_background_robust import FINAL_REGIMES, summarize_runs


def test_final_background_regimes_are_fixed_and_distinct():
    assert len(FINAL_REGIMES) == 5
    assert len({regime.mean_segments for regime in FINAL_REGIMES}) == 5
    assert min(regime.mean_segments for regime in FINAL_REGIMES) > 0


def test_holdout_summary_reports_mean_and_minimum():
    summary = summarize_runs(
        [
            {"exact_accuracy": 0.4, "meioses_accuracy": 0.6},
            {"exact_accuracy": 0.5, "meioses_accuracy": 0.7},
        ]
    )
    assert summary["mean_exact_accuracy"] == 0.45
    assert summary["minimum_exact_accuracy"] == 0.4
    assert summary["mean_meioses_accuracy"] == pytest.approx(0.65)
