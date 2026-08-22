import pytest

from scripts.freeze_single_sweep import validate


def _recovery(q90: float) -> dict:
    metrics = {"quantiles": {"0.9": q90}, "median_absolute_error": 5_000}
    return {
        "held_out_method_gate": {
            "selected": "median_all_three",
            "evaluation_metrics": {"median_all_three": metrics},
        }
    }


def _candidate() -> dict:
    return {
        "interval_selection": {"interval": [965_000, 1_035_000], "final_half_width": 35_000},
        "observed_estimates": {"gap": 1_000_000, "composite": 1_020_000},
    }


def _raisd() -> dict:
    return {
        "maximum": {"position": 990_000},
        "separated_peaks_100kb": [{"mu": 120.0}, {"mu": 30.0}],
    }


def test_freeze_gate_accepts_precise_robust_recovery() -> None:
    result = validate(
        _recovery(60_000),
        _recovery(70_000),
        {"robustness_gate": {"passed": True}},
        _candidate(),
        _raisd(),
    )
    assert result["candidate_hit_score"] == pytest.approx(1 / 1.35)


def test_freeze_gate_rejects_scale_failure() -> None:
    with pytest.raises(ValueError, match="scale stress"):
        validate(
            _recovery(60_000),
            _recovery(70_000),
            {"robustness_gate": {"passed": False}},
            _candidate(),
            _raisd(),
        )
