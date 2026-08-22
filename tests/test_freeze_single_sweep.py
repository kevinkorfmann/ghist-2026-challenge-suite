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


def _candidate(half_width: int = 35_000) -> dict:
    return {
        "interval_selection": {
            "interval": [1_000_000 - half_width, 1_000_000 + half_width],
            "final_half_width": half_width,
        },
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
    with pytest.raises(ValueError, match="exact-scale recovery"):
        validate(
            _recovery(60_000),
            _recovery(70_000),
            {"robustness_gate": {"passed": False}},
            _candidate(),
            _raisd(),
        )


def test_freeze_gate_accepts_calibrated_60kb_half_width() -> None:
    result = validate(
        _recovery(60_000),
        _recovery(70_000),
        {"robustness_gate": {"passed": True}},
        _candidate(60_000),
        _raisd(),
    )
    assert result["candidate_hit_score"] == pytest.approx(1 / 1.6)


def test_freeze_gate_rejects_half_width_above_75kb() -> None:
    with pytest.raises(ValueError, match="precision gate"):
        validate(
            _recovery(60_000),
            _recovery(70_000),
            {"robustness_gate": {"passed": True}},
            _candidate(76_000),
            _raisd(),
        )
