from scripts.summarize_single_sweep_scale_stress import summarize


def _record(scale: int, truth: int, gap: int, composite: int, linkage: int) -> dict:
    return {
        "truth": {"position": truth},
        "estimates": {
            "largest_gap_midpoint": gap,
            "composite": composite,
            "flank_ld_100kb": linkage,
        },
        "slim_scaling_factor": scale,
    }


def test_scale_stress_uses_locked_median_and_passes_good_cells() -> None:
    records = []
    for scale in [1, 2, 5, 10]:
        records.extend(
            [
                _record(scale, 1_000_000, 995_000, 1_010_000, 1_020_000),
                _record(scale, 2_000_000, 2_004_000, 1_980_000, 2_030_000),
            ]
        )
    result = summarize(records)
    assert result["selected_method_locked_before_stress_test"] == "median_all_three"
    assert result["robustness_gate"]["passed"] is True
    assert result["replicates"] == 8


def test_scale_stress_gate_rejects_large_scale_specific_error() -> None:
    records = [
        _record(1, 1_000_000, 1_200_000, 1_200_000, 1_200_000),
        _record(1, 2_000_000, 2_200_000, 2_200_000, 2_200_000),
        _record(10, 1_000_000, 1_001_000, 1_002_000, 1_003_000),
        _record(10, 2_000_000, 2_001_000, 2_002_000, 2_003_000),
    ]
    result = summarize(records)
    assert result["robustness_gate"]["passed"] is False
