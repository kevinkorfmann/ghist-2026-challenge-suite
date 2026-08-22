from scripts.summarize_single_sweep_exact_recovery import summarize


def _record(seed: int, truth: int, offsets: tuple[int, int, int]) -> dict:
    return {
        "seed": seed,
        "truth": {"position": truth},
        "estimates": {
            "largest_gap_midpoint": truth + offsets[0],
            "composite": truth + offsets[1],
            "flank_ld_100kb": truth + offsets[2],
        },
    }


def test_exact_recovery_gate_passes_precise_consensus() -> None:
    records = [
        _record(75_225_000 + index, 1_000_000 + index * 100_000, (-5_000, 5_000, 10_000))
        for index in range(128)
    ]
    result = summarize(records)
    assert result["calibration_consensus_replicates"] == 64
    assert result["evaluation_consensus_replicates"] == 64
    assert result["robustness_gate"]["passed"] is True


def test_exact_recovery_gate_rejects_consistent_offset() -> None:
    records = [
        _record(75_225_000 + index, 1_000_000 + index * 100_000, (150_000, 155_000, 160_000))
        for index in range(128)
    ]
    result = summarize(records)
    assert result["robustness_gate"]["passed"] is False
