from scripts.summarize_single_sweep_exact_recovery import summarize


def _record(truth: int, offsets: tuple[int, int, int]) -> dict:
    return {
        "truth": {"position": truth},
        "estimates": {
            "largest_gap_midpoint": truth + offsets[0],
            "composite": truth + offsets[1],
            "flank_ld_100kb": truth + offsets[2],
        },
    }


def test_exact_recovery_gate_passes_precise_consensus() -> None:
    records = [
        _record(1_000_000 + index * 100_000, (-5_000, 5_000, 10_000))
        for index in range(12)
    ]
    result = summarize(records)
    assert result["consensus_replicates"] == 12
    assert result["robustness_gate"]["passed"] is True


def test_exact_recovery_gate_rejects_consistent_offset() -> None:
    records = [
        _record(1_000_000 + index * 100_000, (70_000, 75_000, 80_000))
        for index in range(12)
    ]
    result = summarize(records)
    assert result["consensus_replicates"] == 12
    assert result["robustness_gate"]["passed"] is False
