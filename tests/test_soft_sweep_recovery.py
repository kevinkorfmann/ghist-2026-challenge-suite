from scripts.summarize_soft_sweep_recovery import summarize


def _record(index: int, offset: int) -> dict:
    truth = 2_000_000 + index * 1_000
    estimates = {
        "composite": truth + offset + 40_000,
        "soft_composite": truth + offset,
        "h12": truth + offset + 2_000,
        "h12_h2": truth + offset - 2_000,
        "flank_ld_100kb": truth + offset + 10_000,
        "largest_gap_midpoint": truth + offset + 80_000,
    }
    return {
        "seed": 75_226_000 + index,
        "truth": {"position": truth},
        "estimates": estimates,
    }


def test_soft_recovery_selects_precise_haplotype_localizer() -> None:
    records = [_record(index, (index % 5 - 2) * 2_000) for index in range(48)]
    result = summarize(records)
    assert result["selected_method"] in {"soft_composite", "median_soft_haplotype"}
    assert result["robustness_gate"]["passed"] is True
    assert result["untouched_evaluation"]["mean_competition_score"] >= 0.4


def test_soft_recovery_rejects_large_untouched_offset() -> None:
    records = [
        _record(index, 0 if index < 24 else 400_000)
        for index in range(48)
    ]
    result = summarize(records)
    assert result["robustness_gate"]["passed"] is False
