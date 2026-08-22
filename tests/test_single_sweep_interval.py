import numpy as np

from scripts.select_single_sweep_interval import choose_interval, competition_score


def _summary(errors: list[int], selection_count: int) -> dict:
    records = []
    selection_indices = []
    evaluation_indices = []
    for index in range(500):
        is_selection = ((index * 17 + 3) % 5) in [0, 1]
        target = selection_indices if is_selection else evaluation_indices
        target.append(index)
        if len(selection_indices) >= selection_count and len(evaluation_indices) >= len(errors) - selection_count:
            break
    indices = selection_indices[:selection_count] + evaluation_indices[: len(errors) - selection_count]
    for error, index in zip(errors, indices, strict=True):
        truth = 1_000_000
        estimate = truth + error
        records.append(
            {
                "seed": 75_222_000 + index,
                "truth": {"position": truth},
                "estimates": {
                    "largest_gap_midpoint": estimate,
                    "composite": estimate,
                    "flank_ld_100kb": estimate,
                },
            }
        )
    return {"held_out_method_gate": {"selected": "median_all_three"}, "records": records}


def test_competition_score_penalizes_interval_length() -> None:
    errors = np.asarray([1_000, 40_000])
    assert np.allclose(competition_score(errors, 10_000), [1 / 1.1, 0])
    assert competition_score(errors, 50_000).mean() < 1.0


def test_consensus_guard_expands_selection_width() -> None:
    summaries = [_summary([1_000, 2_000, 3_000, 4_000], 2)]
    observed = {
        "largest_gap_midpoint": 1_000_000,
        "composite_100kb": 1_010_000,
        "flank_ld_100kb": 970_000,
        "raisd": 990_000,
    }
    result = choose_interval(summaries, observed, np.asarray([5_000, 25_000]))
    assert result["center"] == 1_000_000
    assert result["selection_optimal_half_width"] == 5_000
    assert result["consensus_guard_half_width"] == 32_000
    assert result["final_half_width"] == 32_000
    assert result["interval"] == [968_000, 1_032_000]
