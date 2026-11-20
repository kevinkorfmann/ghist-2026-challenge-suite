import json
from pathlib import Path

import numpy as np
import pytest

from ghist_single_sweep import recovery
from ghist_single_sweep.scan import VariantMatrix
from scripts.freeze_task7_sweep_bgs import validate_matched_recovery
from scripts.summarize_single_sweep_selection_recovery import load_reference_estimates


def test_localizer_uses_scan_peak_without_truth(monkeypatch):
    observed = {}

    def ld(variants, candidates, flank_width):
        observed["candidates"] = candidates
        return {
            "candidate_position": candidates,
            "flank_ld_ratio": np.arange(len(candidates)),
        }

    monkeypatch.setattr(recovery, "flank_ld_scan", ld)
    variants = VariantMatrix(
        "sim",
        4_000_000,
        np.asarray([2_400_000, 2_410_000, 2_450_000]),
        np.zeros((3, 4)),
        0,
    )
    result = recovery.localize_from_scan(
        variants,
        {
            "center": np.asarray([500_000, 2_500_000]),
            "composite_score": np.asarray([0, 1]),
        },
        radius=300_000,
    )
    assert result["localization"]["region"] == [2_200_000, 2_800_000]
    assert result["localization"]["truth_used_for_inference"] is False
    assert observed["candidates"].min() == 2_200_000
    assert observed["candidates"].max() < 2_800_000


def test_legacy_oracle_recovery_and_duplicate_seeds_rejected():
    with pytest.raises(ValueError, match="truth-blind"):
        recovery.validate_blind_records([{"seed": 1}])
    record = {
        "seed": 1,
        "localization": {
            "truth_used_for_inference": False,
            "region_source": "full_scan_composite_maximum",
        },
    }
    with pytest.raises(ValueError, match="distinct"):
        recovery.validate_blind_records([record, record])


def test_task7_cannot_transfer_non_bgs_recovery():
    payload = {
        "inference_method": "median_multiscale_raisd_gap_ld",
        "records": [
            {
                "seed": 1,
                "localization": {
                    "truth_used_for_inference": False,
                    "region_source": "full_scan_composite_maximum",
                },
                "background_selection": False,
                "demography": "Zigzag_1S14",
            }
        ],
    }
    with pytest.raises(ValueError, match="BGS"):
        validate_matched_recovery(payload)
    # Task 7 itself does not mandate Task 8's bottleneck/growth history.
    payload["records"][0]["background_selection"] = True
    payload["records"][0]["demography"] = "constant_10000"
    validate_matched_recovery(payload)
    payload["records"][0].pop("demography")
    with pytest.raises(ValueError, match="demographic history"):
        validate_matched_recovery(payload)


def test_raisd_reference_matches_filename_not_lexical_position(tmp_path):
    paths = [Path("task_10.json"), Path("task_2.json")]
    records = [{"truth": {"position": 100}}, {"truth": {"position": 200}}]
    for path, truth, estimate in zip(paths, (100, 200), (110, 220), strict=True):
        (tmp_path / path).write_text(
            json.dumps({"truth_position": truth, "maximum": {"position": estimate}})
        )
    np.testing.assert_array_equal(
        load_reference_estimates(paths, records, tmp_path), [110, 220]
    )
    with pytest.raises(ValueError, match="different simulation"):
        load_reference_estimates(paths, records[::-1], tmp_path)


def test_region_metrics_respect_half_open_bed_coordinates():
    from scripts.score_sweep_background_validation import region_metrics

    result = region_metrics([[10, 20], [30, 40]], [10, 20])
    assert result["covered_sweeps"] == 1
    assert result["true_regions"] == 1
    assert result["false_regions"] == 1
    assert result["covered_bases"] == 20
