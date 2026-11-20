import numpy as np
import pytest

from scripts.validate_task12_reference_candidate import score
from scripts.train_task12_reference_candidate import validate_protocol


def test_candidate_score_preserves_pair_identity_and_calibration():
    records = [dict(indID1="a", indID2="b", label="siblings"),
               dict(indID1="c", indID2="d", label="avuncular")]
    predictions = [
        dict(indID1="c", indID2="d", relation="siblings", confidence=.8,
             probabilities={"siblings": .8, "avuncular": .2}),
        dict(indID1="a", indID2="b", relation="siblings", confidence=.8,
             probabilities={"siblings": .8, "avuncular": .2}),
    ]
    result = score(predictions, records)
    assert result["accuracy"] == .5
    assert result["meioses_accuracy"] == .5
    assert result["ece_10_equal_width_bins"] == pytest.approx(.3)
    assert result["multiclass_brier_sum"] == pytest.approx(.68)
    assert result["negative_log_likelihood"] == pytest.approx(-np.log(.16) / 2)


def test_protocol_rejects_source_drift(tmp_path):
    protocol = {"source_sha256": {"scripts/train_task12_reference_candidate.py": "wrong"}}
    with pytest.raises(ValueError, match="Locked Task12 source"):
        validate_protocol(protocol, tmp_path)


@pytest.mark.parametrize("probability", [float("nan"), -1.0, 0.4, 2.0])
def test_candidate_score_rejects_invalid_probability_distribution(probability):
    records = [dict(indID1="a", indID2="b", label="siblings")]
    predictions = [dict(indID1="a", indID2="b", relation="siblings", confidence=probability,
                        probabilities={"siblings": probability})]
    with pytest.raises(ValueError, match="Invalid relationship probabilities"):
        score(predictions, records)
