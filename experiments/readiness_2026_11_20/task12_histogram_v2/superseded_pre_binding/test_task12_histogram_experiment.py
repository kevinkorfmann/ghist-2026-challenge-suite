import json

import pytest

from scripts.fit_human_family_relationships import LABELS
from scripts.task12_histogram_protocol import (
    score_predictions,
    sha256,
    validate_pairs,
    validate_protocol,
    write_json,
)
from scripts.train_task12_histogram_candidate import build_training_records
from scripts.validate_task12_histogram_candidate import validate_candidate_lock


def score_fixture():
    rows, predictions = [], []
    for label in LABELS:
        pair = {"indID1": f"{label}1_g1", "indID2": f"{label}1_g2"}
        rows.append({**pair, "label": label})
        predictions.append(
            {
                **pair,
                "relation": label,
                "confidence": 1.0,
                "probabilities": {other: float(other == label) for other in LABELS},
            }
        )
    return rows, predictions


def test_missing_zero_ibd_pairs_remain_in_full_template_denominator():
    rows, predictions = score_fixture()
    result = score_predictions(predictions, rows, expected_per_label=2)
    assert result["accuracy"] == 1
    assert result["template_accuracy_counting_unpredicted_as_wrong"] == 0.5
    assert result["expected_template_pairs"] == len(LABELS) * 2
    assert result["zero_ibd_omitted_by_label"] == dict.fromkeys(LABELS, 1)
    assert result["multiclass_brier_sum"] == result["negative_log_likelihood"] == 0


def test_scorer_accepts_emitted_nonfirst_tied_maximum():
    rows, predictions = score_fixture()
    predictions[0]["probabilities"] = dict.fromkeys(LABELS, 0.0)
    predictions[0]["probabilities"][LABELS[0]] = 0.5
    predictions[0]["probabilities"][LABELS[1]] = 0.5
    predictions[0]["relation"] = LABELS[1]
    predictions[0]["confidence"] = 0.5
    result = score_predictions(predictions, rows)
    assert result["accuracy"] == pytest.approx((len(LABELS) - 1) / len(LABELS))
    assert result["confusion"][LABELS[0]][LABELS[1]] == 1


@pytest.mark.parametrize(
    "failure",
    ["missing", "duplicate", "label", "confidence", "nan", "probability_class"],
)
def test_scorer_rejects_incomplete_or_inconsistent_predictions(failure):
    rows, predictions = score_fixture()
    if failure == "missing":
        predictions.pop()
    elif failure == "duplicate":
        predictions[-1] = predictions[0]
    elif failure == "label":
        predictions[0]["relation"] = LABELS[1]
    elif failure == "confidence":
        predictions[0]["confidence"] = 0.3
    elif failure == "nan":
        predictions[0]["probabilities"][LABELS[0]] = float("nan")
    else:
        predictions[0]["probabilities"].pop(LABELS[0])
    with pytest.raises(ValueError):
        score_predictions(predictions, rows)


def test_pair_identity_and_relationship_contract():
    rows, _ = score_fixture()
    assert len(validate_pairs(rows)) == len(LABELS)
    with pytest.raises(ValueError, match="unique"):
        validate_pairs(rows + rows[:1])
    rows[0]["indID2"] = "avuncular2_g2"
    with pytest.raises(ValueError, match="labels|canonical"):
        validate_pairs(rows)


def test_training_reference_excludes_each_calibration_replicate(monkeypatch):
    from scripts import train_task12_histogram_candidate as module

    clean = {str(i): [{"label": "siblings", "source": i}] for i in range(3)}
    seen = []
    monkeypatch.setattr(
        module, "read_pedsim_with_lengths", lambda path, genome: (clean[path], {})
    )
    monkeypatch.setattr(module, "validate_pairs", lambda rows: {rows[0]["source"]})
    monkeypatch.setattr(module, "augment_empirical", lambda rows, *args: rows)

    def corrected(rows, reference, **kwargs):
        seen.append((rows[0]["source"], {r["source"] for r in reference}))
        return rows

    monkeypatch.setattr(module, "correct_records_v2", corrected)
    settings = {
        "training_seed": 122226,
        "genome_cm": 3300,
        "augmentations": 1,
        "background_profile": [],
        "minimum_cm": 5,
    }
    rows, _ = build_training_records(list(clean), settings)
    assert len(rows) == 6
    assert seen == [(0, {1, 2}), (1, {0, 2}), (2, {0, 1})]


def test_candidate_lock_rejects_tampered_model_before_confirmation(tmp_path):
    protocol = {"source_sha256": {"x": "source"}, "runtime_manifest_sha256": "runtime"}
    protocol_path = tmp_path / "protocol.json"
    write_json(protocol_path, protocol)
    (tmp_path / "model.joblib").write_bytes(b"fitted model")
    write_json(tmp_path / "training.json", {"complete": True})
    lock = {
        "protocol_sha256": sha256(protocol_path),
        "source_sha256": protocol["source_sha256"],
        "runtime_manifest_sha256": "runtime",
        "model_sha256": sha256(tmp_path / "model.joblib"),
        "training_sha256": sha256(tmp_path / "training.json"),
    }
    write_json(tmp_path / "candidate_lock.json", lock)
    assert validate_candidate_lock(tmp_path, protocol_path, protocol) == lock
    (tmp_path / "model.joblib").write_bytes(b"other model")
    with pytest.raises(ValueError, match="lock"):
        validate_candidate_lock(tmp_path, protocol_path, protocol)


def test_source_input_and_runtime_are_locked(monkeypatch, tmp_path):
    from scripts import task12_histogram_protocol as module

    (tmp_path / "source.py").write_text("source")
    (tmp_path / "input.seg").write_text("simulation")
    write_json(tmp_path / "runtime.json", {"python": "frozen"})
    protocol = {
        "source_sha256": {"source.py": sha256(tmp_path / "source.py")},
        "input_sha256": {"input.seg": sha256(tmp_path / "input.seg")},
        "runtime_manifest": "runtime.json",
        "runtime_manifest_sha256": sha256(tmp_path / "runtime.json"),
    }
    monkeypatch.setattr(module, "runtime_manifest", lambda: {"python": "frozen"})
    validate_protocol(protocol, tmp_path, tmp_path)
    monkeypatch.setattr(module, "runtime_manifest", lambda: {"python": "changed"})
    with pytest.raises(ValueError, match="runtime"):
        validate_protocol(protocol, tmp_path, tmp_path)
    (tmp_path / "input.seg").write_text("other simulation")
    with pytest.raises(ValueError, match="input"):
        validate_protocol(protocol, tmp_path, tmp_path)
    (tmp_path / "source.py").write_text("other source")
    with pytest.raises(ValueError, match="source"):
        validate_protocol(protocol, tmp_path, tmp_path)


def test_evidence_writer_refuses_replacement_and_nonfinite_json(tmp_path):
    path = tmp_path / "evidence.json"
    write_json(path, {"value": 1})
    with pytest.raises(FileExistsError):
        write_json(path, {"value": 2})
    assert json.loads(path.read_text()) == {"value": 1}
    with pytest.raises(ValueError):
        write_json(tmp_path / "bad.json", {"value": float("nan")})
    assert not (tmp_path / "bad.json").exists()
