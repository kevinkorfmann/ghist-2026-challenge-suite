import json

import pytest

from scripts.fit_human_family_relationships import LABELS
from scripts.task12_histogram_protocol import (
    compare_runs,
    load_frozen_protocol,
    require_protocol_digest,
    score_predictions,
    sha256,
    validate_pairs,
    validate_protocol,
    validate_seed_registry,
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


def test_positive_seed_gate_uses_exact_counts_for_observed_floating_tie():
    import numpy as np

    historical = [2617, 2546, 2496, 2447, 2476]
    candidate = [2632, 2580, 2478, 2433, 2459]
    assert sum(historical) == sum(candidate)
    assert (
        np.mean([new / 6797 - old / 6797 for old, new in zip(historical, candidate)])
        > 0
    )
    runs = []
    for old, new in zip(historical, candidate):
        models = {}
        for name, count in (("historical", old), ("full_histogram_v2", new)):
            models[name] = {
                "accuracy": count / 6797,
                "meioses_accuracy": 0.5,
                "segment_present_pairs": 6797,
                "confusion": {
                    label: {label: count if index == 0 else 0}
                    for index, label in enumerate(LABELS)
                },
            }
        runs.append({"seed": 2026115012, "models": models})
    settings = {
        "mean_exact_gain_min": -1,
        "positive_mean_gain_fraction_by_seed_min": 1,
        "worst_regime_exact_gain_min": -1,
        "mean_meioses_gain_min": -1,
    }
    result = compare_runs(runs, "historical", settings)
    assert result["mean_exact_gain_by_seed"] == {"2026115012": 0.0}
    assert result["exact_correct_count_gain_by_seed"] == {"2026115012": 0}
    assert not result["gate"]["positive_seed_fraction"]
    runs[0]["models"]["full_histogram_v2"]["confusion"][LABELS[0]][LABELS[0]] += 1
    assert compare_runs(runs, "historical", settings)["gate"]["positive_seed_fraction"]


def test_exact_seed_gate_refuses_different_pair_denominators():
    runs = [
        {
            "seed": 1,
            "models": {
                name: {
                    "accuracy": 0.5,
                    "meioses_accuracy": 0.5,
                    "segment_present_pairs": count,
                }
                for name, count in (("historical", 3), ("full_histogram_v2", 4))
            },
        }
    ]
    with pytest.raises(ValueError, match="denominators"):
        compare_runs(runs, "historical", {})


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
    protocol_digest = sha256(protocol_path)
    assert (
        validate_candidate_lock(tmp_path, protocol_path, protocol, protocol_digest)
        == lock
    )
    (tmp_path / "model.joblib").write_bytes(b"other model")
    with pytest.raises(ValueError, match="lock"):
        validate_candidate_lock(tmp_path, protocol_path, protocol, protocol_digest)


def test_source_input_and_runtime_are_locked(monkeypatch, tmp_path):
    from scripts import task12_histogram_protocol as module

    monkeypatch.setattr(module, "validate_seed_registry", lambda *args: None)

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


def test_protocol_bytes_are_bound_at_parse_and_mutation_rejected(tmp_path):
    path = tmp_path / "protocol.json"
    path.write_text('{"training_seed": 1}')
    parsed, digest = load_frozen_protocol(path)
    assert parsed == {"training_seed": 1}
    require_protocol_digest(path, digest)
    path.write_text('{"training_seed": 2}')
    with pytest.raises(ValueError, match="changed after"):
        require_protocol_digest(path, digest)
    # A rewritten lock cannot bind this already-parsed run to the new config.
    with pytest.raises(ValueError, match="changed after"):
        validate_candidate_lock(tmp_path, path, parsed, digest)


def test_train_main_rejects_mid_fit_protocol_change_before_publication(
    monkeypatch, tmp_path
):
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from scripts import train_task12_histogram_candidate as module

    protocol_path = tmp_path / "protocol.json"
    output = tmp_path / "candidate"
    settings = {
        "artifact_directory": str(output),
        "background_profile": [1.0],
        "calibration": ["a.seg"],
        "metadata": {},
        "training_seed": 1,
    }
    protocol = {
        "training": settings,
        "confirmation": {"comparators": [{"path": "old.joblib"}]},
        "fixed_hyperparameters": {},
        "expected_training_rows": len(LABELS),
        "training_regimes_reference": [],
    }
    protocol_path.write_text(json.dumps(protocol))
    historical = {
        "corrected_payload": {
            "kind": "direct",
            "model": HistGradientBoostingClassifier(),
        },
        "corrected_model_weight": 1.0,
        "background_profile": [1.0],
        "feature_names": module.corrected_feature_names(5),
    }

    class Estimator:
        def set_params(self, **kwargs):
            return self

        def get_params(self):
            return {}

        def fit(self, *args, **kwargs):
            changed = json.loads(protocol_path.read_text())
            changed["training"]["training_seed"] = 999
            protocol_path.write_text(json.dumps(changed))

    monkeypatch.setattr(module, "validate_protocol", lambda *args: None)
    monkeypatch.setattr(module.joblib, "load", lambda *args: historical)
    monkeypatch.setattr(
        module.joblib,
        "dump",
        lambda *args: pytest.fail("Model published after protocol changed"),
    )
    monkeypatch.setattr(module, "clone", lambda *args: Estimator())
    monkeypatch.setattr(module, "validate_simulation_splits", lambda *args: {})
    monkeypatch.setattr(
        module,
        "build_training_records",
        lambda *args: ([{"label": label} for label in LABELS], []),
    )
    monkeypatch.setattr(module, "matrix", lambda *args: np.ones((len(LABELS), 2)))
    monkeypatch.setattr(
        "sys.argv",
        [
            "train",
            "--root",
            str(tmp_path),
            "--protocol",
            str(protocol_path),
            "--out-dir",
            str(output),
        ],
    )
    with pytest.raises(ValueError, match="changed after"):
        module.main()
    assert not (output / "model.joblib").exists()
    assert not (output / "candidate_lock.json").exists()


def seed_registry_fixture(tmp_path):
    log = tmp_path / "a.log"
    log.write_text("Random seed:\t101001\n")
    row = {"path": "a.log", "sha256": sha256(log), "seed": 101001}
    return {
        "training": {"calibration": ["a.seg"]},
        "confirmation": {"seeds": [222], "augmentation_seed": 333, "regimes": [{}, {}]},
        "seed_registry": {
            "calibration_logs": [row],
            "confirmation_seeds": [222],
            "augmentation_seeds": [333, 334],
        },
        "input_sha256": {"a.log": row["sha256"]},
        "prior_confirmation_seeds_excluded": [111],
    }


@pytest.mark.parametrize(
    "failure",
    [
        "log_content",
        "parsed_seed",
        "calibration_collision",
        "augmentation_collision",
        "previous_confirmation",
        "missing_log",
    ],
)
def test_simulation_seed_registry_is_executable_and_disjoint(tmp_path, failure):
    protocol = seed_registry_fixture(tmp_path)
    validate_seed_registry(protocol, tmp_path)
    registry = protocol["seed_registry"]
    if failure == "log_content":
        (tmp_path / "a.log").write_text("Random seed: 101002\n")
    elif failure == "parsed_seed":
        registry["calibration_logs"][0]["seed"] = 101002
    elif failure == "calibration_collision":
        protocol["confirmation"]["seeds"] = registry["confirmation_seeds"] = [101001]
    elif failure == "augmentation_collision":
        protocol["confirmation"]["augmentation_seed"] = 222
        registry["augmentation_seeds"] = [222, 223]
    elif failure == "previous_confirmation":
        protocol["prior_confirmation_seeds_excluded"] = [222]
    else:
        registry["calibration_logs"] = []
    with pytest.raises(ValueError, match="seed|random"):
        validate_seed_registry(protocol, tmp_path)
