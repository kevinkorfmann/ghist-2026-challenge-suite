from __future__ import annotations

import csv
import json
from pathlib import Path

import joblib
import numpy as np
import pytest
import sklearn
from sklearn.ensemble import RandomForestClassifier

from ghist_relatedness import ibd
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.task9_release import input_pairs, load_release, sha256
from scripts import predict_task9_current
from scripts.prepare_task9_release import release_gates


def segment(**changes):
    return dict(indID1="a", indID2="b", chr="1", start_bp="0", end_bp="20",
                type="IBD1", start_cM="0", end_cM="10", length="10", **changes)


@pytest.mark.parametrize("changes", [
    {"start_bp": "nan"}, {"end_bp": "-1"}, {"end_cM": "9"},
    {"length": "inf"}, {"indID1": " a"}, {"type": "HBD"},
])
def test_input_contract_rejects_bad_segments(changes):
    with pytest.raises(ValueError):
        input_pairs([{**segment(), **changes}])


def test_duplicate_reverse_segment_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        input_pairs([segment(), {**segment(), "indID1": "b", "indID2": "a"}])


def test_reordered_bom_csv_and_bad_header(tmp_path):
    path = tmp_path / "segments.csv"
    row = segment()
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, list(reversed(row)))
        writer.writeheader()
        writer.writerow(row)
    assert ibd.read_segment_rows(path) == [row]
    path.write_text("indID1,indID1\na,b\n")
    with pytest.raises(ValueError, match="columns"):
        ibd.read_segment_rows(path)


@pytest.mark.parametrize("text", ["", "1 20 1 1\n1 0 0 0\n", "1 0 0 nan\n", "1 0 2 2\n1 20 1 3\n"])
def test_invalid_map_rejected(tmp_path, text):
    path = tmp_path / "map"
    path.write_text(text)
    with pytest.raises(ValueError):
        ibd.read_sex_map(path)


def test_probability_shape_and_retired_nan_rejected():
    with pytest.raises(ValueError, match="matrix"):
        predict_task9_current.accepted_probabilities(["siblings"], [[.5, .5]])
    with pytest.raises(ValueError, match="matrix"):
        predict_task9_current.accepted_probabilities(["siblings", "half_sibling"], [[1, np.nan]])


def test_final_requires_provenance_before_loading_model(tmp_path, monkeypatch):
    monkeypatch.setattr("sys.argv", ["predict_task9_current", "--model", "absent",
                                    "--observed", "absent", "--sex-map", "absent",
                                    "--dataset", "final", "--out-dir", str(tmp_path)])
    with pytest.raises(ValueError, match="Synapse entity/version/checksum"):
        predict_task9_current.main()


@pytest.mark.parametrize('entity_field', ['entity_id', 'entity'])
def test_final_accepts_download_manifest_identity(tmp_path, monkeypatch, entity_field):
    observed = tmp_path / 'input.csv'
    observed.write_text('input snapshot')
    provenance = tmp_path / 'manifest.json'
    record = {entity_field: 'syn75202631', 'version': 1, 'sha256': sha256(observed)}
    provenance.write_text(json.dumps(record))
    monkeypatch.setattr('sys.argv', ['predict_task9_current', '--model', 'unused',
                                    '--observed', str(observed), '--sex-map', 'unused',
                                    '--dataset', 'final', '--input-provenance', str(provenance),
                                    '--out-dir', str(tmp_path / 'out')])
    def reached_model(*args, **kwargs):
        raise RuntimeError('validated provenance reached model loading')
    monkeypatch.setattr(joblib, 'load', reached_model)
    with pytest.raises(RuntimeError, match='reached model loading'):
        predict_task9_current.main()
    provenance.write_text(json.dumps({**record, 'entity': 'synWrong', 'entity_id': 'syn75202631'}))
    with pytest.raises(ValueError, match='provenance'):
        predict_task9_current.main()


def test_production_checks_coverage_and_preserves_output_on_failure(tmp_path, monkeypatch):
    path = tmp_path / "segments.csv"
    row = segment()
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(row))
        writer.writeheader()
        writer.writerow(row)
    map_path = tmp_path / "map"
    map_path.write_text("1 0 0 0\n1 20 5 15\n")
    records = ibd.pair_features([row], 3350.933402, ibd.read_sex_map(map_path))
    matrix, features = ibd.matrix(records)
    labels = list(IBD_LABEL_MEIOSES)
    model = RandomForestClassifier(n_estimators=2, random_state=1).fit(
        np.repeat(matrix, len(labels), axis=0), labels)
    model_path = tmp_path / "model"
    joblib.dump({"model": model, "features": features}, model_path)
    pairs = tmp_path / "pairs.csv"
    pairs.write_text("indID1,indID2\na,b\n")
    output = tmp_path / "out"
    argv = ["predict_task9_current", "--model", str(model_path), "--observed", str(path),
            "--sex-map", str(map_path), "--pairs", str(pairs), "--dataset", "testing",
            "--synthetic-input", "--out-dir", str(output)]
    monkeypatch.setattr("sys.argv", argv)
    predict_task9_current.main()
    candidate = output / "task9_testing_candidate.csv"
    before = candidate.read_bytes()
    summary = json.loads((output / "summary.json").read_text())
    assert summary["validation"]["input_pair_coverage_verified"]
    pairs.write_text("indID1,indID2\na,b\na,c\n")
    with pytest.raises(ValueError, match="zero-IBD/unlisted"):
        predict_task9_current.main()
    assert candidate.read_bytes() == before


def test_release_pins_model_map_features_and_environment(tmp_path):
    model = tmp_path / "model"
    model.write_bytes(b"model")
    genetic_map = tmp_path / "map"
    genetic_map.write_bytes(b"map")
    manifest = {
        "schema": "ghist.task9.release.v1", "status": "VALIDATED_RELEASE_CANDIDATE",
        "accepted_labels": list(IBD_LABEL_MEIOSES), "pair_domain": "unordered_pairs_present_in_segment_csv",
        "gates": {"all": True}, "environment": {"scikit_learn": sklearn.__version__},
        "model": {"path": "model", "sha256": sha256(model)}, "sex_map_sha256": sha256(genetic_map),
        "feature_contract_sha256": {name: sha256(Path(ibd.__file__).with_name(name)) for name in ("ibd.py", "labels.py")},
        "genome_cm": 3350.933402,
    }
    path = tmp_path / "release.json"
    path.write_text(json.dumps(manifest))
    assert load_release(path, tmp_path, genetic_map)[1] == model
    with pytest.raises(ValueError, match="normalization"):
        load_release(path, tmp_path, genetic_map, 100)
    model.write_bytes(b"different model")
    with pytest.raises(ValueError, match="model checksum"):
        load_release(path, tmp_path, genetic_map)


def test_release_gate_rejects_model_or_provenance_mismatch():
    folder = Path("experiments/scientific_2026_11_20/relatedness")
    training = json.loads((folder / "task9_current_training.json").read_text())
    evidence = [json.loads((folder / f"task9_confirmation_current_{scenario}.json").read_text())
                for scenario in ("interference", "poisson")]
    operational = {"model_sha256": "different", "passed": False, "datasets": []}
    gates = release_gates(training, evidence, operational, evidence[0]["model_sha256"])
    assert gates["both_exact_accuracy_at_least_80pct"]
    assert not gates["same_fixed_model"]
    assert not gates["production_pair_coverage_and_format"]
    training["source_splits"]["calibration"][0]["sha256"] = evidence[0]["segment_sha256"]
    assert not release_gates(training, evidence, operational, evidence[0]["model_sha256"])["distinct_calibration_and_confirmation"]
