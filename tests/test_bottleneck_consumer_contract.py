from __future__ import annotations

import json
import sys

import numpy as np
import pytest

from ghist_bottleneck.optimization import neighborhood, optimizer_contract
from scripts import bottleneck_consumer_guard as guard


def valid_fit(engine="moments", restarts=2):
    contract = optimizer_contract(engine)
    optimizer = {"warnflag": 0, "converged": True}
    point = np.log([.4, .1])
    probe = neighborhood(lambda value: 1. + float(np.sum((value - point) ** 2)), point)
    attempts = [{"restart": origin, "arm": arm, "state": "COMPLETED", "parameters": [.4, .1],
                 "start_log_parameters": (point if arm == "refined" else point + [.001 * origin, 0.]).tolist(),
                 "objective": 1., "optimizer": optimizer, "neighborhood": probe,
                 "theta": 136800.}
                for origin in range(restarts) for arm in contract["epsilons"]]
    return {
        "engine": engine, "optimizer_contract": contract,
        "generations": 2000., "post_decline_fraction": .4,
        "ancestral_effective_size": 10000., "best_restart_converged": True,
        "folded": True, "mutation_rate": 3.42e-8, "sequence_length": 100000000,
        "sample_haplotypes": 114,
        "optimizer": optimizer, "restarts": restarts, "best_restart": 0,
        "best_optimizer_arm": "coarse", "multinomial_log_likelihood": -1.,
        "scaled_time_2N": .1, "theta": 136800.,
        "successful_restarts": restarts, "converged_restarts": restarts,
        "optimizer_attempts": attempts,
        "optimizer_diagnostics": {"restarts_within_0_1_log_likelihood": restarts},
        "optimizer_neighborhood": probe,
    }


def bundle(fits=None, producer="scripts/evaluate_moments_chromosomes.py"):
    fits = [valid_fit() for _ in range(3)] if fits is None else fits
    return {
        "cell_id": np.arange(len(fits)),
        "coalescent_ploidy": np.int64(2),
        "raw_latent": np.asarray([guard.latent(fit) for fit in fits]),
        "estimate": np.asarray([guard.physical(fit) for fit in fits]),
        "target_latent": np.asarray([guard.latent(fit) for fit in fits]),
        **guard.bundle_metadata(fits, producer=producer),
    }


def ensemble_inputs():
    mechanical = bundle()
    learned = {
        "cell_id": np.arange(4), "scenario": np.asarray([0, 0, 0, 1]),
        "target_latent": np.tile(mechanical["target_latent"][0], (4, 1)),
        "summary_calibrated_latent": np.tile(mechanical["target_latent"][0], (4, 1)),
    }
    split = {"split_cells": {"calibration": [0, 3], "selection": [1], "evaluation": [2]},
             "observed_summary": {"generations": 2000., "post_decline_fraction": .4}}
    return mechanical, learned, split, valid_fit()


def test_bundle_retains_full_fit_and_matches_array():
    data = bundle()
    contract, fits, _ = guard.require_bundle(data, producer="scripts/evaluate_moments_chromosomes.py",
                                           engine="moments", latent_field="raw_latent")
    assert contract == optimizer_contract("moments")
    assert len(fits) == len(data["cell_id"]) == 3


@pytest.mark.parametrize("change", ["legacy", "stale_source", "nonconverged", "alias", "agreement", "neighborhood", "engine"])
def test_rejects_unsupported_fit_even_when_expected_self_matches(change):
    fit = valid_fit()
    if change == "legacy":
        fit.pop("optimizer_contract")
    elif change == "stale_source":
        key = next(iter(fit["optimizer_contract"]["source_sha256"]))
        fit["optimizer_contract"]["source_sha256"][key] = "0" * 64
    elif change == "nonconverged":
        fit["optimizer"] = {"warnflag": 1, "converged": False}
        fit["best_restart_converged"] = False
    elif change == "alias":
        fit["best_restart_converged"] = False
    elif change == "agreement":
        fit["optimizer_diagnostics"]["restarts_within_0_1_log_likelihood"] = 1
    elif change == "neighborhood":
        fit["optimizer_neighborhood"]["passes_finite_neighborhood_check"] = False
    else:
        fit["engine"] = "dadi"
    with pytest.raises(ValueError):
        guard.require_fit(fit, expected=fit.get("optimizer_contract"))


@pytest.mark.parametrize("change", ["no_schema", "wrong_source", "missing_fit", "duplicate_cell", "array", "nonfinite"])
def test_bundle_rejects_stale_or_incomplete_evidence(change):
    data = bundle()
    if change == "no_schema":
        data.pop("consumer_schema")
    elif change == "wrong_source":
        context = guard.scalar_json(data, "artifact_context_json")
        context["source_sha256"]["scripts/evaluate_moments_chromosomes.py"] = "bad"
        data["artifact_context_json"] = np.asarray(json.dumps(context))
    elif change == "missing_fit":
        data["fit_records_json"] = np.asarray(json.dumps([valid_fit()]))
    elif change == "duplicate_cell":
        data["cell_id"] = np.asarray([0, 1, 1])
    elif change == "array":
        data["raw_latent"][0, 0] += 1
    else:
        data["raw_latent"][0, 0] = np.nan
    with pytest.raises(ValueError):
        guard.require_bundle(data, producer="scripts/evaluate_moments_chromosomes.py", latent_field="raw_latent")


def test_complete_nominal_ensemble_retains_stress_cells_in_full_partition():
    _, _, splits = guard.require_ensemble_inputs(*ensemble_inputs())
    assert {name: values.tolist() for name, values in splits.items()} == {
        "calibration": [0], "selection": [1], "evaluation": [2]}


@pytest.mark.parametrize("change", [{"folded": False}, {"mutation_rate": 1e-8},
                                   {"sequence_length": 1000000}, {"callable_sequence_length": 99999999},
                                   {"sample_haplotypes": 112}, {"sfs": {"projection": 112}}])
def test_calibration_rejects_changed_data_convention(change):
    mechanical, learned, split, observed = ensemble_inputs()
    observed.update(change)
    with pytest.raises(ValueError, match="folded|exposure|masked|projection"):
        guard.require_ensemble_inputs(mechanical, learned, split, observed)


@pytest.mark.parametrize("change", ["missing_mechanistic", "extra_mechanistic", "missing_learned", "duplicate_learned", "target", "overlap", "unassigned", "missing_split", "nonfinite", "observed"])
def test_ensemble_cannot_silently_use_common_subset(change):
    mechanical, learned, split, observed = ensemble_inputs()
    if change == "missing_mechanistic":
        mechanical = bundle([valid_fit(), valid_fit()])
    elif change == "extra_mechanistic":
        mechanical = bundle([valid_fit() for _ in range(4)])
    elif change == "missing_learned":
        learned = {key: value[:3] for key, value in learned.items()}
    elif change == "duplicate_learned":
        learned["cell_id"][-1] = 0
    elif change == "target":
        learned["target_latent"][0, 0] += 1
    elif change == "overlap":
        split["split_cells"]["selection"].append(0)
    elif change == "unassigned":
        split["split_cells"]["calibration"] = [0]
    elif change == "missing_split":
        del split["split_cells"]["evaluation"]
    elif change == "nonfinite":
        learned["summary_calibrated_latent"][0, 0] = np.inf
    else:
        split["observed_summary"]["post_decline_fraction"] = 1.1
    with pytest.raises(ValueError):
        guard.require_ensemble_inputs(mechanical, learned, split, observed)


def test_input_snapshot_rejects_midrun_changes_and_existing_outputs(tmp_path):
    path = tmp_path / "input"
    path.write_text("original")
    hashes = guard.capture_inputs(path)
    path.write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        guard.require_unchanged(hashes)
    with pytest.raises(FileExistsError):
        guard.require_new_outputs(path)


def test_cell_failure_is_retained_without_calibration(monkeypatch):
    from scripts import evaluate_moments_chromosomes as evaluation
    from scripts import run_multiscale_sfs_recovery as recovery

    def fail(*args, **kwargs):
        raise RuntimeError("optimizer failed")

    monkeypatch.setattr(evaluation, "fit_two_epoch_spectrum", fail)
    assert evaluation._fit_cell((7, np.ones(5), 2, "moments")) == (7, {"failure": "RuntimeError: optimizer failed"})
    monkeypatch.setattr(recovery, "simulate_cell", fail)
    assert recovery.recorded_cell((7,)) == {"cell_id": 7, "failure": "RuntimeError: optimizer failed"}


def test_learned_producer_binds_npz_and_rejects_replacement(tmp_path, monkeypatch):
    from scripts import evaluate_bottleneck_chromosomes as producer

    dataset = tmp_path / "data.npz"
    target = np.asarray([[7., -1.], [8., 0.], [9., 1.]])
    np.savez(dataset, coalescent_ploidy=2, cell_id=np.arange(3), targets=target,
             nuisance=np.zeros((3, 1)), sfs=np.ones((3, 100)), windows=np.ones((3, 2, 2)),
             n_variants=np.ones(3))
    checkpoint = tmp_path / "checkpoint"
    checkpoint.mkdir()
    for name in ("model_0.pt", "normalization.npz", "reference.joblib", "experiment.json"):
        (checkpoint / name).write_text("synthetic software fixture")
    output = tmp_path / "learned.json"
    monkeypatch.setattr(sys, "argv", ["evaluate", "--dataset", str(dataset), "--checkpoint-dir", str(checkpoint), "--out", str(output)])
    monkeypatch.setattr(producer, "stratified_cell_split", lambda *args: tuple(np.asarray([i]) for i in range(3)))
    monkeypatch.setattr(producer, "_network_predictions", lambda *args: (None, None, None, None))
    monkeypatch.setattr(producer, "_aggregate_cells", lambda *args: (target, np.ones_like(target)))
    monkeypatch.setattr(producer, "fit_affine", lambda *args: (np.eye(2), np.zeros(2)))
    monkeypatch.setattr(producer, "_metrics", lambda *args, **kwargs: {})
    monkeypatch.setattr(producer, "predict_blocks", lambda *args, **kwargs: {
        "generations": 2000., "post_decline_fraction": .4})
    producer.main()
    record = json.loads(output.read_text())
    predictions = output.with_suffix(".npz")
    guard.require_learned_binding(predictions, record)
    predictions.write_bytes(b"replaced")
    with pytest.raises(ValueError, match="hash-bound"):
        guard.require_learned_binding(predictions, record)
    with pytest.raises(FileExistsError):
        producer.main()


def test_multiscale_cli_rejects_legacy_before_training(tmp_path, monkeypatch):
    from scripts import evaluate_multiscale_sfs_recovery as evaluation

    shard = tmp_path / "legacy.npz"
    np.savez(shard, coalescent_ploidy=2)
    observed = tmp_path / "observed.json"
    observed.write_text(json.dumps(valid_fit()))
    monkeypatch.setattr(sys, "argv", ["evaluate", "--input-glob", str(shard), "--observed", str(observed),
                                      "--out", str(tmp_path / "report.json"), "--model-out", str(tmp_path / "model.joblib")])
    with pytest.raises(ValueError, match="Legacy"):
        evaluation.main()
    assert not (tmp_path / "model.joblib").exists()


@pytest.mark.parametrize("changed", [{"seed": 200}, {"blocks": 25, "block_length": 4000000.}, {"restarts": 12}])
def test_multiscale_cli_rejects_mixed_simulation_designs(tmp_path, monkeypatch, changed):
    from scripts import evaluate_multiscale_sfs_recovery as evaluation
    from scripts.run_multiscale_sfs_recovery import simulation_design

    baseline = dict(seed=100, blocks=50, block_length=2000000., mutation_rate=3.42e-8, restarts=8)
    for index, settings in enumerate((baseline, {**baseline, **changed})):
        design = simulation_design(**settings)
        fits = [valid_fit(restarts=settings["restarts"]) for _ in range(3)]
        for cell, fit in enumerate(fits):
            fit.update(optimizer_seed=settings["seed"] + 900000 + cell, restarts=settings["restarts"])
        data = bundle(fits, producer="scripts/run_multiscale_sfs_recovery.py")
        data.update(guard.bundle_metadata(fits, producer="scripts/run_multiscale_sfs_recovery.py", simulation_design=design))
        np.savez(tmp_path / f"shard_{index}.npz", **data)
    observed = tmp_path / "observed.json"
    observed.write_text(json.dumps(valid_fit()))
    monkeypatch.setattr(sys, "argv", ["evaluate", "--input-glob", str(tmp_path / "shard_*.npz"),
        "--observed", str(observed), "--out", str(tmp_path / "report.json"),
        "--model-out", str(tmp_path / "model.joblib")])
    with pytest.raises(ValueError, match="different simulation"):
        evaluation.main()
    assert not (tmp_path / "model.joblib").exists()


def test_multiscale_design_rejects_forged_fit_seed():
    from scripts.run_multiscale_sfs_recovery import require_simulation_design, simulation_design

    design = simulation_design(seed=100, blocks=50, block_length=2000000., mutation_rate=3.42e-8, restarts=8)
    fits = [valid_fit(restarts=8)]
    fits[0].update(optimizer_seed=900101, restarts=8)
    data = bundle(fits, producer="scripts/run_multiscale_sfs_recovery.py")
    data.update(guard.bundle_metadata(fits, producer="scripts/run_multiscale_sfs_recovery.py", simulation_design=design))
    _, retained, context = guard.require_bundle(data, producer="scripts/run_multiscale_sfs_recovery.py")
    with pytest.raises(ValueError, match="seed/restarts"):
        require_simulation_design(data, retained, context)


def test_mechanistic_cli_rejects_nonconverged_observed_before_fitting(tmp_path, monkeypatch):
    from scripts import evaluate_moments_chromosomes as evaluation

    observed = valid_fit()
    observed["optimizer"] = {"warnflag": 1, "converged": False}
    observed["best_restart_converged"] = False
    observed["converged_restarts"] = 0
    for attempt in observed["optimizer_attempts"]:
        attempt["optimizer"] = observed["optimizer"]
    path = tmp_path / "observed.json"
    path.write_text(json.dumps(observed))
    dataset, split = tmp_path / "dataset", tmp_path / "split"
    dataset.write_text("not read before rejection")
    split.write_text("not read before rejection")
    monkeypatch.setattr(sys, "argv", ["evaluate", "--dataset", str(dataset), "--split-record", str(split),
                                      "--observed-fit", str(path), "--out", str(tmp_path / "report.json")])
    with pytest.raises(ValueError, match="converge"):
        evaluation.main()


def test_ensemble_cli_requires_bound_full_inputs_and_remains_provisional(tmp_path, monkeypatch):
    from scripts import evaluate_bottleneck_ensemble as evaluation

    mechanical, learned, split, observed = ensemble_inputs()
    learned["coalescent_ploidy"] = np.int64(2)
    learned_path = tmp_path / "learned.npz"
    np.savez(learned_path, **learned)
    split.update(coalescent_ploidy=2, prediction_npz_sha256=guard.sha256(learned_path),
                 evaluator_sha256=guard.sha256(guard.ROOT / "scripts/evaluate_bottleneck_chromosomes.py"),
                 dataset_sha256="0" * 64)
    split_path = tmp_path / "split.json"
    split_path.write_text(json.dumps(split))
    observed_path = tmp_path / "observed.json"
    observed_path.write_text(json.dumps(observed))
    mechanical.update(guard.bundle_metadata([valid_fit() for _ in range(3)],
        producer="scripts/evaluate_moments_chromosomes.py", inputs={
            "dataset_sha256": split["dataset_sha256"],
            "split_record_sha256": guard.sha256(split_path),
            "observed_fit_sha256": guard.sha256(observed_path),
        }))
    mechanical_path = tmp_path / "mechanistic.npz"
    np.savez(mechanical_path, **mechanical)
    output = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["evaluate", "--moments", str(mechanical_path),
        "--learned", str(learned_path), "--split-record", str(split_path),
        "--observed-moments", str(observed_path), "--out", str(output)])
    evaluation.main()
    report = json.loads(output.read_text())
    assert report["common_nominal_cells"] == 3
    assert report["split_sizes"] == {"calibration": 1, "selection": 1, "evaluation": 1}
    assert report["scientifically_validated"] is report["release_ready"] is False
    output.unlink()
    observed["unrelated_note"] = "different source bytes, same valid numerical fit"
    observed_path.write_text(json.dumps(observed))
    with pytest.raises(ValueError, match="inputs differ"):
        evaluation.main()
    assert not output.exists()
