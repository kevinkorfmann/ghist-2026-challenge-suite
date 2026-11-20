import copy

import numpy as np
import pytest

from ghist_bottleneck import optimization as opt


def optimizer_info(converged=True):
    return dict(warnflag=0 if converged else 2, task="CONVERGENCE" if converged else "ABNORMAL",
                funcalls=4, nit=1, grad=np.zeros(2))


def passing_neighborhood(*args):
    point = np.log([.4, .2]) if not args else args[-1]
    # Save a real probe geometry around the requested endpoint.
    return REAL_NEIGHBORHOOD(lambda p: float(np.sum((p - point)**2)), point)


REAL_NEIGHBORHOOD = opt.neighborhood


def test_paired_candidates_preserve_unconverged_winner_and_distinct_ancestry(monkeypatch):
    endpoints = iter([(-2.0, True), (-3.0, False), (-1.5, True), (-1.6, True), (-2.0, True), (-1.7, True)])
    def optimize(*args, **kwargs):
        value, converged = next(endpoints)
        return np.array([value, -2.0]), 0, optimizer_info(converged)
    monkeypatch.setattr(opt, "fmin_l_bfgs_b", optimize)
    monkeypatch.setattr(opt, "neighborhood", passing_neighborhood)
    result = opt.run_optimizer(lambda point: point[0], engine="moments", starts=np.array([[-2., -2.], [-1., -1.]]))
    assert len(result["attempts"]) == 6
    assert len(result["origin_best"]) == 2
    assert result["best"]["arm"] == "fine"
    assert not result["best"]["optimizer"]["converged"]
    assert sum(row["optimizer"]["converged"] for row in result["origin_best"]) == 1
    assert result["arm_summary"]["fine"]["reported_converged"] == 1


def test_failed_refinement_and_probe_never_erase_finite_ancestor(monkeypatch):
    calls = []
    def optimize(*args, **kwargs):
        calls.append(kwargs["epsilon"])
        if len(calls) == 3:
            raise ArithmeticError("refinement failed")
        return np.array([-2., -2.]), 0, optimizer_info()
    def failed_probe(*args):
        raise ArithmeticError("probe failed after finite endpoint")
    monkeypatch.setattr(opt, "fmin_l_bfgs_b", optimize)
    monkeypatch.setattr(opt, "neighborhood", failed_probe)
    result = opt.run_optimizer(lambda point: point[0], engine="moments", starts=np.array([[-2., -2.]]))
    assert result["best"]["arm"] == "coarse"  # original-first exact tie
    assert result["best"]["objective"] == -2
    assert result["best"]["neighborhood"]["status"] == "failed"
    assert result["attempts"][2]["state"] == "FAILED"
    assert len(result["origin_best"]) == 1
    assert result["arm_summary"]["coarse"]["neighborhood_failed"] == 1


def test_missing_coarse_endpoint_is_explicitly_skipped_without_losing_fine(monkeypatch):
    calls = []
    def optimize(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ArithmeticError("coarse failed")
        return np.array([-2., -2.]), 0, optimizer_info()
    monkeypatch.setattr(opt, "fmin_l_bfgs_b", optimize)
    monkeypatch.setattr(opt, "neighborhood", passing_neighborhood)
    result = opt.run_optimizer(lambda point: point[0], engine="dadi", starts=np.array([[-2., -2.]]))
    assert result["best"]["arm"] == "fine"
    assert result["attempts"][-1]["state"] == "SKIPPED"
    assert result["failed_origins"] == []


def test_seeded_starts_keep_historical_warm_start_consumption():
    rng = np.random.default_rng(1729)
    expected = np.log(np.exp(rng.uniform(np.log(opt.LOWER), np.log(opt.UPPER), size=(2, 2))))
    starts = opt.log_starts(3, 1729, initial=(.999, .2))
    np.testing.assert_array_equal(starts[1:], expected)
    np.testing.assert_array_equal(starts[0], np.log([.999, .2]))
    np.testing.assert_array_equal(opt.log_starts(2, 3, starts=expected), expected)
    with pytest.raises(ValueError, match="mutually exclusive"):
        opt.log_starts(2, 3, starts=expected, initial=(.4, .2))


def valid_fit(engine="moments", strategy=opt.DEFAULT_STRATEGY):
    contract = opt.optimizer_contract(engine, strategy)
    optimizer = {"converged": True, "warnflag": 0}
    probe = passing_neighborhood()
    attempts = [{"restart": origin, "arm": arm, "state": "COMPLETED", "parameters": [.4, .2],
                 "objective": 1., "theta": .4, "optimizer": optimizer, "neighborhood": probe,
                 "start_log_parameters": [-2. - origin, -2. - origin]}
                for origin in range(2) for arm in contract["epsilons"]]
    return {"engine": engine, "optimizer_contract": contract, "sample_haplotypes": 10,
            "theta": .4, "mutation_rate": 1e-8, "sequence_length": 1000,
            "ancestral_effective_size": 10000., "generations": 4000.,
            "optimizer": optimizer, "best_restart_converged": True,
            "restarts": 2, "best_restart": 0, "best_optimizer_arm": "coarse",
            "multinomial_log_likelihood": -1., "post_decline_fraction": .4, "scaled_time_2N": .2,
            "successful_restarts": 2, "converged_restarts": 2, "optimizer_attempts": attempts,
            "optimizer_diagnostics": {"restarts_within_0_1_log_likelihood": 2},
            "optimizer_neighborhood": probe}


@pytest.mark.parametrize("corruption", ["missing", "stale_source", "wrong_engine", "nonconverged", "alias", "agreement", "probe", "physical", "theta", "raw_probe", "displacement"])
def test_contract_blocks_unknown_mismatched_or_failed_promotion(corruption):
    fit = valid_fit()
    if corruption == "missing":
        del fit["optimizer_contract"]
    elif corruption == "stale_source":
        fit["optimizer_contract"]["source_sha256"]["ghist_bottleneck/optimization.py"] = "old"
    elif corruption == "wrong_engine":
        fit["engine"] = "dadi"
    elif corruption == "nonconverged":
        fit["optimizer"].update(converged=False, warnflag=2)
        fit["best_restart_converged"] = False
    elif corruption == "alias":
        fit["best_restart_converged"] = False
    elif corruption == "agreement":
        fit["optimizer_diagnostics"]["restarts_within_0_1_log_likelihood"] = 1
    elif corruption == "raw_probe":
        fit["optimizer_neighborhood"]["probes"][0]["objective_change"] = -1.
    elif corruption == "displacement":
        fit["optimizer_neighborhood"]["probes"][0]["actual_log_displacement"] = [0., 0.]
    elif corruption == "physical":
        fit["generations"] *= 2
    elif corruption == "theta":
        fit["theta"] *= 2
    else:
        fit["optimizer_neighborhood"] = {"status": "failed"}
    with pytest.raises(ValueError):
        opt.validate_optimizer_contract(fit, expected=fit.get("optimizer_contract"), require_agreement=True, require_neighborhood=True)


def test_legacy_contract_requires_explicit_development_opt_in():
    fit = valid_fit(strategy=opt.LEGACY_STRATEGY)
    with pytest.raises(ValueError, match="Legacy"):
        opt.validate_optimizer_contract(fit)
    opt.validate_optimizer_contract(fit, allow_legacy=True)


def test_bootstrap_retains_nonconverged_fit_and_rejects_changed_contract():
    from scripts.bootstrap_bottleneck_sfs import record_replica, validate_reference
    reference = valid_fit()
    reference.update(folded=True, sequence_length=1000, theta=.4, ancestral_effective_size=10000)
    opt.validate_optimizer_contract(reference, require_agreement=True, require_neighborhood=True)
    validate_reference(reference, mutation_rate=1e-8, sequence_length=1000, folded=True)
    failed = copy.deepcopy(reference)
    failed["optimizer"].update(converged=False, warnflag=2)
    failed["best_restart_converged"] = False
    failed["converged_restarts"] = 0
    with pytest.raises(ValueError, match="did not converge"):
        validate_reference(failed, mutation_rate=1e-8, sequence_length=1000, folded=True)
    record = record_replica(failed, expected=reference["optimizer_contract"])
    assert record["fit"] is failed and "intervals omitted" in record["failure"]
    changed = copy.deepcopy(reference["optimizer_contract"])
    changed["strategy"] = "different"
    assert "failure" in record_replica(reference, expected=changed)


def test_asymmetric_boundary_probes_do_not_claim_centered_derivative():
    upper = np.log(opt.UPPER)
    point = upper - [.00005, .1]
    probe = opt.neighborhood(lambda x: float(np.sum((x - upper)**2)), point)
    assert probe["axis_slopes"]["0.0001"][0]["centered_slope"] is None
    assert probe["axis_slopes"]["0.0001"][0]["symmetric"] is False
    assert not probe["passes_finite_neighborhood_check"]


def test_one_failed_bootstrap_replica_suppresses_all_intervals(tmp_path, monkeypatch):
    import hashlib
    import json
    import sys
    from concurrent.futures import Future
    from ghist_inference.checkpoint import spectrum_digest
    from scripts import bootstrap_bottleneck_sfs as bootstrap

    blocks = np.array([[0., 4., 2., 1., 0.], [0., 2., 3., 2., 0.]])
    vcf = tmp_path / "input.vcf"
    vcf.write_text("synthetic input identity; extraction is separately tested")
    reference = valid_fit()
    reference.update(folded=True, polarized=False, sample_haplotypes=4, sequence_length=1000, theta=.4,
                     ancestral_effective_size=10000, generations=4000., post_decline_fraction=.4, scaled_time_2N=.2,
                     input_vcf_sha256=hashlib.sha256(vcf.read_bytes()).hexdigest(),
                     sfs={"projection": 4, "spectrum_sha256": spectrum_digest(blocks.sum(axis=0))})
    reference_path = tmp_path / "fit.json"
    reference_path.write_text(json.dumps(reference))
    failed = copy.deepcopy(reference)
    failed["optimizer"].update(converged=False, warnflag=2)
    failed["best_restart_converged"] = False
    failed["converged_restarts"] = 0
    results = iter([reference, failed])

    class ImmediateExecutor:
        def __init__(self, **kwargs):
            pass
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def submit(self, function, argument):
            future = Future()
            future.set_result(next(results))
            return future

    monkeypatch.setattr(bootstrap, "ProcessPoolExecutor", ImmediateExecutor)
    monkeypatch.setattr(bootstrap, "complete_block_spectra_from_vcf", lambda *a, **k: (blocks, {"projection": 4}))
    output = tmp_path / "bootstrap.json"
    monkeypatch.setattr(sys, "argv", ["bootstrap", "--vcf", str(vcf), "--fit", str(reference_path),
                                      "--out", str(output), "--replicates", "2", "--workers", "1",
                                      "--folded", "--mutation-rate", "1e-8", "--sequence-length", "1000", "--block-length", "500"])
    with pytest.raises(RuntimeError, match="1 bootstrap fits failed"):
        bootstrap.main()
    saved = json.loads(output.read_text())
    assert saved["status"] == "failed" and saved["requested_replicates"] == 2
    assert len(saved["completed_replicates"]) == 2
    assert saved["completed_replicates"]["1"]["fit"]["optimizer"]["converged"] is False
    assert "interval_95" not in saved


@pytest.mark.parametrize("engine", ["moments", "dadi"])
@pytest.mark.parametrize("change_input", [False, True])
def test_vcf_wrappers_bind_polarization_and_input_bytes(tmp_path, monkeypatch, engine, change_input):
    import hashlib
    from ghist_bottleneck import sfs_baseline, dadi_reference
    module = sfs_baseline if engine == "moments" else dadi_reference
    wrapper = module.fit_two_epoch if engine == "moments" else module.fit_two_epoch_dadi
    fitter_name = "fit_two_epoch_spectrum" if engine == "moments" else "fit_two_epoch_dadi_spectrum"
    vcf = tmp_path / "input.vcf"
    vcf.write_text("original bytes")
    spectrum = np.array([0., 4., 2., 1., 0.])
    metadata = {"projection": 4}
    monkeypatch.setattr(module, "projected_sfs_from_vcf", lambda *a, **k: (spectrum, metadata))
    monkeypatch.setattr(module, "spectrum_exposure", lambda *a: 1000)
    def fit(*args, **kwargs):
        if change_input:
            vcf.write_text("changed after extraction")
        return {"folded": kwargs["folded"]}
    monkeypatch.setattr(module, fitter_name, fit)
    kwargs = {"folded": True, "sequence_length": 1000}
    if engine == "moments":
        kwargs["polarized"] = False
    if change_input:
        with pytest.raises(ValueError, match="VCF changed"):
            wrapper(vcf, **kwargs)
    else:
        result = wrapper(vcf, **kwargs)
        assert result["polarized"] is False
        assert result["input_vcf_sha256"] == hashlib.sha256(vcf.read_bytes()).hexdigest()


def test_bootstrap_input_and_projection_guards(tmp_path, monkeypatch):
    import sys
    from scripts import bootstrap_bottleneck_sfs as bootstrap
    path = tmp_path / "input.vcf"
    path.write_text("first")
    locked = {path: bootstrap.file_sha256(path)}
    bootstrap.verify_unchanged(locked)
    path.write_text("changed")
    with pytest.raises(ValueError, match="changed during"):
        bootstrap.verify_unchanged(locked)
    with pytest.raises(ValueError, match="sample_haplotypes"):
        bootstrap.validate_input_reference({"sample_haplotypes": 6, "sfs": {"projection": 4}},
                                           path, np.zeros((2, 5)), {"projection": 4}, allow_unverified=True)
    output = tmp_path / "existing.json"
    output.write_text("preserve")
    monkeypatch.setattr(sys, "argv", ["bootstrap", "--vcf", str(path), "--fit", "missing.json", "--out", str(output)])
    with pytest.raises(FileExistsError, match="already exists"):
        bootstrap.main()
    assert output.read_text() == "preserve"


def test_duplicate_original_starts_cannot_create_restart_agreement():
    starts = np.array([[-2., -2.], [-2., -2.]])
    with pytest.raises(ValueError, match="must be distinct"):
        opt.log_starts(2, 1, starts=starts)
    with pytest.raises(ValueError, match="must be distinct"):
        opt.run_optimizer(lambda point: float(np.sum(point**2)), engine="moments", starts=starts)
    fit = valid_fit()
    for row in fit["optimizer_attempts"]:
        if row["arm"] in {"coarse", "fine"}:
            row["start_log_parameters"] = [-2., -2.]
    with pytest.raises(ValueError, match="must be distinct"):
        opt.validate_optimizer_contract(fit, require_agreement=True, require_neighborhood=True)


def test_retained_paired_arms_must_share_the_same_original_start():
    fit = valid_fit()
    fit["optimizer_attempts"][1]["start_log_parameters"] = [-2.1, -2.]
    with pytest.raises(ValueError, match="disagree on their original start"):
        opt.validate_optimizer_contract(fit)
