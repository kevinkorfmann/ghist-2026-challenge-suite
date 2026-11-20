import json

import moments
import numpy as np
import pytest

from ghist_inference.checkpoint import (
    load_spectrum_checkpoint,
    save_spectrum_checkpoint,
)
from ghist_inference.integration import refine_spectrum
from scripts import validate_demography_independently as benchmark


def test_seed_streams_remain_unique_beyond_one_hundred_replicates():
    seeds = [
        seed
        for index in range(11)
        for replicate in range(150)
        for seed in benchmark.simulation_seeds(index, replicate).values()
    ]
    assert len(seeds) == len(set(seeds))
    assert benchmark.simulation_seeds(0, 100) != benchmark.simulation_seeds(1, 0)
    assert benchmark.simulation_seeds(0, 0, "legacy_v1")["ancestry"] == 202611210
    with pytest.raises(ValueError, match="collides"):
        benchmark.simulation_seeds(0, 100, "legacy_v1")
    with pytest.raises(ValueError, match="uint32"):
        benchmark.simulation_seeds(0, 2**32)


def test_refinement_is_finite_and_never_clips_invalid_bins():
    calls = []

    def evaluate(step):
        calls.append(step)
        return moments.Spectrum(np.array([0, -1e-9, 1, 0]))

    with pytest.raises(ValueError, match="7 finite refinements"):
        refine_spectrum(evaluate, dt_fac=0.02)
    assert calls == [0.02 / 2**i for i in range(7)]


def test_refinement_requires_two_valid_agreeing_resolutions():
    calls = []

    def evaluate(step):
        calls.append(step)
        return moments.Spectrum(np.array([0, 2, 1, 0]))

    expected = refine_spectrum(evaluate, dt_fac=0.02)
    assert calls == [0.02, 0.01]
    assert expected[1] == 2


def test_real_nonprogress_proposal_is_rejected_without_native_adaptive_steps(
    monkeypatch,
):
    from ghist_split_migration_bgs.spectrum import isolation_with_migration_spectrum

    original = moments.Spectrum.integrate

    def integrate(self, *args, **kwargs):
        assert kwargs.get("adapt_dt") is False
        return original(self, *args, **kwargs)

    monkeypatch.setattr(moments.Spectrum, "integrate", integrate)
    with pytest.raises(ValueError, match="finite refinements"):
        isolation_with_migration_spectrum(
            [
                37.002346797132866,
                6.218138603544293,
                0.059115084943017114,
                22.980411155782058,
            ],
            (20, 16),
            adapt_dt=True,
        )


def test_spectrum_checkpoint_rejects_changed_configuration_and_content(tmp_path):
    path = tmp_path / "spectrum.npz"
    spectrum = np.array([0, 5, 2, 0])
    metadata = save_spectrum_checkpoint(path, spectrum, {"seed": 1})
    loaded, _ = load_spectrum_checkpoint(path, expected={"seed": 1})
    np.testing.assert_array_equal(loaded, spectrum)
    with pytest.raises(ValueError, match="configuration mismatch"):
        load_spectrum_checkpoint(path, expected={"seed": 2})
    np.savez(path, spectrum=spectrum + 1, metadata=json.dumps(metadata))
    with pytest.raises(ValueError, match="digest mismatch"):
        load_spectrum_checkpoint(path, expected={"seed": 1})


@pytest.mark.parametrize("failure", [KeyboardInterrupt, ValueError])
def test_resume_preserves_simulation_and_first_model_after_interruption(
    tmp_path, monkeypatch, failure
):
    config = next(row for row in benchmark.configurations() if row["task"] == 3)
    simulations, fits = [], []

    def simulate(*args):
        simulations.append(1)
        return np.array([0, 5, 2, 0]), {
            "truth": {
                "ancestral_effective_size": 7000,
                "generations": 9200,
                "post_growth_scaling": 3.7,
            },
            "sites": 7,
            "length": 1000,
            "mutation_rate": 3.36e-8,
            "recombination_rate": 2e-8,
        }

    def fit(spectrum, **settings):
        fits.append(settings["convention"])
        if fits == ["exponential", "sudden"]:
            raise failure("simulate process interruption")
        return {
            "ancestral_effective_size": 7000,
            "generations": 9200,
            "post_growth_scaling": 3.7,
        }

    monkeypatch.setattr(benchmark, "simulate_configuration", simulate)
    monkeypatch.setattr(benchmark, "fit_growth_spectrum", fit)
    item = (config, 0, tmp_path, 1000, 2, "counter_v2")
    if failure is KeyboardInterrupt:
        with pytest.raises(KeyboardInterrupt):
            benchmark.run_one(item)
    else:
        first = json.loads(open(benchmark.run_one(item)).read())
        assert first["status"] == "partial"
    result = json.loads(open(benchmark.run_one(item)).read())
    assert simulations == [1]
    assert fits == ["exponential", "sudden", "sudden"]
    assert result["status"] == "completed"
    assert "1000 bp" in result["scope"]
    if failure is ValueError:
        checkpoint = json.loads(next(tmp_path.glob("*.fit1.json")).read_text())
        assert "ValueError" in checkpoint["previous_failed_attempts"][0]["error"]
    with pytest.raises(ValueError, match="Fit checkpoint configuration mismatch"):
        benchmark.run_one((config, 0, tmp_path, 1000, 3, "counter_v2"))


def test_tree_checkpoint_reuses_only_matching_verified_simulation(tmp_path):
    from types import SimpleNamespace
    import msprime
    from ghist_inference.checkpoint import simulation_tree_checkpoint, file_digest

    calls = []

    def simulate():
        calls.append(1)
        return msprime.sim_ancestry(2, sequence_length=100, random_seed=1)

    args = SimpleNamespace(
        out=tmp_path / "run.json",
        trees=tmp_path / "run.trees",
        resume=False,
        existing_trees_sha256=None,
    )
    simulation_tree_checkpoint(args, {"seed": 1}, simulate)
    with pytest.raises(FileExistsError, match="--resume"):
        simulation_tree_checkpoint(args, {"seed": 1}, simulate)
    args.resume = True
    recovered = simulation_tree_checkpoint(args, {"seed": 1}, simulate)
    assert recovered.num_samples == 4
    assert calls == [1]
    with pytest.raises(ValueError, match="configuration or tree digest mismatch"):
        simulation_tree_checkpoint(args, {"seed": 2}, simulate)
    args.out = tmp_path / "imported.json"
    with pytest.raises(ValueError, match="recorded SHA-256"):
        simulation_tree_checkpoint(args, {"seed": 1}, simulate)
    args.existing_trees_sha256 = file_digest(args.trees)
    simulation_tree_checkpoint(args, {"seed": 1}, simulate)
    assert json.loads(args.out.with_suffix(".simulation.json").read_text())[
        "imported_existing_tree"
    ]


@pytest.mark.parametrize(
    "settings",
    [
        {"max_refinements": 0},
        {"max_refinements": 1.5},
        {"max_refinements": 100000},
        {"max_refinements": True},
        {"relative_tolerance": np.nan},
        {"relative_tolerance": 0},
        {"relative_tolerance": np.inf},
    ],
)
def test_refinement_rejects_invalid_controls(settings):
    with pytest.raises(ValueError):
        refine_spectrum(lambda step: None, dt_fac=0.02, **settings)


def test_refinement_retries_arithmetic_failures_but_not_parameter_errors():
    calls = []

    def evaluate(step):
        calls.append(step)
        if len(calls) == 1:
            raise FloatingPointError("coarse numerical failure")
        return moments.Spectrum(np.array([0, 2, 1, 0]))

    refine_spectrum(evaluate, dt_fac=0.02)
    assert calls == [0.02, 0.01, 0.005]

    def invalid(step):
        raise ValueError("bad demographic parameter")

    with pytest.raises(ValueError, match="bad demographic parameter"):
        refine_spectrum(invalid, dt_fac=0.02)


@pytest.mark.parametrize("task", [2, 4])
def test_exact_bound_bootstrap_start_is_valid_but_boundary_is_reported(
    monkeypatch, task
):
    if task == 2:
        import ghist_secondary_contact.spectrum as module

        model_name = "secondary_contact_spectrum"
        fit = module.fit_secondary_contact_spectrum
        initial = (0.01, 1.0, 0.3, 0.1, 1.0)
    else:
        import ghist_split_migration_bgs.spectrum as module

        model_name = "isolation_with_migration_spectrum"
        fit = module.fit_isolation_with_migration_spectrum
        initial = (0.01, 1.0, 0.3, 1.0)
    model = moments.Spectrum(np.ones((5, 5)))
    monkeypatch.setattr(module, model_name, lambda *args, **kwargs: model)

    def optimizer(start, *args, **kwargs):
        assert kwargs["full_output"] is True
        return (
            np.array(start),
            0.0,
            {"warnflag": 1, "nit": 1, "funcalls": 2, "task": "iteration limit"},
        )

    monkeypatch.setattr(moments.Inference, "optimize_log_lbfgsb", optimizer)
    result = fit(
        np.asarray(model),
        mutation_rate=1e-8,
        sequence_length=10000,
        initial_scaled=initial,
        restarts=1,
    )
    assert result["scaled_parameters"][next(iter(result["scaled_parameters"]))] == 0.01
    assert result["best_restart_converged"] is False
    assert result["converged_restarts"] == 0
    assert 0 in result["optimizer_diagnostics"]["near_boundary_parameter_indices"]


@pytest.mark.parametrize("flag", [False, None, 1, "true"])
def test_freeze_rejects_finite_fit_without_actual_convergence(flag):
    from ghist_inference.freeze import require_optimizer_convergence

    with pytest.raises(ValueError, match="unconverged optimizer"):
        require_optimizer_convergence(
            {
                "generations": 1000.0,
                "multinomial_log_likelihood": -5.0,
                "best_restart_converged": flag,
            }
        )


def test_freeze_requires_current_optimizer_diagnostics():
    from ghist_inference.freeze import (
        require_optimizer_convergence,
        require_projection_convergence,
    )

    with pytest.raises(ValueError, match="regenerate the fit"):
        require_optimizer_convergence({"generations": 1000.0})
    with pytest.raises(ValueError, match="individual projection fits"):
        require_projection_convergence({})
    require_optimizer_convergence({"best_restart_converged": True})


@pytest.mark.parametrize("task", [2, 3, 4])
@pytest.mark.parametrize("convergence", [False, None])
def test_all_demography_freeze_paths_reject_unconverged_or_legacy_candidates(
    tmp_path, monkeypatch, task, convergence
):
    import sys
    from scripts import freeze_secondary_contact, freeze_growth_bgs
    from scripts import freeze_split_migration_bgs

    candidate = {"multinomial_log_likelihood": -1.0, "generations": 1000.0}
    if convergence is not None:
        candidate["best_restart_converged"] = convergence
    if task == 2:
        module = freeze_secondary_contact
        flags = [
            "raw-fit",
            "standard-fit",
            "calibration",
            "resolution",
            "broad-resolution",
            "varying-resolution",
            "stress-comparison",
            "stress-transfer",
            "bootstrap",
            "vcf-audit",
        ]
        sources = {flag: {} for flag in flags}
        sources["raw-fit"] = candidate
        monkeypatch.setattr(module, "require_spectrum_contract", lambda _: None)
    elif task == 3:
        module = freeze_growth_bgs
        flags = [
            "audit",
            "middle-fit",
            "bootstrap",
            "convention",
            "model-comparison",
            "smoke",
            "grid",
            "calibration",
            "stress",
            "final-sealed-sha256",
        ]
        sources = {flag: {} for flag in flags}
        sources["middle-fit"] = {
            "strata": {"noncoding_1kb_100kb": {
                "by_projection": {str(p): candidate for p in (84, 60, 40)}
            }}
        }
        monkeypatch.setattr(module, "require_physical_scaling", lambda _: None)
    else:
        module = freeze_split_migration_bgs
        flags = [
            "baseline",
            "refined",
            "bootstrap",
            "calibration",
            "near-recovery",
            "testing-vcf-sha256",
            "coding-bed-sha256",
            "recombination-map-sha256",
            "final-sealed-sha256",
        ]
        sources = {flag: {} for flag in flags}
        sources["baseline"] = {
            "by_projection": {p: candidate for p in ("50x40", "70x50", "94x68")}
        }
        monkeypatch.setattr(module, "require_physical_scaling", lambda _: None)
    monkeypatch.setattr(module, "load", lambda path: sources[path.name])
    args = ["freeze", "--out", str(tmp_path / "frozen.json")]
    for flag in flags:
        args += ["--" + flag, flag]
    monkeypatch.setattr(sys, "argv", args)
    with pytest.raises(ValueError, match="unconverged optimizer|regenerate the fit"):
        module.main()
    assert not (tmp_path / "frozen.json").exists()


def test_growth_summary_preserves_projection_convergence(tmp_path, monkeypatch):
    import sys
    from scripts import summarize_growth_bgs_fits as summarize

    fit = {
        "generations": 1000,
        "post_growth_scaling": 2,
        "ancestral_effective_size": 10000,
        "modern_effective_size": 20000,
        "best_restart_converged": False,
        "converged_restarts": 0,
    }
    source = tmp_path / "fit.json"
    source.write_text(json.dumps({"projection": 40, "results": {"noncoding": fit}}))
    target = tmp_path / "summary.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["summary", "--fit", str(source), "--job", "1", "--out", str(target)],
    )
    summarize.main()
    summary = json.loads(target.read_text())
    assert summary["strata"]["noncoding"]["by_projection"]["40"] == fit
