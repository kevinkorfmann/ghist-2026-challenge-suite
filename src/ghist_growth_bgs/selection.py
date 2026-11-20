"""Compare growth histories by predicting entire withheld chromosome arms."""

from __future__ import annotations

import numpy as np

from .spectrum import fit_growth_spectrum, growth_spectrum, sudden_growth_spectrum
from ghist_inference.validation import validate_expected_spectrum


def compare_growth_models(
    spectra,
    exposure,
    groups,
    *,
    mutation_rate=3.36e-8,
    projection=None,
    restarts=12,
    workers=1,
    maxiter=500,
    seed=20261120,
):
    import moments

    spectra = np.asarray(spectra, dtype=float)
    exposure = np.asarray(exposure, dtype=float)
    groups = np.asarray(groups)
    if (
        spectra.ndim != 2
        or exposure.ndim != 1
        or spectra.shape[1] < 5
        or len(spectra) != len(exposure)
        or groups.shape != exposure.shape
        or not np.isfinite(spectra).all()
        or np.any(spectra < 0)
        or not np.isfinite(exposure).all()
        or np.any(exposure < 0)
    ):
        raise ValueError(
            "Expected finite nonnegative block spectra, exposure, and matching groups"
        )
    if np.any((exposure == 0) & (spectra.sum(axis=1) > 0)):
        raise ValueError("A block with observed sites cannot have zero exposure")
    unique = np.unique(groups)
    if len(unique) < 3:
        raise ValueError("At least three independent chromosome groups are required")
    models = {"exponential": growth_spectrum, "sudden": sudden_growth_spectrum}
    folds = []
    for fold, group in enumerate(unique):
        holdout = groups == group
        if exposure[holdout].sum() <= 0 or exposure[~holdout].sum() <= 0:
            raise ValueError(
                "Each training and held-out partition needs positive exposure"
            )
        held_out = moments.Spectrum(spectra[holdout].sum(axis=0))
        if projection is not None:
            held_out = held_out.project([projection])
        if held_out.sum() <= 0:
            raise ValueError("Each held-out group needs polymorphic sites")
        record = {"held_out_group": str(group), "models": {}}
        for index, (name, model) in enumerate(models.items()):
            fit = fit_growth_spectrum(
                spectra[~holdout].sum(axis=0),
                sequence_length=exposure[~holdout].sum(),
                mutation_rate=mutation_rate,
                projection=projection,
                restarts=restarts,
                workers=workers,
                maxiter=maxiter,
                seed=seed + 101 * fold + index,
                convention=name,
            )
            expected = model(
                (fit["post_growth_scaling"], fit["scaled_duration_2N"]),
                held_out.sample_sizes,
            )
            validate_expected_spectrum(expected, held_out)
            score = float(moments.Inference.ll_multinom(expected, held_out))
            if not np.isfinite(score):
                raise RuntimeError("Nonfinite held-out predictive likelihood")
            record["models"][name] = {
                "training_fit": fit,
                "held_out_log_likelihood": score,
            }
        folds.append(record)
    scores = {
        name: sum(f["models"][name]["held_out_log_likelihood"] for f in folds)
        for name in models
    }
    selected = max(scores, key=scores.get)
    wins = sum(
        max(f["models"], key=lambda name: f["models"][name]["held_out_log_likelihood"])
        == selected
        for f in folds
    )
    fits = {
        name: fit_growth_spectrum(
            spectra.sum(axis=0),
            sequence_length=exposure.sum(),
            mutation_rate=mutation_rate,
            projection=projection,
            restarts=restarts,
            workers=workers,
            maxiter=maxiter,
            seed=seed + 10007 + index,
            convention=name,
        )
        for index, name in enumerate(models)
    }
    return {
        "protocol": "Leave one chromosome group out; select by summed conditional SFS-shape likelihood (held-out amplitude profiled); refit all groups",
        "public_testing_answers_used": False,
        "folds": folds,
        "predictive_log_likelihood": scores,
        "selected_convention": selected,
        "selected_wins": wins,
        "groups": len(unique),
        "selection_consistent_across_groups": wins == len(unique),
        "full_data_fits": fits,
        "prediction": {
            key: fits[selected][key] for key in ("generations", "post_growth_scaling")
        },
        "scope": "Model comparison and candidate estimate; independent recovery under BGS still required",
    }
