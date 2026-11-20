from __future__ import annotations

import math
import statistics


def require_optimizer_convergence(fit, label="candidate fit"):
    """Require numerical termination evidence, not scientific adequacy."""
    if "best_restart_converged" not in fit:
        raise ValueError(
            f"{label} lacks optimizer convergence diagnostics; regenerate the fit "
            "with the current fitter before freezing (preserve the legacy artifact)"
        )
    if fit["best_restart_converged"] is not True:
        raise ValueError(
            f"{label} selected an unconverged optimizer restart; rerun and resolve "
            "the convergence failure before freezing"
        )


def require_projection_convergence(
    fits, label="candidate projections", *, expected_projections=None
):
    if not isinstance(fits, dict) or not fits:
        raise ValueError(
            f"{label} are missing; regenerate individual projection fits with "
            "optimizer convergence diagnostics before freezing"
        )
    if expected_projections is not None:
        expected = {str(projection) for projection in expected_projections}
        actual = {str(projection) for projection in fits}
        if actual != expected or len(actual) != len(fits):
            raise ValueError(
                f"{label} do not match the required projection set; "
                f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
            )
    for projection, fit in fits.items():
        require_optimizer_convergence(fit, f"{label} {projection}")


def continuous_prediction(fit, parameters, *, nonnegative=()):
    """Export the fitted estimates, without substituting preset targets."""
    if not set(nonnegative).issubset(parameters):
        raise ValueError("Nonnegative parameters must belong to the exported prediction")
    if any(isinstance(fit[name], bool) for name in parameters):
        raise ValueError("Fitted estimates must be numeric, not boolean")
    prediction = {name: float(fit[name]) for name in parameters}
    if any(
        not math.isfinite(value)
        or (value < 0 if name in nonnegative else value <= 0)
        for name, value in prediction.items()
    ):
        raise ValueError("Fitted estimates must be finite and positive, except permitted zero boundaries")
    return prediction


def verified_projection_median(fits, stored, parameters, *, label, expected_projections):
    """Bind an exported median to every required converged projection fit."""
    require_projection_convergence(
        fits, label, expected_projections=expected_projections
    )
    if not isinstance(stored, dict) or not set(parameters).issubset(stored):
        raise ValueError(f"{label} stored median is missing required parameters")
    values = [continuous_prediction(fit, parameters) for fit in fits.values()]
    observed = continuous_prediction(stored, parameters)
    recomputed = {
        parameter: float(statistics.median(row[parameter] for row in values))
        for parameter in parameters
    }
    for parameter in parameters:
        if not math.isclose(observed[parameter], recomputed[parameter], rel_tol=1e-12):
            raise ValueError(
                f"{label} stored median differs from contributing fits: {parameter}; "
                "regenerate the summary before freezing"
            )
    return recomputed


def validated_intervals(intervals, parameters, *, nonnegative=(), minimum_exclusive=None):
    """Validate interval shape and physical domain, without certifying coverage."""
    if set(intervals) != set(parameters):
        raise ValueError("Interval parameters do not match the prediction contract")
    minimum_exclusive = minimum_exclusive or {}
    if not (set(nonnegative) | set(minimum_exclusive)).issubset(parameters):
        raise ValueError("Interval bounds reference unknown parameters")
    result = {}
    for parameter in parameters:
        bounds = intervals[parameter]
        if (
            not isinstance(bounds, (list, tuple))
            or len(bounds) != 2
            or any(isinstance(value, bool) for value in bounds)
        ):
            raise ValueError(f"{parameter} interval must contain two numeric bounds")
        try:
            lower, upper = (float(value) for value in bounds)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{parameter} interval bounds must be numeric") from error
        minimum = minimum_exclusive.get(parameter, 0.0)
        if (
            not all(math.isfinite(value) for value in (lower, upper))
            or lower > upper
            or (lower < 0 if parameter in nonnegative else lower <= minimum)
        ):
            raise ValueError(
                f"{parameter} interval must be finite, ordered, and within its physical domain"
            )
        result[parameter] = [lower, upper]
    return result


def validated_quantile_intervals(records, parameters, **domain):
    """Read the existing 2.5/97.5-percentile schema with explicit required fields."""
    intervals = {}
    for parameter in parameters:
        try:
            quantiles = records[parameter]["quantiles"]
            intervals[parameter] = [quantiles["0.025"], quantiles["0.975"]]
        except (KeyError, TypeError) as error:
            raise ValueError(
                f"{parameter} interval lacks required 0.025/0.975 quantiles"
            ) from error
    return validated_intervals(intervals, parameters, **domain)


def validate_physical_identity(fit, task, *, selected=None):
    """Check retained theta and scaled-to-physical arithmetic, not model adequacy.

    The fitters omit each restart's parameter array from their public diagnostics.
    These checks therefore cannot establish full selected-parameter provenance.
    """
    if task not in (3, 4):
        raise ValueError('Physical arithmetic requires task 3 or 4')
    # Older arithmetic-complete fits may omit restart rows. If rows are supplied,
    # their selected theta is evidence and must agree with the fitted theta.
    if selected is None and 'restart_diagnostics' in fit:
        diagnostics = fit['restart_diagnostics']
        best_id = fit.get('best_restart')
        if (not isinstance(diagnostics, list) or isinstance(best_id, bool)
                or not isinstance(best_id, int)
                or any(not isinstance(row, dict) for row in diagnostics)):
            raise ValueError('Fitted selected-restart theta metadata is malformed')
        candidates = [row for row in diagnostics if row.get('restart') == best_id]
        if len(candidates) != 1:
            raise ValueError('Fitted selected-restart theta metadata is not unique')
        selected = candidates[0]

    def positive(record, names):
        try:
            return continuous_prediction(record, names)
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError('Fitted arithmetic metadata must be complete, finite and positive') from error

    common = positive(fit, ('theta', 'mutation_rate', 'sequence_length'))
    if selected is not None and common['theta'] != positive(selected, ('theta',))['theta']:
        raise ValueError('Fitted theta differs from the selected retained restart')
    denominator = 4 * common['mutation_rate'] * common['sequence_length']
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError('Fitted theta exposure scaling is nonfinite or nonpositive')
    ancestral = common['theta'] / denominator

    def check(name, expected):
        observed = positive(fit, (name,))[name]
        # Arithmetic/serialization tolerance; this is not a recovery threshold.
        if not math.isfinite(expected) or not math.isclose(observed, expected, rel_tol=1e-12, abs_tol=0):
            raise ValueError(f'Fitted {name} disagrees with retained scaled-to-physical arithmetic')

    if task == 3:
        scaled = positive(fit, ('scaled_duration_2N', 'post_growth_scaling'))
        check('ancestral_effective_size', ancestral)
        check('generations', scaled['scaled_duration_2N'] * 2 * ancestral)
        check('modern_effective_size', scaled['post_growth_scaling'] * ancestral)
    else:
        scaled = positive(
            fit.get('scaled_parameters'),
            ('AFR_size_ratio', 'EUR_size_ratio', 'split_duration_2N', 'migration_2N'),
        )
        check('ancestral_population_size', ancestral)
        check('AFR_population_size', scaled['AFR_size_ratio'] * ancestral)
        check('EUR_population_size', scaled['EUR_size_ratio'] * ancestral)
        check('generations_since_split', scaled['split_duration_2N'] * 2 * ancestral)
        if not math.isfinite(ancestral) or ancestral <= 0:
            raise ValueError('Fitted ancestral scaling is nonfinite or nonpositive')
        check('migration_rate', scaled['migration_2N'] / (2 * ancestral))
