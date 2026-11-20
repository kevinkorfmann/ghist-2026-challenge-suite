from __future__ import annotations

import math


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
