"""Bounded numerical refinement without moments' non-progressing retry loop."""

from __future__ import annotations

import numpy as np

from .validation import validate_expected_spectrum


def refine_spectrum(evaluate, *, dt_fac, max_refinements=6, relative_tolerance=1e-5):
    """Require a valid spectrum and agreement at two adjacent time resolutions.

    ``evaluate`` must perform a *non-adaptive* complete model integration. The
    moments 1.7.1 adaptive migration integrator can restore the same negative
    state forever once its timestep floor is reached. Whole-model retries avoid
    that loop, preserve every likelihood bin, and fail explicitly when spatial
    approximation errors persist despite timestep refinement.
    """
    if not np.isfinite(dt_fac) or dt_fac <= 0:
        raise ValueError("integration timestep must be finite and positive")
    if (
        isinstance(max_refinements, bool)
        or not isinstance(max_refinements, (int, np.integer))
        or not 1 <= max_refinements <= 16
    ):
        raise ValueError("max_refinements must be an integer between 1 and 16")
    if not np.isfinite(relative_tolerance) or relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be finite and positive")
    previous = None
    last_reason = "no converged refinement"
    for refinement in range(max_refinements + 1):
        try:
            expected = evaluate(dt_fac / 2**refinement)
        except (ArithmeticError, np.linalg.LinAlgError) as error:
            previous = None
            last_reason = f"{type(error).__name__}: {error}"
            continue
        domain = np.ma.array(
            np.zeros(expected.shape), mask=np.ma.getmaskarray(expected)
        )
        try:
            validate_expected_spectrum(expected, domain)
        except ValueError as error:
            previous = None
            last_reason = str(error)
            continue
        values = np.asarray(expected)[~np.ma.getmaskarray(domain)]
        if previous is not None:
            relative_change = np.abs(values - previous).sum() / max(
                np.abs(values).sum(), np.finfo(float).tiny
            )
            if relative_change <= relative_tolerance:
                return expected
            last_reason = f"relative spectrum change {relative_change:.6g}"
        previous = values.copy()
    raise ValueError(
        f"Spectrum integration failed after {max_refinements + 1} finite "
        f"refinements: {last_reason}"
    )
