from __future__ import annotations

import math


def continuous_prediction(fit, parameters):
    """Export the fitted estimates, without substituting preset targets."""
    prediction = {name: float(fit[name]) for name in parameters}
    if any(not math.isfinite(value) or value <= 0 for value in prediction.values()):
        raise ValueError("Fitted estimates must be finite and positive")
    return prediction
