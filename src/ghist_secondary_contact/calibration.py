from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LogCorrection:
    center: np.ndarray
    scale: np.ndarray
    coefficients: np.ndarray

    def predict(self, estimates: np.ndarray) -> np.ndarray:
        logged = np.log(estimates)
        features = np.column_stack(
            (np.ones(len(estimates)), (logged - self.center) / self.scale)
        )
        return np.exp(logged + features @ self.coefficients)


def competition_scores(predictions: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predictions = np.asarray(predictions, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    if predictions.shape != truth.shape or predictions.ndim != 2:
        raise ValueError("predictions and truth must be matching matrices")
    if np.any(predictions <= 0) or np.any(truth <= 0):
        raise ValueError("predictions and truth must be positive")
    return np.sqrt(np.sum((predictions / truth - 1.0) ** 2, axis=1))


def fit_log_correction(
    method: str,
    estimates: np.ndarray,
    truth: np.ndarray,
) -> LogCorrection | None:
    estimates = np.asarray(estimates, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    if estimates.shape != truth.shape or estimates.ndim != 2:
        raise ValueError("estimates and truth must be matching matrices")
    if np.any(estimates <= 0) or np.any(truth <= 0):
        raise ValueError("estimates and truth must be positive")
    if method == "identity":
        return None
    logged = np.log(estimates)
    residual = np.log(truth) - logged
    center = logged.mean(axis=0)
    scale = logged.std(axis=0)
    scale[scale < 1e-8] = 1.0
    if method == "median_ratio":
        coefficients = np.zeros((estimates.shape[1] + 1, estimates.shape[1]))
        coefficients[0] = np.median(residual, axis=0)
        return LogCorrection(center, scale, coefficients)
    prefix = "ridge_"
    if not method.startswith(prefix):
        raise ValueError(f"unknown correction {method!r}")
    penalty_strength = float(method[len(prefix) :])
    if penalty_strength < 0:
        raise ValueError("ridge penalty must be nonnegative")
    features = np.column_stack((np.ones(len(estimates)), (logged - center) / scale))
    penalty = np.eye(features.shape[1]) * penalty_strength
    penalty[0, 0] = 0.0
    coefficients = np.linalg.solve(
        features.T @ features + penalty,
        features.T @ residual,
    )
    return LogCorrection(center, scale, coefficients)


def apply_log_correction(
    correction: LogCorrection | None,
    estimates: np.ndarray,
) -> np.ndarray:
    estimates = np.asarray(estimates, dtype=np.float64)
    return estimates if correction is None else correction.predict(estimates)
