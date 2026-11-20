from __future__ import annotations

from pathlib import Path

import numpy as np


def read_raisd(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="/", ndmin=2)


def separated_peaks(report: np.ndarray, separation: int) -> list[int]:
    report = np.asarray(report, dtype=float)
    if (
        report.ndim != 2
        or report.shape[1] < 7
        or not np.isfinite(report[:, [0, 6]]).all()
        or np.any(report[:, 6] < 0)
    ):
        raise ValueError(
            "RAiSD report needs finite coordinates and nonnegative scores in at least seven columns"
        )
    if separation <= 0:
        raise ValueError("Peak separation must be positive")
    order = np.argsort(-report[:, 6], kind="stable")
    selected: list[int] = []
    # Occupied separation-sized bins bound neighbor lookup to three bins,
    # avoiding an all-selected-peaks scan for each of millions of windows.
    occupied = {}
    for index in order:
        position = float(report[index, 0])
        cell = int(position // separation)
        if any(
            abs(position - other) < separation
            for neighbor in (cell - 1, cell, cell + 1)
            for other in occupied.get(neighbor, ())
        ):
            continue
        selected.append(int(index))
        occupied.setdefault(cell, []).append(position)
    return selected


def nearest_exonic_position(position: int, intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        raise ValueError("No exons found for the selected chromosome")
    best_position = position
    best_distance = float("inf")
    for start, end in intervals:
        candidate = min(max(position, start), end - 1)
        distance = abs(candidate - position)
        if distance < best_distance:
            best_position = candidate
            best_distance = distance
    return best_position
