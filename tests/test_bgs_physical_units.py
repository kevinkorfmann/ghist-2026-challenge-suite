import math

import numpy as np
import pytest

from ghist_growth_bgs.simulation import (
    SIMULATION_SCALING_SCHEME,
    growth_demography,
    require_physical_scaling,
)
from ghist_split_migration_bgs.simulation import split_migration_demography


def test_growth_target_units_and_backward_trajectory():
    demography = growth_demography(7000, 9200, 3.7, "exponential")
    population = demography.populations[0]
    assert population.initial_size == pytest.approx(25900)
    assert population.growth_rate == pytest.approx(math.log(3.7) / 9200)
    assert population.initial_size * math.exp(
        -population.growth_rate * 9200
    ) == pytest.approx(7000)
    assert demography.events[0].time == 9200
    assert demography.events[0].initial_size == 7000
    assert growth_demography(7000, 9200, 3.7, "sudden").populations[0].growth_rate == 0


def test_split_target_sizes_time_and_symmetric_backward_migration():
    demography = split_migration_demography(11000, 8000, 19000, 7800, 7e-5)
    assert [population.initial_size for population in demography.populations] == [
        11000,
        8000,
        19000,
    ]
    assert demography.events[0].time == 7800
    assert (
        demography.migration_matrix[1, 2] == demography.migration_matrix[2, 1] == 7e-5
    )
    assert np.diag(demography.migration_matrix).sum() == 0


@pytest.mark.parametrize(
    "records",
    [
        [],
        [{}],
        [{"generator": {"slim_scaling_factor": 10}}],
        [{"simulation_scaling_scheme": "legacy"}],
    ],
)
def test_legacy_or_missing_bgs_scaling_cannot_certify_recovery(records):
    with pytest.raises(ValueError, match="physical-unit"):
        require_physical_scaling(records)


def test_verified_scheme_survives_aggregation_but_mixed_records_fail():
    raw = {"generator": {"simulation_scaling_scheme": SIMULATION_SCALING_SCHEME}}
    aggregate = {"simulation_scaling_scheme": require_physical_scaling([raw])}
    assert require_physical_scaling([raw, aggregate]) == SIMULATION_SCALING_SCHEME
    with pytest.raises(ValueError):
        require_physical_scaling([raw, {}])


@pytest.mark.parametrize("value", [float("nan"), float("inf"), 0, -1])
def test_invalid_physical_demography_rejected(value):
    with pytest.raises(ValueError):
        growth_demography(value, 9200, 3.7, "exponential")
    with pytest.raises(ValueError):
        split_migration_demography(11000, 8000, 19000, 7800, value)
