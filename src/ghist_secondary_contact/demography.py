from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import msprime


@dataclass(frozen=True)
class SecondaryContactParameters:
    ancestral_size: float
    population1_size: float
    population2_size: float
    split_time: float
    contact_time: float
    migration_rate: float

    def validate(self) -> None:
        sizes = (self.ancestral_size, self.population1_size, self.population2_size)
        if any(size <= 0 for size in sizes):
            raise ValueError("population sizes must be positive")
        if not 0 < self.contact_time < self.split_time:
            raise ValueError("contact time must be positive and younger than split time")
        if self.migration_rate < 0:
            raise ValueError("migration rate cannot be negative")


@dataclass(frozen=True)
class SizeChange:
    population: str
    time: float
    size: float


def build_demography(
    parameters: SecondaryContactParameters,
    *,
    size_changes: Iterable[SizeChange] = (),
    migration_asymmetry: float = 1.0,
) -> msprime.Demography:
    """Build isolation followed by recent bidirectional migration.

    ``migration_asymmetry`` is the forward-rate ratio between directions while
    preserving their geometric mean. A value of one gives symmetric migration.
    Additional size changes support generator-misspecification stress tests.
    """
    parameters.validate()
    if migration_asymmetry <= 0:
        raise ValueError("migration asymmetry must be positive")
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=parameters.ancestral_size)
    demography.add_population(name="population1", initial_size=parameters.population1_size)
    demography.add_population(name="population2", initial_size=parameters.population2_size)
    root = migration_asymmetry**0.5
    demography.set_migration_rate(
        source="population1",
        dest="population2",
        rate=parameters.migration_rate * root,
    )
    demography.set_migration_rate(
        source="population2",
        dest="population1",
        rate=parameters.migration_rate / root,
    )
    demography.add_migration_rate_change(
        time=parameters.contact_time,
        source="population1",
        dest="population2",
        rate=0,
    )
    demography.add_migration_rate_change(
        time=parameters.contact_time,
        source="population2",
        dest="population1",
        rate=0,
    )
    allowed = {"ancestral", "population1", "population2"}
    for change in size_changes:
        if change.population not in allowed:
            raise ValueError(f"unknown population {change.population!r}")
        if change.time <= 0 or change.size <= 0:
            raise ValueError("size-change times and sizes must be positive")
        if change.population == "ancestral" and change.time <= parameters.split_time:
            raise ValueError("ancestral size changes must predate the split")
        if change.population != "ancestral" and change.time >= parameters.split_time:
            raise ValueError("derived-population size changes must postdate the split")
        demography.add_population_parameters_change(
            time=change.time,
            initial_size=change.size,
            population=change.population,
        )
    demography.add_population_split(
        time=parameters.split_time,
        derived=["population1", "population2"],
        ancestral="ancestral",
    )
    demography.sort_events()
    return demography
