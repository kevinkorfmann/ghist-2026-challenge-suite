"""Target demographic units; stdpopsim owns all computational rescaling."""

import math

import msprime


def split_migration_demography(
    ancestral_size, afr_size, eur_size, generations, migration_rate
):
    if not all(
        math.isfinite(x) and x > 0
        for x in (ancestral_size, afr_size, eur_size, generations, migration_rate)
    ):
        raise ValueError("finite positive demographic parameters required")
    demography = msprime.Demography()
    for name, size in zip(("ANC", "AFR", "EUR"), (ancestral_size, afr_size, eur_size)):
        demography.add_population(name=name, initial_size=size)
    demography.set_symmetric_migration_rate(["AFR", "EUR"], rate=migration_rate)
    demography.add_population_split(
        time=generations, derived=["AFR", "EUR"], ancestral="ANC"
    )
    return demography
