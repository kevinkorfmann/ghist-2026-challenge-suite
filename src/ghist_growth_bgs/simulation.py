"""Physical-unit demographic models for stdpopsim's SLiM engine.

The engine alone applies its computational Q scaling. Pre-multiplying N and
time by Q while dividing mutation by Q changes Nr and Ns if map/DFE rates are
left unchanged, so those older simulations are not matched target-model recovery.
"""

import math

import msprime


SIMULATION_SCALING_SCHEME = "physical_units_stdpopsim_q_v1"


def growth_demography(ancestral_size, generations, ratio, convention):
    if (
        not all(
            math.isfinite(x) and x > 0 for x in (ancestral_size, generations, ratio)
        )
        or ratio <= 1
    ):
        raise ValueError("finite positive sizes/time and growth ratio > 1 required")
    if convention not in {"exponential", "sudden"}:
        raise ValueError("unknown growth convention")
    demography = msprime.Demography()
    demography.add_population(
        name="population",
        initial_size=ancestral_size * ratio,
        growth_rate=math.log(ratio) / generations if convention == "exponential" else 0,
    )
    demography.add_population_parameters_change(
        time=generations,
        initial_size=ancestral_size,
        growth_rate=0,
        population="population",
    )
    return demography


def require_physical_scaling(records):
    """Do not allow historical differently scaled BGS data to certify a freeze."""
    if not records or any(
        record.get(
            "simulation_scaling_scheme",
            record.get("generator", {}).get("simulation_scaling_scheme"),
        )
        != SIMULATION_SCALING_SCHEME
        for record in records
    ):
        raise ValueError(
            "BGS recovery requires regenerated physical-unit simulations; legacy scaling is not matched validation"
        )
    return SIMULATION_SCALING_SCHEME
