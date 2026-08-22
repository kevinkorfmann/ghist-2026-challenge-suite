"""Secondary-contact demographic inference and validation."""

from .demography import SecondaryContactParameters, SizeChange, build_demography
from .spectrum import fit_secondary_contact_spectrum, secondary_contact_spectrum
from .model import SecondaryContactNet
from .simulate import SimulationNuisance, parameter_design, simulate_block
from .features import (
    SecondaryFeatureBundle,
    inverse_transform_targets,
    population_pair_hazards,
    summarize_tree_sequence,
    transform_targets,
)

__all__ = [
    "SecondaryContactParameters",
    "SizeChange",
    "SecondaryContactNet",
    "SecondaryFeatureBundle",
    "SimulationNuisance",
    "build_demography",
    "fit_secondary_contact_spectrum",
    "inverse_transform_targets",
    "population_pair_hazards",
    "parameter_design",
    "secondary_contact_spectrum",
    "summarize_tree_sequence",
    "simulate_block",
    "transform_targets",
]
