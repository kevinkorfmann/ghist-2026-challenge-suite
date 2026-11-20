import json
from pathlib import Path
import msprime
import moments
import numpy as np
from ghist_secondary_contact.demography import (
    SecondaryContactParameters,
    build_demography,
)
from ghist_secondary_contact.simulate import SimulationNuisance, nuisance_size_changes
from ghist_secondary_contact.spectrum import (
    generator_matched_secondary_contact_spectrum,
)
from ghist_inference.validation import validate_expected_spectrum

rows = []
for i, contact in enumerate((2700, 9000)):
    p = SecondaryContactParameters(12000, 18000, 7000, 10500, contact, 6e-5)
    nuisance = SimulationNuisance(
        ancestral_epoch_factors=(0.5, 1.5, 2),
        population1_epoch_factors=(0.8, 0.5),
        population2_epoch_factor=1.5,
        population2_contact_factor=1 / 3,
    )
    demography = build_demography(p, size_changes=nuisance_size_changes(p, nuisance))
    total = np.zeros((13, 11))
    for ts in msprime.sim_ancestry(
        samples={"population1": 6, "population2": 5},
        demography=demography,
        sequence_length=1,
        num_replicates=10000,
        random_seed=202611230 + i,
    ):
        total += ts.allele_frequency_spectrum(
            sample_sets=[ts.samples(population=1), ts.samples(population=2)],
            mode="branch",
            polarised=True,
            span_normalise=False,
        )
    empirical = moments.Spectrum(total / 10000 / (4 * p.ancestral_size))
    model = generator_matched_secondary_contact_spectrum(
        (
            p.population1_size / p.ancestral_size,
            p.population2_size / p.ancestral_size,
            p.split_time / (2 * p.ancestral_size),
            p.contact_time / p.split_time,
            2 * p.ancestral_size * p.migration_rate,
        ),
        [12, 10],
    )
    validate_expected_spectrum(model, empirical)
    a = np.ma.filled(empirical, 0)
    b = np.ma.filled(model, 0)
    rows.append(
        {
            "contact_time": contact,
            "replicates": 10000,
            "seed": 202611230 + i,
            "simulation": "msprime independent unlinked genealogies, branch SFS, diploid ploidy",
            "relative_total_branch_length_error": float(a.sum() / b.sum() - 1),
            "normalized_L1_difference": float(np.abs(a / a.sum() - b / b.sum()).sum()),
            "empirical_fixed_difference_mass": float(a[0, -1] + a[-1, 0]),
            "model_fixed_difference_mass": float(b[0, -1] + b[-1, 0]),
        }
    )
    print(rows[-1], flush=True)
Path(
    "experiments/scientific_2026_11_20/demography/hidden_epoch_consistency.json"
).write_text(json.dumps(rows, indent=2) + "\n")
