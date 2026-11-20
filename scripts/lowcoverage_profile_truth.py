"""Development-only demographic objective profile on independent genotype truth.

A composite-likelihood curve is not a calibrated confidence interval: linked
sites require block uncertainty analysis. This diagnostic isolates identifiability
and numerical behavior from genotype-likelihood and realSFS estimation errors.
"""
import argparse
import json
from pathlib import Path

import moments
import numpy as np
from scipy.optimize import minimize_scalar

from ghist_lowcoverage.pipeline import (expected_spectrum, folded_spectrum, physical_parameters,
                                        sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    simulation_path = args.case / 'synthetic/simulation.json'
    simulation = json.loads(simulation_path.read_text())
    if simulation['task'] != 'bottleneck' or simulation['testing_data_access'] or simulation['final_data_access']:
        raise ValueError('independent synthetic bottleneck only')
    spectrum_path = args.case / 'synthetic/truth_folded_sfs.npy'
    values = np.load(spectrum_path).astype(float)
    values[0] = simulation['length'] - values.sum() + values[0]
    data = folded_spectrum(values, [40])
    fractions = np.geomspace(.001, .9999, 25)
    log_times = np.linspace(np.log(.0001), np.log(10), 49)
    rows = []
    for fraction in fractions:
        def objective(log_time):
            model = expected_spectrum([fraction, np.exp(log_time)], [40], 'bottleneck')
            return -float(moments.Inference.ll_multinom(model, data))
        coarse = np.array([objective(t) for t in log_times])
        candidates = [(float(coarse[i]), float(log_times[i])) for i in range(len(log_times))]
        for i in range(1, len(log_times)-1):
            if coarse[i] <= coarse[i-1] and coarse[i] <= coarse[i+1]:
                optimum = minimize_scalar(objective, bounds=(log_times[i-1], log_times[i+1]),
                                          method='bounded', options={'xatol': 1e-6})
                if np.isfinite(optimum.fun):
                    candidates.append((float(optimum.fun), float(optimum.x)))
        loss, log_time = min(candidates)
        parameters = [float(fraction), float(np.exp(log_time))]
        model = expected_spectrum(parameters, [40], 'bottleneck')
        theta = float(moments.Inference.optimal_sfs_scaling(model, data))
        physical = physical_parameters('bottleneck', parameters, theta, simulation['length'])
        rows.append({'fraction': float(fraction), 'scaled_time': parameters[1], 'theta': theta,
                     'negative_log_likelihood': loss, **physical})
    design = simulation['design']
    true_parameters = [design['fraction'], design['generations']/(2*design['ancestral_size'])]
    true_model = expected_spectrum(true_parameters, [40], 'bottleneck')
    true_theta = 4 * simulation['mutation_rate'] * simulation['length'] * design['ancestral_size']
    best = min(row['negative_log_likelihood'] for row in rows)
    for row in rows:
        row['delta_negative_log_likelihood'] = row['negative_log_likelihood'] - best
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'profile.json', {'source_sha256': sha256(Path(__file__)),
        'simulation_sha256': sha256(simulation_path), 'spectrum_sha256': sha256(spectrum_path),
        'task': simulation['task'], 'length': simulation['length'], 'true_parameters': true_parameters,
        'true_parameters_negative_log_likelihood': -float(moments.Inference.ll_multinom(true_model, data)),
        'true_theta': true_theta, 'expected_segregating_sites_at_truth': float(true_model.sum()*true_theta),
        'observed_segregating_sites': float(data.sum()), 'profile': rows,
        'calibrated_confidence_interval': False, 'scientifically_ready': False,
        'testing_data_access': False, 'final_data_access': False})


if __name__ == '__main__':
    main()
