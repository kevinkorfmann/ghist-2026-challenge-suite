"""Evaluate four frozen existing-case points; never optimize or simulate."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np

from ghist_bottleneck.dadi_reference import physical_parameters
from ghist_inference.checkpoint import atomic_json
from ghist_inference.validation import validate_expected_spectrum, validate_spectrum

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest_path = HERE / 'manifest.json'
    manifest_hash = digest(manifest_path)
    config = json.loads(manifest_path.read_text())
    for relative, expected in config['frozen_sha256'].items():
        if digest(HERE / relative) != expected:
            raise ValueError(f'Frozen input/source changed: {relative}')
    versions = {name: importlib.metadata.version(name) for name in config['runtime_versions']}
    if versions != config['runtime_versions']:
        raise ValueError('Must use identical archived runtime versions')
    output = HERE / 'results'
    output.mkdir(exist_ok=False)
    with np.load(HERE / config['spectrum'], allow_pickle=False) as loaded:
        spectrum = validate_spectrum(loaded['spectrum'], 1)
    record = {'schema': 'bottleneck_stress_fixed_points_v1', 'manifest_sha256': manifest_hash,
              'state': 'RUNNING', 'planned_points': 4, 'completed_points': [],
              'versions': versions, 'truth_used_posthoc': True,
              'scope': config['scope'], 'optimized': False, 'simulated': False}
    atomic_json(output / 'evaluations.json', record)
    arrays = {'observed_spectrum_unfolded_storage': spectrum}
    started = time.monotonic()
    try:
        for point in config['points']:
            engine = point['engine']
            if engine == 'moments':
                import moments as package
                expected = package.Demographics1D.two_epoch(point['parameters'], [len(spectrum) - 1])
            else:
                import dadi as package
                extrapolated = package.Numerics.make_extrap_log_func(package.Demographics1D.two_epoch)
                expected = extrapolated(point['parameters'], [len(spectrum) - 1], config['dadi_grid_points'])
            data = package.Spectrum(spectrum)
            if config['folded']:
                data = data.fold()
                expected = expected.fold()
            expected = validate_expected_spectrum(expected, data)
            likelihood = float(package.Inference.ll_multinom(expected, data))
            theta = float(package.Inference.optimal_sfs_scaling(expected, data))
            if not np.isfinite([likelihood, theta]).all() or theta <= 0:
                raise ValueError('Invalid fixed-point likelihood/scaling')
            fraction, scaled_time = point['parameters']
            generations, fraction, ancestral_ne = physical_parameters(
                fraction, scaled_time, theta,
                mutation_rate=config['mutation_rate'], sequence_length=config['sequence_length'])
            row = {**point, 'conditional_log_likelihood': likelihood, 'negative_objective': -likelihood,
                   'optimal_theta': theta, 'profiled_ancestral_effective_size': ancestral_ne,
                   'profiled_generations': generations,
                   'profiled_current_effective_size': fraction * ancestral_ne,
                   'active_bins': int(np.count_nonzero(~np.ma.getmaskarray(data))),
                   'active_observed_sum': float(data.sum()),
                   'active_expected_unit_theta_sum': float(expected.sum()),
                   'physical_roundtrip_scaled_time': generations / (2 * ancestral_ne),
                   'physical_roundtrip_theta': 4 * config['mutation_rate'] * config['sequence_length'] * ancestral_ne}
            if point['label'] == 'selected':
                archived = json.loads((HERE / point['fit_record']).read_text())['fit']
                row['archived_ll_difference'] = likelihood - archived['multinomial_log_likelihood']
                row['archived_theta_relative_difference'] = theta / archived['theta'] - 1
                row['archived_physical_differences'] = {
                    'generations': generations - archived['generations'],
                    'ancestral_effective_size': ancestral_ne - archived['ancestral_effective_size'],
                    'fraction': fraction - archived['post_decline_fraction']}
            record['completed_points'].append(row)
            key = engine + '_' + point['label']
            arrays[key + '_expected_per_unit_theta'] = np.asarray(expected)
            arrays[key + '_mask'] = np.ma.getmaskarray(expected)
            arrays[key + '_observed_mask'] = np.ma.getmaskarray(data)
            atomic_json(output / 'evaluations.json', record)
        for relative, expected_hash in config['frozen_sha256'].items():
            if digest(HERE / relative) != expected_hash:
                raise ValueError(f'Frozen source/input changed during evaluation: {relative}')
        if digest(manifest_path) != manifest_hash:
            raise ValueError('Manifest changed during evaluation')
        for row in record['completed_points']:
            if row['label'] == 'selected' and (
                abs(row['archived_ll_difference']) > config['selected_comparison_tolerance']['ll_absolute']
                or abs(row['archived_theta_relative_difference']) > config['selected_comparison_tolerance']['theta_relative']
            ):
                raise ValueError('Selected fixed-point result does not reproduce archived fit')
        np.savez_compressed(output / 'spectra.npz', **arrays)
        record['state'] = 'COMPLETED'
    except Exception as error:
        record.update(state='FAILED', failure=f'{type(error).__name__}: {error}')
    record['seconds'] = time.monotonic() - started
    atomic_json(output / 'evaluations.json', record)
    atomic_json(output / 'completion.json', {
        'manifest_sha256': manifest_hash, 'state': record['state'],
        'output_sha256': {p.name: digest(p) for p in output.iterdir() if p.is_file()}})
    print(json.dumps(record, indent=2), flush=True)
    if record['state'] != 'COMPLETED':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
