"""Audit fixed-point provenance and report descriptive likelihood/scaling contrasts."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    for relative, expected in manifest['frozen_sha256'].items():
        assert digest(HERE / relative) == expected, relative
    completion = json.loads((HERE / 'results/completion.json').read_text())
    assert completion['manifest_sha256'] == digest(HERE / 'manifest.json')
    for name, expected in completion['output_sha256'].items():
        assert digest(HERE / 'results' / name) == expected, name
    data = json.loads((HERE / 'results/evaluations.json').read_text())
    assert data['manifest_sha256'] == completion['manifest_sha256']
    assert data['state'] == completion['state'] == 'COMPLETED'
    assert data['versions'] == manifest['runtime_versions']
    rows = {(r['engine'], r['label']): r for r in data['completed_points']}
    assert len(data['completed_points']) == len(rows) == manifest['planned_points'] == 4
    assert set(rows) == {(e, label) for e in ['moments', 'dadi'] for label in ['true_shape', 'selected']}
    true_current = manifest['true_current_effective_size']
    with np.load(HERE / 'results/spectra.npz', allow_pickle=False) as spectra:
        for engine in ['moments', 'dadi']:
            true_mask = spectra[engine + '_true_shape_observed_mask']
            selected_mask = spectra[engine + '_selected_observed_mask']
            assert np.array_equal(true_mask, selected_mask)
            for label in ['true_shape', 'selected']:
                active = ~spectra[engine + '_' + label + '_observed_mask']
                expected = spectra[engine + '_' + label + '_expected_per_unit_theta']
                assert not np.any(spectra[engine + '_' + label + '_mask'][active])
                row = rows[engine, label]
                assert int(active.sum()) == row['active_bins']
                theta = row['active_observed_sum'] / expected[active].sum()
                assert math.isclose(theta, row['optimal_theta'], rel_tol=1e-12)
    for point in manifest['points']:
        row = rows[point['engine'], point['label']]
        assert row['parameters'] == point['parameters']
        assert math.isfinite(row['conditional_log_likelihood'])
        nu, scaled_time = row['parameters']
        n = row['optimal_theta'] / (4 * manifest['mutation_rate'] * manifest['sequence_length'])
        assert math.isclose(n, row['profiled_ancestral_effective_size'], rel_tol=1e-12)
        assert math.isclose(2 * n * scaled_time, row['profiled_generations'], rel_tol=1e-12)
        assert math.isclose(nu * n, row['profiled_current_effective_size'], rel_tol=1e-12)
    comparisons = []
    for engine in ['moments', 'dadi']:
        truth, selected = rows[engine, 'true_shape'], rows[engine, 'selected']
        assert truth['active_bins'] == selected['active_bins']
        assert truth['active_observed_sum'] == selected['active_observed_sum']
        assert abs(selected['archived_ll_difference']) <= 1e-8
        assert abs(selected['archived_theta_relative_difference']) <= 1e-12
        comparisons.append({
            'engine': engine,
            'selected_minus_true_shape_log_likelihood': selected['conditional_log_likelihood'] - truth['conditional_log_likelihood'],
            'true_shape_profiled_theta': truth['optimal_theta'],
            'true_shape_profiled_theta_relative_error': truth['optimal_theta'] / manifest['true_theta'] - 1,
            'selected_profiled_theta': selected['optimal_theta'],
            'selected_current_effective_size': selected['profiled_current_effective_size'],
            'selected_current_effective_size_relative_error': selected['profiled_current_effective_size'] / true_current - 1,
            'true_shape_profiled_current_effective_size': truth['profiled_current_effective_size'],
            'true_shape_profiled_current_effective_size_relative_error': truth['profiled_current_effective_size'] / true_current - 1,
            'selected_physical_relative_errors': {
                'generations': selected['profiled_generations'] / manifest['truth']['generations'] - 1,
                'ancestral_effective_size': selected['profiled_ancestral_effective_size'] / manifest['truth']['ancestral_effective_size'] - 1,
                'post_decline_fraction': selected['parameters'][0] / manifest['truth']['post_decline_fraction'] - 1,
            },
        })
    summary = {
        'schema': 'bottleneck_stress_fixed_point_audit_v1',
        'scope': manifest['scope'], 'manifest_sha256': digest(HERE / 'manifest.json'),
        'audit_source_sha256': digest(Path(__file__)),
        'output_attestation_sha256': digest(HERE / 'results/completion.json'),
        'planned_points': 4, 'completed_attested_points': 4,
        'true_theta': manifest['true_theta'], 'true_current_effective_size': true_current,
        'comparisons': comparisons, 'limits': manifest['interpretation_limits'],
        'scientifically_validated': False, 'release_ready': False,
    }
    (HERE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
