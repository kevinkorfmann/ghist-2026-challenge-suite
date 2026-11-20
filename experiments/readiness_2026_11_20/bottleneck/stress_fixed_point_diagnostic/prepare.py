"""Freeze the authorized existing-case truth and selected endpoints without evaluating them."""
from pathlib import Path
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
REPLAY = PARENT / 'production_guard_replay'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if (HERE / 'manifest.json').exists():
        raise ValueError('Refuse replacing frozen manifest')
    replay = json.loads((REPLAY / 'manifest.json').read_text())
    context_path = PARENT / 'production_stress_biological_context.json'
    context = json.loads(context_path.read_text())
    assert context['final_replay_manifest_sha256'] == digest(REPLAY / 'manifest.json')
    case = next(c for c in replay['cases'] if c['name'] == '100mb_case_16_new_stress')
    inputs = HERE / 'inputs'
    inputs.mkdir(exist_ok=False)
    shutil.copy2(context_path, inputs / 'biological_context.json')
    shutil.copy2(REPLAY / case['input'], inputs / 'spectrum.npz')
    assert digest(inputs / 'spectrum.npz') == case['input_sha256']
    shutil.copytree(REPLAY / 'source_snapshot', HERE / 'source_snapshot',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    truth = context['truth']
    true_parameters = [truth['post_decline_fraction'],
                       truth['generations'] / (2 * truth['ancestral_effective_size'])]
    points = []
    for engine in ['moments', 'dadi']:
        fit_path = REPLAY / 'results' / case['name'] / f'{engine}.json'
        assert digest(fit_path) == context['methods'][engine]['fit_sha256']
        destination = inputs / f'{engine}.json'
        shutil.copy2(fit_path, destination)
        fit = json.loads(fit_path.read_text())['fit']
        points.extend([
            {'engine': engine, 'label': 'true_shape', 'parameters': true_parameters},
            {'engine': engine, 'label': 'selected',
             'parameters': [fit['post_decline_fraction'], fit['scaled_time_2N']],
             'fit_record': str(destination.relative_to(HERE))},
        ])
    frozen = {str(path.relative_to(HERE)): digest(path)
              for path in HERE.rglob('*') if path.is_file() and '__pycache__' not in path.parts}
    manifest = {
        'schema': 'bottleneck_stress_fixed_point_protocol_v1',
        'scope': 'Posthoc four-point diagnostic of one existing development case selected after poor recovery was known; no independent confirmation, optimization, grid search, new biology or holdout access.',
        'case': case['name'], 'truth': truth,
        'true_theta': 4 * replay['mutation_rate'] * case['sequence_length'] * truth['ancestral_effective_size'],
        'true_current_effective_size': truth['ancestral_effective_size'] * truth['post_decline_fraction'],
        'points': points, 'planned_points': 4, 'planned_engines': 2,
        'spectrum': 'inputs/spectrum.npz', 'folded': replay['folded'],
        'mutation_rate': replay['mutation_rate'], 'sequence_length': case['sequence_length'],
        'sample_haplotypes': 114, 'dadi_grid_points': [134, 144, 154],
        'likelihood': 'Archived engine Inference.ll_multinom on exactly the fixed folded active-bin domain; theta profiled analytically at each of two fixed shape points per engine.',
        'physical_transform': 'N_ancestral=theta/(4*mu*L); generations=2*N_ancestral*T; N_current=nu*N_ancestral.',
        'selected_comparison_tolerance': {'ll_absolute': 1e-8, 'theta_relative': 1e-12},
        'runtime_versions': replay['runtime_versions'],
        'replay_manifest_sha256': digest(REPLAY / 'manifest.json'),
        'frozen_sha256': frozen,
        'allocation': {'cpus': 1, 'seconds': 300},
        'interpretation_limits': ['No optimization or global optimum claim', 'No identifiability conclusion from one realization', 'Posthoc development diagnostic, not independent validation', 'No changed production model/settings/thresholds', 'Profiled theta at true shape does not fix ancestral size to generating truth'],
    }
    (HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(digest(HERE / 'manifest.json'))


if __name__ == '__main__':
    main()
