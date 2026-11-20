"""Independent reproductions for projection-summary convergence gates."""
import json
import sys

import pytest

from ghist_inference.freeze import require_projection_convergence
from scripts import summarize_growth_bgs_fits


def test_growth_summary_rejects_duplicate_projection_hiding_failed_fit(tmp_path, monkeypatch):
    arguments = ['summary']
    for index, (time, converged) in enumerate(((1000, False), (3000, True))):
        path = tmp_path / f'fit_{index}.json'
        path.write_text(json.dumps({
            'projection': 40,
            'results': {'noncoding_1kb_100kb': {
                'generations': time, 'post_growth_scaling': 2,
                'ancestral_effective_size': 10000, 'modern_effective_size': 20000,
                'best_restart_converged': converged,
            }},
        }))
        arguments.extend(['--fit', str(path), '--job', str(index + 1)])
    output = tmp_path / 'summary.json'
    monkeypatch.setattr(sys, 'argv', arguments + ['--out', str(output)])
    with pytest.raises(ValueError, match='Duplicate projection'):
        summarize_growth_bgs_fits.main()
    assert not output.exists()


@pytest.mark.parametrize('required', [('84', '60', '40'), ('50x40', '70x50', '94x68')])
def test_freeze_requires_every_planned_projection_and_no_extras(required):
    fit = {'best_restart_converged': True}
    complete = {projection: fit for projection in required}
    require_projection_convergence(complete, expected_projections=required)
    for changed in ({required[0]: fit}, {**complete, 'unexpected': fit}):
        with pytest.raises(ValueError, match='required projection set'):
            require_projection_convergence(changed, expected_projections=required)


@pytest.mark.parametrize('task', [3, 4])
def test_freeze_entrypoints_reject_incomplete_converged_projection_set(tmp_path, monkeypatch, task):
    from scripts import freeze_growth_bgs, freeze_split_migration_bgs

    candidate = {'best_restart_converged': True}
    if task == 3:
        module = freeze_growth_bgs
        flags = ['audit', 'middle-fit', 'bootstrap', 'convention', 'model-comparison',
                 'smoke', 'grid', 'calibration', 'stress', 'final-sealed-sha256']
        sources = {flag: {} for flag in flags}
        sources['middle-fit'] = {'strata': {'noncoding_1kb_100kb': {
            'by_projection': {'40': candidate},
        }}}
    else:
        module = freeze_split_migration_bgs
        flags = ['baseline', 'refined', 'bootstrap', 'calibration', 'near-recovery',
                 'testing-vcf-sha256', 'coding-bed-sha256', 'recombination-map-sha256',
                 'final-sealed-sha256']
        sources = {flag: {} for flag in flags}
        sources['baseline'] = {'by_projection': {'50x40': candidate}}
    monkeypatch.setattr(module, 'load', lambda path: sources[path.name])
    monkeypatch.setattr(module, 'require_physical_scaling', lambda records: None)
    output = tmp_path / 'freeze.json'
    arguments = ['freeze', '--out', str(output)]
    for flag in flags:
        arguments += ['--' + flag, flag]
    monkeypatch.setattr(sys, 'argv', arguments)
    with pytest.raises(ValueError, match='required projection set'):
        module.main()
    assert not output.exists()


def secondary_contact_sources():
    from ghist_secondary_contact.contract import SPECTRUM_CONTRACT
    from scripts.freeze_secondary_contact import PARAMETER_KEYS

    raw = dict.fromkeys(PARAMETER_KEYS, 1.0)
    raw.update(generations_since_split=2., best_restart_converged=True,
               projection=[4, 4], integration_dt_fac=.02, adaptive_integration=False)
    contract = {'spectrum_contract': SPECTRUM_CONTRACT}
    resolution = {**contract, 'selection_valid': True, 'selected_resolution': 'p40',
                  'replicates': [], 'summaries': {}, 'paired_comparisons': {}, 'ranking': []}
    return {
        'raw-fit': dict(raw), 'standard-fit': dict(raw),
        'calibration': {**contract, 'observed_accepted': {k: raw[k] for k in PARAMETER_KEYS},
                        'untouched_evaluation': {'correction_accepted': True},
                        'accepted_prediction_method': 'synthetic', 'extra_calibration_patterns': [],
                        'point_grouped_partitions': {}, 'selected_method': 'synthetic'},
        'resolution': resolution, 'broad-resolution': resolution,
        'varying-resolution': resolution,
        'stress-comparison': {**contract, 'selection_valid': True, 'summaries': {}},
        'stress-transfer': {**contract, 'conditions': {'nominal': {
            'corrected': {'mean_competition_rrmse': .1},
            'raw': {'mean_competition_rrmse': .2},
        }}},
        'bootstrap': {**contract, 'observed_sampling_interval_95': {}},
        'vcf-audit': {'population': {'file_sha256': '0' * 64, 'records': 1,
                                    'skipped_non_biallelic_records': 0,
                                    'ancestral_reference_records': 1}},
    }


def run_secondary_contact_freeze(sources, tmp_path, monkeypatch):
    from scripts import freeze_secondary_contact

    output = tmp_path / 'freeze.json'
    arguments = ['freeze', '--out', str(output)]
    for flag, source in sources.items():
        path = tmp_path / (flag + '.json')
        path.write_text(json.dumps(source))
        arguments += ['--' + flag, str(path)]
    monkeypatch.setattr(sys, 'argv', arguments)
    freeze_secondary_contact.main()
    return json.loads(output.read_text())


@pytest.mark.parametrize('location', ['prediction', 'raw', 'standard', 'corrected_risk', 'raw_risk', 'interval'])
def test_secondary_contact_real_freeze_rejects_nonfinite_values(tmp_path, monkeypatch, location):
    sources = secondary_contact_sources()
    if location == 'prediction':
        sources['calibration']['observed_accepted']['population1_effective_size'] = float('nan')
    elif location in ('raw', 'standard'):
        sources[location + '-fit']['population1_effective_size'] = float('nan')
    elif location.endswith('_risk'):
        role = location.removesuffix('_risk')
        sources['stress-transfer']['conditions']['nominal'][role]['mean_competition_rrmse'] = float('nan')
    else:
        sources['bootstrap']['observed_sampling_interval_95'] = {'bad': float('nan')}
    with pytest.raises(ValueError, match='finite|JSON compliant'):
        run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


@pytest.mark.parametrize('change', [{'population1_effective_size': 0.},
                                  {'migration_rate': -1.},
                                  {'generations_since_contact': 3.}])
def test_secondary_contact_real_freeze_rejects_invalid_physical_domain(tmp_path, monkeypatch, change):
    sources = secondary_contact_sources()
    sources['calibration']['observed_accepted'].update(change)
    with pytest.raises(ValueError, match='finite|precede'):
        run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


@pytest.mark.parametrize('contact_time', [0., 2.])
def test_secondary_contact_real_freeze_preserves_physical_boundaries(tmp_path, monkeypatch, contact_time):
    sources = secondary_contact_sources()
    for source in (sources['raw-fit'], sources['standard-fit'],
                   sources['calibration']['observed_accepted']):
        source.update(generations_since_contact=contact_time, migration_rate=0.)
    result = run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert result['testing_prediction']['migration_rate'] == 0.
    assert result['testing_prediction']['generations_since_contact'] == contact_time
    assert result['validation']['numerical_convergence']['refined_minus_standard_relative']['migration_rate'] == 0.


@pytest.mark.parametrize('status', ['unknown', None])
def test_unknown_fit_checkpoint_status_cannot_call_fitter_or_overwrite(tmp_path, status):
    from ghist_inference.checkpoint import fit_checkpoint

    path = tmp_path / 'fit.json'
    spectrum = [0, 3, 2, 0]
    calls = []

    def fit(values, **settings):
        calls.append(True)
        return {'value': 1}

    fit_checkpoint(path, spectrum=spectrum, fit=fit, settings={}, resume=False)
    record = json.loads(path.read_text())
    record['status'] = status
    path.write_text(json.dumps(record))
    before = path.read_bytes()
    calls.clear()
    with pytest.raises(ValueError, match='unknown status'):
        fit_checkpoint(path, spectrum=spectrum, fit=fit, settings={}, resume=True)
    assert calls == []
    assert path.read_bytes() == before
