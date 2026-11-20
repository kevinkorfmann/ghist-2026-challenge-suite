"""Freeze outputs must describe their retained fits and valid physical intervals."""
from __future__ import annotations

import json
import sys

import pytest

from ghist_growth_bgs.simulation import SIMULATION_SCALING_SCHEME
from scripts import freeze_growth_bgs as growth
from scripts import freeze_split_migration_bgs as split


def sources_for(task, tmp_path):
    from test_demography_bootstrap_contract import bootstrap_fixture, synthetic_fit
    from ghist_inference.bootstrap import summarize_payloads
    verified_bootstrap = summarize_payloads([bootstrap_fixture(tmp_path, task)], task)
    scaling = {'simulation_scaling_scheme': SIMULATION_SCALING_SCHEME}
    if task == 3:
        fit = synthetic_fit(3)
        sources = {
            'audit': {'files': {name: {'sha256': '0' * 64} for name in ['vcf', 'coding_bed', 'recombination_map']},
                      'segregating_biallelic_sites': 1, 'expected_haplotypes': 84,
                      'physical_bp_by_distance_stratum': {'noncoding_1kb_10kb': 1, 'noncoding_10kb_100kb': 1}},
            'middle-fit': {'strata': {'noncoding_1kb_100kb': {
                'by_projection': {str(p): dict(synthetic_fit(3), projection=p) for p in [84, 60, 40]},
                'median': {name: fit[name] for name in growth.PARAMETERS}}}},
            'bootstrap': {'replicates': 2, 'generations': {'quantiles': {'0.025': 500., '0.975': 1500.}},
                          'post_growth_scaling': {'quantiles': {'0.025': 1.5, '0.975': 2.5}}},
            'convention': {}, 'model-comparison': {'selected_convention': 'exponential', 'selection_consistent_across_groups': True},
            'smoke': {**scaling, 'preferred': 'exponential'},
            'grid': {**scaling, 'results': {'synthetic': {'relative_error': {'generations': 0., 'post_growth_scaling': 0.}}}},
            'calibration': {**scaling, 'selected_model': 'identity', 'validation': {}, 'evaluations': {}},
            'stress': {**scaling, 'freeze_gate_passed': True},
        }
        sources['bootstrap'] = verified_bootstrap
        return growth, sources, ['final-sealed-sha256']
    fit = synthetic_fit(4)
    sources = {
        'baseline': {'by_projection': {p: dict(synthetic_fit(4), projection=[int(value) for value in p.split('x')]) for p in ['50x40', '70x50', '94x68']},
                     'projection_median': {name: fit[name] for name in split.PARAMETERS},
                     'metadata': {'ancestral_reference_records': 1, 'haplotypes': [94, 68]}, 'physical_exposure_bp': 2},
        'refined': dict(fit),
        'bootstrap': {'replicates': 2, 'parameters': {
            name: {'quantiles': {'0.025': fit[name] * .5, '0.975': fit[name] * 1.5}}
            for name in split.PARAMETERS}},
        'calibration': {**scaling, 'selected_model': 'identity', 'validation': {}, 'truth_points': [], 'evaluations': {}},
        'near-recovery': {**scaling, 'relative_error': {name: 0. for name in split.PARAMETERS}},
    }
    sources['bootstrap'] = verified_bootstrap
    return split, sources, ['testing-vcf-sha256', 'coding-bed-sha256', 'recombination-map-sha256', 'final-sealed-sha256']


def fields(task, sources):
    if task == 3:
        block = sources['middle-fit']['strata']['noncoding_1kb_100kb']
        return block['by_projection'], block['median'], sources['bootstrap'], 'generations'
    block = sources['baseline']
    return block['by_projection'], block['projection_median'], sources['bootstrap']['parameters'], 'generations_since_split'


def freeze(module, sources, flags, tmp_path, monkeypatch):
    output = tmp_path / 'freeze.json'
    argv = ['freeze', '--out', str(output)]
    for flag in flags:
        argv.extend(['--' + flag, '0' * 64])
    for flag, payload in sources.items():
        path = tmp_path / (flag + '.json')
        path.write_text(json.dumps(payload))
        argv.extend(['--' + flag, str(path)])
    monkeypatch.setattr(sys, 'argv', argv)
    module.main()
    return json.loads(output.read_text())


@pytest.mark.parametrize('task', [3, 4])
def test_freeze_recomputes_median_and_preserves_valid_quantile_schema(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    fits, stored, _, parameter = fields(task, sources)
    projections = list(fits)
    fits[projections[1]][parameter] = 800.
    fits[projections[2]][parameter] = 1200.
    for fit in fits.values():
        if task == 3:
            fit['scaled_duration_2N'] = fit[parameter] / (2 * fit['ancestral_effective_size'])
        else:
            fit['scaled_parameters']['split_duration_2N'] = fit[parameter] / (2 * fit['ancestral_population_size'])
    # Only serialization-level rounding is tolerated; export the computed median.
    stored[parameter] = 1000. + 1e-10
    output = freeze(module, sources, flags, tmp_path, monkeypatch)
    assert output['testing_prediction'][parameter] == 1000.
    assert len(output['validation']['projection_fits']) == 3


@pytest.mark.parametrize('task', [3, 4])
def test_inconsistent_cached_median_cannot_pass_converged_projection_gate(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    _, stored, _, parameter = fields(task, sources)
    stored[parameter] = 2000.
    with pytest.raises(ValueError, match='stored median differs'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -1., 0., True, None, 'invalid', 2000.])
def test_invalid_interval_cannot_create_freeze_output(task, bad, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    _, _, intervals, parameter = fields(task, sources)
    intervals[parameter]['quantiles']['0.025'] = bad
    with pytest.raises(ValueError, match='interval|quantiles'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()
    assert not list(tmp_path.glob('.freeze.json.*'))


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('omission', ['parameter', 'quantile', 'malformed_quantiles'])
def test_required_interval_fields_fail_explicitly(task, omission, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    _, _, intervals, parameter = fields(task, sources)
    if omission == 'parameter':
        del intervals[parameter]
    elif omission == 'quantile':
        del intervals[parameter]['quantiles']['0.975']
    else:
        intervals[parameter]['quantiles'] = []
    with pytest.raises(ValueError, match='quantiles'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


def test_growth_interval_must_respect_existing_growth_ratio_domain(tmp_path, monkeypatch):
    module, sources, flags = sources_for(3, tmp_path)
    sources['bootstrap']['post_growth_scaling']['quantiles']['0.025'] = 1.
    with pytest.raises(ValueError, match='physical domain|quantiles'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


def test_zero_migration_interval_boundary_is_valid_without_widening_fit_policy():
    from ghist_inference.freeze import validated_intervals
    assert validated_intervals({'migration_rate': [0., .01]}, ('migration_rate',),
                               nonnegative=('migration_rate',)) == {'migration_rate': [0., .01]}


@pytest.mark.parametrize('task', [3, 4])
def test_nonfinite_embedded_evidence_does_not_replace_existing_artifact(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    sources['calibration']['evaluations']['invalid'] = float('nan')
    target = tmp_path / 'freeze.json'
    target.write_text('preserved old artifact')
    with pytest.raises(ValueError, match='JSON compliant'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert target.read_text() == 'preserved old artifact'


def test_every_growth_projection_must_respect_model_domain(tmp_path, monkeypatch):
    module, sources, flags = sources_for(3, tmp_path)
    fits, _, _, _ = fields(3, sources)
    fits['40']['post_growth_scaling'] = 0.9
    with pytest.raises(ValueError, match='scaling greater than one'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


@pytest.mark.parametrize('task', [3, 4])
def test_freeze_rejects_physical_time_changed_without_scaled_endpoint(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    parameter = 'generations' if task == 3 else 'generations_since_split'
    record = sources['bootstrap']['replica_records'][0]
    # The same retained record also appears under task_records. All selected
    # theta, scaled parameters, likelihood, and convergence evidence stay fixed.
    record['fit'][parameter] *= 10
    with pytest.raises(ValueError, match='missing, duplicated, failed, or unconverged'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()
    assert not list(tmp_path.glob('.freeze.json.*'))


@pytest.mark.parametrize('task', [3, 4])
def test_observed_time_and_matching_median_cannot_hide_inconsistent_scaled_fit(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    fits, stored, _, parameter = fields(task, sources)
    # Establish an existing valid artifact before exercising the actual gate.
    freeze(module, sources, flags, tmp_path, monkeypatch)
    before = (tmp_path / 'freeze.json').read_bytes()
    changed = list(fits) if task == 3 else ['70x50', '94x68']
    for projection in changed:
        fits[projection][parameter] *= 10
    stored[parameter] *= 10
    # Task4's 50x40 reference and refined fit are deliberately untouched.
    with pytest.raises(ValueError, match='scaled-to-physical arithmetic'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert (tmp_path / 'freeze.json').read_bytes() == before
    assert not list(tmp_path.glob('freeze.json.*'))


@pytest.mark.parametrize('nested', [False, True])
@pytest.mark.parametrize('corrupt', [False, True])
def test_refined_split_arithmetic_checked_in_both_historical_layouts(nested, corrupt, tmp_path, monkeypatch):
    module, sources, flags = sources_for(4, tmp_path)
    refined = sources['refined']
    if corrupt:
        # This is below the existing 1% refinement threshold but violates the
        # retained scaled endpoint; numerical agreement alone must not pass it.
        refined['generations_since_split'] *= 1.0001
    if nested:
        sources['refined'] = {'fit': refined}
    if corrupt:
        with pytest.raises(ValueError, match='scaled-to-physical arithmetic'):
            freeze(module, sources, flags, tmp_path, monkeypatch)
        assert not (tmp_path / 'freeze.json').exists()
    else:
        result = freeze(module, sources, flags, tmp_path, monkeypatch)
        assert result['testing_prediction']['generations_since_split'] == 1000.


@pytest.mark.parametrize('task', [3, 4])
def test_arithmetic_complete_observed_fit_without_restart_rows_remains_supported(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    fits, stored, _, _ = fields(task, sources)
    observed = list(fits.values()) + ([sources['refined']] if task == 4 else [])
    for fit in observed:
        del fit['restart_diagnostics']
        del fit['best_restart']
    # A componentwise median is an aggregate, not a complete model fit. It has
    # no theta/ancestral/scaled fields and must not receive fit arithmetic checks.
    assert set(stored) == set(module.PARAMETERS)
    result = freeze(module, sources, flags, tmp_path, monkeypatch)
    assert result['testing_prediction'] == stored


@pytest.mark.parametrize('task', [3, 4])
def test_observed_available_selected_theta_cannot_disagree(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    fits, _, _, _ = fields(task, sources)
    next(iter(fits.values()))['restart_diagnostics'][0]['theta'] *= 2
    with pytest.raises(ValueError, match='theta differs'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()


@pytest.mark.parametrize('task', [3, 4])
def test_partial_observed_arithmetic_metadata_is_rejected_not_skipped(task, tmp_path, monkeypatch):
    module, sources, flags = sources_for(task, tmp_path)
    fits, _, _, _ = fields(task, sources)
    del next(iter(fits.values()))['theta']
    with pytest.raises(ValueError, match='arithmetic metadata must be complete'):
        freeze(module, sources, flags, tmp_path, monkeypatch)
    assert not (tmp_path / 'freeze.json').exists()
