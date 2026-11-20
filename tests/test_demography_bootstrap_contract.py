"""Declared linked-bootstrap plans cannot silently drop or duplicate replicas."""
from __future__ import annotations

from copy import deepcopy
import json
import sys

import numpy as np
import pytest

from ghist_inference.bootstrap import (
    make_contract, replica_record, require_summary, summarize_payloads,
)
from scripts import bootstrap_growth_bgs_blocks as growth_producer
from scripts import bootstrap_split_migration_bgs_blocks as split_producer
from scripts import summarize_growth_bgs_bootstrap as growth_summary
from scripts import summarize_split_migration_bgs_bootstrap as split_summary


def synthetic_fit(task, replica=0, converged=True):
    settings = ({'projection': 84, 'mutation_rate': 3.36e-8, 'restarts': 2, 'convention': 'exponential'}
                if task == 3 else {'projection': [50, 40], 'mutation_rate': 1.015e-8, 'restarts': 2})
    fit = {**settings, 'generations': 1000. + replica * 100., 'post_growth_scaling': 2. + replica * .1,
           'AFR_population_size': 12000. + replica * 100., 'EUR_population_size': 9000. + replica * 100.,
           'generations_since_split': 1000. + replica * 100., 'migration_rate': 1e-4 + replica * 1e-5,
           'best_restart': 0, 'best_restart_converged': converged, 'multinomial_log_likelihood': -1.,
           'restart_diagnostics': [
               {'ok': True, 'restart': 0, 'likelihood': -1., 'optimizer': {'converged': converged, 'warnflag': 0 if converged else 2}},
               {'ok': True, 'restart': 1, 'likelihood': -2., 'optimizer': {'converged': True, 'warnflag': 0}},
           ]}
    ancestral = 10000.
    fit.update(sequence_length=400., theta=4 * settings['mutation_rate'] * 400. * ancestral)
    for row in fit['restart_diagnostics']:
        row['theta'] = fit['theta']
    if task == 3:
        fit.update(ancestral_effective_size=ancestral,
                   modern_effective_size=fit['post_growth_scaling'] * ancestral,
                   scaled_duration_2N=fit['generations'] / (2 * ancestral))
    else:
        fit.update(ancestral_population_size=ancestral, scaled_parameters={
            'AFR_size_ratio': fit['AFR_population_size'] / ancestral,
            'EUR_size_ratio': fit['EUR_population_size'] / ancestral,
            'split_duration_2N': fit['generations_since_split'] / (2 * ancestral),
            'migration_2N': fit['migration_rate'] * 2 * ancestral,
        })
    return fit


def bootstrap_fixture(tmp_path, task, planned=2):
    blocks = tmp_path / f'synthetic_blocks_{task}.npz'
    if not blocks.exists():
        np.savez_compressed(blocks, spectra=np.ones((4, 5) if task == 3 else (4, 5, 5)),
                            exposure_bp=np.full(4, 100.), contig=np.array(['a', 'a', 'b', 'b']),
                            start=np.array([0, 100, 0, 100]), end=np.array([100, 200, 100, 200]), block_size=100)
    initial = tmp_path / 'synthetic_initial.json'
    initial.write_text(json.dumps({'fit': synthetic_fit(4)}))
    settings = {key: synthetic_fit(task)[key] for key in ['projection', 'mutation_rate', 'restarts']}
    if task == 3:
        settings.update(convention='exponential', maxiter=400, seed=75194000,
                        fit_seed_stride=1009, resampling='within_contig_complete_blocks')
    else:
        settings.update(maxiter=600, resampling_seed=75198000, fit_seed=75199000,
                        resampling='within_contig_all_stored_blocks')
    contract = make_contract(task=task, planned_replicates=planned, blocks=blocks, settings=settings,
                             initial_fit=initial if task == 4 else None)
    rows = [replica_record(i, [0, 1, 2, 3], synthetic_fit(task, i), contract) for i in range(planned)]
    payload = {'bootstrap_contract': contract, 'state': 'COMPLETED',
               'assigned_replicate_ids': list(range(planned)), 'results': rows}
    return payload


@pytest.mark.parametrize('task', [3, 4])
def test_complete_declared_plan_recomputes_intervals(task, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    summary = summarize_payloads([payload], task)
    intervals = require_summary(summary, task)
    assert summary['replicates'] == summary['planned_replicates'] == 2
    assert summary['completed_replicates'] == 2
    assert summary['interval_coverage_validated'] is False
    key = 'generations' if task == 3 else 'generations_since_split'
    assert intervals[key] == [1002.5, 1097.5]


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('defect', ['prefix', 'duplicate', 'missing', 'unconverged', 'failed', 'running', 'wrong_id', 'status_spoof'])
def test_bad_plans_preserve_records_and_suppress_all_intervals(task, defect, tmp_path):
    payload = bootstrap_fixture(tmp_path, task, planned=4)
    if defect in ('prefix', 'missing'):
        payload['results'] = payload['results'][:2] if defect == 'prefix' else payload['results'][1:]
    elif defect == 'duplicate':
        payload['results'][-1] = deepcopy(payload['results'][0])
    elif defect == 'wrong_id':
        payload['results'][0]['replicate'] = 90
    elif defect == 'failed':
        payload['results'][0] = {'replicate': 0, 'state': 'FAILED', 'failure': 'preserved synthetic exception'}
    elif defect == 'running':
        payload['state'] = 'RUNNING'
    else:
        payload['results'][0]['fit']['best_restart_converged'] = defect == 'status_spoof'
        payload['results'][0]['fit']['restart_diagnostics'][0]['optimizer'] = {'converged': False, 'warnflag': 2}
    summary = summarize_payloads([payload], task)
    assert summary['state'] == 'FAILED'
    assert summary['replica_records'] == payload['results']
    assert 'parameters' not in summary and 'generations' not in summary
    with pytest.raises(ValueError, match='missing, duplicated, failed, or unconverged'):
        require_summary(summary, task)


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('change', ['contract', 'count', 'quantile', 'fit', 'legacy'])
def test_freeze_guard_rejects_unverifiable_or_changed_summary(task, change, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    summary = summarize_payloads([payload], task)
    if change == 'contract':
        summary['bootstrap_contract']['source_sha256'] = {}
    elif change == 'count':
        summary['planned_replicates'] = 1
    elif change == 'quantile':
        stats = summary if task == 3 else summary['parameters']
        next(stats[key] for key in ['generations', 'generations_since_split'] if key in stats)['quantiles']['0.025'] += 1.
    elif change == 'fit':
        summary['replica_records'] = deepcopy(summary['replica_records'])
        summary['replica_records'][0]['fit']['best_restart_converged'] = False
    else:
        del summary['schema']
    with pytest.raises(ValueError):
        require_summary(summary, task)


@pytest.mark.parametrize('task', [3, 4])
def test_mixed_protocol_shards_are_rejected(task, tmp_path):
    one = bootstrap_fixture(tmp_path, task)
    two = deepcopy(one)
    two['bootstrap_contract']['settings']['restarts'] += 1
    with pytest.raises(ValueError, match='different declared protocols'):
        summarize_payloads([one, two], task)


@pytest.mark.parametrize('task', [3, 4])
def test_actual_summarizer_writes_failed_denominator_without_quantiles(task, tmp_path, monkeypatch):
    payload = bootstrap_fixture(tmp_path, task, planned=4)
    payload['assigned_replicate_ids'] = [0, 1]
    payload['results'] = payload['results'][:2]
    (tmp_path / 'task_0.json').write_text(json.dumps(payload))
    out = tmp_path / 'summary.json'
    if task == 3:
        module = growth_summary
        argv = ['summary', '--input-dir', str(tmp_path), '--out', str(out), '--job', '1']
    else:
        module = split_summary
        blocks = tmp_path / 'block_summary.json'
        blocks.write_text('{}')
        argv = ['summary', '--task-dir', str(tmp_path), '--out', str(out), '--block-summary', str(blocks)]
    monkeypatch.setattr(sys, 'argv', argv)
    with pytest.raises(RuntimeError, match='interval output suppressed'):
        module.main()
    result = json.loads(out.read_text())
    assert result['planned_replicates'] == 4
    assert result['missing_replicate_ids'] == [2, 3]
    assert result['state'] == 'FAILED'
    assert 'parameters' not in result and 'generations' not in result


@pytest.mark.parametrize('task', [3, 4])
def test_producer_retains_unconverged_fit_and_never_overwrites(task, tmp_path, monkeypatch):
    payload = bootstrap_fixture(tmp_path, task)
    contract = payload['bootstrap_contract']
    out = tmp_path / 'task_0.json'
    argv = ['bootstrap', '--blocks', contract['inputs']['blocks']['path'], '--out', str(out),
            '--planned-replicates', '2', '--restarts', '2']
    if task == 3:
        module = growth_producer
        argv += ['--replicate-start', '0', '--replicates', '2']
        calls = []
        def fit(*args, **kwargs):
            calls.append(kwargs)
            return synthetic_fit(3, len(calls)-1, converged=len(calls) != 1)
        monkeypatch.setattr(module, 'fit_growth_spectrum', fit)
    else:
        module = split_producer
        argv += ['--initial-fit', contract['inputs']['initial_fit']['path'], '--replicate', '0']
        monkeypatch.setattr(module, 'fit_isolation_with_migration_spectrum', lambda *a, **k: synthetic_fit(4, converged=False))
    monkeypatch.setattr(sys, 'argv', argv)
    with pytest.raises(RuntimeError, match='attempted.*retained'):
        module.main()
    result = json.loads(out.read_text())
    assert result['state'] == 'FAILED'
    assert result['results'][0]['fit']['best_restart_converged'] is False
    assert result['bootstrap_contract']['planned_replicates'] == 2
    if task == 3:
        assert result['results'][1]['state'] == 'COMPLETED'
    previous = out.read_bytes()
    with pytest.raises(FileExistsError):
        module.main()
    assert out.read_bytes() == previous


@pytest.mark.parametrize('task', [3, 4])
def test_selected_converged_flag_cannot_hide_better_unconverged_restart(task, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    fit = payload['results'][0]['fit']
    fit['restart_diagnostics'][1].update(likelihood=0., optimizer={'converged': False, 'warnflag': 2})
    summary = summarize_payloads([payload], task)
    assert summary['state'] == 'FAILED'
    assert summary['failed_replicate_ids'] == [0]
    with pytest.raises(ValueError, match='unconverged'):
        require_summary(summary, task)


@pytest.mark.parametrize('task', [3, 4])
def test_contract_missing_planned_denominator_cannot_be_inferred_from_present_rows(task, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    del payload['bootstrap_contract']['planned_replicates']
    with pytest.raises(ValueError, match='planned_replicates'):
        summarize_payloads([payload], task)


@pytest.mark.parametrize('task', [3, 4])
def test_producer_preserves_exception_as_failed_artifact(task, tmp_path, monkeypatch):
    payload = bootstrap_fixture(tmp_path, task)
    contract = payload['bootstrap_contract']
    out = tmp_path / 'failed_task.json'
    argv = ['bootstrap', '--blocks', contract['inputs']['blocks']['path'], '--out', str(out),
            '--planned-replicates', '2', '--restarts', '2']
    def broken_fit(*args, **kwargs):
        raise RuntimeError('synthetic native failure')
    if task == 3:
        module = growth_producer
        argv += ['--replicate-start', '0', '--replicates', '2']
        monkeypatch.setattr(module, 'fit_growth_spectrum', broken_fit)
    else:
        module = split_producer
        argv += ['--initial-fit', contract['inputs']['initial_fit']['path'], '--replicate', '0']
        monkeypatch.setattr(module, 'fit_isolation_with_migration_spectrum', broken_fit)
    monkeypatch.setattr(sys, 'argv', argv)
    with pytest.raises(RuntimeError, match='attempted.*retained'):
        module.main()
    record = json.loads(out.read_text())
    assert record['state'] == 'FAILED'
    assert all('synthetic native failure' in row['failure'] for row in record['results'])


@pytest.mark.parametrize('task', [3, 4])
def test_producer_detects_input_change_without_losing_completed_fits(task, tmp_path, monkeypatch):
    from pathlib import Path
    payload = bootstrap_fixture(tmp_path, task)
    contract = payload['bootstrap_contract']
    block_path = Path(contract['inputs']['blocks']['path'])
    original_hash = contract['inputs']['blocks']['sha256']
    out = tmp_path / 'changed_task.json'
    argv = ['bootstrap', '--blocks', str(block_path), '--out', str(out),
            '--planned-replicates', '2', '--restarts', '2']
    calls = []
    def modifying_fit(*args, **kwargs):
        calls.append(True)
        if len(calls) == 1:
            block_path.write_bytes(block_path.read_bytes() + b'changed-after-extraction')
        return synthetic_fit(task, len(calls) - 1)
    if task == 3:
        module = growth_producer
        argv += ['--replicate-start', '0', '--replicates', '2']
        monkeypatch.setattr(module, 'fit_growth_spectrum', modifying_fit)
    else:
        module = split_producer
        argv += ['--initial-fit', contract['inputs']['initial_fit']['path'], '--replicate', '0']
        monkeypatch.setattr(module, 'fit_isolation_with_migration_spectrum', modifying_fit)
    monkeypatch.setattr(sys, 'argv', argv)
    with pytest.raises(RuntimeError, match='attempted.*retained'):
        module.main()
    record = json.loads(out.read_text())
    assert record['state'] == 'FAILED'
    assert 'input changed' in record['failure']
    assert record['bootstrap_contract']['inputs']['blocks']['sha256'] == original_hash
    assert all(row['state'] == 'COMPLETED' and 'fit' in row for row in record['results'])
    summary = summarize_payloads([record], task)
    assert summary['state'] == 'FAILED' and 'parameters' not in summary and 'generations' not in summary


@pytest.mark.parametrize('task', [3, 4])
def test_valid_bootstrap_cannot_attach_to_different_observed_fit_protocol(task, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    summary = summarize_payloads([payload], task)
    key = '84' if task == 3 else '50x40'
    fit = synthetic_fit(task)
    require_summary(summary, task, observed_fits={key: fit})
    fit['mutation_rate'] *= 2
    with pytest.raises(ValueError, match='mutation-rate identity differs'):
        require_summary(summary, task, observed_fits={key: fit})


@pytest.mark.parametrize('task', [3, 4])
def test_outer_summary_contract_cannot_override_retained_task_plan(task, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    summary = summarize_payloads([payload], task)
    summary['bootstrap_contract'] = deepcopy(summary['bootstrap_contract'])
    summary['bootstrap_contract']['planned_replicates'] = 3
    summary['bootstrap_contract']['planned_replicate_ids'] = [0, 1, 2]
    with pytest.raises(ValueError, match='differs from retained task contracts'):
        require_summary(summary, task)


@pytest.mark.parametrize('task,field', [
    (3, 'generations'), (3, 'post_growth_scaling'), (3, 'modern_effective_size'),
    (3, 'ancestral_effective_size'), (3, 'scaled_duration_2N'),
    (4, 'AFR_population_size'), (4, 'EUR_population_size'), (4, 'generations_since_split'),
    (4, 'migration_rate'), (4, 'ancestral_population_size'),
])
def test_inconsistent_physical_fit_is_retained_without_intervals(task, field, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    payload['results'][0]['fit'][field] *= 10
    summary = summarize_payloads([payload], task)
    assert summary['state'] == 'FAILED'
    assert summary['failed_replicate_ids'] == [0]
    assert summary['replica_records'] == payload['results']
    assert 'parameters' not in summary and 'generations' not in summary
    assert any('scaled-to-physical arithmetic' in error for error in summary['errors'])
    with pytest.raises(ValueError, match='missing, duplicated, failed, or unconverged'):
        require_summary(summary, task)


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('field', ['theta', 'selected_theta', 'sequence_length'])
@pytest.mark.parametrize('bad', [None, True, 0., -1., float('nan'), float('inf')])
def test_invalid_theta_or_exposure_cannot_authorize_intervals(task, field, bad, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    fit = payload['results'][0]['fit']
    if field == 'selected_theta':
        fit['restart_diagnostics'][0]['theta'] = bad
    else:
        fit[field] = bad
    summary = summarize_payloads([payload], task)
    assert summary['state'] == 'FAILED'
    assert summary['failed_replicate_ids'] == [0]
    assert 'parameters' not in summary and 'generations' not in summary


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('field', ['theta', 'selected_theta', 'sequence_length'])
def test_finite_theta_or_exposure_mismatch_cannot_change_intervals(task, field, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    fit = payload['results'][0]['fit']
    if field == 'selected_theta':
        fit['restart_diagnostics'][0]['theta'] *= 2
    else:
        fit[field] *= 2
    assert summarize_payloads([payload], task)['state'] == 'FAILED'


@pytest.mark.parametrize('task', [3, 4])
@pytest.mark.parametrize('length', [5e-324, 1e308])
def test_unrepresentable_theta_exposure_arithmetic_is_failed(task, length, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    fit = payload['results'][0]['fit']
    fit['sequence_length'] = length
    if length > 1:
        fit['mutation_rate'] = 1e308
    assert summarize_payloads([payload], task)['state'] == 'FAILED'


@pytest.mark.parametrize('task,field', [
    (3, 'scaled_duration_2N'), (3, 'post_growth_scaling'), (3, 'ancestral_effective_size'),
    (4, 'ancestral_population_size'), (4, 'AFR_size_ratio'), (4, 'EUR_size_ratio'),
    (4, 'split_duration_2N'), (4, 'migration_2N'),
])
@pytest.mark.parametrize('bad', [True, 0., float('nan')])
def test_invalid_scaled_parameters_or_ancestral_size_fail(task, field, bad, tmp_path):
    payload = bootstrap_fixture(tmp_path, task)
    fit = payload['results'][0]['fit']
    (fit['scaled_parameters'] if task == 4 and field != 'ancestral_population_size' else fit)[field] = bad
    assert summarize_payloads([payload], task)['state'] == 'FAILED'
