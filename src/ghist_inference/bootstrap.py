"""Declared linked-bootstrap denominators and retained-fit integrity, not coverage."""
from __future__ import annotations

import importlib.metadata
import json
import math
from pathlib import Path

import numpy as np

from .checkpoint import file_digest
from .freeze import (
    continuous_prediction, require_optimizer_convergence, validated_quantile_intervals,
    validate_physical_identity,
)

CONTRACT = 'demography_linked_bootstrap_v1'
PARAMETERS = {
    3: ('generations', 'post_growth_scaling'),
    4: ('AFR_population_size', 'EUR_population_size', 'generations_since_split', 'migration_rate'),
}
ROOT = Path(__file__).resolve().parents[2]


def source_hashes(task):
    stem = 'growth_bgs' if task == 3 else 'split_migration_bgs'
    paths = [
        'src/ghist_inference/bootstrap.py', 'src/ghist_inference/freeze.py',
        'src/ghist_inference/checkpoint.py', 'src/ghist_inference/validation.py',
        'src/ghist_inference/integration.py', f'src/ghist_{stem}/spectrum.py',
        f'scripts/bootstrap_{stem}_blocks.py', f'slurm/{stem}_block_bootstrap.sbatch',
    ]
    return {path: file_digest(ROOT / path) for path in paths}


def positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def make_contract(*, task, planned_replicates, blocks, settings, initial_fit=None):
    if task not in PARAMETERS:
        raise ValueError('Unsupported bootstrap task')
    planned_replicates = positive_integer(planned_replicates, 'planned_replicates')
    sources = {'blocks': {'path': str(blocks), 'sha256': file_digest(blocks)}}
    if initial_fit is not None:
        sources['initial_fit'] = {'path': str(initial_fit), 'sha256': file_digest(initial_fit)}
    return {
        'schema': CONTRACT, 'task': task, 'planned_replicates': planned_replicates,
        'planned_replicate_ids': list(range(planned_replicates)),
        'inputs': sources, 'settings': settings, 'source_sha256': source_hashes(task),
        'runtime_versions': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'moments-popgen']},
        'scope': 'Declared complete linked-block bootstrap; integrity checks do not validate interval coverage.',
    }


def validate_contract(contract, task):
    if not isinstance(contract, dict) or contract.get('schema') != CONTRACT or contract.get('task') != task:
        raise ValueError('Missing current bootstrap contract; regenerate legacy bootstrap evidence')
    total = positive_integer(contract.get('planned_replicates'), 'planned_replicates')
    planned_ids = contract.get('planned_replicate_ids')
    if (not isinstance(planned_ids, list)
            or any(isinstance(value, bool) or not isinstance(value, int) for value in planned_ids)
            or planned_ids != list(range(total))):
        raise ValueError('Bootstrap planned IDs must match its declared total')
    if contract.get('source_sha256') != source_hashes(task):
        raise ValueError('Bootstrap source identity differs; regenerate with the current producer')
    if not isinstance(contract.get('settings'), dict) or not isinstance(contract.get('runtime_versions'), dict):
        raise ValueError('Bootstrap settings/runtime identity missing')
    settings = contract['settings']
    for name in ['restarts', 'maxiter']:
        positive_integer(settings.get(name), name)
    rate = settings.get('mutation_rate')
    if isinstance(rate, bool) or not isinstance(rate, (int, float)) or not math.isfinite(rate) or rate <= 0:
        raise ValueError('Bootstrap mutation rate must be finite and positive')
    if task == 3:
        positive_integer(settings.get('projection'), 'projection')
        if settings['projection'] < 2 or settings.get('convention') != 'exponential':
            raise ValueError('Bootstrap growth projection/convention differs from its producer')
        seed = settings.get('seed')
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError('Bootstrap seed must be a nonnegative integer')
        if settings.get('fit_seed_stride') != 1009 or settings.get('resampling') != 'within_contig_complete_blocks':
            raise ValueError('Bootstrap growth resampling/seed protocol differs from its producer')
    elif (settings.get('projection') != [50, 40] or settings['mutation_rate'] != 1.015e-8
          or settings['maxiter'] != 600 or settings.get('resampling_seed') != 75198000
          or settings.get('fit_seed') != 75199000
          or settings.get('resampling') != 'within_contig_all_stored_blocks'):
        raise ValueError('Bootstrap split settings differ from its declared producer')
    versions = contract['runtime_versions']
    if set(versions) != {'numpy', 'scipy', 'moments-popgen'} or not all(isinstance(v, str) and v for v in versions.values()):
        raise ValueError('Bootstrap runtime version identity is incomplete')
    required_inputs = {'blocks'} if task == 3 else {'blocks', 'initial_fit'}
    if set(contract.get('inputs', {})) != required_inputs:
        raise ValueError('Bootstrap input identity missing')
    for identity in contract['inputs'].values():
        value = identity.get('sha256', '')
        if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
            raise ValueError('Bootstrap input SHA256 missing or malformed')
    return total


def verify_inputs(contract):
    for name, identity in contract['inputs'].items():
        if file_digest(identity['path']) != identity['sha256']:
            raise ValueError(f'Bootstrap {name} input changed during execution')
    validate_contract(contract, contract['task'])


def validate_fit(fit, contract):
    task = contract['task']
    require_optimizer_convergence(fit, 'bootstrap replica')
    values = continuous_prediction(fit, PARAMETERS[task])
    if task == 3 and values['post_growth_scaling'] <= 1:
        raise ValueError('Bootstrap growth ratio must be greater than one')
    diagnostics = fit.get('restart_diagnostics', [])
    requested = contract['settings']['restarts']
    ids = [row.get('restart') for row in diagnostics]
    if (len(ids) != requested or any(isinstance(value, bool) or not isinstance(value, int) for value in ids)
            or set(ids) != set(range(requested))):
        raise ValueError('Bootstrap retained restart denominator is incomplete or duplicated')
    best_id = fit.get('best_restart')
    if isinstance(best_id, bool) or not isinstance(best_id, int):
        raise ValueError('Bootstrap selected restart identity is invalid')
    selected = [row for row in diagnostics if row.get('restart') == best_id]
    if len(selected) != 1 or selected[0].get('ok') is not True:
        raise ValueError('Bootstrap selected restart lacks unique retained diagnostics')
    optimizer = selected[0].get('optimizer', {})
    if optimizer.get('converged') is not True or type(optimizer.get('warnflag')) is not int or optimizer['warnflag'] != 0:
        raise ValueError('Bootstrap selected restart convergence disagrees with retained diagnostics')
    if not math.isfinite(float(fit.get('multinomial_log_likelihood', float('nan')))):
        raise ValueError('Bootstrap selected likelihood is nonfinite')
    valid = [row for row in diagnostics if row.get('ok') is True]
    if any(not math.isfinite(float(row.get('likelihood', float('nan')))) for row in valid):
        raise ValueError('Bootstrap retained likelihood is nonfinite')
    best = max(valid, key=lambda row: row['likelihood'])
    if best['restart'] != best_id or selected[0]['likelihood'] != fit['multinomial_log_likelihood']:
        raise ValueError('Bootstrap selected fit does not match the best retained likelihood')
    validate_physical_identity(fit, task, selected=selected[0])
    for name in ['projection', 'mutation_rate', 'restarts'] + (['convention'] if task == 3 else []):
        if fit.get(name) != contract['settings'].get(name):
            raise ValueError(f'Bootstrap fitted {name} differs from the declared protocol')
    json.dumps(fit, allow_nan=False)
    return values


def replica_record(replicate, selected_blocks, fit, contract):
    """Retain finite but unconverged fits, without permitting their intervals."""
    record = {'replicate': replicate, 'selected_blocks': selected_blocks, 'fit': fit}
    try:
        validate_fit(fit, contract)
        record['state'] = 'COMPLETED'
    except (ValueError, TypeError, KeyError) as error:
        record.update(state='FAILED', failure=f'{type(error).__name__}: {error}')
        try:
            json.dumps(fit, allow_nan=False)
        except (ValueError, TypeError):
            record['invalid_fit_repr'] = repr(record.pop('fit'))
    return record


def summarize_records(records, contract):
    task = contract['task']
    total = validate_contract(contract, task)
    expected_ids = contract['planned_replicate_ids']
    by_id = {}
    errors = []
    for row in records:
        identity = row.get('replicate')
        if isinstance(identity, bool) or not isinstance(identity, int) or identity not in expected_ids:
            errors.append(f'Unexpected or malformed replicate ID: {identity!r}')
            continue
        by_id.setdefault(identity, []).append(row)
    missing = sorted(set(expected_ids) - set(by_id))
    duplicates = sorted(identity for identity, rows in by_id.items() if len(rows) != 1)
    failed = []
    for identity, rows in by_id.items():
        for row in rows:
            try:
                if row.get('state') != 'COMPLETED':
                    raise ValueError(f"replica state is {row.get('state')!r}")
                validate_fit(row.get('fit', {}), contract)
            except (ValueError, TypeError, KeyError) as error:
                failed.append(identity)
                errors.append(f'Replica {identity}: {error}')
    complete = not (missing or duplicates or failed or errors)
    result = {
        'schema': CONTRACT + '_summary', 'bootstrap_contract': contract,
        'state': 'COMPLETED' if complete else 'FAILED',
        'planned_replicates': total, 'replicates': len(records),
        'completed_replicates': len([identity for identity, rows in by_id.items()
                                     if len(rows) == 1 and identity not in failed]),
        'missing_replicate_ids': missing, 'duplicate_replicate_ids': duplicates,
        'failed_replicate_ids': sorted(set(failed)), 'errors': errors,
        'replica_records': records, 'interval_coverage_validated': False,
    }
    if complete:
        statistics = {}
        for parameter in PARAMETERS[task]:
            values = np.asarray([by_id[identity][0]['fit'][parameter] for identity in expected_ids], dtype=float)
            statistics[parameter] = {
                'mean': float(values.mean()),
                'standard_deviation': float(values.std(ddof=1)) if total > 1 else None,
                'quantiles': dict(zip(['0.025', '0.25', '0.5', '0.75', '0.975'],
                                      np.quantile(values, [.025, .25, .5, .75, .975]).tolist())),
                'minimum': float(values.min()), 'maximum': float(values.max()),
            }
        if task == 3:
            result.update(statistics)
        else:
            result['parameters'] = statistics
    return result



def summarize_payloads(payloads, task):
    if not payloads:
        raise ValueError('No bootstrap task records found')
    contract = payloads[0].get('bootstrap_contract')
    validate_contract(contract, task)
    records = []
    task_errors = []
    for index, payload in enumerate(payloads):
        if payload.get('bootstrap_contract') != contract:
            raise ValueError('Bootstrap tasks have different declared protocols or inputs')
        assigned = payload.get('assigned_replicate_ids')
        if (not isinstance(assigned, list) or not assigned
                or any(isinstance(value, bool) or not isinstance(value, int) for value in assigned)
                or len(set(assigned)) != len(assigned)
                or not set(assigned).issubset(contract['planned_replicate_ids'])):
            raise ValueError('Bootstrap task assignment is malformed or outside the plan')
        rows = payload.get('results')
        if not isinstance(rows, list):
            raise ValueError('Bootstrap task lacks retained results')
        records.extend(rows)
        if payload.get('state') != 'COMPLETED':
            task_errors.append(f"Task {index} state is {payload.get('state')!r}")
        row_ids = [row.get('replicate') for row in rows]
        if any(isinstance(value, bool) or not isinstance(value, int) for value in row_ids) or sorted(row_ids) != sorted(assigned):
            task_errors.append(f'Task {index} results differ from its declared assignment')
    result = summarize_records(records, contract)
    result['task_records'] = payloads
    result['task_errors'] = task_errors
    if task_errors:
        result['state'] = 'FAILED'
        for parameter in PARAMETERS[task]:
            result.pop(parameter, None)
        result.pop('parameters', None)
    return result


def require_summary(summary, task, *, observed_fits=None):
    """Recompute the complete planned denominator and intervals from retained fits."""
    if not isinstance(summary, dict) or summary.get('schema') != CONTRACT + '_summary':
        raise ValueError('Legacy bootstrap summary lacks verifiable planned IDs and fits; regenerate it')
    contract = summary.get('bootstrap_contract')
    validate_contract(contract, task)
    if observed_fits is not None:
        settings = contract['settings']
        primary = 84 if task == 3 else [50, 40]
        if settings['projection'] != primary:
            raise ValueError('Bootstrap projection differs from the declared freeze reference')
        for projection, fit in observed_fits.items():
            expected_projection = int(projection) if task == 3 else [int(value) for value in projection.split('x')]
            if fit.get('projection') != expected_projection or fit.get('mutation_rate') != settings['mutation_rate']:
                raise ValueError('Bootstrap/observed fit projection or mutation-rate identity differs')
            if task == 3 and fit.get('convention') != settings['convention']:
                raise ValueError('Bootstrap/observed growth convention differs')
    tasks = summary.get('task_records')
    if not isinstance(tasks, list):
        raise ValueError('Bootstrap summary lacks retained task records')
    recomputed = summarize_payloads(tasks, task)
    if contract != recomputed['bootstrap_contract']:
        raise ValueError('Bootstrap summary contract differs from retained task contracts')
    if summary.get('replica_records') != recomputed['replica_records']:
        raise ValueError('Bootstrap retained replicas disagree with retained task records')
    if recomputed['state'] != 'COMPLETED':
        raise ValueError('Bootstrap planned replicas are missing, duplicated, failed, or unconverged')
    for key in ['state', 'planned_replicates', 'replicates', 'completed_replicates',
                'missing_replicate_ids', 'duplicate_replicate_ids', 'failed_replicate_ids', 'errors',
                'interval_coverage_validated', 'task_errors']:
        if summary.get(key) != recomputed[key]:
            raise ValueError(f'Bootstrap summary {key} disagrees with retained replicas')
    recorded = summary if task == 3 else summary.get('parameters', {})
    expected = recomputed if task == 3 else recomputed['parameters']
    for parameter in PARAMETERS[task]:
        if recorded.get(parameter) != expected[parameter]:
            raise ValueError(f'Bootstrap {parameter} quantiles/statistics disagree with retained replicas')
    return validated_quantile_intervals(
        expected, PARAMETERS[task],
        minimum_exclusive={'post_growth_scaling': 1.0} if task == 3 else None,
        nonnegative=('migration_rate',) if task == 4 else (),
    )
