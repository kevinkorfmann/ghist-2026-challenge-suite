"""Summarize every planned case, retaining failures and incomplete cases."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_inference.checkpoint import atomic_json, load_spectrum_checkpoint
from scripts.validate_bottleneck_science import panel_work


def summarize(directory, *, allow_legacy_contracts=False):
    config = json.loads((directory / 'protocol.json').read_text())
    protocol_hash = hashlib.sha256((directory / 'protocol.json').read_bytes()).hexdigest()
    cases = []
    groups = {}
    for index, (point, replica, scenario, length, seed, restarts) in enumerate(panel_work(
        2, config['sequence_length'], config['seed'], config['restarts'],
    )):
        case = dict(case=index, point=point, replicate=replica, scenario=scenario, methods={})
        for method in (['moments', 'dadi'] if scenario == 'matched' else ['moments']):
            file = directory / f'case_{index:02d}' / f'{method}.json'
            result = json.loads(file.read_text()) if file.exists() else {'incomplete': True}
            expected = dict(case=index, point=point, replicate=replica, scenario=scenario,
                            ancestry_seed=seed, mutation_seed=seed+1, optimizer_seed=seed+2,
                            protocol_sha256=protocol_hash)
            if file.exists() and result.get('identity') != expected:
                raise ValueError(f'Fit identity or protocol mismatch: {file}')
            if file.exists():
                _, metadata = load_spectrum_checkpoint(file.with_name('spectrum.npz'), expected=expected)
                if config.get('schema') != 'bottleneck_target_length_v1' and result.get('spectrum_sha256') != metadata['spectrum_sha256']:
                    raise ValueError(f'Fit spectrum mismatch: {file}')
                truth = dict(zip(('generations', 'post_decline_fraction', 'ancestral_effective_size'),
                                 config['parameter_points'][point], strict=True))
                if result.get('truth') != truth or metadata.get('truth') != truth:
                    raise ValueError(f'Fit or spectrum truth differs from locked protocol: {file}')
                if 'fit' in result:
                    errors = {key: result['fit'][key] / value - 1 for key, value in truth.items()}
                    if not np.isfinite(list(errors.values())).all():
                        raise ValueError(f'Nonfinite fit: {file}')
                    if any(not np.isclose(result.get('relative_error', {}).get(key, np.nan), value,
                                          rtol=1e-12, atol=1e-12) for key, value in errors.items()):
                        raise ValueError(f'Stored error differs from recomputed error: {file}')
                    result = {**result, 'relative_error': errors}
            case['methods'][method] = result
            groups.setdefault(f'{scenario}/{method}', []).append(result)
        if index in (0, 12):
            file = directory / f'case_{index:02d}' / 'vcf_contract.json'
            case['vcf_contract'] = json.loads(file.read_text()) if file.exists() else {'passed': False}
            if file.exists():
                _, metadata = load_spectrum_checkpoint(file.with_name('spectrum.npz'), expected=expected)
                check = case['vcf_contract']
                if check.get('identity') != expected or check.get('spectrum_sha256') != metadata['spectrum_sha256']:
                    if not allow_legacy_contracts or config.get('schema') != 'bottleneck_target_length_v1':
                        raise ValueError(f'Unverified VCF contract identity: {file}')
                    case['vcf_contract'] = dict(passed=False, verification_status='unverified_legacy_identity',
                                               original_record=check)
                elif check.get('contig_length') != length:
                    raise ValueError(f'Wrong VCF contract contig length: {file}')
        cases.append(case)
    summary = {}
    for group, rows in groups.items():
        finite = [row for row in rows if 'fit' in row]
        summary[group] = dict(
            planned=len(rows), fitted=len(finite), failed=sum('failure' in row for row in rows),
            incomplete=sum('incomplete' in row for row in rows),
            optimizer_converged=sum(row['fit']['optimizer']['converged'] for row in finite),
            error_scope='Completed cases only; cannot represent full panel until complete',
            absolute_relative_error={
                key: dict(zip(('median', 'p90', 'maximum'), np.quantile([
                    abs(row['relative_error'][key]) for row in finite
                ], [.5, .9, 1]).tolist()))
                for key in ('generations', 'post_decline_fraction', 'ancestral_effective_size')
            } if finite else {},
        )
    return dict(
        protocol=config, cases=cases, summary=summary,
        complete=all(not row['incomplete'] and not row['failed'] for row in summary.values()),
        vcf_checks_passed=all(case.get('vcf_contract', {'passed': True})['passed'] for case in cases),
        scientific_status='diagnostic evidence; not a full-prior accuracy or uncertainty-coverage certification',
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--allow-legacy-contracts', action='store_true',
                        help='Retain original v1 VCF checks as unverified; never counts them as passing')
    args = parser.parse_args()
    result = summarize(args.directory, allow_legacy_contracts=args.allow_legacy_contracts)
    atomic_json(args.out, result)
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    main()
