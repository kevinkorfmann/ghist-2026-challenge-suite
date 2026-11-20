"""Report diagnostic recovery without calibrating or selecting an estimator."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import moments
import numpy as np


def summarize(path: Path) -> dict:
    record = json.loads(path.read_text())
    if not record['complete']:
        raise ValueError('Recovery panel is incomplete')
    groups = {}
    cases = []
    names = ['generations', 'post_decline_fraction', 'ancestral_effective_size']
    for case in record['cases']:
        data = moments.Spectrum(case['spectrum']).fold()
        null = moments.Demographics1D.snm(data.sample_sizes).fold()
        null_ll = float(moments.Inference.ll_multinom(null, data))
        truth = case['truth']
        expected = moments.Demographics1D.two_epoch(
            [truth['post_decline_fraction'], truth['generations']/(2*truth['ancestral_effective_size'])],
            data.sample_sizes,
        ).fold()
        truth_ll = float(moments.Inference.ll_multinom(expected, data))
        comparisons = {}
        for method, result in case['methods'].items():
            group = groups.setdefault(f"{case['scenario']}/{method}", dict(errors=[], failures=0))
            if 'failure' in result:
                group['failures'] += 1
                continue
            group['errors'].append([abs(result['relative_error'][key]) for key in names])
            comparisons[method] = dict(
                log_likelihood_gain_over_constant=result['fit']['multinomial_log_likelihood']-null_ll,
                fitted_minus_two_epoch_at_target_parameters_log_likelihood=result['fit']['multinomial_log_likelihood']-truth_ll,
            )
        cases.append({key:case[key] for key in ['point','replicate','scenario','truth']} | {'likelihood_comparisons':comparisons})
    summary = {}
    for name, group in groups.items():
        errors = np.asarray(group['errors'])
        summary[name] = dict(successes=len(errors), failures=group['failures'])
        if len(errors):
            summary[name].update(
                median_absolute_relative_error=dict(zip(names, np.median(errors,axis=0))),
                maximum_absolute_relative_error=dict(zip(names, np.max(errors,axis=0))),
            )
    return dict(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Descriptive diagnostic panel; no calibration, model selection, confidence intervals, or hypothesis-test p-values',
                summary=summary, cases=cases)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    args.out.write_text(json.dumps(summarize(args.input),indent=2)+'\n')
