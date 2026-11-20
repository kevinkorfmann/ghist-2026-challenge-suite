"""Assess every predeclared cohort, preserving failures and fixed numeric gates."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from calibrate_phase import metrics
from ghist_relatedness.labels import HUMAN_LABELS
from structural_confidence import reliability

from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs


def numerical_gates(score, baseline):
    return {'all_labels': score['complete_label_coverage'],
        'exact_accuracy': score['exact_accuracy'] >= .8,
        'macro_recall': score['macro_recall'] >= .7,
        'minimum_label_recall': score['minimum_label_recall'] >= .5,
        'gain_over_baseline': score['exact_accuracy'] - baseline['exact_accuracy'] >= .1,
        'no_cohort_regression': all(value >= baseline['cohort_accuracy'][cohort]
            for cohort, value in score['cohort_accuracy'].items())}


def assess(root, protocol_path, accounting_path, stage, output):
    protocol = json.loads(protocol_path.read_text())
    accounting = json.loads(accounting_path.read_text())
    indices = list(range(8 if stage == 'small' else 10))
    labels = list(HUMAN_LABELS.values())
    variants = {name: [] for name in ['structural', 'pair_only', 'baseline']}
    inputs = {str(path): sha256(path) for path in [protocol_path, accounting_path]}
    failures, runtime, roots = [], {}, []
    for index in indices:
        directory = root / f'cohort_{index:03d}'
        required = [directory / name for name in ['result.json', 'identifiers.json', 'truth.json',
            'structural_predictions.json', 'pair_only_predictions.json', 'baseline_predictions.json',
            'structural_tree.json', 'missingness_contract/result.json']]
        if any(not path.is_file() for path in required):
            failures.append({'cohort_index': index, 'reason': 'Missing required artifacts',
                'missing': [str(path) for path in required if not path.is_file()]})
            continue
        inputs.update({str(path): sha256(path) for path in required})
        record = json.loads(required[0].read_text())
        verify_inputs(record['output_sha256'])
        roster = sorted(json.loads(required[1].read_text()).values())
        truth = json.loads(required[2].read_text())
        expected_pairs = set(itertools.combinations(roster, 2))
        actual_pairs = {tuple(sorted([row['indID1'], row['indID2']])) for row in truth}
        if len(roster) != protocol['cohorts'][index]['people'] or len(set(roster)) != len(roster):
            raise ValueError('Declared roster count changed')
        if actual_pairs != expected_pairs or len(truth) != len(expected_pairs):
            raise ValueError('Truth omits/duplicates full roster pairs')
        for name in variants:
            path = directory / f'{name}_predictions.json'
            variants[name].append({'cohort_id': f'nested_{index:03d}', 'truth': truth,
                'predictions': json.loads(path.read_text())})
        audit = accounting.get(str(index))
        if audit is None or audit['job_id'] != record['slurm_job_id']:
            failures.append({'cohort_index': index, 'reason': 'Missing/mismatched final Slurm accounting'})
        else:
            runtime[str(index)] = {'completed': audit['state'] == 'COMPLETED',
                'wall_time': audit['elapsed_seconds'] <= (5400 if index >= 8 else 2700),
                'memory': audit['max_rss_kib'] > 0 and audit['max_rss_kib'] <= (44 if index >= 8 else 22) * 1024**2}
        contract = json.loads(required[-1].read_text())
        if contract['state'] != 'MISSINGNESS_REJECTION_VERIFIED' or contract['biological_robustness_qualified']:
            failures.append({'cohort_index': index, 'reason': 'Missingness rejection contract failed'})
        roots.append({'cohort_index': index, **json.loads((directory / 'structural_tree.json').read_text())})
    output.mkdir(parents=True, exist_ok=False)
    for name, cohorts in variants.items():
        save_json(output / f'{name}_scorer_input.json', {'labels': labels, 'cohorts': cohorts})
    if failures or len(variants['structural']) != len(indices):
        save_json(output / 'result.json', {'state': 'INCOMPLETE_QUALIFICATION_FAILED', 'stage': stage,
            'failures': failures, 'passed': False, 'runtime': runtime, 'input_sha256': inputs,
            'missing_cohorts_excluded_from_gate_scoring': False, 'submission_ready': False})
        return
    scores = {name: metrics(cohorts, labels) for name, cohorts in variants.items()}
    small = {name: metrics(cohorts[:8], labels) for name, cohorts in variants.items()}
    gates = {'small_nested': numerical_gates(small['structural'], small['baseline'])}
    reliability_report = reliability(variants['structural'], labels)
    small_reliability = reliability(variants['structural'][:8], labels)
    gates['small_nested']['ece'] = small_reliability['all']['ece_10bins'] <= .1
    if stage == 'complete':
        for index in [8, 9]:
            by_variant = {name: metrics([cohorts[index]], labels) for name, cohorts in variants.items()}
            gates[f'large_{index}'] = numerical_gates(by_variant['structural'], by_variant['baseline'])
            gates[f'large_{index}']['ece'] = reliability_report[f'cohort:nested_{index:03d}']['ece_10bins'] <= .1
    gates['all'] = {'ece': reliability_report['all']['ece_10bins'] <= .1,
        'complete_predeclared_cohorts': len(variants['structural']) == len(indices),
        'runtime_and_memory': len(runtime) == len(indices) and all(all(row.values()) for row in runtime.values())}
    save_json(output / 'reliability.json', {'all': reliability_report, 'small_nested': small_reliability})
    save_json(output / 'root_diagnostics.json', roots)
    verify_inputs(inputs)
    save_json(output / 'result.json', {'state': 'ASSESSED_FIXED_QUALIFICATION_GATES', 'stage': stage,
        'scores': scores, 'small_scores': small, 'gates': gates,
        'passed': all(all(group.values()) for group in gates.values()),
        'runtime': runtime, 'input_sha256': inputs, 'protocol_sha256': sha256(protocol_path),
        'submission_ready': False, 'requires_independent_rescore': True,
        'confidence_interpretation': 'Global correctness estimate only; zero individual resolution; no root posterior',
        'domain': protocol['domain']})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--accounting', type=Path, required=True)
    parser.add_argument('--stage', choices=['small', 'complete'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assess(args.root, args.protocol, args.accounting, args.stage, args.output)


if __name__ == '__main__':
    main()
