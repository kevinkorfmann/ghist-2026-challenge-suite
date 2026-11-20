"""Apply the prospectively fixed pragmatic thinning gates, with full denominators."""
from __future__ import annotations
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path
from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs


def measure(cohorts, labels):
    totals, correct = Counter(), Counter()
    for cohort in cohorts:
        truth = {(min(r['indID1'],r['indID2']), max(r['indID1'],r['indID2'])): r['relation'] for r in cohort['truth']}
        predictions = {(min(r['indID1'],r['indID2']), max(r['indID1'],r['indID2'])): r['relation'] for r in cohort['predictions']}
        if len(truth) != len(cohort['truth']) or len(predictions) != len(cohort['predictions']) or set(truth) != set(predictions):
            raise ValueError('Missing, duplicate or extra pairs')
        for key, label in truth.items():
            totals[label] += 1
            correct[label] += predictions[key] == label
    recalls = {label: correct[label] / totals[label] for label in labels}
    return {'pairs': sum(totals.values()), 'correct': sum(correct.values()), 'accuracy': sum(correct.values())/sum(totals.values()), 'macro_recall': sum(recalls.values())/len(labels), 'per_label_recall': recalls}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    protocol_path = args.source / 'pragmatic_missingness_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    intact_path = args.root / 'integer_edge_decoder_v1/pair_only_scorer_input.json'
    intact = json.loads(intact_path.read_text())
    cohorts, bound = [], {str(intact_path): sha256(intact_path), str(protocol_path): sha256(protocol_path)}
    for index, original in zip(protocol['development_cohorts'], intact['cohorts'], strict=True):
        folder = args.root / 'site_thinning_diagnostic_v1' / f'cohort_{index:03d}'
        result = json.loads((folder / 'result.json').read_text())
        verify_inputs(result['input_sha256'])
        verify_inputs(result['output_sha256'])
        if result['index'] != index or result['protocol_sha256'] != sha256(protocol_path):
            raise ValueError('Wrong cohort or policy')
        for filename in ['result.json','samples.json','predictions.json','thinning.json']:
            bound[str(folder / filename)] = sha256(folder / filename)
        predictions = json.loads((folder / 'predictions.json').read_text())
        samples = json.loads((folder / 'samples.json').read_text())
        expected = set(itertools.combinations(sorted(samples), 2))
        if expected != {tuple(sorted([r['indID1'],r['indID2']])) for r in original['truth']}:
            raise ValueError('Intact comparator differs from exact roster')
        cohorts.append({'cohort_id': original['cohort_id'], 'truth': original['truth'], 'predictions': predictions})
    before, after = [measure(rows, intact['labels']) for rows in [intact['cohorts'], cohorts]]
    gates = {'complete_all_pairs': before['pairs'] == after['pairs'] == 585,
        'exact_loss_at_most_0_03': before['accuracy']-after['accuracy'] <= .03,
        'macro_loss_at_most_0_03': before['macro_recall']-after['macro_recall'] <= .03,
        'no_previously_adequate_class_below_0_5': all(after['per_label_recall'][label] >= .5 for label, value in before['per_label_recall'].items() if value >= .5)}
    args.output.mkdir(parents=True, exist_ok=False)
    save_json(args.output / 'scorer_input.json', {'labels': intact['labels'], 'cohorts': cohorts})
    save_json(args.output / 'result.json', {'protocol_sha256': sha256(protocol_path), 'intact': before, 'thinned': after, 'gates': gates, 'all_gates_pass': all(gates.values()), 'scope': protocol['limitation'], 'strict_qualification_pass': False, 'input_sha256': bound})


if __name__ == '__main__':
    main()
