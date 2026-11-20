"""Paired exact/detected feature diagnostic; never fits or confirms a candidate."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np

from scripts.develop_task10_detected_ibd import (
    genome_length,
    phase_pair_features,
    read_map,
    save_json,
    score_predictions,
    sha256,
    verify_inputs,
)


def decompressed_hash(path):
    digest = hashlib.sha256()
    with gzip.open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment-root', type=Path, required=True)
    parser.add_argument('--asset-root', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use Slurm for genome-scale diagnostics')
    root, assets = args.experiment_root.resolve(), args.asset_root.resolve()
    cohort = root / f'phase_development_cohort_{args.index}'
    output = root / f'paired_phase_diagnostic/cohort_{args.index}'
    output.mkdir(parents=True, exist_ok=False)
    paths = [cohort / name for name in ['commands.json', 'identifiers.json', 'truth.json',
        'family.vcf.gz', 'family-everyone.fam', 'founders/founders.vcf.gz']]
    model_path = root / 'phase_training/phase_model.joblib'
    control_path = root / 'phase_training/phase_free_model.joblib'
    detected_path = root / f'phase_predictions/cohort_{args.index}/phase/features.json'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    paths.extend([model_path, control_path, detected_path, map_path, Path(__file__).resolve()])
    hashes = {str(path): sha256(path) for path in paths}
    original = json.loads((cohort / 'commands.json').read_text())[0]['command']
    command = list(original)
    command[command.index('-o') + 1] = str(output / 'family')
    command.append('--seg_haps')
    save_json(output / 'protocol.json', {'stage': 'PAIRED_DEVELOPMENT_DIAGNOSTIC', 'command': command,
        'hashes': hashes, 'feature_swaps': ['raw_', 'state_', 'phase_', 'ibd2_'], 'submission_ready': False})
    with (output / 'pedsim.log').open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    genotype_hash = decompressed_hash(cohort / 'family.vcf.gz')
    if decompressed_hash(output / 'family.vcf.gz') != genotype_hash:
        raise ValueError('Haplotype-resolved rerun changed genotypes; paired diagnostic invalid')
    if sha256(output / 'family-everyone.fam') != sha256(cohort / 'family-everyone.fam'):
        raise ValueError('Haplotype-resolved rerun changed parentage')
    identifiers = json.loads((cohort / 'identifiers.json').read_text())
    grouped = defaultdict(lambda: defaultdict(list))
    with (output / 'family.seg').open() as handle:
        for line in handle:
            a, h1, b, h2, chrom, start, end, state, *_ = line.split()
            if a not in identifiers or b not in identifiers or a == b:
                continue
            if state != 'IBD1':
                raise ValueError('Unexpected phased ancestral segment state')
            a, b, h1, h2 = identifiers[a], identifiers[b], int(h1) + 1, int(h2) + 1
            if a > b:
                a, b, h1, h2 = b, a, h2, h1
            if int(end) > int(start):
                grouped[(a, b)][chrom].append((int(start), int(end), h1, h2))
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    exact = [{'indID1': a, 'indID2': b, **phase_pair_features(grouped.get((a, b), {}), genetic_map, genome_cm)}
        for a, b in itertools.combinations(sorted(identifiers.values()), 2)]
    detected = json.loads(detected_path.read_text())
    if [(row['indID1'], row['indID2']) for row in exact] != [(row['indID1'], row['indID2']) for row in detected]:
        raise ValueError('Paired feature order differs')
    save_json(output / 'exact_features.json', exact)
    truth = json.loads((cohort / 'truth.json').read_text())
    scores = {}
    for variant, path in [('phase', model_path), ('phase_free', control_path)]:
        bundle = joblib.load(path)
        classifier, names = bundle['model'], bundle['feature_names']
        policies = {'exact': names, 'detected': [], **{f'swap_{prefix}': [name for name in names if name.startswith(prefix)]
            for prefix in ['raw_', 'state_', 'phase_', 'ibd2_']}}
        for policy, swapped in policies.items():
            swapped = set(swapped)
            matrix = [[exact_row[name] if name in swapped else detected_row[name] for name in names]
                for exact_row, detected_row in zip(exact, detected, strict=True)]
            probabilities = classifier.predict_proba(matrix)
            predictions = [{'indID1': row['indID1'], 'indID2': row['indID2'],
                'relation': str(classifier.classes_[int(np.argmax(values))]),
                'probabilities': dict(zip(map(str, classifier.classes_), map(float, values), strict=True))}
                for row, values in zip(exact, probabilities, strict=True)]
            scores[f'{variant}_{policy}'] = score_predictions(truth, predictions)
            save_json(output / f'{variant}_{policy}_predictions.json', predictions)
        save_json(output / f'{variant}_importance.json', dict(sorted(zip(names, map(float, classifier.feature_importances_), strict=True),
            key=lambda entry: -entry[1])))
    verify_inputs(hashes)
    save_json(output / 'result.json', {'scores': scores, 'identical_uncompressed_genotype_sha256': genotype_hash,
        'identical_parentage_sha256': sha256(cohort / 'family-everyone.fam'),
        'submission_ready': False, 'submission_count': 0})


if __name__ == '__main__':
    main()
