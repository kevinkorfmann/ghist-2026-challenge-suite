"""Frozen pair-feature parity and 229-person inference runtime, not qualification."""
from __future__ import annotations

import itertools
import json
import os
import resource
import time
from pathlib import Path

import joblib
import numpy as np
from structural_pair_candidate import infer, pair_features, pair_predictions

from scripts.develop_task10_detected_ibd import (
    genome_length,
    read_exact_phase,
    read_haplotype_segments,
    read_map,
    save_json,
    sha256,
    verify_inputs,
)


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Use Slurm')
    root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
    source = Path(__file__).resolve().parent
    frozen = root / 'detected-panel-source'
    manifest = json.loads((source / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    verify_inputs(json.loads((frozen / 'source_manifest.json').read_text()))
    protocol = json.loads((source / 'structural_pair_runtime_protocol.json').read_text())
    evaluation = root / 'detected_panel_v2/model_evaluation'
    model_path = evaluation / 'pair_only_model.joblib'
    freeze_path = evaluation / 'model_calibration_freeze.json'
    freeze = json.loads(freeze_path.read_text())
    if sha256(model_path) != freeze['model_hashes']['pair_only']:
        raise ValueError('Frozen model changed')
    bundle = joblib.load(model_path)
    temperature = freeze['temperatures']['pair_only']['temperature']
    map_path = frozen / 'assets/external/ped-sim-map/refined_average.plink.map'
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    output = root / 'structural_pair_runtime_v1'
    output.mkdir(exist_ok=False)
    inputs = {str(path): sha256(path) for path in [model_path, freeze_path, map_path]}
    checks = []
    for index in protocol['training_parity_cohorts']:
        cohort = root / f'detected_panel_v2/cohort_{index:03d}'
        paths = [cohort / name for name in ['identifiers.json', 'features.json', 'detected.ibd.gz']]
        inputs.update({str(path): sha256(path) for path in paths})
        samples = sorted(json.loads(paths[0].read_text()).values())
        grouped, _ = read_haplotype_segments(paths[2])
        segments = {}
        for (a, b, chrom), rows in grouped.items():
            segments.setdefault((a, b), {})[chrom] = rows
        old = json.loads(paths[1].read_text())
        new = pair_features(samples, segments, genetic_map, genome_cm)
        matrices = [np.asarray([[row[name] for name in bundle['feature_names']] for row in rows], dtype=np.float32)
            for rows in [old, new]]
        if not np.array_equal(*matrices):
            raise ValueError('Direct pair-only feature inputs differ')
        original, direct = [pair_predictions(rows, bundle, temperature) for rows in [old, new]]
        if original != direct:
            raise ValueError('Direct pair-only model probabilities differ')
        checks.append({'training_cohort': index, 'pairs': len(new),
            'float32_inputs_identical': True, 'serial_probabilities_identical': True})
    # Verify integration and serial reduction preserve all reused development labels.
    previous_path = root / 'integer_edge_decoder_v1/pair_only_scorer_input.json'
    inputs[str(previous_path)] = sha256(previous_path)
    previous = json.loads(previous_path.read_text())
    integrations = []
    for cohort in previous['cohorts']:
        index = int(cohort['cohort_id'].split('_')[-1])
        path = root / f'detected_panel_v2/cohort_{index:03d}/features.json'
        inputs[str(path)] = sha256(path)
        predictions = pair_predictions(json.loads(path.read_text()), bundle, temperature)
        rows, _ = infer(predictions)
        if rows != cohort['predictions']:
            raise ValueError('Structural labels changed with integrated serial prediction')
        integrations.append({'cohort_id': cohort['cohort_id'], 'pairs': len(rows), 'labels_identical': True})
    save_json(output / 'parity.json', {'training': checks, 'reused_development_integration': integrations,
        'input_sha256': inputs, 'truth_used_for_inference': False})
    data = root / 'scaling_preflight'
    paths = [data / name for name in ['roster.json', 'printed_roster.json', 'profile.seg']]
    inputs.update({str(path): sha256(path) for path in paths})
    started = time.monotonic()
    samples = sorted(json.loads(paths[0].read_text()))
    if len(samples) != protocol['people']:
        raise ValueError('Unexpected runtime roster')
    expected = dict.fromkeys(itertools.combinations(sorted(json.loads(paths[1].read_text())), 2))
    segments = read_exact_phase(paths[2], expected)
    selected = set(samples)
    segments = {pair: value for pair, value in segments.items() if set(pair) <= selected}
    loaded = time.monotonic()
    features = pair_features(samples, segments, genetic_map, genome_cm)
    featured = time.monotonic()
    predictions = pair_predictions(features, bundle, temperature)
    predicted = time.monotonic()
    try:
        rows, tree = infer(predictions)
    except Exception as error:
        save_json(output / 'result.json', {'state': 'INFERENCE_FAILED', 'error_type': type(error).__name__,
            'error': str(error), 'people': len(samples), 'pairs': len(predictions),
            'total_seconds': time.monotonic() - started, 'biological_qualification': False,
            'input_sha256': inputs, 'protocol_sha256': sha256(source / 'structural_pair_runtime_protocol.json')})
        raise
    inferred = time.monotonic()
    save_json(output / 'predictions.json', rows)
    save_json(output / 'tree.json', tree)
    verify_inputs(inputs)
    verify_inputs(manifest)
    save_json(output / 'result.json', {'state': 'COMPLETE_RUNTIME_PREFLIGHT', 'people': len(samples),
        'pairs': len(rows), 'load_seconds': loaded - started, 'feature_seconds': featured - loaded,
        'pair_predict_seconds': predicted - featured, 'tree_seconds': inferred - predicted,
        'total_seconds': inferred - started, 'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'biological_qualification': False, 'data': 'Retained exact ancestral development segments; no new genomes',
        'input_sha256': inputs, 'source_manifest_sha256': sha256(source / 'source_manifest.json'),
        'protocol_sha256': sha256(source / 'structural_pair_runtime_protocol.json')})


if __name__ == '__main__':
    main()
