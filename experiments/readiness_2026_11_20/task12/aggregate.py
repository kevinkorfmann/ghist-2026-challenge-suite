"""Summarize fixed holdout records; no model selection or fitting."""
import hashlib
import json
from pathlib import Path
from statistics import mean

path = Path(__file__).with_name('holdout.json')
report = json.loads(path.read_text())
runs = report['runs']
seeds = sorted({row['seed'] for row in runs})
metrics = ('accuracy', 'balanced_accuracy', 'meioses_accuracy',
           'ece_10_equal_width_bins', 'multiclass_brier_sum', 'negative_log_likelihood')
summary = {
    'holdout_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
    'independent_seeds': seeds,
    'regimes_per_seed': {seed: sum(row['seed'] == seed for row in runs) for seed in seeds},
    'paired_comparisons': len(runs),
    'validated_csv_outputs': 2 * len(runs),
    'unique_simulated_input_pairs': sum(next(row['models']['candidate']['pairs'] for row in runs if row['seed'] == seed) for seed in seeds),
    'all_csv_formats_valid': all(row['models'][model]['validation']['format_valid'] for row in runs for model in ('historical', 'candidate')),
    'all_simulated_input_pairs_covered': all(row['models'][model]['validation']['input_pair_coverage_verified'] for row in runs for model in ('historical', 'candidate')),
    'models': {model: {metric: mean(row['models'][model][metric] for row in runs) for metric in metrics} for model in ('historical', 'candidate')},
    'per_seed': [
        {'seed': seed, 'pairs': next(row['models']['candidate']['pairs'] for row in runs if row['seed'] == seed),
         'models': {model: {metric: mean(row['models'][model][metric] for row in runs if row['seed'] == seed) for metric in metrics} for model in ('historical', 'candidate')}}
        for seed in seeds
    ],
    'mean_exact_gain': report['mean_exact_gain'],
    'mean_exact_gain_by_seed': report['mean_exact_gain_by_seed'],
    'mean_meioses_gain': report['mean_meioses_gain'],
    'minimum_regime_exact_gain': min(row['models']['candidate']['accuracy'] - row['models']['historical']['accuracy'] for row in runs),
    'gate': report['gate'],
    'averaging': 'Equal weight per fixed regime and independent seed; regimes within a seed are repeated measurements',
    'confidence_is_uncalibrated': True,
    'component_only': True,
    'full_vcf_pipeline_validated': False,
    'candidate_promoted': False,
}
Path(__file__).with_name('aggregate.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
