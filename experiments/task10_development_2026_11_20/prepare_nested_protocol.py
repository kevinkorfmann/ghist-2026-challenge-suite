"""Freeze a future qualification design before any v2 development outcomes."""
import json
from pathlib import Path

import numpy as np
from nested_qualification_design import design_summary, large_nested_shape, nested_shapes

from scripts.develop_task10_detected_ibd import derived_seed, save_json, sha256


def main():
    directory = Path('/Users/kevinkorfmann/Documents/ChatGPT/ghist/experiments/task10_development_2026_11_20')
    cohorts = []
    for index in range(10):
        shape = index // 2 if index < 8 else 4
        leaves = (nested_shapes() + [large_nested_shape()])[shape]
        row = design_summary(leaves, f'nested{index:03d}', 229 if index >= 8 else None)
        seed = int(np.random.SeedSequence([2026192010, index]).generate_state(1, dtype=np.uint32)[0])
        row.update(index=index, shape=shape, seed=seed, stage='UNRUN_INDEPENDENT_QUALIFICATION',
            derived_seeds={'meioses': derived_seed(seed, 0, 3), 'chromosomes': {
                str(c): {'ancestry': derived_seed(seed, c, 1), 'mutation': derived_seed(seed, c, 2),
                    'dropout_contract_only': derived_seed(seed, c, 4)} for c in range(1, 23)}})
        cohorts.append(row)
    previous = json.loads((directory / 'detected_panel_protocol_v2.json').read_text())['cohorts']
    old_seeds = [row['seed'] for row in previous] + [2026132010, 2026133010, 2026134010, 2026135010,
        2026142010, 2026143010, 2026144010, 2026145010]
    streams = [derived_seed(seed, c, s) for seed in old_seeds + [row['seed'] for row in cohorts]
        for c in range(23) for s in range(1, 5)]
    if len(streams) != len(set(streams)) or 2026180010 in streams:
        raise ValueError('Random stream collision')
    numerical = {'exact_accuracy_min': .80, 'macro_recall_min': .70, 'minimum_label_recall_min': .50,
        'exact_accuracy_gain_over_original_baseline_min': .10}
    protocol = {
        'status': 'PROSPECTIVE_DESIGN_ONLY_NOT_DISPATCHED',
        'candidate': 'detected-multiway-family-forest-v2',
        'candidate_source_manifest_sha256': '7d4c00c9cfd9d197fce722e5de566078336ff0e6c2d7a39a84504e28a3bd9489',
        'design_source_sha256': sha256(directory / 'nested_qualification_design.py'),
        'launch_requires': ['All frozen v2 primary context development gates pass',
            'Model and scalar-temperature SHA binding before new genotype generation',
            '229-person exact-segment preflight completes within fixed feature time limit'],
        'cohorts': cohorts,
        'gates': {'small_nested_aggregate': {**numerical, 'no_individual_cohort_regression_from_baseline': True},
            'each_229_person_family': {**numerical, 'all_26106_pairs_required': True,
                'whole_pipeline_wall_minutes_max': 90, 'peak_memory_GiB_max': 44},
            'all': {'complete_labels': 14, 'missing_extra_duplicate_pairs': 0,
                'finite_normalized_probabilities': True, 'calibrated_toplabel_ECE_10bins_max': .10}},
        'comparators': 'Original frozen baseline and v2 pair-only ablation on identical complete cohorts; no retraining or recalibration. Context is the sole qualifying candidate.',
        'genotype_model': 'Same frozen v2 founder and detector configuration; independent founders per family; 22 autosomes; genuinely phased complete calls; no reads or phase errors. Exact segments retained only for diagnosis.',
        'domain': 'Outbred nested pedigrees, unrelated spouses, accepted two-MRCA labels. Does not qualify inbreeding, multiple ancestral paths, ancestry shifts, sequencing errors, phase errors or missing genotypes.',
        'missingness_contract': {'rates': [.05, .20], 'expected': 'Reject missing alleles without publishing prepared VCF. Rejection is a contract check, not robustness. No imputation or site deletion.',
            'mask_rule': 'Independent whole-genotype masks from per-chromosome stage-4 seed, canonical sorted sample order and genomic record order. Reuse same uniforms at both thresholds.'},
        'resources_proposed': {'small_families': 8, 'small_limit_minutes_each': 45, 'small_cpus_each': 4,
            'large_families': 2, 'large_limit_minutes_each': 90, 'large_cpus_each': 8,
            'hard_CPU_hours': 48, 'small_concurrency_max': 4, 'large_concurrency_max': 2,
            'estimate': 'Must be estimated from measured development runtime before dispatch; these are caps.'},
        'seed_audit': {'old_and_new_cartesian_streams': len(streams), 'unique': True, 'old_confirmation_unused': True},
        'freeze_rule': 'Roster and gates precede observation of v2 development outcomes. Runner, source manifest and model hashes must additionally freeze before dispatch. Failure cannot tune this candidate and then be relabeled confirmation.',
    }
    save_json(directory / 'nested_qualification_protocol.json', protocol)
    (directory / 'nested_qualification_structures.def').write_text(''.join(row['definition'] for row in cohorts[::2]))
    print(json.dumps({'protocol_sha256': sha256(directory / 'nested_qualification_protocol.json'),
        'cohorts': len(cohorts), 'pairs': sum(row['pairs'] for row in cohorts), 'streams': len(streams)}))


if __name__ == '__main__':
    main()
