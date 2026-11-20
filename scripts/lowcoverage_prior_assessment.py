"""Summarize retained synthetic evidence before new pragmatic cap500 outcomes."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import mean, median

from ghist_lowcoverage.pipeline import sha256


def summarize(report):
    rows = report['cohorts']
    keys = ['generations_since_split', 'pop1_population_size', 'pop2_population_size']
    primary = [abs(row['relative_errors']['seed1'][key]) for row in rows.values() for key in keys]
    return {
        'cohort_count': len(rows),
        'primary_target_count': len(primary),
        'primary_median_absolute_relative_error': median(primary),
        'all_runs_maximum_absolute_relative_error': max(abs(value) for row in rows.values()
            for estimates in row['relative_errors'].values() for value in estimates.values()),
        'primary_signed_error_by_target': {key: {
            'mean': mean(row['relative_errors']['seed1'][key] for row in rows.values()),
            'minimum': min(row['relative_errors']['seed1'][key] for row in rows.values()),
            'maximum': max(row['relative_errors']['seed1'][key] for row in rows.values())} for key in keys},
        'maximum_target_sensitivity': max(abs(value) for row in rows.values()
            for comparison in row['comparisons'].values() for value in comparison['target_relative_change'].values()),
        'maximum_frequency_L1': max(comparison['frequency_L1'] for row in rows.values()
            for comparison in row['comparisons'].values()),
        'maximum_segregating_relative_L1': max(comparison['segregating_relative_L1'] for row in rows.values()
            for comparison in row['comparisons'].values()),
        'minimum_callable_fraction': min(value for row in rows.values() for value in row['callable_fraction'].values()),
        'failed_gates_by_case': {name: row['failed_gates'] for name, row in rows.items()},
        'noise_stratum_primary_median_absolute_relative_errors': report['stratum_medians'],
        'all_registered_gates_passed': report['all_registered_gates_passed'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120' and not args.output.exists()
    names = ['task14_confirmation_gate_audit_v1.json', 'task14_stabilized_development_gate_audit_v1.json',
             'task13_refinement_gate_audit_v1.json']
    original, stabilized, bottleneck = [json.loads((args.root / name).read_text()) for name in names]
    result = {
        'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Previously inspected synthetic evidence only; no cap500 pragmatic results or challenge data read.',
        'sources_sha256': {name: sha256(args.root / name) for name in names},
        'summarizer_sha256': sha256(Path(__file__)),
        'task14_original_fixed_window100_seed1': summarize(original),
        'task14_fixed_window500_cap200_seed1': summarize(stabilized),
        'task13_fullbatch_refinement': {
            'all_registered_diagnostic_gates_passed': bottleneck['all_registered_diagnostic_gates_passed'],
            'failed_gates_by_initialization': {name: row['failed_gates'] for name, row in bottleneck['initializations'].items()},
            'age_error_range': [min(row['relative_errors'][run]['generations'] for row in bottleneck['initializations'].values()
                for run in ['refined', 'continued']), max(row['relative_errors'][run]['generations']
                for row in bottleneck['initializations'].values() for run in ['refined', 'continued'])],
            'fraction_error_range': [min(row['relative_errors'][run]['post_decline_fraction'] for row in bottleneck['initializations'].values()
                for run in ['refined', 'continued']), max(row['relative_errors'][run]['post_decline_fraction']
                for row in bottleneck['initializations'].values() for run in ['refined', 'continued'])],
        },
        'assessment': {
            'strictly_qualified_candidates': [],
            'broadest_measured_task14_candidate': 'Original fixed window100 seed1 has six prespecified cohorts and <=11.34% error across all runs; 1.08-3.93% target sensitivity fails strict 1% but remains within the separately authorized pragmatic 5% rule.',
            'cap200_tradeoff': 'Full window reduces seed sensitivity below 0.04% in the two regime0 cases; five of six runs hit cap200 above 1e-6 and no new regime1/2 validation was performed.',
            'cap500_plan': 'Separate prospectively frozen same engine, two existing development cases, no fresh confirmation. Root reviews all unchanged strict gates and separate pragmatic usefulness criteria before any final inference.',
            'task13': 'No current supported recovery claim for decline fraction: full-batch refinement agrees across optimizers yet overestimates fraction by about51.7% in the single full-width development case, and stopping/spectrum stability still fail.',
            'bias_interpretation': 'Signed error summaries describe these small synthetic designs only, not calibrated bias corrections or universal uncertainty bounds.',
            'estimate_selection': 'Fixed seed1 primary; no seed selection, truth-based correction, or final-data tuning.'},
        'scientifically_ready': False,
        'final_data_access': False,
        'submitted': False,
    }
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
