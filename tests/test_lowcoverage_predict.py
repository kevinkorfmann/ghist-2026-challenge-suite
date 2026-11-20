import copy

import pytest

from scripts.lowcoverage_predict_window import numerical_gates


def valid_comparison():
    labels = ['seed1', 'seed2', 'continued']
    return {'estimates': dict.fromkeys(labels, {}),
            'callable_fraction': dict.fromkeys(labels, .99),
            'early_stop_observed': dict.fromkeys(labels, True),
            'optimizer_agreement': {name: {'passed': True} for name in labels},
            'comparisons': {name: {'frequency_L1': .0009, 'segregating_relative_L1': .049,
                                   'target_relative_change': dict.fromkeys([
                                       'pop1_population_size', 'pop2_population_size', 'generations_since_split'], -.009)}
                            for name in ['seed2', 'continued']}}


def test_all_frozen_numerical_gates_required():
    comparison = valid_comparison()
    assert all(numerical_gates(comparison).values())
    fields = [('frequency_L1', .001, 'frequency_L1'),
              ('segregating_relative_L1', .05, 'segregating_relative_L1')]
    for field, value, gate in fields:
        changed = copy.deepcopy(comparison)
        changed['comparisons']['seed2'][field] = value
        assert not numerical_gates(changed)[gate]
    comparison['comparisons']['continued']['target_relative_change']['generations_since_split'] = -.01
    assert not numerical_gates(comparison)['target_stability']


@pytest.mark.parametrize('field,value,gate', [
    ('callable_fraction', .989999, 'exposure'),
    ('early_stop_observed', False, 'observed_stopping'),
    ('optimizer_agreement', {'passed': False}, 'two_optimizer_agreement')])
def test_alternate_seed_failure_blocks_draft(field, value, gate):
    comparison = valid_comparison()
    comparison[field]['seed2'] = value
    assert not numerical_gates(comparison)[gate]


def test_missing_estimator_run_rejected():
    comparison = valid_comparison()
    del comparison['estimates']['seed2']
    with pytest.raises(ValueError, match='all three'):
        numerical_gates(comparison)
