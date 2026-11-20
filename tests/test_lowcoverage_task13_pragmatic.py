"""Root authorization is checked before any explicit final input is accessed."""
import pytest

from scripts.lowcoverage_predict_task13_pragmatic import FREEZE, check_decision


def decision():
    return {'authorize_final_inference': True, 'task': 'bottleneck',
            'strict_qualification': False, 'freeze_sha256': FREEZE,
            'wrapper_sha256': 'wrapper', 'primary_estimate': 'continued',
            'development_selection_disclosed': True}


def test_exact_decision():
    check_decision(decision(), 'wrapper')


@pytest.mark.parametrize('field,value', [('primary_estimate', 'seed1'),
    ('task', 'split_isolation'), ('development_selection_disclosed', False),
    ('strict_qualification', True), ('authorize_final_inference', False)])
def test_incompatible_decisions_rejected(field, value):
    record = decision()
    record[field] = value
    with pytest.raises(ValueError):
        check_decision(record, 'wrapper')
