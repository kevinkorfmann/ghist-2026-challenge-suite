import json
from pathlib import Path
import pytest
import structural_fallback_policy as policy


def fixture(name):
    return json.loads((Path(__file__).parent/name).read_text())


def test_actual_no_root_fixture_returns_all_frozen_pair_rows_and_order_invariance():
    rows=fixture('fallback_no_root_fixture.json')
    with pytest.raises(ValueError,match='Initial NJ tree has no accepted root'):
        policy.structural_infer(rows)
    prediction,tree=policy.infer_with_fallback(rows)
    assert prediction==policy.canonical_complete(rows)[0]
    reverse=[dict(row,indID1=row['indID2'],indID2=row['indID1']) for row in reversed(rows)]
    assert policy.infer_with_fallback(reverse)==(prediction,tree)
    assert len(prediction)==120 and tree['pragmatic_policy']['used_fallback']
    assert tree['pragmatic_policy']['pairs']==120 and len(tree['pragmatic_policy']['sample_order'])==16


def test_actual_valid_fixture_preserves_exact_structural_rows():
    rows=fixture('fallback_valid_fixture.json')
    original,original_tree=policy.structural_infer(rows)
    prediction,tree=policy.infer_with_fallback(rows)
    assert prediction==policy.canonical_complete(original)[0]
    assert not tree['pragmatic_policy']['used_fallback']
    assert {k:v for k,v in tree.items() if k!='pragmatic_policy'}==original_tree


@pytest.mark.parametrize('error',[ValueError(policy.NO_ROOT+' extra'),ValueError('bad model'),RuntimeError(policy.NO_ROOT)])
def test_any_other_exception_is_not_fallback(monkeypatch,error):
    def fail(rows):
        raise error
    monkeypatch.setattr(policy,'structural_infer',fail)
    with pytest.raises(type(error),match=str(error)):
        policy.infer_with_fallback(fixture('fallback_no_root_fixture.json'))


def test_missing_or_duplicate_pairs_rejected_before_inference():
    rows=fixture('fallback_no_root_fixture.json')
    for broken in [rows[:-1],rows+[rows[0]]]:
        with pytest.raises(ValueError,match='complete unique'):
            policy.infer_with_fallback(broken)
