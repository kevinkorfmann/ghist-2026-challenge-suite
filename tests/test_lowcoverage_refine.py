import numpy as np
import pytest

from scripts.lowcoverage_refine_window import initialization_counts


def test_initialization_only_floors_allowed_zero_bins_without_refolding():
    values = np.zeros(41)
    values[0], values[1], values[20] = 9990, 9, 1
    original = values.copy()
    result = initialization_counts(values)
    np.testing.assert_equal(values, original)
    np.testing.assert_equal(result[21:], 0)
    np.testing.assert_equal(result[[0, 1, 20]], [9990, 9, 1])
    assert np.all(result[2:20] == 1e-8)


@pytest.mark.parametrize('bad', [np.ones(40), np.ones(41), np.full(41, np.nan), np.full(41, -1)])
def test_initialization_rejects_invalid_or_unfolded_arrays(bad):
    with pytest.raises(ValueError):
        initialization_counts(bad)
