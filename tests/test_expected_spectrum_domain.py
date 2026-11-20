import numpy as np
import pytest

from ghist_inference.validation import validate_expected_spectrum


@pytest.mark.parametrize("bad", [-0.01, float("nan"), float("inf"), 0])
def test_numerical_failure_cannot_remove_positive_observed_bin(bad):
    data = np.ma.array([0, 20, 5, 0], mask=[True, False, False, True])
    expected = np.ma.array([0, bad, 1, 0], mask=data.mask)
    with pytest.raises(ValueError, match="likelihood bins"):
        validate_expected_spectrum(expected, data)


def test_additional_model_mask_cannot_remove_observed_bin():
    data = np.ma.array([0, 20, 5, 0], mask=[True, False, False, True])
    model = np.ma.array([0, 2, 1, 0], mask=[True, True, False, True])
    with pytest.raises(ValueError):
        validate_expected_spectrum(model, data)


def test_monomorphic_mask_and_zero_data_zero_model_bin_are_allowed():
    data = np.ma.array([0, 20, 0, 0], mask=[True, False, False, True])
    model = np.ma.array([float("nan"), 2, 0, float("nan")], mask=data.mask)
    assert validate_expected_spectrum(model, data) is model
