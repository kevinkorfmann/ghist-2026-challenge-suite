"""Strict count/shape checks for the experimental independent SFS estimator."""
import numpy as np
import pytest

from ghist_lowcoverage.winsfs import read_winsfs


def test_read_count_spectrum_keeps_invariant_exposure(tmp_path):
    path = tmp_path / 'sfs'
    path.write_text('#SHAPE=<3/3>\n900 10 20 30 40 0 0 0 0\n')
    spectrum = read_winsfs(path, [2, 2])
    assert spectrum.sum() == 1000 and spectrum[0] == 900


@pytest.mark.parametrize('text', ['1 2 3', '#SHAPE=<3>\n1 2', '#SHAPE=<2>\n1 2',
                                  '#SHAPE=<3>\n1 nan 2', '#SHAPE=<3>\n1 -1 2'])
def test_read_spectrum_rejects_malformed_shape_or_counts(tmp_path, text):
    path = tmp_path / 'sfs'
    path.write_text(text)
    with pytest.raises(ValueError):
        read_winsfs(path, [2])


def test_folded_window_spectrum_preserves_joint_fixed_differences(tmp_path):
    from ghist_lowcoverage.pipeline import folded_spectrum
    path = tmp_path / 'sfs'
    path.write_text('#SHAPE=<3/3>\n900 10 30 20 40 0 30 0 0\n')
    spectrum = folded_spectrum(read_winsfs(path, [2, 2]), [2, 2])
    assert not spectrum.mask[0, 2] and not spectrum.mask[2, 0]
    np.testing.assert_allclose([spectrum[0, 2], spectrum[2, 0]], [30, 30])
