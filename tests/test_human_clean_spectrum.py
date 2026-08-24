from scripts.fit_human_clean_spectrum import spectrum_feature_names


def test_spectrum_features_keep_counts_and_bins_without_ordered_shape():
    names = spectrum_feature_names(5)
    assert "segments_ge_5" in names
    assert "length_bin_1" in names
    assert "max_cm" not in names
    assert not any(name.startswith("ordered_length_") for name in names)
