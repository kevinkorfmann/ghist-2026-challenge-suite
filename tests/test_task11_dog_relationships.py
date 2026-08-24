from scripts.fit_task11_dog_relationships import (
    feature_names,
    merge_intervals,
    relationship_label,
    subtract_intervals,
    summarize_segments,
)


def test_dog_relationship_label_parsing():
    assert relationship_label("full_1st_cousin17_g3-b1-i1") == "full_1st_cousin"
    assert relationship_label("siblings2_g2-b2-i1") == "siblings"


def test_dog_segment_summary_separates_ibd_states():
    values = summarize_segments(
        [("1", 12.0, False), ("1", 8.0, True), ("2", 3.0, False)]
    )
    result = dict(zip(feature_names(), values))
    assert result["segments"] == 3
    assert result["chromosomes"] == 2
    assert result["total_ibd1"] == 15
    assert result["total_ibd2"] == 8
    assert result["total_sharing"] == 31


def test_hbd_mask_interval_operations():
    assert merge_intervals([(8, 12), (1, 4), (3, 6)]) == [(1, 6), (8, 12)]
    assert subtract_intervals(0, 20, [(3, 6), (10, 14)]) == [
        (0, 3),
        (6, 10),
        (14, 20),
    ]
