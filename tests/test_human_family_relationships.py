from scripts.fit_human_family_relationships import (
    feature_names,
    relationship_label,
    summarize_raw_lengths,
)


def test_human_relationship_label_parsing():
    assert (
        relationship_label("full_1st_cousin_twice_removed9_g5-b1-i1")
        == "full_1st_cousin_twice_removed"
    )


def test_human_raw_length_summary():
    result = summarize_raw_lengths([("1", 7.0), ("2", 25.0)], genome_cm=100.0)
    assert result["segments_ge_5"] == 2
    assert result["segments_ge_10"] == 1
    assert result["relatedness_ge_20"] == 0.125
    names = feature_names(20)
    assert names[:7] == [
        "max_cm", "segments_ge_20", "total_cm_ge_20", "relatedness_ge_20",
        "segments_ge_40", "total_cm_ge_40", "relatedness_ge_40",
    ]
    assert "length_bin_3" in names
    assert names[-1] == "ordered_length_20"
