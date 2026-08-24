from scripts.fit_human_family_relationships import (
    feature_names,
    individual_background_correct,
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


def test_individual_background_correction_preserves_reference_level():
    observed = []
    reference = []
    effects = {"a": 0.01, "b": -0.01, "c": 0.0}
    for left, right in (("a", "b"), ("a", "c"), ("b", "c")):
        base = 0.02 + effects[left] + effects[right]
        row = {"indID1": left, "indID2": right, "max_cm": 20.0}
        ref = {"indID1": left, "indID2": right, "max_cm": 20.0}
        for threshold in (5, 10, 20, 40):
            row[f"segments_ge_{threshold}"] = base * 100
            row[f"relatedness_ge_{threshold}"] = base
            row[f"total_cm_ge_{threshold}"] = base * 200
            ref[f"segments_ge_{threshold}"] = 1.0
            ref[f"relatedness_ge_{threshold}"] = 0.01
            ref[f"total_cm_ge_{threshold}"] = 2.0
        observed.append(row)
        reference.append(ref)
    individual_background_correct(observed, reference, genome_cm=100.0)
    corrected = [row["relatedness_ge_20"] for row in observed]
    assert max(corrected) <= 0.011
