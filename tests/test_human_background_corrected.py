from scripts.fit_human_background_corrected import correct_records
from scripts.fit_human_family_relationships import summarize_raw_lengths


def test_corrected_domain_transform_reduces_added_background():
    reference = []
    observed = []
    for index in range(8):
        clean = {"indID1": f"a{index}", "indID2": f"b{index}", "label": "siblings"}
        clean.update(summarize_raw_lengths([("1", 60.0)], genome_cm=100.0))
        noisy = {"indID1": clean["indID1"], "indID2": clean["indID2"], "label": "siblings"}
        noisy.update(
            summarize_raw_lengths(
                [("1", 22.0), ("1", 28.0), ("1", 60.0)], genome_cm=100.0
            )
        )
        reference.append(clean)
        observed.append(noisy)
    corrected = correct_records(observed, reference, minimum_cm=20, genome_cm=100.0)
    assert corrected[0]["segments_ge_20"] == 1.0
    assert corrected[0]["ordered_length_1"] == 60.0
