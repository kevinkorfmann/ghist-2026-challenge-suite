from ghist_relatedness.labels import accepted_label, canonical_tuple, meioses_count


def test_human_removed_relationships() -> None:
    assert accepted_label((3, 1, 2), "human") == "avuncular_once_removed"
    assert accepted_label((5, 2, 2), "human") == "full_1st_cousin_thrice_removed"


def test_dog_half_relationships() -> None:
    assert accepted_label((2, 1, 1), "dog") == "half_avuncular"
    assert accepted_label((2, 2, 1), "dog") == "half_1st_cousin"


def test_meioses_count() -> None:
    assert canonical_tuple((4, 2, 2)) == (2, 4, 2)
    assert meioses_count((4, 2, 2)) == 6
