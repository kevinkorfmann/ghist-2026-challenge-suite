from scripts.consensus_task11_dog_pedigrees import (
    bonsai_label,
    choose_from_probabilities,
)


def test_accepted_dog_topology_mapping():
    row = {"relationship_tuple": [2, 1, 1]}
    assert bonsai_label(row) == ("half_avuncular", 3)


def test_unsupported_topology_still_constrains_meioses():
    row = {"relationship_tuple": [0, 4, 1]}
    assert bonsai_label(row) == (None, 4)
    label, confidence = choose_from_probabilities(
        {
            "siblings": 0.50,
            "avuncular_once_removed": 0.20,
            "full_1st_cousin": 0.10,
            "half_1st_cousin": 0.05,
        },
        4,
    )
    assert label == "avuncular_once_removed"
    assert confidence > 0.5
