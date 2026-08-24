from scripts.consensus_human_family_pedigree import (
    choose_relation,
    topology_evidence,
)


def test_human_topology_mapping():
    assert topology_evidence({"relationship_tuple": [2, 3, 2]}) == (
        "full_1st_cousin_once_removed",
        5,
    )


def test_primary_and_guarded_topology_policies():
    probabilities = {
        "siblings": 0.7,
        "avuncular_once_removed": 0.3,
    }
    assert choose_relation(
        probabilities, "avuncular_once_removed", 4, None
    )[0] == "avuncular_once_removed"
    assert choose_relation(probabilities, "avuncular_once_removed", 4, 2.0)[0] == "siblings"
    assert (
        choose_relation(probabilities, "avuncular_once_removed", 4, 3.0)[0]
        == "avuncular_once_removed"
    )
