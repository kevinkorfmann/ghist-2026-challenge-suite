import copy
import json
from itertools import combinations

import pytest

from scripts.consensus_task11_dog_pedigrees import (
    MEIOSES,
    bonsai_label,
    build_consensus,
    choose_from_probabilities,
    pair_key,
)
from scripts.consensus_task11_guarded_candidate import build_guarded_candidate


def probabilities(**values):
    return {label: values.get(label, 0.0) for label in MEIOSES}


def fixture(*, members=("a", "b"), family_id=0):
    pairs = list(combinations(members, 2))
    runs = [
        [
            {"indID1": a, "indID2": b, "probabilities": probabilities(siblings=1.0)}
            for a, b in pairs
        ]
        for _ in range(3)
    ]
    topology = {
        "families": [{"family": family_id, "members": list(members)}],
        "relationships": [
            {"indID1": a, "indID2": b, "relationship_tuple": [1, 1, 2]}
            for a, b in pairs
        ],
    }
    return runs, topology


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


def test_propagation_updates_confidence_and_supports_nonsequential_family_ids():
    runs, topology = fixture(members=("a", "b", "c"), family_id=17)
    topology["relationships"][1]["relationship_tuple"] = [1, 1, 1]
    topology["relationships"][2]["relationship_tuple"] = [2, 0, 1]
    rows, summary = build_consensus(runs, topology)
    propagated = next(row for row in rows if pair_key(row) == ("b", "c"))
    assert propagated["relation"] == "half_sibling"
    assert propagated["source"] == "full-sibling component propagation"
    assert propagated["confidence"] == 0.0
    assert summary["full_sibling_relation_vector_propagations"] == 1
    assert summary["final_labels_with_zero_pair_probability"] == 2
    assert summary["submission_ready"] is False


def test_relation_vector_refinement_updates_confidence():
    runs, topology = fixture()
    runs = runs[:2]
    runs[0][0]["probabilities"] = probabilities(siblings=0.8, half_sibling=0.2)
    runs[1][0]["probabilities"] = probabilities(siblings=0.3, half_sibling=0.7)
    topology["relationships"][0]["relationship_tuple"] = [2, 0, 1]
    rows, summary = build_consensus(runs, topology)
    assert rows[0]["relation"] == "half_sibling"
    assert rows[0]["source"] == "family relation-vector refinement"
    assert rows[0]["confidence"] == pytest.approx(0.45)
    assert summary["degree_two_relation_vector_refinements"] == 1


def test_confidence_is_unconditional_support_not_topology_renormalization():
    runs, topology = fixture()
    for run in runs:
        run[0]["probabilities"] = probabilities(siblings=0.9, avuncular=0.1)
    topology["relationships"][0]["relationship_tuple"] = [3, 0, 1]
    rows, _ = build_consensus(runs, topology)
    assert rows[0]["relation"] == "avuncular"
    assert rows[0]["confidence"] == pytest.approx(0.1)


def test_raw_sensitivity_is_not_hidden_by_topology_override():
    runs, topology = fixture()
    runs[1][0]["probabilities"] = probabilities(half_sibling=1.0)
    rows, summary = build_consensus(runs, topology)
    assert rows[0]["sensitivity_labels"] == ["siblings", "half_sibling", "siblings"]
    assert rows[0]["topology_conditioned_sensitivity_labels"] == ["siblings"] * 3
    assert summary["hbd_sensitivity_stability"] == 0.0
    assert summary["topology_conditioned_hbd_sensitivity_stability"] == 1.0
    assert summary["fallback_hbd_sensitivity_stability"] is None
    json.dumps(summary, allow_nan=False)


@pytest.mark.parametrize("meioses", [None, 2, 99])
def test_probability_ties_do_not_depend_on_mapping_order(meioses):
    tied = probabilities(siblings=0.5, half_sibling=0.5)
    reversed_keys = dict(reversed(list(tied.items())))
    assert choose_from_probabilities(tied, meioses) == ("siblings", 0.5)
    assert choose_from_probabilities(reversed_keys, meioses) == ("siblings", 0.5)


def test_consensus_tied_sensitivity_is_invariant_to_json_key_order():
    runs, topology = fixture()
    for run in runs:
        run[0]["probabilities"] = probabilities(siblings=0.5, half_sibling=0.5)
    expected = build_consensus(runs, topology)
    runs[1][0]["probabilities"] = dict(
        reversed(list(runs[1][0]["probabilities"].items()))
    )
    assert build_consensus(runs, topology) == expected
    assert expected[1]["hbd_sensitivity_stability"] == 1.0
    # The separate development policy must still abstain on the same exact tie.
    rows, _ = build_guarded_candidate(runs, topology)
    assert "pair_probability_tie" in rows[0]["abstention_reasons"]


@pytest.mark.parametrize(
    "relationship",
    [
        [float("nan"), 1, 2],
        [1, float("nan"), 2],
        [float("-inf"), 1, 2],
        [float("inf"), float("inf"), 2],
        [1, 1, None],
        [float("nan"), float("inf"), None],
        [float("-inf"), float("-inf"), None],
        [None, None, None],
        ["inf", "inf", None],
    ],
)
def test_corrupt_topology_cannot_trigger_unconstrained_consensus(relationship):
    runs, topology = fixture()
    topology["relationships"][0]["relationship_tuple"] = relationship
    with pytest.raises(ValueError, match="Bonsai"):
        build_consensus(runs, topology)


def test_explicit_bonsai_unrelated_sentinel_is_preserved():
    runs, topology = fixture()
    unrelated = [float("inf"), float("inf"), None]
    topology["relationships"][0]["relationship_tuple"] = unrelated
    assert bonsai_label({"relationship_tuple": unrelated}) == (None, None)
    rows, summary = build_consensus(runs, topology)
    assert rows[0]["relation"] == "siblings"
    assert rows[0]["topology_meioses"] is None
    assert summary["pair_model_fallbacks"] == 1
    json.dumps((rows, summary), allow_nan=False)


@pytest.mark.parametrize("component", ["prediction", "topology"])
def test_reversed_duplicate_pairs_are_rejected(component):
    runs, topology = fixture()
    rows = runs[0] if component == "prediction" else topology["relationships"]
    duplicate = dict(rows[0], indID1="b", indID2="a")
    rows.append(duplicate)
    with pytest.raises(ValueError, match="Duplicate unordered pair"):
        build_consensus(runs, topology)


@pytest.mark.parametrize("component", ["prediction", "topology", "partition"])
def test_missing_pairs_and_repeated_partition_samples_are_rejected(component):
    runs, topology = fixture(members=("a", "b", "c"))
    if component == "prediction":
        for run in runs:
            run.pop()
    elif component == "topology":
        topology["relationships"].pop()
    else:
        topology["families"].append({"family": 1, "members": ["a"]})
    with pytest.raises(ValueError, match="exactly cover|Duplicate sample"):
        build_consensus(runs, topology)


def test_different_hbd_pair_sets_are_rejected():
    runs, topology = fixture(members=("a", "b", "c"))
    runs[1].pop()
    with pytest.raises(ValueError, match="different pairs"):
        build_consensus(runs, topology)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -0.1, 1.1, True, "1.0"])
def test_invalid_probabilities_are_rejected(value):
    runs, topology = fixture()
    runs[0][0]["probabilities"]["siblings"] = value
    with pytest.raises(ValueError, match="finite numbers"):
        build_consensus(runs, topology)


@pytest.mark.parametrize(
    "invalid",
    [probabilities(), {"siblings": 1.0}, {**probabilities(siblings=1), "unknown": 0}],
)
def test_incomplete_or_unnormalized_probabilities_are_rejected(invalid):
    runs, topology = fixture()
    runs[0][0]["probabilities"] = invalid
    with pytest.raises(ValueError, match="sum to one|every accepted|unknown labels"):
        build_consensus(runs, topology)


def test_zero_mass_meiosis_constraint_cannot_invent_a_label():
    with pytest.raises(ValueError, match="No probability support"):
        choose_from_probabilities(probabilities(siblings=1), 4)


@pytest.mark.parametrize(
    "relationship",
    [[1.9, 1, 2], [-1, 1, 2], [1, 1, float("nan")], [True, 1, 2], [1, 1]],
)
def test_invalid_topology_is_not_silently_coerced(relationship):
    with pytest.raises(ValueError):
        bonsai_label({"relationship_tuple": relationship})


def test_guarded_candidate_preserves_raw_proposal_and_abstains_on_conflict():
    runs, topology = fixture()
    topology["relationships"][0]["relationship_tuple"] = [1, 1, 1]
    original = copy.deepcopy((runs, topology))
    rows, summary = build_guarded_candidate(runs, topology)
    assert rows[0]["proposed_relation"] == "siblings"
    assert rows[0]["relation"] is None
    assert rows[0]["abstention_reasons"] == ["topology_disagreement"]
    assert summary["provisional_pairs"] == 0
    assert summary["submission_ready"] is False
    assert (runs, topology) == original


def test_guarded_candidate_agreement_is_only_provisional():
    runs, topology = fixture()
    rows, summary = build_guarded_candidate(runs, topology)
    assert rows[0]["relation"] == "siblings"
    assert rows[0]["decision"] == "provisional"
    assert summary["submission_ready"] is False
    assert summary["submission_csv_exported"] is False


@pytest.mark.parametrize("case", ["tie", "unstable", "unmapped", "one_run"])
def test_guarded_candidate_abstains_when_evidence_is_ambiguous(case):
    runs, topology = fixture()
    if case == "tie":
        for run in runs:
            run[0]["probabilities"] = probabilities(siblings=0.5, half_sibling=0.5)
    elif case == "unstable":
        runs[0][0]["probabilities"] = probabilities(half_sibling=1.0)
    elif case == "unmapped":
        topology["relationships"][0]["relationship_tuple"] = [0, 2, 1]
    else:
        runs = runs[:1]
    rows, _ = build_guarded_candidate(runs, topology)
    assert rows[0]["decision"] == "abstain"
    assert rows[0]["relation"] is None
