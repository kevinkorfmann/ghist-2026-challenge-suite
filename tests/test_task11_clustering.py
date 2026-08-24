import numpy as np
import pandas as pd

from scripts.cluster_task11_families import affinity_matrix, canonical_labels


def test_affinity_is_symmetric() -> None:
    pairs = pd.DataFrame(
        [{"indID1": "id_0", "indID2": "id_1", "relatedness_ge_10": 0.25}]
    )
    result = affinity_matrix(pairs, ["id_0", "id_1"], "relatedness_ge_10")
    assert np.array_equal(result, np.asarray([[0.25, 0.25], [0.25, 0.25]]))


def test_canonical_labels_follow_smallest_identifier() -> None:
    labels = canonical_labels(np.asarray([4, 2, 4, 2]), ["id_0", "id_1", "id_2", "id_3"])
    assert labels.tolist() == [0, 1, 0, 1]
