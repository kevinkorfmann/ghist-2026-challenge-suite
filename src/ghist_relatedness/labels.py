from __future__ import annotations


HUMAN_LABELS = {
    (1, 1, 2): "siblings",
    (1, 2, 2): "avuncular",
    (1, 3, 2): "avuncular_once_removed",
    (1, 4, 2): "avuncular_twice_removed",
    (2, 2, 2): "full_1st_cousin",
    (2, 3, 2): "full_1st_cousin_once_removed",
    (2, 4, 2): "full_1st_cousin_twice_removed",
    (2, 5, 2): "full_1st_cousin_thrice_removed",
    (3, 3, 2): "full_2nd_cousin",
    (3, 4, 2): "full_2nd_cousin_once_removed",
    (3, 5, 2): "full_2nd_cousin_twice_removed",
    (4, 4, 2): "full_3rd_cousin",
    (4, 5, 2): "full_3rd_cousin_once_removed",
    (5, 5, 2): "full_4th_cousin",
}

DOG_LABELS = {
    (1, 1, 2): "siblings",
    (1, 1, 1): "half_sibling",
    (1, 2, 2): "avuncular",
    (1, 3, 2): "avuncular_once_removed",
    (1, 4, 2): "avuncular_twice_removed",
    (2, 2, 2): "full_1st_cousin",
    (2, 2, 1): "half_1st_cousin",
    (1, 2, 1): "half_avuncular",
}


def canonical_tuple(relationship: tuple[int, int, int]) -> tuple[int, int, int]:
    up, down, ancestors = relationship
    return min(up, down), max(up, down), ancestors


def accepted_label(
    relationship: tuple[int, int, int], species: str
) -> str | None:
    labels = HUMAN_LABELS if species == "human" else DOG_LABELS
    return labels.get(canonical_tuple(relationship))


def meioses_count(relationship: tuple[int, int, int]) -> int:
    return relationship[0] + relationship[1]
