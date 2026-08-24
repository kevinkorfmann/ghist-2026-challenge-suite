from scripts.cluster_human_family_branches import find, union


def test_human_branch_union_find():
    parent = {sample: sample for sample in ("a", "b", "c")}
    union(parent, "a", "b")
    assert find(parent, "a") == find(parent, "b")
    assert find(parent, "a") != find(parent, "c")
