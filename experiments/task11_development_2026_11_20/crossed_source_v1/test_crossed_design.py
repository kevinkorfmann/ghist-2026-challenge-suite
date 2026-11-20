"""Prospective nuisance-independence and unchanged-mechanics checks."""
import ast
import copy
import importlib.util
import inspect
import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest
from scripts import develop_task11_pedigree as dev

ROOT = Path(__file__).resolve().parents[1]
PLAN = json.loads((ROOT / "crossed_development_plan_v1.json").read_text())
SPEC = importlib.util.spec_from_file_location("crossed_generator", Path(__file__).with_name("generate_crossed_genotypes.py"))
generator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(generator)


def test_full_conditional_crossing_in_every_split_and_label():
    for split, replicates in (("train", 2), ("calibration", 1), ("development", 1)):
        cells = defaultdict(Counter)
        configs = [c for c in PLAN["cohorts"] if c["split"] == split]
        assert len(configs) == 54 * replicates
        for config in configs:
            generator.validate_assignment(config)
            for label, kind in zip(dev.LABELS, config["ancestry_kinds"], strict=True):
                cells[(label, config["population_size"], config["missing_rate"], config["complex_topology"])][kind] += 1
        assert len(cells) == 8 * 3 * 3 * 2
        assert all(dict(counts) == {kind: replicates for kind in ("unrelated", "first_cousins", "siblings")} for counts in cells.values())


def test_balanced_oof_folds_and_no_cohort_leakage():
    train = [c for c in PLAN["cohorts"] if c["split"] == "train"]
    for fold in range(4):
        held = [c for c in train if c["oof_fold"] == fold]
        fitted = [c for c in train if c["oof_fold"] != fold]
        assert len(held) == 27 and len(fitted) == 81
        assert not {c["seed"] for c in held} & {c["seed"] for c in fitted}
        assert len({(c["population_size"], c["missing_rate"], c["ancestry_shift"]) for c in held}) == 27
        for key in ("population_size", "missing_rate", "ancestry_shift"):
            assert sorted(Counter(c[key] for c in held).values()) == [9, 9, 9]
        assert sorted(Counter(c["complex_topology"] for c in held).values()) == [13, 14]


def test_training_third_branches_balanced_within_nuisance_cells():
    cells = defaultdict(list)
    for config in PLAN["cohorts"]:
        if config["split"] == "train":
            cells[(config["population_size"], config["missing_rate"], config["ancestry_shift"], config["complex_topology"])].append(config)
    for cell, replicas in cells.items():
        assert len(replicas) == 2
        for index, label in enumerate(dev.LABELS):
            count = sum(config["third_branches"][index] for config in replicas)
            assert count == (2 if cell[-1] and label.startswith("avuncular") else 1)


def test_all_new_streams_and_future_stress_confirmation_are_disjoint():
    all_configs = PLAN["cohorts"] + [c for values in PLAN["reserved_future_cohorts"].values() for c in values]
    seeds = [c["seed"] for c in all_configs]
    assert len(seeds) == len(set(seeds)) == 232
    assert not set(seeds) & set(PLAN["old_seed_union"])
    assert PLAN["gates"] == dev.GATES
    for split in ("stress", "confirmation"):
        configs = PLAN["reserved_future_cohorts"][split]
        assert len(configs) == 8 and sum(c["complex_topology"] for c in configs) == 4
        for label_index in range(8):
            assert sorted(Counter(c["ancestry_kinds"][label_index] for c in configs).values()) == [2, 3, 3]


def test_assignment_is_explicit_and_does_not_depend_on_reindexing():
    config = copy.deepcopy(PLAN["cohorts"][0])
    config["index"] = 987654
    generator.validate_assignment(config)
    body = inspect.getsource(generator.generate_genotypes)
    assert 'config["index"]' not in body
    assert 'config["ancestry_kinds"][index]' in body
    assert 'config["third_branches"][index]' in body


@pytest.mark.parametrize("field,value", [("ancestry_kinds", []), ("third_branches", [1] * 8), ("max_branch_size", 5)])
def test_invalid_explicit_assignments_fail(field, value):
    config = copy.deepcopy(PLAN["cohorts"][0])
    config[field] = value
    with pytest.raises(ValueError):
        generator.validate_assignment(config)


def test_genotype_mechanics_ast_diff_is_only_explicit_assignment_and_validation():
    expected = inspect.getsource(dev.generate_genotypes)
    expected = expected.replace('    import msprime\n', '    validate_assignment(config)\n    import msprime\n')
    expected = expected.replace('    kinds = ("unrelated", "first_cousins", "siblings")\n', '')
    expected = expected.replace('third = bool((config["index"] + index) % 2) or (complex_topology and label.startswith("avuncular"))', 'third = config["third_branches"][index]')
    expected = expected.replace('kinds[(index + config["index"]) % 3]', 'config["ancestry_kinds"][index]')
    actual = inspect.getsource(generator.generate_genotypes)
    assert ast.dump(ast.parse(actual)) == ast.dump(ast.parse(expected))
