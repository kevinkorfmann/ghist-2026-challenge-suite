"""Exercise the comparison metric without executing the compute-stage module."""
import ast
from pathlib import Path

import pytest

source = ast.parse(Path(__file__).with_name("evaluate_pragmatic_on_crossed_v2.py").read_text())
function = next(node for node in source.body if isinstance(node, ast.FunctionDef) and node.name == "weighted_jaccard")
namespace = {}
exec(compile(ast.Module(body=[function], type_ignores=[]), "comparison_metric", "exec"), namespace)  # noqa: S102 - trusted local function AST only
weighted_jaccard = namespace["weighted_jaccard"]


def row(a, b, label):
    return {"indID1": a, "indID2": b, "relation": label}


def test_extra_pairs_are_false_positives_and_missing_pairs_are_false_negatives():
    truth = [row("a", "b", "x"), row("c", "d", "x"), row("e", "f", "y")]
    # One correct x, one missing x, one y misclassified x, and one extra x.
    predictions = [row("b", "a", "x"), row("f", "e", "x"), row("g", "h", "x")]
    # x IoU=1/(1TP+1FN+2FP)=1/4, weighted by2/3; y IoU=0.
    assert weighted_jaccard([(None, {"pairs": truth})], [predictions]) == pytest.approx(1 / 6)


def test_cohort_namespaces_are_distinct_and_perfect_is_one():
    first = [row("a", "b", "x")]
    second = [row("a", "b", "y")]
    assert weighted_jaccard([(None, {"pairs": first}), (None, {"pairs": second})], [first, second]) == 1


def test_duplicate_prediction_is_rejected():
    truth = [row("a", "b", "x")]
    with pytest.raises(ValueError):
        weighted_jaccard([(None, {"pairs": truth})], [truth + truth])


def test_empty_truth_is_rejected():
    with pytest.raises(ValueError):
        weighted_jaccard([(None, {"pairs": []})], [[]])
