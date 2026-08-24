import ast
from pathlib import Path


def test_task9_writer_uses_official_headers() -> None:
    tree = ast.parse(Path("scripts/fit_task9_ibd_classifier.py").read_text())
    header_tuples = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Tuple)
        and all(isinstance(item, ast.Constant) for item in node.elts)
        and node.elts
        and node.elts[0].value == "indID1"
    ]
    values = [tuple(item.value for item in node.elts) for node in header_tuples]
    assert ("indID1", "indID2", "meioses_count", "relation") in values
