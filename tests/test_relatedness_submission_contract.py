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


def test_task9_freezer_rejects_retired_labels_without_overwriting_candidate(tmp_path, monkeypatch):
    import pytest
    from scripts import freeze_task9_ibd
    raw = tmp_path / 'raw.csv'
    candidate = tmp_path / 'candidate.csv'
    raw.write_text('indID1,indID2,relation\na,b,half_sibling\n')
    candidate.write_text('previous validated artifact\n')
    monkeypatch.setattr(freeze_task9_ibd, 'RAW_PREDICTIONS', raw)
    monkeypatch.setattr(freeze_task9_ibd, 'CANDIDATE', candidate)
    with pytest.raises(ValueError, match='Unaccepted current Task 9'):
        freeze_task9_ibd.write_candidate()
    assert candidate.read_text() == 'previous validated artifact\n'
