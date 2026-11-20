"""Verify the fourteen archived submitted files using the recorded SHA-256 values."""
from pathlib import Path
import hashlib
import json


def main():
    root = Path(__file__).resolve().parents[1]
    folder = root / 'experiments/final_submission_2026_11_20'
    audit = json.loads((folder / 'all_tasks_scored_audit.json').read_text())
    if len(audit['tasks']) != 14:
        raise ValueError('Expected fourteen task records')
    for task in audit['tasks']:
        relative = task['candidate']['path'].split('/ghist/', 1)[-1]
        candidate = root / relative
        if not candidate.resolve().is_relative_to(root):
            raise ValueError('Candidate path leaves repository')
        digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if digest != task['candidate']['sha256']:
            raise ValueError(f"Task {task['task']}: candidate digest mismatch")
        receipt = json.loads((folder / f"task{task['task']}" / 'submission_receipt.json').read_text())
        if receipt['status'] != 'SCORED' or str(receipt['submission_id']) != str(task['submission_id']):
            raise ValueError(f"Task {task['task']}: inconsistent receipt")
    print('PASS: all 14 submitted files match their recorded hashes and scored receipts')


if __name__ == '__main__':
    main()
