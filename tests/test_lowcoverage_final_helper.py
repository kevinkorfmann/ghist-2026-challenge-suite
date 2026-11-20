"""Independent YAML checks use synthetic values and the verified official template."""
from pathlib import Path

import pytest
import yaml

from scripts.lowcoverage_audit_final_pragmatic import audit_draft


TEMPLATE = (Path(__file__).resolve().parents[1] / 'experiments/lowcoverage_2026_11_20'
            / 'official_templates/split_isolation_low_coverage.submission.yaml')
TARGETS = {'pop1_population_size': 4000.0, 'pop2_population_size': 6000.0,
           'generations_since_split': 1000.0}


def draft(tmp_path):
    document = yaml.safe_load(TEMPLATE.read_text())
    document['parameters'] = dict(TARGETS)
    path = tmp_path / 'draft.yaml'
    path.write_text(yaml.safe_dump(document))
    return path


def test_exact_official_draft(tmp_path):
    audit_draft(draft(tmp_path), TEMPLATE, TARGETS, 'split_isolation')


@pytest.mark.parametrize('value', [True, float('inf'), 0, -1, '4000', 4001])
def test_malformed_or_changed_estimate_rejected(tmp_path, value):
    path = draft(tmp_path)
    document = yaml.safe_load(path.read_text())
    document['parameters']['pop1_population_size'] = value
    path.write_text(yaml.safe_dump(document))
    with pytest.raises(ValueError):
        audit_draft(path, TEMPLATE, TARGETS, 'split_isolation')
