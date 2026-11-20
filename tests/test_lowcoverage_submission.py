from pathlib import Path

import pytest
import yaml

from ghist_lowcoverage.submission import format_submission_draft

TEMPLATES = Path(__file__).resolve().parents[1] / 'experiments/lowcoverage_2026_11_20/official_templates'


@pytest.mark.parametrize('task,name,targets', [
    ('bottleneck', 'bottleneck.submission.yaml', {'generations': 2000, 'post_decline_fraction': .2}),
    ('split_isolation', 'split_isolation_low_coverage.submission.yaml', {
        'pop1_population_size': 5000, 'pop2_population_size': 15000, 'generations_since_split': 5000}),
])
def test_official_template_fill_retains_units(task, name, targets):
    template = TEMPLATES / name
    original = yaml.safe_load(template.read_text())
    result = yaml.safe_load(format_submission_draft(task, targets, template))
    assert result['parameters'] == targets
    assert {k: v for k, v in result.items() if k != 'parameters'} == {k: v for k, v in original.items() if k != 'parameters'}


@pytest.mark.parametrize('value', [0, -1, float('nan'), float('inf'), True, '2000'])
def test_invalid_parameter_rejected(value):
    with pytest.raises(ValueError, match='finite positive'):
        format_submission_draft('bottleneck', {'generations': value, 'post_decline_fraction': .2}, TEMPLATES / 'bottleneck.submission.yaml')


def test_template_tamper_rejected(tmp_path):
    path = tmp_path / 'template.yaml'
    path.write_text((TEMPLATES / 'bottleneck.submission.yaml').read_text() + '\n')
    with pytest.raises(ValueError, match='checksum'):
        format_submission_draft('bottleneck', {'generations': 2000, 'post_decline_fraction': .2}, path)


def test_schema_rejected():
    with pytest.raises(ValueError, match='schema'):
        format_submission_draft('bottleneck', {'generations': 2000, 'fraction': .2}, TEMPLATES / 'bottleneck.submission.yaml')
