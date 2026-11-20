"""Real freeze regressions for malformed secondary-contact sampling intervals."""
import importlib.util
import json
from pathlib import Path

import pytest

from scripts.freeze_secondary_contact import INTERVAL_KEYS

_spec = importlib.util.spec_from_file_location(
    "interval_fixture", Path(__file__).with_name("test_demography_freeze_cross_review.py")
)
_fixture = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fixture)


@pytest.mark.parametrize("parameter", INTERVAL_KEYS)
@pytest.mark.parametrize("bounds", [
    [2., 1.], [-3., -1.], [float("nan"), 3.], [1., float("inf")],
    [True, 3.], [1.], [1., 2., 3.], "1,2", [None, 3.],
])
def test_malformed_interval_cannot_replace_existing_freeze(tmp_path, monkeypatch, parameter, bounds):
    sources = _fixture.secondary_contact_sources()
    sources["bootstrap"]["observed_sampling_interval_95"][parameter] = bounds
    output = tmp_path / "freeze.json"
    original = b'{"preserved": true}\n'
    output.write_bytes(original)
    with pytest.raises(ValueError):
        _fixture.run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert output.read_bytes() == original
    assert not list(tmp_path.glob("freeze.json.*"))


@pytest.mark.parametrize("change", ["empty", "missing", "extra"])
def test_interval_requires_exact_producer_parameter_set(tmp_path, monkeypatch, change):
    sources = _fixture.secondary_contact_sources()
    intervals = sources["bootstrap"]["observed_sampling_interval_95"]
    if change == "empty":
        intervals.clear()
    elif change == "missing":
        intervals.pop(INTERVAL_KEYS[0])
    else:
        intervals["unexpected"] = [1., 2.]
    with pytest.raises(ValueError, match="Interval parameters"):
        _fixture.run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert not (tmp_path / "freeze.json").exists()


@pytest.mark.parametrize("parameter", INTERVAL_KEYS)
def test_interval_zero_boundary_matches_parameter_domain(tmp_path, monkeypatch, parameter):
    sources = _fixture.secondary_contact_sources()
    sources["bootstrap"]["observed_sampling_interval_95"][parameter] = [0., 0.]
    if parameter in ("generations_since_migration", "migration_rate"):
        result = _fixture.run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
        assert result["conservative_sampling_interval_95"][parameter] == [0., 0.]
    else:
        with pytest.raises(ValueError, match="physical domain"):
            _fixture.run_secondary_contact_freeze(sources, tmp_path, monkeypatch)


def test_valid_interval_schema_and_values_are_preserved(tmp_path, monkeypatch):
    sources = _fixture.secondary_contact_sources()
    expected = sources["bootstrap"]["observed_sampling_interval_95"]
    result = _fixture.run_secondary_contact_freeze(sources, tmp_path, monkeypatch)
    assert result["conservative_sampling_interval_95"] == expected
    assert set(result["conservative_sampling_interval_95"]) == set(INTERVAL_KEYS)
    json.loads((tmp_path / "freeze.json").read_text(), parse_constant=lambda _: pytest.fail("nonfinite JSON"))
