"""Queue-delay guards must refuse every stage after its registered deadline."""
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path

import pytest

from scripts import lowcoverage_clock_revision_guard as guard


@pytest.mark.parametrize('stage', ['development', 'prepare', 'confirmation', 'final'])
def test_stage_start_deadline_fails_closed(tmp_path, monkeypatch, stage):
    source = Path(__file__).resolve().parents[1] / 'experiments/lowcoverage_2026_11_20/task14_extension_protocol_v2.json'
    (tmp_path / source.name).write_bytes(source.read_bytes())
    plan = json.loads(source.read_text())
    deadline = datetime.fromisoformat(plan['stage_latest_start_utc'][stage].replace('Z', '+00:00'))

    class ControlledClock(datetime):
        current = deadline - timedelta(seconds=1)

        @classmethod
        def now(cls, tz=None):
            assert tz == timezone.utc
            return cls.current

    monkeypatch.setattr(guard, 'datetime', ControlledClock)
    assert guard.check_clock(tmp_path, stage)['actual_start_utc'] == ControlledClock.current.isoformat()
    ControlledClock.current = deadline + timedelta(seconds=1)
    with pytest.raises(RuntimeError, match='actual start exceeds registered deadline'):
        guard.check_clock(tmp_path, stage)
