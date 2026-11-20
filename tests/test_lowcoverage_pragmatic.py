"""Pragmatic release still requires explicit bound root approval before input access."""
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import sys

import pytest

from scripts import lowcoverage_pragmatic_guard as guard
from scripts import lowcoverage_predict_pragmatic as predict


@pytest.mark.parametrize('stage', ['development', 'final'])
def test_actual_start_cutoff(tmp_path, monkeypatch, stage):
    deadline = datetime(2026, 11, 20, 20, 45, tzinfo=timezone.utc)
    (tmp_path / 'task14_pragmatic_protocol_v1.json').write_text(json.dumps({
        'stage_latest_start_utc': {stage: deadline.isoformat()}}))

    class Clock(datetime):
        current = deadline - timedelta(seconds=1)

        @classmethod
        def now(cls, tz=None):
            assert tz == timezone.utc
            return cls.current

    monkeypatch.setattr(guard, 'datetime', Clock)
    assert guard.check_clock(tmp_path, stage)['actual_start_utc'] == Clock.current.isoformat()
    Clock.current = deadline + timedelta(seconds=1)
    with pytest.raises(RuntimeError, match='actual start exceeds'):
        guard.check_clock(tmp_path, stage)


@pytest.mark.parametrize('complete, authorized', [(False, True), (True, False), (True, True)])
def test_explicit_decision_and_both_complete_precede_input_access(tmp_path, monkeypatch, complete, authorized):
    report_path = tmp_path / 'task14_pragmatic_development_gate_audit_v1.json'
    report_path.write_text(json.dumps({'freeze_sha256': 'fixture', 'mode': 'development',
        'all_cohorts_completed': complete, 'primary_target_count': 6}))
    decision_path = tmp_path / 'decision.json'
    decision_path.write_text(json.dumps({'authorize_final_inference': authorized,
        'strict_qualification': False, 'freeze_sha256': 'fixture',
        'wrapper_sha256': predict.sha256(Path(predict.__file__)),
        'development_report_sha256': predict.sha256(report_path), 'primary_estimate': 'seed1'}))
    monkeypatch.setattr(sys, 'argv', ['predict', '--frozen-root', str(tmp_path),
        '--freeze-sha', 'fixture', '--root-decision', str(decision_path),
        '--decision-sha', predict.sha256(decision_path), '--directory', str(tmp_path / 'bam'),
        '--reference', str(tmp_path / 'ref'), '--output', str(tmp_path / 'out')])
    monkeypatch.setattr(predict, 'verify', lambda *args: None)
    monkeypatch.setattr(predict, 'check_clock', lambda *args: {})
    accesses = []

    def bind(*args):
        accesses.append(True)
        raise RuntimeError('input-access sentinel')

    monkeypatch.setattr(predict, 'bind_inputs', bind)
    if complete and authorized:
        with pytest.raises(RuntimeError, match='input-access sentinel'):
            predict.main()
        assert accesses == [True]
    else:
        with pytest.raises(ValueError, match='explicit root pragmatic decision'):
            predict.main()
        assert not accesses
