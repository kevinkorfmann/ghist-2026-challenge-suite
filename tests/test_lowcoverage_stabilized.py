import json

import pytest

from ghist_lowcoverage.stabilized import run_stabilized_comparison


@pytest.mark.parametrize('last_epoch, expected_stop', [(20, True), (200, False)])
def test_full_window_policy_and_cap_are_applied_to_all_runs(tmp_path, monkeypatch, last_epoch, expected_stop):
    import ghist_lowcoverage.stabilized as module

    tool = tmp_path / 'winsfs'
    tool.write_text('pinned-test-tool')
    indices = []
    for i in range(2):
        index = tmp_path / f'population_{i}.saf.idx'
        for suffix in ['.idx', '.gz', '.pos.gz']:
            (tmp_path / f'population_{i}.saf{suffix}').write_text('synthetic SAF fixture')
        indices.append(index)
    calls = []

    def run(command, *, stdout, stderr=None, check):
        calls.append(command)
        stdout.write('#SHAPE=<21/21>\n' + ' '.join(map(str, [999000, 1000] + [0] * 439)) + '\n')
        if stderr is not None:
            stderr.write(f'Finished epoch {last_epoch}\n')

    monkeypatch.setattr(module.subprocess, 'run', run)
    monkeypatch.setattr(module, 'fit_spectrum', lambda *args: {
        'targets': {'generations_since_split': 1000., 'pop1_population_size': 4000., 'pop2_population_size': 4000.},
        'optimizer_agreement': {'passed': True}})
    output = tmp_path / 'output'
    result = run_stabilized_comparison({'sample_sizes': [20, 20], 'task': 'split_isolation', 'length': 1_000_000},
        indices, output, tool)
    inference = [command for command in calls if '--seed' in command]
    assert len(inference) == 3
    for command in inference:
        assert command[command.index('--blocks') + 1] == '500'
        assert command[command.index('--window-size') + 1] == '500'
        assert command[command.index('--tolerance') + 1] == '0.000001'
        assert command[command.index('--max-epochs') + 1] == '200'
    assert [command[command.index('--seed') + 1] for command in inference] == ['761001', '761002', '761001']
    assert '--initial' not in inference[0] and '--initial' not in inference[1]
    assert inference[2][inference[2].index('--initial') + 1] == str(output / 'seed1.unfolded.sfs')
    assert set(result['early_stop_observed'].values()) == {expected_stop}
    assert set(result['callable_fraction'].values()) == {1.0}
    protocol = json.loads((output / 'protocol.json').read_text())
    assert protocol['window_size'] == protocol['blocks'] == 500 and protocol['tolerance'] == 1e-6
