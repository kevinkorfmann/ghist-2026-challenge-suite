"""Read-only review of actual freeze entrypoints; synthetic arithmetic only."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile

root = Path(sys.argv[1]).resolve()
sys.path[:0] = [str(root / 'src'), str(root), str(root / 'tests')]
import pytest
from test_demography_freeze_integrity import fields, freeze, sources_for

results = []
for task in [3, 4]:
    scenarios = ['consistent', 'bootstrap_time', 'observed_time']
    if task == 4:
        scenarios += ['consistent_nested_refined', 'refined_flat', 'refined_nested']
    for scenario in scenarios:
        folder = Path(tempfile.mkdtemp(prefix='ghist-v3-independent-'))
        module, sources, flags = sources_for(task, folder)
        if scenario in ['consistent_nested_refined', 'refined_nested']:
            sources['refined'] = {'fit': sources['refined']}
        with pytest.MonkeyPatch.context() as monkeypatch, contextlib.redirect_stdout(io.StringIO()):
            valid = freeze(module, sources, flags, folder, monkeypatch)
        destination = folder / 'freeze.json'
        before = destination.read_bytes()
        parameter = 'generations' if task == 3 else 'generations_since_split'
        if scenario.startswith('consistent'):
            results.append({'task': task, 'scenario': scenario, 'accepted': True,
                            'timing': valid['testing_prediction'][parameter]})
            continue
        if scenario == 'bootstrap_time':
            sources['bootstrap']['replica_records'][0]['fit'][parameter] *= 10
        elif scenario == 'observed_time':
            fits, median, _, _ = fields(task, sources)
            for projection in (list(fits) if task == 3 else ['70x50', '94x68']):
                fits[projection][parameter] *= 10
            median[parameter] *= 10
        else:
            refined = sources['refined'].get('fit', sources['refined'])
            refined[parameter] *= 1.0001
        try:
            with pytest.MonkeyPatch.context() as monkeypatch, contextlib.redirect_stdout(io.StringIO()):
                freeze(module, sources, flags, folder, monkeypatch)
        except ValueError as error:
            rejected = str(error)
        else:
            raise AssertionError(f'Corruption accepted: task={task}, scenario={scenario}')
        assert destination.read_bytes() == before
        assert not [path for path in folder.iterdir() if '.freeze.json.' in path.name]
        results.append({'task': task, 'scenario': scenario, 'rejected': rejected,
                        'existing_output_preserved_sha256': hashlib.sha256(before).hexdigest(),
                        'partial_output_files': 0})
print(json.dumps(results, indent=2))
