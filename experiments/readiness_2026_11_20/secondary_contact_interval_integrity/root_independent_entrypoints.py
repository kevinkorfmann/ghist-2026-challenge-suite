"""Bounded synthetic producer/freezer contract checks; no model evaluations."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile

root = Path(sys.argv[1]).resolve()
sys.path[:0] = [str(root / 'src'), str(root)]

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

fixed = module('fixed_task2', root / 'scripts/freeze_secondary_contact.py')
baseline = module('baseline_task2', root / 'baseline_freeze_secondary_contact.py')
fixture = module('fixture_task2', root / 'tests/test_demography_freeze_cross_review.py')
producer = module('producer_task2', root / 'scripts/summarize_secondary_contact_generator_bootstrap.py')
from ghist_secondary_contact.contract import SPECTRUM_CONTRACT

assert fixed.INTERVAL_KEYS == producer.PARAMETER_NAMES
work = Path(tempfile.mkdtemp(prefix='ghist-task2-synthetic-inputs-'))

def invoke(mod, argv):
    previous = sys.argv
    sys.argv = argv
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mod.main()
    finally:
        sys.argv = previous

checks = []

def compare_valid(name, sources):
    case = work / name
    case.mkdir()
    output = case / 'freeze.json'
    argv = ['freeze', '--out', str(output)]
    for flag, source in sources.items():
        path = case / (flag + '.json')
        path.write_text(json.dumps(source))
        argv += ['--' + flag, str(path)]
    invoke(baseline, argv)
    before = output.read_bytes()
    invoke(fixed, argv)
    after = output.read_bytes()
    assert before == after
    assert not list(case.glob('freeze.json.*'))
    checks.append({'case': name, 'baseline_and_patch_byte_identical': True,
                   'output_sha256': hashlib.sha256(after).hexdigest(),
                   'intervals': json.loads(after)['conservative_sampling_interval_95']})

compare_valid('ordinary_complete', fixture.secondary_contact_sources())
for contact in (0.0, 2.0):
    sources = fixture.secondary_contact_sources()
    for fit in (sources['raw-fit'], sources['standard-fit'], sources['calibration']['observed_accepted']):
        fit.update(generations_since_contact=contact, migration_rate=0.0)
    sources['bootstrap']['observed_sampling_interval_95']['generations_since_migration'] = [0.0, 0.0]
    sources['bootstrap']['observed_sampling_interval_95']['migration_rate'] = [0.0, 0.0]
    compare_valid('zero_migration_contact_' + str(contact), sources)

calibration = work / 'synthetic_calibration'
calibration.mkdir()
truth = [10000., 20000., 3000., 1500., 0.0001]
for index, ratio in enumerate((0.9, 1.0, 1.1)):
    values = [x * ratio for x in truth]
    record = dict(zip(fixed.PARAMETER_KEYS, values))
    record.update(truth=truth, spectrum_contract=SPECTRUM_CONTRACT)
    (calibration / f'fit{index}.json').write_text(json.dumps(record))
produced = work / 'produced.json'
invoke(producer, ['producer', '--inputs', str(calibration / 'fit*.json'), '--out', str(produced), '--observed', *map(str, truth)])
sources = fixture.secondary_contact_sources()
sources['bootstrap'] = json.loads(produced.read_text())
assert set(sources['bootstrap']['observed_sampling_interval_95']) == set(fixed.INTERVAL_KEYS)
compare_valid('actual_producer_output', sources)

for name, mutate in (
    ('empty', lambda s: s['bootstrap'].__setitem__('observed_sampling_interval_95', {})),
    ('reversed', lambda s: s['bootstrap']['observed_sampling_interval_95'].__setitem__('migration_rate', [2., 1.])),
    ('negative_population', lambda s: s['bootstrap']['observed_sampling_interval_95'].__setitem__('mainland_population_size', [-3., -1.])),
    ('late_serialization_nan', lambda s: s['calibration']['untouched_evaluation'].__setitem__('synthetic_nonfinite_payload', float('nan'))),
):
    case = work / name
    case.mkdir()
    sources = fixture.secondary_contact_sources()
    mutate(sources)
    output = case / 'freeze.json'
    original = b'{"preserved": true}\n'
    output.write_bytes(original)
    argv = ['freeze', '--out', str(output)]
    for flag, source in sources.items():
        path = case / (flag + '.json')
        path.write_text(json.dumps(source))
        argv += ['--' + flag, str(path)]
    try:
        invoke(fixed, argv)
    except ValueError as error:
        detail = str(error)
    else:
        raise AssertionError(f'{name} accepted')
    assert output.read_bytes() == original
    assert not list(case.glob('freeze.json.*'))
    checks.append({'case': name, 'rejected': True, 'existing_output_unchanged': True,
                   'no_temporary_files_correct_prefix': True, 'error': detail})

print(json.dumps({'scope': 'synthetic arithmetic and actual CLI checks only; no fitting or model evaluation',
                  'snapshot_root': str(root), 'synthetic_inputs_root': str(work),
                  'checks': checks}, indent=2, allow_nan=False))
