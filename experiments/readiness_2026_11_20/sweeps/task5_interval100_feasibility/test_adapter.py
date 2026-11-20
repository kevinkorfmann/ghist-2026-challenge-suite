import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('adapter', Path(__file__).with_name('adapter.py'))
ADAPTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ADAPTER)


class TransformationTests(unittest.TestCase):
    def test_only_bookkeeping_path_and_passive_logging_change(self):
        original = 'initialize() {\n' + ADAPTER.DECLARATION + '\n    defineConstant("trees_file", "/tmp/checkpoint");\n}\nfunction(void)save(void){sim.treeSeqOutput(trees_file);}\nfunction(void)restore(void){sim.readFromPopulationFile(trees_file);}\n'
        changed = ADAPTER.transform_script(original, '/retained/raw.trees', '/retained/phase.tsv')
        core, diagnostics = changed.split('// Passive phase logging:', 1)
        restored = core.replace(ADAPTER.FIXED_DECLARATION, ADAPTER.DECLARATION).replace('/retained/raw.trees', '/tmp/checkpoint')
        self.assertEqual(restored.rstrip(), original.rstrip())
        self.assertIn('community.tick % 100', diagnostics)
        self.assertIn('sim.getValue("n_restores")', diagnostics)
        self.assertNotIn('runif', diagnostics)
        self.assertNotIn('sample(', diagnostics)
        self.assertNotIn('setSubpopulationSize', diagnostics)

    def test_fixed_interval_is_exactly_100(self):
        self.assertEqual(ADAPTER.FIXED_DECLARATION,
                         'initializeTreeSeq(simplificationInterval=100, timeUnit="generations");')

    def test_adapter_diff_is_only_interval_and_diagnostic_names(self):
        here = Path(__file__).parent
        previous = (here / 'reference' / 'adapter.py').read_text()
        expected = previous.replace('simplificationInterval=20,', 'simplificationInterval=100,')
        expected = expected.replace('generated_interval20.slim', 'generated_interval100.slim')
        expected = expected.replace("'interval20_sha256'", "'interval100_sha256'")
        self.assertEqual((here / 'adapter.py').read_text(), expected)

    def test_missing_or_repeated_declaration_fails(self):
        for text in ['', ADAPTER.DECLARATION, ADAPTER.DECLARATION * 2]:
            with self.assertRaises(ValueError):
                ADAPTER.transform_script(text, '/a', '/b')

    def test_unsafe_eidos_path_fails(self):
        original = ADAPTER.DECLARATION + '\n    defineConstant("trees_file", "/tmp/checkpoint");\n'
        for path in ['/a"b', '/a\\b', '/a\nb']:
            with self.assertRaises(ValueError):
                ADAPTER.transform_script(original, path, '/b')

    def test_original_loader_receives_identical_final_raw_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            retained = Path(directory) / 'retained.trees'
            original = Path(directory) / 'temporary.trees'
            retained.write_bytes(b'complete raw tree payload')
            ADAPTER.restore_loader_input(retained, original)
            self.assertEqual(retained.read_bytes(), original.read_bytes())
            with self.assertRaises(ValueError):
                ADAPTER.restore_loader_input(retained, retained)

    def test_missing_final_tree_cannot_silently_pass_loader(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                ADAPTER.restore_loader_input(Path(directory) / 'missing', Path(directory) / 'output')

    def test_final_verification_error_is_retained_and_raised(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'status.json'
            with self.assertRaisesRegex(RuntimeError, 'verification failed'):
                ADAPTER.finalize_status({'state': 'COMPLETED'}, path, ValueError('source changed'), [])
            status = json.loads(path.read_text())
            self.assertEqual(status['state'], 'FAILED')
            self.assertFalse(status['frozen_inputs_unchanged'])
            self.assertEqual(status['provenance_error'], 'source changed')

    def test_monitor_exception_is_retained_and_prevents_complete_status(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'resources.jsonl'
            errors = []
            with patch.object(ADAPTER, 'resource_row', side_effect=PermissionError('cgroup unreadable')):
                ADAPTER.monitor_resources(path, [], {'value': 'test'}, None, errors)
            self.assertIn('cgroup unreadable', errors[0])
            self.assertEqual(json.loads(path.with_name('monitor_failure.json').read_text())['errors'], errors)
            status_path = Path(directory) / 'status.json'
            with self.assertRaisesRegex(RuntimeError, 'monitor failed'):
                ADAPTER.finalize_status({'state': 'COMPLETED'}, status_path, None, errors)
            status = json.loads(status_path.read_text())
            self.assertEqual(status['state'], 'FAILED')
            self.assertFalse(status['resource_monitor_complete'])

    def test_missing_memory_counters_fail_before_simulation(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'No readable memory'):
                ADAPTER.resource_row([Path(directory)], {'value': 'test'})


if __name__ == '__main__':
    unittest.main()
