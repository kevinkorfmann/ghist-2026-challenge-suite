import importlib.util
from pathlib import Path
import tempfile
import unittest

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


if __name__ == '__main__':
    unittest.main()
