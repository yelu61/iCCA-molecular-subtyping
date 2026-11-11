"""A renderer writing its session input must not change the frozen input bundle."""
from pathlib import Path
from types import SimpleNamespace
from unittest import TestCase, main
from unittest.mock import patch
import hashlib
import importlib.util
import json
import tempfile

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('icca_entry', REPO / 'run.py')
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)


class FrozenInputIsolationTest(TestCase):
    def test_renderer_write_is_confined_to_session(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            code = root / 'code'
            data = root / 'data'
            output = root / 'output'
            original = data / 'source_data/project/table.csv'
            original.parent.mkdir(parents=True)
            original.write_bytes(b'value\n1\n')
            (code / 'config').mkdir(parents=True)
            (code / '_implementation/project').mkdir(parents=True)
            (code / 'config/input_manifest.json').write_text(json.dumps([
                {'path': 'source_data/project/table.csv',
                 'sha256': hashlib.sha256(original.read_bytes()).hexdigest()}
            ]))
            (code / 'config/figure_jobs.json').write_text('{"jobs": []}')

            def renderer(command, **kwargs):
                project = Path(command[command.index('--original-root') + 1])
                (project / 'table.csv').write_bytes(b'value\n2\n')
                return SimpleNamespace(returncode=0)

            arguments = ['run.py', 'figures', '--data-root', str(data),
                         '--output', str(output), '--r-home', str(root / 'R')]
            with patch.object(entry, 'REPO', code), patch('sys.argv', arguments), \
                    patch.object(entry.subprocess, 'run', side_effect=renderer):
                with self.assertRaises(SystemExit) as result:
                    entry.main()
            self.assertEqual(result.exception.code, 0)
            self.assertEqual(original.read_bytes(), b'value\n1\n')
            self.assertEqual((output / 'session/source_data/project/table.csv').read_bytes(),
                             b'value\n2\n')


if __name__ == '__main__':
    main()
