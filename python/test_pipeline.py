"""Two-run shell integration: real 05/06/generator; data preparation/R are stubs.

No Space-Track credentials or production outputs are touched. This tests ordering,
backup/replacement and downstream hand-offs, not live data access or R statistics.
Run from root: python -B python/test_pipeline.py.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


class PipelineTests(unittest.TestCase):
    def test_two_complete_orchestrated_runs(self):
        bash = shutil.which('bash')
        self.assertIsNotNone(bash, 'Bash is required for the shell pipeline test')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            shutil.copytree(ROOT / 'python', root / 'python', ignore=shutil.ignore_patterns('__pycache__'))
            shutil.copy2(ROOT / 'run_all.sh', root / 'run_all.sh')
            processed = root / 'data/processed'
            processed.mkdir(parents=True)
            shutil.copy2(ROOT / 'data/processed/params.csv', processed / 'params.csv')
            conjunctions = pd.read_csv(ROOT / 'data/processed/conjunctions.csv')
            conjunctions[conjunctions.conj_id.isin([183, 911, 1198])].to_csv(processed / 'conjunctions.csv', index=False)
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in processed.iterdir()}
            for name in ['01_download.py', '02_clean.py', '03_residuals.py', '04_screen_conjunctions.py']:
                (root / 'python' / name).write_text('print("Stubbed unavailable upstream data stage")\n')
            bin_dir = root / 'bin'
            bin_dir.mkdir()
            python_runner = bin_dir / 'pipeline-python'
            python_runner.write_text(f'#!{sys.executable}\n' + '''import json, subprocess, sys
from pathlib import Path
with Path('calls.jsonl').open('a') as f: f.write(json.dumps(sys.argv[1:])+'\\n')
subprocess.run([sys.executable, '-B', *sys.argv[1:]], check=True)
''')
            r_runner = bin_dir / 'Rscript'
            r_runner.write_text(f'#!{sys.executable}\n' + '''import json, sys
from pathlib import Path
import pandas as pd
with Path('calls.jsonl').open('a') as f: f.write(json.dumps(sys.argv[1:])+'\\n')
if sys.argv[1] in ('r/12_poisson_fleet_model.R', 'r/14_plots.R'):
    data = pd.read_csv('data/processed/pc_results.csv')
    assert len(data) == 15
    assert set(data.method) == {'exact', 'mc', 'is'}
    assert Path('data/processed/mc_convergence.csv').exists()
''')
            for runner in [python_runner, r_runner]:
                runner.chmod(0o755)
            env = dict(os.environ, PYTHON=str(python_runner), PATH=str(bin_dir)+os.pathsep+os.environ['PATH'])
            expected = ['python/01_download.py', 'python/02_clean.py', 'python/03_residuals.py',
                        'r/10_error_model.R', 'python/04_screen_conjunctions.py', 'python/05_exact_pc.py',
                        'python/06_monte_carlo.py', 'python/generate_results.py',
                        'r/11_gaussian_vs_learned_tests.R', 'r/12_poisson_fleet_model.R',
                        'r/13_calibration.R', 'r/14_plots.R']
            old_bytes, old_tables = {}, {}
            for run in range(2):
                completed = subprocess.run([bash, 'run_all.sh'], cwd=root, env=env,
                                           capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                for name, count in [('pc_results.csv', 15), ('pc_mc.csv', 6), ('pc_is.csv', 6),
                                    ('mc_convergence.csv', 42), ('pc_radius_sensitivity.csv', 18)]:
                    path = processed / name
                    table = pd.read_csv(path)
                    self.assertEqual(len(table), count)
                    if run:
                        for column in table.columns:
                            if column != 'runtime_s':
                                np.testing.assert_array_equal(table[column], old_tables[name][column])
                        backups = list((processed / 'backups').glob('*/' + name))
                        self.assertTrue(any(p.read_bytes() == old_bytes[name] for p in backups))
                    old_bytes[name], old_tables[name] = path.read_bytes(), table
            calls = [json.loads(line)[0] for line in (root / 'calls.jsonl').read_text().splitlines()]
            self.assertEqual(calls, expected * 2)
            for name, digest in before.items():
                self.assertEqual(hashlib.sha256((processed / name).read_bytes()).hexdigest(), digest)


if __name__ == '__main__':
    unittest.main()
