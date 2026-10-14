"""Integration checks: .venv/bin/python -B python/test_result_generation.py."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('generator', Path(__file__).with_name('generate_results.py'))
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)


class IntegrationTests(unittest.TestCase):
    def test_generated_log_diagnostics_and_801(self):
        logs = pd.read_csv(g.PROCESSED / 'is_log_diagnostics.csv')
        results = pd.read_csv(g.PROCESSED / 'pc_is.csv')
        self.assertEqual(tuple(logs.columns), g.LOG_COLUMNS)
        self.assertEqual(len(logs), len(results))
        self.assertFalse(logs.duplicated(['conj_id', 'error_model']).any())
        self.assertTrue(np.isfinite(logs[['log_pc', 'contribution_ess']]).all().all())
        joined = results.merge(logs, on=['conj_id', 'error_model', 'n_samples'], validate='one_to_one')
        np.testing.assert_allclose(joined.pc, np.exp(joined.log_pc), rtol=1e-8, atol=0)
        self.assertTrue((joined.numerical_underflow == (joined.pc == 0)).all())
        row = logs[(logs.conj_id == 801) & (logs.error_model == 'gaussian')].iloc[0]
        self.assertAlmostEqual(row.log_pc, -745.0924876472135, places=7)
        self.assertAlmostEqual(row.contribution_ess, 7270.4022062, places=4)

    def test_extreme_gaussian_contribution_diagnostics(self):
        moments = g.importance.Moments()
        moments.update(np.array([1e-200, 2e-200]))
        summary = moments.summary()
        self.assertAlmostEqual(summary['pc'] / 1e-200, 1.5)
        self.assertTrue(np.isnan(summary['contribution_ess']))

    def test_existing_inputs_and_mc_intervals(self):
        exact = pd.read_csv(g.PROCESSED / 'pc_exact.csv')
        mc = pd.read_csv(g.PROCESSED / 'pc_mc.csv')
        ids = set(pd.read_csv(g.PROCESSED / 'conjunctions.csv').conj_id)
        g.validate_results(exact, ids, [('exact', 'gaussian')])
        g.validate_results(mc, ids, [('mc', model) for model in g.mc.MODELS])
        g.validate_mc(mc)

    def test_reject_invalid_schema_duplicates_coverage(self):
        mean, lo, hi = g.mc.jeffreys_summary(0, 100)
        frame = pd.DataFrame([(1, 'mc', 'gaussian', mean, lo, hi, 100, .1)], columns=g.mc.RESULT_COLUMNS)
        g.validate_results(frame, {1}, [('mc', 'gaussian')])
        for invalid in [frame.iloc[:, ::-1], pd.concat([frame, frame]), frame.assign(pc=2)]:
            with self.assertRaises(ValueError):
                g.validate_results(invalid, {1}, [('mc', 'gaussian')])
        with self.assertRaises(ValueError):
            g.validate_results(frame, {1, 2}, [('mc', 'gaussian')])

    def test_bound_validates_existing_is_and_confidence_interval(self):
        row = pd.read_csv(g.PROCESSED / 'conjunctions.csv').query('conj_id == 183').iloc[0].to_dict()
        parameters = g.mc.load_parameters()
        for model in g.mc.MODELS:
            bound = g.contribution_bound(row, parameters[model], model)
            d = g.importance.estimate_pc(row, parameters[model], model, 1000, 100, np.random.default_rng(42))
            self.assertLessEqual(d['max_contribution'], bound)
            lo, hi = g.bayes.weighted_interval(d, bound)
            self.assertTrue(0 <= lo <= hi <= 1)

    def test_explicit_backups_and_no_implicit_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            frame = pd.DataFrame({'value': [1]})
            tables = {name: frame for name in g.OUTPUTS if name.endswith('.csv')}
            g.save_outputs(folder, tables, {'seed': 42})
            path = Path(folder) / 'pc_results.csv'
            original = path.read_bytes()
            with self.assertRaises(FileExistsError):
                g.save_outputs(folder, tables, {})
            self.assertEqual(path.read_bytes(), original)
            g.save_outputs(folder, tables, {'seed': 43}, True)
            backups = list((Path(folder) / 'backups').glob('*/pc_results.csv'))
            self.assertEqual(len(backups), 1)
            self.assertEqual(backups[0].read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
