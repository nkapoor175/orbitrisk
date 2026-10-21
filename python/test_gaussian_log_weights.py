"""Regression checks; run .venv/bin/python -B python/test_gaussian_log_weights.py."""
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from scipy.special import logsumexp

spec = importlib.util.spec_from_file_location('generation', Path(__file__).with_name('generate_results.py'))
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)


class LogWeightTests(unittest.TestCase):
    def test_extremely_small_contributions_and_ess(self):
        logs = np.array([-10000., -10001., -10002., -10003.])
        moments = g.importance.LogMoments()
        moments.update(logs[:2]); moments.update(logs[2:])
        d = moments.summary()
        self.assertEqual(d['pc'], 0)
        self.assertTrue(d['numerical_underflow'])
        self.assertFalse(d['zero_contributions'])
        self.assertAlmostEqual(d['log_pc'], logsumexp(logs) - np.log(4), places=10)
        self.assertAlmostEqual(d['contribution_ess'], np.exp(2 * logsumexp(logs) - logsumexp(2 * logs)), places=9)
        low, high = g.bayes.weighted_interval(d, .01)
        self.assertEqual(low, 0)
        self.assertGreater(high, 0)

    def test_excluded_events_are_distinct_from_underflow(self):
        moments = g.importance.LogMoments()
        moments.update(np.full(10, -np.inf))
        d = moments.summary()
        self.assertTrue(d['zero_contributions'])
        self.assertFalse(d['numerical_underflow'])
        moments.update(np.array([-1000., -1001.]))
        self.assertTrue(np.isfinite(moments.summary()['log_pc']))

    def test_801_chunk_averaging_regression(self):
        row = pd.read_csv(g.PROCESSED / 'conjunctions.csv').query('conj_id == 801').iloc[0].to_dict()
        parameters = g.mc.load_parameters()['gaussian']
        results = []
        for chunk in [31, 1000, 10000]:
            rng = np.random.default_rng(np.random.SeedSequence([42, 801, 0, 7]))
            d = g.importance.estimate_pc(row, parameters, 'gaussian', 10000, chunk, rng)
            self.assertAlmostEqual(d['log_pc'], -745.0924876472135, places=9)
            self.assertEqual(d['pc'], np.nextafter(0., 1.))
            self.assertGreater(d['contribution_ess'], 7000)
            self.assertFalse(d['numerical_underflow'])
            results.append(d)
        np.testing.assert_allclose([d['contribution_ess'] for d in results], results[0]['contribution_ess'], rtol=1e-10)

    def test_scaled_variance_and_normal_weights(self):
        logs = np.log(np.array([1., 2., 3., 4.])) - 400
        moments = g.importance.LogMoments(); moments.update(logs)
        d = moments.summary()
        self.assertGreater(d['standard_error'], 0)
        self.assertAlmostEqual(d['standard_error'] / np.exp(-400), np.std([1., 2., 3., 4.], ddof=1)/2)
        self.assertEqual(d['contribution_variance'], 0) # unavoidable final variance representation

    def test_student_results_and_csv_contracts_unchanged(self):
        manifest = json.loads((g.PROCESSED / 'result_generation.json').read_text())
        params = g.mc.load_parameters()['student_t']
        row = pd.read_csv(g.PROCESSED / 'conjunctions.csv').query('conj_id == 183').iloc[0].to_dict()
        rng = np.random.default_rng(np.random.SeedSequence([42, 183, 1, 7]))
        d = g.importance.estimate_pc(row, params, 'student_t', 10000, 1000, rng)
        old = next(x for x in manifest['focused_is_diagnostics'] if x['conj_id']==183 and x['model']=='student_t')
        for key in ['pc', 'contribution_variance', 'contribution_ess', 'max_contribution']:
            np.testing.assert_allclose(d[key], old[key], rtol=1e-12, atol=0,
                                       err_msg=f"Student-t diagnostic changed: {key}")
        for name in ['pc_is.csv', 'pc_results.csv']:
            self.assertEqual(tuple(pd.read_csv(g.PROCESSED/name, nrows=0).columns), g.mc.RESULT_COLUMNS)
        self.assertEqual(tuple(pd.read_csv(g.PROCESSED/'mc_convergence.csv', nrows=0).columns), g.CONVERGENCE_COLUMNS)


if __name__ == '__main__':
    unittest.main()
