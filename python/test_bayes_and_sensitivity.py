"""Run: .venv/bin/python -B python/test_bayes_and_sensitivity.py."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import beta


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


module = load("bayes_sensitivity", "08_bayes_and_sensitivity.py")
mc_tests = load("mc_fixture", "test_monte_carlo.py")


class UncertaintyTests(unittest.TestCase):
    def test_jeffreys_zero_few_many(self):
        for hits in [0, 2, 500, 1000]:
            mean, low, high = module.mc_interval(hits, 1000)
            self.assertAlmostEqual(mean, (hits + .5) / 1001)
            np.testing.assert_allclose([low, high], beta.ppf([.025, .975], hits + .5, 1000 - hits + .5))
            self.assertTrue(0 < low < high < 1)
        _, low1, high1 = module.mc_interval(0, 100)
        _, low2, high2 = module.mc_interval(0, 1000)
        self.assertGreater(high1, high2)
        self.assertGreater(low1, low2)
        with self.assertRaises(ValueError):
            module.mc_interval(.5, 100)

    def test_weighted_interval_and_zero_contributions(self):
        values = np.array([0., .2, .4, .8] * 1000)
        moments = module.importance.Moments()
        moments.update(values)
        low, high = module.weighted_interval(moments.summary(), 1.)
        width = np.sqrt(np.log(40) / (2 * len(values)))
        np.testing.assert_allclose([low, high], [values.mean() - width, values.mean() + width])
        low, high = module.weighted_interval(dict(n_samples=100, pc=0., max_contribution=0.), 1.)
        self.assertEqual(low, 0)
        self.assertGreater(high, 0)
        self.assertEqual(module.weighted_interval(dict(n_samples=100, pc=0.), 100), (0, 1))
        with self.assertRaises(ValueError):
            module.weighted_interval(moments.summary(), .1)

    def test_weighted_interval_frequentist_coverage_known_distribution(self):
        # Bounded non-Bernoulli weights with known mean .1.
        rng = np.random.default_rng(42)
        covered = 0
        for _ in range(200):
            values = rng.uniform(0, .2, 100)
            interval = module.weighted_interval(dict(n_samples=100, pc=values.mean()), .2)
            covered += interval[0] <= .1 <= interval[1]
        self.assertGreaterEqual(covered, 190)

    def test_nested_radii_and_input_preservation(self):
        row = dict(mc_tests.conjunction(), conj_id=1)
        original = copy.deepcopy(row)
        params = (np.full(3, .1), np.zeros(3), np.array([7.14, 1.0003947, 4.45]))
        saved_params = [p.copy() for p in params]
        for method in ["mc", "is"]:
            for model in module.mc.MODELS:
                frame, diagnostics = module.sensitivity(row, params, model, [50, 100, 200], 20_000, 7919, 42, method)
                self.assertEqual(tuple(frame.columns), module.RADIUS_COLUMNS)
                self.assertTrue((np.diff(frame.pc) >= 0).all())
                self.assertEqual(row, original)
                for a, b in zip(params, saved_params):
                    np.testing.assert_array_equal(a, b)
                for diagnostic in diagnostics:
                    self.assertLessEqual(diagnostic["interval_low"], diagnostic["interval_high"])
                    if method == "is":
                        self.assertLessEqual(diagnostic["max_contribution"], diagnostic["contribution_bound"])

    def test_gaussian_sensitivity_known_probability(self):
        row = dict(mc_tests.conjunction(), conj_id=1)
        params = (np.full(3, .1), np.zeros(3), np.full(3, np.nan))
        radii = np.array([50., 100., 200.])
        truth = -np.expm1(-(radii / 1000) ** 2 / .04)
        for method in ["mc", "is"]:
            frame, diagnostics = module.sensitivity(row, params, "gaussian", radii, 100_000, 7919, 7, method)
            np.testing.assert_allclose(frame.pc, truth, atol=.005, rtol=0)
            for expected, diagnostic in zip(truth, diagnostics):
                self.assertLessEqual(diagnostic["interval_low"], expected)
                self.assertGreaterEqual(diagnostic["interval_high"], expected)

    def test_seed_and_chunk_reproducibility(self):
        row = dict(mc_tests.conjunction(), conj_id=1)
        params = (np.full(3, .1), np.zeros(3), np.ones(3))
        for method in ["mc", "is"]:
            for model in module.mc.MODELS:
                a, _ = module.sensitivity(row, params, model, [50, 100], 1001, 31, 42, method)
                b, _ = module.sensitivity(row, params, model, [50, 100], 1001, 1001, 42, method)
                np.testing.assert_allclose(a.pc, b.pc, rtol=1e-13, atol=1e-15)

    def test_cli_output_contract_and_preservation(self):
        files = [module.PROCESSED / n for n in ["params.csv", "conjunctions.csv", "pc_exact.csv", "pc_mc.csv"]]
        hashes = {p: hashlib.sha256(p.read_bytes()).digest() for p in files if p.exists()}
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "sensitivity.csv"
            command = [sys.executable, "-B", str(Path(__file__).with_name("08_bayes_and_sensitivity.py")),
                       "--limit", "1", "--radii-m", "1", "5", "--n-samples", "100", "--chunk-size", "31",
                       "--output", str(output)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            frame = module.pd.read_csv(output)
            self.assertEqual(tuple(frame.columns), module.RADIUS_COLUMNS)
            self.assertIn(2.5, frame.radius_m.to_list())
            original = output.read_bytes()
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertEqual(original, output.read_bytes())
            for p, digest in hashes.items():
                self.assertEqual(hashlib.sha256(p.read_bytes()).digest(), digest)
            with self.assertRaises(ValueError):
                module.write_sensitivity(frame, Path(folder) / "pc_results.csv")


if __name__ == "__main__":
    unittest.main()
