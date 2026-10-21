"""Run from root: .venv/bin/python -B python/test_importance_sampling.py."""
import importlib.util
from pathlib import Path
import unittest
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.stats import t


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


is_module = load("importance", "07_importance_sampling.py")
exact = load("exact", "05_exact_pc.py")
mc_test = load("mc_test", "test_monte_carlo.py")


class ImportanceTests(unittest.TestCase):
    def estimate(self, row, params, model, n=100_000, chunk=7919, seed=42):
        return is_module.estimate_pc(row, params, model, n, chunk, np.random.default_rng(seed))

    def test_uniform_disk(self):
        points = is_module.sample_disk(np.random.default_rng(1), 100_000, 2)
        radius2 = np.sum(points ** 2, axis=1)
        self.assertTrue((radius2 <= 4).all())
        self.assertLess(abs(radius2.mean() - 2), .02)
        np.testing.assert_allclose(points.mean(axis=0), [0, 0], atol=.015)

    def test_online_variance(self):
        values = np.array([.01, .02, .03, .04, .1])
        moments = is_module.Moments()
        moments.update(values[:2])
        moments.update(values[2:])
        result = moments.summary()
        self.assertAlmostEqual(result["pc"], values.mean())
        self.assertAlmostEqual(result["contribution_variance"], values.var(ddof=1))
        self.assertAlmostEqual(result["contribution_ess"], values.sum() ** 2 / np.square(values).sum())

    def test_cli_contract_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "pc_is.csv"
            command = [sys.executable, "-B", str(Path(__file__).with_name("07_importance_sampling.py")),
                       "--limit", "1", "--n-samples", "100", "--chunk-size", "31", "--output", str(output)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            data = pd.read_csv(output)
            self.assertEqual(tuple(data.columns), is_module.RESULT_COLUMNS)
            self.assertEqual(set(data.error_model), set(is_module.MODELS))
            self.assertTrue(data[["ci_low", "ci_high"]].isna().all().all())
            original = output.read_bytes()
            refused = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(refused.returncode, 0)
            self.assertEqual(output.read_bytes(), original)

    def test_gaussian_known_probability_and_variance(self):
        row = mc_test.conjunction()
        params = (np.full(3, .1), np.zeros(3), np.full(3, np.nan))
        result = self.estimate(row, params, "gaussian")
        truth = -np.expm1(-.1 ** 2 / (4 * .1 ** 2))
        self.assertLess(abs(result["pc"] - truth), 5 * result["standard_error"])
        self.assertLess(result["contribution_variance"], truth * (1 - truth))
        self.assertAlmostEqual(result["estimator_variance"], result["standard_error"] ** 2)
        self.assertGreater(result["contribution_ess"], 1)

    def test_student_weights_jacobian_and_target_density(self):
        matrix = np.array([[2., 0, 0, 0, 0, 0], [0, 3., 0, 0, 0, 0]])
        pivots, remaining, jacobian = is_module.select_pivots(matrix)
        self.assertAlmostEqual(jacobian, 6)
        points = np.array([[0., 0.], [.2, -.3]])
        scales, nu = np.ones(6), np.ones(6)
        weights = is_module.student_contributions(points, np.zeros((2, 4)), np.zeros(2),
                                                  matrix, scales, nu, pivots, remaining, 1)
        want = np.pi / 6 * t.pdf(points[:, 0] / 2, 1) * t.pdf(points[:, 1] / 3, 1)
        np.testing.assert_allclose(weights, want, rtol=1e-14)

    def test_student_cauchy_known_probability(self):
        # Geometry projects independent radial and cross components. Differences
        # of independent Cauchy(s) variables are Cauchy(2s); along drops out.
        row = mc_test.conjunction()
        params = (np.full(3, .1), np.zeros(3), np.ones(3))
        radius, scale = .1, .2
        truth = quad(lambda x: t.pdf(x / scale, 1) / scale *
                     (2 * t.cdf(np.sqrt(max(0, radius ** 2 - x ** 2)) / scale, 1) - 1),
                     -radius, radius, epsabs=1e-12)[0]
        result = self.estimate(row, params, "student_t", n=300_000)
        self.assertLess(abs(result["pc"] - truth), 5 * result["standard_error"])

    def test_reproducibility_chunks_and_checkpoints(self):
        row = mc_test.conjunction()
        for model in is_module.MODELS:
            params = (np.full(3, .1), np.ones(3), np.array([7.14, 1.0003947, 4.45]))
            a = self.estimate(row, params, model, n=1001, chunk=31)
            b = self.estimate(row, params, model, n=1001, chunk=1001)
            self.assertAlmostEqual(a["pc"], b["pc"], places=14)
            self.assertAlmostEqual(a["estimator_variance"], b["estimator_variance"], places=14)
        checkpoints = []
        is_module.estimate_pc(row, params, "student_t", 100, 31,
                              np.random.default_rng(1), checkpoints.append)
        self.assertEqual([d["n_samples"] for d in checkpoints], [31, 62, 93, 100])

    def test_actual_gaussian_against_recomputed_exact(self):
        conjunctions = pd.read_csv(is_module.PROCESSED / "conjunctions.csv")
        params = is_module.load_parameters()["gaussian"]
        sig = {c: (params[0][i], params[1][i]) for i, c in enumerate(is_module._mc.COMPONENTS)}
        for conj_id in [183, 911, 1198]:
            selected = conjunctions[conjunctions.conj_id == conj_id]
            truth = exact.compute_pc(selected, sig)[0][0]
            result = self.estimate(selected.iloc[0].to_dict(), params, "gaussian")
            self.assertLess(abs(result["pc"] - truth), 5 * result["standard_error"])

    def test_saved_exact_matches_current_inputs(self):
        conjunctions = pd.read_csv(is_module.PROCESSED / "conjunctions.csv")
        saved = pd.read_csv(is_module.PROCESSED / "pc_exact.csv").set_index("conj_id")
        params = is_module.load_parameters()["gaussian"]
        sig = {c: (params[0][i], params[1][i]) for i, c in enumerate(is_module._mc.COMPONENTS)}
        for conj_id in [183, 911, 1198]:
            with self.subTest(conj_id=conj_id):
                truth = exact.compute_pc(conjunctions[conjunctions.conj_id == conj_id], sig)[0][0]
                np.testing.assert_allclose(truth, saved.loc[conj_id, "pc"], rtol=1e-6, atol=0)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            is_module.select_pivots(np.zeros((2, 6)))
        with self.assertRaises(ValueError):
            self.estimate(mc_test.conjunction(), None, "invalid")


if __name__ == "__main__":
    unittest.main()
