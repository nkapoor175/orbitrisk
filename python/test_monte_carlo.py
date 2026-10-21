"""Focused MC checks. Run: python -B python/test_monte_carlo.py."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.stats import beta

spec = importlib.util.spec_from_file_location("monte_carlo", Path(__file__).with_name("06_monte_carlo.py"))
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)


def conjunction():
    row = dict(enc_x_km=0.0, enc_y_km=0.0, radius_m=100.0,
               age1_days=0.01, age2_days=2.0)
    for prefix, values, unit in [("r1", [7000, 0, 0], "km"),
                                 ("r2", [7000, 0, 0], "km"),
                                 ("v1", [0, 7, 0], "kms"),
                                 ("v2", [0, -7, 0], "kms")]:
        row.update({f"{prefix}_{axis}_{unit}": value for axis, value in zip("xyz", values)})
    return row


class MonteCarloTests(unittest.TestCase):
    def test_jeffreys_zero_and_all_hits(self):
        for hits in [0, 17]:
            mean, low, high = mc.jeffreys_summary(hits, 17)
            self.assertEqual(mean, (hits + 0.5) / 18)
            np.testing.assert_allclose([low, high], beta.ppf([0.025, 0.975], hits + .5, 17 - hits + .5))
            self.assertTrue(0 < low < high < 1)

    def test_individual_ages_and_clamp(self):
        params = (np.array([1., 2., 3.]), np.ones(3), np.ones(3))
        np.testing.assert_allclose(mc.scales_at_age(params, .01), [.1, .2, .3])
        np.testing.assert_allclose(mc.scales_at_age(params, 2), [2, 4, 6])

    def test_gaussian_centered_circle_and_chunk_reproducibility(self):
        row = conjunction()
        params = (np.full(3, .1), np.zeros(3), np.full(3, np.nan))
        n = 200_000
        hits = mc.count_collisions(row, params, "gaussian", n, 7919, np.random.default_rng(42))
        repeat = mc.count_collisions(row, params, "gaussian", n, n, np.random.default_rng(42))
        self.assertEqual(hits, repeat)
        # Two independent isotropic objects: combined sigma^2 = 2 * .1^2.
        exact = -np.expm1(-.1 ** 2 / (2 * 2 * .1 ** 2))
        self.assertLess(abs(hits / n - exact), 5 * np.sqrt(exact * (1 - exact) / n))

    def test_student_t_direct_scales_projection_and_inclusive_boundary(self):
        row = conjunction()
        row["radius_m"] = 2000.0
        params = (np.ones(3), np.ones(3), np.array([7.14, 1.0003947, 4.45]))

        class FixedDraws:
            sizes = []

            def standard_t(self, nu, size):
                np.testing.assert_array_equal(nu, params[2])
                self.sizes.append(size[0])
                errors = np.zeros(size)
                errors[:, 1, 0] = 1.0  # object 2 radial scale = 2 km
                return errors

        rng = FixedDraws()
        self.assertEqual(mc.count_collisions(row, params, "student_t", 7, 3, rng), 7)
        self.assertEqual(rng.sizes, [3, 3, 1])
        row["radius_m"] = 1999.0
        self.assertEqual(mc.count_collisions(row, params, "student_t", 7, 3, rng), 0)

    def test_fitted_parameters(self):
        models = mc.load_parameters()
        self.assertEqual(set(models), set(mc.MODELS))
        self.assertLess(models["student_t"][2][1], 2)  # infinite along-track variance


if __name__ == "__main__":
    unittest.main()
