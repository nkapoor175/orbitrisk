"""
07_importance_sampling.py  --  OWNER: Parthvi

Purpose : Importance-sampled Pc for rare events; compare variance and runtime with plain Monte Carlo.
Reads   : data/processed/conjunctions.csv, data/processed/params.csv
Writes  : data/processed/pc_is.csv
Proposal: uniform collision disk. Gaussian contributions are area * target
density. Student-t uses an exact change of variables in the six independent
RTN errors: draw four components from scale*t(nu), solve two for the disk
position, and include their target densities and the inverse Jacobian.
No Student-t covariance or multivariate-t approximation is used.

pc is the unbiased mean contribution, not a posterior mean. ci_low/ci_high
remain blank: the contract calls these credible intervals, and no Bayesian
interval for weighted samples is asserted. Diagnostics expose sample variance,
standard error, contribution ESS, and cumulative convergence via a callback;
the CLI prints them without inventing additional CSV columns.
Existing output files are never overwritten.
Run from the repo root:  python python/07_importance_sampling.py
"""
import argparse
import importlib.util
from itertools import combinations
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal, t
from scipy.special import logsumexp

from common import PROCESSED, RANDOM_SEED
from orbit_utils import encounter_basis, rtn_matrix, project_covariance_to_plane

# Reuse Parthvi's validated parameter ordering, age clamp, and output columns.
_spec = importlib.util.spec_from_file_location("orbitrisk_mc", Path(__file__).with_name("06_monte_carlo.py"))
_mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mc)
load_parameters = _mc.load_parameters
scales_at_age = _mc.scales_at_age
MODELS = _mc.MODELS
RESULT_COLUMNS = _mc.RESULT_COLUMNS


def sample_disk(rng, size, radius):
    """Uniform area sampling: r = R sqrt(U), theta = 2 pi V."""
    if size <= 0 or not np.isfinite(radius) or radius <= 0:
        raise ValueError("Positive size and finite positive radius required")
    uv = rng.random((size, 2))
    angle = 2 * np.pi * uv[:, 1]
    radial = radius * np.sqrt(uv[:, 0])
    return np.column_stack([radial * np.cos(angle), radial * np.sin(angle)])


def select_pivots(matrix):
    """Choose the best-conditioned area projection from two RTN components."""
    pairs = list(combinations(range(6), 2))
    determinants = [abs(np.linalg.det(matrix[:, pair])) for pair in pairs]
    pivots = np.array(pairs[int(np.argmax(determinants))])
    jacobian = max(determinants)
    if not np.isfinite(jacobian) or jacobian <= 1e-12:
        raise ValueError("Degenerate RTN-to-encounter projection")
    remaining = np.array([i for i in range(6) if i not in pivots])
    return pivots, remaining, jacobian


def student_contributions(points, free_errors, mean, matrix, scales, nu,
                          pivots, remaining, radius):
    """Exact target/proposal weights after solving the two pivot errors.

    x = mean + A_p e_p + A_f e_f. With e_f drawn from its target and
    x uniform in the disk, q(e) = product(f_f(e_f))*|det(A_p)|/area.
    Thus I(collision)*p(e)/q(e) = area*product(f_p(e_p))/|det(A_p)|.
    Proposal support covers every error vector that causes a collision.
    """
    pivot_errors = np.linalg.solve(matrix[:, pivots],
                                  (points - mean - free_errors @ matrix[:, remaining].T).T).T
    log_density = t.logpdf(pivot_errors / scales[pivots], df=nu[pivots]) - np.log(scales[pivots])
    log_weight = (np.log(np.pi) + 2 * np.log(radius)
                  - np.log(abs(np.linalg.det(matrix[:, pivots]))) + log_density.sum(axis=1))
    return np.exp(log_weight)


class LogMoments:
    """Gaussian contributions accumulated before float64 conversion.

    -inf denotes an excluded collision event, not density underflow. Finite
    log_pc with pc==0 explicitly identifies final representation underflow.
    Scaled moments retain variance/SE even when ordinary squared weights vanish.
    """
    def __init__(self):
        self.log_sum = -np.inf
        self.log_sum2 = -np.inf
        self.scale = -np.inf
        self.scaled = Moments()

    def update(self, logs):
        logs = np.asarray(logs, float)
        if not len(logs) or np.isnan(logs).any() or np.isposinf(logs).any():
            raise FloatingPointError("Invalid Gaussian log contributions")
        self.log_sum = np.logaddexp(self.log_sum, logsumexp(logs))
        self.log_sum2 = np.logaddexp(self.log_sum2, logsumexp(2 * logs))
        peak = float(logs.max())
        if peak > self.scale:
            factor = np.exp(self.scale - peak)
            self.scaled.mean *= factor
            self.scaled.m2 *= factor ** 2
            self.scaled.maximum *= factor
            self.scale = peak
        values = np.exp(logs - self.scale) if np.isfinite(self.scale) else np.zeros(len(logs))
        self.scaled.update(values)

    def summary(self):
        n = self.scaled.n
        log_pc = float(self.log_sum - np.log(n))
        pc = float(np.exp(log_pc))
        variance = self.scaled.m2 / (n - 1) if n > 1 else np.nan
        log_variance = (2 * self.scale + np.log(variance)) if variance > 0 else (-np.inf if variance == 0 else np.nan)
        return dict(n_samples=n, pc=pc, log_pc=log_pc,
                    numerical_underflow=bool(pc == 0 and np.isfinite(log_pc)),
                    zero_contributions=bool(np.isneginf(log_pc)),
                    contribution_variance=float(np.exp(log_variance)),
                    estimator_variance=float(np.exp(log_variance - np.log(n))),
                    standard_error=float(np.exp(.5 * (log_variance - np.log(n)))),
                    log_contribution_variance=float(log_variance),
                    contribution_ess=float(min(n, np.exp(2 * self.log_sum - self.log_sum2))) if np.isfinite(log_pc) else 0.,
                    max_contribution=float(np.exp(self.scale)) if np.isfinite(self.scale) else 0.)


class Moments:
    """Stable online moments of iid weighted contributions (not Bernoulli)."""
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0
        self.maximum = 0.0

    def update(self, values):
        if not np.isfinite(values).all() or (values < 0).any():
            raise FloatingPointError("Nonfinite or negative importance contribution")
        size, mean = len(values), float(values.mean())
        delta = mean - self.mean
        total = self.n + size
        self.m2 += float(np.sum((values - mean) ** 2)) + delta ** 2 * self.n * size / total
        self.mean += delta * size / total
        self.n = total
        self.maximum = max(self.maximum, float(values.max()))

    def summary(self):
        variance = self.m2 / (self.n - 1) if self.n > 1 else np.nan
        # Contribution ESS is a concentration diagnostic, not a binomial N.
        if self.mean == 0:
            ess = 0.0
        elif self.mean < np.sqrt(np.finfo(float).tiny):
            # Squared moments are not representable at this scale.
            ess = np.nan
        else:
            ess = self.n / (1 + self.m2 / self.n / self.mean ** 2)
        return dict(n_samples=self.n, pc=self.mean, contribution_variance=variance,
                    estimator_variance=variance / self.n,
                    standard_error=np.sqrt(variance / self.n), contribution_ess=ess,
                    max_contribution=self.maximum)


def estimate_pc(row, parameters, model, n_samples, chunk_size, rng, callback=None):
    """Unbiased disk-proposal estimator; callback receives cumulative diagnostics.

    Both models have bounded contributions (finite density maxima and fixed
    nonzero Jacobian), even when RTN errors have infinite variance. The standard
    error estimates Var(contribution)/N. Normal intervals would be asymptotic,
    not Bayesian; interval validation/integration is deliberately deferred.
    """
    if model not in MODELS or n_samples < 2 or chunk_size <= 0:
        raise ValueError("Known model, at least two samples, positive chunk size required")
    started = time.perf_counter()
    def state(prefix, unit):
        return np.array([row[f"{prefix}_{axis}_{unit}"] for axis in "xyz"], float)
    r1, v1 = state("r1", "km"), state("v1", "kms")
    r2, v2 = state("r2", "km"), state("v2", "kms")
    mean = np.array([row["enc_x_km"], row["enc_y_km"]], float)
    radius = float(row["radius_m"]) / 1000
    if not np.isfinite(np.concatenate([r1, v1, r2, v2, mean])).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError("Finite states/miss vector and positive finite radius required")
    s1 = scales_at_age(parameters, row["age1_days"])
    s2 = scales_at_age(parameters, row["age2_days"])
    if model == "gaussian":
        cov = project_covariance_to_plane(r1, v1, s1, r2, v2, s2)
        density = multivariate_normal(mean=mean, cov=cov)
    else:
        x, y, _ = encounter_basis(r1, v1, r2, v2)
        plane = np.vstack([x, y])
        matrix = np.column_stack([-plane @ rtn_matrix(r1, v1), plane @ rtn_matrix(r2, v2)])
        pivots, remaining, _ = select_pivots(matrix)
        scales, nu = np.concatenate([s1, s2]), np.tile(parameters[2], 2)
    # Separate streams ensure chunk boundaries do not change paired disk/free draws.
    disk_rng, free_rng = [np.random.default_rng(int(seed)) for seed in
                          rng.integers(0, 2**63, size=2)]
    moments = LogMoments() if model == "gaussian" else Moments()
    for start in range(0, n_samples, chunk_size):
        size = min(chunk_size, n_samples - start)
        points = sample_disk(disk_rng, size, radius)
        if model == "gaussian":
            values = np.log(np.pi) + 2 * np.log(radius) + density.logpdf(points)
            if not np.isfinite(values).all():
                raise FloatingPointError("Gaussian log density is not finite")
        else:
            free = scales[remaining] * free_rng.standard_t(nu[remaining], size=(size, 4))
            values = student_contributions(points, free, mean, matrix, scales, nu,
                                           pivots, remaining, radius)
        moments.update(values)
        if callback is not None:
            callback(dict(moments.summary(), runtime_s=time.perf_counter() - started))
    return dict(moments.summary(), runtime_s=time.perf_counter() - started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--model", choices=("both",) + MODELS, default="both")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, default=PROCESSED / "pc_is.csv")
    args = parser.parse_args()
    if args.n_samples < 2 or args.chunk_size <= 0 or args.seed < 0 or (args.limit is not None and args.limit <= 0):
        parser.error("At least two samples, positive chunk/limit and nonnegative seed required")
    protected = [PROCESSED / name for name in
                 ("pc_exact.csv", "pc_mc.csv", "pc_results.csv", "mc_convergence.csv")]
    if args.output.exists() or args.output.resolve() in [p.resolve() for p in protected]:
        parser.error("Output exists or is reserved; choose a new --output path")
    parameters = load_parameters()
    conjunctions = pd.read_csv(PROCESSED / "conjunctions.csv")
    if conjunctions.empty or conjunctions.conj_id.duplicated().any():
        raise ValueError("Nonempty conjunction table with unique IDs required")
    if args.limit is not None:
        conjunctions = conjunctions.head(args.limit)
    results = []
    for row in conjunctions.to_dict("records"):
        conj_id = int(row["conj_id"])
        if conj_id < 0 or conj_id != row["conj_id"]:
            raise ValueError("conj_id must be a nonnegative integer")
        for model in (MODELS if args.model == "both" else (args.model,)):
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, conj_id, MODELS.index(model), 7]))
            diagnostics = estimate_pc(row, parameters[model], model, args.n_samples, args.chunk_size, rng)
            print(f"{conj_id} {model}: {diagnostics}", flush=True)
            results.append((conj_id, "is", model, diagnostics["pc"], np.nan, np.nan,
                            args.n_samples, diagnostics["runtime_s"]))
    with args.output.open("x") as output:
        pd.DataFrame(results, columns=RESULT_COLUMNS).to_csv(output, index=False, float_format="%.10e")
    print(f"Wrote {len(results)} rows to {args.output}; credible intervals intentionally blank")


if __name__ == "__main__":
    main()
