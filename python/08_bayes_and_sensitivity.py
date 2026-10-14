"""
08_bayes_and_sensitivity.py  --  OWNER: Parthvi

Purpose : MC credible intervals, weighted-IS confidence intervals, radius sensitivity.
Reads   : data/processed/conjunctions.csv, params.csv (fitted parameters only).
Writes  : pc_radius_sensitivity.csv; integrated pc_results.csv is deferred.

MC: iid Bernoulli indicators, Jeffreys Beta(1/2,1/2) prior, equal-tailed
95% posterior credible interval and posterior mean. Zero hits do not imply
zero probability. Coverage is Bayesian, not a finite-sample frequentist claim.

IS: iid nonnegative contributions bounded by a known density/Jacobian bound
B, with Hoeffding radius B*sqrt(log(2/alpha)/(2*N)). This is a finite-sample
95% frequentist confidence interval, not a Bayesian credible interval. It
remains conservative for zero contributions and may be uninformative [0,1].
The bound must be analytic, never the observed maximum or contribution ESS.
All intervals condition on fitted parameters and geometry: they omit fitting
uncertainty, model misspecification, and floating-point error. Extreme target
densities can underflow; a computed zero is not proof of zero risk.

Sensitivity uses shared samples and nested events. IS uses ONE largest-disk
proposal, so smaller-disk weights retain the largest disk area and Jacobian.
Radii are correlated experiments; individual intervals are marginal, not
simultaneous confidence bands. Only four contracted columns are exported;
intervals and diagnostics are returned through the Python API/printed by CLI.
Existing files are never overwritten.
Run from the repo root:  python python/08_bayes_and_sensitivity.py
"""
import argparse
import importlib.util
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal, t

from common import PROCESSED, RANDOM_SEED
from orbit_utils import encounter_basis, rtn_matrix, project_covariance_to_plane


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mc = _load("radius_mc", "06_monte_carlo.py")
importance = _load("radius_is", "07_importance_sampling.py")
RADIUS_COLUMNS = ("conj_id", "radius_m", "error_model", "pc")


def mc_interval(hits, n_samples):
    """Jeffreys posterior mean and 95% interval from raw integer counts only."""
    if (not isinstance(hits, (int, np.integer)) or
            not isinstance(n_samples, (int, np.integer))):
        raise ValueError("MC requires integer hit and sample counts")
    return mc.jeffreys_summary(hits, n_samples)


def weighted_interval(diagnostics, bound, confidence=.95):
    """Hoeffding CI for an iid, unbiased, bounded weighted mean.

    Diagnostics are Moments.summary()/estimate_pc() output from script 07.
    No conversion to binomial counts, self-normalization, or bootstrap of
    insufficient aggregate statistics. Known bound B applies to each draw.
    """
    n, mean = diagnostics["n_samples"], diagnostics["pc"]
    if (not isinstance(n, (int, np.integer)) or n < 1 or
            not np.isfinite(mean) or mean < 0 or
            not np.isfinite(bound) or bound <= 0 or not 0 < confidence < 1):
        raise ValueError("Invalid diagnostics, contribution bound, or confidence")
    if mean > bound * (1 + 1e-12) or diagnostics.get("max_contribution", 0) > bound * (1 + 1e-12):
        raise ValueError("Contributions exceed the declared bound")
    half_width = bound * np.sqrt(np.log(2 / (1 - confidence)) / (2 * n))
    # Conservative allowance for accumulation/conversion roundoff, including
    # final subnormal rounding. This does not certify scipy density evaluation
    # or covariance conditioning; those numerical-model errors remain excluded.
    rounding = 64 * np.finfo(float).eps * (n + 1) * bound + np.nextafter(0., 1.)
    half_width += rounding
    low = max(0., float(np.nextafter(mean - half_width, -np.inf)))
    high = min(1., float(np.nextafter(mean + half_width, np.inf)))
    if low > high:
        # Probability constraint contradicts the confidence set on this rare
        # realization. Returning the whole parameter space is conservative.
        return 0., 1.
    return low, high


def _geometry(row, parameters):
    def state(prefix, unit):
        return np.array([row[f"{prefix}_{axis}_{unit}"] for axis in "xyz"], float)
    r1, v1 = state("r1", "km"), state("v1", "kms")
    r2, v2 = state("r2", "km"), state("v2", "kms")
    mean = np.array([row["enc_x_km"], row["enc_y_km"]], float)
    if not np.isfinite(np.concatenate([r1, v1, r2, v2, mean])).all():
        raise ValueError("Nonfinite state or miss vector")
    s1 = mc.scales_at_age(parameters, row["age1_days"])
    s2 = mc.scales_at_age(parameters, row["age2_days"])
    return r1, v1, r2, v2, mean, s1, s2


def sensitivity(row, parameters, model, radii_m, n_samples, chunk_size, seed,
                method="is"):
    """Return four-column results and uncertainty diagnostics, without mutation.

    Geometry/densities are prepared once and only one chunk is held in memory.
    This extends the sampling loop from scripts 06/07 for multiple nested events
    while reusing their scales, proposal, weight, moments and interval helpers.
    """
    radii = np.unique(np.asarray(radii_m, float))
    if (radii.ndim != 1 or not len(radii) or not np.isfinite(radii).all() or
            (radii <= 0).any() or model not in mc.MODELS or method not in ("mc", "is") or
            n_samples < 2 or chunk_size <= 0 or seed < 0):
        raise ValueError("Invalid radii, model, method, sample/chunk count or seed")
    conj_id = int(row["conj_id"])
    if conj_id < 0 or conj_id != row["conj_id"]:
        raise ValueError("conj_id must be a nonnegative integer")
    started = time.perf_counter()
    r1, v1, r2, v2, mean, s1, s2 = _geometry(row, parameters)
    maximum_radius = radii[-1] / 1000
    if model == "gaussian":
        cov = project_covariance_to_plane(r1, v1, s1, r2, v2, s2)
        transform = np.linalg.cholesky(cov)
        density = multivariate_normal(mean=mean, cov=cov)
        bound = np.exp(np.log(np.pi) + 2 * np.log(maximum_radius) +
                       density.logpdf(mean))
    else:
        x, y, _ = encounter_basis(r1, v1, r2, v2)
        plane = np.vstack([x, y])
        matrix = np.column_stack([-plane @ rtn_matrix(r1, v1), plane @ rtn_matrix(r2, v2)])
        pivots, remaining, jacobian = importance.select_pivots(matrix)
        scales, nu = np.concatenate([s1, s2]), np.tile(parameters[2], 2)
        bound = np.exp(np.log(np.pi) + 2 * np.log(maximum_radius) - np.log(jacobian) +
                       (t.logpdf(0, nu[pivots]) - np.log(scales[pivots])).sum())
    if method == "is" and (not np.isfinite(bound) or bound <= 0):
        raise FloatingPointError("Contribution bound underflowed or overflowed")
    rng = np.random.default_rng(np.random.SeedSequence([seed, conj_id, mc.MODELS.index(model), 8]))
    disk_rng, free_rng = [np.random.default_rng(int(s)) for s in rng.integers(0, 2**63, size=2)]
    moments = [(importance.LogMoments() if model == "gaussian" else importance.Moments()) for _ in radii]
    hits = np.zeros(len(radii), dtype=np.int64)
    for start in range(0, n_samples, chunk_size):
        size = min(chunk_size, n_samples - start)
        if method == "mc":
            if model == "gaussian":
                positions = mean + rng.standard_normal((size, 2)) @ transform.T
            else:
                errors = rng.standard_t(parameters[2], size=(size, 2, 3))
                positions = mean + (errors.reshape(size, 6) * scales) @ matrix.T
        else:
            positions = importance.sample_disk(disk_rng, size, maximum_radius)
            if model == "gaussian":
                weights = np.log(np.pi) + 2 * np.log(maximum_radius) + density.logpdf(positions)
                if not np.isfinite(weights).all():
                    raise FloatingPointError("Gaussian log density is not finite")
            else:
                free = scales[remaining] * free_rng.standard_t(nu[remaining], size=(size, 4))
                weights = importance.student_contributions(positions, free, mean, matrix, scales, nu,
                                                           pivots, remaining, maximum_radius)
        distances = np.hypot(positions[:, 0], positions[:, 1])
        for i, radius in enumerate(radii):
            inside = distances <= radius / 1000
            if method == "mc":
                hits[i] += np.count_nonzero(inside)
            else:
                moments[i].update(np.where(inside, weights, -np.inf if model == "gaussian" else 0.))
    results, diagnostics = [], []
    for i, radius in enumerate(radii):
        if method == "mc":
            pc, low, high = mc_interval(int(hits[i]), n_samples)
            summary = dict(hits=int(hits[i]), n_samples=n_samples, pc=pc,
                           interval_kind="jeffreys_credible")
        else:
            summary = moments[i].summary()
            pc = summary["pc"]
            low, high = weighted_interval(summary, bound)
            summary.update(interval_kind="hoeffding_confidence", contribution_bound=bound)
        summary.update(radius_m=float(radius), interval_low=low, interval_high=high,
                       runtime_s=time.perf_counter() - started)
        diagnostics.append(summary)
        results.append((conj_id, float(radius), model, pc))
    return pd.DataFrame(results, columns=RADIUS_COLUMNS), diagnostics


def write_sensitivity(results, output):
    """Exclusive creation of only the contracted radius-sensitivity CSV."""
    output = Path(output)
    reserved = ("params.csv", "conjunctions.csv", "pc_exact.csv", "pc_mc.csv", "pc_is.csv",
                "pc_results.csv", "mc_convergence.csv")
    if output.name in reserved or output.exists():
        raise ValueError("Output exists or is reserved")
    if tuple(results.columns) != RADIUS_COLUMNS:
        raise ValueError("Incorrect radius-sensitivity columns")
    with output.open("x") as stream:
        results.to_csv(stream, index=False, float_format="%.10e")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radii-m", nargs="+", type=float, default=[1., 2.5, 5.])
    parser.add_argument("--method", choices=("mc", "is"), default="is")
    parser.add_argument("--model", choices=("both",) + mc.MODELS, default="both")
    parser.add_argument("--n-samples", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, default=PROCESSED / "pc_radius_sensitivity.csv")
    args = parser.parse_args()
    if (args.output.exists() or args.output.name in
            ("params.csv", "conjunctions.csv", "pc_exact.csv", "pc_mc.csv", "pc_is.csv",
             "pc_results.csv", "mc_convergence.csv")):
        parser.error("Output exists or is reserved")
    if args.n_samples < 2 or args.chunk_size <= 0 or args.seed < 0 or (args.limit is not None and args.limit <= 0):
        parser.error("Invalid sample/chunk/limit count or seed")
    if not np.isfinite(args.radii_m).all() or min(args.radii_m) <= 0:
        parser.error("Radii must be finite and positive")
    parameters = mc.load_parameters()
    conjunctions = pd.read_csv(PROCESSED / "conjunctions.csv")
    if conjunctions.empty or conjunctions.conj_id.duplicated().any():
        raise ValueError("Nonempty conjunctions with unique IDs required")
    if args.limit is not None:
        conjunctions = conjunctions.head(args.limit)
    # Always include the baseline and each selected row's current radius.
    radii = np.unique(args.radii_m + [2.5] + conjunctions.radius_m.tolist())
    frames = []
    for row in conjunctions.to_dict("records"):
        for model in (mc.MODELS if args.model == "both" else (args.model,)):
            frame, diagnostics = sensitivity(row, parameters[model], model, radii,
                                              args.n_samples, args.chunk_size, args.seed, args.method)
            frames.append(frame)
            print(f"{row['conj_id']} {model} {args.method}: {diagnostics}", flush=True)
    write_sensitivity(pd.concat(frames, ignore_index=True), args.output)
    print(f"Wrote radius sensitivity to {args.output}; integrated results deferred")


if __name__ == "__main__":
    main()
