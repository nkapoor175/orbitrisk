"""
06_monte_carlo.py  --  OWNER: Parthvi

Purpose : Plain Monte Carlo Pc under the Gaussian and the learned Student-t error model, with runtime.
Reads   : data/processed/conjunctions.csv, data/processed/params.csv
Writes  : data/processed/pc_mc.csv
Pc is the Jeffreys posterior mean; intervals are equal-tailed 95% Beta
credible intervals. Output uses the existing pc_results.csv column contract.
Existing output files are never overwritten. Runtime includes geometry and
simulation for each row. Student-t scales are not standard deviations.
Run from the repo root:  python python/06_monte_carlo.py
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from common import PROCESSED, RANDOM_SEED, MIN_PAIR_AGE_DAYS
from orbit_utils import encounter_basis, rtn_matrix, project_covariance_to_plane

COMPONENTS = ("radial", "along", "cross")
MODELS = ("gaussian", "student_t")
RESULT_COLUMNS = ("conj_id", "method", "error_model", "pc", "ci_low",
                  "ci_high", "n_samples", "runtime_s")


def load_parameters():
    """Read only fitted params.csv; require one row per component/model."""
    table = pd.read_csv(PROCESSED / "params.csv")
    models = {}
    for model in MODELS:
        rows = []
        for component in COMPONENTS:
            match = table[(table.model == model) & (table.component == component)]
            if len(match) != 1:
                raise ValueError(f"Expected one {model}/{component} parameter row")
            rows.append(match.iloc[0])
        a = np.array([r.a for r in rows], dtype=float)
        b = np.array([r.b for r in rows], dtype=float)
        nu = np.array([r.nu for r in rows], dtype=float)
        if not (np.isfinite(a).all() and (a > 0).all() and np.isfinite(b).all()):
            raise ValueError(f"Invalid scales or exponents for {model}")
        if model == "student_t" and not (np.isfinite(nu).all() and (nu > 0).all()):
            raise ValueError("Student-t degrees of freedom must be finite and positive")
        models[model] = (a, b, nu)
    return models


def scales_at_age(parameters, age_days):
    """Gaussian sigmas or Student-t scales, in radial/along/cross order."""
    a, b, _ = parameters
    if not np.isfinite(age_days):
        raise ValueError("TLE age must be finite")
    scales = a * max(float(age_days), MIN_PAIR_AGE_DAYS) ** b
    if not (np.isfinite(scales).all() and (scales > 0).all()):
        raise ValueError("Age-dependent scales must be finite and positive")
    return scales


def jeffreys_summary(hits, n_samples):
    """Posterior mean and equal-tailed interval for independent Bernoulli draws."""
    if n_samples <= 0 or not 0 <= hits <= n_samples:
        raise ValueError("Invalid collision count or sample count")
    alpha, beta_shape = hits + 0.5, n_samples - hits + 0.5
    low, high = beta.ppf([0.025, 0.975], alpha, beta_shape)
    return alpha / (n_samples + 1), float(low), float(high)


def count_collisions(row, parameters, model, n_samples, chunk_size, rng):
    """Sample relative position; return raw hits without retaining sample arrays.

    Student-t components and objects are independent. Each object's RTN error
    is rotated to TEME, then the relative error (object 2 minus object 1) is
    projected into the encounter plane. This is not a 2D multivariate t.
    """
    if model not in MODELS or n_samples <= 0 or chunk_size <= 0:
        raise ValueError("Unknown model or nonpositive sample/chunk count")
    def state(prefix, unit):
        return np.array([row[f"{prefix}_{axis}_{unit}"] for axis in "xyz"], float)

    r1, v1 = state("r1", "km"), state("v1", "kms")
    r2, v2 = state("r2", "km"), state("v2", "kms")
    mean = np.array([row["enc_x_km"], row["enc_y_km"]], float)
    radius = float(row["radius_m"]) / 1000.0
    if not (np.isfinite(np.concatenate([r1, v1, r2, v2, mean])).all()
            and np.isfinite(radius) and radius > 0):
        raise ValueError("States and miss vector must be finite; radius must be positive")
    s1 = scales_at_age(parameters, row["age1_days"])
    s2 = scales_at_age(parameters, row["age2_days"])
    if model == "gaussian":
        cov = project_covariance_to_plane(r1, v1, s1, r2, v2, s2)
        transform = np.linalg.cholesky(cov)
    else:
        x, y, _ = encounter_basis(r1, v1, r2, v2)
        plane = np.vstack([x, y])
        t1, t2 = plane @ rtn_matrix(r1, v1), plane @ rtn_matrix(r2, v2)
        if not (np.isfinite(t1).all() and np.isfinite(t2).all()):
            raise ValueError("Invalid encounter geometry")

    hits = 0
    for start in range(0, n_samples, chunk_size):
        size = min(chunk_size, n_samples - start)
        if model == "gaussian":
            positions = mean + rng.standard_normal((size, 2)) @ transform.T
        else:
            nu = parameters[2]
            # One independent univariate t per object and RTN component.
            errors = rng.standard_t(nu, size=(size, 2, 3))
            positions = mean + (errors[:, 1, :] * s2) @ t2.T - (errors[:, 0, :] * s1) @ t1.T
        hits += int(np.count_nonzero(np.hypot(positions[:, 0], positions[:, 1]) <= radius))
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--model", choices=("both",) + MODELS, default="both")
    parser.add_argument("--limit", type=int, help="First N conjunctions for a quick run")
    parser.add_argument("--output", type=Path, default=PROCESSED / "pc_mc.csv")
    args = parser.parse_args()
    if args.n_samples <= 0 or args.chunk_size <= 0 or args.seed < 0 or (args.limit is not None and args.limit <= 0):
        parser.error("Sample/chunk/limit counts must be positive; seed must be nonnegative")
    if args.output.exists():
        parser.error(f"Output already exists: {args.output}; choose a new --output path")
    models = load_parameters()
    conjunctions = pd.read_csv(PROCESSED / "conjunctions.csv")
    if conjunctions.empty or conjunctions.conj_id.duplicated().any():
        raise ValueError("Conjunctions must be nonempty with unique conj_id values")
    if args.limit is not None:
        conjunctions = conjunctions.head(args.limit)
    selected = MODELS if args.model == "both" else (args.model,)
    results = []
    for row in conjunctions.to_dict("records"):
        conj_id = int(row["conj_id"])
        if conj_id < 0 or conj_id != row["conj_id"]:
            raise ValueError("conj_id must be a nonnegative integer")
        for model in selected:
            # Stable streams across subsets and model selection; no shared draws.
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, conj_id, MODELS.index(model)]))
            started = time.perf_counter()
            hits = count_collisions(row, models[model], model, args.n_samples, args.chunk_size, rng)
            pc, low, high = jeffreys_summary(hits, args.n_samples)
            results.append((conj_id, "mc", model, pc, low, high, args.n_samples,
                            time.perf_counter() - started))
        if len(results) % (100 * len(selected)) == 0:
            print(f"Simulated {len(results) // len(selected)} conjunctions", flush=True)
    # Exclusive creation also protects against an output appearing during simulation.
    with args.output.open("x") as output:
        pd.DataFrame(results, columns=RESULT_COLUMNS).to_csv(output, index=False, float_format="%.10e")
    print(f"Wrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
