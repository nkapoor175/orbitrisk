"""
05_exact_pc.py  --  OWNER: Navika

Purpose : Exact collision probability (Pc) of every conjunction under the GAUSSIAN position-error model.

          In the encounter plane the relative position of the two objects at closest approach is
          modelled as a 2D Gaussian: mean = the miss vector (enc_x_km, enc_y_km), covariance = the combined
          position covariance of both objects. A collision happens if the relative position falls inside the
          hard-body circle (radius = combined radius of the two objects). So

              Pc = integral over the circle  |p| <= R  of  N(p; miss vector, covariance) dp

          This 2D integral is computed numerically in polar coordinates (Gauss-Legendre in radius, equally
          spaced points in angle). For a smooth density this is accurate to many digits, so the result is
          treated as the "exact" value that Monte Carlo (script 06) and importance sampling (script 07) are
          validated against.

          How the covariance is built for each conjunction:
            sigma_j(age) = a_j * age^b_j  for j = radial, along-track, cross-track  (fitted in R by Vishi,
            Gaussian rows of params.csv), evaluated at each object's own TLE age (age1_days, age2_days).
            The two diagonal covariances are rotated into the encounter plane with
            project_covariance_to_plane() from orbit_utils.py.

          Notes for the report (limitations):
            * Only the Gaussian model has a closed-form integral. The learned Student-t model is handled by
              Monte Carlo / importance sampling in scripts 06 and 07.
            * Ages below MIN_PAIR_AGE_DAYS are raised to it (the power law was not fitted below that age).
              Ages above MAX_PAIR_AGE_DAYS are extrapolated and counted in the printed summary.
            * Straight-line relative motion and a constant covariance during the encounter are assumed
              (standard short-encounter assumption).

Reads   : data/processed/conjunctions.csv, data/processed/params.csv
Writes  : data/processed/pc_exact.csv   (columns as in pc_results.csv in docs/DATA_CONTRACT.md;
          ci_low, ci_high, n_samples are blank because the value is not an estimate;
          runtime_s = total run time divided equally between the rows)
Status  : DONE (first version). Until Vishi's params.csv exists, test with the stand-in file made by
          python/make_placeholder_params.py:   python python/05_exact_pc.py --params data/processed/params_placeholder.csv
Run from the repo root:  python python/05_exact_pc.py
"""
import argparse
import sys
import time

import numpy as np
import pandas as pd

from common import PROCESSED, MIN_PAIR_AGE_DAYS, MAX_PAIR_AGE_DAYS, ALERT_THRESHOLD
from orbit_utils import project_covariance_to_plane

COMPONENTS = ["radial", "along", "cross"]
CHUNK_POINTS = 4_000_000          # max number of integration points held in memory at once


def pc_gaussian_circle(mx, my, cov, radius_km, n_r=64, n_theta=256):
    """Probability that a 2D Gaussian N((mx, my), cov) falls inside the circle |p| <= radius_km.

    mx, my : arrays (n,) mean of the relative position in the encounter plane (km)
    cov    : array (n, 2, 2) covariance (km^2)
    radius_km : scalar or array (n,)
    Returns an array (n,) of probabilities.
    """
    mx = np.atleast_1d(np.asarray(mx, float))
    my = np.atleast_1d(np.asarray(my, float))
    cov = np.asarray(cov, float).reshape(-1, 2, 2)
    n = len(mx)
    radius = np.broadcast_to(np.asarray(radius_km, float), (n,))

    xr, wr = np.polynomial.legendre.leggauss(n_r)
    r_unit = 0.5 * (xr + 1.0)                      # radial nodes on [0, 1]
    w_r = 0.5 * wr
    theta = 2.0 * np.pi * (np.arange(n_theta) + 0.5) / n_theta
    cos_t, sin_t = np.cos(theta)[None, None, :], np.sin(theta)[None, None, :]
    w_theta = 2.0 * np.pi / n_theta

    c00, c01, c11 = cov[:, 0, 0], cov[:, 0, 1], cov[:, 1, 1]
    det = c00 * c11 - c01 ** 2
    i00, i01, i11 = c11 / det, -c01 / det, c00 / det
    norm = 1.0 / (2.0 * np.pi * np.sqrt(det))

    out = np.empty(n)
    step = max(1, CHUNK_POINTS // (n_r * n_theta))
    for s in range(0, n, step):
        e = min(n, s + step)
        sl = slice(s, e)
        r = radius[sl, None] * r_unit[None, :]                       # (m, n_r) radii in km
        dx = r[:, :, None] * cos_t - mx[sl, None, None]
        dy = r[:, :, None] * sin_t - my[sl, None, None]
        q = i00[sl, None, None] * dx ** 2 + 2.0 * i01[sl, None, None] * dx * dy + i11[sl, None, None] * dy ** 2
        dens = np.exp(-0.5 * q) * norm[sl, None, None]
        # area element: r dr dtheta, with dr = radius * d(unit radius)
        integrand = dens.sum(axis=2) * w_theta * r * (radius[sl, None] * w_r[None, :])
        out[sl] = integrand.sum(axis=1)
    return np.clip(out, 0.0, 1.0)


def load_sigma_model(params_path, model="gaussian"):
    """Read params.csv and return {component: (a, b)} for the chosen model."""
    p = pd.read_csv(params_path)
    p = p[p["model"] == model]
    out = {}
    for comp in COMPONENTS:
        row = p[p["component"] == comp]
        if row.empty:
            sys.exit(f"{params_path} has no '{model}' row for component '{comp}'.")
        out[comp] = (float(row["a"].iloc[0]), float(row["b"].iloc[0]))
    return out


def sigma_at_age(sig, age_days):
    """(n, 3) array of sigma_radial, sigma_along, sigma_cross in km for an array of TLE ages."""
    age = np.maximum(np.asarray(age_days, float), MIN_PAIR_AGE_DAYS)
    return np.column_stack([sig[c][0] * age ** sig[c][1] for c in COMPONENTS])


def conjunction_covariances(conj, sig):
    """Combined 2x2 encounter-plane covariance (n, 2, 2) in km^2 for every conjunction."""
    s1 = sigma_at_age(sig, conj["age1_days"].to_numpy())
    s2 = sigma_at_age(sig, conj["age2_days"].to_numpy())
    r1 = conj[["r1_x_km", "r1_y_km", "r1_z_km"]].to_numpy()
    v1 = conj[["v1_x_kms", "v1_y_kms", "v1_z_kms"]].to_numpy()
    r2 = conj[["r2_x_km", "r2_y_km", "r2_z_km"]].to_numpy()
    v2 = conj[["v2_x_kms", "v2_y_kms", "v2_z_kms"]].to_numpy()
    cov = np.empty((len(conj), 2, 2))
    for i in range(len(conj)):
        cov[i] = project_covariance_to_plane(r1[i], v1[i], s1[i], r2[i], v2[i], s2[i])
    return cov


def compute_pc(conj, sig):
    """Exact Gaussian Pc for every row of a conjunctions table. Returns (pc array, covariance array)."""
    cov = conjunction_covariances(conj, sig)
    radius_km = conj["radius_m"].to_numpy() / 1000.0
    # make the integration grid fine enough for the narrowest covariance direction
    tr, det = cov[:, 0, 0] + cov[:, 1, 1], cov[:, 0, 0] * cov[:, 1, 1] - cov[:, 0, 1] ** 2
    lam_min = tr / 2.0 - np.sqrt(np.maximum((tr / 2.0) ** 2 - det, 0.0))
    ratio = float(np.max(radius_km / np.sqrt(lam_min)))
    n_r = int(np.clip(np.ceil(6.0 * ratio), 64, 1024))
    pc = pc_gaussian_circle(conj["enc_x_km"].to_numpy(), conj["enc_y_km"].to_numpy(), cov, radius_km,
                            n_r=n_r, n_theta=4 * n_r)
    return pc, cov


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default=str(PROCESSED / "params.csv"),
                        help="fitted error model (default: data/processed/params.csv)")
    parser.add_argument("--limit", type=int, default=None, help="only the first N conjunctions (quick test)")
    args = parser.parse_args()

    conj_path = PROCESSED / "conjunctions.csv"
    if not conj_path.exists():
        sys.exit("data/processed/conjunctions.csv not found. Run python python/04_screen_conjunctions.py first.")
    from pathlib import Path
    params_path = Path(args.params)
    if not params_path.exists():
        sys.exit(f"{params_path} not found. Vishi's r/10_error_model.R writes data/processed/params.csv. To test "
                 "now, run python python/make_placeholder_params.py and add --params "
                 "data/processed/params_placeholder.csv")

    conj = pd.read_csv(conj_path)
    if args.limit:
        conj = conj.head(args.limit)
    sig = load_sigma_model(params_path)
    print("Gaussian error model (sigma = a * age^b, km):")
    for c in COMPONENTS:
        print(f"  {c:7s} a = {sig[c][0]:.4f}   b = {sig[c][1]:.3f}")

    t0 = time.perf_counter()
    pc, cov = compute_pc(conj, sig)
    elapsed = time.perf_counter() - t0

    res = pd.DataFrame({"conj_id": conj["conj_id"], "method": "exact", "error_model": "gaussian", "pc": pc,
                        "ci_low": np.nan, "ci_high": np.nan, "n_samples": np.nan,
                        "runtime_s": elapsed / len(conj)})
    res.to_csv(PROCESSED / "pc_exact.csv", index=False, float_format="%.6e")
    print(f"\nWrote data/processed/pc_exact.csv: {len(res)} conjunctions in {elapsed:.1f} s.")

    ages = np.concatenate([conj["age1_days"], conj["age2_days"]])
    print(f"TLE ages outside the fitted range {MIN_PAIR_AGE_DAYS}-{MAX_PAIR_AGE_DAYS} days: "
          f"{int((ages < MIN_PAIR_AGE_DAYS).sum())} below (raised), {int((ages > MAX_PAIR_AGE_DAYS).sum())} above (extrapolated).")
    print(f"Pc: min {pc.min():.2e}   median {np.median(pc):.2e}   max {pc.max():.2e}")
    print(f"Conjunctions with Pc >= {ALERT_THRESHOLD:g} (alerts): {int((pc >= ALERT_THRESHOLD).sum())}")
    top = conj.assign(pc=pc).sort_values("pc", ascending=False).head(5)
    print("\nFive highest-Pc conjunctions:")
    print(top[["conj_id", "norad_id_1", "norad_id_2", "tca_utc", "miss_distance_km", "pc"]].to_string(index=False))


if __name__ == "__main__":
    main()
