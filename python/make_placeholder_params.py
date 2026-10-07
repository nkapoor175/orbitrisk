"""
make_placeholder_params.py  --  TEMPORARY STAND-IN (owner: Navika)

Vishi's r/10_error_model.R will write the real data/processed/params.csv. Until it exists, this script makes a quick
stand-in so that scripts 05 to 08 can be built and tested: for each direction it takes the robust spread
(1.4826 x median absolute deviation) of the residuals in age bins and fits  sigma = a * age^b  by a straight line on
log-log axes. It does NOT replace Vishi's model (no Student-t fit, no AIC). Do not use it for report results.

Reads   : data/processed/residuals.csv
Writes  : data/processed/params_placeholder.csv   (same columns as params.csv; NOT params.csv, so it can never
          be mistaken for, or overwrite, the real file)
Run from the repo root:  python python/make_placeholder_params.py
"""
import sys

import numpy as np
import pandas as pd

from common import PROCESSED

N_BINS = 8
COLS = {"radial": "d_radial_km", "along": "d_along_km", "cross": "d_cross_km"}


def main():
    src = PROCESSED / "residuals.csv"
    if not src.exists():
        sys.exit("data/processed/residuals.csv not found. Run python python/03_residuals.py first.")
    res = pd.read_csv(src)
    bins = pd.qcut(res["age_days"], N_BINS, duplicates="drop")
    rows = []
    for comp, col in COLS.items():
        ages, sigmas = [], []
        for _, g in res.groupby(bins, observed=True):
            x = g[col].to_numpy()
            ages.append(g["age_days"].median())
            sigmas.append(1.4826 * np.median(np.abs(x - np.median(x))))
        b, log_a = np.polyfit(np.log(ages), np.log(sigmas), 1)
        a = float(np.exp(log_a))
        print(f"  {comp:7s} a = {a:.4f} km   b = {b:.3f}")
        rows.append({"component": comp, "model": "gaussian", "a": a, "b": b, "nu": np.nan, "aic": np.nan,
                     "n_obs": len(res)})
        rows.append({"component": comp, "model": "student_t", "a": a * 0.8, "b": b, "nu": 4.0, "aic": np.nan,
                     "n_obs": len(res)})
    out = pd.DataFrame(rows)
    out.to_csv(PROCESSED / "params_placeholder.csv", index=False, float_format="%.6g")
    print("Wrote data/processed/params_placeholder.csv (STAND-IN, not Vishi's fit).")


if __name__ == "__main__":
    main()
