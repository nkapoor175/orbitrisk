"""
03_residuals.py  --  OWNER: Navika

Purpose : Build the position-error ("residual") data set that the learned error model is fitted on.

          For an object with TLEs at epochs t1 < t2 < t3 < ...:
            - take an OLDER TLE (epoch tk) and propagate it forward with SGP4 to the epoch tj of a NEWER TLE;
            - compare that predicted position with the position given by the newer TLE at its own epoch;
            - rotate the difference into radial / along-track / cross-track directions.
          age_days = tj - tk is the "TLE age": how old the prediction is.

          To cover a range of ages, each older TLE is paired with up to MAX_PAIRS_PER_TLE newer TLEs that lie
          between MIN_PAIR_AGE_DAYS and MAX_PAIR_AGE_DAYS later. Pairs are skipped if a manoeuvre happened in
          between (flag from 02_clean.py) or SGP4 reports an error.

          Notes for the report (limitations):
            * The newer TLE is used as the "truth" proxy, but it has its own error, so the residual contains
              the errors of BOTH TLEs. Small-age residuals are therefore slightly overestimated.
            * Pairs that share a TLE are not independent. When fitting, Vishi can thin the data (for example one
              pair per object per age bin) or use robust/cluster-aware methods.
Reads   : data/processed/tle_clean.csv
Writes  : data/processed/residuals.csv   (columns fixed in docs/DATA_CONTRACT.md)
Status  : DONE (first version)
Run from the repo root:  python python/03_residuals.py
"""
import sys

import numpy as np
import pandas as pd

from common import PROCESSED, RANDOM_SEED, MIN_PAIR_AGE_DAYS, MAX_PAIR_AGE_DAYS, MAX_PAIRS_PER_TLE, RESIDUAL_MAX_KM
from orbit_utils import satrec_from_lines, to_rtn


RNG = np.random.default_rng(RANDOM_SEED)


def residuals_for_object(g):
    """All residual rows for one object (g = its TLEs sorted by epoch)."""
    sats = [satrec_from_lines(a, b) for a, b in zip(g["line1"], g["line2"])]
    jd = np.array([s.jdsatepoch for s in sats])
    fr = np.array([s.jdsatepochF for s in sats])

    # position and velocity of every TLE at its OWN epoch (the "truth" proxy)
    own_ok = np.zeros(len(sats), dtype=bool)
    own_r = np.zeros((len(sats), 3))
    own_v = np.zeros((len(sats), 3))
    for i, s in enumerate(sats):
        e, r, v = s.sgp4(jd[i], fr[i])
        if e == 0:
            own_ok[i] = True
            own_r[i], own_v[i] = r, v

    jump_cum = np.cumsum(g["jump"].to_numpy())
    epochs = g["epoch_utc"].to_numpy()
    alt = g["altitude_km"].to_numpy()
    bstar = g["bstar"].to_numpy()

    out = []
    for k in range(len(sats) - 1):
        if not own_ok[k]:
            continue
        age = (jd[k + 1:] - jd[k]) + (fr[k + 1:] - fr[k])
        idx = np.where((age >= MIN_PAIR_AGE_DAYS) & (age <= MAX_PAIR_AGE_DAYS))[0] + k + 1
        idx = idx[own_ok[idx] & (jump_cum[idx] - jump_cum[k] == 0)]      # no manoeuvre in between
        if len(idx) == 0:
            continue
        if len(idx) > MAX_PAIRS_PER_TLE:                                  # random choice -> ages spread evenly
            idx = np.sort(RNG.choice(idx, size=MAX_PAIRS_PER_TLE, replace=False))

        err, r_pred, _ = sats[k].sgp4_array(jd[idx], fr[idx])             # older TLE propagated to newer epochs
        good = err == 0
        if not good.any():
            continue
        idx, r_pred = idx[good], r_pred[good]
        d = own_r[idx] - r_pred                                           # newer minus propagated
        d_rad, d_along, d_cross = to_rtn(d, own_r[idx], own_v[idx])
        age_days = (jd[idx] - jd[k]) + (fr[idx] - fr[k])
        for n, j in enumerate(idx):
            out.append((epochs[k], epochs[j], age_days[n], d_rad[n], d_along[n], d_cross[n], alt[j], bstar[k]))
    return out


def main():
    src = PROCESSED / "tle_clean.csv"
    if not src.exists():
        sys.exit("data/processed/tle_clean.csv not found. Run python python/02_clean.py first.")
    df = pd.read_csv(src, parse_dates=["epoch_utc"])
    df = df.sort_values(["norad_id", "epoch_utc"]).reset_index(drop=True)

    frames = []
    groups = df.groupby("norad_id")
    n_obj = groups.ngroups
    for count, (norad, g) in enumerate(groups, start=1):
        g = g.reset_index(drop=True)
        rows = residuals_for_object(g)
        if rows:
            part = pd.DataFrame(rows, columns=["epoch_old_utc", "epoch_new_utc", "age_days", "d_radial_km",
                                               "d_along_km", "d_cross_km", "altitude_km", "bstar"])
            part.insert(0, "object_type", g["object_type"].iloc[0])
            part.insert(0, "object_name", g["object_name"].iloc[0])
            part.insert(0, "norad_id", norad)
            frames.append(part)
        if count % 50 == 0 or count == n_obj:
            print(f"  processed {count}/{n_obj} objects")

    if not frames:
        sys.exit("No TLE pairs found. Check MIN/MAX_PAIR_AGE_DAYS in python/common.py or the downloaded data.")
    res = pd.concat(frames, ignore_index=True)

    big = (res[["d_radial_km", "d_along_km", "d_cross_km"]].abs() > RESIDUAL_MAX_KM).any(axis=1)
    print(f"Dropped {int(big.sum())} pairs with an error component above {RESIDUAL_MAX_KM} km (bad TLEs).")
    res = res[~big].reset_index(drop=True)

    res = res.round({"age_days": 5, "d_radial_km": 5, "d_along_km": 5, "d_cross_km": 5, "altitude_km": 3, "bstar": 8})
    res["epoch_old_utc"] = pd.to_datetime(res["epoch_old_utc"]).dt.strftime("%Y-%m-%dT%H:%M:%S")
    res["epoch_new_utc"] = pd.to_datetime(res["epoch_new_utc"]).dt.strftime("%Y-%m-%dT%H:%M:%S")
    res.to_csv(PROCESSED / "residuals.csv", index=False)
    print(f"Wrote data/processed/residuals.csv: {len(res)} residuals from {res['norad_id'].nunique()} objects.")

    # quick sanity check: errors should grow with TLE age
    res["age_bin"] = pd.cut(res["age_days"], [0, 1, 2, 3, 5, 7.0001])
    summary = res.groupby("age_bin", observed=True)[["d_radial_km", "d_along_km", "d_cross_km"]].agg(
        lambda x: np.median(np.abs(x))).round(3)
    summary.insert(0, "n", res.groupby("age_bin", observed=True).size())
    print("\nMedian absolute error (km) by TLE age (days):")
    print(summary.to_string())


if __name__ == "__main__":
    main()
