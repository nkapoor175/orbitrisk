"""
make_synthetic_tles.py  --  TEST DATA ONLY (owner: Navika)

Creates a FAKE tle_history.csv so the pipeline (02_clean.py, 03_residuals.py, and the R scripts) can be tested
before Space-Track data is available. The satellites are invented; the "TLE errors" are random noise.

            *** NEVER use this data for results in the report. ***

Writes  : data/raw/tle_history.csv   (refuses to overwrite a real download unless you pass --force)
Run from the repo root:  python python/make_synthetic_tles.py [--objects 60] [--days 60] [--force]
"""
import argparse
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from sgp4.api import Satrec, WGS72
from sgp4.exporter import export_tle

from common import RAW, RANDOM_SEED, SHELL_MIN_ALT_KM, SHELL_MAX_ALT_KM
from orbit_utils import MU_KM3_S2, RE_KM

JD_1949 = 2433281.5   # sgp4init wants epoch as days since 1949-12-31 00:00 UT


def make_object(rng, norad, start, days):
    alt = rng.uniform(SHELL_MIN_ALT_KM, SHELL_MAX_ALT_KM)
    a = RE_KM + alt
    n = np.sqrt(MU_KM3_S2 / a ** 3) * 60.0                  # rad/min
    inc = np.radians(rng.uniform(80, 100))
    raan0, argp0, m0 = rng.uniform(0, 2 * np.pi, 3)
    bstar = 10 ** rng.uniform(-5, -3.3)
    ecc = rng.uniform(0.0005, 0.01)
    jd0 = 2451545.0 + (start - datetime(2000, 1, 1, 12, tzinfo=timezone.utc)).total_seconds() / 86400.0
    base = Satrec()
    base.sgp4init(WGS72, "i", norad, jd0 - JD_1949, bstar, 0.0, 0.0, ecc, argp0, inc, m0, n, raan0)
    rows = []
    t = 0.0
    while t < days:
        epoch = start + timedelta(days=t)
        jd = jd0 + t
        tmin = t * 1440.0
        # mean elements drift with SGP4's own secular rates, so the synthetic TLEs are consistent with SGP4
        m = (m0 + base.mdot * tmin + np.radians(rng.normal(0, 0.01))) % (2 * np.pi)     # random "fit error"
        argp = (argp0 + base.argpdot * tmin) % (2 * np.pi)
        raan = (raan0 + base.nodedot * tmin) % (2 * np.pi)
        n_eff = n * (1 + rng.normal(0, 4e-7))                                           # error in mean motion
        sat = Satrec()
        sat.sgp4init(WGS72, "i", norad, jd - JD_1949, bstar, 0.0, 0.0, ecc, argp, inc, m, n_eff, raan)
        l1, l2 = export_tle(sat)
        rows.append({"norad_id": norad, "epoch_utc": epoch.replace(tzinfo=None), "bstar": bstar,
                     "inclination_deg": np.degrees(inc), "line1": l1, "line2": l2})
        t += rng.uniform(0.3, 1.2)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--objects", type=int, default=60)
    parser.add_argument("--days", type=int, default=60)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    out = RAW / "tle_history.csv"
    if out.exists() and not args.force:
        sys.exit("data/raw/tle_history.csv already exists. Use --force to overwrite it with FAKE data.")

    rng = np.random.default_rng(RANDOM_SEED)
    start = datetime.now(timezone.utc) - timedelta(days=args.days)
    rows = []
    for i in range(args.objects):
        rows.extend(make_object(rng, 90000 + i, start, args.days))
    df = pd.DataFrame(rows)
    df.insert(1, "object_name", "SYNTHETIC " + df["norad_id"].astype(str))
    df.insert(2, "object_type", np.where(df["norad_id"] % 5 < 3, "debris", "satellite"))
    df.to_csv(out, index=False)
    print(f"Wrote FAKE data/raw/tle_history.csv: {len(df)} TLEs for {args.objects} synthetic objects. Test use only.")


if __name__ == "__main__":
    main()
