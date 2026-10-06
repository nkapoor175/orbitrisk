"""
02_clean.py  --  OWNER: Navika

Purpose : Clean the raw TLE history:
            - drop rows with missing TLE lines and duplicate TLEs,
            - compute semi-major axis and mean altitude from each TLE with SGP4,
            - drop rows outside the LEO altitude band,
            - drop objects with too few TLEs,
            - flag manoeuvre jumps (a sudden change in semi-major axis), so 03_residuals.py never pairs
              TLEs across a manoeuvre.
Reads   : data/raw/tle_history.csv
Writes  : data/processed/tle_clean.csv
Status  : DONE (first version)
Run from the repo root:  python python/02_clean.py
"""
import sys

import numpy as np
import pandas as pd

from common import RAW, PROCESSED, LEO_MIN_ALT_KM, LEO_MAX_ALT_KM, MIN_TLES_PER_OBJECT, MANEUVER_DA_KM
from orbit_utils import satrec_from_lines, sma_km, RE_KM


def main():
    src = RAW / "tle_history.csv"
    if not src.exists():
        sys.exit("data/raw/tle_history.csv not found. Run python python/01_download.py first.")
    df = pd.read_csv(src, parse_dates=["epoch_utc"])
    n0, objs0 = len(df), df["norad_id"].nunique()
    print(f"Loaded {n0} TLEs for {objs0} objects.")

    df = df.dropna(subset=["line1", "line2", "object_type"])
    df = df.drop_duplicates(subset=["norad_id", "epoch_utc"], keep="last")
    df = df.sort_values(["norad_id", "epoch_utc"]).reset_index(drop=True)

    # semi-major axis and altitude from the TLE itself
    sma = np.full(len(df), np.nan)
    for i, (l1, l2) in enumerate(zip(df["line1"], df["line2"])):
        try:
            sma[i] = sma_km(satrec_from_lines(l1, l2))
        except Exception:
            pass                                    # unparsable TLE -> NaN, dropped below
    df["sma_km"] = sma
    df["altitude_km"] = df["sma_km"] - RE_KM
    df = df.dropna(subset=["sma_km"])

    before = len(df)
    df = df[(df["altitude_km"] >= LEO_MIN_ALT_KM) & (df["altitude_km"] <= LEO_MAX_ALT_KM)]
    print(f"Dropped {before - len(df)} TLEs outside {LEO_MIN_ALT_KM}-{LEO_MAX_ALT_KM} km altitude.")

    counts = df.groupby("norad_id")["norad_id"].transform("size")
    before_objs = df["norad_id"].nunique()
    df = df[counts >= MIN_TLES_PER_OBJECT]
    print(f"Dropped {before_objs - df['norad_id'].nunique()} objects with fewer than {MIN_TLES_PER_OBJECT} TLEs.")

    # manoeuvre flag: semi-major axis jumps by more than MANEUVER_DA_KM compared with the previous TLE
    df = df.sort_values(["norad_id", "epoch_utc"]).reset_index(drop=True)
    d_sma = df.groupby("norad_id")["sma_km"].diff().abs()
    df["jump"] = (d_sma > MANEUVER_DA_KM).astype(int)
    print(f"Flagged {int(df['jump'].sum())} manoeuvre jumps (|change in semi-major axis| > {MANEUVER_DA_KM} km).")

    cols = ["norad_id", "object_name", "object_type", "epoch_utc", "bstar", "inclination_deg",
            "sma_km", "altitude_km", "jump", "line1", "line2"]
    df = df[cols]
    df.to_csv(PROCESSED / "tle_clean.csv", index=False)
    per_type = df.drop_duplicates("norad_id")["object_type"].value_counts().to_dict()
    print(f"Kept {len(df)} TLEs for {df['norad_id'].nunique()} objects {per_type}.")
    print("Wrote data/processed/tle_clean.csv")


if __name__ == "__main__":
    main()
