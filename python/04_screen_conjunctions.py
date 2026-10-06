"""
04_screen_conjunctions.py  --  OWNER: Navika

Purpose : Find close approaches ("conjunctions") between satellites and debris in the next SCREEN_DAYS days.

          Method (two stages, so it is fast enough for ~5000 objects):
            1. COARSE: propagate every object with SGP4 every SCREEN_STEP_S seconds. At each step a KD-tree finds
               satellite-debris pairs that are close enough to possibly meet within that step. For those pairs the
               time and distance of closest approach are estimated with straight-line relative motion.
            2. FINE: every candidate is re-propagated on a 1-second grid around the estimated time, and the exact
               closest approach is found with a parabola through the distance-squared samples.
          Candidates of the same pair that are less than 2 minutes apart are merged (same encounter).

          Only satellite (object 1) vs debris (object 2) pairs are kept. Pairs that move almost together
          (relative speed below MIN_REL_SPEED_KMS) are dropped: they are not short encounters, and the 2D
          collision-probability formula does not apply to them.

          TLE age: the screening starts at the newest TLE epoch in the snapshot (SCREEN_START_UTC=None), so the
          result is reproducible. age1_days / age2_days = (time of closest approach) - (epoch of that object's TLE).
          Conjunctions where either age exceeds MAX_PAIR_AGE_DAYS are dropped, because the error model of 03/10
          is only fitted for ages up to that value.
Reads   : data/raw/catalog_snapshot.csv
Writes  : data/processed/conjunctions.csv   (columns fixed in docs/DATA_CONTRACT.md)
Status  : DONE (first version)
Run from the repo root:  python python/04_screen_conjunctions.py
Options                  --days N   screen only N days (quick test, for example --days 1)
"""
import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sgp4.api import SatrecArray

from common import (RAW, PROCESSED, SCREEN_DAYS, SCREEN_START_UTC, SCREEN_STEP_S, SCREEN_DISTANCE_KM,
                    SCREEN_MARGIN_KM, MIN_REL_SPEED_KMS, MAX_PAIR_AGE_DAYS, HARD_BODY_RADIUS_M)
from orbit_utils import satrec_from_lines, encounter_basis, RE_KM

MAX_REL_SPEED_KMS = 16.0        # head-on at LEO speeds, used to size the coarse search radius
CHUNK_STEPS = 120               # coarse time steps propagated at once (memory)
FINE_HALF_WINDOW_S = 45         # fine grid is +/- this around the coarse estimate
MERGE_GAP_S = 120.0             # candidates of one pair closer in time than this are one encounter
JD_UNIX_EPOCH = 2440587.5


def to_jd_pair(unix_s):
    """Unix seconds (array) -> (jd, fr) pair accepted by SGP4."""
    days = np.asarray(unix_s, float) / 86400.0
    jd_whole = np.floor(days)
    return JD_UNIX_EPOCH + jd_whole, days - jd_whole


def epoch_unix(ts):
    return pd.Timestamp(ts).tz_localize(None).to_datetime64().astype("datetime64[ns]").astype("int64") / 1e9


def coarse_candidates(sa, is_sat, is_deb, t0, n_steps, step_s):
    """Stage 1. Returns arrays (i, j, t_est_unix) of candidate satellite/debris pairs (indices into sa)."""
    sat_idx = np.where(is_sat)[0]
    deb_idx = np.where(is_deb)[0]
    radius = MAX_REL_SPEED_KMS * step_s / 2.0 + SCREEN_DISTANCE_KM + SCREEN_MARGIN_KM
    keep_lin = SCREEN_DISTANCE_KM + SCREEN_MARGIN_KM
    out_i, out_j, out_t = [], [], []

    t_start = time.time()
    for c0 in range(0, n_steps, CHUNK_STEPS):
        steps = np.arange(c0, min(c0 + CHUNK_STEPS, n_steps))
        t_unix = t0 + steps * step_s
        jd, fr = to_jd_pair(t_unix)
        err, r, v = sa.sgp4(jd, fr)                   # shapes (n_obj, n_steps_chunk, 3)
        r = np.where(np.isnan(r), 1e9, r)             # failed propagation -> far away
        v = np.where(np.isnan(v), 0.0, v)
        for k, tk in enumerate(t_unix):
            ts = cKDTree(r[sat_idx, k, :])
            td = cKDTree(r[deb_idx, k, :])
            sdm = ts.sparse_distance_matrix(td, radius, output_type="ndarray")
            if sdm.size == 0:
                continue
            i = sat_idx[sdm["i"]]
            j = deb_idx[sdm["j"]]
            d = r[j, k, :] - r[i, k, :]
            w = v[j, k, :] - v[i, k, :]
            ww = (w * w).sum(axis=1)
            ok = ww > 1e-12
            tau = np.zeros(len(ww))
            tau[ok] = -(d[ok] * w[ok]).sum(axis=1) / ww[ok]
            tau = np.clip(tau, -step_s / 2.0, step_s / 2.0)    # each step owns the window around itself
            dmin = np.linalg.norm(d + w * tau[:, None], axis=1)
            sel = ok & (dmin < keep_lin)
            if sel.any():
                out_i.append(i[sel])
                out_j.append(j[sel])
                out_t.append(tk + tau[sel])
        done = min(c0 + CHUNK_STEPS, n_steps)
        if (c0 // CHUNK_STEPS) % 10 == 0 or done == n_steps:
            el = time.time() - t_start
            print(f"  coarse screening: {done}/{n_steps} steps ({el:.0f} s elapsed)")
    if not out_i:
        return np.array([], int), np.array([], int), np.array([], float)
    return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_t)


def merge_candidates(i, j, t):
    """Keep one candidate per encounter: sort by pair and time, start a new group when the gap exceeds MERGE_GAP_S."""
    order = np.lexsort((t, j, i))
    i, j, t = i[order], j[order], t[order]
    new_group = np.ones(len(t), bool)
    same_pair = (i[1:] == i[:-1]) & (j[1:] == j[:-1]) & ((t[1:] - t[:-1]) < MERGE_GAP_S)
    new_group[1:] = ~same_pair
    return i[new_group], j[new_group], t[new_group]


def refine(sats, i, j, t_est):
    """Stage 2. For each candidate find the exact time, distance and states at closest approach."""
    grid = np.arange(-FINE_HALF_WINDOW_S, FINE_HALF_WINDOW_S + 1, 1.0)
    rows = []
    for n in range(len(i)):
        a, b = sats[i[n]], sats[j[n]]
        tt = t_est[n] + grid
        jd, fr = to_jd_pair(tt)
        ea, ra, va = a.sgp4_array(jd, fr)
        eb, rb, vb = b.sgp4_array(jd, fr)
        if (ea != 0).any() or (eb != 0).any():
            continue
        d2 = ((rb - ra) ** 2).sum(axis=1)
        m = int(np.argmin(d2))
        if m == 0 or m == len(grid) - 1:
            continue                                  # minimum at the edge: not a clean encounter
        y0, y1, y2 = d2[m - 1], d2[m], d2[m + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom > 0 else 0.0
        shift = float(np.clip(shift, -1.0, 1.0))
        t_tca = tt[m] + shift
        jd1, fr1 = to_jd_pair(np.array([t_tca]))
        _, r1, v1 = a.sgp4_array(jd1, fr1)
        _, r2, v2 = b.sgp4_array(jd1, fr1)
        rows.append((i[n], j[n], t_tca, r1[0], v1[0], r2[0], v2[0]))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--days", type=float, default=SCREEN_DAYS, help="days to screen (default from common.py)")
    args = parser.parse_args()

    src = RAW / "catalog_snapshot.csv"
    if not src.exists():
        sys.exit("data/raw/catalog_snapshot.csv not found. Run python python/01_download.py first.")
    cat = pd.read_csv(src, parse_dates=["epoch_utc"]).dropna(subset=["line1", "line2", "object_type"])
    cat = cat.reset_index(drop=True)
    print(f"Loaded {len(cat)} objects: {cat['object_type'].value_counts().to_dict()}")

    sats = [satrec_from_lines(a, b) for a, b in zip(cat["line1"], cat["line2"])]
    sa = SatrecArray(sats)
    is_sat = (cat["object_type"] == "satellite").to_numpy()
    is_deb = (cat["object_type"] == "debris").to_numpy()
    epoch_s = np.array([epoch_unix(t) for t in cat["epoch_utc"]])

    t0 = epoch_unix(SCREEN_START_UTC) if SCREEN_START_UTC else float(epoch_s.max())
    n_steps = int(np.ceil(args.days * 86400.0 / SCREEN_STEP_S)) + 1
    start_txt = pd.to_datetime(t0, unit="s")
    print(f"Screening {args.days:g} days from {start_txt} UTC in {SCREEN_STEP_S} s steps ({n_steps} steps) ...")

    i, j, t_est = coarse_candidates(sa, is_sat, is_deb, t0, n_steps, SCREEN_STEP_S)
    print(f"Coarse stage: {len(i)} raw candidates.")
    if len(i) == 0:
        sys.exit("No candidates found. Check SCREEN_DISTANCE_KM or the catalogue.")
    i, j, t_est = merge_candidates(i, j, t_est)
    print(f"After merging the same encounter: {len(i)} candidates. Refining ...")

    t_refine = time.time()
    rows = refine(sats, i, j, t_est)
    print(f"Fine stage done in {time.time() - t_refine:.0f} s: {len(rows)} refined encounters.")

    out = []
    for (ii, jj, t_tca, r1, v1, r2, v2) in rows:
        miss_vec = r2 - r1
        miss = float(np.linalg.norm(miss_vec))
        w = v2 - v1
        rel_speed = float(np.linalg.norm(w))
        if miss > SCREEN_DISTANCE_KM or rel_speed < MIN_REL_SPEED_KMS:
            continue
        age1 = (t_tca - epoch_s[ii]) / 86400.0
        age2 = (t_tca - epoch_s[jj]) / 86400.0
        if age1 > MAX_PAIR_AGE_DAYS or age2 > MAX_PAIR_AGE_DAYS or age1 < 0 or age2 < 0:
            continue
        x_hat, y_hat, _ = encounter_basis(r1, v1, r2, v2)
        out.append({
            "norad_id_1": int(cat.at[ii, "norad_id"]),
            "norad_id_2": int(cat.at[jj, "norad_id"]),
            "tca_utc": pd.to_datetime(t_tca, unit="s").strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3],
            "miss_distance_km": miss,
            "rel_speed_kms": rel_speed,
            "enc_x_km": float(np.dot(miss_vec, x_hat)),
            "enc_y_km": float(np.dot(miss_vec, y_hat)),
            "age1_days": age1,
            "age2_days": age2,
            "altitude_km": float(0.5 * (np.linalg.norm(r1) + np.linalg.norm(r2)) - RE_KM),
            "radius_m": HARD_BODY_RADIUS_M["satellite"] + HARD_BODY_RADIUS_M["debris"],
            "r1_x_km": r1[0], "r1_y_km": r1[1], "r1_z_km": r1[2],
            "v1_x_kms": v1[0], "v1_y_kms": v1[1], "v1_z_kms": v1[2],
            "r2_x_km": r2[0], "r2_y_km": r2[1], "r2_z_km": r2[2],
            "v2_x_kms": v2[0], "v2_y_kms": v2[1], "v2_z_kms": v2[2],
        })
    if not out:
        sys.exit("No conjunctions left after filtering.")
    df = pd.DataFrame(out).sort_values("tca_utc").reset_index(drop=True)
    df.insert(0, "conj_id", np.arange(1, len(df) + 1))
    df = df.round({c: 5 for c in df.columns if c.endswith(("_km", "_kms", "_days", "_m"))})
    df.to_csv(PROCESSED / "conjunctions.csv", index=False)

    print(f"\nWrote data/processed/conjunctions.csv: {len(df)} conjunctions with miss distance < {SCREEN_DISTANCE_KM} km.")
    print(f"  distinct satellites: {df['norad_id_1'].nunique()}, distinct debris objects: {df['norad_id_2'].nunique()}")
    print(f"  miss distance (km): median {df['miss_distance_km'].median():.2f}, min {df['miss_distance_km'].min():.3f}")
    print(f"  relative speed (km/s): median {df['rel_speed_kms'].median():.1f}")
    print(f"  TLE age at closest approach (days): median {df[['age1_days', 'age2_days']].max(axis=1).median():.1f}")
    print(f"  closest 5 encounters:\n{df.nsmallest(5, 'miss_distance_km')[['conj_id', 'norad_id_1', 'norad_id_2', 'tca_utc', 'miss_distance_km', 'rel_speed_kms']].to_string(index=False)}")


if __name__ == "__main__":
    main()
