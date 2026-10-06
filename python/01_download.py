"""
01_download.py  --  OWNER: Navika

Purpose : Download TLE data from Space-Track in two steps:
            1. a snapshot of every payload and debris object in our orbital shell (latest TLE each),
               used later for conjunction screening;
            2. the TLE HISTORY of a random sample of those objects (N_OBJECTS), used for the error model.
Reads   : SPACETRACK_USER and SPACETRACK_PASS from the .env file (or environment variables)
Writes  : data/raw/catalog_snapshot.csv   latest TLE of every object in the shell
          data/raw/tle_history.csv        TLE history of the sampled objects
          data/raw/sampled_ids.csv        which objects were sampled
          data/raw/cache/                 raw Space-Track answers (so re-runs do not query again)
Status  : DONE (first version)

Run from the repo root:   python python/01_download.py
Options                   --refresh        ignore the cache and query Space-Track again
                          --snapshot-only  only step 1 (quick login/connection test)

Space-Track asks users not to repeat the same query and limits requests to 30 per minute and 300 per hour,
so this script pauses between requests and caches every answer.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from common import (ROOT, RAW, RANDOM_SEED, N_OBJECTS, HISTORY_DAYS, SHELL_MIN_ALT_KM, SHELL_MAX_ALT_KM,
                    SNAPSHOT_MAX_EPOCH_AGE_DAYS, DEBRIS_FRACTION, EXTRA_NORAD_IDS, IDS_PER_REQUEST,
                    SPACETRACK_DELAY_S)
from orbit_utils import mean_motion_range_for_shell

BASE = "https://www.space-track.org"
LOGIN_URL = BASE + "/ajaxauth/login"
QUERY_URL = BASE + "/basicspacedata/query/"
CACHE = RAW / "cache"

SNAPSHOT_FIELDS = "NORAD_CAT_ID,OBJECT_NAME,OBJECT_TYPE,EPOCH,BSTAR,INCLINATION,MEAN_MOTION,TLE_LINE1,TLE_LINE2"
HISTORY_FIELDS = "NORAD_CAT_ID,EPOCH,BSTAR,INCLINATION,TLE_LINE1,TLE_LINE2"
TYPE_MAP = {"PAYLOAD": "satellite", "DEBRIS": "debris"}


# --------------------------------------------------------------------------------------------
# Login and queries
# --------------------------------------------------------------------------------------------
def load_credentials():
    """Read the Space-Track login from .env (or from environment variables)."""
    values = {}
    env_file = ROOT / ".env"
    if env_file.exists():
        for line in env_file.read_text(encoding="utf-8-sig").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, val = line.split("=", 1)
            values[key.strip()] = val.strip().strip('"').strip("'")
    user = os.environ.get("SPACETRACK_USER") or values.get("SPACETRACK_USER")
    password = os.environ.get("SPACETRACK_PASS") or values.get("SPACETRACK_PASS")
    if not user or not password or user == "your_email_here":
        sys.exit("No Space-Track login found. Copy .env.example to .env and fill in SPACETRACK_USER and "
                 "SPACETRACK_PASS (your Space-Track email and password).")
    return user, password


def login(session):
    user, password = load_credentials()
    r = session.post(LOGIN_URL, data={"identity": user, "password": password}, timeout=60)
    if r.status_code != 200 or "Failed" in r.text:
        sys.exit("Space-Track login failed. Check SPACETRACK_USER / SPACETRACK_PASS in .env, and that your "
                 "account is confirmed (log in once on the website).")
    time.sleep(SPACETRACK_DELAY_S)
    print("Logged in to Space-Track.")


def run_query(session, path, refresh=False):
    """Run one Space-Track query (with caching, politeness delay and retries). Returns a list of dicts."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE / (hashlib.sha1(path.encode("utf-8")).hexdigest()[:16] + ".json")
    if cache_file.exists() and not refresh:
        return json.loads(cache_file.read_text(encoding="utf-8"))

    url = QUERY_URL + quote(path, safe="/,-_.~")
    for attempt in range(4):
        r = session.get(url, timeout=300)
        time.sleep(SPACETRACK_DELAY_S)
        if r.status_code == 200:
            try:
                data = r.json()
            except ValueError:
                sys.exit("Space-Track did not return JSON. First 200 characters of the answer:\n" + r.text[:200])
            if isinstance(data, dict):
                sys.exit(f"Space-Track returned an error: {data}")
            cache_file.write_text(json.dumps(data), encoding="utf-8")
            return data
        if r.status_code in (429, 500, 502, 503, 504):
            wait = 65 * (attempt + 1)
            print(f"Space-Track answered {r.status_code}; waiting {wait} s before retrying ...")
            time.sleep(wait)
            continue
        sys.exit(f"Space-Track error {r.status_code}: {r.text[:200]}")
    sys.exit("Gave up after 4 attempts. Wait a few minutes and run the script again (cached parts are kept).")


# --------------------------------------------------------------------------------------------
# Step 1: snapshot of the shell
# --------------------------------------------------------------------------------------------
def snapshot_query_path():
    n_lo, n_hi = mean_motion_range_for_shell(SHELL_MIN_ALT_KM, SHELL_MAX_ALT_KM)
    return (f"class/gp/DECAY_DATE/null-val/EPOCH/>now-{SNAPSHOT_MAX_EPOCH_AGE_DAYS}"
            f"/MEAN_MOTION/{n_lo:.4f}--{n_hi:.4f}/OBJECT_TYPE/PAYLOAD,DEBRIS"
            f"/orderby/NORAD_CAT_ID asc/predicates/{SNAPSHOT_FIELDS}/format/json")


def tidy_snapshot(rows):
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("The snapshot query returned no objects. Check the shell settings in python/common.py.")
    out = pd.DataFrame({
        "norad_id": df["NORAD_CAT_ID"].astype(int),
        "object_name": df["OBJECT_NAME"],
        "object_type": df["OBJECT_TYPE"].map(TYPE_MAP),
        "epoch_utc": pd.to_datetime(df["EPOCH"]),
        "bstar": pd.to_numeric(df["BSTAR"], errors="coerce"),
        "inclination_deg": pd.to_numeric(df["INCLINATION"], errors="coerce"),
        "mean_motion_rev_day": pd.to_numeric(df["MEAN_MOTION"], errors="coerce"),
        "line1": df["TLE_LINE1"],
        "line2": df["TLE_LINE2"],
    })
    out = out.dropna(subset=["object_type", "line1", "line2"])
    out = out.sort_values("epoch_utc").drop_duplicates("norad_id", keep="last").sort_values("norad_id")
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------------------------
# Step 2: sample objects and download their history
# --------------------------------------------------------------------------------------------
def choose_sample(snapshot):
    rng = np.random.default_rng(RANDOM_SEED)
    n_debris = int(round(N_OBJECTS * DEBRIS_FRACTION))
    n_sat = N_OBJECTS - n_debris
    chosen = []
    for obj_type, n_wanted in (("debris", n_debris), ("satellite", n_sat)):
        pool = snapshot.loc[snapshot["object_type"] == obj_type, "norad_id"].to_numpy()
        n_take = min(n_wanted, len(pool))
        if n_take < n_wanted:
            print(f"Only {len(pool)} {obj_type} objects available (wanted {n_wanted}); taking all of them.")
        chosen.extend(rng.choice(pool, size=n_take, replace=False).tolist())
    for extra in EXTRA_NORAD_IDS:
        if extra not in chosen:
            chosen.append(int(extra))
    return sorted(set(int(x) for x in chosen))


def history_query_path(ids):
    id_list = ",".join(str(i) for i in ids)
    return (f"class/gp_history/NORAD_CAT_ID/{id_list}/EPOCH/>now-{HISTORY_DAYS}"
            f"/orderby/EPOCH asc/predicates/{HISTORY_FIELDS}/format/json")


def tidy_history(rows, snapshot):
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("The history queries returned no TLEs. Is your account allowed to query gp_history?")
    out = pd.DataFrame({
        "norad_id": df["NORAD_CAT_ID"].astype(int),
        "epoch_utc": pd.to_datetime(df["EPOCH"]),
        "bstar": pd.to_numeric(df["BSTAR"], errors="coerce"),
        "inclination_deg": pd.to_numeric(df["INCLINATION"], errors="coerce"),
        "line1": df["TLE_LINE1"],
        "line2": df["TLE_LINE2"],
    })
    names = snapshot.set_index("norad_id")[["object_name", "object_type"]]
    out = out.join(names, on="norad_id")
    out = out[["norad_id", "object_name", "object_type", "epoch_utc", "bstar", "inclination_deg", "line1", "line2"]]
    out = out.sort_values(["norad_id", "epoch_utc"]).reset_index(drop=True)
    return out


def main():
    parser = argparse.ArgumentParser(description="Download TLE data from Space-Track.")
    parser.add_argument("--refresh", action="store_true", help="ignore the cache and query again")
    parser.add_argument("--snapshot-only", action="store_true", help="only download the shell snapshot")
    args = parser.parse_args()

    session = requests.Session()
    login(session)

    # Step 1: snapshot
    print(f"Step 1: snapshot of payloads and debris with mean altitude {SHELL_MIN_ALT_KM}-{SHELL_MAX_ALT_KM} km ...")
    snapshot = tidy_snapshot(run_query(session, snapshot_query_path(), refresh=args.refresh))
    snapshot.to_csv(RAW / "catalog_snapshot.csv", index=False)
    counts = snapshot["object_type"].value_counts().to_dict()
    print(f"  {len(snapshot)} objects in the shell: {counts}")
    print("  wrote data/raw/catalog_snapshot.csv")
    if args.snapshot_only:
        return

    # Step 2: history of a sample
    ids = choose_sample(snapshot)
    pd.DataFrame({"norad_id": ids}).to_csv(RAW / "sampled_ids.csv", index=False)
    batches = [ids[i:i + IDS_PER_REQUEST] for i in range(0, len(ids), IDS_PER_REQUEST)]
    print(f"Step 2: downloading {HISTORY_DAYS} days of TLE history for {len(ids)} objects "
          f"in {len(batches)} requests (about {len(batches) * SPACETRACK_DELAY_S / 60:.1f} min if not cached) ...")
    rows = []
    for batch in tqdm(batches, unit="request"):
        rows.extend(run_query(session, history_query_path(batch), refresh=args.refresh))
    history = tidy_history(rows, snapshot)
    history.to_csv(RAW / "tle_history.csv", index=False)
    per_obj = history.groupby("norad_id").size()
    print(f"  {len(history)} TLEs for {per_obj.shape[0]} objects "
          f"(median {int(per_obj.median())} TLEs per object, min {int(per_obj.min())})")
    print("  wrote data/raw/tle_history.csv")


if __name__ == "__main__":
    main()
