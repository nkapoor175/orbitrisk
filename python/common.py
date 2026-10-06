"""Shared paths and project-wide settings. Every Python script imports from here.

Change a number HERE, never inside a script, so the whole team runs with the same values.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
for _p in (RAW, PROCESSED, FIGURES):
    _p.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 42

# --- Which objects to download (01_download.py) -------------------------------------------
# We work in ONE orbital shell so that close approaches actually happen in a small data set.
SHELL_MIN_ALT_KM = 700            # mean altitude band of the shell
SHELL_MAX_ALT_KM = 900
SNAPSHOT_MAX_EPOCH_AGE_DAYS = 14  # only objects with a TLE newer than this (still in orbit)
N_OBJECTS = 300                   # objects whose TLE HISTORY is downloaded (for the error model)
HISTORY_DAYS = 120                # days of history per object
DEBRIS_FRACTION = 0.6             # share of the N_OBJECTS that are debris (rest are payloads/satellites)
EXTRA_NORAD_IDS = []              # force-include objects, e.g. [19460] (put NORAD numbers here)
IDS_PER_REQUEST = 25              # objects per Space-Track history request
SPACETRACK_DELAY_S = 3.0          # pause after every request (Space-Track limits: 30/min, 300/hour)

# --- Cleaning (02_clean.py) ---------------------------------------------------------------
LEO_MIN_ALT_KM = 300              # rows outside this altitude band are dropped
LEO_MAX_ALT_KM = 1200
MIN_TLES_PER_OBJECT = 5           # objects with fewer TLEs are dropped
MANEUVER_DA_KM = 0.5              # a jump in semi-major axis bigger than this between TLEs = manoeuvre

# --- Residuals (03_residuals.py) ----------------------------------------------------------
MIN_PAIR_AGE_DAYS = 0.1           # pair an older TLE with a newer one at least this much later ...
MAX_PAIR_AGE_DAYS = 7.0           # ... and at most this much later (the "TLE age" range we model)
MAX_PAIRS_PER_TLE = 2             # newer TLEs paired with each older TLE (spread across the age range); keeps the CSV small
RESIDUAL_MAX_KM = 100.0           # drop pairs with any error component bigger than this (bad TLE)

# --- Conjunction screening (04_screen_conjunctions.py) ------------------------------------
SCREEN_DAYS = 7                   # look-ahead window
SCREEN_START_UTC = None           # None = start at the newest TLE epoch in the snapshot (reproducible)
SCREEN_STEP_S = 30                # coarse time step of the screening
SCREEN_DISTANCE_KM = 5.0          # keep close approaches below this miss distance
SCREEN_MARGIN_KM = 5.0            # extra margin for the coarse (straight-line) estimate, refined afterwards
MIN_REL_SPEED_KMS = 0.5           # slower pairs are co-orbiting objects, not real short encounters: dropped

# --- Assumptions to justify in the report and vary in the sensitivity study ---------------
HARD_BODY_RADIUS_M = {"satellite": 2.0, "debris": 0.5}
ALERT_THRESHOLD = 1e-4            # Pc above this triggers an alert (confirm the usual value in the literature review)


def not_implemented(script, owner):
    print(f"{script}: not implemented yet (owner: {owner})")
