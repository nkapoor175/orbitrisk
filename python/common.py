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

# --- Data size (kept small so we hit the Oct 12 date) ---
N_OBJECTS = 300            # how many objects to download
HISTORY_DAYS = 120         # days of TLE history per object
LEO_MIN_ALT_KM = 300       # low-Earth-orbit band kept after cleaning
LEO_MAX_ALT_KM = 1200

# --- Conjunction screening ---
SCREEN_DAYS = 7            # look-ahead window
SCREEN_DISTANCE_KM = 10.0  # keep close approaches below this miss distance

# --- Assumptions to justify in the report and vary in the sensitivity study ---
HARD_BODY_RADIUS_M = {"satellite": 2.0, "debris": 0.5}
ALERT_THRESHOLD = 1e-4     # Pc above this triggers an alert (confirm the usual value in the literature review)


def not_implemented(script, owner):
    print(f"{script}: not implemented yet (owner: {owner})")
