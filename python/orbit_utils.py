"""Small orbit helpers shared by the Python scripts (SGP4 parsing, altitude, radial/along/cross frame)."""
import numpy as np
from sgp4.api import Satrec

MU_KM3_S2 = 398600.4418   # Earth gravitational parameter
RE_KM = 6378.137          # Earth equatorial radius


def satrec_from_lines(line1, line2):
    """Parse a two-line element set into an SGP4 satellite record."""
    return Satrec.twoline2rv(line1, line2)


def sma_km(sat):
    """Semi-major axis in km from the TLE mean motion."""
    n_rad_s = sat.no_kozai / 60.0           # sgp4 stores rad/min
    return (MU_KM3_S2 / n_rad_s ** 2) ** (1.0 / 3.0)


def mean_altitude_km(sat):
    """Mean altitude above the equatorial radius (semi-major axis minus Earth radius)."""
    return sma_km(sat) - RE_KM


def mean_motion_range_for_shell(min_alt_km, max_alt_km):
    """Mean motion (rev/day) range matching a mean-altitude band: (lowest, highest)."""
    def rev_per_day(alt_km):
        a = RE_KM + alt_km
        n_rad_s = np.sqrt(MU_KM3_S2 / a ** 3)
        return n_rad_s * 86400.0 / (2.0 * np.pi)
    return rev_per_day(max_alt_km), rev_per_day(min_alt_km)


def to_rtn(d, r, v):
    """Rotate vectors d (n x 3) into the frame of positions r and velocities v (both n x 3).

    Returns three arrays: radial, along-track, cross-track components.
    """
    r_hat = r / np.linalg.norm(r, axis=1, keepdims=True)
    h = np.cross(r, v)
    n_hat = h / np.linalg.norm(h, axis=1, keepdims=True)    # cross-track (orbit normal)
    t_hat = np.cross(n_hat, r_hat)                          # along-track
    return (d * r_hat).sum(axis=1), (d * t_hat).sum(axis=1), (d * n_hat).sum(axis=1)
