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


# ---------------------------------------------------------------------------------------------
# Encounter-plane helpers (used by 04_screen_conjunctions.py and by the Pc scripts)
# ---------------------------------------------------------------------------------------------
def encounter_basis(r1, v1, r2, v2):
    """Right-handed unit vectors (x_hat, y_hat, z_hat) of the encounter frame, each a length-3 array.

    z_hat = direction of the relative velocity (object 2 minus object 1); the encounter plane is
    perpendicular to it. x_hat = the part of object 1's orbit normal that lies in that plane; y_hat = z x x.
    Same convention for every conjunction, so enc_x_km / enc_y_km in conjunctions.csv mean the same thing.
    """
    w = np.asarray(v2, float) - np.asarray(v1, float)
    z = w / np.linalg.norm(w)
    h1 = np.cross(r1, v1)
    n1 = h1 / np.linalg.norm(h1)
    x = n1 - np.dot(n1, z) * z
    if np.linalg.norm(x) < 1e-6:                      # degenerate: use radial direction instead
        r_hat = np.asarray(r1, float) / np.linalg.norm(r1)
        x = r_hat - np.dot(r_hat, z) * z
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return x, y, z


def rtn_matrix(r, v):
    """3x3 matrix whose columns are the radial, along-track and cross-track unit vectors of an orbit state."""
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    r_hat = r / np.linalg.norm(r)
    n_hat = np.cross(r, v)
    n_hat = n_hat / np.linalg.norm(n_hat)
    t_hat = np.cross(n_hat, r_hat)
    return np.column_stack([r_hat, t_hat, n_hat])


def project_covariance_to_plane(r1, v1, sigma1, r2, v2, sigma2):
    """Combined 2x2 position covariance (km^2) in the encounter plane.

    sigma1, sigma2 = (sigma_radial, sigma_along, sigma_cross) in km for object 1 and object 2 (independent
    errors, each diagonal in its own radial/along/cross frame). Returns the 2x2 matrix in (x_hat, y_hat)
    coordinates, so it can be used directly with enc_x_km and enc_y_km from conjunctions.csv.
    """
    x, y, _ = encounter_basis(r1, v1, r2, v2)
    p = np.vstack([x, y])                             # 2 x 3
    cov = np.zeros((3, 3))
    for r, v, s in ((r1, v1, sigma1), (r2, v2, sigma2)):
        m = rtn_matrix(r, v)
        cov += m @ np.diag(np.square(s)) @ m.T
    return p @ cov @ p.T
