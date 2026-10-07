"""
test_exact_pc.py  --  checks that the exact Pc integral in 05_exact_pc.py is correct (owner: Navika)

Run from the repo root:  python python/test_exact_pc.py
Prints PASS for every check and stops with an error if one fails.
"""
import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("exact_pc", Path(__file__).with_name("05_exact_pc.py"))
exact_pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact_pc)
pc_circle = exact_pc.pc_gaussian_circle


def cov_of(sx, sy, rho=0.0):
    return np.array([[[sx ** 2, rho * sx * sy], [rho * sx * sy, sy ** 2]]])


def check(name, ok, detail=""):
    print(("PASS  " if ok else "FAIL  ") + name + (f"   {detail}" if detail else ""))
    assert ok, name


# 1. Zero miss distance, round covariance: closed form  Pc = 1 - exp(-R^2 / (2 sigma^2))
for R, s in [(0.0025, 0.05), (0.1, 0.05), (0.2, 0.05), (0.01, 0.01)]:
    got = pc_circle([0.0], [0.0], cov_of(s, s), R)[0]
    want = 1.0 - np.exp(-R ** 2 / (2 * s ** 2))
    check(f"centred circle R={R}, sigma={s}", abs(got - want) < 1e-9 * max(1, want) + 1e-12, f"{got:.6e} vs {want:.6e}")

# 2. Small circle far from the centre: Pc ~ area x density at the miss vector
mx, my, sx, sy, rho, R = 0.08, -0.05, 0.1, 0.06, 0.3, 0.0025
c = cov_of(sx, sy, rho)
d = np.array([mx, my])
dens = np.exp(-0.5 * d @ np.linalg.inv(c[0]) @ d) / (2 * np.pi * np.sqrt(np.linalg.det(c[0])))
got = pc_circle([mx], [my], c, R)[0]
want = np.pi * R ** 2 * dens
check("small circle ~ area x density", abs(got / want - 1) < 1e-3, f"{got:.6e} vs {want:.6e}")

# 3. Correlated, off-centre, larger circle: compare with Monte Carlo (5 million draws)
rng = np.random.default_rng(1)
mx, my, sx, sy, rho, R = 0.05, 0.03, 0.04, 0.09, -0.5, 0.06
c = cov_of(sx, sy, rho)
draw = rng.multivariate_normal([mx, my], c[0], size=5_000_000)
mc = float((np.hypot(draw[:, 0], draw[:, 1]) <= R).mean())
got = pc_circle([mx], [my], c, R)[0]
se = np.sqrt(mc * (1 - mc) / 5_000_000)
check("matches Monte Carlo", abs(got - mc) < 4 * se, f"exact {got:.5f}  MC {mc:.5f} +- {se:.5f}")

# 4. Very large circle gives 1, far-away circle gives about 0, all values stay in [0, 1]
check("huge circle -> 1", abs(pc_circle([0.0], [0.0], cov_of(0.05, 0.05), 5.0)[0] - 1.0) < 1e-9)
far = pc_circle([20.0], [0.0], cov_of(0.05, 0.05), 0.0025)[0]
check("far-away circle -> 0", 0.0 <= far < 1e-12, f"{far:.2e}")

# 5. Circle bigger than sigma, off-centre (needs the finer grid that compute_pc chooses): check against MC
mx, my, sx, sy, R = 0.0, 0.1, 0.01, 0.01, 0.1
c = cov_of(sx, sy)
draw = rng.multivariate_normal([mx, my], c[0], size=2_000_000)
mc = float((np.hypot(draw[:, 0], draw[:, 1]) <= R).mean())
got = pc_circle([mx], [my], c, R, n_r=256, n_theta=1024)[0]
se = np.sqrt(mc * (1 - mc) / 2_000_000)
check("circle >> sigma, fine grid", abs(got - mc) < 4 * se, f"exact {got:.5f}  MC {mc:.5f} +- {se:.5f}")

# 6. Vectorised call gives the same answer as one-at-a-time calls
cs = np.concatenate([cov_of(0.05, 0.08, 0.2), cov_of(0.2, 0.1, -0.4), cov_of(0.03, 0.03)])
mxs, mys = np.array([0.02, 0.1, 0.0]), np.array([0.0, -0.05, 0.01])
batch = pc_circle(mxs, mys, cs, 0.0025)
single = np.array([pc_circle(mxs[i:i + 1], mys[i:i + 1], cs[i:i + 1], 0.0025)[0] for i in range(3)])
check("batch == single", np.allclose(batch, single, rtol=1e-12, atol=0))

print("\nAll checks passed.")
