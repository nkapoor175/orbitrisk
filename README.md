# OrbitRisk

Data-driven Bayesian Monte Carlo estimation of satellite-debris collision probability.
BAMAT207 Probability and Statistics course project. Track 4: Data Simulation (custom topic). Deadline: 25 Oct 2026.

**Idea in one line:** public TLE data has no position-uncertainty (covariance) information, so we learn the uncertainty from TLE history itself,
then estimate collision probability with importance-sampled Monte Carlo and report Bayesian credible intervals.

## Team and ownership

| Person | Owns | Scripts |
| --- | --- | --- |
| Navika | Data pipeline, conjunction screening, exact 2D Pc integral, this repo | `python/01` to `python/05` |
| Parthvi | Plain Monte Carlo, importance sampling, Bayesian intervals, radius sensitivity | `python/06` to `python/08` |
| Vishi | Learned error model, statistical tests, ROC, Poisson fit, calibration, plots (R) | `r/10` to `r/14` |

## Setup (once)

```bash
# 1. Python packages
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 2. R packages
Rscript install.R

# 3. Space-Track login (free account) for TLE history
cp .env.example .env               # then edit .env with your username and password
```

## How to run

Always run from the repository root (the folder that contains `data/`, `python/` and `r/`).

```bash
bash run_all.sh                    # whole pipeline, in the right order
```

Or one script at a time, in this order:

| Step | Command | Reads | Writes |
| --- | --- | --- | --- |
| 1 | `python python/01_download.py` | Space-Track | `data/raw/tle_history.csv` |
| 2 | `python python/02_clean.py` | raw TLEs | `data/processed/tle_clean.csv` |
| 3 | `python python/03_residuals.py` | clean TLEs | `data/processed/residuals.csv` |
| 4 | `Rscript r/10_error_model.R` | `residuals.csv` | `data/processed/params.csv` |
| 5 | `python python/04_screen_conjunctions.py` | clean TLEs | `data/processed/conjunctions.csv` |
| 6 | `python python/05_exact_pc.py` | conjunctions, params | `data/processed/pc_exact.csv` |
| 7 | `python python/06_monte_carlo.py` | conjunctions, params | `data/processed/pc_mc.csv` |
| 8 | `python python/07_importance_sampling.py` | conjunctions, params | `data/processed/pc_is.csv` |
| 9 | `python python/08_bayes_and_sensitivity.py` | MC results | `data/processed/pc_results.csv`, `pc_radius_sensitivity.csv` |
| 10 | `Rscript r/11_gaussian_vs_learned_tests.R` | `pc_results.csv` | tables and `figures/` |
| 11 | `Rscript r/12_poisson_fleet_model.R` | `conjunctions.csv` | tables and `figures/` |
| 12 | `Rscript r/13_calibration.R` | `residuals.csv`, `params.csv` | tables and `figures/` |
| 13 | `Rscript r/14_plots.R` | everything above | final `figures/` |

The exact column names of every CSV are fixed in [docs/DATA_CONTRACT.md](docs/DATA_CONTRACT.md). Do not change a column name without telling the people who read that file.

## Folder structure

```text
orbitrisk/
  README.md
  requirements.txt        Python packages
  install.R               R packages
  run_all.sh              runs everything in order
  .env.example            template for Space-Track login (copy to .env, never commit it)
  docs/
    DATA_CONTRACT.md      exact columns of every CSV hand-off
    CONTRIBUTING.md       git workflow for the team
  data/
    raw/                  downloaded TLEs (not committed, too big)
    processed/            hand-off CSVs between teammates (committed, small)
  python/                 01_download.py ... 08_bayes_and_sensitivity.py, common.py
  r/                      10_error_model.R ... 14_plots.R, common.R
  figures/                every figure used in the report
  report/                 IEEE paper, plagiarism report
```

## Shared settings

All project-wide numbers (number of objects, history length, altitude band, hard-body radius, alert threshold, random seed) live in
`python/common.py` and `r/common.R`. Change them there, never inside a script, so all three of us run with the same values.
Hard-body radii and the alert threshold are assumptions: justify them in the report and vary them in the sensitivity study.

## References to check before citing

Vallado, Crawford, Hujsak and Kelso, "Revisiting Spacetrack Report #3" (SGP4); Foster and Estes, NASA/JSC-25898, 1992; F. K. Chan, "Spacecraft Collision Probability", 2008.
Confirm every citation yourselves before it goes in the report.
