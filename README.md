# OrbitRisk

Data-driven estimation of satellite-debris collision probability from public TLE data.
BAMAT207 Probability and Statistics course project (Prof. Pallavi Mishra), VIT Vellore. Track 4: Data Simulation (custom topic).

**Idea in one line:** public TLE data has no position-uncertainty (covariance) information, so we learn the uncertainty from TLE history itself,
then estimate collision probability exactly (Gaussian) and with Monte Carlo and importance sampling (Student-t), with Bayesian credible intervals.

**Main results** (full details in [report/OrbitRisk_IEEE_Report.pdf](report/OrbitRisk_IEEE_Report.pdf)):

- TLE position error grows about linearly with TLE age; along-track error is about 3× radial and cross-track.
- Errors are heavy-tailed: Student-t (ν = 7.1 radial, 1.0 along-track, 4.5 cross-track) beats the Gaussian on AIC and gives calibrated 95% intervals on held-out data.
- 8,172 satellite-debris conjunctions in one week; maximum Gaussian Pc 4.3×10⁻⁴; 0.22% chance of at least one collision.
- Under the Student-t model the total risk barely changes, but the risky encounters change: 38 exceed 10⁻⁵ instead of 18, and a Gaussian alert rule misses 23 of them (ROC AUC 0.995).

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

On Windows, run it from **Git Bash** with the virtual environment's Python and R on the path, for example
`PYTHON=.venv/Scripts/python.exe bash run_all.sh`. The full run takes about **13.5 hours** on a laptop, mostly the
10⁶-sample Monte Carlo checks and the R simulation study. All results are already committed in `data/processed/` and `figures/`,
so nothing has to be rerun to read them. Running it again backs up the previous results in `data/processed/backups/` (not committed).

Or one script at a time, in this order:

| Step | Command | Reads | Writes |
| --- | --- | --- | --- |
| 1 | `python python/01_download.py` | Space-Track | `data/raw/catalog_snapshot.csv`, `data/raw/tle_history.csv` |
| 2 | `python python/02_clean.py` | `tle_history.csv` | `data/processed/tle_clean.csv` |
| 3 | `python python/03_residuals.py` | `tle_clean.csv` | `data/processed/residuals.csv` |
| 4 | `Rscript r/10_error_model.R` | `residuals.csv` | `data/processed/params.csv`, `error_model_fit.csv` |
| 5 | `python python/04_screen_conjunctions.py` | `catalog_snapshot.csv` | `data/processed/conjunctions.csv` |
| 6 | `python python/05_exact_pc.py` | conjunctions, params | `data/processed/pc_exact.csv` |
| 7 | `python python/06_monte_carlo.py --backup-existing` | conjunctions, params | `data/processed/pc_mc.csv` |
| 8 | `python python/generate_results.py --is-samples 100000 --checkpoints 1000 10000 100000 --backup-existing` | conjunctions, params, `pc_exact.csv`, `pc_mc.csv` | `pc_is.csv`, `pc_results.csv`, `pc_radius_sensitivity.csv`, `mc_convergence.csv` |
| 9 | `Rscript r/11_gaussian_vs_learned_tests.R` | `residuals.csv`, `error_model_fit.csv` | `model_comparison_tests.csv`, `type_error_sim.csv`, `figures/` |
| 10 | `Rscript r/12_poisson_fleet_model.R` | conjunctions, `pc_results.csv` | `data/processed/fleet_risk/` |
| 11 | `Rscript r/13_calibration.R` | `residuals.csv` | `calibration_results.csv`, `figures/` |
| 12 | `Rscript r/14_plots.R` | everything above | `roc_summary.csv`, final `figures/` |

Step 8 runs the estimators of `07_importance_sampling.py` and `08_bayes_and_sensitivity.py` for every conjunction and combines
all results; those two scripts can also be run on their own. The exact column names of every CSV are fixed in
[docs/DATA_CONTRACT.md](docs/DATA_CONTRACT.md), and the Monte Carlo settings are described in
[docs/MONTE_CARLO_METHODOLOGY.md](docs/MONTE_CARLO_METHODOLOGY.md).

## Tests

Run from the repository root; each prints `OK` (or `All checks passed.`):

```bash
python python/test_exact_pc.py
python python/test_monte_carlo.py
python python/test_importance_sampling.py
python python/test_gaussian_log_weights.py
python python/test_bayes_and_sensitivity.py
python python/test_result_generation.py
python python/test_pipeline.py        # runs run_all.sh twice on a small copy; needs Git Bash on Windows
```

## Report

`report/OrbitRisk_IEEE_Report.docx` (and `.pdf`) is the IEEE-format paper. Its extra figures come from
`python report/make_report_figures.py`, and `report/build/build_full.py` rebuilds the document from the committed results.

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
  python/                 01_download.py ... 08_bayes_and_sensitivity.py, generate_results.py, common.py, orbit_utils.py,
                          test_*.py (tests), make_synthetic_tles.py (fake data, for testing only)
  r/                      10_error_model.R ... 14_plots.R, common.R
  figures/                figures produced by the pipeline
  report/                 IEEE paper (.docx, .pdf), its build scripts and extra figures
```

## Shared settings

All project-wide numbers (number of objects, history length, altitude band, hard-body radius, alert threshold, random seed) live in
`python/common.py` and `r/common.R`. Change them there, never inside a script, so all three of us run with the same values.
Hard-body radii and the alert threshold are assumptions: justify them in the report and vary them in the sensitivity study.

## Data source

TLE data from Space-Track.org (18th Space Defense Squadron), https://www.space-track.org, API classes `gp` and `gp_history`
(free account needed). The full reference list is in the report.
