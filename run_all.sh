#!/usr/bin/env bash
# Runs the whole pipeline in dependency order. Run from the repo root:  bash run_all.sh
set -e
# Override PYTHON for a specific interpreter, including a Windows venv under Bash.
PYTHON="${PYTHON:-python}"

echo "== Stage 1: data (Navika) =="
"$PYTHON" python/01_download.py
"$PYTHON" python/02_clean.py
"$PYTHON" python/03_residuals.py

echo "== Stage 2: learned error model (Vishi, R) =="
Rscript r/10_error_model.R

echo "== Stage 3: conjunctions and probabilities (Navika, Parthvi) =="
"$PYTHON" python/04_screen_conjunctions.py
"$PYTHON" python/05_exact_pc.py
"$PYTHON" python/06_monte_carlo.py --backup-existing
# Generates IS, radius sensitivity, integrated results and convergence before R.
# Scripts 07/08 provide its estimators; running them separately duplicates work.
"$PYTHON" python/generate_results.py --backup-existing

echo "== Stage 4: statistics and figures (Vishi, R) =="
Rscript r/11_gaussian_vs_learned_tests.R
Rscript r/12_poisson_fleet_model.R
Rscript r/13_calibration.R
Rscript r/14_plots.R

echo "Done. Figures are in figures/, tables in data/processed/."
