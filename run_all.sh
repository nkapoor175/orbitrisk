#!/usr/bin/env bash
# Runs the whole pipeline in dependency order. Run from the repo root:  bash run_all.sh
set -e

echo "== Stage 1: data (Navika) =="
python python/01_download.py
python python/02_clean.py
python python/03_residuals.py

echo "== Stage 2: learned error model (Vishi, R) =="
Rscript r/10_error_model.R

echo "== Stage 3: conjunctions and probabilities (Navika, Parthvi) =="
python python/04_screen_conjunctions.py
python python/05_exact_pc.py
python python/06_monte_carlo.py
python python/07_importance_sampling.py
python python/08_bayes_and_sensitivity.py

echo "== Stage 4: statistics and figures (Vishi, R) =="
Rscript r/11_gaussian_vs_learned_tests.R
Rscript r/12_poisson_fleet_model.R
Rscript r/13_calibration.R
Rscript r/14_plots.R

echo "Done. Figures are in figures/, tables in data/processed/."
