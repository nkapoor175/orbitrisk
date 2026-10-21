# Probability hand-offs

Run from the repository root:

```sh
python -B python/generate_results.py --is-samples 10000 --chunk-size 1000 --seed 42 --focus-ids 183 911 1198 --checkpoints 1000 2000 5000 10000 --radii-m 1 2.5 5 --radius-samples 100000 --backup-existing
```

The generator preserves existing `pc_mc.csv` and corrected `pc_exact.csv`, reads
only fitted `params.csv`, and writes all-catalogue IS and integrated results.
It checks the corrected exact input against recomputation on the focused IDs.
There are no exact Student-t rows. Existing MC has 100,000 draws per event/model.
Its original runtimes are retained. To independently reproduce MC estimates:

```sh
python -B python/06_monte_carlo.py --n-samples 100000 --chunk-size 1000 --seed 42 --output data/processed/pc_mc_replay.csv
```

All output paths must be unused. For repeat generation, pass a fresh
`--output-dir`, or explicitly pass `--backup-existing`; the latter copies existing
generated outputs into a timestamped `backups/` directory before replacement.
Input MC/exact results are never replaced. Wall-clock runtime varies between runs;
probability estimates are reproducible with unchanged inputs, code and numerical
libraries. `result_generation.json` records seeds, input/code hashes, sample counts,
NumPy version and focused diagnostics. Original MC seed provenance is not embedded
in its CSV; replay is checked on the focused cases, not asserted for other rows.

MC estimates are Jeffreys posterior means `(hits + 0.5)/(N + 1)`, with equal-tailed
95% Beta credible intervals. Zero hits give a positive posterior mean and a positive
upper endpoint, not evidence of no risk. The current 100,000-draw catalogue has many
zero-hit cases. Conjunction 183's seed-42 Gaussian interval misses exact Pc due to
an audited sampling fluctuation; results are not edited to force agreement.

IS is the unbiased mean of likelihood-weighted collision contributions. Its 95%
Hoeffding intervals are frequentist confidence intervals, not Bayesian credible
intervals. They use analytic density/Jacobian upper bounds; observed maxima or ESS
are not substituted for bounds. Intervals may be wide even when empirical standard
errors are small. Full-catalogue IS uses 10,000 draws per conjunction/model, chosen
after a small runtime benchmark; this is not a claim of uniform relative precision.

Detailed convergence and radius sensitivity focus on conjunctions 183, 911 and
1198, for both models. Convergence uses cumulative IS callbacks and deterministic
MC prefixes at 1,000, 2,000, 5,000 and 10,000 draws, plus MC's existing 100,000-draw
endpoint. Radius sensitivity uses 100,000 IS draws with one shared maximum-radius
disk proposal for 1, 2.5 and 5 metres. Nested events guarantee nondecreasing estimates.
The radius CSV has no method column; this run is exclusively IS. No full-catalogue
radius sweep is implied. Each radius run has a separate deterministic stream.

The R plotting reader accepts the generated grouping and numeric columns. ROC
uses matched Gaussian/Student-t IS rows; its Student-t alerts are model-derived
labels, not observed collisions. The fleet R script prioritizes corrected exact
Gaussian results. Neither script consumes interval endpoints. No R scripts are changed.
Gaussian IS now accumulates log contributions using chunked log-sum-exp. Its
ESS uses log sums; scaled moments retain standard errors when squared ordinary
weights underflow. The Student-t accumulation path is unchanged. The API reports
`log_pc` (natural logarithm), `numerical_underflow` and `zero_contributions`.
A finite `log_pc` with float64 `pc == 0` means representation underflow, not zero
collision probability. No artificial floor is substituted. Gaussian nonfinite
log densities raise an error. Radius masks may legitimately exclude every draw.

The regenerated `is_log_diagnostics.csv` has 16,344 rows and columns
`conj_id,error_model,n_samples,log_pc,numerical_underflow,contribution_ess`.
Gaussian log probabilities come from stable internal log sums. Student-t log
probabilities come from the unchanged positive internal weighted mean, before
CSV rounding. This separate diagnostic preserves the four existing CSV schemas.

Hoeffding intervals retain the analytic contribution bound and receive an outward
rounding allowance of `64 * eps * (N + 1) * bound` plus the smallest subnormal.
This conservatively allows for accumulation/conversion roundoff, but does not
certify scipy density calculations or ill-conditioned covariance calculations.
No confidence claim uses zero empirical variance as proof of zero uncertainty.

## Regenerated results

All 8,172 events have exact Gaussian, Gaussian/Student-t MC and Gaussian/Student-t
IS rows: `pc_results.csv` has 40,860 rows; `pc_is.csv` has 16,344. There are no
exact Student-t rows. `mc_convergence.csv` has 54 rows, and the radius CSV has
18 rows. Input hashes match the pre-generation dataset hashes. Catalogue outputs
and pre-existing error-growth plots were copied into timestamped
`data/processed/backups/` directories before intentional replacement.

There are 475 flagged Gaussian float64 underflows and no Student-t underflows.
Their positive estimates retain finite log probabilities. Conjunction 801 now
has `pc = 4.9406564584e-324`, `log_pc = -745.0924876472`, and ESS about 7270.4022.
These are conversions/diagnostics of the same estimator, not arbitrary floors.
Finite logs range down to approximately -127435.9247; values at that scale cannot
be exported as positive float64 probabilities. The 8,172 Student-t point estimates
match the backed-up catalogue exactly. Focused replay on IDs 183, 801, 911 and
1198, using chunk size 7919 instead of 1000, matches the regenerated estimates.

Reported per-row IS runtimes sum to 36.954 seconds; preserved MC runtimes sum to
136.454 seconds. This compares different sample budgets (10,000 versus 100,000)
and runs at different times. It is not a matched-budget efficiency experiment.
The three focused Gaussian IS estimates at 10,000 draws are approximately
4.248968e-4 (183), 4.372344e-5 (911), and 4.495934e-5 (1198), all containing the
corrected exact values in their Hoeffding intervals. For 183, Gaussian IS stays
near 4.25e-4 across checkpoints; MC moves from a 0.001499 posterior mean at 1,000
draws to 0.000615 at 100,000. The latter remains the audited high realization,
not an estimator bug. Neither MC nor IS sample paths are guaranteed to converge
monotonically with sample count. Radius paths are monotonic because events are nested.

Both models' focused convergence and the catalogue probabilities are available
for comparison; three focused cases do not establish universal superiority of IS.
Student-t contribution ESS varies substantially: for example, about 64 for 183
at 10,000 proposal samples. Its Hoeffding interval is correspondingly conservative;
there is no exact Student-t probability to use as ground truth.

The existing R workflow was run without modifying R code:

```sh
.venv/bin/python -B - <<'PY'
from datetime import datetime, timezone
from pathlib import Path
import shutil
backup = Path('data/processed/backups') / ('plots-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
backup.mkdir(parents=True)
for source in [Path('figures'), Path('data/processed/roc_summary.csv'), Path('data/processed/fleet_risk')]:
    if source.is_dir():
        shutil.copytree(source, backup / source.name)
    elif source.exists():
        shutil.copy2(source, backup / source.name)
PY
Rscript r/12_poisson_fleet_model.R
Rscript r/14_plots.R
```

It refreshed `figures/final_error_growth_{radial,along,cross}.png` and generated
`figures/final_roc_curve.png` and `figures/mc_convergence.png`. Fleet tables are
in `data/processed/fleet_risk/`, and ROC summary is `data/processed/roc_summary.csv`.
Only one exact-Gaussian event exceeds 1e-4, so an inter-event gap model/plot cannot
be fitted. Fleet collision probability is about 0.002177828 under the R workflow's
assumptions. ROC reports AUC 1 with one model-derived positive among 8,172 events;
this is not evidence of perfect real-world detection. The existing convergence
plot draws 12 grouped curves and omits its legend above ten groups; use CSV keys
to identify individual traces.

Before rerunning R, preserve any existing plots and summary outputs that it will
replace. The generator's backup flag protects its own six outputs, not R outputs.

## Repeatable pipeline and portability

Activate the project's virtual environment, then run `bash run_all.sh`. The shell
runner uses `python` from the active environment; `PYTHON` can select a specific
interpreter. On Windows activate `.venv/Scripts/Activate.ps1` in PowerShell or
`source .venv/Scripts/activate` in Git Bash; the Bash runner still requires Bash
and Rscript on PATH. No POSIX-only interpreter path is hard-coded in the generator.

Stage 3 runs exact Gaussian quadrature, script 06 with `--backup-existing`, then
`generate_results.py --backup-existing`. The generator calls the IS and sensitivity
implementations, so standalone 07/08 are not also run. All hand-offs therefore
exist before the R consumers start. MC backups and generator backups are retained
in timestamped directories; no silent overwrite is enabled by default in script 06.
The upstream data stages retain their existing refresh/cache behavior.

`python -B python/test_pipeline.py` runs the complete shell orchestration twice in
an isolated workspace using actual exact/MC/IS/integration calculations on three
real conjunctions. It asserts order, hand-offs, deterministic values and byte-for-
byte backups. Unavailable data preparation and R stages are explicitly stubbed;
it does not claim to validate Space-Track access or R statistics. A live invocation
in this checkout stops at stage 1 because no Space-Track login or raw TLE archive
is available. Two live download-to-plots runs remain to be verified on a configured
team machine. No current input datasets or catalogue outputs were changed by these tests.

Never sum the MC posterior-mean Pc column for fleet risk. At N=100,000, each
zero-hit row contributes about 5e-6 from the prior, producing a misleading fleet
sum across thousands of events. Fleet analyses use exact Gaussian or unbiased IS
estimates. The existing R fleet workflow selects exact Gaussian when available;
the pipeline always supplies these rows. Do not use its fallback on an MC-only file.

## Additional million-sample checks

Each model/method was rerun in memory for IDs 183, 911 and 1198 with N=1,000,000,
chunk size 10,000, seed 42 and the existing stream conventions: MC uses
`SeedSequence([42, conj_id, model_index])`; IS appends `7`. These checks do not
replace the reproducible catalogue or its original sampling budgets.

| ID | Gaussian MC posterior mean | MC 95% credible interval | Exact Gaussian |
| --- | --- | --- | --- |
| 183 | 4.094996e-4 | [3.708004e-4, 4.500912e-4] | 4.248835e-4 |
| 911 | 5.749994e-5 | [4.360668e-5, 7.328491e-5] | 4.372667e-5 |
| 1198 | 5.049995e-5 | [3.754196e-5, 6.534934e-5] | 4.495986e-5 |

All three exact values lie in these MC intervals and in the million-sample
Gaussian IS confidence intervals. Gaussian IS estimates are 4.248841e-4,
4.372588e-5 and 4.495964e-5, respectively.

| ID | Student-t IS estimate | IS 95% Hoeffding confidence interval |
| --- | --- | --- |
| 183 | 1.278490e-4 | [8.227726e-5, 1.734207e-4] |
| 911 | 1.022119e-5 | [8.576890e-6, 1.186549e-5] |
| 1198 | 9.404635e-5 | [9.245858e-5, 9.563412e-5] |

The larger IS budget reduces the Hoeffding widths by about a factor of ten
relative to 10,000 draws; uncertainty from fitted parameters remains excluded.

Intervals condition on fixed fitted parameters, ages and encounter geometry. They
exclude parameter-fitting uncertainty, model misspecification, and numerical error.
Student-t components use their fitted scales directly, with no finite-variance
assumption. Tiny Gaussian densities may underflow to zero. Weighted point estimates
can theoretically exceed one; the generator rejects such output rather than silently
clipping it. Checkpoint and radius intervals are marginal, not simultaneous bands.
