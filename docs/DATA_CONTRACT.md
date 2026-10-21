# Data contract

Every hand-off between teammates is a CSV in `data/processed/` with exactly these columns. Units are in the column name.
If you need to add a column, add it at the END and tell the readers. Never rename or reorder existing columns.

## residuals.csv  (written by `python/03_residuals.py`, Navika; read by `r/10`, `r/13`)

One row per pair of successive TLEs of the same object: the older TLE is propagated to the epoch of the newer one with SGP4 and compared.

| Column | Meaning |
| --- | --- |
| norad_id | NORAD catalogue number |
| object_name | Name from the catalogue |
| object_type | `satellite` or `debris` |
| epoch_old_utc | Epoch of the older TLE (ISO 8601) |
| epoch_new_utc | Epoch of the newer TLE (ISO 8601) |
| age_days | Days between the two epochs (the TLE age) |
| d_radial_km | Radial error, newer minus propagated |
| d_along_km | Along-track error |
| d_cross_km | Cross-track error |
| altitude_km | Mean altitude at the newer epoch |
| bstar | BSTAR drag term of the older TLE |

## params.csv  (written by `r/10_error_model.R`, Vishi; read by `python/05` to `python/08`)

Fitted error model: `sigma(age) = a * age_days ^ b` per direction, for both the learned Student-t model and the Gaussian baseline.

| Column | Meaning |
| --- | --- |
| component | `radial`, `along` or `cross` |
| model | `student_t` or `gaussian` |
| a | Scale at age = 1 day, km |
| b | Growth exponent |
| nu | Degrees of freedom (blank for `gaussian`) |
| aic | AIC of the fit |
| n_obs | Number of residuals used |

## conjunctions.csv  (written by `python/04_screen_conjunctions.py`, Navika; read by Parthvi and Vishi)

One row per satellite-debris close approach found by screening (miss distance below `SCREEN_DISTANCE_KM`).
Object 1 is always the satellite, object 2 the debris. Positions are in the SGP4 TEME frame.

| Column | Meaning |
| --- | --- |
| conj_id | Unique integer id |
| norad_id_1 | The satellite |
| norad_id_2 | The debris |
| tca_utc | Time of closest approach (ISO 8601) |
| miss_distance_km | Distance at closest approach |
| rel_speed_kms | Relative speed at closest approach |
| enc_x_km | Miss vector (object 2 minus object 1) along x_hat of the encounter plane |
| enc_y_km | Miss vector along y_hat of the encounter plane |
| age1_days | TLE age of the satellite at TCA (days since its TLE epoch) |
| age2_days | TLE age of the debris at TCA |
| altitude_km | Altitude at TCA |
| radius_m | Combined hard-body radius used (satellite plus debris) |
| r1_x_km, r1_y_km, r1_z_km | Satellite position at TCA |
| v1_x_kms, v1_y_kms, v1_z_kms | Satellite velocity at TCA |
| r2_x_km, r2_y_km, r2_z_km | Debris position at TCA |
| v2_x_kms, v2_y_kms, v2_z_kms | Debris velocity at TCA |

Encounter plane convention: `z_hat` = direction of the relative velocity (debris minus satellite), the plane is perpendicular to it,
`x_hat` = the part of the satellite's orbit normal that lies in the plane, `y_hat` = z_hat cross x_hat.
Do not rebuild this yourself: use `encounter_basis()` and `project_covariance_to_plane()` from `python/orbit_utils.py`.
To get the combined 2x2 covariance for a conjunction, give `project_covariance_to_plane()` the two states and the
(radial, along-track, cross-track) sigmas of each object from `params.csv` at that object's TLE age.

## pc_results.csv  (written by `python/08_bayes_and_sensitivity.py`, Parthvi; read by Vishi)

One row per conjunction, method and error model.

| Column | Meaning |
| --- | --- |
| conj_id | Matches `conjunctions.csv` |
| method | `exact`, `mc` or `is` |
| error_model | `gaussian` or `student_t` |
| pc | Collision probability estimate (posterior mean for Bayesian rows) |
| ci_low | Lower end of the 95% Jeffreys credible interval for `mc`, or Hoeffding confidence interval for `is`; blank for `exact` |
| ci_high | Upper end of the corresponding 95% interval; blank for `exact` |
| n_samples | Number of Monte Carlo samples (blank for `exact`) |
| runtime_s | Wall-clock seconds |

## catalog_snapshot.csv  (written by `python/01_download.py`, Navika, in data/raw/; read by `python/04`)

Latest TLE of every payload and debris object in the orbital shell.

| Column | Meaning |
| --- | --- |
| norad_id | NORAD catalogue number |
| object_name | Name from the catalogue |
| object_type | `satellite` or `debris` |
| epoch_utc | TLE epoch (UTC) |
| bstar | BSTAR drag term |
| inclination_deg | Inclination |
| mean_motion_rev_day | Mean motion, revolutions per day |
| line1, line2 | The two TLE lines |

## Other files

`pc_radius_sensitivity.csv` (Parthvi): `conj_id, radius_m, error_model, pc`.

`pc_mc.csv` and `pc_is.csv` use the same eight-column order as `pc_results.csv`.
The integration generator supplies Hoeffding confidence endpoints in `pc_is.csv`;
the standalone script 07 leaves endpoints blank when it does not calculate intervals.

`mc_convergence.csv` (Parthvi): `conj_id, method, error_model, n_samples, pc, ci_low, ci_high`.
One row per conjunction/method/model/cumulative sample checkpoint. Interval meanings
are the same as above; exact quadrature has no sample-count convergence rows.

See `docs/MONTE_CARLO_METHODOLOGY.md` for generation settings and limitations.

`is_log_diagnostics.csv` (Parthvi): `conj_id, error_model, n_samples, log_pc, numerical_underflow, contribution_ess`.
One row per IS event/model. `log_pc` is the natural logarithm of the internal
estimate (Gaussian log-sum-exp, not a logarithm of rounded CSV data).
`numerical_underflow` marks a positive estimate whose float64 conversion is zero.
ESS describes concentration of weighted contributions, not binomial sample size.
