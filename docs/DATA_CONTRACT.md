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

One row per close approach found by screening.

| Column | Meaning |
| --- | --- |
| conj_id | Unique integer id |
| norad_id_1 | First object (the satellite) |
| norad_id_2 | Second object (the debris) |
| tca_utc | Time of closest approach (ISO 8601) |
| miss_distance_km | Distance at closest approach |
| rel_speed_kms | Relative speed at closest approach |
| enc_x_km | Miss vector component 1 in the encounter plane |
| enc_y_km | Miss vector component 2 in the encounter plane |
| age1_days | TLE age of object 1 at TCA |
| age2_days | TLE age of object 2 at TCA |
| altitude_km | Altitude at TCA |
| radius_m | Combined hard-body radius used |

## pc_results.csv  (written by `python/08_bayes_and_sensitivity.py`, Parthvi; read by Vishi)

One row per conjunction, method and error model.

| Column | Meaning |
| --- | --- |
| conj_id | Matches `conjunctions.csv` |
| method | `exact`, `mc` or `is` |
| error_model | `gaussian` or `student_t` |
| pc | Collision probability estimate (posterior mean for Bayesian rows) |
| ci_low | Lower end of the 95% credible interval |
| ci_high | Upper end of the 95% credible interval |
| n_samples | Number of Monte Carlo samples (blank for `exact`) |
| runtime_s | Wall-clock seconds |

## Other files

`pc_radius_sensitivity.csv` (Parthvi): `conj_id, radius_m, error_model, pc`.
