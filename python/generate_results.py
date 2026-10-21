"""Generate contracted hand-offs without changing estimators or input results.

Run from root with python -B python/generate_results.py --help, using the
activated virtual environment (Windows: .venv/Scripts/Activate.ps1).
Existing outputs require --backup-existing; timestamped backups are preserved.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal, t

from common import PROCESSED, RANDOM_SEED
from orbit_utils import encounter_basis, rtn_matrix, project_covariance_to_plane

spec = importlib.util.spec_from_file_location('integration_bayes', Path(__file__).with_name('08_bayes_and_sensitivity.py'))
bayes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bayes)
mc, importance = bayes.mc, bayes.importance
CONVERGENCE_COLUMNS = ('conj_id', 'method', 'error_model', 'n_samples', 'pc', 'ci_low', 'ci_high')
LOG_COLUMNS = ('conj_id', 'error_model', 'n_samples', 'log_pc', 'numerical_underflow', 'contribution_ess')
OUTPUTS = ('pc_is.csv', 'pc_results.csv', 'pc_radius_sensitivity.csv', 'mc_convergence.csv', 'is_log_diagnostics.csv', 'result_generation.json')


def contribution_bound(row, parameters, model):
    """Analytic upper bound for script 07's single-radius contributions."""
    r1, v1, r2, v2, mean, s1, s2 = bayes._geometry(row, parameters)
    radius = float(row['radius_m']) / 1000
    if model == 'gaussian':
        cov = project_covariance_to_plane(r1, v1, s1, r2, v2, s2)
        peak = multivariate_normal(mean=mean, cov=cov).logpdf(mean)
    else:
        x, y, _ = encounter_basis(r1, v1, r2, v2)
        plane = np.vstack([x, y])
        matrix = np.column_stack([-plane @ rtn_matrix(r1, v1), plane @ rtn_matrix(r2, v2)])
        pivots, _, jacobian = importance.select_pivots(matrix)
        scales = np.concatenate([s1, s2])
        nu = np.tile(parameters[2], 2)
        peak = (t.logpdf(0, nu[pivots]) - np.log(scales[pivots])).sum() - np.log(jacobian)
    bound = float(np.exp(np.log(np.pi) + 2 * np.log(radius) + peak))
    if not np.isfinite(bound) or bound <= 0:
        raise FloatingPointError('Invalid contribution bound')
    return bound


def validate_results(frame, ids, combinations):
    if tuple(frame.columns) != mc.RESULT_COLUMNS:
        raise ValueError('Incorrect result schema/order')
    if frame.duplicated(['conj_id', 'method', 'error_model']).any():
        raise ValueError('Duplicate result keys')
    if set(map(tuple, frame[['method', 'error_model']].drop_duplicates().to_numpy())) != set(combinations):
        raise ValueError('Unexpected method/model coverage')
    for method, model in combinations:
        selected = frame[(frame.method == method) & (frame.error_model == model)]
        if set(selected.conj_id) != set(ids):
            raise ValueError('Incomplete conjunction coverage')
    if not np.isfinite(frame[['pc', 'runtime_s']]).all().all() or not frame.pc.between(0, 1).all() or (frame.runtime_s < 0).any():
        raise ValueError('Invalid probabilities/runtimes')
    exact = frame.method == 'exact'
    if frame.loc[exact, ['ci_low', 'ci_high', 'n_samples']].notna().any().any():
        raise ValueError('Exact interval/sample fields must be blank')
    simulation = frame.loc[~exact]
    if not (np.isfinite(simulation[['ci_low', 'ci_high', 'n_samples']]).all().all()
            and (simulation.n_samples > 0).all()
            and (simulation.n_samples == simulation.n_samples.round()).all()
            and simulation.ci_low.between(0, 1).all()
            and simulation.ci_high.between(0, 1).all()
            and (simulation.ci_low <= simulation.ci_high).all()):
        raise ValueError('Invalid simulation counts/intervals')


def validate_mc(frame):
    n = frame.n_samples.to_numpy()
    counts = frame.pc.to_numpy() * (n + 1) - .5
    hits = np.rint(counts).astype(int)
    if not np.allclose(counts, hits, atol=1e-6, rtol=0) or (hits < 0).any() or (hits > n).any():
        raise ValueError('MC probabilities are not Jeffreys posterior means')
    # Recalculate once per distinct count/sample pair, not per catalogue row.
    for (h, samples), indexes in frame.assign(hits=hits).groupby(['hits', 'n_samples']).groups.items():
        mean, low, high = mc.jeffreys_summary(int(h), int(samples))
        np.testing.assert_allclose(frame.loc[indexes, ['pc', 'ci_low', 'ci_high']],
                                   np.tile([mean, low, high], (len(indexes), 1)), rtol=1e-8, atol=1e-14)


def save_outputs(directory, tables, metadata, backup_existing=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / name for name in OUTPUTS]
    existing = [p for p in paths if p.exists()]
    if existing and not backup_existing:
        raise FileExistsError('Outputs exist; use a new directory or --backup-existing')
    if existing:
        backup = directory / 'backups' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        backup.mkdir(parents=True)
        for path in existing:
            shutil.copy2(path, backup / path.name)
        for path in existing:
            path.unlink()
    for name, table in tables.items():
        with (directory / name).open('x') as stream:
            table.to_csv(stream, index=False, float_format='%.10e')
    with (directory / 'result_generation.json').open('x') as stream:
        json.dump(metadata, stream, indent=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--is-samples', type=int, default=10_000)
    parser.add_argument('--chunk-size', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=RANDOM_SEED)
    parser.add_argument('--focus-ids', type=int, nargs='+', default=[183, 911, 1198])
    parser.add_argument('--checkpoints', type=int, nargs='+', default=[1000, 2000, 5000, 10000])
    parser.add_argument('--radii-m', type=float, nargs='+', default=[1., 2.5, 5.])
    parser.add_argument('--radius-samples', type=int, default=100_000)
    parser.add_argument('--output-dir', type=Path, default=PROCESSED)
    parser.add_argument('--backup-existing', action='store_true')
    args = parser.parse_args()
    if args.is_samples < 2 or args.radius_samples < 2 or args.chunk_size <= 0 or args.seed < 0:
        parser.error('Invalid sample/chunk count or seed')
    checkpoints = sorted(set(args.checkpoints + [args.is_samples]))
    if min(checkpoints) < 2 or max(checkpoints) > args.is_samples or any(n % args.chunk_size for n in checkpoints if n != args.is_samples):
        parser.error('Checkpoints must align with chunks and not exceed IS samples')
    if not np.isfinite(args.radii_m).all() or min(args.radii_m) <= 0:
        parser.error('Radii must be finite and positive')
    if not args.backup_existing and any((args.output_dir / n).exists() for n in OUTPUTS):
        parser.error('Outputs exist; choose a new output directory or --backup-existing')
    inputs = {name: PROCESSED / name for name in ('params.csv', 'conjunctions.csv', 'pc_exact.csv', 'pc_mc.csv')}
    hashes = {n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in inputs.items()}
    conjunctions = pd.read_csv(inputs['conjunctions.csv'])
    exact, saved_mc = pd.read_csv(inputs['pc_exact.csv']), pd.read_csv(inputs['pc_mc.csv'])
    ids = set(conjunctions.conj_id)
    if conjunctions.conj_id.duplicated().any() or not set(args.focus_ids) <= ids:
        raise ValueError('Invalid conjunction IDs')
    validate_results(exact, ids, [('exact', 'gaussian')])
    validate_results(saved_mc, ids, [('mc', m) for m in mc.MODELS])
    validate_mc(saved_mc)
    parameters = mc.load_parameters()
    # Guard the known reference cases against stale exact input, without rewriting it.
    exact_module = bayes._load('integration_exact', '05_exact_pc.py')
    sig = exact_module.load_sigma_model(inputs['params.csv'])
    selected = conjunctions[conjunctions.conj_id.isin(args.focus_ids)]
    recomputed, _ = exact_module.compute_pc(selected, sig)
    np.testing.assert_allclose(recomputed, exact.set_index('conj_id').loc[selected.conj_id, 'pc'], rtol=1e-6, atol=0)
    rows, convergence, sensitivity_frames = [], [], []
    diagnostics_records, log_records = [], []
    for index, row in enumerate(conjunctions.to_dict('records')):
        cid = int(row['conj_id'])
        for model_index, model in enumerate(mc.MODELS):
            started = time.perf_counter()
            bound = contribution_bound(row, parameters[model], model)
            def checkpoint(d):
                if cid in args.focus_ids and d['n_samples'] in checkpoints:
                    low, high = bayes.weighted_interval(d, bound)
                    convergence.append((cid, 'is', model, d['n_samples'], d['pc'], low, high))
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, cid, model_index, 7]))
            diagnostic = importance.estimate_pc(row, parameters[model], model, args.is_samples, args.chunk_size, rng, checkpoint)
            # Gaussian log_pc comes from log-sum-exp, before float conversion.
            # Student-t retains its unchanged, positive internal weighted mean.
            log_pc = diagnostic.get('log_pc', float(np.log(diagnostic['pc'])) if diagnostic['pc'] > 0 else -np.inf)
            log_records.append((cid, model, args.is_samples, log_pc,
                                diagnostic.get('numerical_underflow', False), diagnostic['contribution_ess']))
            low, high = bayes.weighted_interval(diagnostic, bound)
            elapsed = time.perf_counter() - started
            rows.append((cid, 'is', model, diagnostic['pc'], low, high, args.is_samples, elapsed))
            if cid in args.focus_ids:
                diagnostics_records.append(dict(conj_id=cid, model=model, bound=bound, **diagnostic))
                original = saved_mc[(saved_mc.conj_id == cid) & (saved_mc.error_model == model)].iloc[0]
                final_n = int(original.n_samples)
                for n in sorted(set(checkpoints + [final_n])):
                    rng_mc = np.random.default_rng(np.random.SeedSequence([args.seed, cid, model_index]))
                    h = mc.count_collisions(row, parameters[model], model, n, args.chunk_size, rng_mc)
                    mean, lo, hi = mc.jeffreys_summary(h, n)
                    convergence.append((cid, 'mc', model, n, mean, lo, hi))
                    if n == final_n:
                        np.testing.assert_allclose([mean, lo, hi], original[['pc', 'ci_low', 'ci_high']].to_numpy(float), rtol=1e-8)
                radii = sorted(set(args.radii_m + [2.5, float(row['radius_m'])]))
                frame, _ = bayes.sensitivity(row, parameters[model], model, radii, args.radius_samples, args.chunk_size, args.seed, 'is')
                sensitivity_frames.append(frame)
        if (index + 1) % 500 == 0:
            print(f'IS completed {index + 1}/{len(conjunctions)} conjunctions', flush=True)
    is_results = pd.DataFrame(rows, columns=mc.RESULT_COLUMNS)
    validate_results(is_results, ids, [('is', m) for m in mc.MODELS])
    combined = pd.concat([exact, saved_mc, is_results], ignore_index=True).sort_values(['conj_id', 'method', 'error_model'])
    validate_results(combined, ids, [('exact', 'gaussian')] + [(method, m) for method in ('mc', 'is') for m in mc.MODELS])
    tables = {'pc_is.csv': is_results, 'pc_results.csv': combined,
              'is_log_diagnostics.csv': pd.DataFrame(log_records, columns=LOG_COLUMNS),
              'pc_radius_sensitivity.csv': pd.concat(sensitivity_frames, ignore_index=True),
              'mc_convergence.csv': pd.DataFrame(convergence, columns=CONVERGENCE_COLUMNS).sort_values(['conj_id', 'method', 'error_model', 'n_samples'])}
    if any(hashlib.sha256(p.read_bytes()).hexdigest() != hashes[n] for n, p in inputs.items()):
        raise RuntimeError('Inputs changed during generation')
    metadata = dict(seed=args.seed, is_samples=args.is_samples, radius_samples=args.radius_samples,
                    chunk_size=args.chunk_size, focus_ids=args.focus_ids, checkpoints=checkpoints,
                    radii_m=args.radii_m, input_sha256=hashes, numpy_version=np.__version__,
                    mc_source='Existing pc_mc.csv; seed replay verified on focused IDs only; original runtime preserved',
                    script_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')},
                    focused_is_diagnostics=diagnostics_records,
                    rows={n: len(d) for n, d in tables.items()})
    save_outputs(args.output_dir, tables, metadata, args.backup_existing)
    print(metadata['rows'], flush=True)


if __name__ == '__main__':
    main()
