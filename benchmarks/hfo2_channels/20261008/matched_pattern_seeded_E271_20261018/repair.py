"""One fixed-seed numerical repair, without changing the E270 experiment."""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

REPOSITORY = Path(__file__).resolve().parents[4]
if (REPOSITORY/'vcneb').is_dir():
    sys.path.insert(0, str(REPOSITORY))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb import CalculatorModeEvaluator, ModePlane, relax_orthogonal_at_q, validate_path_geometry
from vcneb.mode_surface import _orthogonal_directions

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-matched-pattern-E271-r1')
BASE = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-matched-pattern-E270-r1')
BASE_PIN = 'f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca'
BASE_JOB = '28947338'
POLICY = dict(source_probe_step_A=.01, new_attempts_per_arm=12, maximum_new_attempts=72,
    successful_arms_reused=True, physical_inputs_unchanged=True, baseline_labels_unchanged=True,
    numerical_initializer='transported atomic/cell block-Jacobi', eigenvalue_clipping=False,
    maximum_study_allocations=2, prediction_or_TS_or_stability_certificate=False)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def prepare(case, output):
    source = Path(case).with_name('matched_pattern_release_E270_20261018')
    numeric = load_module(source/'numerical_seed.py', 'E271_inert_seed')
    matrices, metadata = [], []
    for i in range(3):
        matrix, report = numeric.seed_for(source, i, probe_index=0)
        matrices.append(matrix.tolist()); metadata.append(report)
    data = dict(policy=POLICY, base_bundle_pin=BASE_PIN, base_job=BASE_JOB,
        source_builder_sha256=sha256(source/'numerical_seed.py'),
        source_curvature_native_sha256='a182e4c386889880c77c2069e939cd65665724d98d3838e251a077204cad1c3b',
        matrices=matrices, metadata=metadata, new_DFT_calls=0)
    with Path(output).open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')
    return dict(status='prepared_without_DFT_or_submission', new_DFT_calls=0)


def baseline_plan(rows):
    """No continuation of live/failed-calculator jobs and no rerun of success."""
    if len(rows) != 3:
        raise ValueError('exactly the three original T neighborhoods required')
    plan = []
    for i, arms in enumerate(rows):
        if [r['treatment'] for r in arms] != ['A', 'J']:
            raise ValueError('all original A/J outcomes required')
        for r in arms:
            candidate = r['candidate']
            if not 0 <= r['new_attempts'] <= 12:
                raise ValueError('original attempt cap differs')
            if candidate is not None:
                if r['status'] != 'gradient_candidate_stability_unchecked':
                    raise ValueError('unknown candidate status')
            elif r['status'] not in ('attempt_budget_exhausted', 'released_gradient_target_not_reached',
                                     'not_attempted_atomic_candidate_unavailable'):
                raise ValueError('unknown original numerical failure')
        if arms[0]['candidate'] is None and (arms[1]['candidate'] is not None
                or arms[1]['status'] != 'not_attempted_atomic_candidate_unavailable'
                or arms[1]['new_attempts'] != 0):
            raise ValueError('J cannot precede an available A')
        if arms[0]['candidate'] is not None and arms[1]['status'] == 'not_attempted_atomic_candidate_unavailable':
            raise ValueError('an available A cannot justify skipping original J')
        plan.append(dict(index=i, repair_A=arms[0]['candidate'] is None,
                         repair_J=arms[1]['candidate'] is None))
    return plan


def release_space(data, atomic, reference):
    c = np.asarray(data['covectors_39_by_4']); q = np.asarray(data['orthonormal_39_by_4'])
    if (c.shape != (39, 4) or q.shape != c.shape or np.any(c[36:]) or np.any(q[36:])
            or not np.allclose(q.T@q, np.eye(4), atol=1e-12, rtol=0)
            or not np.allclose(q@(q.T@c), c, atol=1e-12, rtol=0)):
        raise ValueError('same four purely atomic physical constraints required')
    translations = np.zeros((39, 3)); translations[:36] = np.tile(np.eye(3), (12, 1))
    n = 36 if atomic else 39
    plane = ModePlane(np.asarray(reference)[:n], q[:n, :2], np.ones(n), np.eye(2),
        ('fixed physical span 1', 'fixed physical span 2'), 'angstrom', 'E270-four-pattern-controls')
    frozen = np.column_stack([q[:n, 2:], translations[:n]])
    if _orthogonal_directions(plane, frozen).shape[1] != (29 if atomic else 32):
        raise ValueError('original release rank changed')
    return plane, frozen, c, translations


def seeded_arm(study, chart, data, evaluate, matrix, arm, atomic_candidate=None):
    if arm not in ('A', 'J'):
        raise ValueError('only original A/J treatments admitted')
    atomic = arm == 'A'
    if not atomic and atomic_candidate is None:
        return dict(treatment=arm, status='not_attempted_atomic_candidate_unavailable',
                    candidate=None, new_attempts=0, origin='E271_seeded_repair')
    reference = np.zeros(39) if atomic else np.asarray(atomic_candidate['coordinates'])
    plane, frozen, c, translations = release_space(data, atomic, reference)
    def full(x):
        return np.r_[x, np.zeros(3)] if atomic else np.asarray(x)
    def guard(x):
        x = full(x)
        if np.linalg.norm(x) > study.POLICY['maximum_coordinate_norm_A']:
            return False
        try:
            validate_path_geometry([chart.reference_atoms, chart.to_atoms(x)], cell_scale=study.SCALE,
                minimum_distance=1.6, maximum_deformation=.25)
            return True
        except (ValueError, np.linalg.LinAlgError):
            return False
    evaluate.stage_start = evaluate.fresh
    try:
        point = relax_orthogonal_at_q(plane, [0., 0.],
            lambda x: study.pullback(evaluate, full(x), atomic), frozen_directions=frozen,
            gradient_tolerance=.03, optimizer='safeguarded_bfgs', trial_validator=guard,
            initial_trust_radius=.05, max_iterations=200,
            initial_hessian=matrix[:36, :36] if atomic else matrix)
        x = full(point.coordinates)
        if (np.max(np.abs(c.T@x)) > 1e-10 or np.max(np.abs(translations.T@x)) > 1e-10
                or atomic and np.any(x[36:])):
            raise ValueError('repair changed retained physical variables or fixed cell')
        candidate = dict(coordinates=x.tolist(), energy_eV_cell=point.energy,
            released_gradient_norm_eV_A=point.orthogonal_gradient_norm,
            physical_constraint_delta=(c.T@x).tolist(), stationary_TS_or_stability_certified=False)
        status = 'gradient_candidate_stability_unchecked'
    except study.AttemptCap:
        candidate, status = None, 'attempt_budget_exhausted'
    except RuntimeError as exc:
        if str(exc) != 'no starting point reached the orthogonal-gradient tolerance':
            raise
        candidate, status = None, 'released_gradient_target_not_reached'
    return dict(treatment=arm, status=status, candidate=candidate,
        new_attempts=evaluate.fresh-evaluate.stage_start, origin='E271_seeded_repair')


def replay_candidate(study, chart, candidate, directory):
    x = np.asarray(candidate['coordinates'])
    if not np.any(x):
        raise ValueError('a reused fresh A must identify its native calculation')
    matches = [p for p in directory.glob('evaluations/eval-*/result.json')
               if json.loads(p.read_bytes())['coordinates'] == candidate['coordinates']]
    if len(matches) != 1:
        raise ValueError('reused A requires exactly one complete native record')
    path = matches[0]; native = json.loads(path.read_bytes()); atoms = chart.to_atoms(x)
    atoms.calc = SinglePointCalculator(atoms, **audited_results(path.parent/'scf_000000'))
    validation = study.result_validator(0, atoms, path.parent)
    energy, gradient = chart.energy_and_gradient(x, atoms)
    if (native['status'] != 'complete' or validation != native['validation']
            or energy != candidate['energy_eV_cell'] or energy != native['enthalpy_eV']
            or not np.allclose(gradient, native['gradient'], atol=2e-12, rtol=0)):
        raise ValueError('reused native A input, E/F/stress or candidate changed')
    return energy, gradient


def validate(pin):
    bundle = ROOT/'bundle'
    if sha256(bundle/'bundle_files.json') != pin:
        raise ValueError('E271 bundle pin changed')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('immutable repair bundle changed')
    data = json.loads((bundle/'prepared.json').read_bytes())
    if data['policy'] != POLICY or data['base_bundle_pin'] != BASE_PIN or data['base_job'] != BASE_JOB:
        raise ValueError('predeclared numerical repair differs')
    study = load_module(BASE/'bundle/study.py', 'immutable_E270')
    study.validate(BASE_PIN)
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', BASE_JOB,
                                   '--format=JobIDRaw,State,ExitCode'], text=True)
    if [r.split('|')[:3] for r in raw.strip().splitlines()] != [[BASE_JOB, 'COMPLETED', '0:0']]:
        raise ValueError('original E270 must complete normally before repair admission')
    complete = json.loads((BASE/'results/complete.json').read_bytes())
    if complete['job_id'] != BASE_JOB or complete['bundle_pin'] != BASE_PIN:
        raise ValueError('terminal baseline identity differs')
    originals = [json.loads((BASE/'results'/str(i)/'treatments.json').read_bytes())['treatments'] for i in range(3)]
    plan = baseline_plan(originals)
    if not any(r['repair_A'] or r['repair_J'] for r in plan):
        raise ValueError('no unavailable baseline arm justifies new DFT')
    return study, data, originals, plan


def preflight(pin):
    study, data, originals, plan = validate(pin)
    checks = []
    for i in range(3):
        seed = BASE/'bundle/seeds'/str(i); chart = study.chart_for(seed)
        study.cached_frozen(seed, chart)
        meta = json.loads((seed/'seed.json').read_bytes()); h = np.asarray(data['matrices'][i])
        if h.shape != (39, 39) or not np.isfinite(h).all() or not np.allclose(h, h.T, atol=1e-12, rtol=0):
            raise ValueError('finite symmetric physical-coordinate numerical seed required')
        for atomic in (True, False):
            plane, frozen, _, _ = release_space(meta, atomic, np.zeros(39))
            u = _orthogonal_directions(plane, frozen)
            np.linalg.cholesky(u.T@(h[:36, :36] if atomic else h)@u)
        if originals[i][0]['candidate'] is not None:
            replay_candidate(study, chart, originals[i][0]['candidate'], BASE/'results'/str(i))
        checks.append(dict(index=i, positive_restrictions_without_clipping=True))
    result = dict(bundle_pin=pin, baseline_plan=plan, checks=checks, new_DFT_calls=0,
        original_outcome_sha256={f'{i}/{arm}': sha256(BASE/'results'/str(i)/f'treatment_{arm}.json')
                                for i in range(3) for arm in ('A', 'J')})
    study.save(ROOT/'preflight.json', result)
    return result


def run(pin):
    study, data, originals, plan = validate(pin)
    accepted = json.loads((ROOT/'accepted.json').read_bytes())
    proof = json.loads((ROOT/'preflight.json').read_bytes())
    if (proof['bundle_pin'] != pin or proof['baseline_plan'] != plan or proof['new_DFT_calls'] != 0
            or any(sha256(BASE/'results'/str(i)/f'treatment_{arm}.json') != proof['original_outcome_sha256'][f'{i}/{arm}']
                   for i in range(3) for arm in ('A', 'J'))
            or os.environ.get('RUN_DFT') != '1' or os.environ.get('SLURM_JOB_ID') != accepted['job_id']
            or accepted['bundle_pin'] != pin or os.environ.get('SLURM_JOB_NAME') != 'hfo2-E271-seeded-T'
            or os.environ.get('SLURM_JOB_PARTITION') != 'hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES') != '1' or os.environ.get('SLURM_NTASKS') != '1'
            or os.environ.get('SLURM_CPUS_PER_TASK') != '32' or len(study.queue()) > 2
            or any(os.environ.get(n) != '1' for n in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'))):
        raise ValueError('fixed baseline, accepted32CPU job and two-study-allocation cap required')
    output = ROOT/'results'; output.mkdir(exist_ok=False)
    for i in range(3):
        seed = BASE/'bundle/seeds'/str(i); chart = study.chart_for(seed)
        energy, gradient = study.cached_frozen(seed, chart)
        evaluate = CalculatorModeEvaluator(chart,
            study.make_factory(parameters={'source_directory': str(study.ASSETS)}, command=study.COMMAND),
            workdir=output/str(i)/'evaluations', calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6, result_validator=study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1', max_evaluation_attempts=24)
        cache = study.ExactCache(evaluate, energy, gradient)
        meta = json.loads((seed/'seed.json').read_bytes()); h = np.asarray(data['matrices'][i])
        rows = []
        for arm in ('A', 'J'):
            old = originals[i][0 if arm == 'A' else 1]
            if old['candidate'] is not None:
                row = dict(old, origin='E270_reused', new_attempts=0,
                           original_new_attempts=old['new_attempts'])
                if arm == 'A':
                    x = np.asarray(old['candidate']['coordinates'])
                    cache.cache[x.tobytes()] = replay_candidate(study, chart, old['candidate'], BASE/'results'/str(i))
            else:
                row = seeded_arm(study, chart, meta, cache, h, arm,
                                 None if arm == 'A' else rows[0]['candidate'])
            rows.append(row)
            (output/str(i)).mkdir(parents=True, exist_ok=True)
            study.save(output/str(i)/f'treatment_{arm}.json', row)
            print(json.dumps(dict(index=i, **row)), flush=True)
        study.save(output/str(i)/'treatments.json', dict(index=i, frozen_energy_eV_cell=energy,
            original_treatments=originals[i], treatments=rows, total_new_attempts=cache.fresh,
            original_cost_not_erased=True, old_path_labels_unchanged=True))
    study.save(output/'complete.json', dict(job_id=accepted['job_id'], bundle_pin=pin,
        status='bounded_numerical_repair_finished_not_TS_or_stability_certificate'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'preflight', 'run'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--bundle-pin')
    args = parser.parse_args()
    result = prepare(Path(__file__).resolve().parent, args.output) if args.mode == 'prepare' else (
        preflight(args.bundle_pin) if args.mode == 'preflight' else run(args.bundle_pin))
    if result is not None:
        print(json.dumps(result), flush=True)
