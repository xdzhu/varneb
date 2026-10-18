"""Finite signed preserving controls; exact cached F, no parameter changes."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from scripts.audit_hfo2_static_replica import sha256
from vcneb import CalculatorModeEvaluator, relax_orthogonal_at_q, validate_path_geometry

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-preserving-branches-E272-r1')
BASE = ROOT.with_name('20261018-varneb-matched-pattern-E270-r1')
HELPER = ROOT.with_name('20261018-varneb-matched-pattern-E271-r1')
PINS = {BASE: 'f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca',
        HELPER: 'bcdb328562237d84522ddc5aa597b291ef5ba622716907e6692cd530cbda32f6'}


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def save(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')


def check_files(root, pin):
    bundle = root/'bundle'
    if sha256(bundle/'bundle_files.json') != pin: raise ValueError('immutable bundle pin differs')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('immutable source/input changed')


def validate(pin):
    check_files(ROOT, pin)
    for root, source_pin in PINS.items(): check_files(root, source_pin)
    study = load(BASE/'bundle/study.py', 'E272_original_native')
    helper = load(HELPER/'bundle/repair.py', 'E272_release_space')
    study.validate(PINS[BASE])
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', '28947723',
                                   '--format=JobIDRaw,State,ExitCode'], text=True).strip()
    if raw.split('|')[:3] != ['28947723', 'COMPLETED', '0:0']:
        raise ValueError('finish original T repair before P admission')
    config = json.loads((ROOT/'bundle/prepared/config.json').read_bytes())
    policy = config['policy']
    if (config['new_DFT_calls'] != 0 or config['execution_admitted']
            or policy['prepared_strains'] != [0., .01] or policy['midpoint_admitted']
            or policy['signed_atomic_starts_A'] != [-.1, .1]
            or policy['maximum_attempts_per_arm'] != 12 or policy['maximum_prospective_attempts'] != 104
            or policy['signed_direction_probe_step_A'] != .01 or policy['gradient_tolerance_eV_A'] != .03
            or policy['numerical_seed'] != 'unit branch coordinate plus on-site Jacobi complement'
            or policy['numerical_branch_curvature_eV_A2'] != 1. or policy['eigenvalue_clipping']
            or policy['retained'] != ['Qx', 'Qz', 'Qy', 'polar_proxy_z']):
        raise ValueError('finite same-pattern physical protocol differs')
    return study, helper, config


def guard(study, chart, x):
    if np.linalg.norm(x) > .5: return False
    try:
        validate_path_geometry([chart.reference_atoms, chart.to_atoms(x)], cell_scale=study.SCALE,
                              minimum_distance=1.6, maximum_deformation=.25)
        return True
    except (ValueError, np.linalg.LinAlgError):
        return False


def arm(study, helper, chart, data, cache, hessian, start, name):
    if start is None:
        return dict(treatment=name, status='not_attempted_atomic_candidate_unavailable', candidate=None, new_attempts=0)
    atomic = name == 'A'
    plane, frozen, c, translations = helper.release_space(data, atomic, start)
    def full(x): return np.r_[x, np.zeros(3)] if atomic else np.asarray(x)
    cache.stage_start = cache.fresh
    try:
        point = relax_orthogonal_at_q(plane, [0., 0.],
            lambda x: study.pullback(cache, full(x), atomic), frozen_directions=frozen,
            gradient_tolerance=.03, optimizer='safeguarded_bfgs', initial_trust_radius=.05,
            max_iterations=200, initial_hessian=hessian[:36, :36] if atomic else hessian,
            trial_validator=lambda x: guard(study, chart, full(x)))
        x = full(point.coordinates)
        if (np.max(np.abs(c.T@x)) > 1e-10 or np.max(np.abs(translations.T@x)) > 1e-10
                or atomic and np.any(x[36:])):
            raise ValueError('physical pattern, translation or A cell changed')
        candidate = dict(coordinates=x.tolist(), energy_eV_cell=point.energy,
            released_gradient_norm_eV_A=point.orthogonal_gradient_norm,
            physical_constraint_delta=(c.T@x).tolist(), stationary_TS_or_stability_certified=False)
        status = 'gradient_candidate_stability_unchecked'
    except study.AttemptCap:
        candidate, status = None, 'attempt_budget_exhausted'
    except RuntimeError as exc:
        if str(exc) != 'no starting point reached the orthogonal-gradient tolerance': raise
        candidate, status = None, 'released_gradient_target_not_reached'
    return dict(treatment=name, status=status, candidate=candidate, new_attempts=cache.fresh-cache.stage_start)


def directional_check(study, chart, cache, candidate, direction):
    if candidate is None:
        return dict(status='unavailable_A_candidate', new_attempts=0, whole_stability_certified=False)
    x = np.asarray(candidate['coordinates']); probes = [x+s*.01*direction for s in (-1, 1)]
    if not all(guard(study, chart, p) for p in probes):
        return dict(status='original_geometry_guard_prevented_checks', new_attempts=0, whole_stability_certified=False)
    cache.stage_start = cache.fresh
    values = [cache(p) for p in probes]
    curvature_g = float((values[1][1]-values[0][1])@direction/.02)
    curvature_e = float((values[1][0]-2*candidate['energy_eV_cell']+values[0][0])/.01**2)
    return dict(status='known_direction_positive_at_one_step' if min(curvature_g, curvature_e) > 0
                else 'known_direction_nonpositive_or_unresolved',
        gradient_secant_eV_A2=curvature_g, energy_secant_eV_A2=curvature_e,
        step_A=.01, new_attempts=cache.fresh-cache.stage_start, whole_stability_certified=False,
        coordinates=[p.tolist() for p in probes])


def preflight(pin):
    study, helper, config = validate(pin)
    directory = ROOT/'preflight'; directory.mkdir(exist_ok=False)
    rows = []
    for index, job, strain in ((0, '28727382', 0.), (2, '28812952', .01)):
        seed = ROOT/'bundle/prepared/seeds'/str(index)
        data = json.loads((seed/'seed.json').read_bytes()); chart = study.chart_for(seed)
        if data['source_job'] != job or data['substrate_strain'] != strain:
            raise ValueError('original P source changed')
        energy, gradient = study.cached_frozen(seed, chart)
        prep = json.loads((seed/'branch_preparation.json').read_bytes()); h = np.asarray(prep['numerical_hessian_39'])
        if h.shape != (39, 39) or not np.isfinite(h).all() or not np.allclose(h, h.T, atol=1e-12, rtol=0):
            raise ValueError('finite physical-coordinate numerical initializer required')
        np.linalg.cholesky(.5*(h+h.T))
        starts = np.asarray(prep['signed_starts'])
        if starts.shape != (2, 39) or np.any(starts[:, 36:]) or not np.allclose(starts[0], -starts[1], atol=1e-12, rtol=0):
            raise ValueError('two signed purely atomic starts required')
        for sign, x in zip(('minus', 'plus'), starts):
            plane, frozen, c, translations = helper.release_space(data, True, x)
            if (not guard(study, chart, x) or not np.isclose(np.linalg.norm(x), .1, atol=1e-12, rtol=0)
                    or np.max(np.abs(c.T@x)) > 1e-12 or np.max(np.abs(translations.T@x)) > 1e-12):
                raise ValueError('signed fixed-pattern start changed')
            folder = directory/f'{index}_{sign}'; folder.mkdir()
            study.write_structure(folder/'STRU', chart.to_atoms(x))
            study.geometry_roundtrip(folder/'STRU', chart.to_atoms(x))
        rows.append(dict(index=index, frozen_energy_eV_cell=energy, frozen_gradient_norm_eV_A=float(np.linalg.norm(gradient))))
    result = dict(bundle_pin=pin, rows=rows, new_DFT_calls=0, original_T_job='28947723', midpoint_admitted=False)
    save(ROOT/'preflight.json', result)
    return result


def run(pin):
    study, helper, config = validate(pin)
    accepted = json.loads((ROOT/'accepted.json').read_bytes()); proof = json.loads((ROOT/'preflight.json').read_bytes())
    if (accepted['bundle_pin'] != pin or proof['bundle_pin'] != pin or proof['new_DFT_calls'] != 0
            or len(proof['rows']) != 2 or len(study.queue()) > 2
            or os.environ.get('RUN_DFT') != '1' or os.environ.get('SLURM_JOB_ID') != accepted['job_id']
            or os.environ.get('SLURM_JOB_NAME') != 'hfo2-E272-preserving'
            or os.environ.get('SLURM_JOB_PARTITION') != 'hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES') != '1' or os.environ.get('SLURM_NTASKS') != '1'
            or os.environ.get('SLURM_CPUS_PER_TASK') != '32'
            or any(os.environ.get(n) != '1' for n in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'))):
        raise ValueError('accepted32CPU allocation, native preflight and two-study cap required')
    output = ROOT/'results'; output.mkdir(exist_ok=False); total = 0
    for index in (0, 2):
        seed = ROOT/'bundle/prepared/seeds'/str(index); chart = study.chart_for(seed)
        data = json.loads((seed/'seed.json').read_bytes()); prep = json.loads((seed/'branch_preparation.json').read_bytes())
        e0, g0 = study.cached_frozen(seed, chart); h = np.asarray(prep['numerical_hessian_39'])
        direction = np.asarray(prep['signed_starts'][1])/.1
        for sign, start in zip(('minus', 'plus'), prep['signed_starts']):
            folder = output/str(index)/sign; folder.mkdir(parents=True)
            evaluate = CalculatorModeEvaluator(chart,
                study.make_factory(parameters={'source_directory': str(study.ASSETS)}, command=study.COMMAND),
                workdir=folder/'evaluations', calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
                minimum_distance_A=1.6, result_validator=study.result_validator,
                validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1', max_evaluation_attempts=26)
            cache = study.ExactCache(evaluate, e0, g0)
            a = arm(study, helper, chart, data, cache, h, np.asarray(start), 'A')
            save(folder/'treatment_A.json', a)
            check = directional_check(study, chart, cache, a['candidate'], direction)
            save(folder/'known_direction_check.json', check)
            j = arm(study, helper, chart, data, cache, h,
                    None if a['candidate'] is None else np.asarray(a['candidate']['coordinates']), 'J')
            save(folder/'treatment_J.json', j)
            total += cache.fresh
            result = dict(index=index, strain=data['substrate_strain'], sign=sign,
                frozen_energy_eV_cell=e0, treatments=[a, j], known_direction_check=check,
                total_new_attempts=cache.fresh, old_path_labels_unchanged=True, stationary_or_stability_certified=False)
            if cache.fresh > 26 or total > 104: raise ValueError('registered branch budget exceeded')
            save(folder/'treatments.json', result)
            print(json.dumps(dict(index=index, sign=sign, new_attempts=cache.fresh,
                A_status=a['status'], J_status=j['status'], known_direction_status=check['status'])), flush=True)
    save(output/'complete.json', dict(job_id=accepted['job_id'], bundle_pin=pin, total_new_attempts=total,
        status='bounded_P_extreme_strain_signed_controls_finished_not_TS_or_stability_certificate'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'run')); parser.add_argument('--bundle-pin', required=True)
    args = parser.parse_args(); result = preflight(args.bundle_pin) if args.mode == 'preflight' else run(args.bundle_pin)
    if result is not None: print(json.dumps(result))
