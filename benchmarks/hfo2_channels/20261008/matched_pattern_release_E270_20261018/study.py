"""Bounded same-pattern F/A/J interventions, using the immutable E268 runtime.

Only T's three existing sampled peaks are admitted in this first batch. This
does not certify TSs, replace path labels, fit a model or alter DFT inputs.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile

REPOSITORY = Path(__file__).resolve().parents[4]
if (REPOSITORY/'vcneb').is_dir():
    sys.path.insert(0, str(REPOSITORY))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from examples.hfo2_fixed_input_factory import CONTRACT, make_factory, read_fixed_hfo2_stru
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.hfo2_G3_local_curvature import result_validator, ASSETS, COMMAND
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb import CalculatorModeEvaluator, ModePlane, clamped_plane_vcneb_boundary, relax_orthogonal_at_q, validate_path_geometry
from vcneb.clamped_mode_coordinates import ClampedModeCoordinates

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-matched-pattern-E270-r1')
OLD = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
OLD_PIN = '689a26fed0f47dede6d7c2f2230d3181fe5a16bbc5eee0356da6f033d89987cf'
BASE = Path('benchmarks/hfo2_channels/20261008')
REGISTER = BASE/'midpoint_reverse_E269_20261018/six_neighborhoods.json'
SCALE = 5.12968067458423
POLICY = dict(channel='T', strains=[0., .005, .01], retained=['Qx', 'Qz', 'Qy', 'polar_proxy_z'],
    new_attempts_per_treatment=12, maximum_new_attempts=72, gradient_tolerance_eV_A=.03,
    initial_trust_radius_A=.05, maximum_coordinate_norm_A=.5, minimum_distance_A=1.6,
    maximum_deformation=.25, optimizer='safeguarded_bfgs', maximum_iterations=200,
    atomic_release_rank=29, joint_release_rank=32, ordinary_NEB_force_target_eV_A=.10,
    pressure_GPa=0., climb=False, maximum_allocations=2, nodes=1, MPI_ranks=32,
    walltime='06:00:00', source_E268_labels_unchanged=True,
    same_coordinate_definitions_not_equal_values_across_strains=True,
    independent_prediction_or_TS_or_stable_surface_certified=False)


def save(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def chart_for(seed):
    atoms = read_fixed_hfo2_stru(seed/'raw/STRU')
    boundary = clamped_plane_vcneb_boundary(12, atoms.cell.array, allow_tilt=True)
    return ClampedModeCoordinates(atoms, boundary=boundary, cell_scale_A=SCALE)


def cached_frozen(seed, chart):
    data = json.loads((seed/'seed.json').read_bytes())
    for name, digest in data['native_files_sha256'].items():
        if sha256(seed/'raw'/name) != digest:
            raise ValueError('native frozen geometry/input/log changed')
    if any(data['physical_input_sha256'][k] != h for k, h in CONTRACT.items()):
        raise ValueError('frozen source must use the original six physical inputs')
    raw = audited_results(seed/'raw')
    if any(not np.allclose(raw[k], data['results'][k], atol=1e-12, rtol=0) for k in raw):
        raise ValueError('frozen native E/F/stress replay differs')
    atoms = chart.reference_atoms
    atoms.calc = SinglePointCalculator(atoms, **raw)
    return chart.energy_and_gradient(np.zeros(39), atoms)


def prepare(repository, output):
    from ase.io import read
    from scripts.plot_hfo2_mechanism import rotated_t_triplet
    repository, output = Path(repository), Path(output)
    register = json.loads((repository/REGISTER).read_bytes())
    selected = [r for r in register['rows'] if r['channel'] == 'T']
    if [r['substrate_strain'] for r in selected] != POLICY['strains']:
        raise ValueError('all three original T neighborhoods required')
    output.mkdir(parents=True, exist_ok=False)
    t_path = repository/BASE/'reference_variants/T.vasp'
    orientation_path = repository/BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    t = read(t_path); _, patterns, _ = rotated_t_triplet(t)
    orientation = json.loads(orientation_path.read_bytes())['orientation']
    rotation = np.asarray(orientation['cartesian_rotation_rows'])
    if not np.array_equal(rotation, np.eye(3)[[2, 0, 1]]):
        raise ValueError('only the original fixed G2 frame is admitted')
    records = []
    for index, row in enumerate(selected):
        observation_path = repository/row['observation']
        if sha256(observation_path) != row['observation_sha256']:
            raise ValueError('registered complete sampled peak changed')
        observation = json.loads(observation_path.read_bytes())
        peak = row['sampled_peak_image']; raw_record = observation['raw_image_evaluations'][peak]
        if (peak != int(np.argmax(observation['relative_enthalpy_meV_fu']))
                or observation['climb'] or observation['pressure_GPa'] != 0):
            raise ValueError('original ordinary zero-pressure sampled peak required')
        seed = output/'seeds'/str(index); (seed/'raw/OUT.ABACUS').mkdir(parents=True)
        source = observation_path.parent/'raw'/f'image_{peak:04d}'
        hashes = {n: raw_record['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU')}
        hashes['OUT.ABACUS/running_scf.log'] = raw_record['raw_log_sha256']
        for name, digest in hashes.items():
            if sha256(source/name) != digest:
                raise ValueError('original native cache bytes differ')
            shutil.copyfile(source/name, seed/'raw'/name)
        atoms = read_fixed_hfo2_stru(seed/'raw/STRU')
        transform = np.linalg.inv(atoms.cell.array)@rotation@t.cell.array
        covectors = np.zeros((39, 4))
        for j, mode in enumerate((0, 2, 1)):
            covectors[:36, j] = (patterns[mode]@transform.T).ravel()
        weights = np.r_[np.full(4, -.25), np.full(8, .125)]
        covectors[:36, 3] = (weights[:, None]*transform[:, 2]).ravel()
        orthonormal = np.linalg.qr(covectors[:36], mode='reduced')[0]
        orthonormal = np.vstack([orthonormal, np.zeros((3, 4))])
        data = dict(index=index, channel='T', substrate_strain=row['substrate_strain'],
            source_job=row['source_job'], observation=row['observation'],
            observation_sha256=row['observation_sha256'], native_files_sha256=hashes,
            physical_input_sha256=raw_record['input_sha256'],
            results=dict(energy=raw_record['energy_eV_cell'], forces=raw_record['forces_eV_A'],
                         stress=raw_record['stress_ASE_voigt_eV_A3']),
            sampled_peak_meV_fu=row['sampled_peak_meV_fu'],
            covectors_39_by_4=covectors.tolist(), orthonormal_39_by_4=orthonormal.tolist())
        save(seed/'seed.json', data)
        chart = chart_for(seed); e, g = cached_frozen(seed, chart)
        records.append(dict(index=index, source_job=row['source_job'], substrate_strain=row['substrate_strain'],
                            frozen_energy_eV_cell=e, frozen_gradient_norm_eV_A=float(np.linalg.norm(g))))
    config = dict(policy=POLICY, physical_input_sha256=CONTRACT, runtime_root=OLD.as_posix(),
        runtime_pin=OLD_PIN, source_registration_sha256=sha256(repository/REGISTER),
        T_reference_sha256=sha256(t_path), orientation_sha256=sha256(orientation_path), frozen=records,
        new_DFT_calls=0, preserving_admitted=False)
    save(output/'config.json', config)
    return config


class AttemptCap(RuntimeError):
    pass


class ExactCache:
    """No geometric tolerance reuse: known coordinates and fresh calls share one chart."""
    def __init__(self, delegate, energy, gradient):
        self.delegate = delegate
        self.cache = {np.zeros(39).tobytes(): (float(energy), np.array(gradient, copy=True))}
        self.fresh = 0; self.hits = 0; self.stage_start = 0

    def __call__(self, x):
        x = np.asarray(x, dtype=float); key = x.tobytes()
        if key not in self.cache:
            if self.fresh-self.stage_start >= POLICY['new_attempts_per_treatment']:
                raise AttemptCap('registered per-treatment attempt cap reached; no calculator call')
            self.fresh += 1  # Failed attempts also consume the budget.
            energy, gradient = self.delegate(x)
            self.cache[key] = (float(energy), np.array(gradient, copy=True))
        else:
            self.hits += 1
        energy, gradient = self.cache[key]
        return energy, gradient.copy()


def treatments(chart, seed_data, evaluate, *, record=None):
    """Same physical pattern values within each F/A/J triplet, no new optimizer."""
    c = np.asarray(seed_data['covectors_39_by_4'])
    q = np.asarray(seed_data['orthonormal_39_by_4'])
    if (c.shape != (39, 4) or q.shape != c.shape or not np.isfinite(c).all()
            or np.linalg.matrix_rank(c) != 4
            or not np.allclose(q.T@q, np.eye(4), atol=1e-12, rtol=0)
            or not np.allclose(q@(q.T@c), c, atol=1e-12, rtol=0)
            or np.any(c[36:]) or np.any(q[36:])):
        raise ValueError('four finite purely atomic physical constraints and their frozen span required')
    translation = np.zeros((39, 3)); translation[:36] = np.tile(np.eye(3), (12, 1))
    if not np.allclose(translation.T@c, 0, atol=1e-12, rtol=0):
        raise ValueError('physical constraints must exclude rigid translation')
    def guard(full):
        if np.linalg.norm(full) > POLICY['maximum_coordinate_norm_A']:
            return False
        try:
            validate_path_geometry([chart.reference_atoms, chart.to_atoms(full)], cell_scale=SCALE,
                minimum_distance=POLICY['minimum_distance_A'], maximum_deformation=POLICY['maximum_deformation'])
            return True
        except (ValueError, np.linalg.LinAlgError):
            return False
    rows = []
    for treatment in ('A', 'J'):
        if treatment == 'J' and rows[-1]['candidate'] is None:
            rows.append(dict(treatment='J', status='not_attempted_atomic_candidate_unavailable', candidate=None,
                             new_attempts=0))
            if record is not None: record(rows[-1])
            break
        atomic = treatment == 'A'
        reference = np.zeros(36) if atomic else np.asarray(rows[-1]['candidate']['coordinates'])
        axes = q[:36, :2] if atomic else q[:, :2]
        frozen = np.column_stack([q[:36, 2:], translation[:36]]) if atomic else np.column_stack([q[:, 2:], translation])
        plane = ModePlane(reference, axes, np.ones(len(reference)), np.eye(2),
            ('fixed physical-constraint span 1', 'fixed physical-constraint span 2'), 'angstrom', 'E270-four-pattern-controls')
        def full(x):
            return np.r_[x, np.zeros(3)] if atomic else np.asarray(x)
        evaluate.stage_start = evaluate.fresh
        try:
            point = relax_orthogonal_at_q(plane, [0., 0.], lambda x: pullback(evaluate, full(x), atomic),
                frozen_directions=frozen, gradient_tolerance=POLICY['gradient_tolerance_eV_A'],
                optimizer='safeguarded_bfgs', trial_validator=lambda x: guard(full(x)),
                initial_trust_radius=POLICY['initial_trust_radius_A'], max_iterations=POLICY['maximum_iterations'])
            x = full(point.coordinates)
            if (not np.allclose(c.T@x, 0, atol=1e-10, rtol=0)
                    or not np.allclose(translation.T@x, 0, atol=1e-10, rtol=0)
                    or atomic and np.any(x[36:])):
                raise ValueError('accepted candidate changed retained physical coordinates or fixed cell')
            candidate = dict(coordinates=x.tolist(), energy_eV_cell=point.energy,
                released_gradient_norm_eV_A=point.orthogonal_gradient_norm,
                physical_constraint_delta=(c.T@x).tolist(), stationary_TS_or_stability_certified=False)
            status = 'gradient_candidate_stability_unchecked'
        except AttemptCap:
            candidate, status = None, 'attempt_budget_exhausted'
        except RuntimeError as exc:
            if str(exc) != 'no starting point reached the orthogonal-gradient tolerance':
                raise
            candidate, status = None, 'released_gradient_target_not_reached'
        rows.append(dict(treatment=treatment, status=status, candidate=candidate,
                         new_attempts=evaluate.fresh-evaluate.stage_start))
        if record is not None: record(rows[-1])
    return rows


def pullback(evaluate, x, atomic):
    energy, gradient = evaluate(x)
    return energy, gradient[:36] if atomic else gradient


def validate(pin):
    bundle = ROOT/'bundle'
    if sha256(bundle/'bundle_files.json') != pin:
        raise ValueError('external E270 bundle pin changed')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('E270 immutable input drift')
    config = json.loads((bundle/'config.json').read_bytes())
    if config['policy'] != POLICY or config['physical_input_sha256'] != CONTRACT or config['preserving_admitted']:
        raise ValueError('finite T-only intervention or physical settings changed')
    if config['runtime_root'] != OLD.as_posix() or config['runtime_pin'] != OLD_PIN:
        raise ValueError('fixed native runtime differs')
    files = json.loads((OLD/'bundle/bundle_files.json').read_bytes())
    if sha256(OLD/'bundle/bundle_files.json') != OLD_PIN or sha256(OLD/'bundle/endpoint_worker.py') != files['endpoint_worker.py']:
        raise ValueError('original runtime checker drift')
    spec = importlib.util.spec_from_file_location('E268_immutable_checker', OLD/'bundle/endpoint_worker.py')
    checker = importlib.util.module_from_spec(spec); spec.loader.exec_module(checker)
    checker.validate_bundle(OLD, OLD_PIN)
    return config


def queue():
    raw = subprocess.check_output(['squeue', '-u', 'iai806', '-h', '-o', '%i|%j|%T|%P|%C'], text=True)
    return [s for s in raw.splitlines() if s.split('|')[1].lower().startswith(('hfo2', 'varneb', 'vcneb'))]


def preflight(pin):
    validate(pin)
    output = ROOT/'preflight'; output.mkdir(exist_ok=False)
    rows = []
    for i in range(3):
        seed = ROOT/'bundle/seeds'/str(i); chart = chart_for(seed)
        e, g = cached_frozen(seed, chart)
        directory = output/str(i); directory.mkdir()
        for name in ('INPUT', 'KPT'):
            shutil.copyfile(ASSETS/name, directory/name)
        write_structure(directory/'STRU', chart.reference_atoms)
        geometry_roundtrip(directory/'STRU', chart.reference_atoms)
        rows.append(dict(index=i, cached_energy_eV_cell=e, cached_gradient_norm_eV_A=float(np.linalg.norm(g))))
    result = dict(status='three_T_native_seeds_ready_no_DFT', bundle_pin=pin, rows=rows, new_DFT_calls=0)
    save(ROOT/'preflight.json', result)
    return result


def run(pin):
    validate(pin)
    accepted = json.loads((ROOT/'accepted.json').read_bytes())
    if (os.environ.get('RUN_DFT') != '1' or os.environ.get('SLURM_JOB_ID') != accepted['job_id']
            or os.environ.get('SLURM_JOB_NAME') != 'hfo2-E270-pattern-T'
            or os.environ.get('SLURM_JOB_PARTITION') != 'hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES') != '1' or os.environ.get('SLURM_NTASKS') != '1'
            or os.environ.get('SLURM_CPUS_PER_TASK') != '32'
            or any(os.environ.get(n) != '1' for n in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'))
            or len(queue()) > 2):
        raise ValueError('original accepted single32CPU allocation and two-study-job cap required')
    output = ROOT/'results'; output.mkdir(exist_ok=False)
    for i in range(3):
        seed = ROOT/'bundle/seeds'/str(i); chart = chart_for(seed)
        energy, gradient = cached_frozen(seed, chart)
        evaluate = CalculatorModeEvaluator(chart, make_factory(parameters={'source_directory': str(ASSETS)}, command=COMMAND),
            workdir=output/str(i)/'evaluations', calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6, result_validator=result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1', max_evaluation_attempts=24)
        cache = ExactCache(evaluate, energy, gradient)
        result = treatments(chart, json.loads((seed/'seed.json').read_bytes()), cache,
            record=lambda row: save(output/str(i)/('treatment_'+row['treatment']+'.json'), row))
        save(output/str(i)/'treatments.json', dict(index=i, frozen_energy_eV_cell=energy,
            frozen_gradient_eV_A=gradient.tolist(), treatments=result, total_new_attempts=cache.fresh,
            exact_cache_hits=cache.hits, old_path_labels_unchanged=True))
        print(json.dumps(dict(index=i, treatments=result, new_attempts=cache.fresh)), flush=True)
    save(output/'complete.json', dict(status='bounded_T_treatments_finished_not_stability_or_TS_certificate',
                                      job_id=accepted['job_id'], bundle_pin=pin))


def dispatch(mode, pin):
    validate(pin)
    proof = json.loads((ROOT/'preflight.json').read_bytes())
    if proof['bundle_pin'] != pin or len(proof['rows']) != 3 or proof['new_DFT_calls'] != 0:
        raise ValueError('all three native zero-DFT seed checks required')
    if mode == 'submit-held':
        if len(queue()) >= 2 or (ROOT/'intent.json').exists():
            raise ValueError('no free study slot or existing submission intent; never duplicate')
        cmd = ['sbatch', '--hold', '--parsable', '--job-name=hfo2-E270-pattern-T',
            f'--export=ALL,RUN_DFT=1,BUNDLE_SHA256={pin}', f'--output={ROOT}/job-%j.out', str(ROOT/'bundle/run.slurm')]
        save(ROOT/'intent.json', dict(command=cmd))
        reply = subprocess.check_output(cmd, text=True).strip(); job = reply.split(';')[0]
        if not job.isdigit():
            raise ValueError('ambiguous submission reply; inspect accounting, never resubmit')
        save(ROOT/'held.json', dict(job_id=job, reply=reply, bundle_pin=pin))
        return dict(job_id=job, state='held')
    held = json.loads((ROOT/'held.json').read_bytes()); job = held['job_id']
    raw = subprocess.check_output(['scontrol', 'show', 'job', '-o', job], text=True)
    fields = dict(re.findall(r'(\S+)=(\S+)', raw))
    expected = dict(JobId=job, JobName='hfo2-E270-pattern-T', JobState='PENDING', Reason='JobHeldUser',
        Partition='hfacnormal01', NumCPUs='32', NumTasks='1', TimeLimit='06:00:00')
    expected['CPUs/Task'] = '32'
    if (held['bundle_pin'] != pin or len(queue()) > 2 or fields.get('NumNodes') not in ('1', '1-1')
            or any(fields.get(k) != v for k, v in expected.items())):
        raise ValueError('held resource contract differs; leave held')
    save(ROOT/'accepted.json', dict(**held, raw_scontrol=raw))
    save(ROOT/'release_intent.json', dict(job_id=job))
    reply = subprocess.check_output(['scontrol', 'release', job], text=True).strip()
    result = dict(job_id=job, reply=reply, bundle_pin=pin)
    save(ROOT/'released.json', result)
    return result


def package(case, output):
    output.mkdir(parents=True, exist_ok=False); bundle = output/'bundle'; bundle.mkdir()
    for name in ('study.py', 'run.slurm', 'protocol.md'):
        if b'\r' in (case/name).read_bytes():
            raise ValueError('execution sources must be LF')
        shutil.copyfile(case/name, bundle/name)
    shutil.copyfile(case/'prepared/config.json', bundle/'config.json')
    shutil.copytree(case/'prepared/seeds', bundle/'seeds')
    files = {p.relative_to(bundle).as_posix(): sha256(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    save(bundle/'bundle_files.json', files)
    archive = output/'bundle.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar:
        tar.add(bundle, arcname='bundle')
    result = dict(bundle_pin=sha256(bundle/'bundle_files.json'), archive_sha256=sha256(archive), files=len(files), new_DFT_calls=0)
    save(output/'package_receipt.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'package', 'preflight', 'submit-held', 'release', 'run'))
    parser.add_argument('--repository', type=Path, default=REPOSITORY)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--bundle-pin')
    args = parser.parse_args()
    if args.mode == 'prepare':
        result = prepare(args.repository, args.output)
    elif args.mode == 'package':
        result = package(Path(__file__).resolve().parent, args.output)
    elif args.mode == 'preflight':
        result = preflight(args.bundle_pin)
    elif args.mode == 'run':
        result = run(args.bundle_pin)
    else:
        result = dispatch(args.mode, args.bundle_pin)
    if result is not None:
        print(json.dumps(result), flush=True)
