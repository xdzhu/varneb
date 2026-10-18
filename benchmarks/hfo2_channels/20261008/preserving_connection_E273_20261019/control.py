"""Exact-cache admission for one finite free-interior connection run."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from ase.io import read
from examples.hfo2_fixed_input_factory import FixedHfo2Calculator, make_clamped_resume_cached_factory
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb import VCNEB, clamped_plane_vcneb_boundary, validate_path_geometry, validate_periodic_path_lift

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-connection-E273-r1')
OLD = ROOT.with_name('20261018-varneb-physical-E268-r1')
PARENT = ROOT.with_name('20261018-varneb-preserving-branches-E272-r1')
T = ROOT.with_name('20261018-varneb-matched-pattern-E270-r1')
PARENT_PIN = 'd12af56eac84edf28e536c69b524d12681f876db1918c889fbf4d89eb59e06c2'
T_PIN = 'f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca'
SCALE = 5.12968067458423


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def pinned_files(bundle, pin):
    if sha256(bundle/'bundle_files.json') != pin:
        raise ValueError('executed bundle identity differs')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('immutable source/input changed')


def require_parent_terminal(raw):
    rows = [row.split('|') for row in raw.strip().splitlines() if row.split('|')[0] == '28948425']
    if len(rows) != 1 or rows[0][1:3] != ['COMPLETED', '0:0']:
        raise ValueError('wait for the original E272 whole job, not a live/timeout promotion')


def validate_held(raw, job):
    fields = dict(re.findall(r'(\S+)=(\S+)', raw))
    expected = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser',
        Partition='hfacnormal01', NumCPUs='32', NumTasks='1', TimeLimit='06:00:00',
        JobName='hfo2-E273-connection')
    expected['CPUs/Task'] = '32'
    if fields.get('NumNodes') not in ('1', '1-1') or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError('held allocation is not the declared one-node32CPU6h job')


def validate(pin):
    bundle = ROOT/'bundle'; pinned_files(bundle, pin); pinned_files(PARENT/'bundle', PARENT_PIN)
    accepted = json.loads((PARENT/'accepted.json').read_bytes())
    complete = json.loads((PARENT/'results/complete.json').read_bytes())
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', '28948425',
                                  '--format=JobIDRaw,State,ExitCode'], text=True)
    require_parent_terminal(raw)
    if (accepted['job_id'] != '28948425' or accepted['bundle_pin'] != PARENT_PIN
            or complete['job_id'] != '28948425' or complete['bundle_pin'] != PARENT_PIN
            or not 0 < complete['total_new_attempts'] <= 104):
        raise ValueError('original finite parent result required')
    pinned_files(T/'bundle', T_PIN)
    spec = importlib.util.spec_from_file_location('E273_original_runtime_check', T/'bundle/study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    study.validate(T_PIN)  # Also checks the original E268375files, binary and six assets.
    seed = bundle/'prepared'; meta = json.loads((seed/'manifest.json').read_bytes())
    if (meta['branch'] != 'zero_minus_J' or meta['free_interiors'] != list(range(1, 8))
            or meta['retained_pattern_constraints_in_NEB'] != [] or meta['fmax_target_eV_A'] != .10
            or meta['maximum_updates'] != 24 or meta['maximum_new_interior_SCFs'] != 168
            or meta['walltime'] != '06:00:00' or meta['maximum_study_allocations'] != 2
            or meta['pressure_GPa'] != 0 or meta['climb'] or meta['execution_admitted']
            or meta['new_DFT_calls'] != 0 or meta['scheduler_mutations'] != 0
            or meta['protocol_sha256'] != sha256(bundle/'protocol.md')):
        raise ValueError('single finite connection protocol differs')
    for name, digest in meta['files_sha256'].items():
        if sha256(seed/name) != digest:
            raise ValueError('prepared connection seed changed')
    return study


def preflight(pin):
    study = validate(pin); seed = ROOT/'bundle/prepared'
    output = ROOT/'preflight'; output.mkdir(exist_ok=False)
    images = read(seed/'seed.traj', ':')
    if len(images) != 9 or any(a.get_chemical_symbols() != ['Hf']*4+['O']*8 for a in images):
        raise ValueError('one ordered nine-image Hf4O8 chain required')
    boundary = clamped_plane_vcneb_boundary(12, read(seed/'substrate.vasp').cell.array, allow_tilt=True)
    boundary.validate_images(images); validate_periodic_path_lift(images)
    validate_path_geometry(images, cell_scale=SCALE, minimum_distance=1.6, maximum_deformation=.25)
    parameters = json.loads((seed/'factory_parameters.json').read_bytes())
    factory = make_clamped_resume_cached_factory(parameters=parameters, command=study.COMMAND)
    previous = FixedHfo2Calculator.calculate
    def prohibit(*args, **kwargs):
        raise RuntimeError('new DFT is forbidden during exact-cache preflight')
    FixedHfo2Calculator.calculate = prohibit
    decks = []
    try:
        for i, image in enumerate(images):
            folder = output/f'image_{i:04d}'; folder.mkdir()
            image.calc = factory(i, image, folder/'calculator')
            image.get_potential_energy(); image.get_forces(); image.get_stress()
            write_structure(folder/'STRU', image); geometry_roundtrip(folder/'STRU', image)
            for name in ('INPUT', 'KPT'): shutil.copyfile(study.ASSETS/name, folder/name)
            decks.append(dict(image=i, STRU_sha256=sha256(folder/'STRU'), exact_cache_checked=True))
        chain = VCNEB(images, cell_scale=SCALE, k=.2, climb=False, pressure=0., **boundary.vcneb_kwargs(images))
        residual = float(np.linalg.norm(chain.get_forces(), axis=1).max())
        expected = json.loads((seed/'connection_replay.json').read_bytes())['replaced_central_image']
        if (abs(residual-expected['replayed_fmax_eV_A']) > 7e-7
                or not np.allclose(chain.enthalpies, expected['energy_eV_cell'], atol=1e-10, rtol=0)):
            raise ValueError('native resumed polygon is not the prepared E/F/stress chain')
    finally:
        FixedHfo2Calculator.calculate = previous
    result = dict(status='nine_exact_native_caches_ready', bundle_pin=pin, decks=decks,
        initial_fmax_eV_A=residual, new_DFT_calls=0, parent_job='28948425',
        free_interiors=list(range(1, 8)), fixed_endpoints=[0, 8], physical_parameters_changed=False)
    save(ROOT/'preflight.json', result)
    return result


def run_check(pin):
    study = validate(pin); accepted = json.loads((ROOT/'accepted.json').read_bytes())
    proof = json.loads((ROOT/'preflight.json').read_bytes())
    expected = dict(RUN_DFT='1', SLURM_JOB_ID=accepted['job_id'], SLURM_JOB_NAME='hfo2-E273-connection',
        SLURM_JOB_PARTITION='hfacnormal01', SLURM_JOB_NUM_NODES='1', SLURM_NTASKS='1', SLURM_CPUS_PER_TASK='32',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    if (accepted['bundle_pin'] != pin or proof['bundle_pin'] != pin or proof['new_DFT_calls'] != 0
            or len(proof['decks']) != 9 or len(study.queue()) > 2
            or any(os.environ.get(k) != v for k, v in expected.items())):
        raise ValueError('only the exact accepted32CPU allocation and verified cache seed may run')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'run-check')); parser.add_argument('--bundle-pin', required=True)
    args = parser.parse_args()
    if args.mode == 'preflight': print(json.dumps(preflight(args.bundle_pin)))
    else: run_check(args.bundle_pin)
