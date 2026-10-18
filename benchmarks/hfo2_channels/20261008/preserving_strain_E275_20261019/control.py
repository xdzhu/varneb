"""Native zero-SCF admission and exact held dispatch for two E275 bands."""
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
from examples.hfo2_fixed_input_factory import CONTRACT, FixedHfo2Calculator, make_clamped_seed_cached_factory
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb import clamped_plane_vcneb_boundary, validate_path_geometry, validate_periodic_path_lift

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-strain-E275-r1')
PREVIOUS = ROOT.with_name('20261018-varneb-midpoint-reverse-E269-r1')
PREVIOUS_PIN = '381ae9787cca91b1435b121dafd74ee072d1aa3bd22335f05a7e7ba42bfc7ae1'
SOURCE = ROOT.with_name('20261019-varneb-preserving-connection-E273-r1')
LABELS = ('strain_0050', 'strain_0100')


def save(path, data):
    with Path(path).open('x', encoding='utf-8', newline='\n') as h:
        json.dump(data, h, indent=2, allow_nan=False); h.write('\n')


def command(args): return subprocess.check_output(args, text=True).strip()


def own_queue():
    raw = command(['squeue', '-u', 'iai806', '-h', '-o', '%i|%j|%T|%P|%C'])
    return [r for r in raw.splitlines() if r.split('|')[1].lower().startswith(('hfo2', 'varneb', 'vcneb'))]


def validate(pin):
    bundle = ROOT/'bundle'
    if sha256(bundle/'bundle_files.json') != pin: raise ValueError('external bundle pin differs')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('executed input/source changed')
    old = PREVIOUS/'bundle'
    if sha256(old/'bundle_files.json') != PREVIOUS_PIN: raise ValueError('previous pinned source changed')
    files = json.loads((old/'bundle_files.json').read_bytes())
    if sha256(old/'control.py') != files['control.py']: raise ValueError('previous validator changed')
    spec = importlib.util.spec_from_file_location('E275_original_validator', old/'control.py')
    previous = importlib.util.module_from_spec(spec); spec.loader.exec_module(previous)
    worker = previous.validate(PREVIOUS_PIN)  # Immutable E268375 files, binary and all six assets.
    raw = command(['sacct', '-n', '-X', '-P', '-j', '28949435', '--format=JobIDRaw,State,ExitCode'])
    rows = [r.split('|') for r in raw.splitlines() if r.split('|')[0] == '28949435']
    if len(rows) != 1 or rows[0][1:3] != ['COMPLETED', '0:0']:
        raise ValueError('the exact source job must be completed, not live or timed out')
    meta = json.loads((bundle/'seeds/preparation.json').read_bytes())
    if (meta['protocol_sha256'] != sha256(bundle/'protocol.md') or meta['maximum_updates'] != 60
            or meta['maximum_new_interior_SCFs_per_chain'] != 427 or meta['walltime'] != '18:00:00'
            or meta['new_DFT_calls'] != 0 or meta['scheduler_mutations'] != 0):
        raise ValueError('finite protocol changed')
    for label, strain in zip(LABELS, (.005, .01)):
        seed = bundle/'seeds'/label; row = json.loads((seed/'manifest.json').read_bytes())
        if (row['label'] != label or row['target_strain'] != strain or row['source_strain'] != 0.
                or row['source_job'] != '28949435' or row['source_step'] != 9
                or row['physical_input_sha256'] != CONTRACT or row['free_interiors'] != list(range(1, 8))
                or row['retained_pattern_constraints_in_NEB'] != [] or row['independent_forecast_test']
                or row['cell_scale_A'] != 5.12968067458423
                or sha256(SOURCE/'path/band/vcneb_summary.json') != row['source_summary_sha256']):
            raise ValueError('same source, boundary and seven free interiors are required')
        for name, digest in row['files_sha256'].items():
            if sha256(seed/name) != digest: raise ValueError('seed changed')
        parameters = json.loads((seed/'factory_parameters.json').read_bytes())
        if (set(parameters) != {'source_directory', 'seed_cache_records'}
                or parameters['source_directory'] != str(worker.PHYSICS)
                or len(parameters['seed_cache_records']) != 9
                or parameters['seed_cache_records'][1:8] != [None]*7
                or any(parameters['seed_cache_records'][i] is None for i in (0, 8))):
            raise ValueError('only exact endpoint caches, no physics overrides/interior reuse')
    return worker


def preflight(pin):
    worker = validate(pin); output = ROOT/'preflight'; output.mkdir(exist_ok=False)
    old = FixedHfo2Calculator.calculate
    def prohibit(*args, **kwargs): raise RuntimeError('preflight forbids all fresh DFT')
    FixedHfo2Calculator.calculate = prohibit
    decks = []
    try:
        for label in LABELS:
            seed = ROOT/'bundle/seeds'/label; images = read(seed/'seed.traj', ':')
            if len(images) != 9 or any(a.get_chemical_symbols() != ['Hf']*4+['O']*8 for a in images):
                raise ValueError('ordered Hf4O8 nine-image source required')
            boundary = clamped_plane_vcneb_boundary(12, images[0].cell.array, allow_tilt=True)
            boundary.validate_images(images); validate_periodic_path_lift(images)
            validate_path_geometry(images, cell_scale=5.12968067458423,
                                   minimum_distance=1.6, maximum_deformation=.25)
            factory = make_clamped_seed_cached_factory(
                parameters=json.loads((seed/'factory_parameters.json').read_bytes()),
                command=f'mpirun -np 32 {worker.BIN}')
            for i, image in enumerate(images):
                folder = output/label/f'image_{i:04d}'; folder.mkdir(parents=True)
                image.calc = factory(i, image, folder/'calculator')
                if i in (0, 8): image.get_potential_energy(); image.get_forces(); image.get_stress()
                write_structure(folder/'STRU', image); geometry_roundtrip(folder/'STRU', image)
                for n in ('INPUT', 'KPT'): shutil.copyfile(worker.PHYSICS/n, folder/n)
                decks.append(dict(label=label, image=i, exact_endpoint_cache_checked=i in (0, 8),
                                  STRU_sha256=sha256(folder/'STRU')))
    finally: FixedHfo2Calculator.calculate = old
    result = dict(status='two_native_strain_transfers_ready', bundle_sha256=pin,
        new_DFT_calls=0, decks=decks, source_job='28949435', physics_changed=False)
    save(ROOT/'preflight.json', result); return result


def validate_held(raw, index, job):
    fields = dict(re.findall(r'(\S+)=(\S+)', raw))
    expected = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser', Partition='hfacnormal01',
        NumCPUs='32', NumTasks='1', TimeLimit='18:00:00', JobName=f'hfo2-E275-strain-{index}')
    expected['CPUs/Task'] = '32'
    if fields.get('NumNodes') not in ('1', '1-1') or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError('only exact one-node32CPU18h held identity may be released')


def dispatch(mode, pin):
    validate(pin); proof = json.loads((ROOT/'preflight.json').read_bytes())
    if proof['bundle_sha256'] != pin or proof['new_DFT_calls'] != 0 or len(proof['decks']) != 18:
        raise ValueError('complete native zero-SCF preflight required')
    if mode == 'submit-held':
        if own_queue() or any((ROOT/f'{n}_{i}.json').exists()
                              for n in ('intent', 'held', 'accepted', 'released') for i in (0, 1)):
            raise ValueError('existing study jobs/intents: inspect original handles, never duplicate')
        for i, label in enumerate(LABELS):
            if len(own_queue()) >= 2: raise ValueError('two study allocation cap')
            cmd = ['sbatch', '--hold', '--parsable', f'--job-name=hfo2-E275-strain-{i}',
                f'--export=ALL,RUN_DFT=1,CHANNEL_ID={i},BUNDLE_SHA256={pin}',
                f'--output={ROOT}/path-%j.out', str(ROOT/'bundle/paths.slurm')]
            save(ROOT/f'intent_{i}.json', dict(command=cmd, bundle_sha256=pin))
            reply = command(cmd); job = reply.split(';')[0]
            if not job.isdigit(): raise ValueError('ambiguous submission; inspect, never repeat')
            row = dict(job_id=job, label=label, bundle_sha256=pin, reply=reply,
                raw_scontrol=command(['scontrol', 'show', 'job', '-o', job]))
            save(ROOT/f'held_{i}.json', row); print(json.dumps(row), flush=True)
    else:
        held = [json.loads((ROOT/f'held_{i}.json').read_bytes()) for i in (0, 1)]
        if {r.split('|')[0] for r in own_queue()} != {r['job_id'] for r in held}:
            raise ValueError('only the original two held study jobs may be released')
        checked = []
        for i, row in enumerate(held):
            if row['bundle_sha256'] != pin or (ROOT/f'accepted_{i}.json').exists():
                raise ValueError('pin changed or already accepted; inspect original jobs')
            raw = command(['scontrol', 'show', 'job', '-o', row['job_id']])
            validate_held(raw, i, row['job_id']); checked.append((row, raw))
        for i, (row, raw) in enumerate(checked):
            save(ROOT/f'accepted_{i}.json', dict(**row, validated_scontrol=raw))
        for i, (row, _) in enumerate(checked):
            save(ROOT/f'release_intent_{i}.json', dict(job_id=row['job_id']))
            reply = command(['scontrol', 'release', row['job_id']])
            save(ROOT/f'released_{i}.json', dict(job_id=row['job_id'], reply=reply,
                raw_scontrol=command(['scontrol', 'show', 'job', '-o', row['job_id']])))


def run_check(pin, index):
    validate(pin); row = json.loads((ROOT/f'accepted_{index}.json').read_bytes())
    proof = json.loads((ROOT/'preflight.json').read_bytes())
    expected = dict(RUN_DFT='1', SLURM_JOB_ID=row['job_id'], SLURM_JOB_NAME=f'hfo2-E275-strain-{index}',
        SLURM_JOB_PARTITION='hfacnormal01', SLURM_JOB_NUM_NODES='1', SLURM_NTASKS='1', SLURM_CPUS_PER_TASK='32',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', I_MPI_HYDRA_BOOTSTRAP='fork')
    if (row['bundle_sha256'] != pin or proof['bundle_sha256'] != pin or proof['new_DFT_calls'] != 0
            or len(proof['decks']) != 18 or len(own_queue()) > 2
            or any(os.environ.get(k) != v for k, v in expected.items())):
        raise ValueError('only the exact accepted allocation with native preflight may run')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('preflight', 'submit-held', 'release', 'run-check'))
    p.add_argument('--bundle-sha256', required=True); p.add_argument('--channel-id', type=int, choices=(0, 1))
    args = p.parse_args()
    if args.mode == 'preflight': print(json.dumps(preflight(args.bundle_sha256)))
    elif args.mode == 'run-check': run_check(args.bundle_sha256, args.channel_id)
    else: dispatch(args.mode, args.bundle_sha256)
