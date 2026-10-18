"""Finite E269 dispatch using the unchanged, hash-validated E268 runtime."""
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

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-midpoint-reverse-E269-r1')
OLD = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
OLD_PIN = '689a26fed0f47dede6d7c2f2230d3181fe5a16bbc5eee0356da6f033d89987cf'
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def command(args):
    return subprocess.check_output(args, text=True).strip()


def own_queue():
    raw = command(['squeue', '-u', 'iai806', '-h', '-o', '%i|%j|%T|%P|%C'])
    return [r for r in raw.splitlines() if r.split('|')[1].lower().startswith(('hfo2', 'varneb', 'vcneb'))]


def validate(pin):
    bundle = ROOT/'bundle'
    if sha256(bundle/'bundle_files.json') != pin:
        raise ValueError('external bundle pin changed')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or sha256(bundle/relative) != digest:
            raise ValueError('new bundle changed')
    # Read-only reuse: the original checker validates its original source archive,
    # all 375 source files, six physical assets and the actual ABACUS binary.
    if sha256(OLD/'bundle/bundle_files.json') != OLD_PIN:
        raise ValueError('original runtime manifest changed')
    old_files = json.loads((OLD/'bundle/bundle_files.json').read_bytes())
    checker = OLD/'bundle/endpoint_worker.py'
    if sha256(checker) != old_files['endpoint_worker.py']:
        raise ValueError('original runtime checker changed')
    spec = importlib.util.spec_from_file_location('E268_fixed_runtime', checker)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    worker.validate_bundle(OLD, OLD_PIN)
    config = json.loads((bundle/'config.json').read_bytes())
    if (config['channels'] != list(CHANNELS) or config['physical_input_sha256'] != CONTRACT
            or config['target_strain'] != .005 or config['source_strain'] != .01
            or config['runtime_root'] != str(OLD) or config['runtime_bundle_sha256'] != OLD_PIN):
        raise ValueError('fixed continuation definition changed')
    for channel in CHANNELS:
        seed = bundle/'seeds'/channel
        manifest = json.loads((seed/'manifest.json').read_bytes())
        if (manifest['channel'] != channel or manifest['target_strain'] != .005
                or manifest['source_strain'] != .01 or manifest['physical_input_sha256'] != CONTRACT
                or manifest['cell_scale_A'] != 5.12968067458423
                or any(sha256(seed/n) != h for n, h in manifest['files_sha256'].items())):
            raise ValueError('prepared continuation seed changed')
        old_parameters = json.loads((OLD/'path_bundle/seeds'/channel/'factory_parameters.json').read_bytes())
        parameters = json.loads((seed/'factory_parameters.json').read_bytes())
        if parameters != old_parameters:
            raise ValueError('same midpoint endpoint caches required; physical overrides prohibited')
    return worker


def preflight(pin):
    worker = validate(pin)
    out = ROOT/'preflight'; out.mkdir(exist_ok=False)
    old_calculate = FixedHfo2Calculator.calculate
    def prohibit(*args, **kwargs):
        raise RuntimeError('preflight prohibits any fresh DFT')
    FixedHfo2Calculator.calculate = prohibit
    records = []
    try:
        for channel in CHANNELS:
            seed = ROOT/'bundle/seeds'/channel
            images = read(seed/'seed.traj', index=':')
            if len(images) != 9:
                raise ValueError('seven interior images and two endpoints required')
            boundary = clamped_plane_vcneb_boundary(12, images[0].cell.array, allow_tilt=True)
            boundary.validate_images(images); validate_periodic_path_lift(images)
            validate_path_geometry(images, cell_scale=5.12968067458423,
                                   minimum_distance=1.6, maximum_deformation=.25)
            factory = make_clamped_seed_cached_factory(
                parameters=json.loads((seed/'factory_parameters.json').read_bytes()),
                command=f'mpirun -np 32 {worker.BIN}')
            for i, image in enumerate(images):
                folder = out/channel/f'image_{i:04d}'; folder.mkdir(parents=True)
                image.calc = factory(i, image, folder/'calculator')
                if i in (0, 8):
                    image.get_potential_energy(); image.get_forces(); image.get_stress()
                write_structure(folder/'STRU', image); geometry_roundtrip(folder/'STRU', image)
                for name in ('INPUT', 'KPT'):
                    shutil.copyfile(worker.PHYSICS/name, folder/name)
                records.append(dict(channel=channel, image=i, STRU_sha256=sha256(folder/'STRU'),
                                    exact_endpoint_cache_checked=i in (0, 8)))
    finally:
        FixedHfo2Calculator.calculate = old_calculate
    result = dict(status='two_native_reverse_continuations_ready', bundle_sha256=pin,
        decks=records, new_DFT_calls=0, input_parameters_changed=False,
        independent_forecast_test=False)
    save(ROOT/'preflight.json', result)
    return result


def validate_held(raw, index, job):
    fields = dict(re.findall(r'(\S+)=(\S+)', raw))
    expected = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser',
        Partition='hfacnormal01', NumCPUs='32', NumTasks='1',
        TimeLimit='1-00:00:00', JobName=f'hfo2-E269-midpoint-{index}')
    expected['CPUs/Task'] = '32'
    if fields.get('NumNodes') not in ('1', '1-1') or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError('held job differs from one-node32CPU contract; leave held')


def dispatch(mode, pin):
    validate(pin)
    proof = json.loads((ROOT/'preflight.json').read_bytes())
    if proof['bundle_sha256'] != pin or len(proof['decks']) != 18 or proof['new_DFT_calls'] != 0:
        raise ValueError('all18 native decks required')
    if mode == 'submit-held':
        if own_queue() or any((ROOT/f'{name}_{i}.json').exists()
                              for name in ('intent', 'held', 'accepted', 'released') for i in (0, 1)):
            raise ValueError('existing jobs or intents: inspect original handles, never resubmit')
        for i, channel in enumerate(CHANNELS):
            if len(own_queue()) >= 2:
                raise ValueError('two-allocation cap')
            cmd = ['sbatch', '--hold', '--parsable', f'--job-name=hfo2-E269-midpoint-{i}',
                f'--export=ALL,RUN_DFT=1,CHANNEL_ID={i},BUNDLE_SHA256={pin}',
                f'--output={ROOT}/path-%j.out', str(ROOT/'bundle/paths.slurm')]
            save(ROOT/f'intent_{i}.json', dict(command=cmd, bundle_sha256=pin))
            reply = command(cmd); job = reply.split(';')[0]
            if not job.isdigit():
                raise ValueError('ambiguous reply: do not resubmit')
            row = dict(job_id=job, channel=channel, bundle_sha256=pin, reply=reply,
                       raw_scontrol=command(['scontrol', 'show', 'job', '-o', job]))
            save(ROOT/f'held_{i}.json', row); print(json.dumps(row), flush=True)
    else:
        held = [json.loads((ROOT/f'held_{i}.json').read_bytes()) for i in (0, 1)]
        if {r.split('|')[0] for r in own_queue()} != {r['job_id'] for r in held}:
            raise ValueError('only the exact two original held jobs may be released')
        accepted = []
        for i, row in enumerate(held):
            if row['bundle_sha256'] != pin or (ROOT/f'accepted_{i}.json').exists():
                raise ValueError('pin differs or already accepted; inspect original jobs')
            raw = command(['scontrol', 'show', 'job', '-o', row['job_id']])
            validate_held(raw, i, row['job_id']); accepted.append((row, raw))
        for i, (row, raw) in enumerate(accepted):
            save(ROOT/f'accepted_{i}.json', dict(**row, validated_scontrol=raw))
        for i, (row, _) in enumerate(accepted):
            save(ROOT/f'release_intent_{i}.json', dict(job_id=row['job_id']))
            reply = command(['scontrol', 'release', row['job_id']])
            report = dict(job_id=row['job_id'], reply=reply,
                raw_scontrol=command(['scontrol', 'show', 'job', '-o', row['job_id']]))
            save(ROOT/f'released_{i}.json', report); print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'submit-held', 'release', 'run-check'))
    parser.add_argument('--bundle-sha256', required=True)
    parser.add_argument('--channel-id', type=int, choices=(0, 1))
    args = parser.parse_args()
    if args.mode == 'preflight':
        print(json.dumps(preflight(args.bundle_sha256)))
    elif args.mode == 'run-check':
        validate(args.bundle_sha256)
        row = json.loads((ROOT/f'accepted_{args.channel_id}.json').read_bytes())
        expected = dict(SLURM_JOB_ID=row['job_id'], SLURM_JOB_NAME=f'hfo2-E269-midpoint-{args.channel_id}',
            SLURM_JOB_PARTITION='hfacnormal01', SLURM_JOB_NUM_NODES='1', SLURM_NTASKS='1',
            SLURM_CPUS_PER_TASK='32')
        if (row['bundle_sha256'] != args.bundle_sha256
                or any(os.environ.get(k) != v for k, v in expected.items())):
            raise ValueError('only the accepted exact32CPU allocation may run')
    else:
        dispatch(args.mode, args.bundle_sha256)
