"""Bounded two-path preflight and exact-handle dispatch after forecast freeze."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import FixedHfo2Calculator, make_clamped_seed_cached_factory
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from vcneb import clamped_plane_vcneb_boundary, validate_path_geometry, validate_periodic_path_lift
from vcneb.forecast_freeze import load_frozen_forecast_panel

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')
ENDPOINT_PIN = '689a26fed0f47dede6d7c2f2230d3181fe5a16bbc5eee0356da6f033d89987cf'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2); f.write('\n')


def command(args): return subprocess.check_output(args, text=True).strip()


def queue():
    raw = command(['squeue', '-u', 'iai806', '-h', '-o', '%i|%j|%T|%P|%C'])
    return [r for r in raw.splitlines() if r.split('|')[1].lower().startswith(('hfo2', 'varneb', 'vcneb'))]


def validate(pin):
    bundle = ROOT/'path_bundle'
    if sha(bundle/'bundle_files.json') != pin: raise ValueError('external path bundle pin changed')
    files = json.loads((bundle/'bundle_files.json').read_bytes())
    for name, digest in files.items():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts or sha(bundle/p) != digest:
            raise ValueError('path bundle changed')
    # Reuse the original fixed-source/physical validation without editing it.
    old = ROOT/'bundle/endpoint_worker.py'
    old_files = json.loads((ROOT/'bundle/bundle_files.json').read_bytes())
    if sha(ROOT/'bundle/bundle_files.json') != ENDPOINT_PIN or sha(old) != old_files['endpoint_worker.py']:
        raise ValueError('original endpoint checker drift')
    s = importlib.util.spec_from_file_location('E268_original_worker', old)
    w = importlib.util.module_from_spec(s); s.loader.exec_module(w)
    w.validate_bundle(ROOT, ENDPOINT_PIN)
    config = json.loads((bundle/'config.json').read_bytes())
    forecasts = load_frozen_forecast_panel(bundle/'forecasts', expected_manifest_sha256=config['forecast_manifest_sha256'])
    if (len(forecasts) != 138 or {r.case_id for r in forecasts} != set(CHANNELS)
            or config['channels'] != list(CHANNELS)):
        raise ValueError('complete two-path forecast panel required')
    for channel in CHANNELS:
        seed = bundle/'seeds'/channel
        m = json.loads((seed/'manifest.json').read_bytes())
        if (m['forecast_manifest_sha256'] != config['forecast_manifest_sha256']
                or m['channel'] != channel or m['strain'] != .005
                or m['recipe']['steps'] != 80 or m['recipe']['fmax'] != .10
                or m['recipe']['climb'] or m['cell_scale_A'] != 5.12968067458423
                or any(sha(seed/name) != h for name, h in m['files_sha256'].items())):
            raise ValueError('registered path geometry or fixed recipe changed')
    return config, w


def preflight(pin):
    config, worker = validate(pin)
    out = ROOT/'path_preflight'; out.mkdir(exist_ok=False)
    old_calculate = FixedHfo2Calculator.calculate
    def prohibit(*args, **kwargs): raise RuntimeError('preflight prohibits fresh DFT')
    FixedHfo2Calculator.calculate = prohibit
    records = []
    try:
        for channel in CHANNELS:
            seed = ROOT/'path_bundle/seeds'/channel
            m = json.loads((seed/'manifest.json').read_bytes())
            images = read(seed/'seed.traj', index=':')
            boundary = clamped_plane_vcneb_boundary(12, images[0].cell.array, allow_tilt=True)
            boundary.validate_images(images); validate_periodic_path_lift(images)
            validate_path_geometry(images, cell_scale=m['cell_scale_A'], minimum_distance=1.6, maximum_deformation=.25)
            factory = make_clamped_seed_cached_factory(parameters=json.loads((seed/'factory_parameters.json').read_bytes()),
                                                       command=f'mpirun -np 32 {worker.BIN}')
            for i, image in enumerate(images):
                folder = out/channel/f'image_{i:04d}'; folder.mkdir(parents=True)
                image.calc = factory(i, image, folder/'calculator')
                if i in (0, 8):
                    image.get_potential_energy(); image.get_forces(); image.get_stress()
                write_structure(folder/'STRU', image); geometry_roundtrip(folder/'STRU', image)
                for name in ('INPUT', 'KPT'): shutil.copyfile(worker.PHYSICS/name, folder/name)
                records.append(dict(channel=channel, image=i, STRU_sha256=sha(folder/'STRU'),
                                    exact_endpoint_cache_checked=i in (0, 8)))
    finally:
        FixedHfo2Calculator.calculate = old_calculate
    result = dict(status='two_native_nine_image_inputs_ready', path_bundle_sha256=pin,
        forecast_manifest_sha256=config['forecast_manifest_sha256'], decks=records,
        new_DFT_calls=0, target_interior_energies_read=False, physical_parameters_changed=False)
    save(ROOT/'path_preflight.json', result); return result


def dispatch(mode, pin):
    config, _ = validate(pin)
    proof = json.loads((ROOT/'path_preflight.json').read_bytes())
    if proof['path_bundle_sha256'] != pin or len(proof['decks']) != 18 or proof['new_DFT_calls'] != 0:
        raise ValueError('all18 native inputs required')
    if mode == 'submit-held':
        if queue() or any((ROOT/f'path_{name}_{i}.json').exists() for i in (0, 1)
                          for name in ('intent', 'held', 'accepted', 'released')):
            raise ValueError('existing allocation or intent; reconcile before new submission')
        for i, channel in enumerate(CHANNELS):
            if len(queue()) >= 2: raise ValueError('two allocation cap')
            cmd = ['sbatch', '--hold', '--parsable', f'--job-name=hfo2-E268-path-{i}',
                   f'--export=ALL,RUN_DFT=1,CHANNEL_ID={i},PATH_BUNDLE_SHA256={pin}',
                   f'--output={ROOT}/path-%j.out', str(ROOT/'path_bundle/paths.slurm')]
            save(ROOT/f'path_intent_{i}.json', dict(command=cmd, path_bundle_sha256=pin))
            reply = command(cmd); job = reply.split(';')[0]
            if not job.isdigit(): raise ValueError('ambiguous job reply; never resubmit')
            row = dict(job_id=job, channel=channel, path_bundle_sha256=pin,
                       forecast_manifest_sha256=config['forecast_manifest_sha256'], reply=reply,
                       raw_scontrol=command(['scontrol', 'show', 'job', '-o', job]))
            save(ROOT/f'path_held_{i}.json', row); print(json.dumps(row), flush=True)
    else:
        if len(queue()) != 2: raise ValueError('exactly the two held allocations required')
        accepted = []
        for i in (0, 1):
            if (ROOT/f'path_accepted_{i}.json').exists(): raise ValueError('already accepted; inspect original handles')
            row = json.loads((ROOT/f'path_held_{i}.json').read_bytes())
            raw = command(['scontrol', 'show', 'job', '-o', row['job_id']])
            fields = dict(re.findall(r'(\S+)=(\S+)', raw))
            expected = dict(JobId=row['job_id'], JobState='PENDING', Reason='JobHeldUser',
                Partition='hfacnormal01', NumCPUs='32', NumTasks='1', TimeLimit='1-00:00:00',
                JobName=f'hfo2-E268-path-{i}')
            expected['CPUs/Task'] = '32'
            if fields.get('NumNodes') not in ('1', '1-1') or any(fields.get(k) != v for k, v in expected.items()):
                raise ValueError('held path resources differ; leave held')
            accepted.append((i, row, raw))
        for i, row, raw in accepted:
            save(ROOT/f'path_accepted_{i}.json', dict(**row, validated_scontrol=raw))
        for i, row, _ in accepted:
            save(ROOT/f'path_release_intent_{i}.json', dict(job_id=row['job_id']))
            reply = command(['scontrol', 'release', row['job_id']])
            result = dict(job_id=row['job_id'], reply=reply,
                          raw_scontrol=command(['scontrol', 'show', 'job', '-o', row['job_id']]))
            save(ROOT/f'path_released_{i}.json', result); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('preflight', 'submit-held', 'release', 'run-check'))
    p.add_argument('--bundle-sha256', required=True); p.add_argument('--channel-id', type=int, choices=(0, 1))
    a = p.parse_args()
    if a.mode == 'preflight': print(json.dumps(preflight(a.bundle_sha256)))
    elif a.mode == 'run-check':
        validate(a.bundle_sha256)
        row = json.loads((ROOT/f'path_accepted_{a.channel_id}.json').read_bytes())
        expected = dict(SLURM_JOB_ID=row['job_id'], SLURM_JOB_NAME=f'hfo2-E268-path-{a.channel_id}',
                        SLURM_JOB_PARTITION='hfacnormal01', SLURM_JOB_NUM_NODES='1',
                        SLURM_NTASKS='1', SLURM_CPUS_PER_TASK='32')
        if row['path_bundle_sha256'] != a.bundle_sha256 or any(os.environ.get(k) != v for k, v in expected.items()):
            raise ValueError('only the accepted original32CPU path allocation may run')
    else: dispatch(a.mode, a.bundle_sha256)
