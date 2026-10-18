"""Collect only terminal E269 paths, retaining both initializations separately.

Uses the frozen material runtime's native readers. No calculator, optimizer,
submission, source overwrite or re-scoring of the E268 prediction is invoked.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.audit_hfo2_clamped_terminal import completed_step, validate_terminal
from scripts.export_hfo2_clamped_observation import export, REGISTERED_FACTORIES

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-midpoint-reverse-E269-r1')
PIN = '381ae9787cca91b1435b121dafd74ee072d1aa3bd22335f05a7e7ba42bfc7ae1'
JOBS = ('28946180', '28946181')
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')


def save(path, value):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def accounting(job):
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True)
    rows = [r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] == job]
    if len(rows) != 1 or len(rows[0]) < 10 or rows[0][6] != '1':
        raise ValueError('complete accounting for the exact one-node handle required')
    return rows[0]


def require_completed(row):
    if len(row) < 10 or row[1:7] != ['COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1']:
        raise ValueError('exact original one-node32CPU job must be completed0:0')


def audit(index, output):
    if output.exists():
        raise FileExistsError(output)
    accepted = json.loads((ROOT/f'accepted_{index}.json').read_bytes())
    channel, job = CHANNELS[index], JOBS[index]
    if (accepted['job_id'] != job or accepted['channel'] != channel
            or accepted['bundle_sha256'] != PIN):
        raise ValueError('only the original E269 handles and fixed bundle are accepted')
    script = ROOT/'bundle/paths.slurm'
    files = json.loads((ROOT/'bundle/bundle_files.json').read_bytes())
    if (sha256(ROOT/'bundle/bundle_files.json') != PIN
            or sha256(script) != files['paths.slurm']):
        raise ValueError('executed bundle changed')
    row = accounting(job)
    require_completed(row)
    work = ROOT/'paths'/channel/'band'
    # A live handle, missing summary or failed run must not become a completed
    # result. Timeout/failure diagnosis remains with the original handle.
    summary = json.loads((work/'vcneb_summary.json').read_bytes())
    validate_terminal(job, '|'.join(row[:6]), summary, (work/'vcneb_failure.json').exists(),
                      outcome='converged' if summary.get('converged') else 'step-cap')
    step = completed_step(work/'vcneb.opt.log')
    calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
    if step > 80 or len(calls) != 7*(step+1):
        raise ValueError('registered update/call count differs; inspect without relaunching')
    records = []
    for path in calls:
        call = json.loads(path.read_bytes()); folder = path.parent
        if (any(call['input_sha256'][n] != h or sha256(folder/n) != h for n, h in CONTRACT.items())
                or sha256(folder/'STRU') != call['input_sha256']['STRU']
                or sha256(folder/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']):
            raise ValueError('native input or output changed')
        parsed = audited_results(folder)
        if any(not np.allclose(parsed[k], call['results'][k], atol=1e-12, rtol=0) for k in parsed):
            raise ValueError('complete native E/F/stress does not replay')
        seconds = call['elapsed_seconds']
        if not np.isfinite(seconds) or seconds <= 0:
            raise ValueError('invalid SCF timing')
        records.append(dict(directory=str(folder), call_audit_sha256=sha256(path),
            input_sha256=call['input_sha256'], raw_log_sha256=call['raw_log_sha256'], elapsed_seconds=seconds))
    REGISTERED_FACTORIES[sha256(script)] = 'examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory'
    output.mkdir()
    observation = export(work, step, output/f'step_{step:04d}', job, production_script=script)
    if bool(observation['ordinary_residual_pass']) != bool(summary['converged']):
        raise ValueError('replayed ordinary criterion disagrees with terminal summary')
    for image in observation['raw_image_evaluations']:
        source = Path(image['raw_source'])
        destination = output/f'step_{step:04d}'/'raw'/f"image_{image['image_index']:04d}"
        (destination/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
            shutil.copyfile(source/name, destination/name)
    for i, call in enumerate(records):
        source = Path(call['directory']); destination = output/'all_calls'/f'{i:04d}'
        (destination/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT', 'KPT', 'STRU', 'call_audit.json', 'OUT.ABACUS/running_scf.log'):
            shutil.copyfile(source/name, destination/name)
    for name in ('vcneb_summary.json', 'vcneb.opt.log'):
        shutil.copyfile(work/name, output/name)
    result = dict(status='ordinary_reverse_midpoint_audited' if summary['converged'] else 'bounded_step_cap_audited',
        job_id=job, channel=channel, scheduler_row=row, step=step,
        ordinary_converged=bool(summary['converged']), final_fmax_eV_A=observation['replayed_fmax_eV_A'],
        sampled_peak_eV_fu=observation['sampled_forward_barrier_meV_fu']/1000,
        observation_sha256=sha256(output/f'step_{step:04d}'/'observation.json'),
        new_SCF_calls=len(calls), SCF_core_hours=sum(c['elapsed_seconds'] for c in records)*32/3600,
        allocation_core_hours=float(row[7])*32/3600, calls=records,
        all_native_calls_reparsed=True, all_six_assets_checked=True, genuine_MPI32_checked=True,
        bundle_sha256=PIN, initializer_strain=.01, target_strain=.005,
        original_E268_prediction_score_unchanged=True, independent_prediction_test=False,
        physical_barrier_error_bound=None, stationary_TS_certified=False, new_DFT_calls=0)
    save(output/'audit.json', result)
    return result


def collect():
    output = ROOT/'completed_paths'; archive = ROOT/'completed_paths.tar.gz'
    if output.exists() or archive.exists():
        raise FileExistsError('export exists: inspect it, never overwrite or resubmit')
    # Verify all original runtime/code/physical bytes once, outside any DFT job.
    spec = importlib.util.spec_from_file_location('E269_frozen_control', ROOT/'bundle/control.py')
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control)
    control.validate(PIN)
    for i, job in enumerate(JOBS):
        row = accounting(job)
        require_completed(row)
    output.mkdir()
    results = [audit(i, output/name) for i, name in enumerate(('preserving', 'T'))]
    files = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'files.json', files)
    with tarfile.open(archive, 'x:gz') as tar:
        tar.add(output, arcname=output.name)
    receipt = dict(archive=str(archive), archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
        exported_files=len(files), results=[{k: v for k, v in r.items() if k != 'calls'} for r in results],
        new_DFT_calls=0, scheduler_mutations=0)
    save(ROOT/'collection_receipt.json', receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(collect(), allow_nan=False))
