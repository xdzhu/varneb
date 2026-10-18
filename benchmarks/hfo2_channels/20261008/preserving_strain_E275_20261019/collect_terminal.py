"""Collect one original E275 terminal; never launch, retry or replace runtime."""
from __future__ import annotations
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
from scripts.audit_hfo2_clamped_terminal import completed_step, fresh_interior_call_count, validate_terminal
from scripts.export_hfo2_clamped_observation import export, REGISTERED_FACTORIES

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-strain-E275-r1')
PIN = '00ea6f9f82e50aeea60cfd560ecdcc76cec68b3308a0683e8bdf0fe922dfe8e3'
JOBS = {'strain_0050': '28951226', 'strain_0100': '28951227'}
FACTORY = 'examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory'


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as h:
        json.dump(value, h, indent=2, allow_nan=False); h.write('\n')


def terminal_row(raw, job):
    rows = [r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] == job]
    if (len(rows) != 1 or len(rows[0]) < 10
            or rows[0][1:7] != ['COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1']):
        raise ValueError('only the exact completed original one-node32CPU handle is a terminal result')
    return rows[0]


def admit(job, row, summary, step, failure_exists, preflight, producer):
    REGISTERED_FACTORIES[producer] = FACTORY
    validate_terminal(job, '|'.join(row[:6]), summary, failure_exists,
                      outcome='converged' if summary.get('converged') else 'step-cap')
    if not 0 <= step <= 60 or not summary['converged'] and step != 60:
        raise ValueError('ordinary convergence or the original60-update cap required')
    initial, expected = fresh_interior_call_count(step, preflight, producer)
    if initial != 7 or expected != 7*(step+1):
        raise ValueError('exact endpoints only; every interior evaluation is fresh')
    return expected


def copy_native(source, destination):
    (destination/'OUT.ABACUS').mkdir(parents=True)
    for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
        shutil.copyfile(source/name, destination/name)
    for name in ('call_audit.json', 'abacus.out', 'abacus.err'):
        if (source/name).is_file(): shutil.copyfile(source/name, destination/name)


def collect(label):
    job = JOBS[label]; index = tuple(JOBS).index(label)
    accepted = json.loads((ROOT/f'accepted_{index}.json').read_bytes())
    if (accepted['job_id'] != job or accepted['label'] != label or accepted['bundle_sha256'] != PIN):
        raise ValueError('only the two originally accepted handles and bundle may be collected')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True)
    row = terminal_row(raw, job)  # Fail before exports if live, missing, timeout or failed.
    work = ROOT/'paths'/label/'band'
    script = ROOT/'bundle/paths.slurm'; producer = sha256(script)
    manifest = json.loads((ROOT/'bundle/bundle_files.json').read_bytes())
    if sha256(ROOT/'bundle/bundle_files.json') != PIN or producer != manifest['paths.slurm']:
        raise ValueError('executed source pin changed')
    summary = json.loads((work/'vcneb_summary.json').read_bytes())
    preflight = json.loads((work/'vcneb_preflight.json').read_bytes())
    step = completed_step(work/'vcneb.opt.log')
    expected = admit(job, row, summary, step, (work/'vcneb_failure.json').exists(), preflight, producer)
    spec = importlib.util.spec_from_file_location('E275_terminal_validator', ROOT/'bundle/control.py')
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control); control.validate(PIN)
    destination = ROOT/'exports'/label; output = destination/'native_terminal'
    if destination.exists(): raise FileExistsError('original export exists; inspect it without overwriting')
    records = []; counts = {f'image_{i:04d}': 0 for i in range(1, 8)}
    for path in sorted(work.glob('image_*/scf_*/call_audit.json')):
        folder = path.parent; image = folder.parent.name
        if image not in counts: raise ValueError('fixed endpoints must not be recalculated')
        counts[image] += 1; call = json.loads(path.read_bytes())
        if (any(call['input_sha256'][n] != h or sha256(folder/n) != h for n, h in CONTRACT.items())
                or sha256(folder/'STRU') != call['input_sha256']['STRU']
                or sha256(folder/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']):
            raise ValueError('native six inputs, structure or raw log changed')
        parsed = audited_results(folder)
        if any(not np.allclose(parsed[k], call['results'][k], atol=1e-12, rtol=0) for k in parsed):
            raise ValueError('native E/F/stress does not replay')
        if not np.isfinite(call['elapsed_seconds']) or call['elapsed_seconds'] <= 0:
            raise ValueError('native elapsed time must be positive')
        records.append(dict(relative_directory=folder.relative_to(work).as_posix(),
            call_audit_sha256=sha256(path), input_sha256=call['input_sha256'],
            raw_log_sha256=call['raw_log_sha256'], seconds=call['elapsed_seconds']))
    if len(records) != expected or any(n != step+1 for n in counts.values()):
        raise ValueError('all actual per-image calls must match the fixed cache/update policy')
    output.mkdir(parents=True)
    obs = export(work, step, output/'observation', job, production_script=script)
    if obs['ordinary_residual_pass'] != summary['converged']:
        raise ValueError('native residual replay and terminal outcome disagree')
    for image in obs['raw_image_evaluations']:
        copy_native(Path(image['raw_source']), output/'observation/raw'/f"image_{image['image_index']:04d}")
    for call in records:
        copy_native(work/call['relative_directory'], output/'all_calls'/call['relative_directory'])
    for name in ('vcneb_summary.json', 'vcneb_preflight.json', 'vcneb.opt.log'):
        shutil.copyfile(work/name, output/name)
    receipt = dict(status='ordinary_converged_terminal' if summary['converged'] else 'bounded_unconverged_terminal',
        job_id=job, label=label, target_strain=.005 if index == 0 else .01, bundle_pin=PIN,
        producer_sha256=producer, sidecar_source_sha256=sha256(Path(__file__)),
        scheduler_row_at_export=row, step=step, ordinary_converged=summary['converged'],
        new_SCF_calls=len(records), calls=records, all_new_native_calls_six_assets_MPI_EFS_replayed_on_HF=True,
        native_proxy_core_hours=sum(c['seconds'] for c in records)*32/3600,
        allocated_core_hours=float(row[7])*32/3600, prior_path_release_costs_additional=True,
        fmax_eV_A=obs['replayed_fmax_eV_A'], sampled_peak_meV_fu=obs['sampled_forward_barrier_meV_fu'],
        old_labels_and_scores_unchanged=True, stationary_TS_or_sampling_error_bound_certified=False,
        new_DFT_calls=0, scheduler_mutations=0)
    save(output/'receipt.json', receipt)
    files = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'files.json', files)
    archive = destination/'native_terminal.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar: tar.add(output, arcname='native_terminal')
    result = {k: v for k, v in receipt.items() if k != 'calls'}
    result.update(archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
                  exported_files=len(files), all_runtime_calls_exported=True)
    save(destination/'collection_receipt.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--label', choices=tuple(JOBS), required=True)
    print(json.dumps(collect(p.parse_args().label)))
