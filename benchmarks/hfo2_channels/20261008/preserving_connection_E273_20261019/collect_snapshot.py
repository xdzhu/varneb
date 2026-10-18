"""Read-only export of one complete E273 frame, not a new calculation.

Registers only the pinned producer in this sidecar's reader process. The
running source, optimizer and original labels remain unmodified.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_clamped_observation import export, REGISTERED_FACTORIES

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-connection-E273-r1')
PIN = 'a452d5b125d6d95c729ea60607524ce5d35370d5aab57ac6f2b1b7869d61991a'
PRODUCER_SHA = 'df72663a76e50866dd54b637a3ee15612364f6eba1b711443f947e05b29847da'
JOB = '28949435'


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def require_original_handle(raw):
    rows = [r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] == JOB]
    if (len(rows) != 1 or len(rows[0]) < 10 or rows[0][1] not in ('RUNNING', 'COMPLETED')
            or rows[0][2:7] != ['0:0', '32', 'hfacnormal01', 'iai806', '1']):
        raise ValueError('only original live/successful one-node32CPU handle admitted')
    return rows[0]


def completed_steps(text):
    steps = []
    for line in text.splitlines():
        columns = line.split()
        if len(columns) == 5 and columns[0] == 'CheckedFIRE:':
            step = int(columns[1])
            if not 0 <= step <= 24 or step != len(steps):
                raise ValueError('original consecutive bounded optimizer steps required')
            float(columns[3]); float(columns[4]); steps.append(step)
    if not steps:
        raise ValueError('no complete force-evaluated frame is recorded yet')
    return steps


def collect(step):
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', JOB,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True)
    accounting = require_original_handle(raw)
    work = ROOT/'path/band'
    if step not in completed_steps((work/'vcneb.opt.log').read_text()):
        raise ValueError('requested frame not yet complete; keep original job running')
    parent = ROOT/'snapshot_exports'; parent.mkdir(exist_ok=True)
    output = parent/f'step_{step:04d}'; archive = parent/f'step_{step:04d}.tar.gz'
    if output.exists() or archive.exists():
        raise FileExistsError('snapshot export exists; inspect it without overwriting')
    accepted = json.loads((ROOT/'accepted.json').read_bytes())
    if accepted['job_id'] != JOB or accepted['bundle_pin'] != PIN:
        raise ValueError('original admitted E273 source and handle required')
    spec = importlib.util.spec_from_file_location('E273_snapshot_fixed_control', ROOT/'bundle/control.py')
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control)
    control.validate(PIN)
    producer = ROOT/'bundle/run.slurm'
    if sha256(producer) != PRODUCER_SHA:
        raise ValueError('registered finite producer differs')
    # Scoped to this read-only postprocessor, never written to the runtime.
    REGISTERED_FACTORIES[PRODUCER_SHA] = 'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory'
    output.mkdir()
    observation = export(work, step, output/'observation', JOB, production_script=producer)
    fresh = []
    for record in observation['raw_image_evaluations']:
        source = Path(record['raw_source'])
        destination = output/'observation/raw'/f"image_{record['image_index']:04d}"
        (destination/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
            shutil.copyfile(source/name, destination/name)
        for name in ('call_audit.json', 'abacus.out', 'abacus.err'):
            if (source/name).is_file(): shutil.copyfile(source/name, destination/name)
        if (any(sha256(destination/n) != record['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
                or sha256(destination/'OUT.ABACUS/running_scf.log') != record['raw_log_sha256']):
            raise ValueError('selected native bytes changed during copying')
        if source.is_relative_to(work):
            audit = json.loads((source/'call_audit.json').read_bytes())
            fresh.append(dict(image=record['image_index'],source=str(source),seconds=audit['elapsed_seconds']))
    if len(fresh) != (0 if step == 0 else 7):
        raise ValueError('the selected frame must preserve endpoints and exact initial caches')
    summary = None
    if accounting[1] == 'COMPLETED':
        summary = json.loads((work/'vcneb_summary.json').read_bytes())
        if (work/'vcneb_failure.json').exists(): raise ValueError('completed scheduler job has a native failure')
        if step == max(completed_steps((work/'vcneb.opt.log').read_text())):
            if bool(summary['converged']) != observation['ordinary_residual_pass']:
                raise ValueError('last-frame ordinary residual disagrees with native terminal summary')
    receipt = dict(status='complete_evaluated_frame_not_automatically_a_final_path',job_id=JOB,step=step,
        scheduler_row_at_export=accounting,bundle_pin=PIN,producer_sha256=PRODUCER_SHA,
        sidecar_source_sha256=sha256(Path(__file__)),
        ordinary_residual_pass_at_frame=observation['ordinary_residual_pass'],
        fmax_eV_A=observation['replayed_fmax_eV_A'],sampled_peak_meV_fu=observation['sampled_forward_barrier_meV_fu'],
        selected_fresh_calls=fresh,selected_frame_native_proxy_core_hours=sum(r['seconds'] for r in fresh)*32/3600,
        earlier_calls_and_preparation_costs_additional=True,terminal_summary_if_available=summary,
        all_six_assets_checked_on_HF=True,new_DFT_calls=0,scheduler_mutations=0,
        original_labels_and_forecast_scores_unchanged=True,stationary_TS_or_barrier_error_bound_certified=False)
    save(output/'receipt.json',receipt)
    files = {p.relative_to(output).as_posix():sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'files.json',files)
    with tarfile.open(archive,'x:gz') as tar: tar.add(output,arcname=output.name)
    result = dict(archive_sha256=sha256(archive),archive_bytes=archive.stat().st_size,exported_files=len(files),
        job_id=JOB,step=step,new_DFT_calls=0,scheduler_mutations=0,whole_path_completion_claimed=False)
    result['sidecar_source_sha256']=sha256(Path(__file__))
    save(parent/f'step_{step:04d}_receipt.json',result)
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--step',type=int,required=True)
    print(json.dumps(collect(parser.parse_args().step)))
