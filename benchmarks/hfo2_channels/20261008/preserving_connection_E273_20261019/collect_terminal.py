"""Collect the original bounded E273 terminal, including every fresh SCF.

No scheduler writes, DFT calls, source replacement or automatic continuation.
Ordinary convergence and an exhausted step cap remain different outcomes.
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

import numpy as np
from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_clamped_terminal import completed_step,fresh_interior_call_count,validate_terminal
from scripts.audit_hfo2_static_replica import audited_results,sha256
from scripts.export_hfo2_clamped_observation import export,REGISTERED_FACTORIES

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-connection-E273-r1')
PIN='a452d5b125d6d95c729ea60607524ce5d35370d5aab57ac6f2b1b7869d61991a'
PRODUCER_SHA='df72663a76e50866dd54b637a3ee15612364f6eba1b711443f947e05b29847da'
JOB='28949435'


def save(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(value,handle,indent=2,allow_nan=False);handle.write('\n')


def terminal_row(raw):
    rows=[r.split('|') for r in raw.strip().splitlines() if r.split('|')[0]==JOB]
    if (len(rows)!=1 or len(rows[0])<10
            or rows[0][1:7]!=['COMPLETED','0:0','32','hfacnormal01','iai806','1']):
        raise ValueError('exact completed original E273 one-node32CPU job required')
    return rows[0]


def admit(row,summary,step,failure_exists,preflight):
    REGISTERED_FACTORIES[PRODUCER_SHA]='examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory'
    validate_terminal(JOB,'|'.join(row[:6]),summary,failure_exists,
                      outcome='converged' if summary.get('converged') else 'step-cap')
    if not 0<=step<=24 or not summary['converged'] and step!=24:
        raise ValueError('ordinary convergence or the original24-update bound required')
    initial,expected=fresh_interior_call_count(step,preflight,PRODUCER_SHA)
    if initial!=0 or expected!=7*step:
        raise ValueError('all nine exact seed caches, no new endpoint evaluations required')
    return expected


def collect():
    raw=subprocess.check_output(['sacct','-n','-X','-P','-j',JOB,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'],text=True)
    row=terminal_row(raw);work=ROOT/'path/band'
    summary=json.loads((work/'vcneb_summary.json').read_bytes());step=completed_step(work/'vcneb.opt.log')
    preflight=json.loads((work/'vcneb_preflight.json').read_bytes())
    expected=admit(row,summary,step,(work/'vcneb_failure.json').exists(),preflight)
    accepted=json.loads((ROOT/'accepted.json').read_bytes())
    if accepted['job_id']!=JOB or accepted['bundle_pin']!=PIN:
        raise ValueError('only the original accepted source and handle admitted')
    spec=importlib.util.spec_from_file_location('E273_terminal_fixed_runtime',ROOT/'bundle/control.py')
    control=importlib.util.module_from_spec(spec);spec.loader.exec_module(control);control.validate(PIN)
    producer=ROOT/'bundle/run.slurm'
    if sha256(producer)!=PRODUCER_SHA:raise ValueError('immutable producer changed')
    output=ROOT/'native_terminal';archive=ROOT/'native_terminal.tar.gz'
    if output.exists() or archive.exists():raise FileExistsError('terminal export exists; inspect it without overwriting')
    calls=[]
    for path in sorted(work.glob('image_*/scf_*/call_audit.json')):
        folder=path.parent;image=folder.parent.name
        if image not in [f'image_{i:04d}' for i in range(1,8)]:
            raise ValueError('fixed endpoints must never be recalculated')
        call=json.loads(path.read_bytes())
        if (any(call['input_sha256'][n]!=h or sha256(folder/n)!=h for n,h in CONTRACT.items())
                or sha256(folder/'STRU')!=call['input_sha256']['STRU']
                or sha256(folder/'OUT.ABACUS/running_scf.log')!=call['raw_log_sha256']):
            raise ValueError('immutable native physical input or raw output differs')
        parsed=audited_results(folder)
        if any(not np.allclose(parsed[k],call['results'][k],atol=1e-12,rtol=0) for k in parsed):
            raise ValueError('native SCF energy/forces/stress do not replay')
        seconds=call['elapsed_seconds']
        if not np.isfinite(seconds) or seconds<=0:raise ValueError('native SCF timing must be positive')
        calls.append(dict(directory=str(folder),relative_directory=folder.relative_to(work).as_posix(),
            call_audit_sha256=sha256(path),input_sha256=call['input_sha256'],raw_log_sha256=call['raw_log_sha256'],seconds=seconds))
    if len(calls)!=expected:raise ValueError('all fresh interior calls must match the native cache policy')
    output.mkdir()
    observation=export(work,step,output/'observation',JOB,production_script=producer)
    if observation['ordinary_residual_pass']!=summary['converged']:
        raise ValueError('independent native force replay and terminal outcome disagree')
    selected=[]
    for record in observation['raw_image_evaluations']:
        source=Path(record['raw_source']);destination=output/'observation/raw'/f"image_{record['image_index']:04d}"
        (destination/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT','KPT','STRU','OUT.ABACUS/running_scf.log','call_audit.json','abacus.out','abacus.err'):
            if (source/name).is_file():shutil.copyfile(source/name,destination/name)
        if source.is_relative_to(work):
            call=json.loads((source/'call_audit.json').read_bytes())
            selected.append(dict(image=record['image_index'],source=str(source),seconds=call['elapsed_seconds']))
    for call in calls:
        source=Path(call['directory']);destination=output/'all_calls'/call['relative_directory']
        (destination/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT','KPT','STRU','call_audit.json','OUT.ABACUS/running_scf.log','abacus.out','abacus.err'):
            if (source/name).is_file():shutil.copyfile(source/name,destination/name)
    for name in ('vcneb_summary.json','vcneb.opt.log'):shutil.copyfile(work/name,output/name)
    receipt=dict(status='ordinary_converged_terminal' if summary['converged'] else 'bounded_unconverged_terminal',
        job_id=JOB,bundle_pin=PIN,producer_sha256=PRODUCER_SHA,sidecar_source_sha256=sha256(Path(__file__)),
        scheduler_row_at_export=row,step=step,ordinary_converged=summary['converged'],new_SCF_calls=len(calls),calls=calls,
        native_proxy_core_hours=sum(c['seconds'] for c in calls)*32/3600,allocated_core_hours=float(row[7])*32/3600,
        selected_fresh_calls=selected,selected_frame_native_proxy_core_hours=sum(c['seconds'] for c in selected)*32/3600,
        all_new_native_calls_six_assets_MPI_EFS_replayed_on_HF=True,original_preparation_costs_additional=True,
        fmax_eV_A=observation['replayed_fmax_eV_A'],sampled_peak_meV_fu=observation['sampled_forward_barrier_meV_fu'],
        old_labels_and_forecast_scores_unchanged=True,stationary_TS_or_sampling_bound_certified=False,
        new_DFT_calls=0,scheduler_mutations=0)
    save(output/'receipt.json',receipt)
    files={p.relative_to(output).as_posix():sha256(p) for p in sorted(output.rglob('*')) if p.is_file()};save(output/'files.json',files)
    with tarfile.open(archive,'x:gz') as tar:tar.add(output,arcname='native_terminal')
    result=dict(archive_sha256=sha256(archive),archive_bytes=archive.stat().st_size,exported_files=len(files),
        job_id=JOB,step=step,ordinary_converged=summary['converged'],all_runtime_calls_exported=True,
        new_SCF_calls=len(calls),native_proxy_core_hours=receipt['native_proxy_core_hours'],
        allocated_core_hours=receipt['allocated_core_hours'],new_DFT_calls=0,scheduler_mutations=0)
    save(ROOT/'native_terminal_receipt.json',result)
    return result


if __name__=='__main__':print(json.dumps(collect()))
