"""Exactly four same-parameter frozen SCFs; no relaxation or continuation."""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

import numpy as np
from examples.hfo2_fixed_input_factory import CONTRACT,make_factory,read_fixed_hfo2_stru,same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results,sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_clamped_endpoints import geometry_guard
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-midpoint-factorial-E274-r1')
PARENT=ROOT.with_name('20261018-varneb-midpoint-reverse-E269-r1')
PARENT_PIN='381ae9787cca91b1435b121dafd74ee072d1aa3bd22335f05a7e7ba42bfc7ae1'
KEYS=[('preserving','AB'),('preserving','BA'),('T','AB'),('T','BA')]


def save(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(value,handle,indent=2,allow_nan=False);handle.write('\n')


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def queue():
    raw=subprocess.check_output(['squeue','-h','-u','iai806','-o','%i|%j|%T'],text=True)
    return [r for r in raw.splitlines() if r.split('|')[1].lower().startswith(('hfo2','varneb','vcneb'))]


def completed_rows(jobs):
    raw=subprocess.check_output(['sacct','-n','-X','-P','-j',','.join(jobs),
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'],text=True)
    rows=[r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] in jobs]
    if ({r[0] for r in rows}!=set(jobs) or len(rows)!=len(jobs)
            or any(len(r)<10 or r[1:7]!=['COMPLETED','0:0','32','hfacnormal01','iai806','1'] for r in rows)):
        raise ValueError('exact original completed one-node32CPU handles required')
    return rows


def validate(pin):
    bundle=ROOT/'bundle'
    if sha256(bundle/'bundle_files.json')!=pin:raise ValueError('external bundle pin differs')
    for name,digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        p=Path(name)
        if p.is_absolute() or '..' in p.parts or sha256(bundle/p)!=digest:
            raise ValueError('fixed source or input differs')
    if sha256(PARENT/'bundle/bundle_files.json')!=PARENT_PIN:raise ValueError('parent pin differs')
    files=json.loads((PARENT/'bundle/bundle_files.json').read_bytes())
    if sha256(PARENT/'bundle/control.py')!=files['control.py']:raise ValueError('parent checker differs')
    parent=load_module(PARENT/'bundle/control.py','E274_parent_check')
    runtime=parent.validate(PARENT_PIN)  # Original375files, binary and six physical assets.
    completed_rows(['28938772','28938773','28946180','28946181'])
    prep=json.loads((bundle/'prepared/preparation.json').read_bytes())
    if (prep['upper_stage']!='ordinary_terminal' or prep['physical_input_sha256']!=CONTRACT
            or prep['target_strain']!=.005 or prep['pressure_GPa']!=0 or prep['formula_units']!=4
            or prep['execution_admitted'] or prep['additional_SCFs_required']!=4
            or prep['new_DFT_calls']!=0 or prep['jobs_submitted']!=0
            or [(r['channel'],r['corner']) for r in prep['decks']]!=KEYS):
        raise ValueError('exact four terminal-source frozen crossovers required')
    for deck in prep['decks']:
        raw=bundle/'prepared'/deck['path']
        if (not raw.resolve().is_relative_to((bundle/'prepared').resolve())
                or any(sha256(raw/n)!=h for n,h in deck['input_sha256'].items())
                or any(deck['input_sha256'][n]!=CONTRACT[n] for n in ('INPUT','KPT'))):
            raise ValueError('prepared crossover input bytes differ')
    return runtime,prep


def preflight(pin):
    runtime,prep=validate(pin);out=ROOT/'preflight';out.mkdir(exist_ok=False)
    for channel in prep['channels']:
        for source in channel['sources'].values():
            raw=Path(source['source_static_directory'])
            if (any(sha256(raw/n)!=h for n,h in source['source_input_sha256'].items())
                    or sha256(raw/'OUT.ABACUS/running_scf.log')!=source['source_log_sha256']
                    or audited_results(raw)['energy']!=source['cached_energy_eV_cell']):
                raise ValueError('native AA/BB cache input, output or energy differs')
    decks=[]
    for deck in prep['decks']:
        raw=ROOT/'bundle/prepared'/deck['path'];atoms=read_fixed_hfo2_stru(raw/'STRU')
        geometry_roundtrip(raw/'STRU',atoms);geometry_guard([atoms])
        target=out/deck['path'];target.mkdir(parents=True)
        for name in ('INPUT','KPT'):shutil.copyfile(raw/name,target/name)
        # Match the immutable calculator's actual writer on this native runtime,
        # rather than conflating an analysis lift with a serialized STRU hash.
        write_structure(target/'STRU',atoms);geometry_roundtrip(target/'STRU',atoms)
        native=read_fixed_hfo2_stru(target/'STRU')
        delta=native.get_scaled_positions(wrap=False)-atoms.get_scaled_positions(wrap=False)
        shifts=np.rint(delta).astype(int)
        if not same_ordered_geometry(atoms,native) or not np.allclose(delta-shifts,0,atol=1e-12,rtol=0):
            raise ValueError('native serialization must preserve ordered periodic geometry')
        decks.append(dict(channel=deck['channel'],corner=deck['corner'],
            input_sha256={n:sha256(target/n) for n in ('INPUT','KPT','STRU')},
            prepared_STRU_sha256=deck['input_sha256']['STRU'],native_integer_coordinate_shifts=shifts.tolist()))
    result=dict(bundle_pin=pin,status='four_native_frozen_crossovers_ready',decks=decks,new_DFT_calls=0,
        exact_existing_diagonal_caches_replayed=True,physical_parameters_changed=False)
    save(ROOT/'preflight.json',result);return result


def validate_held(raw,job):
    fields=dict(re.findall(r'(\S+)=(\S+)',raw))
    expected=dict(JobId=str(job),JobState='PENDING',Reason='JobHeldUser',Partition='hfacnormal01',
        NumCPUs='32',NumTasks='1',TimeLimit='01:00:00',JobName='hfo2-E274-crossovers')
    expected['CPUs/Task']='32'
    if fields.get('NumNodes') not in ('1','1-1') or any(fields.get(k)!=v for k,v in expected.items()):
        raise ValueError('held resources must remain one node/32CPU/one hour')


def dispatch(mode,pin):
    validate(pin);proof=json.loads((ROOT/'preflight.json').read_bytes())
    if proof['bundle_pin']!=pin or proof['new_DFT_calls']!=0 or len(proof['decks'])!=4:
        raise ValueError('native four-deck preflight required')
    if mode=='submit-held':
        if len(queue())>=2 or any((ROOT/f'{n}.json').exists() for n in ('intent','held','accepted','released')):
            raise ValueError('existing intent/handle or two-study cap: inspect, do not resubmit')
        save(ROOT/'intent.json',dict(bundle_pin=pin,maximum_new_SCFs=4))
        raw=subprocess.check_output(['sbatch','--parsable','--hold','--job-name=hfo2-E274-crossovers',
            f'--export=ALL,RUN_DFT=1,BUNDLE_SHA256={pin}',f'--chdir={ROOT}',
            f'--output={ROOT}/slurm-%j.out',str(ROOT/'bundle/run.slurm')],text=True).strip()
        job=raw.split(';')[0]
        if not job.isdigit():raise ValueError('unknown submission outcome: inspect intent/accounting')
        save(ROOT/'held.json',dict(job_id=job,bundle_pin=pin,submission_response=raw))
    else:
        held=json.loads((ROOT/'held.json').read_bytes());job=held['job_id']
        if held['bundle_pin']!=pin or (ROOT/'released.json').exists() or len(queue())>2:
            raise ValueError('original held handle only')
    resource=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    validate_held(resource,job)
    if mode=='release':
        save(ROOT/'accepted.json',dict(job_id=job,bundle_pin=pin,held_resource_text=resource))
        subprocess.check_call(['scontrol','release',job])
        save(ROOT/'released.json',dict(job_id=job,bundle_pin=pin))
    return dict(job_id=job,bundle_pin=pin,mode=mode,resource_text=resource)


def require_resources(environment,job):
    expected=dict(RUN_DFT='1',SLURM_JOB_ID=job,SLURM_JOB_NAME='hfo2-E274-crossovers',
        SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',SLURM_NTASKS='1',
        SLURM_CPUS_PER_TASK='32',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if any(environment.get(k)!=v for k,v in expected.items()):
        raise ValueError('accepted one-node32CPU compute allocation required; never run on login')


def run(pin):
    runtime,prep=validate(pin);accepted=json.loads((ROOT/'accepted.json').read_bytes())
    proof=json.loads((ROOT/'preflight.json').read_bytes())
    require_resources(os.environ,accepted['job_id'])
    if (accepted['bundle_pin']!=pin or proof['bundle_pin']!=pin or proof['new_DFT_calls']!=0
            or len(proof['decks'])!=4 or len(queue())>2):
        raise ValueError('only the exact accepted study handle may run')
    output=ROOT/'results';output.mkdir(exist_ok=False);records=[]
    factory=make_factory(parameters={'source_directory':str(runtime.PHYSICS)},command=f'mpirun -np 32 {runtime.BIN}')
    for deck in prep['decks']:
        atoms=read_fixed_hfo2_stru(ROOT/'bundle/prepared'/deck['path']/'STRU')
        native_proof=next(r for r in proof['decks'] if (r['channel'],r['corner'])==(deck['channel'],deck['corner']))
        folder=output/deck['path'];atoms.calc=factory(0,atoms,folder)
        atoms.get_potential_energy();atoms.get_forces();atoms.get_stress()
        calls=list(folder.glob('scf_*/call_audit.json'))
        if len(calls)!=1:raise ValueError('exactly one frozen SCF per crossover required')
        raw=calls[0].parent;call=json.loads(calls[0].read_bytes())
        if (call['input_sha256']['STRU']!=native_proof['input_sha256']['STRU']
                or not same_ordered_geometry(atoms,read_fixed_hfo2_stru(raw/'STRU'))):
            raise ValueError('frozen ordered crossover geometry changed')
        record=dict(channel=deck['channel'],corner=deck['corner'],raw_source=str(raw),
            input_sha256=call['input_sha256'],raw_log_sha256=call['raw_log_sha256'],results=call['results'],
            seconds=call['elapsed_seconds']);records.append(record);print(json.dumps(record),flush=True)
    result=dict(job_id=accepted['job_id'],bundle_pin=pin,new_SCF_calls=4,records=records,
        native_proxy_core_hours=sum(r['seconds'] for r in records)*32/3600,
        relaxation_or_new_NEB_performed=False,old_labels_and_forecasts_unchanged=True)
    save(output/'complete.json',result);return result


def collect(pin):
    _,prep=validate(pin);accepted=json.loads((ROOT/'accepted.json').read_bytes())
    row=completed_rows([accepted['job_id']])[0]
    complete=json.loads((ROOT/'results/complete.json').read_bytes())
    if (complete['job_id']!=accepted['job_id'] or complete['bundle_pin']!=pin or accepted['bundle_pin']!=pin
            or complete['new_SCF_calls']!=4 or [(r['channel'],r['corner']) for r in complete['records']]!=KEYS):
        raise ValueError('complete original four-call terminal required')
    if len(list((ROOT/'results').glob('*/*/scf_*/call_audit.json')))!=4:
        raise ValueError('exactly four actual native crossover calls required')
    output=ROOT/'native_terminal';output.mkdir(exist_ok=False)
    for record in complete['records']:
        raw=Path(record['raw_source']);parsed=audited_results(raw)
        call=json.loads((raw/'call_audit.json').read_bytes())
        if (any(sha256(raw/n)!=h for n,h in record['input_sha256'].items())
                or any(record['input_sha256'][n]!=h for n,h in CONTRACT.items())
                or sha256(raw/'OUT.ABACUS/running_scf.log')!=record['raw_log_sha256']
                or any(not np.allclose(parsed[k],record['results'][k],atol=1e-12,rtol=0) for k in parsed)
                or not np.isfinite(record['seconds']) or record['seconds']<=0
                or record['seconds']!=call['elapsed_seconds']):
            raise ValueError('complete native crossover E/F/stress or six input bytes differ')
        target=output/record['channel']/record['corner'];(target/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT','KPT','STRU','call_audit.json','input_sha256.json','OUT.ABACUS/running_scf.log','abacus.out','abacus.err'):
            shutil.copyfile(raw/name,target/name)
    save(output/'complete.json',complete);save(output/'preparation.json',prep)
    shutil.copyfile(ROOT/'preflight.json',output/'preflight.json')
    save(output/'scheduler.json',dict(row=row,allocated_core_hours=float(row[7])*32/3600))
    files={p.relative_to(output).as_posix():sha256(p) for p in sorted(output.rglob('*')) if p.is_file()};save(output/'files.json',files)
    archive=ROOT/'native_terminal.tar.gz'
    with tarfile.open(archive,'x:gz') as tar:tar.add(output,arcname='native_terminal')
    receipt=dict(archive_sha256=sha256(archive),archive_bytes=archive.stat().st_size,exported_files=len(files),
        job_id=accepted['job_id'],bundle_pin=pin,control_source_sha256=sha256(Path(__file__)),
        native_preflight_sha256=sha256(ROOT/'preflight.json'),new_SCF_calls=4,all_native_EFS_six_assets_replayed_on_HF=True,
        allocated_core_hours=float(row[7])*32/3600,native_proxy_core_hours=complete['native_proxy_core_hours'],
        new_DFT_calls_for_analysis=0,scheduler_mutations=0)
    save(ROOT/'native_terminal_receipt.json',receipt);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('preflight','submit-held','release','run','collect'));p.add_argument('--bundle-pin',required=True)
    a=p.parse_args()
    result=(dispatch(a.mode,a.bundle_pin) if a.mode in ('submit-held','release')
        else {'preflight':preflight,'run':run,'collect':collect}[a.mode](a.bundle_pin))
    print(json.dumps(result,allow_nan=False))
