"""One audited same-chain held submission/release; never rewrite dependencies."""
from datetime import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-flip-E069-r1'
SOURCE = RUNS/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
PARENT = '28722810'
OTHER = '28692775'
WAIT = '28692776'
ALLOWED = {PARENT,OTHER,WAIT,'28709788','28709789','28709790','28709791'}
OLD_WORK = RUNS/'20261011-varneb-clamped-flip-E067-r1/PO_flip_T_pattern_preserving/band'
SCRIPT = ROOT/'hf_hfo2_clamped_flip_continue_E069_20261011.slurm'
SCRIPT_SHA = '61199684c543a016a546f02cbfe1eb7ee659c6310db4f7ad0c0340e714647379'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
RULE_SOURCE = RUNS/'20261011-varneb-clamped-flip-E067-r1/submit_held_once.py'
RULE_SHA = 'e74ba9ba9090ec0a46cb56d046c4931f7889e4a9e4855a38fe59544b12708927'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv,check=True,capture_output=True,text=True,timeout=45).stdout.strip()


def rules():
    if sha(RULE_SOURCE) != RULE_SHA:
        raise ValueError('immutable reviewed identity parser changed')
    spec = importlib.util.spec_from_file_location('E069_read_only_rules',RULE_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spare_slot(body, *, accepted=None):
    active = set()
    for line in body.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed queue')
        job,name,state,partition,cpus = parts
        if not name.startswith('hfo2'):
            continue
        if job not in ALLOWED|({accepted} if accepted else set()) or partition != 'hfacnormal01' or cpus != '32':
            raise ValueError('unexpected study allocation')
        if state in ('RUNNING','COMPLETING','CONFIGURING'):
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    if active != {OTHER}:
        raise ValueError('exactly registered PO-to-T must occupy the first slot')
    return sorted(active)


def verify_audit(reviewed_sha):
    if sha(ROOT/'audit_receipt.json') != reviewed_sha:
        raise ValueError('reviewed material receipt changed')
    receipt = json.loads((ROOT/'audit_receipt.json').read_text())
    if (receipt['status'] != 'audited_step_cap_not_converged_same_chain_seed'
            or receipt['job_id'] != PARENT or receipt['ordinary_converged'] is not False
            or receipt['terminal_step'] != 20 or receipt['fresh_interior_SCFs'] != 140
            or receipt['current_frame_exact_caches_checked'] != 9
            or receipt['new_DFT_calls'] != 0 or receipt['scheduler_mutations'] != 0
            or not receipt['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
            or receipt['physical_inputs_or_running_source_changed'] is not False
            or receipt['holdout_generated_or_read'] is not False
            or receipt['terminal_summary_sha256'] != sha(OLD_WORK/'vcneb_summary.json')
            or receipt['executed_auditor_sha256'] != sha(ROOT/'audit_source/scripts/audit_hfo2_clamped_terminal.py')
            or receipt['executed_exporter_sha256'] != sha(ROOT/'audit_source/scripts/export_hfo2_clamped_observation.py')
            or receipt['executed_residual_analyzer_sha256'] != sha(ROOT/'audit_source/scripts/analyze_hfo2_clamped_residual.py')
            or len(receipt['runtime_code_sha256']) != 334
            or any(sha(SOURCE/n) != h for n,h in receipt['runtime_code_sha256'].items())
            or sha(SCRIPT) != SCRIPT_SHA or sha(BINARY) != BINARY_SHA
            or (ROOT/'PO_flip_T_pattern_preserving/band').exists()):
        raise ValueError('actual audited material/runtime/script namespace differs')
    seed = ROOT/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    required = dict(status='audited_G2_geometry_continuation_prepared',source_job_id=PARENT,source_step=20,
        current_frame_exact_caches=9,new_independent_chains=0,new_DFT_calls=0,physical_inputs_changed=False,
        holdout_generated_or_read=False,optimizer_state_restored=False,pressure_GPa=0,
        ordinary_fmax_eV_A=.10,climb=False,n_total_images=9,n_active_images=7,n_fixed_endpoints=2)
    names = {'initial.vasp','final.vasp','seed.traj','substrate.vasp','factory_parameters.json'}
    if (sha(seed/'manifest.json') != receipt['seed_manifest_sha256']
            or any(manifest.get(k) != v for k,v in required.items())
            or set(manifest['files_sha256']) != names
            or any(sha(seed/n) != h for n,h in manifest['files_sha256'].items())):
        raise ValueError('exact latest complete same-chain seed changed')
    return receipt


def dispatch_once(reviewed_sha):
    journal_path = ROOT/'dispatch_journal.jsonl'
    if journal_path.exists() or (ROOT/'dispatch_receipt.json').exists():
        raise FileExistsError('existing intent/handle; reconcile without resubmission or release retry')
    receipt = verify_audit(reviewed_sha)
    parser = rules()
    accounting = run(['sacct','-X','-n','-P','-j',PARENT,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if accounting != PARENT+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('registered capped parent must remain terminal')
    waiter = parser.held_waiter(run(['scontrol','show','job','-o',WAIT]),WAIT,required=set(),allowed=set())
    queue = ['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']
    active = spare_slot(run(queue))
    run(['bash','-n',str(SCRIPT)])
    with journal_path.open('x') as journal:
        def record(event,**data):
            journal.write(json.dumps(dict(event=event,check_CST=datetime.now().isoformat(),**data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        argv = ['sbatch','--parsable','--hold','--job-name=hfo2-G2-flip-E069',
            '--export=ALL,RUN_DFT=1,SEED_ROOT='+str(ROOT/'seed')+',WORKDIR='+str(ROOT/'PO_flip_T_pattern_preserving/band')
            +',SEED_MANIFEST_SHA256='+receipt['seed_manifest_sha256'],
            '--output='+str(ROOT/'flip.slurm.out'),'--error='+str(ROOT/'flip.slurm.err'),str(SCRIPT)]
        record('one_held_submission_intent',argv=argv,audit_receipt_sha256=reviewed_sha,
               parent_accounting=accounting,active_jobs=active,unchanged_held_successor=waiter)
        try:
            handle = run(argv)
            if not handle.isdecimal():
                raise ValueError('ambiguous accepted handle; never resubmit')
            record('one_accepted_handle',job_id=handle)
            row = parser.fields(run(['scontrol','show','job','-o',handle]))
            expected = dict(JobId=handle,JobName='hfo2-G2-flip-E069',JobState='PENDING',Reason='JobHeldUser',
                Priority='0',Partition='hfacnormal01',NumCPUs='32',Command=str(SCRIPT),
                StdOut=str(ROOT/'flip.slurm.out'),StdErr=str(ROOT/'flip.slurm.err'),Dependency='(null)')
            if (any(row.get(k) != v for k,v in expected.items()) or not row.get('UserId','').startswith('iai806(')):
                raise ValueError('exact accepted held continuation not corroborated')
            spare_slot(run(queue),accepted=handle)
            parser.held_waiter(run(['scontrol','show','job','-o',WAIT]),WAIT,required=set(),allowed=set())
            record('one_release_intent',job_id=handle,accepted_native_state=row)
            reply = run(['scontrol','release',handle])
            record('one_existing_handle_released',job_id=handle,reply=reply)
            result = dict(status='one_same_chain_continuation_submitted_and_released_not_start_proof',
                check_CST=datetime.now().isoformat(),job_id=handle,parent_job=PARENT,
                source_step=20,cumulative_source_updates=50,max_new_updates=40,
                max_new_fresh_SCFs_without_geometry_retries=280,total_submissions=1,release_writes=1,
                dependency_writes=0,new_independent_chains=0,audit_receipt_sha256=reviewed_sha,
                seed_manifest_sha256=receipt['seed_manifest_sha256'],script_sha256=sha(SCRIPT),
                executed_helper_sha256=sha(Path(__file__)),journal_sha256=sha(journal_path),
                physical_inputs_running_source_or_holdout_changed=False,
                other_successor_still_held=WAIT,execution_started_claimed=False)
            with (ROOT/'dispatch_receipt.json').open('x') as stream:
                json.dump(result,stream,indent=2,allow_nan=False)
        except Exception as exc:
            record('reconcile_existing_intent_never_repeat_writes',exception=type(exc).__name__,message=str(exc))
            raise
    print(json.dumps(result))
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-sha256',required=True)
    dispatch_once(parser.parse_args().audit_sha256)
