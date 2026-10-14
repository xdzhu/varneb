"""One reviewed child40 continuation; never automatic rolling submission."""
import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261014-varneb-clamped-reversing-continue-E146-r1'
PARENT = '28764704'
ORIGINAL_PARENT = '28709791'
CHANNEL = 'PO_flip_T_pattern_reversing'
ORIGINAL = RUNS/'20261012-varneb-clamped-flip-continue-E103-r1'
TERMINAL = ORIGINAL/ORIGINAL_PARENT/'child_terminal'
AUDIT_SHA = '6cf7e71447da77aa9f869e7590f70368f4510cb60326017d728c277c06f3c8c1'
SEED_SHA = 'c176730bcc50a9782ecfa8483e97c949b59b2c83dbcb43b284304df8ff342cdc'
LEGACY = ORIGINAL/'continue_flips_HF.py'
LEGACY_SHA = '6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7'
SOURCE = RUNS/'20261010-varneb-nested-controls-E055-r1/source'
SCRIPT = ROOT/'hf_hfo2_clamped_continue_E146_20261014.slurm'
SCRIPT_SHA = 'd79f538fd0381349d7149cd277ec5756162847e4d9c673038ba4de25ad3f6f9c'
LOCK = RUNS/'20261011-varneb-clamped-held-flips-E100-r1/admission.lock'
EVENTS = ['one_held_submission_intent','one_accepted_handle','one_release_intent',
          'one_existing_handle_released','released_slot_independently_corroborated']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


def legacy():
    if sha(LEGACY) != LEGACY_SHA:
        raise ValueError('unchanged pinned E103 reader required')
    spec = importlib.util.spec_from_file_location('E146_reused_E103', LEGACY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def paths():
    return dict(work=ROOT/CHANNEL/'band', stdout=ROOT/(CHANNEL+'.slurm.out'),
                stderr=ROOT/(CHANNEL+'.slurm.err'))


def environment():
    return dict(RUN_DFT='1', SOURCE_ROOT=str(SOURCE), SEED_ROOT=str(TERMINAL/'seed'),
        WORKDIR=str(paths()['work']), PARENT_JOB_ID=PARENT,
        PARENT_AUDIT_PATH=str(TERMINAL/'audit_receipt.json'), PARENT_AUDIT_SHA256=AUDIT_SHA,
        SEED_MANIFEST_SHA256=SEED_SHA)


def verify_parent():
    """Require the actual cap and exact ancestry before replaying startup checks."""
    old = legacy()
    if old.accepted(ORIGINAL_PARENT)['job_id'] != PARENT:
        raise ValueError('original accepted child differs')
    body = run(['sacct','-X','-n','-P','-j',PARENT,
        '--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000'])
    old.identity(body, PARENT, ORIGINAL_PARENT, child=True)
    if (sha(TERMINAL/'audit_receipt.json') != AUDIT_SHA
            or sha(TERMINAL/'seed/manifest.json') != SEED_SHA or sha(SCRIPT) != SCRIPT_SHA
            or sha(old.BINARY) != old.BINARY_SHA or paths()['work'].exists()):
        raise ValueError('reviewed parent/seed/wrapper/binary or fresh namespace differs')
    receipt = json.loads((TERMINAL/'audit_receipt.json').read_text())
    old_work = old.paths(ORIGINAL_PARENT)['work']
    if sha(old_work/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
        raise ValueError('terminal source summary changed')
    _, _, support = old.load_auditor()
    for key, name in [('executed_auditor_sha256','audit_hfo2_clamped_terminal.py'),
                      ('executed_exporter_sha256','export_hfo2_clamped_observation.py'),
                      ('executed_residual_analyzer_sha256','analyze_hfo2_clamped_residual.py')]:
        if receipt[key] != support[name]:
            raise ValueError('original audit support differs')
    body = SCRIPT.read_text().split(" 'import hashlib,json,pathlib,sys\n",1)[1]
    body = 'import hashlib,json,pathlib,sys\n'+body.split("' "+chr(34)+"\u0024{seed_root}/manifest.json"+chr(34),1)[0]
    before = sys.argv
    try:
        sys.argv = ['reviewed-startup',str(TERMINAL/'seed/manifest.json'),SEED_SHA,
                    str(TERMINAL/'audit_receipt.json'),AUDIT_SHA,PARENT,str(SOURCE)]
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(body,str(SCRIPT),'exec'),{})
    finally:
        sys.argv = before
    return receipt


def accepted():
    result = json.loads((ROOT/'dispatch_receipt.json').read_text())
    journal = ROOT/'dispatch_journal.jsonl'
    records = [json.loads(line) for line in journal.read_text().splitlines()]
    if (result.get('parent_job') != PARENT or result.get('source_step') != 40
            or result.get('status') != 'one_same_chain_continuation_submitted_released_and_slot_corroborated'
            or not str(result.get('job_id','')).isdecimal()
            or result.get('total_submissions') != 1 or result.get('release_writes') != 1
            or result.get('new_independent_chains') != 0 or result.get('script_sha256') != SCRIPT_SHA
            or result.get('executed_helper_sha256') != sha(Path(__file__))
            or result.get('reused_E103_driver_sha256') != LEGACY_SHA
            or result.get('journal_sha256') != sha(journal)
            or [r['event'] for r in records] != EVENTS
            or any(r.get('job_id') != result['job_id'] for r in records[1:])):
        raise ValueError('accepted E146 evidence differs; never reset/replay')
    return result


def slots(body, *, spare, extra=None):
    known = dict(extra or {})
    if (ROOT/'dispatch_receipt.json').exists():
        child = accepted()['job_id']
        known[child] = 'hfo2-G2-E146-'+PARENT
    return legacy().slots(body, ORIGINAL_PARENT, spare=spare, extra=known)


def dispatch():
    journal_path = ROOT/'dispatch_journal.jsonl'
    if journal_path.exists() or (ROOT/'dispatch_receipt.json').exists() or LOCK.exists():
        raise FileExistsError('existing intent or lock; reconcile, never resubmit')
    verify_parent()
    accounting = run(['sacct','-X','-n','-P','-j',PARENT,
                      '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if accounting != PARENT+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('actual reviewed successful cap required')
    queue = ['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C|%R']
    slots(run(queue), spare=True)
    run(['bash','-n',str(SCRIPT)])
    name = 'hfo2-G2-E146-'+PARENT
    p = paths()
    argv = ['sbatch','--parsable','--hold','--job-name='+name,
        '--export=ALL,'+','.join(k+'='+v for k,v in environment().items()),
        '--output='+str(p['stdout']),'--error='+str(p['stderr']),str(SCRIPT)]
    with LOCK.open('x') as stream:
        stream.write(json.dumps(dict(parent_job=PARENT,helper_sha256=sha(Path(__file__)))))
        stream.flush(); os.fsync(stream.fileno())
    with journal_path.open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event,check_UTC=datetime.now(timezone.utc).isoformat(),**data))+'\n')
            journal.flush(); os.fsync(journal.fileno())
        try:
            committed = slots(run(queue),spare=True)
            record(EVENTS[0],argv=argv,parent_accounting=accounting,
                   audit_receipt_sha256=AUDIT_SHA,committed_slots=committed)
            handle = run(argv)
            if not handle.isdecimal():
                raise ValueError('ambiguous submission reply; never repeat')
            record(EVENTS[1],job_id=handle)
            row = dict(t.split('=',1) for t in shlex.split(run(['scontrol','show','job','-o',handle])) if '=' in t)
            expected = dict(JobId=handle,JobName=name,JobState='PENDING',Reason='JobHeldUser',
                Priority='0',Partition='hfacnormal01',NumCPUs='32',NumTasks='1',Command=str(SCRIPT),
                Dependency='(null)',StdOut=str(p['stdout']),StdErr=str(p['stderr']),
                TimeLimit='12:00:00',Restarts='0',**{'CPUs/Task':'32'})
            def same_identity(native):
                return (all(native.get(k)==v for k,v in expected.items() if k not in ('JobState','Reason','Priority'))
                        and native.get('NumNodes') in ('1','1-1')
                        and native.get('UserId','').startswith('iai806('))
            if not same_identity(row) or any(row.get(k)!=expected[k] for k in ('JobState','Reason','Priority')):
                raise ValueError('exact held native allocation not corroborated')
            slots(run(queue),spare=True,extra={handle:name})
            record(EVENTS[2],job_id=handle,accepted_native_state=row)
            reply = run(['scontrol','release',handle])
            record(EVENTS[3],job_id=handle,reply=reply)
            after = dict(t.split('=',1) for t in shlex.split(run(['scontrol','show','job','-o',handle])) if '=' in t)
            if (not same_identity(after) or after.get('JobState') not in ('PENDING','RUNNING','COMPLETING','CONFIGURING')
                    or after.get('Reason')=='JobHeldUser' or after.get('Priority')=='0'):
                raise ValueError('release not independently corroborated; retain lock')
            committed = slots(run(queue),spare=False,extra={handle:name})
            if handle not in committed:
                raise ValueError('released job not a committed slot; retain lock')
            record(EVENTS[4],job_id=handle,native_control=after,committed_slots=committed)
        except Exception as exc:
            record('reconcile_never_repeat_write',error=type(exc).__name__,message=str(exc))
            raise
    result = dict(status='one_same_chain_continuation_submitted_released_and_slot_corroborated',
        job_id=handle,parent_job=PARENT,source_step=40,max_new_updates=40,wall_hours=12,
        total_submissions=1,release_writes=1,dependency_writes=0,new_independent_chains=0,
        audit_receipt_sha256=AUDIT_SHA,seed_manifest_sha256=SEED_SHA,script_sha256=SCRIPT_SHA,
        executed_helper_sha256=sha(Path(__file__)),reused_E103_driver_sha256=LEGACY_SHA,
        journal_sha256=sha(journal_path),physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,FIRE_optimizer_state_restored=False,execution_started_claimed=False)
    with (ROOT/'dispatch_receipt.json').open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    LOCK.unlink()
    return result


def audit_child():
    child = accepted()['job_id']
    output = ROOT/'child_terminal'
    if output.exists():
        raise FileExistsError('inspect existing child audit; never overwrite')
    native = run(['sacct','-X','-n','-P','-j',child,
        '--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000'])
    fields = native.split('|'); p = paths()
    if (len(fields)!=9 or fields[:6]!=[child,'COMPLETED','iai806','32','hfacnormal01','hfo2-G2-E146-'+PARENT]
            or fields[6:8]!=[str(p['stdout']),str(p['stderr'])]):
        raise ValueError('actual exact successful native terminal identity required')
    argv = shlex.split(fields[8]); exports = [v[len('--export='):] for v in argv if v.startswith('--export=')]
    if (not argv or argv[0]!='sbatch' or argv[-1]!=str(SCRIPT) or len(exports)!=1
            or exports[0].split(',')[0]!='ALL' or len(exports[0].split(','))!=len(environment())+1
            or dict(v.split('=',1) for v in exports[0].split(',')[1:])!=environment()):
        raise ValueError('actual submitted ancestry differs')
    old = legacy(); auditor, registration, support = old.load_auditor()
    if sha(SCRIPT)!=SCRIPT_SHA:
        raise ValueError('reviewed wrapper changed')
    from scripts.export_hfo2_clamped_observation import REGISTERED_FACTORIES
    REGISTERED_FACTORIES[SCRIPT_SHA] = 'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory'
    summary = json.loads((p['work']/'vcneb_summary.json').read_text())
    step = auditor.completed_step(p['work']/'vcneb.opt.log')
    capped = summary.get('status')=='max_steps_reached'
    if not 0<=step<=40 or (capped and step!=40):
        raise ValueError('unexpected finite segment history')
    receipt = auditor.audit(p['work'],child,SCRIPT,SOURCE,registration['archive'],registration['archive_sha'],
                            output,outcome='step-cap' if capped else 'converged',seed_latest=capped)
    if len(receipt['runtime_code_sha256'])!=338:
        raise ValueError('original runtime count differs')
    with (output/'driver_receipt.json').open('x') as stream:
        json.dump(dict(job_id=child,parent_job=PARENT,executed_driver_sha256=sha(Path(__file__)),
            reused_E103_driver_sha256=LEGACY_SHA,terminal_audit_receipt_sha256=sha(output/'audit_receipt.json'),
            support_sha256=support,new_DFT_calls=0,scheduler_writes=0,new_independent_chains=0),stream,indent=2)
    return dict(status=receipt['status'],job_id=child,terminal_step=step,
                terminal_audit_sha256=sha(output/'audit_receipt.json'))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('preflight','dispatch','audit-child'))
    action = parser.parse_args().action
    result = (dispatch() if action=='dispatch' else audit_child() if action=='audit-child' else
              dict(status='parent_ancestry_and_startup_ready_not_submitted',terminal_step=verify_parent()['terminal_step']))
    print(json.dumps(result,indent=2,allow_nan=False))
