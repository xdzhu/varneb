"""Mock-only one-slot dispatch checks; no DFT or scheduler interaction."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

REPO = Path(__file__).resolve().parents[1]
PATH = REPO/'benchmarks/hfo2_channels/20261008/clamped_reversing_terminal_E083_20261011/release_existing_plus1_once.py'
spec = importlib.util.spec_from_file_location('E083_test', PATH)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
spec = importlib.util.spec_from_file_location('E067_fields_only',
    REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/submit_held_once.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


def report():
    return dict(job_id=helper.DONE, ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7, terminal_step=5, fresh_interior_SCFs=42,
        terminal_observation=dict(ordinary_residual_pass=True, replayed_fmax_eV_A=.092))


def queue():
    return '\n'.join([helper.LIVE+'|hfo2-G2-flip-E069|RUNNING|hfacnormal01|32',
        *[job+'|hfo2-G2-E060-'+channel+'|PENDING|hfacnormal01|32'
          for job, channel in helper.CHANNELS.items()]])


def held_row(job, pilot):
    channel = helper.CHANNELS[job]
    dep = '(null)' if job in {'28709788', '28709789'} else 'afterok:28709788(unfulfilled),afterok:28709789(unfulfilled)'
    return dict(JobId=job, JobName='hfo2-G2-E060-'+channel,
        JobState='PENDING', Priority='0', Reason='JobHeldUser',
        UserId='iai806(16284)', Partition='hfacnormal01', NumCPUs='32',
        NumTasks='1', **{'CPUs/Task': '32'}, Command=pilot,
        StdOut=str(helper.WAIT/(channel+'.slurm.out')),
        StdErr=str(helper.WAIT/(channel+'.slurm.err')), TimeLimit='04:00:00', Dependency=dep)


def row_body(row):
    return ' '.join(key+'='+value for key, value in row.items())


def test_positive_native_parser_and_spare_slot():
    helper.verified_audit(report())
    assert helper.one_used_slot(queue()+'\n999|pc-unrelated|RUNNING|other|128') == [helper.LIVE]
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    for job in helper.CHANNELS:
        helper.held_training(row_body(held_row(job, rules.PILOT)), job, rules)


@pytest.mark.parametrize('field,value', [('ordinary_converged',False), ('job_id','123'),
    ('new_DFT_calls',1), ('scheduler_mutations',1), ('holdout_generated_or_read',True),
    ('initial_fresh_interior_SCFs',0), ('fresh_interior_SCFs',35), ('terminal_step',True)])
def test_missing_or_miscounted_starter_evidence_rejected(field, value):
    with pytest.raises(ValueError): helper.verified_audit({**report(),field:value})


@pytest.mark.parametrize('force', [.1001, True, None, float('nan'), -1])
def test_threshold_not_rounded_or_boolean(force):
    value = report()
    value['terminal_observation']['replayed_fmax_eV_A'] = force
    with pytest.raises(ValueError): helper.verified_audit(value)


@pytest.mark.parametrize('change', ['missing_live','extra_active','unknown','missing_waiter','duplicate'])
def test_native_capacity_rejects_extra_or_missing_jobs(change):
    body = queue()
    if change == 'missing_live': body = '\n'.join(body.splitlines()[1:])
    if change == 'extra_active': body = body.replace('PO_to_M|PENDING','PO_to_M|RUNNING')
    if change == 'unknown': body += '\n999|hfo2-unknown|RUNNING|hfacnormal01|32'
    if change == 'missing_waiter': body = '\n'.join(body.splitlines()[:-1])
    if change == 'duplicate': body += '\n'+body.splitlines()[0]
    with pytest.raises(ValueError): helper.one_used_slot(body)


@pytest.mark.parametrize('field,value', [('Priority','16'), ('NumTasks','32'),
    ('NumCPUs','64'), ('CPUs/Task','1'), ('StdErr','wrong'), ('TimeLimit','12:00:00'),
    ('Dependency','afterok:28723655'), ('UserId','other(1)')])
def test_exact_pending_identity_required(field, value):
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    row = {**held_row(helper.TARGET, rules.PILOT), field:value}
    with pytest.raises(ValueError): helper.held_training(row_body(row),helper.TARGET,rules)


def test_second_wave_dependencies_cannot_disappear():
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    row = {**held_row('28709790',rules.PILOT), 'Dependency':'(null)'}
    with pytest.raises(ValueError): helper.held_training(row_body(row),'28709790',rules)


@pytest.fixture
def transaction(tmp_path, monkeypatch):
    root, wait, source, seed = [tmp_path/name for name in ('audit','wait','source','seed')]
    for path in (root,wait,source,seed): path.mkdir()
    monkeypatch.setattr(helper,'ROOT',root)
    monkeypatch.setattr(helper,'WAIT',wait)
    common_file, t_audit = tmp_path/'common.py', tmp_path/'T_audit.json'
    common_file.write_text('# synthetic nonexecuted dependency\n')
    t_audit.write_text('{}\n')
    monkeypatch.setattr(helper,'COMMON',common_file)
    monkeypatch.setattr(helper,'COMMON_SHA',helper.sha(common_file))
    monkeypatch.setattr(helper,'T_AUDIT',t_audit)
    monkeypatch.setattr(helper,'T_AUDIT_SHA',helper.sha(t_audit))
    terminal=tmp_path/'terminal'
    terminal.mkdir()
    (terminal/'vcneb_summary.json').write_text('{}\n')
    value={**report(),'workdir':str(terminal),
        'terminal_summary_sha256':helper.sha(terminal/'vcneb_summary.json')}
    (root/'audit_receipt.json').write_text(json.dumps(value))
    script,binary=source/'production.slurm',source/'binary'
    script.write_text('# synthetic only\n')
    binary.write_text('not executable\n')
    (seed/'manifest.json').write_text('{}\n')
    (seed/'seed.traj').write_bytes(b'synthetic fixture, not physical data')
    proof=dict(condition='strain_p0100',source_root=str(source),pilot_sha256=helper.sha(script),
        ABACUS_binary_sha256=helper.sha(binary),source_code_sha256={},reports=[dict(
        channel=helper.CHANNELS[helper.TARGET],condition='strain_p0100',seed_root=str(seed),
        manifest_sha256=helper.sha(seed/'manifest.json'),seed_files_sha256={'seed.traj':helper.sha(seed/'seed.traj')})])
    (wait/'preflight.json').write_text(json.dumps(proof))
    monkeypatch.setattr(helper,'PREFLIGHT_SHA',helper.sha(wait/'preflight.json'))
    calls=[]
    state=dict(ambiguous=False,bad_waiter=False,queue_race=False,queue_reads=0)
    def run(argv):
        calls.append(argv)
        if argv[0]=='sacct':
            job=argv[argv.index('-j')+1]
            status='RUNNING' if job==helper.LIVE else 'COMPLETED'
            return job+'|'+status+'|0:0|32|hfacnormal01|iai806'
        if argv[0]=='squeue':
            state['queue_reads']+=1
            body=queue()
            if state['queue_race'] and state['queue_reads']>1:
                body=body.replace('PO_to_M|PENDING','PO_to_M|RUNNING')
            return body
        if argv[:2]==['scontrol','release']:
            if state['ambiguous']: raise TimeoutError('possibly accepted native write')
            return ''
        job=argv[-1]
        row=held_row(job,str(script))
        if state['bad_waiter'] and job=='28709791': row['Priority']='16'
        return row_body(row)
    rules=SimpleNamespace(run=run,PILOT=str(script),BINARY=binary,
        fields=parser.fields,dependencies=parser.dependencies)
    common=SimpleNamespace(load_rules=lambda:rules,
        pending_seed_preflight=lambda path:dict(mock_only=True,new_DFT_calls=0))
    monkeypatch.setattr(helper,'load_common',lambda:common)
    return root,calls,state,helper.sha(root/'audit_receipt.json')


def test_mock_admission_only_releases_registered_target_once(transaction):
    root,calls,state,digest=transaction
    receipt=helper.release_once(digest)
    assert receipt['release_writes']==1 and receipt['dependency_writes']==0
    assert [c for c in calls if c[:2]==['scontrol','release']]==[['scontrol','release',helper.TARGET]]
    assert not any(c[0]=='sbatch' or c[:2]==['scontrol','update'] for c in calls)
    with pytest.raises(FileExistsError): helper.release_once(digest)


@pytest.mark.parametrize('failure',['bad_waiter','queue_race'])
def test_protected_waiters_and_second_capacity_read_block_write(transaction,failure):
    root,calls,state,digest=transaction
    state[failure]=True
    with pytest.raises(ValueError): helper.release_once(digest)
    assert not any(c[:2]==['scontrol','release'] for c in calls)


def test_ambiguous_release_recorded_not_replayed(transaction):
    root,calls,state,digest=transaction
    state['ambiguous']=True
    with pytest.raises(TimeoutError): helper.release_once(digest)
    assert 'reconcile_do_not_repeat' in (root/'admission_journal.jsonl').read_text()
    with pytest.raises(FileExistsError): helper.release_once(digest)
    assert sum(c[:2]==['scontrol','release'] for c in calls)==1


def test_registered_real_training_preflight_hash_and_condition():
    p=REPO/'benchmarks/hfo2_channels/20261008/clamped_G2_training_queue_E060_20261010/preflight.json'
    assert helper.sha(p)==helper.PREFLIGHT_SHA
    proof=json.loads(p.read_text(encoding='utf-8'))
    assert proof['condition']=='strain_p0100' and len(proof['source_code_sha256'])==338
    assert [r['channel'] for r in proof['reports']]==list(helper.CHANNELS.values())
