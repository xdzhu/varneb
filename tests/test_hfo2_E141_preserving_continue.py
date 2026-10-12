"""Review/admission fixtures plus the real audited child40; zero DFT."""
import importlib.util
import json
import pathlib
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_preserving_continue_E141_20261012'
SCRIPT = REPO/'cluster/hf_hfo2_clamped_continue_E141_20261012.slurm'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manager():
    return load(CASE/'continue_preserving_HF.py','E141_fixture')


def recipe():
    m = load(REPO/'tests/test_hfo2_rolling_continue_recipe.py','E141_recipe_fixture')
    m.SCRIPT = SCRIPT
    m.PARENTS = {'28737465':(40,0,338,m.RUNS+'20261010-varneb-nested-controls-E055-r1/source',
        m.RUNS+'20261012-varneb-clamped-flip-continue-E103-r1/28709790/PO_flip_T_pattern_preserving/band')}
    return m


def test_exact_new_embedded_startup_with_synthetic_child40(tmp_path,monkeypatch,capsys):
    m = recipe()
    m.run_validator(*m.fixture(tmp_path,'28737465'),'28737465',monkeypatch)
    assert float(capsys.readouterr().out)==5.12968067458423


@pytest.mark.parametrize('field,value',[
    ('terminal_step',10),('ordinary_converged',True),('fresh_interior_SCFs',77),
    ('initial_fresh_interior_SCFs',7),('current_frame_exact_caches_checked',8),
    ('holdout_generated_or_read',True),('physical_inputs_or_running_source_changed',True),
    ('FIRE_state_restored',True),('new_DFT_calls',1),('scheduler_mutations',1),
    ('scheduler_row','28737465|RUNNING|0:0|32|hfacnormal01|iai806')])
def test_startup_does_not_masquerade_child40_as_starter(tmp_path,monkeypatch,field,value):
    m = recipe(); seed,audit,runtime,source = m.fixture(tmp_path,'28737465')
    r = json.loads(audit.read_text()); r[field]=value; m.write_json(audit,r)
    with pytest.raises(ValueError):
        m.run_validator(seed,audit,runtime,source,'28737465',monkeypatch)


def test_original_runtime_command_and_actual_parent_evidence():
    m = manager()
    old = REPO/'cluster/hf_hfo2_clamped_continue_E103_20261012.slurm'
    assert SCRIPT.read_bytes().split(b'module purge',1)[1]==old.read_bytes().split(b'module purge',1)[1]
    assert m.sha(SCRIPT)==m.SCRIPT_SHA
    cap = REPO/'benchmarks/hfo2_channels/20261008/clamped_preserving_cap_E140_20261012/completed_HF'
    assert m.sha(cap/'audit_receipt.json')==m.AUDIT_SHA
    assert m.sha(cap/'seed/manifest.json')==m.SEED_SHA
    assert m.sha(REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012/continue_flips_HF.py')==m.LEGACY_SHA
    record = json.loads((cap/'audit_receipt.json').read_text())
    assert record['terminal_step']==40 and record['job_id']==m.PARENT
    assert record['fresh_interior_SCFs']==280 and not record['ordinary_converged']
    assert m.environment()['PARENT_JOB_ID']=='28737465'


@pytest.mark.parametrize('change',[None,'live','audit','seed','binary','summary','namespace','support'])
def test_actual_manager_parent_verifier_executes_exact_startup(tmp_path,monkeypatch,change):
    from types import SimpleNamespace
    m=manager(); f=recipe()
    seed,audit,runtime,source=f.fixture(tmp_path,'28737465')
    m.ROOT=tmp_path/'new-study'; m.TERMINAL=audit.parent; m.SOURCE=source; m.SCRIPT=SCRIPT
    # Production is POSIX. Route only its immutable source lookup, as in the
    # already-tested startup fixture; don't change the actual embedded code.
    target_seed=m.TERMINAL/'seed'; shutil.copytree(seed.parent,target_seed)
    m.SEED_SHA=f.digest(target_seed/'manifest.json')
    native=tmp_path/'old-native'; native.mkdir()
    (native/'vcneb_summary.json').write_text('{"synthetic_not_DFT":true}')
    record=json.loads(audit.read_text())
    record['terminal_summary_sha256']=f.digest(native/'vcneb_summary.json')
    support={name:'a'*64 for name in ('audit_hfo2_clamped_terminal.py',
        'export_hfo2_clamped_observation.py','analyze_hfo2_clamped_residual.py')}
    for key,name in [('executed_auditor_sha256','audit_hfo2_clamped_terminal.py'),
                     ('executed_exporter_sha256','export_hfo2_clamped_observation.py'),
                     ('executed_residual_analyzer_sha256','analyze_hfo2_clamped_residual.py')]:
        record[key]=support[name]
    f.write_json(audit,record); m.AUDIT_SHA=f.digest(audit)
    binary=tmp_path/'binary'; binary.write_text('fixture_not_executable')
    commands=[]; identities=[]
    def identity(body,job,parent,*,child):
        identities.append((body,job,parent,child))
        if 'RUNNING' in body: raise ValueError('native live fixture')
    old=SimpleNamespace(accepted=lambda parent:{'job_id':'28737465'},identity=identity,
        BINARY=binary,BINARY_SHA=f.digest(binary),paths=lambda parent:{'work':native},
        load_auditor=lambda:(None,None,support))
    monkeypatch.setattr(m,'legacy',lambda:old)
    def run(argv):
        commands.append(argv)
        assert argv[0]=='sacct'
        return 'native '+('RUNNING' if change=='live' else 'COMPLETED')
    monkeypatch.setattr(m,'run',run)
    if change=='audit': audit.write_text('changed')
    elif change=='seed': (target_seed/'manifest.json').write_text('changed')
    elif change=='binary': binary.write_text('changed')
    elif change=='summary': (native/'vcneb_summary.json').write_text('changed')
    elif change=='namespace': m.paths()['work'].mkdir(parents=True)
    elif change=='support': support['audit_hfo2_clamped_terminal.py']='b'*64
    concrete=type(runtime)
    with monkeypatch.context() as context:
        context.setattr(pathlib,'Path',lambda p:runtime if str(p)==source else concrete(p))
        if change is None: assert m.verify_parent()['terminal_step']==40
        else:
            with pytest.raises(ValueError): m.verify_parent()
    assert identities==[('native '+('RUNNING' if change=='live' else 'COMPLETED'),'28737465','28709790',True)]
    assert len(commands)==1 and commands[0][0]=='sacct'


def old_slot_reader(tmp_path):
    old_case = REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012'
    old = load(old_case/'continue_flips_HF.py','E141_old_slots')
    old.ROOT = tmp_path/'old'
    for parent,source in [('28709790',old_case),('28709791',
            REPO/'benchmarks/hfo2_channels/20261008/clamped_reversing_continue_E113_20261012')]:
        out = old.ROOT/parent; out.mkdir(parents=True)
        for name in ('dispatch_receipt.json','dispatch_journal.jsonl'):
            shutil.copyfile(source/name,out/name)
    return old


def test_native_slot_ceiling_and_held_vs_runnable_pending(tmp_path,monkeypatch):
    m = manager(); m.ROOT=tmp_path/'study'
    old = old_slot_reader(tmp_path); monkeypatch.setattr(m,'legacy',lambda:old)
    running = '28764704|hfo2-G2-E103-28709791|RUNNING|hfacnormal01|32|node6'
    extra={'99900001':'hfo2-G2-E141-28737465'}
    held = running+'\n99900001|hfo2-G2-E141-28737465|PENDING|hfacnormal01|32|(JobHeldUser)'
    assert m.slots(held,spare=True,extra=extra)==['28764704']
    with pytest.raises(ValueError):
        m.slots(held.replace('(JobHeldUser)','(Priority)'),spare=True,extra=extra)
    assert m.slots(held.replace('(JobHeldUser)','(Priority)'),spare=False,extra=extra)==['28764704','99900001']
    with pytest.raises(ValueError):
        m.slots(running+'\n99900002|hfo2-unknown|RUNNING|hfacnormal01|32|node1',spare=False)
    with pytest.raises(ValueError):
        m.slots(running+'\n'+running,spare=False)


def fixture(tmp_path,monkeypatch,mode='success'):
    m = manager(); m.ROOT=tmp_path/'study'; m.ROOT.mkdir()
    m.LOCK=tmp_path/'shared.lock'; m.SCRIPT=tmp_path/'wrapper.slurm'
    m.SCRIPT.write_text('# fixture, not a production wrapper\n')
    old = old_slot_reader(tmp_path); monkeypatch.setattr(m,'legacy',lambda:old)
    monkeypatch.setattr(m,'verify_parent',lambda:dict(terminal_step=40))
    calls=[]; released=[]
    running='28764704|hfo2-G2-E103-28709791|RUNNING|hfacnormal01|32|node6'
    def run(argv):
        calls.append(argv)
        if argv[0]=='sacct':
            return '28737465|'+('RUNNING' if mode=='live' else 'COMPLETED')+'|0:0|32|hfacnormal01|iai806'
        if argv[0]=='bash': return ''
        if argv[0]=='squeue':
            if mode=='capacity' and sum(c[0]=='squeue' for c in calls)>=2:
                return running+'\n28737465|hfo2-G2-E103-28709790|RUNNING|hfacnormal01|32|node25'
            if any(c[0]=='sbatch' for c in calls):
                return running+'\n99900001|hfo2-G2-E141-28737465|PENDING|hfacnormal01|32|'+('(Priority)' if released else '(JobHeldUser)')
            return running
        if argv[0]=='sbatch':
            records=[json.loads(line) for line in (m.ROOT/'dispatch_journal.jsonl').read_text().splitlines()]
            assert records[-1]['event']==m.EVENTS[0]
            return '99900001;ambiguous' if mode=='ambiguous' else '99900001'
        if argv[:2]==['scontrol','release']:
            if mode=='timeout': raise subprocess.TimeoutExpired(argv,45)
            if mode!='stale': released.append(True)
            return ''
        assert argv[:3]==['scontrol','show','job']
        p=m.paths()
        row=dict(JobId='99900001',JobName='hfo2-G2-E141-28737465',JobState='PENDING',
            Reason='Priority' if released else 'JobHeldUser',Priority='1' if released else '0',
            Partition='hfacnormal01',NumCPUs='16' if mode=='wrong_cpu' else '32',NumTasks='1',
            Command=str(m.SCRIPT),Dependency='(null)',StdOut=str(p['stdout']),StdErr=str(p['stderr']),
            TimeLimit='12:00:00',Restarts='0',NumNodes='2' if mode=='wrong_nodes' else '1',
            UserId='iai806(1000)',**{'CPUs/Task':'32'})
        return ' '.join(k+'='+shlex.quote(v) for k,v in row.items())
    monkeypatch.setattr(m,'run',run)
    return m,calls


def test_exact_single_submission_release_and_no_replay(tmp_path,monkeypatch):
    m,calls=fixture(tmp_path,monkeypatch)
    result=m.dispatch()
    assert result['source_step']==40 and result['parent_job']=='28737465'
    assert result['total_submissions']==result['release_writes']==1
    assert result['dependency_writes']==result['new_independent_chains']==0
    assert not result['FIRE_optimizer_state_restored'] and not result['execution_started_claimed']
    assert m.accepted()==result and not m.LOCK.exists()
    before=len(calls)
    with pytest.raises(FileExistsError): m.dispatch()
    assert len(calls)==before
    records=[json.loads(line) for line in (m.ROOT/'dispatch_journal.jsonl').read_text().splitlines()]
    assert [r['event'] for r in records]==m.EVENTS
    assert records[-1]['committed_slots']==['28764704','99900001']


@pytest.mark.parametrize('mode,submits,releases',[('live',0,0),('capacity',0,0),
    ('ambiguous',1,0),('wrong_cpu',1,0),('wrong_nodes',1,0),('timeout',1,1),('stale',1,1)])
def test_faults_retain_intent_and_never_retry(tmp_path,monkeypatch,mode,submits,releases):
    m,calls=fixture(tmp_path,monkeypatch,mode)
    with pytest.raises((ValueError,subprocess.TimeoutExpired)): m.dispatch()
    assert sum(c[0]=='sbatch' for c in calls)==submits
    assert sum(c[:2]==['scontrol','release'] for c in calls)==releases
    assert not (m.ROOT/'dispatch_receipt.json').exists()
    if mode!='live':
        assert m.LOCK.exists() and (m.ROOT/'dispatch_journal.jsonl').exists()
        before=len(calls)
        with pytest.raises(FileExistsError): m.dispatch()
        assert len(calls)==before


def test_actual_live_child_refused_before_support_imports(tmp_path,monkeypatch):
    m=manager(); m.ROOT=tmp_path/'study'; m.ROOT.mkdir()
    monkeypatch.setattr(m,'accepted',lambda:dict(job_id='99900001'))
    monkeypatch.setattr(m,'run',lambda argv:'99900001|RUNNING|iai806|32|hfacnormal01|name|out|err|submit')
    monkeypatch.setattr(m,'legacy',lambda:pytest.fail('live child must reject before imports'))
    with pytest.raises(ValueError): m.audit_child()
    assert not (m.ROOT/'child_terminal').exists()


def test_actual_one_admission_and_durable_journal_bytes():
    m=manager()
    receipt=json.loads((CASE/'dispatch_receipt.json').read_text())
    assert m.sha(CASE/'dispatch_receipt.json')=='827b5530bc2543a6c172483b1ad30e9ac294f019cc79085c4464c78a7b84221f'
    assert receipt['job_id']=='28812952' and receipt['source_step']==40
    assert receipt['parent_job']=='28737465'
    assert receipt['total_submissions']==receipt['release_writes']==1
    assert receipt['dependency_writes']==receipt['new_independent_chains']==0
    assert receipt['executed_helper_sha256']==m.sha(CASE/'continue_preserving_HF.py')
    assert receipt['script_sha256']==m.sha(SCRIPT)
    journal=CASE/'dispatch_journal.jsonl'
    assert m.sha(journal)==receipt['journal_sha256']=='015f77d3e30fb8b574ca2aeaabccef3440b64d7cf93ac45bfd4751f4d7886ffa'
    records=[json.loads(line) for line in journal.read_text().splitlines()]
    assert [r['event'] for r in records]==m.EVENTS
    assert records[0]['committed_slots']==['28764704']
    assert records[-1]['committed_slots']==['28764704','28812952']
    assert records[2]['accepted_native_state']['Reason']=='JobHeldUser'
    assert records[-1]['native_control']['Reason']!='JobHeldUser'


def test_actual_startup_is_cache_and_MPI_evidence_not_terminal_success():
    m=manager()
    r=json.loads((CASE/'startup_observation.json').read_text())
    assert r['source_job_id']=='28812952'
    assert '28812952|RUNNING|0:0|32|hfacnormal01|iai806|2026-10-12T09:49:29' in r['native_accounting']
    assert r['committed_study_slots']==['28764704','28812952']
    assert r['actual_wrapper_sha256']==m.SCRIPT_SHA
    assert r['actual_helper_sha256']==m.sha(CASE/'continue_preserving_HF.py')
    assert r['seed_cache_count']==r['seed_cache_audit_files']==9 and r['none_cache_count']==0
    assert r['first_native_DSIZE']==['32']
    assert r['completed_fresh_call_audits']==1
    assert 'CheckedFIRE:    0 09:50:15    -9782.379267        0.114147' in r['optimizer_log_raw']
    assert r['stderr_bytes']==0 and not r['terminal_summary_exists'] and not r['failure_exists']
    assert r['no_new_DFT_initiated_by_reader'] and r['no_scheduler_writes_by_reader']


def test_one_actual_fresh_SCF_was_reparsed_not_the_whole_live_run():
    r=json.loads((CASE/'first_fresh_SCF_replay.json').read_text())
    assert r['source_job_id']=='28812952' and r['actual_native_DSIZE']==32
    assert r['six_physical_bytes_checked'] and r['raw_SCFs_reparsed']==1
    assert r['elapsed_seconds']>0
    assert len(r['raw_log_sha256'])==len(r['call_audit_sha256'])==64
    assert r['new_DFT_initiated_by_reader']==0
