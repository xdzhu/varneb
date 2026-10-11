"""Pinned child audit wiring; synthetic accounting is not native convergence."""
import importlib.util
import json
from pathlib import Path
import shlex
import subprocess
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_child_audit_E093_20261011'


def module():
    spec = importlib.util.spec_from_file_location('E093_audit_child', CASE/'audit_child_HF.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def native_row(m, job, state='RUNNING', *, change=None):
    r = m.REGISTRY[job]
    p = m.paths(r)
    env = dict(RUN_DFT='1', WORKDIR=str(p['work']), SOURCE_ROOT=str(r['source']),
        PARENT_JOB_ID=r['parent'], PARENT_AUDIT_PATH=str(p['parent_audit']),
        PARENT_AUDIT_SHA256=r['parent_audit_sha'], SEED_ROOT=str(p['seed']),
        SEED_MANIFEST_SHA256=r['seed_sha'])
    if change in env:
        env[change] = 'different'
    if change == 'physical_override':
        env['ECUTWFC'] = '120'
    submit = 'sbatch '+shlex.quote('--export=ALL,'+','.join(k+'='+v for k,v in env.items()))+' '+shlex.quote(str(m.SCRIPT))
    if change == 'duplicate_export': submit += ' --export=ALL'
    fields = [job, state, 'iai806', '32', 'hfacnormal01', 'hfo2-G2-E091-'+r['parent'], str(p['stdout']), str(p['stderr']), submit]
    if isinstance(change, int): fields[change] = 'different'
    return '|'.join(fields)


@pytest.mark.parametrize('job,count', [('28727368',338), ('28727382',334)])
def test_exact_children_and_complete_matching_support(job, count):
    m = module()
    assert set(m.REGISTRY) == {'28727368','28727382'}
    assert m.REGISTRY[job]['count'] == count
    assert m.paths(m.REGISTRY[job])['work'].name == 'band'
    assert len(m.SUPPORT_SHA) == 4
    for name,digest in m.SUPPORT_SHA.items():
        assert m.sha(REPO/'scripts'/name) == digest
    assert m.sha(REPO/'cluster/hf_hfo2_clamped_continue_E090_20261011.slurm') == m.SCRIPT_SHA
    assert 'a8b23b9355' in m.SUPPORT_SHA['export_hfo2_clamped_observation.py']


@pytest.mark.parametrize('change', [0,1,2,3,4,5,6,7,'RUN_DFT','WORKDIR','SOURCE_ROOT',
    'PARENT_JOB_ID','PARENT_AUDIT_PATH','PARENT_AUDIT_SHA256','SEED_ROOT','SEED_MANIFEST_SHA256',
    'physical_override','duplicate_export'])
def test_wrong_native_identity_or_environment_rejected(change):
    m = module()
    with pytest.raises(ValueError):
        m.validate_native_identity(native_row(m,'28727368',change=change),'28727368',m.REGISTRY['28727368'],preflight_only=True)


@pytest.mark.parametrize('state', ['PENDING','CONFIGURING','FAILED','TIMEOUT','CANCELLED','RUNNING'])
def test_live_failed_or_waiting_job_never_gets_terminal_audit(state):
    m = module()
    with pytest.raises(ValueError, match='terminal state'):
        m.validate_native_identity(native_row(m,'28727368',state),'28727368',m.REGISTRY['28727368'])


def test_only_preflight_accepts_untransformed_running_row():
    m = module()
    body = native_row(m,'28727368')
    assert m.validate_native_identity(body,'28727368',m.REGISTRY['28727368'],preflight_only=True) == 'RUNNING'
    assert '|RUNNING|' in body
    with pytest.raises(ValueError):
        m.validate_native_identity(body+'\n'+body,'28727368',m.REGISTRY['28727368'],preflight_only=True)


@pytest.mark.parametrize('job', ['28709789','28723655','28709790','unknown'])
def test_old_parents_waiters_and_unknown_reject_before_io(job, monkeypatch):
    m = module()
    def forbidden(*args,**kwargs): raise AssertionError('native IO must not happen')
    monkeypatch.setattr(m.subprocess,'run',forbidden)
    with pytest.raises(ValueError,match='accepted E091 children'):
        m.audit_child(job)


def fixture(tmp_path, monkeypatch, state='RUNNING'):
    m = module()
    monkeypatch.setattr(m,'AUDIT_ROOT',tmp_path/'audit-only')
    monkeypatch.setattr(m,'PRODUCTION',tmp_path/'production')
    r = m.REGISTRY['28727368']
    work = m.paths(r)['work']
    work.mkdir(parents=True)
    (work/'vcneb_preflight.json').write_text(json.dumps({'fixture_not_DFT':True}))
    rows = []
    def query(args,**kwargs):
        assert args[:6] == ['sacct','-X','-n','-P','-j','28727368']
        rows.append(args)
        return subprocess.CompletedProcess(args,0,stdout=native_row(m,'28727368',state))
    monkeypatch.setattr(m.subprocess,'run',query)
    return m,r,work,rows


def test_running_native_state_rejects_before_summary_or_import(tmp_path, monkeypatch):
    m,r,work,rows = fixture(tmp_path,monkeypatch)
    def forbidden(*args): raise AssertionError('do not import/read terminal data for a live job')
    monkeypatch.setattr(m,'load_auditor',forbidden)
    with pytest.raises(ValueError,match='terminal state'):
        m.audit_child('28727368')
    assert len(rows) == 1
    assert not m.AUDIT_ROOT.exists()
    assert not (work/'vcneb_summary.json').exists()


def test_preflight_counts_zero_initial_calls_but_creates_no_terminal_proof(tmp_path, monkeypatch):
    m,r,work,rows = fixture(tmp_path,monkeypatch)
    checked=[]
    def verify(*args):
        assert args == (r['source'],r['archive'],r['archive_sha'])
        checked.append('runtime')
        return {str(i):'pinned' for i in range(r['count'])}
    def count(step,preflight,digest):
        assert (step,preflight,digest) == (40,{'fixture_not_DFT':True},m.SCRIPT_SHA)
        checked.append('cache policy')
        return 0,280
    monkeypatch.setattr(m,'load_auditor',lambda registration:SimpleNamespace(verify_runtime=verify,fresh_interior_call_count=count))
    result=m.audit_child('28727368',preflight_only=True)
    assert checked == ['runtime','cache policy']
    assert result['native_state'] == 'RUNNING'
    assert result['all_native_SCFs_or_final_convergence_audited'] is False
    assert result['new_DFT_calls'] == result['scheduler_writes'] == result['output_files_created'] == 0
    assert len(rows) == 1 and not m.AUDIT_ROOT.exists()


def test_existing_audit_is_not_overwritten_or_queried(tmp_path, monkeypatch):
    m,r,work,rows=fixture(tmp_path,monkeypatch,'COMPLETED')
    (m.AUDIT_ROOT/'28727368/terminal').mkdir(parents=True)
    with pytest.raises(FileExistsError):
        m.audit_child('28727368')
    assert rows == []


@pytest.mark.parametrize('step', [39,41,-1])
def test_unexpected_cap_rejects_before_material_auditor(tmp_path,monkeypatch,step):
    m,r,work,rows=fixture(tmp_path,monkeypatch,'COMPLETED')
    (work/'vcneb_summary.json').write_text(json.dumps({'status':'max_steps_reached'}))
    def forbidden(*args,**kwargs): raise AssertionError('unexpected cap must not be audited/admitted')
    monkeypatch.setattr(m,'load_auditor',lambda registration:SimpleNamespace(completed_step=lambda path:step,audit=forbidden))
    with pytest.raises(ValueError,match='finite segment'):
        m.audit_child('28727368')
    assert not m.AUDIT_ROOT.exists()


@pytest.mark.parametrize('status,step,capped', [('converged',8,False),('max_steps_reached',40,True)])
def test_terminal_delegation_keeps_actual_child_and_old_runtime(tmp_path,monkeypatch,status,step,capped):
    m,r,work,rows=fixture(tmp_path,monkeypatch,'COMPLETED')
    (work/'vcneb_summary.json').write_text(json.dumps({'status':status}))
    calls=[]
    def actual_audit(*args,**kwargs):
        assert args[:6] == (work,'28727368',m.SCRIPT,r['source'],r['archive'],r['archive_sha'])
        assert args[6] == m.AUDIT_ROOT/'28727368/terminal'
        assert kwargs == dict(outcome='step-cap' if capped else 'converged',seed_latest=capped)
        calls.append('terminal audit, synthetic not DFT')
        output=args[6]
        output.mkdir(parents=True)
        receipt=dict(runtime_code_sha256={str(i):'pinned' for i in range(r['count'])},
            ordinary_converged=not capped,terminal_step=step)
        (output/'audit_receipt.json').write_text(json.dumps(receipt))
        return receipt
    monkeypatch.setattr(m,'load_auditor',lambda registration:SimpleNamespace(completed_step=lambda path:step,audit=actual_audit))
    receipt=m.audit_child('28727368')
    assert len(rows) == len(calls) == 1
    assert receipt['ordinary_converged'] is not capped
    record=json.loads((m.AUDIT_ROOT/'28727368/terminal/driver_receipt.json').read_text())
    assert (record['job_id'],record['parent_job']) == ('28727368','28709789')
    assert record['new_DFT_calls'] == record['scheduler_writes'] == record['new_independent_chains'] == 0
