"""E091 terminal-audit/admission wiring; synthetic scheduler tests, no DFT."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/clamped_rolling_native_E091_20261011'


def module(name):
    spec = importlib.util.spec_from_file_location('E091_'+name, CASE/(name+'.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.mark.parametrize('parent,cap,initial,count,channel', [
    ('28723655',40,0,334,'PO_flip_T_pattern_preserving'),
    ('28709789',10,7,338,'PO_to_M'),
])
def test_original_two_runtime_and_cap_registry(parent,cap,initial,count,channel):
    audit = module('audit_native')
    assert set(audit.REGISTRY) == {'28723655','28709789'}
    r = audit.REGISTRY[parent]
    assert (r['cap'],r['initial'],r['count'],r['channel']) == (cap,initial,count,channel)
    assert r['work'].name == 'band'
    assert 'E054' in str(r['source']) if parent == '28723655' else 'E055' in str(r['source'])
    assert len(r['script_sha']) == len(r['archive_sha']) == 64
    assert len(audit.SUPPORT_SHA) == 3
    assert r['audit_subdir'] == ('terminal-r2' if parent=='28723655' else 'terminal')
    assert audit.sha(CASE.parents[3]/'scripts/prepare_hfo2_clamped_resume.py')==audit.PREPARE_SHA


def test_unknown_parent_rejected_before_native_file_access():
    with pytest.raises(ValueError,match='registered G2 parents'):
        module('audit_native').audit_parent('28709790')


def dispatcher(monkeypatch):
    monkeypatch.setitem(sys.modules,'audit_native',module('audit_native'))
    return module('dispatch_once')


def test_two_live_slots_or_unregistered_study_never_admitted(monkeypatch):
    d = dispatcher(monkeypatch)
    allowed = {'28723655','28709789','99900001'}
    assert d.spare_slot('', '28723655',allowed) == []
    assert d.spare_slot('28709789|hfo2-M|RUNNING|hfacnormal01|32','28723655',allowed) == ['28709789']
    for body in ('28723655|hfo2-flip|RUNNING|hfacnormal01|32',
                 '28709789|hfo2-M|RUNNING|hfacnormal01|32\n99900001|hfo2-next|CONFIGURING|hfacnormal01|32',
                 '99999999|hfo2-unknown|PENDING|hfacnormal01|32',
                 '28709789|hfo2-M|RUNNING|wrong|32',
                 '28709789|hfo2-M|RUNNING|hfacnormal01|16',
                 '28709789|hfo2-M|FAILED|hfacnormal01|32','malformed'):
        with pytest.raises(ValueError):
            d.spare_slot(body,'28723655',allowed)


def scheduler_fixture(tmp_path,monkeypatch,parent,mode='success'):
    d = dispatcher(monkeypatch)
    monkeypatch.setattr(d,'ROOT',tmp_path)
    monkeypatch.setattr(d,'SCRIPT',tmp_path/'wrapper.slurm')
    d.SCRIPT.write_text('# synthetic scheduler fixture; never run DFT\n')
    root = tmp_path/parent
    root.mkdir()
    calls = []
    monkeypatch.setattr(d,'verify_audit',lambda p,h:{'seed_manifest_sha256':'a'*64})
    def response(argv):
        calls.append(argv)
        if argv[0]=='sacct':
            return parent+'|'+('RUNNING' if mode=='parent_live' else 'COMPLETED')+'|0:0|32|hfacnormal01|iai806'
        if argv[0]=='squeue':
            if mode=='second_queue' and sum(a[0]=='squeue' for a in calls)==2:
                return '99999999|hfo2-unregistered|RUNNING|hfacnormal01|32'
            return ''
        if argv[0]=='bash':
            return ''
        if argv[0]=='sbatch':
            records=[json.loads(l) for l in (root/'dispatch_journal.jsonl').read_text().splitlines()]
            assert [r['event'] for r in records] == ['one_held_submission_intent']
            return '99900001;cluster' if mode=='ambiguous' else '99900001'
        if argv[:2]==['scontrol','release']:
            records=[json.loads(l) for l in (root/'dispatch_journal.jsonl').read_text().splitlines()]
            assert records[-1]['event']=='one_release_intent'
            if mode=='release_timeout':
                raise subprocess.TimeoutExpired(argv,45,output=b'uncertain reply',stderr=b'transport interrupted')
            return ''
        job = argv[-1]
        if job in d.WAITERS:
            channel = d.WAITERS[job]
            row = dict(JobId=job,JobName='hfo2-G2-E060-'+channel,UserId='iai806(16284)',JobState='PENDING',
                Partition='hfacnormal01',NumCPUs='32',Priority='0',Reason='JobHeldUser',Dependency='(null)',
                Command=str(d.REGISTRY['28709789']['script']),StdOut=str(d.WAIT_ROOT/(channel+'.slurm.out')),
                StdErr=str(d.WAIT_ROOT/(channel+'.slurm.err')),LastSchedEval=str(len(calls)))
            if mode=='held_original':row['Priority']='1'
        else:
            channel=d.REGISTRY[parent]['channel']
            row = dict(JobId=job,JobName='hfo2-G2-E091-'+parent,UserId='iai806(16284)',JobState='PENDING',
                Partition='hfacnormal01',NumCPUs='16' if mode=='accepted_identity' else '32',Priority='0',
                Reason='JobHeldUser',Dependency='(null)',Command=str(d.SCRIPT),
                StdOut=str(root/(channel+'.slurm.out')),StdErr=str(root/(channel+'.slurm.err')))
        return ' '.join(k+'='+v for k,v in row.items())
    monkeypatch.setattr(d,'run',response)
    return d,root,calls


@pytest.mark.parametrize('parent',['28723655','28709789'])
def test_exact_one_submission_release_and_durable_intents(tmp_path,monkeypatch,parent):
    d,root,calls=scheduler_fixture(tmp_path,monkeypatch,parent)
    r=d.dispatch_once(parent,'b'*64)
    assert r['new_independent_chains']==0 and not r['execution_started_claimed']
    assert r['total_submissions']==r['release_writes']==1
    assert r['dependency_writes']==0 and not r['FIRE_optimizer_state_restored']
    submit=next(a for a in calls if a[0]=='sbatch')
    export=next(a for a in submit if a.startswith('--export='))
    for text in ('PARENT_JOB_ID='+parent,'SOURCE_ROOT='+str(d.REGISTRY[parent]['source']),
                 'PARENT_AUDIT_SHA256='+'b'*64,'SEED_MANIFEST_SHA256='+'a'*64):
        assert text in export
    count=len(calls)
    with pytest.raises(FileExistsError):d.dispatch_once(parent,'b'*64)
    assert len(calls)==count
    assert d.admitted_other(next(iter(set(d.REGISTRY)-{parent}))) >= {r['job_id']}


@pytest.mark.parametrize('mode,expected_submit,expected_release',[
    ('parent_live',0,0),('held_original',0,0),('ambiguous',1,0),
    ('accepted_identity',1,0),('second_queue',1,0),('release_timeout',1,1),
])
def test_failure_intents_are_retained_never_replayed(tmp_path,monkeypatch,mode,expected_submit,expected_release):
    parent='28723655'
    d,root,calls=scheduler_fixture(tmp_path,monkeypatch,parent,mode)
    with pytest.raises((ValueError,subprocess.TimeoutExpired)):
        d.dispatch_once(parent,'b'*64)
    assert sum(a[0]=='sbatch' for a in calls)==expected_submit
    assert sum(a[:2]==['scontrol','release'] for a in calls)==expected_release
    journal=root/'dispatch_journal.jsonl'
    if expected_submit:
        records=[json.loads(l) for l in journal.read_text().splitlines()]
        assert records[-1]['event']=='reconcile_existing_intent_never_repeat_writes'
        if mode=='release_timeout':assert records[-1]['stdout']=='uncertain reply'
        count=len(calls)
        with pytest.raises(FileExistsError):d.dispatch_once(parent,'b'*64)
        assert len(calls)==count
    else:
        assert not journal.exists()
