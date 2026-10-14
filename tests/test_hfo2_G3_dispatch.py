"""Admission faults preserve handles; all subprocesses here are explicit mocks."""
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


@pytest.fixture
def helper(tmp_path,monkeypatch):
    path=Path(__file__).parents[1]/'benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/dispatch_HF.py'
    spec=importlib.util.spec_from_file_location('E150_dispatch',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module,'verify',lambda _: ({'manifest_sha256':'a'*64,'runtime_manifest_sha256':'b'*64},
        tmp_path/'source',tmp_path/'manifest.json',tmp_path/'wrapper.slurm'))
    return module


@pytest.mark.parametrize('fault',['none','submission_timeout','release_timeout','wrong_resources','wrong_node_range','still_held'])
def test_two_slots_no_replay_and_native_resources(helper,monkeypatch,fault):
    calls=[];handles={};released=set()
    def execute(argv):
        calls.append(argv)
        if argv[0]=='squeue':
            return '\n'.join(f"{job}|{name}|PENDING|32|({'Priority' if job in released and fault!='still_held' else 'JobHeldUser'})" for job,name in handles.items())
        if argv[0]=='bash':return ''
        if argv[0]=='sbatch':
            if fault=='submission_timeout': raise subprocess.TimeoutExpired(argv,45)
            job=str(101+len(handles));name=next(v.split('=',1)[1] for v in argv if v.startswith('--job-name='))
            handles[job]=name;return job
        if argv[:3]==['scontrol','show','job']:
            job=argv[-1];cpu=16 if fault=='wrong_resources' else 32
            nodes='1-2' if fault=='wrong_node_range' else '1-1'
            return f"JobId={job} JobName={handles[job]} Partition=hfacnormal01 NumNodes={nodes} NumCPUs={cpu} NumTasks=1 CPUs/Task=32 Command={helper.ROOT/'wrapper.slurm'} "
        if argv[:2]==['scontrol','release']:
            if fault=='release_timeout':raise subprocess.TimeoutExpired(argv,45)
            released.add(argv[-1])
            return ''
        raise AssertionError(argv)
    monkeypatch.setattr(helper,'run',execute)
    if fault=='none':
        result=helper.dispatch('c'*64)
        assert result['handles']=={'T':'101','preserving':'102'}
        assert result['total_planned_static_attempts']==296
        assert result['physical_holdout_generated_or_read'] is False
    else:
        with pytest.raises((ValueError,subprocess.TimeoutExpired)):helper.dispatch('c'*64)
        assert not (helper.ROOT/'dispatch_receipt.json').exists()
    before=len(calls)
    with pytest.raises(FileExistsError):helper.dispatch('c'*64)
    assert len(calls)==before
    assert len([a for a in calls if a[0]=='sbatch'])<=2
    events=[json.loads(x)['event'] for x in (helper.ROOT/'dispatch_journal.jsonl').read_text().splitlines()]
    if fault=='submission_timeout':assert events==['one_submission_intent']
    if fault=='wrong_resources':assert 'one_release_intent' not in events


def test_unknown_existing_job_refuses_any_submit(helper,monkeypatch):
    calls=[]
    monkeypatch.setattr(helper,'run',lambda a: calls.append(a) or '999|unrelated|RUNNING|32|node1')
    with pytest.raises(ValueError,match='unknown/missing'):helper.dispatch('c'*64)
    assert all(a[0]!='sbatch' for a in calls)
    assert not (helper.ROOT/'dispatch_journal.jsonl').exists()
