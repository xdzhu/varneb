"""Recovery releases only the two accepted handles; it cannot submit."""
import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize('fault',['none','wrong_nodes','wrong_CPUs','not_held','release_timeout','wrong_handle'])
def test_exact_existing_handle_release_and_no_replay(tmp_path,monkeypatch,fault):
    path=Path(__file__).parents[1]/'benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/release_existing_HF.py'
    spec=importlib.util.spec_from_file_location('E150_release',path)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    monkeypatch.setattr(helper,'ROOT',tmp_path)
    records=[]
    for c,j in helper.HANDLES.items():
        records.extend([dict(event='one_submission_intent',channel=c),dict(event='one_handle_accepted',channel=c,job_id=j)])
    if fault=='wrong_handle': records[-1]['job_id']='999'
    (tmp_path/'dispatch_journal.jsonl').write_text('\n'.join(json.dumps(r) for r in records))
    calls=[];released=set();wrapper=tmp_path/'wrapper.slurm'
    def execute(argv):
        calls.append(argv)
        assert argv[0]!='sbatch'
        if argv[:3]==['scontrol','show','job']:
            j=argv[-1];c=next(c for c,jj in helper.HANDLES.items() if jj==j)
            nodes='1-2' if fault=='wrong_nodes' else '1-1'
            cpus='16' if fault=='wrong_CPUs' else '32'
            reason='Priority' if fault=='not_held' else 'JobHeldUser'
            values=dict(JobId=j,JobName='hfo2-G3-E150-'+c,Partition='hfacnormal01',JobState='PENDING',
                Reason=reason,NumCPUs=cpus,NumTasks='1',NumNodes=nodes,Command=str(wrapper),
                WorkDir=str(tmp_path),StdOut=str(tmp_path/(c+'.slurm.out')),StdErr=str(tmp_path/(c+'.slurm.err')))
            values['CPUs/Task']='32'
            return ' '.join(k+'='+v for k,v in values.items())
        if argv[:2]==['scontrol','release']:
            if fault=='release_timeout':raise subprocess.TimeoutExpired(argv,45)
            released.add(argv[-1]);return ''
        raise AssertionError(argv)
    def queue(accepted):
        return '\n'.join(f"{j}|{n}|PENDING|32|({'Priority' if j in released else 'JobHeldUser'})" for j,n in accepted.items())
    legacy=SimpleNamespace(verify=lambda _: ({},tmp_path/'source',tmp_path/'manifest.json',wrapper),run=execute,queue=queue)
    monkeypatch.setattr(helper,'original',lambda:legacy)
    if fault=='none':
        result=helper.release()
        assert result['additional_submissions']==0 and result['total_release_writes']==2
        assert released==set(helper.HANDLES.values())
        count=len(calls)
        with pytest.raises(FileExistsError):helper.release()
        assert len(calls)==count
    else:
        with pytest.raises((ValueError,subprocess.TimeoutExpired)):helper.release()
        assert not (tmp_path/'reconciled_dispatch_receipt.json').exists()
        if fault=='release_timeout':
            count=len(calls)
            with pytest.raises(FileExistsError):helper.release()
            assert len(calls)==count
        else:
            assert not released
