"""Analyze both original E269 terminals from one compressed native archive.

Reparse every real SCF, retain both initialization outcomes and never replace
the independently frozen E268 labels. No calculator or scheduler is invoked.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results,sha256

CASE=Path(__file__).resolve().parent


def restore(destination):
    archive=CASE/'completed_paths.tar.gz'
    receipt=json.loads((CASE/'collection_receipt.json').read_bytes())
    if (sha256(archive)!=receipt['archive_sha256'] or receipt['new_DFT_calls']!=0
            or receipt['scheduler_mutations']!=0):
        raise ValueError('the complete original E269 terminal archive is required')
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            p=Path(member.name)
            if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='completed_paths'
                    or not(member.isdir() or member.isfile())
                    or not(destination/p).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe completed native archive member')
        tar.extractall(destination)
    folder=destination/'completed_paths'
    files=json.loads((folder/'files.json').read_bytes())
    if len(files)!=receipt['exported_files']:raise ValueError('complete payload count differs')
    for name,digest in files.items():
        p=folder/name
        if not p.resolve().is_relative_to(folder.resolve()) or sha256(p)!=digest:
            raise ValueError('archived native bytes changed')
    return folder,receipt


def replay_calls(folder):
    audit=json.loads((folder/'audit.json').read_bytes())
    expected=7*(audit['step']+1)
    if (not audit['ordinary_converged'] or len(audit['calls'])!=expected
            or audit['new_SCF_calls']!=expected
            or audit['scheduler_row'][1:7]!=['COMPLETED','0:0','32','hfacnormal01','iai806','1']):
        raise ValueError('exact completed ordinary result and all actual calls required')
    seconds=[];image_counts={f'image_{i:04d}':0 for i in range(1,8)}
    for i,record in enumerate(audit['calls']):
        raw=folder/'all_calls'/f'{i:04d}'
        source_image=Path(record['directory']).parent.name
        if source_image not in image_counts:raise ValueError('fixed endpoint recalculation is forbidden')
        image_counts[source_image]+=1
        call=json.loads((raw/'call_audit.json').read_bytes())
        if (sha256(raw/'call_audit.json')!=record['call_audit_sha256']
                or call['input_sha256']!=record['input_sha256']
                or any(call['input_sha256'][n]!=h for n,h in CONTRACT.items())
                or any(sha256(raw/n)!=call['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log')!=record['raw_log_sha256']
                or call['raw_log_sha256']!=record['raw_log_sha256']):
            raise ValueError('original native input/output differs')
        parsed=audited_results(raw)
        if any(not np.allclose(parsed[k],call['results'][k],atol=1e-12,rtol=0) for k in parsed):
            raise ValueError('native intermediate energy/force/stress differs')
        elapsed=call['elapsed_seconds']
        if not np.isfinite(elapsed) or elapsed<=0 or elapsed!=record['elapsed_seconds']:
            raise ValueError('positive unchanged native timing required')
        seconds.append(elapsed)
    if (any(n!=audit['step']+1 for n in image_counts.values())
            or len(list((folder/'all_calls').glob('*/call_audit.json')))!=expected
            or not np.isclose(sum(seconds)*32/3600,audit['SCF_core_hours'],atol=1e-10,rtol=0)):
        raise ValueError('complete per-image call count or cost differs')
    return dict(job_id=audit['job_id'],step=audit['step'],calls=expected,
        all_native_calls_reparsed_locally=True,native_proxy_core_hours=audit['SCF_core_hours'],
        allocated_core_hours=audit['allocation_core_hours'],full_pseudo_orbital_bytes_rechecked_locally=False)


def analyze(output,*,draw=False):
    spec=importlib.util.spec_from_file_location('E269_completed_comparison',CASE/'compare_paths.py')
    comparison=importlib.util.module_from_spec(spec);spec.loader.exec_module(comparison)
    with tempfile.TemporaryDirectory(prefix='varneb_E269_completed_') as tmp:
        exported,receipt=restore(Path(tmp))
        replays={c:replay_calls(exported/c) for c in ('preserving','T')}
        data=comparison.compare(ROOT,exported)
        data['native_cohort']=dict(archive_sha256=receipt['archive_sha256'],channels=replays,
            all_native_calls_reparsed_locally=True,new_DFT_calls=0)
        comparison.write_result(data,output,draw=draw)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--draw',action='store_true')
    a=p.parse_args();result=analyze(a.output,draw=a.draw)
    print(json.dumps(dict(status=result['status'],outcomes=result['outcomes'],paired=result['paired']),allow_nan=False))
