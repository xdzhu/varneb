#!/bin/bash
# One read-only scheduler observation and restricted portable evidence export.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
from datetime import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import tarfile

root = Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E060-r1')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
preflight = json.loads((root/'preflight.json').read_text())
receipt = json.loads((root/'submission_receipt.json').read_text())
assert receipt['preflight_sha256'] == sha(root/'preflight.json')
assert receipt['submission_journal_sha256'] == sha(root/'submission_journal.jsonl')
driver = root/'executed_queue_training_waves.py'
assert sha(driver) == preflight['queue_helper_sha256']
spec = importlib.util.spec_from_file_location('executed_E060',driver)
wave = importlib.util.module_from_spec(spec);spec.loader.exec_module(wave)
assert sha(wave.SOURCE/wave.PILOT) == wave.PILOT_SHA
assert [r['channel'] for r in receipt['jobs']] == list(wave.CHANNELS)
def run(argv):return subprocess.run(argv,check=True,capture_output=True,text=True,timeout=45).stdout
ids = [j['job_id'] for j in receipt['jobs']]
assert len(set(ids+list(wave.KNOWN))) == 8
records = []
for j in receipt['jobs']:
    raw = run(['scontrol','show','job','-o',j['job_id']])
    f = dict(t.split('=',1) for t in raw.split() if '=' in t)
    assert f['JobId'] == j['job_id'] and f['JobName'] == 'hfo2-G2-E060-'+j['channel']
    assert f['UserId'].startswith('iai806(') and f['Partition'] == 'hfacnormal01'
    assert f['NumCPUs'] == f['CPUs/Task'] == '32' and f['NumTasks'] == '1'
    assert f['NumNodes'] in ('1','1-1') and f['TimeLimit'] == '04:00:00'
    assert f['Command'] == (wave.SOURCE/wave.PILOT).as_posix()
    expected = set(j['dependency'].split(':')[1:])
    actual = [re.fullmatch(r'afterok:(\d+)\(unfulfilled\)',t) for t in f['Dependency'].split(',')]
    assert len(actual) == 2 and all(m is not None for m in actual)
    assert {m.group(1) for m in actual} == expected
    assert f['JobState'] == 'PENDING' and f['Reason'] == 'Dependency'
    records.append(dict(channel=j['channel'],job_id=j['job_id'],actual_fields=f,
                        control_raw=raw,resource_and_dependency_checks_passed=True))
all_ids = list(wave.KNOWN)+ids
acct = run(['sacct','-X','-n','-j',','.join(all_ids),'-P',
            '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
rows = [s.split('|') for s in acct.splitlines() if s.strip()]
assert len(rows) == 8 and {r[0] for r in rows} == set(all_ids)
assert all(r[2] == '0:0' and r[4:] == ['hfacnormal01','iai806'] for r in rows)
queue = run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C'])
study = [[s.strip() for s in t.split('|')] for t in queue.splitlines() if '|hfo2' in t]
assert len(study) == 8 and {r[0] for r in study} == set(all_ids)
assert all(r[3:] == ['hfacnormal01','32'] for r in study)
assert sum(r[2] != 'PENDING' for r in study) == 2
proof = dict(check_CST=datetime.now().isoformat(),jobs=records,accounting_raw=acct,
    accounting_rows=rows,queue_raw=queue,active_study_allocations=2,pending_study_starts=6,
    submission_receipt_sha256=sha(root/'submission_receipt.json'),source_pilot_sha256=wave.PILOT_SHA,
    source_helper_sha256=sha(driver),all_eight_registered_G2_chains_have_handles=True,
    study_mutated=False,recurring_monitor_created=False,physical_parameters_changed=False,
    job_completion_is_not_NEB_convergence=True)
with (root/'queued_jobs_verification.json').open('x') as f:
    json.dump(proof,f,indent=2);f.write('\n')
files = ['preflight.json','submission_receipt.json','submission_journal.jsonl',
         'executed_queue_training_waves.py','queued_jobs_verification.json']
for report in preflight['reports']:
    assert report['channel'] in wave.CHANNELS
    for item in report['endpoint_caches']:
        rel = f"cache_preflight/{report['channel']}/image_{item['image_index']:04d}/seed_cache_audit.json"
        assert item['image_index'] in (0,8) and sha(root/rel) == item['audit_sha256']
        files.append(rel)
assert len(files) == 13 and all((root/p).is_file() for p in files)
with (root/'audit_export.tar').open('xb') as out:
    with tarfile.open(fileobj=out,mode='w') as bundle:
        for p in files:bundle.add(root/p,arcname=p,recursive=False)
print(json.dumps(dict(check_CST=proof['check_CST'],active=2,pending=6,
    jobs=[dict(channel=j['channel'],job_id=j['job_id'],state=j['actual_fields']['JobState']) for j in records],
    portable_export_files=len(files),export_sha256=sha(root/'audit_export.tar'),licensed_assets_exported=False)))
PY
# Protect the heredoc terminator from Windows stdin's trailing CRLF.
