#!/bin/bash
# One actual read-only scheduler receipt; no retry, mutation or recurring watch.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
from datetime import datetime
import hashlib,json,subprocess
from pathlib import Path
root=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E058-r1')
receipt=json.loads((root/'submission_receipt.json').read_text())
source=Path(json.loads((root/'preflight.json').read_text())['source_root'])
reports=[]
for j in receipt['jobs']:
 argv=['scontrol','show','job','-o',j['job_id']]
 raw=subprocess.run(argv,check=True,text=True,capture_output=True,timeout=30).stdout
 fields=dict(item.split('=',1) for item in raw.split() if '=' in item)
 assert fields['JobId']==j['job_id'] and fields['JobName']=='hfo2-G2-E058-'+j['channel']
 assert fields['JobState']=='PENDING' and fields['Reason']=='Dependency'
 assert fields['Dependency']=='afterok:28661019(unfulfilled),afterok:28661020(unfulfilled)'
 assert fields['Partition']=='hfacnormal01' and fields['UserId'].startswith('iai806(')
 assert fields['NumCPUs']==fields['CPUs/Task']=='32' and fields['NumTasks']=='1'
 assert fields['NumNodes']=='1-1' and fields['TimeLimit']=='04:00:00'
 assert fields['Command']==str(source/'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm')
 reports.append(dict(channel=j['channel'],job_id=j['job_id'],command=argv,raw_scontrol=raw,
  actual_fields=fields,resource_and_dependency_checks_passed=True))
argv=['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C|%R']
raw=subprocess.run(argv,check=True,text=True,capture_output=True,timeout=30).stdout
rows=[[v.strip() for v in line.split('|')] for line in raw.splitlines() if line.strip()]
study=[r for r in rows if r[1].startswith('hfo2')]
assert {r[0] for r in study}=={'28661019','28661020','28692775','28692776'}
assert sum(r[2] in ('RUNNING','COMPLETING') for r in study)==2
report=dict(check_CST=datetime.now().isoformat(),status='actual_two_pending_dependency_jobs_verified',
 jobs=reports,queue_command=argv,queue_raw=raw,active_study_allocations=2,pending_study_starts=2,
 new_DFT_jobs_running=0,study_mutated=False,submission_receipt_sha256=hashlib.sha256(
 (root/'submission_receipt.json').read_bytes()).hexdigest(),
 job_completion_is_not_NEB_convergence=True,recurring_monitor_created=False)
with (root/'queued_jobs_verification.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(dict(status=report['status'],check_CST=report['check_CST'],
 active_study_allocations=2,pending_study_starts=2)))
PY
# Existing production trees and unrelated jobs remain unchanged.
