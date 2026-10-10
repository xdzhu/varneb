#!/bin/bash
# One read-only snapshot, not a monitor, submission or raw-snapshot TS audit.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
from datetime import datetime
import hashlib,json,subprocess
from pathlib import Path
root=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1')
out=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-remaining-E056-r1/live_jobs_snapshot.json')
if out.exists():raise FileExistsError('fresh observation only')
production=root/'source-fixed/cluster/hf_hfo2_clamped_chain_resume_20261010.slurm'
script_sha=hashlib.sha256(production.read_bytes()).hexdigest()
assert script_sha=='79d9990b42a8da90fe915550332d8060761d0ec108551035e284b69f7768ef02'
jobs={'PO_flip_T_pattern_preserving':'28661019','PO_to_M':'28661020'}
cmd=['sacct','-j',','.join(jobs.values()),'-X','--noheader','--parsable2','--format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,Partition,NodeList']
actual=subprocess.run(cmd,check=True,text=True,capture_output=True,timeout=30).stdout
rows={r.split('|')[0]:r.strip().split('|') for r in actual.splitlines() if r.strip()}
assert set(rows)==set(jobs.values())
reports=[]
for channel,job in jobs.items():
 band=root/channel/'band'
 log=band/'vcneb.opt.log';body=log.read_text();parsed=[]
 for line in body.splitlines():
  fields=line.split()
  if len(fields)==5 and fields[1].isdigit():
   parsed.append(dict(step=int(fields[1]),time_CST=fields[2],highest_E_eV_cell=float(fields[3]),fmax_eV_A=float(fields[4])))
 if not parsed:raise ValueError('no complete optimizer row')
 audits=[]
 for i in range(1,8):
  for p in sorted((band/f'image_{i:04d}').glob('scf_*/call_audit.json')):
   a=json.loads(p.read_text());audits.append(float(a['elapsed_seconds']))
 latest=parsed[-1];remaining=max(0,20-latest['step'])
 mean=sum(audits)/len(audits) if audits else None
 failure=band/'vcneb_failure.json'
 reports.append(dict(channel=channel,job_id=job,scheduler_row=rows[job],latest_complete_log_row=latest,
  raw_optimizer_log_sha256_at_read=hashlib.sha256(body.encode()).hexdigest(),tail=body.splitlines()[-6:],
  completed_new_interior_transport_records=len(audits),recorded_mean_SCF_seconds=mean,
  estimated_remaining_segment_hours_if_recent_speed=remaining*7*mean/3600 if mean else None,
  estimate_is_segment_cap_not_convergence=True,ordinary_threshold=.10,
  current_force_log_is_not_full_raw_snapshot_audit=True,failure_report_path=str(failure),failure_report_present=failure.exists()))
report=dict(check_CST=datetime.now().isoformat(),source_root=str(root),production_script_sha256=script_sha,
 scheduler_command=cmd,scheduler_raw_rows=actual,jobs=reports,
 jobs_mutated_or_submitted=False,live_source_overwritten=False,new_DFT_calls=0,holdout_generated_or_read=False)
with out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(report))
PY
# Preserve production runs and earlier observations.
