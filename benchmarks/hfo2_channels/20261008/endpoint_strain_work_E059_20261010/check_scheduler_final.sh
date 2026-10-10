#!/bin/bash
# Single authoritative four-handle observation, no mutation or monitor.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
from datetime import datetime
import json,subprocess
from pathlib import Path
root=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-endpoint-work-E059-r1')
jobs=['28661019','28661020','28692775','28692776']
argv=['sacct','-X','-n','-j',','.join(jobs),'-P',
 '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,ElapsedRaw']
raw=subprocess.run(argv,check=True,text=True,capture_output=True,timeout=30).stdout
rows=[line.split('|') for line in raw.splitlines() if line.strip()]
assert {r[0] for r in rows}==set(jobs)
queue=subprocess.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C|%R'],
 check=True,text=True,capture_output=True,timeout=30).stdout
receipt=dict(check_CST=datetime.now().isoformat(),jobs=jobs,accounting_command=argv,
 accounting_raw=raw,accounting_rows=rows,queue_raw=queue,job_or_dependency_changed=False,
 new_submissions=0,recurring_monitor_created=False)
with (root/'final_scheduler_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps(receipt))
PY
# Retain the existing parent and child handles even during observation failure.
