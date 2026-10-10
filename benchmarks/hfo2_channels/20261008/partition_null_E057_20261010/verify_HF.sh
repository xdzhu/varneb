#!/bin/bash
# Fresh zero-DFT archive replay, never mutate live production directories.
set -euo pipefail
root=/public/home/iai806/abacus/agent-runs/20261010-varneb-partition-null-E057-r1
archive=/public/home/iai806/abacus/agent-runs/varneb-E057-test-source-23c7b3a.tar
[[ ! -e "$root" ]] || exit 2
[[ $(sha256sum "$archive" | cut -d' ' -f1) == f782806fce4cb3b62a8945a77916d585917390b043cc043ab53e61e0a620b4fc ]] || exit 2
mkdir -p "$root/source"
tar -xf "$archive" -C "$root/source"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$root/source"
cd "$root/source"
/public/home/iai806/.conda/envs/icu/bin/python - "$root" "$archive" <<'PY'
from datetime import datetime
import hashlib,json,shutil,subprocess,sys,tarfile
from pathlib import Path
import numpy as np
import ase
from scripts.check_response_partition_invariance import run_check
root,archive=map(Path,sys.argv[1:]);source=Path.cwd()
def forbidden(*args,**kwargs):raise RuntimeError('No external executable/DFT/scheduler launch in this check')
subprocess.run=forbidden
result=run_check()
with (root/'analytic_HF.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
sources=dict(result['source_SHA256'])
for name in ('scripts/check_nested_response.py','scripts/relax_clamped_ase_endpoint.py'):
 sources[name]=hashlib.sha256((source/name).read_bytes()).hexdigest()
with tarfile.open(archive) as bundle:
 for name,expected in sources.items():
  assert hashlib.sha256(bundle.extractfile(name).read()).hexdigest()==expected
  target=root/'executed_source'/name;target.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(source/name,target)
receipt=dict(check_CST=datetime.now().isoformat(),tested_staged_tree='23c7b3a0040aa529109b39f38bc18bdf4efffdd9',
 tested_archive_SHA256=hashlib.sha256(archive.read_bytes()).hexdigest(),python=sys.version,
 numpy=np.__version__,ASE=ase.__version__,execution_source=str(source),
 source_SHA256=sources,external_executable_prohibited=True,HF_pytest_claimed=False,
 environment_installed_or_upgraded=False,new_DFT_calls=0,new_job_submissions=0,
 production_sources_or_jobs_mutated=False,held_out_condition_read=False,
 analytic_stationary_points=result['analytic']['stationary_points'],
 analytic_paired_points=result['analytic']['paired_points'],
 maximum_errors=result['analytic']['maximum_errors'],full_Hessian_budget_registered=False)
with (root/'verification_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps(receipt))
with tarfile.open(root/'results.tar','w') as output:
 for path in (root/'analytic_HF.json',root/'verification_receipt.json',root/'executed_source'):
  output.add(path,arcname=path.relative_to(root))
print('results_tar_SHA256='+hashlib.sha256((root/'results.tar').read_bytes()).hexdigest())
PY
# A missing observation is not grounds to repeat or overwrite this namespace.
