#!/bin/bash
# Fresh-source zero-DFT replay; never submit, install or modify a production tree.
set -euo pipefail
root=/public/home/iai806/abacus/agent-runs/20261010-varneb-endpoint-work-E059-r1
archive=/public/home/iai806/abacus/agent-runs/varneb-E059-tested.tar
[[ ! -e "$root" ]] || { echo 'fresh replay namespace required' >&2; exit 2; }
[[ $(sha256sum "$archive" | cut -d' ' -f1) == 7c9d7ffad8d6df2d36515b46d209e88b5419be21de8b9a3a93a9e542aaa4c793 ]] || exit 2
mkdir -p "$root/source"
tar -xf "$archive" -C "$root/source"
export PYTHONPATH="$root/source" PYTHONDONTWRITEBYTECODE=1
cd "$root/source"
/public/home/iai806/.conda/envs/icu/bin/python - "$root" "$archive" <<'PY'
from datetime import datetime
import hashlib,json,subprocess,sys,tarfile
from pathlib import Path
import numpy as np
import ase
from scripts.analyze_hfo2_endpoint_strain_work import analyse
from scripts.audit_hfo2_static_replica import sha256,audited_results
from examples.hfo2_fixed_input_factory import CONTRACT
root,archive=map(Path,sys.argv[1:]);source=root/'source'
names=['vcneb/strain_work.py','scripts/analyze_hfo2_endpoint_strain_work.py',
 'scripts/prepare_hfo2_clamped_chains.py','scripts/audit_hfo2_static_replica.py',
 'examples/hfo2_fixed_input_factory.py','vcneb/abacus.py',
 'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md','paper/VARNEB_JCTC/METHODS_DRAFT.md']
with tarfile.open(archive) as bundle:
 for name in names:
  assert (source/name).read_bytes()==bundle.extractfile(name).read()
def no_launch(*a,**k):raise RuntimeError('zero-DFT replay prohibits executable launch')
subprocess.run=no_launch
case=source/'benchmarks/hfo2_channels/20261008'
result=analyse(case)
reference=json.loads((case/'endpoint_strain_work_E059_20261010/endpoint_work.json').read_text())
def numeric_flat(value):
 if isinstance(value,dict):return sum((numeric_flat(v) for v in value.values()),[])
 if isinstance(value,list):return sum((numeric_flat(v) for v in value),[])
 if isinstance(value,(float,int)) and not isinstance(value,bool):return [float(value)]
 return []
assert result['status']==reference['status']
delta=float(np.max(np.abs(np.asarray(numeric_flat(result))-numeric_flat({k:v for k,v in reference.items() if k!='source_sha256'}))))
assert delta<1e-9
physical=[]
for phase,endpoints in result['endpoints'].items():
 for endpoint in endpoints:
  cache=endpoint['source_audit']['cache'];directory=Path(cache['directory'])
  hashes={n:sha256(directory/n) for n in cache['input_sha256']}
  assert hashes==cache['input_sha256'] and all(hashes[n]==h for n,h in CONTRACT.items())
  assert sha256(directory/'OUT.ABACUS/running_scf.log')==cache['raw_log_sha256']
  actual=audited_results(directory)
  archived=audited_results(case/endpoint['terminal_export_relative'])
  assert all(np.allclose(actual[k],archived[k],atol=1e-12,rtol=0) for k in actual)
  physical.append(dict(phase=phase,strain=endpoint['strain'],raw_directory=str(directory),
   original_six_physical_bytes_and_STRU_log_checked=True,native_32MPI_raw_EFS_checked=True))
with (root/'endpoint_work_HF.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
receipt=dict(status='original_HF_environment_and_ten_raw_endpoints_replay_passed',
 check_CST=datetime.now().isoformat(),tested_staged_tree='76dd34f6812bb6cb4dda8cef3f323837b4fade4b',
 archive_sha256=sha256(archive),source_files_byte_checked=len(names),
 source_sha256={n:sha256(source/n) for n in names},python=sys.version,ASE_version=ase.__version__,
 numpy_version=np.__version__,maximum_numeric_replay_difference=delta,numeric_replay_tolerance=1e-9,
 numeric_tolerance_is_NOT_DFT_error_bound=True,actual_terminal_checks=physical,
 new_DFT_calls=0,external_executable_launch_prohibited=True,HF_pytest_run=False,
 environment_installed_or_upgraded=False,production_source_or_job_changed=False,
 physical_parameters_changed=False,holdout_generated_or_read=False,
 material_barrier_or_prediction_advantage_certified=False)
with (root/'replay_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({k:receipt[k] for k in ('status','check_CST','source_files_byte_checked','maximum_numeric_replay_difference','new_DFT_calls')}))
PY
# Keep all earlier raw endpoints, live bands and receipts read-only.
