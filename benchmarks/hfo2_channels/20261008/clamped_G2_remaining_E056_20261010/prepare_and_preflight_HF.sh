#!/bin/bash
# Six remaining preregistered training seeds; no allocation or DFT launch.
set -euo pipefail
root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-remaining-E056-r1
source_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1/source
archive=/public/home/iai806/abacus/agent-runs/varneb-E055-clean-86d2637.tar
[[ ! -e "$root" && -d "$source_root" ]] || exit 2
[[ $(sha256sum "$archive" | cut -d' ' -f1) == cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c ]] || exit 2
module purge
set +u
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
set -u
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$source_root"
mkdir -p "$root"
cd "$source_root"
bash -n cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm
/public/home/iai806/.conda/envs/icu/bin/python - "$root" "$source_root" "$archive" <<'PY'
from datetime import datetime
import hashlib,json,shutil,sys,tarfile
from pathlib import Path
import numpy as np
import ase
from ase.io import read
import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_chains import prepare,CHANNEL_FINAL
from vcneb.material_runner import main as geometry_preflight
from vcneb.periodic_path import validate_periodic_path_lift
root,source,archive=map(Path,sys.argv[1:])
# Verify the existing tested runtime rather than replacing its source.
checked={}
with tarfile.open(archive) as bundle:
 for entry in bundle.getmembers():
  if entry.isfile() and entry.name.endswith('.py') and entry.name.startswith(('vcneb/','scripts/','examples/')):
   raw=bundle.extractfile(entry).read()
   if raw!=(source/entry.name).read_bytes():raise ValueError('immutable tested source changed: '+entry.name)
   checked[entry.name]=hashlib.sha256(raw).hexdigest()
assert checked and sha256(source/'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm')=='7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed'
# Initialize NumPy's documented test helper before the strict executable guard;
# its CPU feature discovery is not an electronic-structure calculation.
np.testing.assert_allclose([1.],[1.])
def no_external_launch(*a,**k):raise RuntimeError('zero-DFT preparation prohibits external executable launch')
transport.subprocess.run=no_external_launch
channels=[('strain_0000','PO_to_T'),('strain_0000','PO_flip_T_pattern_reversing')]
channels += [('strain_p0100',c) for c in CHANNEL_FINAL]
reports=[]
case=source/'benchmarks/hfo2_channels/20261008'
command='mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus'
for condition,channel in channels:
 seed=root/'seeds'/condition/channel
 m=prepare(case,condition,channel,seed)
 assert m['physical_contract_sha256']==transport.CONTRACT and m['strain'] in (0.,.01)
 assert m['new_DFT_calls']==0 and not m['holdout_generated'] and not m['climb']
 assert all(sha256(seed/n)==h for n,h in m['files_sha256'].items())
 images=read(seed/'seed.traj',index=':')
 assert len(images)==9 and all(a.calc is None for a in images)
 validate_periodic_path_lift(images)
 p=json.loads((seed/'factory_parameters.json').read_text())
 factory=transport.make_clamped_seed_cached_factory(parameters=p,command=command)
 endpoints=[]
 for i,a in enumerate(images):
  assert a.get_chemical_symbols()==['Hf']*4+['O']*8
  assert np.allclose(a.cell.array[:2],np.asarray(m['reference_cell_A'])[:2],atol=1e-10,rtol=0.)
  directory=root/'cache_preflight'/condition/channel/f'image_{i:04d}'
  a.calc=factory(i,a,directory)
  if i in (0,8):
   energy=a.get_potential_energy();forces=a.get_forces();stress=a.get_stress()
   assert np.isfinite(energy) and forces.shape==(12,3) and stress.shape==(6,) and a.calc.next_call==0
   endpoints.append(dict(image_index=i,energy_eV_cell=energy,
    cache_audit_sha256=sha256(directory/'seed_cache_audit.json'),
    actual_raw_source=p['seed_cache_records'][i]['directory'],new_DFT_calls=0))
  else:assert a.calc.results=={} and a.calc.next_call==0
 work=root/'geometry_preflight'/condition/channel
 geometry_preflight(['--initial',str(seed/'initial.vasp'),'--final',str(seed/'final.vasp'),
  '--initial-chain',str(seed/'seed.traj'),'--n-images','9','--workdir',str(work),
  '--factory','examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory',
  '--parameters',str(seed/'factory_parameters.json'),'--command',command,
  '--clamped-plane-reference',str(seed/'substrate.vasp'),'--clamped-allow-tilt','true','--cell-scale',str(m['cell_scale_A']),
  '--pressure-gpa','0','--cell-mode','full','--optimizer','FIRE','--fmax','.10','--steps','10','--maxstep','.02','--k','.2',
  '--no-climb','--mapping','identity','--mic','--no-align-cells','--cell-interpolation','linear',
  '--minimum-distance','1.6','--maximum-deformation','.25','--maximum-cell-step','.02','--candidate-step-retries','4',
  '--require-continuous-periodic-lift','--summary-format','brief','--validate-only'])
 g=json.loads((work/'vcneb_preflight.json').read_text())
 assert g['mechanical_boundary']['kind']=='clamped_plane' and g['mechanical_boundary']['cell_dofs']==3
 assert not g['climbing_image_requested'] and g['requires_stress'] and g['n_interior_images']==7
 assert g['calculator_validation']=='not_instantiated_validate_only'
 shutil.copyfile(work/'vcneb_preflight.json',seed/'vcneb_preflight_HF.json')
 reports.append(dict(condition=condition,channel=channel,strain=m['strain'],manifest_sha256=sha256(seed/'manifest.json'),
  geometry_preflight_sha256=sha256(seed/'vcneb_preflight_HF.json'),cached_endpoints=endpoints,
  all_six_original_physical_bytes_and_native_raw_EFS_checked=True,fresh_internal_images=7,
  DFT_calls=0,ordinary_threshold=.10,CI=False,readiness='seed_and_cache_preflight_passed_NOT_submitted'))
# A small exact source subset supplements the full tested archive identity.
for relative in ['scripts/prepare_hfo2_clamped_chains.py','cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm']:
 target=root/'executed_source'/relative;target.parent.mkdir(parents=True,exist_ok=True)
 shutil.copyfile(source/relative,target)
receipt=dict(status='six_remaining_registered_training_seeds_ready_NOT_submitted',check_CST=datetime.now().isoformat(),
 tested_source_commit='86d2637ebcd6c502df43b91ca03977fd50120716',tested_source_archive_sha256=sha256(archive),
 tested_source_full_local_regression_passed=1459,tested_source_full_local_skipped=2,
 execution_source_directory=str(source),execution_source_code_files_byte_checked=len(checked),
 execution_source_code_sha256=checked,python=sys.version,ASE_version=ase.__version__,numpy_version=np.__version__,
 prepared_channels=reports,prepared_chains=6,total_registered_G2_matrix_chains=8,
 existing_independent_G2_chains=2,new_independent_DFT_chains_started=0,new_DFT_calls=0,
 physical_parameters_changed=False,old_source_or_jobs_mutated=False,
 held_out_condition_generated_or_labels_read=False,external_executable_launch_prohibited_during_preparation=True,
 limitations=['Geometry starts and exact terminal EFS caches are not optimized paths, barriers, TS or forecasts',
  'Explicit scheduler reconciliation and study concurrency<=2 required before any actual submission',
  'Original ten endpoint representations are reused, no relaxation rerun and no new material/reference matrix'])
with (root/'preparation_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
for p in root.rglob('*'):
 if p.is_file() and p.suffix in ('.upf','.orb'):raise ValueError('licensed asset should not enter the evidence bundle')
print(json.dumps({k:receipt[k] for k in ('status','check_CST','prepared_chains','new_DFT_calls')}))
PY
tar -cf "$root/results.tar" -C "$root" seeds cache_preflight geometry_preflight executed_source preparation_receipt.json
sha256sum "$root/results.tar"
# No submission, monitor, live production update, DFT or held-out labels.
