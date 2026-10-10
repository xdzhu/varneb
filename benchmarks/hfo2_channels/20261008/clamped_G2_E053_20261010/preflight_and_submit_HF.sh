set -euo pipefail
runtime=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1
source_root="$runtime/source"
seed_base="$runtime/seeds/strain_0000"
[[ -d "$source_root" && ! -e "$runtime/submission_handles.txt" ]] || exit 2
[[ "$(sha256sum "$runtime/source.tar" | cut -d' ' -f1)" == e78042abe89ee43083874f45048c87ae2aac77e2523f1b81c80b3f0574348809 ]] || exit 2
[[ "$(sha256sum "$runtime/seeds.tar" | cut -d' ' -f1)" == 19649f95d406137cbdd78aa9d1a6439f8adb8c32f1829c260f88b5192ae20584 ]] || exit 2
[[ "${FULL_REGRESSION_PASSED:-0}" == 1 ]] || exit 2
module purge
# Site env.sh reads SLURM_JOBID even during a login-node, zero-DFT preflight.
# Leave scheduler variables truthful/unset; do not fabricate an allocation.
set +u
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
set -u
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$source_root"
cd "$source_root"
bash -n cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
import numpy as np
from ase.io import read
import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import sha256
from vcneb.periodic_path import validate_periodic_path_lift
runtime=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1')
def no_external_launch(*a,**k):
    raise RuntimeError('zero-DFT preflight prohibits external executable launch')
transport.subprocess.run=no_external_launch
reports=[]
for channel in ['PO_flip_T_pattern_preserving','PO_to_M']:
    root=runtime/'seeds/strain_0000'/channel
    m=json.loads((root/'manifest.json').read_text())
    assert m['channel']==channel and m['strain']==m['pressure_GPa']==0
    assert m['physical_contract_sha256']==transport.CONTRACT and not m['climb'] and not m['holdout_generated']
    assert m['n_total_images']==9 and m['n_internal_images']==7 and m['n_fixed_endpoints']==2
    assert all(sha256(root/n)==h for n,h in m['files_sha256'].items())
    images=read(root/'seed.traj',index=':')
    assert len(images)==9 and all(a.calc is None for a in images)
    validate_periodic_path_lift(images)
    p=json.loads((root/'factory_parameters.json').read_text())
    factory=transport.make_clamped_seed_cached_factory(parameters=p,command='mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
    energies=[]
    for i,a in enumerate(images):
        assert a.get_chemical_symbols()==['Hf']*4+['O']*8
        # numpy.testing lazily launches a CPU-feature probe on this HF NumPy.
        # Keep the zero-external-launch guard intact, use the same numeric test.
        assert np.allclose(a.cell.array[:2],np.asarray(m['reference_cell_A'])[:2],atol=1e-10,rtol=0)
        a.calc=factory(i,a,runtime/'cache_preflight'/channel/f'image_{i:04d}')
        if i in (0,8):
            e=a.get_potential_energy(); f=a.get_forces(); s=a.get_stress()
            assert np.isfinite(e) and f.shape==(12,3) and s.shape==(6,) and a.calc.next_call==0
            energies.append(e)
        else:
            assert a.calc.results=={}
    reports.append({'channel':channel,'manifest_sha256':sha256(root/'manifest.json'),'new_DFT_calls':0,'full_physical_bytes_and_pinned_raw_hashes_checked':True,'cached_clamped_endpoint_energies_eV_cell':energies,'fresh_interior_images':7,'external_launch_prohibited':True})
output=runtime/'cache_preflight.json'
assert not output.exists()
output.write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports))
PY
for channel in PO_flip_T_pattern_preserving PO_to_M; do
 seed_root="$seed_base/$channel"
 /public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
  --initial "$seed_root/initial.vasp" --final "$seed_root/final.vasp" \
  --initial-chain "$seed_root/seed.traj" --n-images 9 --workdir "$runtime/geometry_preflight/$channel" \
  --clamped-plane-reference "$seed_root/substrate.vasp" --clamped-allow-tilt true --cell-scale 5.12968067458423 \
  --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 10 --maxstep .02 --k .2 \
  --no-climb --mapping identity --mic --no-align-cells --cell-interpolation linear \
  --minimum-distance 1.6 --maximum-deformation .25 --maximum-cell-step .02 --candidate-step-retries 4 \
  --require-continuous-periodic-lift --summary-format brief --validate-only
done
# Unknown/pending/live study handles are never repeated. Other user jobs untouched.
study_count=$(squeue -h -u iai806 -o '%j' | awk '$1 ~ /^hfo2/ {n++} END {print n+0}')
[[ "$study_count" == 0 ]] || { echo 'prior study allocation still present; do not submit' >&2; exit 2; }
for channel in PO_flip_T_pattern_preserving PO_to_M; do
 [[ ! -e "$runtime/$channel" ]] || exit 2
 job=$(sbatch --parsable --job-name="hfo2-G2-${channel}" \
  --export="ALL,RUN_DFT=1,SOURCE_ROOT=$source_root,SEED_ROOT=$seed_base/$channel,WORKDIR=$runtime/$channel/band" \
  --output="$runtime/${channel}.slurm.out" --error="$runtime/${channel}.slurm.err" \
  cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm)
 [[ "$job" =~ ^[0-9]+$ ]] || { echo 'ambiguous submission; reconcile scheduler, do not repeat' >&2; exit 2; }
 printf '%s|%s\n' "$channel" "$job" | tee -a "$runtime/submission_handles.txt"
done
date +%FT%T%z
squeue -h -u iai806 -o '%i|%j|%T|%M|%C|%R'
# End two explicitly reviewed first G2 segments; no automatic continuation.
