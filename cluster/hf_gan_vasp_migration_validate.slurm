#!/bin/bash
# Fresh environment validation only. No endpoint relaxation or CI.
#SBATCH -J varneb-gan-migrate
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/slurm-%x-%j.err
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf
python=/public/home/iai806/.conda/envs/icu/bin/python
case_root=${repo}/gan_b4_b1_vasp_pbe_paw_qian
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export VCNEB_GIT_REVISION=d8af06dbcb72283ae618abb8fe94459749084d9f-dirty-hf-migration
cd "${repo}"
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${case_root}/vcneb_input/initial" \
  --structures tests/fixtures/vasp_bravais/gan_image15.vasp tests/fixtures/vasp_bravais/gan_image24.vasp \
  --workdir "${case_root}/hf_validation/failing_inputs_full_scf" \
  --symprec 1e-4 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" --timeout 600 --full-scf --run
"${python}" scripts/validate_vasp_migration_gate.py \
  --reference "${case_root}/contract_full_scf_symprec1e4_20260918/probe_summary.json" \
  --candidate "${case_root}/hf_validation/failing_inputs_full_scf/probe_summary.json" \
  --output "${case_root}/hf_validation/migration_gate.json"
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${case_root}/vcneb_input/initial" \
  --trajectory "${case_root}/recovery_isym_minus1_symprec1e4_20260918/vcneb_b4_to_b1_tetragonal_n29_cu17_serial/vcneb.traj" \
  --n-images 29 --workdir "${case_root}/hf_validation/complete_chain_initialization" \
  --symprec 1e-4 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" --run
for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
    --workdir "${case_root}/hf_static_${endpoint}" --n-images 29 --fmax 0.10 --steps 0 \
    --k 0.20 --pressure-gpa 45.7 --image-workers 0 --no-climb --mic \
    --cell-interpolation linear --mapping auto --align-translation \
    --minimum-distance 1.4 --maximum-deformation 0.50 --static-only --static-endpoint "${endpoint}" \
    --vasp-isym -1 --vasp-symprec 1e-4 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
done
