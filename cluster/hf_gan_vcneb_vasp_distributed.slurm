#!/bin/bash
# 27 interiors, nine 32-MPI workers, three waves per chain evaluation.
# Each 97-core node has 96 worker cores and one spare core for manager overhead.
#SBATCH -J varneb-gan-vcneb
#SBATCH -p hfacnormal01
#SBATCH -N 3
#SBATCH --ntasks=291
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/slurm-%x-%j.err
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
[[ "${SLURM_NTASKS}" -ge 291 && "${SLURM_JOB_NUM_NODES}" == 3 ]] || exit 2
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf
python=/public/home/iai806/.conda/envs/icu/bin/python
case_root=${repo}/gan_b4_b1_vasp_pbe_paw_qian
workdir=${case_root}/vcneb_hf_n29_w9_mpi32
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export VCNEB_GIT_REVISION=d8af06dbcb72283ae618abb8fe94459749084d9f-dirty-hf-migration
cd "${repo}"
"${python}" scripts/validate_vasp_migration_gate.py \
  --reference "${case_root}/contract_full_scf_symprec1e4_20260918/probe_summary.json" \
  --candidate "${case_root}/hf_validation/failing_inputs_full_scf/probe_summary.json" \
  --output "${case_root}/hf_validation/production_migration_gate.json"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
  --workdir "${workdir}" --n-images 29 --fmax 0.10 --steps 300 --k 0.20 \
  --pressure-gpa 45.7 --optimizer FIRE --image-workers 9 --no-climb \
  --initial-static-summary "${case_root}/hf_static_initial/vasp_static_summary.json" \
  --final-static-summary "${case_root}/hf_static_final/vasp_static_summary.json" \
  --resume --resume-trajectory "${case_root}/recovery_isym_minus1_symprec1e4_20260918/vcneb_b4_to_b1_tetragonal_n29_cu17_serial/vcneb.traj" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --mic --cell-interpolation linear --mapping auto --align-translation \
  --minimum-distance 1.4 --maximum-deformation 0.50 --vasp-isym -1 --vasp-symprec 1e-4 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${workdir}" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target 0.10
