#!/bin/bash
# Eighteen interior images; six concurrent 32-MPI workers; three batches/evaluation.
#SBATCH -J hfo2-prl-vcneb
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=194
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 194 && "${SLURM_JOB_NUM_NODES}" == 2 ]] || exit 2
case "${CASE_NAME:-}" in t_to_po|po_to_m) ;; *) exit 2 ;; esac

repo=/public/home/iai806/abacus/agent-runs/20260919-prl2023-hfo2-fig2a
python=/public/home/iai806/.conda/envs/icu/bin/python
case_root="${repo}/cases/${CASE_NAME}"
initial="${case_root}/endpoint_relax/initial"
final="${case_root}/endpoint_relax/final"
path="${case_root}/initial_corrected.traj"
workdir="${case_root}/vcneb_n20_fmax010_no_ci"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-hfo2-standard-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"

"${python}" examples/run_vcneb_vasp.py \
  --initial "${initial}" --final "${final}" --initial-trajectory "${path}" \
  --workdir "${workdir}" --n-images 20 --fmax 0.10 --steps 300 --k 0.10 \
  --pressure-gpa 0 --optimizer FIRE --image-workers 6 --image-retries 0 --no-climb \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --minimum-distance 1.40 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${workdir}" \
  --max-min-distance 1.40 --max-deformation 0.50 --fmax-target 0.10
