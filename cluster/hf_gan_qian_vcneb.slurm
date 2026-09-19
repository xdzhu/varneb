#!/bin/bash
# 27 interior images; nine concurrent 32-MPI workers; no endpoint workers.
#SBATCH -J varneb-gan-qian-vcneb
#SBATCH -p hfacnormal01
#SBATCH -N 3
#SBATCH --ntasks=291
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 291 && "${SLURM_JOB_NUM_NODES}" == 3 ]] || exit 2
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite
python=/public/home/iai806/.conda/envs/icu/bin/python
case "${CASE_NAME}" in
  b4_b1_hexagonal) pressure=45.7 ;;
  b3_b1_diagonal) pressure=45.0 ;;
  *) exit 2 ;;
esac
case_root="${repo}/cases/${CASE_NAME}"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/input/initial" --final "${case_root}/input/final" \
  --initial-trajectory "${case_root}/initial.traj" --workdir "${case_root}/vcneb_n29_w9_mpi32" \
  --n-images 29 --fmax 0.10 --steps 300 --k 0.20 --pressure-gpa "${pressure}" \
  --optimizer FIRE --image-workers 9 --no-climb \
  --initial-static-summary "${case_root}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${case_root}/static_final/vasp_static_summary.json" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --minimum-distance 1.4 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${case_root}/vcneb_n29_w9_mpi32" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target 0.10
