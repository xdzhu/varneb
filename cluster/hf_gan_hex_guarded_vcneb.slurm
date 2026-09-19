#!/bin/bash
# Immutable standalone source; do not overwrite the active B3 source tree.
#SBATCH -J gan-hex-guarded-vcneb
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
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-lattice-guard
original=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite
case_root="${original}/cases/b4_b1_hexagonal"
work="${repo}/vcneb_hex_n29_w9_mpi32"
python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/input/initial" --final "${case_root}/input/final" \
  --resume --resume-trajectory "${case_root}/vcneb_n29_w9_mpi32/vcneb.traj" \
  --workdir "${work}" --n-images 29 --fmax 0.10 --steps 298 --k 0.20 --pressure-gpa 45.7 \
  --optimizer FIRE --image-workers 9 --no-climb \
  --initial-static-summary "${case_root}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${case_root}/static_final/vasp_static_summary.json" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --native-lattice-probe "${original}/diagnostics/hex_image10_basis/vasp_lattice_probe" \
  --candidate-step-retries 8 --candidate-step-retry-factor 0.5 \
  --minimum-distance 1.4 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${work}" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target 0.10
