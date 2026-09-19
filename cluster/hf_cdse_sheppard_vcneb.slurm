#!/bin/bash
#SBATCH -J cdse-vcneb
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
[[ "${ROUTE}" == cell_mapping || "${ROUTE}" == atomic_mapping ]] || exit 2
repo=${SUITE_REPO:-/public/home/iai806/abacus/agent-runs/20260918-varneb-cdse-sheppard}
python=/public/home/iai806/.conda/envs/icu/bin/python
probe=${LATTICE_PROBE:-/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite/diagnostics/hex_image10_basis/vasp_lattice_probe}
root="${repo}/cases/${ROUTE}"
work="${root}/vcneb_n17_w6_mpi32"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"
"${python}" -c 'import json,sys; assert json.load(open(sys.argv[1]))["all_passed"] is True' "${root}/endpoint_static_gate.json"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${root}/input/initial" --final "${root}/input/final" \
  --initial-trajectory "${root}/initial.traj" --workdir "${work}" \
  --n-images 17 --fmax 0.10 --steps 300 --k 0.20 --pressure-gpa 0.0 \
  --optimizer FIRE --image-workers 6 --no-climb \
  --initial-static-summary "${root}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${root}/static_final/vasp_static_summary.json" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --native-lattice-probe "${probe}" --candidate-step-retries 8 --candidate-step-retry-factor 0.5 \
  --minimum-distance 1.8 --maximum-deformation 0.8 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${work}" \
  --max-min-distance 1.8 --max-deformation 0.8 --fmax-target 0.10
