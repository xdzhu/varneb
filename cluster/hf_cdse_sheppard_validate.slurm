#!/bin/bash
#SBATCH -J cdse-path-check
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
[[ "${ROUTE}" == cell_mapping || "${ROUTE}" == atomic_mapping ]] || exit 2
repo=${SUITE_REPO:-/public/home/iai806/abacus/agent-runs/20260918-varneb-cdse-sheppard}
python=/public/home/iai806/.conda/envs/icu/bin/python
probe=${LATTICE_PROBE:-/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite/diagnostics/hex_image10_basis/vasp_lattice_probe}
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"
candidate_args=()
if [[ -n "${ENDPOINT_CANDIDATES:-}" ]]; then
  candidate_args=(--endpoint-candidates "${ENDPOINT_CANDIDATES}")
fi
"${python}" scripts/prepare_cdse_sheppard_suite.py --stage path --suite-root cases --route "${ROUTE}" "${candidate_args[@]}"
root="${repo}/cases/${ROUTE}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${root}/input/initial" --final "${root}/input/final" \
  --initial-trajectory "${root}/initial.traj" --workdir "${root}/preflight" \
  --n-images 17 --pressure-gpa 0.0 --no-climb --validate-only \
  --native-lattice-probe "${probe}" --minimum-distance 1.8 --maximum-deformation 0.8 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${root}/input/initial" --trajectory "${root}/initial.traj" --n-images 17 \
  --workdir "${root}/initialization_probe" --symprec 1e-4 --ncores 32 \
  --vasp-bin "${VASP_ROOT}/bin/vasp_std" --timeout 180 --run
for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${root}/input/initial" --final "${root}/input/final" \
    --initial-trajectory "${root}/initial.traj" --workdir "${root}/static_${endpoint}" \
    --n-images 17 --pressure-gpa 0.0 --no-climb --static-only --static-endpoint "${endpoint}" \
    --native-lattice-probe "${probe}" --minimum-distance 1.8 --maximum-deformation 0.8 \
    --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
done
"${python}" scripts/verify_cdse_endpoint_statics.py --case-root "${root}"
