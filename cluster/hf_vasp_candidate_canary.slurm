#!/bin/bash
# Full-SCF test of one preserved candidate image rejected by the lattice gate.
# The result is diagnostic only and is never inserted into a production cache.
#SBATCH -J vasp-candidate-canary
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
: "${SOURCE_REPO:?}" "${CANDIDATE_TRAJECTORY:?}" "${SOURCE_INPUT:?}" \
  "${CANARY_WORKDIR:?}" "${N_IMAGES:?}" "${IMAGE_INDEX:?}" "${VASP_SYMPREC:?}"
[[ ! -e "${CANARY_WORKDIR}" ]] || { echo "canary workdir exists" >&2; exit 3; }

python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${SOURCE_REPO}"
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${SOURCE_INPUT}" --trajectory "${CANDIDATE_TRAJECTORY}" \
  --n-images "${N_IMAGES}" --trajectory-index "${IMAGE_INDEX}" \
  --workdir "${CANARY_WORKDIR}" --symprec "${VASP_SYMPREC}" \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" \
  --timeout 1800 --full-scf --run
