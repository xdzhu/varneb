#!/bin/bash
# One immutable endpoint static SCF under the exact VASP path contract.
#SBATCH -J vasp-endpoint-static
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
: "${SOURCE_REPO:?}" "${INITIAL_DIR:?}" "${FINAL_DIR:?}" "${WORKDIR:?}" \
  "${N_IMAGES:?}" "${STATIC_ENDPOINT:?}" "${VASP_SYMPREC:?}" \
  "${MIN_DISTANCE:?}" "${MAX_DEFORMATION:?}"
case "${STATIC_ENDPOINT}" in initial|final) ;; *) exit 2 ;; esac
[[ ! -e "${WORKDIR}" ]] || { echo "endpoint workdir exists" >&2; exit 3; }

python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${SOURCE_REPO}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${INITIAL_DIR}" --final "${FINAL_DIR}" \
  --workdir "${WORKDIR}" --n-images "${N_IMAGES}" \
  --static-only --static-endpoint "${STATIC_ENDPOINT}" \
  --vasp-isym -1 --vasp-symprec "${VASP_SYMPREC}" \
  --minimum-distance "${MIN_DISTANCE}" --maximum-deformation "${MAX_DEFORMATION}" \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
