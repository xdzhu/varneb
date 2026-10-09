#!/bin/bash
# E044: existing interiors only, output-only SCF and fixed-charge quadrature.
#SBATCH -J hfo2-path-berry-E044
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-13%2
#SBATCH --time=00:25:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set immutable tested archive}
batch_root=${BATCH_ROOT:?set audited fresh property batch}
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
/public/home/iai806/.conda/envs/icu/bin/python -m scripts.hfo2_switching_path_polarization \
  run --repository "${source_root}" --root "${batch_root}" --index "${SLURM_ARRAY_TASK_ID}"
