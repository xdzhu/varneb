#!/bin/bash
# Two serial fixed-contract SCFs. No endpoints, band restart, CI or G2.
#SBATCH -J hfo2-G1-reversing-peaks
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
batch_root=${BATCH_ROOT:?set reviewed two-point root}
source_root=${SOURCE_ROOT:?set immutable tested source}
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
mpirun -np 32 /bin/bash -c 'printf "rank=%s host=%s\n" "${PMI_RANK:-unset}" "$(hostname)"' > "${batch_root}/mpi_affinity.txt"
for index in 0 1; do
  /public/home/iai806/.conda/envs/icu/bin/python -m scripts.prepare_hfo2_reversing_peak_sampling \
    run-point --root "${batch_root}" --index "${index}"
done
/public/home/iai806/.conda/envs/icu/bin/python -m scripts.prepare_hfo2_reversing_peak_sampling summarize --root "${batch_root}"
