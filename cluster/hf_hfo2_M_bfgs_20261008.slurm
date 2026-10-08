#!/bin/bash
# New M endpoint only; T/PO endpoints are not reoptimized.
#SBATCH -J hfo2-M-bfgs
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set immutable source archive root}
endpoint_root=${ENDPOINT_ROOT:?set fresh endpoint workdir}
parameters=${PARAMETERS:?set reviewed immutable source-directory JSON}
[[ ! -e "${endpoint_root}" ]] || { echo "refusing existing endpoint run" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
/public/home/iai806/.conda/envs/icu/bin/python -m scripts.relax_ase_endpoint \
  --structure benchmarks/hfo2_channels/20261008/author_path_registration/M_seed_registered.vasp \
  --workdir "${endpoint_root}" --factory examples.hfo2_fixed_input_factory:make_factory \
  --parameters "${parameters}" \
  --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
  --pressure-gpa 0 --fmax 0.03 --stress-kbar 2 --maxstep 0.02 --steps 80
