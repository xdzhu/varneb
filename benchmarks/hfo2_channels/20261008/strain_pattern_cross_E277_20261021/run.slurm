#!/bin/bash
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && -n "${BUNDLE_SHA256:-}" && "${CONDITION:-}" =~ ^(zero_q1|plus1_q0)$ ]] || exit 2
root=/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1
runtime=/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1
[[ ! -e "${root}/results/${CONDITION}" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${runtime}/source"
cd "${runtime}/source"
/public/home/iai806/.conda/envs/icu/bin/python "${root}/bundle/worker.py" run --bundle-pin "${BUNDLE_SHA256}" --condition "${CONDITION}"
