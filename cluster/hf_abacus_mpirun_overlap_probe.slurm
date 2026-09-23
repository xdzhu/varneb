#!/bin/bash
# Validate two concurrent Intel-MPI mpirun workers in one Slurm allocation.
# Each Hydra bootstrap proxy should own four distinct CPUs.
#SBATCH -J varneb-aba-overlap
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=4
#SBATCH --time=00:05:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }
workdir=${WORKDIR:?set a fresh WORKDIR}
[[ ! -e "${workdir}" ]] || { echo "refusing existing WORKDIR" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
export OMP_NUM_THREADS=1
mkdir -p "${workdir}"
printf 'start %s\n' "$(date --iso-8601=seconds)" > "${workdir}/timing.txt"
pids=()
for worker in A B; do
  (
    printf 'start %s %s\n' "${worker}" "$(date --iso-8601=seconds)" >> "${workdir}/timing.txt"
    mpirun -np 4 /bin/bash -c 'printf "rank=%s host=%s " "${PMI_RANK:-unset}" "$(hostname)"; grep Cpus_allowed_list /proc/self/status; sleep 20' \
      > "${workdir}/${worker}.out" 2> "${workdir}/${worker}.err"
    printf 'done %s %s\n' "${worker}" "$(date --iso-8601=seconds)" >> "${workdir}/timing.txt"
  ) &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "${pid}"; done
printf 'end %s\n' "$(date --iso-8601=seconds)" >> "${workdir}/timing.txt"
