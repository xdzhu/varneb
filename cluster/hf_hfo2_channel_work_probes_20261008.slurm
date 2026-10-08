#!/bin/bash
# Eight immutable geometry probes under the historical 100-Ry contract.
#SBATCH -J hfo2-channel-work
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-7%3
#SBATCH --time=00:30:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
probe_root=${PROBE_ROOT:?set fresh prepared probe root}
source_root=${SOURCE_ROOT:?set immutable source archive root}
case_id=$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
inputs=${probe_root}/points/${case_id}
workdir=${probe_root}/calculations/${case_id}
[[ ! -e "${workdir}" ]] || { echo "refusing existing workdir" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p "${workdir}"
for name in INPUT KPT STRU Hf.upf O.upf Hf_gga_10au_100Ry_4s2p2d1f.orb O_gga_10au_100Ry_2s2p1d.orb; do
  cp -p "${inputs}/${name}" "${workdir}/${name}"
done
cd "${workdir}"
sha256sum INPUT KPT STRU *.upf *.orb > input_sha256.txt
mpirun -np 32 /bin/bash -c 'allowed=$(awk "/Cpus_allowed_list/ {print \$2}" /proc/self/status); printf "rank=%s host=%s allowed=%s\n" "${PMI_RANK:-unset}" "$(hostname)" "${allowed}"' > mpi_affinity.txt
date --iso-8601=seconds > started_at.txt
/usr/bin/time -p mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus > abacus.out 2> abacus.err
date --iso-8601=seconds > finished_at.txt
/public/home/iai806/.conda/envs/icu/bin/python -m scripts.score_hfo2_channel_work_probes \
  --root "${probe_root}" --index "${SLURM_ARRAY_TASK_ID}"
