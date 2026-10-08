#!/bin/bash
# Eight independent fixed-cell curvature holdouts, no fitting after DFT.
#SBATCH -J hfo2-reduce-test
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-7%2
#SBATCH --time=00:30:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
batch_root=${BATCH_ROOT:?set prepared holdouts}
source_root=${SOURCE_ROOT:?set immutable source archive}
case_id=$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
inputs=${batch_root}/points/${case_id}
workdir=${batch_root}/calculations/${case_id}
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
mpirun -np 32 /bin/bash -c 'printf "rank=%s host=%s\n" "${PMI_RANK:-unset}" "$(hostname)"' > mpi_affinity.txt
/usr/bin/time -p mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus > abacus.out 2> abacus.err
/public/home/iai806/.conda/envs/icu/bin/python -m scripts.score_hfo2_reduction_holdout \
  --root "${batch_root}" --index "${SLURM_ARRAY_TASK_ID}"
