#!/bin/bash
# Inert: no submission; admit only after an original allocation completes.
#SBATCH -J hfo2-G3-E166-frozen
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
 && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?fresh pinned source required}
runtime_sha=${RUNTIME_MANIFEST_SHA256:?exact Python manifest required}
config=${FROZEN_CONFIG:?exact frozen-control config required}
config_sha=${FROZEN_CONFIG_SHA256:?exact config SHA required}
workdir=${WORKDIR:?fresh frozen_R output required}
[[ ! -e "${workdir}" ]] || exit 2
/public/home/iai806/.conda/envs/icu/bin/python -c '
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1]);p=root/"runtime_files.json"
if hashlib.sha256(p.read_bytes()).hexdigest()!=sys.argv[2]: raise ValueError("runtime manifest drift")
files=json.loads(p.read_bytes())
if not files or any(pathlib.PurePosixPath(n).is_absolute() or ".." in pathlib.PurePosixPath(n).parts or not n.endswith(".py") or hashlib.sha256((root/n).read_bytes()).hexdigest()!=h for n,h in files.items()): raise ValueError("runtime byte drift")
' "${source_root}" "${runtime_sha}"
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${source_root}"
cd "${source_root}"
/public/home/iai806/.conda/envs/icu/bin/python \
 benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/worker.py run \
 --config "${config}" --config-sha256 "${config_sha}" --workdir "${workdir}"
